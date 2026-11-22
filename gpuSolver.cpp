#include "gpuSolver.h"
#include "cellPhysics.h"

#include <cuda_runtime.h>

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <execution>
#include <functional>
#include <map>
#include <numeric>
#include <ranges>
#include <stdexcept>
#include <utility>

namespace {

typedef DESCRIPTOR<T> D;
typedef CellPhysics<T, DESCRIPTOR> Physics;
typedef AtomicAcceleratedLattice3D<T, DESCRIPTOR> Block;
typedef GpuSolver::VelocitySource Source;
typedef Array<T, 3> Vec;

bool isVelocityTag(int model)
{
    return model == CollisionModel::Boundary_RegularizedVelocity_0_1__BGK
        || model == CollisionModel::Boundary_RegularizedVelocity_0_M1__BGK
        || model == CollisionModel::Boundary_RegularizedVelocity_1_1__BGK
        || model == CollisionModel::Boundary_RegularizedVelocity_1_M1__BGK
        || model == CollisionModel::Boundary_RegularizedVelocity_2_1__BGK
        || model == CollisionModel::Boundary_RegularizedVelocity_2_M1__BGK;
}

bool isPressureTag(int model)
{
    return model == CollisionModel::Boundary_RegularizedDensity_0_1__TRT
        || model == CollisionModel::Boundary_RegularizedDensity_0_M1__TRT
        || model == CollisionModel::Boundary_RegularizedDensity_1_1__TRT
        || model == CollisionModel::Boundary_RegularizedDensity_1_M1__TRT
        || model == CollisionModel::Boundary_RegularizedDensity_2_1__TRT
        || model == CollisionModel::Boundary_RegularizedDensity_2_M1__TRT;
}

bool isSupportedTag(int model)
{
    return model == CollisionModel::NoDynamics || model == CollisionModel::BounceBack
        || model == CollisionModel::BGK || isVelocityTag(model) || isPressureTag(model);
}

// HemoFlow's collision for each placeholder tag (gpuSolver.h). The relaxation rate comes from
// the omega function: constant (Newtonian) or local (Carreau-Yasuda).
template <class OmegaFunction>
struct HemoCollision {
    T omega;   // Stored omega of the bulk dynamics (Newtonian, or Carreau at infinite shear)
    OmegaFunction omegaOf;

    // Opening nodes collide like fluid: on the CPU, the velocity nodes have no automatic
    // preparation and the pressure nodes an empty one; both are rebuilt after streaming.
    void operator()(int model, Array<T, D::numPop>& f, Array<T, D::ExternalField::numScalars>& ext,
                    Array<T, GPUconst<T, DESCRIPTOR>::maxStaticScalars>, T*, plint) const
    {
        if (model == CollisionModel::BGK || isVelocityTag(model) || isPressureTag(model)) {
            T rhoBar;
            Physics::collide(f, &ext[0], omega, omegaOf, rhoBar);
        }
        else if (model == CollisionModel::BounceBack) {
            for (int iPop = 1; iPop <= D::q / 2; ++iPop)
                std::swap(f[iPop], f[iPop + D::q / 2]);
        }
        // NoDynamics: nothing
    }
};

void pullExternal(Block const* block, plint cell, T* ext)
{
    for (int k = 0; k < D::ExternalField::numScalars; ++k)
        ext[k] = block->pullExt(cell, k);
}

// Velocity of a pressure node from its known populations (DensityDirichletBoundaryDynamics::computeVelocity)
template <int direction, int orientation>
Vec pressureNodeVelocity(Block const* block, plint cell, T rhoBar)
{
    Array<T, D::q> f;
    block->pullPop(cell, f);
    T ext[D::ExternalField::numScalars];
    pullExternal(block, cell, ext);
    Vec u;
    Physics::template pressureNodeVelocity<direction, orientation>(f, ext, rhoBar, u);
    return u;
}

// Velocity of a neighbour of a pressure node, as its CPU dynamics computes it
Vec sourceVelocity(Block const* block, T const* scalars, Vec const* buffer, Source const& source)
{
    Vec u((T)0, (T)0, (T)0);
    if (source.kind == Source::FLUID) {
        Array<T, D::q> f;
        block->pullPop(source.cell, f);
        T ext[D::ExternalField::numScalars];
        pullExternal(block, source.cell, ext);
        Physics::forcedVelocity(f, ext, u);
    }
    else if (source.kind == Source::VELOCITY_NODE) {
        u = Vec(scalars[source.aux], scalars[source.aux + 1], scalars[source.aux + 2]);
    }
    else if (source.kind == Source::PRESSURE_NODE) {
        u = buffer[source.aux];
    }
    return u;   // ZERO: walls and outside cells
}

// Rebuild the velocity nodes of one opening in one block from the imposed velocity
template <int direction, int orientation, typename Node>
void completeVelocity(Block* block, T const* scalars, Node const* nodes, std::size_t n, T omega)
{
    std::for_each(std::execution::par_unseq, nodes, nodes + n, [=](Node const& node) {
        Array<T, D::q> f;
        block->pullPop(node.cell, f);
        T ext[D::ExternalField::numScalars];
        pullExternal(block, node.cell, ext);
        Vec u(scalars[node.index], scalars[node.index + 1], scalars[node.index + 2]);
        Physics::template completeVelocityNode<direction, orientation>(f, ext, u, omega);
        block->pushPop(node.cell, f);
    });
}

// Finite-difference pressure nodes of one opening in one block, pass 1: their velocities
template <int direction, int orientation, typename Node>
void computePressureVelocities(Block* block, T const* scalars, Node const* nodes, std::size_t n, Vec* buffer)
{
    std::for_each(std::execution::par_unseq, nodes, nodes + n, [=](Node const& node) {
        buffer[node.buffer] = pressureNodeVelocity<direction, orientation>(block, node.cell, scalars[node.index]);
    });
}

// Pass 2: rebuild the nodes from the velocity gradient (PlaneFdBoundaryFunctional3D::processCell)
template <int direction, int orientation, typename Node>
void completePressure(Block* block, T const* scalars, Node const* nodes, std::size_t n, Vec const* buffer, T omega)
{
    std::for_each(std::execution::par_unseq, nodes, nodes + n, [=](Node const& node) {
        if (!node.rebuild)
            return;
        Vec u0 = buffer[node.buffer];
        Vec uInward = sourceVelocity(block, scalars, buffer, node.inward);
        Vec uPlus[3], uMinus[3];
        for (int a = 0; a < 3; ++a) {
            uPlus[a] = sourceVelocity(block, scalars, buffer, node.plus[a]);
            uMinus[a] = sourceVelocity(block, scalars, buffer, node.minus[a]);
        }
        Vec grad[3];
        Physics::template pressureNodeGradient<direction, orientation>(u0, uInward, uPlus, uMinus, grad);
        Array<T, D::q> f;
        block->pullPop(node.cell, f);
        T rho = D::fullRho(scalars[node.index]);
        Physics::completePressureNode(f, rho, u0, grad, omega);
        block->pushPop(node.cell, f);
    });
}

// Porous force in the porous fluid cells of one block (PorousForceFunctional): from the velocity
// after streaming, with the force of the previous step in the half-force correction
template <typename PorousCell>
void porousForce(Block* block, PorousCell const* cells, std::size_t n, T linCoeff, T quadCoeff)
{
    std::for_each(std::execution::par_unseq, cells, cells + n, [=](PorousCell const& c) {
        Array<T, D::q> f;
        block->pullPop(c.cell, f);
        T ext[D::ExternalField::numScalars];
        pullExternal(block, c.cell, ext);
        Vec u;
        Physics::forcedVelocity(f, ext, u);
        Physics::porousForce(u, c.mLin, c.mQuad, linCoeff, quadCoeff, ext);
        for (int iD = 0; iD < D::d; ++iD)
            block->pushExt(c.cell, D::ExternalField::forceBeginsAt + iD, ext[D::ExternalField::forceBeginsAt + iD]);
    });
}

// All populations of envelope copies of opening nodes from their owner cells
void copyCells(Block const* from, Block* to, GpuSolver::CellPair const* pairs, std::size_t n)
{
    std::for_each(std::execution::par_unseq, pairs, pairs + n, [=](GpuSolver::CellPair const& p) {
        Array<T, D::q> f;
        from->pullPop(p.from, f);
        to->pushPop(p.to, f);
    });
}

// The same between processes: the populations of owner cells into a message buffer, and from
// a received buffer into the envelope copies
void packCells(GpuSolver::RemoteCell const* cells, std::size_t n, T* buffer)
{
    auto indices = std::views::iota(std::size_t(0), n);
    std::for_each(std::execution::par_unseq, indices.begin(), indices.end(), [=](std::size_t j) {
        Array<T, D::q> f;
        cells[j].block->pullPop(cells[j].cell, f);
        for (int k = 0; k < D::q; ++k)
            buffer[j * D::q + k] = f[k];
    });
}

void unpackCells(GpuSolver::RemoteCell const* cells, std::size_t n, T const* buffer)
{
    auto indices = std::views::iota(std::size_t(0), n);
    std::for_each(std::execution::par_unseq, indices.begin(), indices.end(), [=](std::size_t j) {
        Array<T, D::q> f;
        for (int k = 0; k < D::q; ++k)
            f[k] = buffer[j * D::q + k];
        cells[j].block->pushPop(cells[j].cell, f);
    });
}

// Face index 0..5 of (direction, orientation), for the dispatch to the templates
int faceIndex(int direction, int orientation)
{
    return direction * 2 + (orientation > 0 ? 1 : 0);
}

}  // namespace


bool queryGpu(std::string& name, std::size_t& freeBytes, std::size_t& totalBytes)
{
    int device = 0;
    cudaDeviceProp properties;
    if (cudaGetDevice(&device) != cudaSuccess || cudaGetDeviceProperties(&properties, device) != cudaSuccess)
        return false;
    if (cudaMemGetInfo(&freeBytes, &totalBytes) != cudaSuccess)
        return false;
    name = properties.name;
    return true;
}

std::vector<GpuDevice> visibleGpus()
{
    std::vector<GpuDevice> devices;
    int count = 0;
    if (cudaGetDeviceCount(&count) != cudaSuccess)
        return devices;
    for (int d = 0; d < count; ++d) {
        cudaDeviceProp properties;
        if (cudaGetDeviceProperties(&properties, d) != cudaSuccess) {
            devices.push_back(GpuDevice {"unknown-" + std::to_string(d), "unknown"});   // Keeps the indices
            continue;
        }
        char uuid[33], bus[32];
        for (int k = 0; k < 16; ++k)
            std::snprintf(uuid + 2 * k, 3, "%02x", (unsigned char)properties.uuid.bytes[k]);
        std::snprintf(bus, sizeof(bus), "%04x:%02x:%02x.0", properties.pciDomainID, properties.pciBusID, properties.pciDeviceID);
        devices.push_back(GpuDevice {uuid, bus});
    }
    return devices;
}

bool selectGpu(int device)
{
    return cudaSetDevice(device) == cudaSuccess;
}

void releaseGpu()
{
    cudaDeviceReset();
}

std::size_t gpuBytesPerCell()
{
    return 2 * D::numPop * sizeof(T) + D::ExternalField::numScalars * sizeof(T) + sizeof(int) + sizeof(std::uint32_t)
         + sizeof(plint);
}

std::size_t gpuBytesPerPorousCell()
{
    return sizeof(GpuSolver::PorousCell);
}


GpuSolver::GpuSolver(MultiBlockLattice3D<T, DESCRIPTOR> const& placeholder, std::vector<OpeningHandler*> const& openings,
                     SimPar const& sim, GpuPorousMedium const& porous)
    : lattice(nullptr), linCoeff(porous.linCoeff), quadCoeff(porous.quadCoeff), omega(sim.omega),
      nonNewtonian(sim.nonNewtonian), carreau(sim.carreau), domain(placeholder.getBoundingBox())
{
    // The copy would keep pointers to the placeholder's data processors and the blocks they modify
    if (!placeholder.getStoredProcessors().empty())
        throw std::runtime_error("GPU lattice: the placeholder lattice must not have data processors.");
    lattice = new AcceleratedLattice3D<T, DESCRIPTOR>(placeholder);
    auto fail = [this](std::string const& message) {
        freeMessageBuffers();
        delete lattice;
        lattice = nullptr;
        throw std::runtime_error("GPU lattice: " + message);
    };

    // Every cell must carry a tag that the collision functor handles
    for (auto const& entry : lattice->getBlockLattices()) {
        Block const& block = *entry.second;
        int const* models = block.getCollisionMatrix();
        for (plint i = 0; i < block.getN(); ++i)
            if (!isSupportedTag(models[i]))
                fail("a cell has the collision model '" + toCollisionString(models[i]) +
                     "', which the GPU collision of HemoFlow does not handle.");
    }

    // Opening nodes per opening and block, envelope copies included
    std::map<std::pair<plint, plint>, plint> pressureBuffer;   // (block id, cell) -> index in the velocity buffer
    std::vector<std::vector<Index3D> > pressureLocal;          // Local coordinates of the pressure nodes, per group
    for (std::size_t o = 0; o < openings.size(); ++o) {
        OpeningHandler const& opening = *openings[o];
        bool velocity = opening.isVelocityType();
        if (!velocity && opening.getOpeningType() != OPENING_PRESSURE)
            fail("opening " + std::to_string(opening.getLabel()) + " has a type that the GPU does not handle.");
        for (auto const& entry : lattice->getBlockLattices()) {
            Block& block = *entry.second;
            Dot3D location = block.getLocation();
            plint nx = block.getNx(), ny = block.getNy(), nz = block.getNz();
            NodeGroup<VelocityNode> vGroup;
            NodeGroup<PressureNode> pGroup;
            std::vector<Index3D> local;
            for (Index3D const& node : opening.getNodes()) {
                plint iX = node.x - location.x, iY = node.y - location.y, iZ = node.z - location.z;
                if (iX < 0 || iY < 0 || iZ < 0 || iX >= nx || iY >= ny || iZ >= nz)
                    continue;
                plint i = iZ + nz * (iY + ny * iX);
                int model = block.getCollisionMatrix()[i];
                T* scalars = nullptr;
                plint index = -1;
                if (velocity) {
                    vec3d u;
                    if (!opening.getUnitProfile(node, u))
                        continue;
                    if (!isVelocityTag(model))
                        fail("a node of opening " + std::to_string(opening.getLabel()) + " has no velocity boundary condition.");
                    block.getDynamicScalar(i, scalars, index);
                    vGroup.dynamicScalars = scalars;
                    vGroup.nodes.push_back(VelocityNode {i, index, u.x, u.y, u.z});
                }
                else {
                    if (!isPressureTag(model))
                        fail("a node of opening " + std::to_string(opening.getLabel()) + " has no pressure boundary condition.");
                    block.getDynamicScalar(i, scalars, index);
                    pGroup.dynamicScalars = scalars;
                    plint buffer = (plint)pressureBuffer.size();
                    pressureBuffer[std::make_pair(entry.first, i)] = buffer;
                    PressureNode p;
                    p.cell = i;
                    p.index = index;
                    p.buffer = buffer;
                    pGroup.nodes.push_back(p);
                    local.push_back(Index3D {(int)iX, (int)iY, (int)iZ});
                }
            }
            if (!vGroup.nodes.empty()) {
                vGroup.opening = (int)o;
                vGroup.direction = opening.getFaceAxis();
                vGroup.orientation = -opening.getFaceSign();
                vGroup.block = &block;
                velocityNodes.push_back(std::move(vGroup));
            }
            if (!pGroup.nodes.empty()) {
                pGroup.opening = (int)o;
                pGroup.direction = opening.getFaceAxis();
                pGroup.orientation = -opening.getFaceSign();
                pGroup.block = &block;
                pressureNodes.push_back(std::move(pGroup));
                pressureLocal.push_back(std::move(local));
            }
        }
    }
    pressureVelocities.resize(pressureBuffer.size());

    // Neighbours of the pressure nodes: how their velocity is obtained
    for (std::size_t g = 0; g < pressureNodes.size(); ++g) {
        NodeGroup<PressureNode>& group = pressureNodes[g];
        Block& block = *group.block;
        plint id = -1;
        for (auto const& entry : lattice->getBlockLattices())
            if (entry.second == group.block)
                id = entry.first;
        plint n[3] = {block.getNx(), block.getNy(), block.getNz()};
        int const* models = block.getCollisionMatrix();
        auto source = [&](Index3D const& c, int axis, int shift, bool& inside) {
            plint p[3] = {c.x, c.y, c.z};
            p[axis] += shift;
            Source s {Source::ZERO, -1, -1};
            if (p[0] < 0 || p[1] < 0 || p[2] < 0 || p[0] >= n[0] || p[1] >= n[1] || p[2] >= n[2]) {
                inside = false;
                return s;
            }
            plint j = p[2] + n[2] * (p[1] + n[1] * p[0]);
            int model = models[j];
            if (model == CollisionModel::BGK) {
                s = Source {Source::FLUID, j, -1};
            }
            else if (isVelocityTag(model)) {
                T* scalars = nullptr;
                plint index = -1;
                block.getDynamicScalar(j, scalars, index);
                s = Source {Source::VELOCITY_NODE, j, index};
            }
            else if (isPressureTag(model)) {
                auto it = pressureBuffer.find(std::make_pair(id, j));
                if (it == pressureBuffer.end())
                    fail("a pressure node has a neighbour with a pressure boundary condition that belongs to no opening.");
                s = Source {Source::PRESSURE_NODE, j, it->second};
            }
            return s;
        };
        for (std::size_t k = 0; k < group.nodes.size(); ++k) {
            PressureNode& node = group.nodes[k];
            Index3D const& c = pressureLocal[g][k];
            bool inside = true;
            node.inward = source(c, group.direction, -group.orientation, inside);
            for (int a = 0; a < 3; ++a) {
                if (a == group.direction) {
                    node.plus[a] = node.minus[a] = Source {Source::ZERO, -1, -1};
                    continue;
                }
                node.plus[a] = source(c, a, 1, inside);
                node.minus[a] = source(c, a, -1, inside);
            }
            node.rebuild = inside ? 1 : 0;
        }
    }

    // Envelope copies of the opening nodes: every block whose domain (bulk and envelope) contains
    // a node, other than the block whose bulk contains it (the owner), gets the owner's
    // populations after each exchange. The pairs come from the global block structure in the
    // same order on all processes, so that the cells of a message between two processes match
    // one by one. None with a single block.
    {
        MultiBlockManagement3D const& management = lattice->getMultiBlockManagement();
        std::map<plint, Box3D> const& bulks = management.getSparseBlockStructure().getBulks();
        ThreadAttribution const& attribution = management.getThreadAttribution();
        int const me = global::mpi().getRank();
        auto localIndex = [this](plint id, Index3D const& node, Block*& block, plint& i) {
            block = lattice->getBlockLattices().at(id);
            Dot3D location = block->getLocation();
            plint iX = node.x - location.x, iY = node.y - location.y, iZ = node.z - location.z;
            if (iX < 0 || iY < 0 || iZ < 0 || iX >= block->getNx() || iY >= block->getNy() || iZ >= block->getNz())
                throw std::runtime_error("an opening node lies outside the block that should hold it.");
            i = iZ + block->getNz() * (iY + block->getNy() * iX);
        };
        std::map<std::pair<plint, plint>, std::size_t> copyGroup;   // (owner id, copy id) -> nodeCopies index
        std::map<int, std::size_t> sendIndex, receiveIndex;          // Other process -> nodeSends/nodeReceives index
        try {
            for (OpeningHandler const* opening : openings) {
                for (Index3D const& node : opening->getNodes()) {
                    plint owner = -1;
                    for (auto const& bulk : bulks)
                        if (contained(node.x, node.y, node.z, bulk.second)) {
                            owner = bulk.first;
                            break;
                        }
                    if (owner < 0)
                        continue;
                    int ownerRank = attribution.getMpiProcess(owner);
                    for (auto const& bulk : bulks) {
                        plint id = bulk.first;
                        if (id == owner || !contained(node.x, node.y, node.z, management.getEnvelope(id)))
                            continue;
                        int copyRank = attribution.getMpiProcess(id);
                        Block *ownerBlock = nullptr, *copyBlock = nullptr;
                        plint ownerCell = -1, copyCell = -1;
                        if (ownerRank == me && copyRank == me) {
                            localIndex(owner, node, ownerBlock, ownerCell);
                            localIndex(id, node, copyBlock, copyCell);
                            auto key = std::make_pair(owner, id);
                            auto it = copyGroup.find(key);
                            if (it == copyGroup.end()) {
                                it = copyGroup.insert(std::make_pair(key, nodeCopies.size())).first;
                                nodeCopies.push_back(EnvelopeCopies {ownerBlock, copyBlock, {}});
                            }
                            nodeCopies[it->second].cells.push_back(CellPair {ownerCell, copyCell});
                        }
                        else if (ownerRank == me) {
                            localIndex(owner, node, ownerBlock, ownerCell);
                            auto it = sendIndex.find(copyRank);
                            if (it == sendIndex.end()) {
                                it = sendIndex.insert(std::make_pair(copyRank, nodeSends.size())).first;
                                nodeSends.push_back(RemoteCopies {copyRank, {}, nullptr});
                            }
                            nodeSends[it->second].cells.push_back(RemoteCell {ownerBlock, ownerCell});
                        }
                        else if (copyRank == me) {
                            localIndex(id, node, copyBlock, copyCell);
                            auto it = receiveIndex.find(ownerRank);
                            if (it == receiveIndex.end()) {
                                it = receiveIndex.insert(std::make_pair(ownerRank, nodeReceives.size())).first;
                                nodeReceives.push_back(RemoteCopies {ownerRank, {}, nullptr});
                            }
                            nodeReceives[it->second].cells.push_back(RemoteCell {copyBlock, copyCell});
                        }
                    }
                }
            }
            for (std::vector<RemoteCopies>* list : {&nodeSends, &nodeReceives})
                for (RemoteCopies& r : *list)
                    if (cudaMallocHost((void**)&r.buffer, sizeof(T) * D::q * r.cells.size()) != cudaSuccess) {
                        r.buffer = nullptr;
                        throw std::runtime_error("no pinned host memory for the opening-node messages.");
                    }
        }
        catch (std::runtime_error const& error) {
            fail(error.what());
        }
    }

    // Porous fluid cells: nonzero multipliers and the fluid tag. As on the CPU, opening nodes in
    // the porous medium get no force (PorousForceFunctional skips boundary dynamics); porous voxels
    // are fluid or opening voxels (GeometryFile::isPorous), so walls and outside cells do not occur.
    if (porous.multipliers) {
        MultiBlockManagement3D const& management = lattice->getMultiBlockManagement();
        MultiBlockManagement3D const& fieldManagement = porous.multipliers->getMultiBlockManagement();
        for (auto const& entry : lattice->getBlockLattices()) {
            Block& block = *entry.second;
            Box3D bulk = management.getBulk(entry.first), cells;
            if (!(fieldManagement.getBulk(entry.first) == bulk))
                fail("the multiplier field of the flow diverter has another block structure than the lattice.");
            if (!intersect(bulk, porous.box, cells))
                continue;
            TensorField3D<T, 2> const& field = porous.multipliers->getComponent(entry.first);
            Dot3D location = block.getLocation(), fieldLocation = field.getLocation();
            plint ny = block.getNy(), nz = block.getNz();
            int const* models = block.getCollisionMatrix();
            PorousGroup group;
            group.block = &block;
            for (plint x = cells.x0; x <= cells.x1; ++x)
                for (plint y = cells.y0; y <= cells.y1; ++y)
                    for (plint z = cells.z0; z <= cells.z1; ++z) {
                        Array<T, 2> const& m = field.get(x - fieldLocation.x, y - fieldLocation.y, z - fieldLocation.z);
                        if (m[0] == 0.0 && m[1] == 0.0)
                            continue;
                        plint i = (z - location.z) + nz * ((y - location.y) + ny * (x - location.x));
                        if (models[i] == CollisionModel::BGK)
                            group.cells.push_back(PorousCell {i, m[0], m[1]});
                        else if (!isVelocityTag(models[i]) && !isPressureTag(models[i]))
                            fail("a porous voxel is neither fluid nor an opening node.");
                    }
            if (!group.cells.empty())
                porousCells.push_back(std::move(group));
        }
    }

    imposeOpenings(openings);
}

GpuSolver::~GpuSolver()
{
    freeMessageBuffers();
    delete lattice;
    if (staging != nullptr)
        cudaFreeHost(staging);
}

void GpuSolver::collideAndStream()
{
    if (nonNewtonian)
        lattice->externalCollideAndStream(HemoCollision<Physics::CarreauOmega> {omega, Physics::CarreauOmega {carreau}});
    else
        lattice->externalCollideAndStream(HemoCollision<Physics::ConstantOmega> {omega, Physics::ConstantOmega {omega}});
    completeVelocityNodes();
    completePressureNodes();
    applyPorousForce();
    lattice->executeInternalProcessors();   // Envelope exchange
    copyNodeEnvelopes();
    lattice->incrementTime();
}

void GpuSolver::completeVelocityNodes()
{
    for (NodeGroup<VelocityNode> const& group : velocityNodes) {
        VelocityNode const* nodes = group.nodes.data();
        std::size_t n = group.nodes.size();
        switch (faceIndex(group.direction, group.orientation)) {
        case 0: completeVelocity<0, -1>(group.block, group.dynamicScalars, nodes, n, omega); break;
        case 1: completeVelocity<0, 1>(group.block, group.dynamicScalars, nodes, n, omega); break;
        case 2: completeVelocity<1, -1>(group.block, group.dynamicScalars, nodes, n, omega); break;
        case 3: completeVelocity<1, 1>(group.block, group.dynamicScalars, nodes, n, omega); break;
        case 4: completeVelocity<2, -1>(group.block, group.dynamicScalars, nodes, n, omega); break;
        default: completeVelocity<2, 1>(group.block, group.dynamicScalars, nodes, n, omega); break;
        }
    }
}

void GpuSolver::completePressureNodes()
{
    Vec* buffer = pressureVelocities.data();
    for (NodeGroup<PressureNode> const& group : pressureNodes) {
        PressureNode const* nodes = group.nodes.data();
        std::size_t n = group.nodes.size();
        switch (faceIndex(group.direction, group.orientation)) {
        case 0: computePressureVelocities<0, -1>(group.block, group.dynamicScalars, nodes, n, buffer); break;
        case 1: computePressureVelocities<0, 1>(group.block, group.dynamicScalars, nodes, n, buffer); break;
        case 2: computePressureVelocities<1, -1>(group.block, group.dynamicScalars, nodes, n, buffer); break;
        case 3: computePressureVelocities<1, 1>(group.block, group.dynamicScalars, nodes, n, buffer); break;
        case 4: computePressureVelocities<2, -1>(group.block, group.dynamicScalars, nodes, n, buffer); break;
        default: computePressureVelocities<2, 1>(group.block, group.dynamicScalars, nodes, n, buffer); break;
        }
    }
    for (NodeGroup<PressureNode> const& group : pressureNodes) {
        PressureNode const* nodes = group.nodes.data();
        std::size_t n = group.nodes.size();
        switch (faceIndex(group.direction, group.orientation)) {
        case 0: completePressure<0, -1>(group.block, group.dynamicScalars, nodes, n, buffer, omega); break;
        case 1: completePressure<0, 1>(group.block, group.dynamicScalars, nodes, n, buffer, omega); break;
        case 2: completePressure<1, -1>(group.block, group.dynamicScalars, nodes, n, buffer, omega); break;
        case 3: completePressure<1, 1>(group.block, group.dynamicScalars, nodes, n, buffer, omega); break;
        case 4: completePressure<2, -1>(group.block, group.dynamicScalars, nodes, n, buffer, omega); break;
        default: completePressure<2, 1>(group.block, group.dynamicScalars, nodes, n, buffer, omega); break;
        }
    }
}

void GpuSolver::copyNodeEnvelopes()
{
    for (EnvelopeCopies const& group : nodeCopies)
        copyCells(group.from, group.to, group.cells.data(), group.cells.size());
    if (nodeSends.empty() && nodeReceives.empty())
        return;
    // Every pair of processes that shares opening nodes exchanges one message each way at most
    MPI_Comm communicator = global::mpi().getGlobalCommunicator();
    int const tag = 7101;
    std::vector<MPI_Request> requests(nodeSends.size() + nodeReceives.size());
    std::size_t r = 0;
    for (RemoteCopies& in : nodeReceives)
        MPI_Irecv(in.buffer, (int)(sizeof(T) * D::q * in.cells.size()), MPI_BYTE, in.rank, tag, communicator, &requests[r++]);
    for (RemoteCopies& out : nodeSends) {
        packCells(out.cells.data(), out.cells.size(), out.buffer);
        MPI_Isend(out.buffer, (int)(sizeof(T) * D::q * out.cells.size()), MPI_BYTE, out.rank, tag, communicator, &requests[r++]);
    }
    MPI_Waitall((int)requests.size(), requests.data(), MPI_STATUSES_IGNORE);
    for (RemoteCopies const& in : nodeReceives)
        unpackCells(in.cells.data(), in.cells.size(), in.buffer);
}

void GpuSolver::freeMessageBuffers()
{
    for (std::vector<RemoteCopies>* list : {&nodeSends, &nodeReceives})
        for (RemoteCopies& r : *list)
            if (r.buffer != nullptr) {
                cudaFreeHost(r.buffer);
                r.buffer = nullptr;
            }
}

void GpuSolver::applyPorousForce()
{
    for (PorousGroup const& group : porousCells)
        porousForce(group.block, group.cells.data(), group.cells.size(), linCoeff, quadCoeff);
}

// Velocity openings: u = profile * scale, rounded to T as ImposeOpeningValues3D does on the CPU.
// Pressure openings: rhoBar of the imposed density, as StoreDensityDynamics::defineDensity stores it.
void GpuSolver::imposeOpenings(std::vector<OpeningHandler*> const& openings)
{
    for (NodeGroup<VelocityNode> const& group : velocityNodes) {
        double scale = openings[group.opening]->getVelocityScale();
        T* scalars = group.dynamicScalars;
        VelocityNode const* nodes = group.nodes.data();
        std::for_each(std::execution::par_unseq, nodes, nodes + group.nodes.size(), [=](VelocityNode const& node) {
            scalars[node.index] = (T)(node.ux * scale);
            scalars[node.index + 1] = (T)(node.uy * scale);
            scalars[node.index + 2] = (T)(node.uz * scale);
        });
    }
    for (NodeGroup<PressureNode> const& group : pressureNodes) {
        T rhoBar = D::rhoBar(openings[group.opening]->getDensity());
        T* scalars = group.dynamicScalars;
        PressureNode const* nodes = group.nodes.data();
        std::for_each(std::execution::par_unseq, nodes, nodes + group.nodes.size(),
                      [=](PressureNode const& node) { scalars[node.index] = rhoBar; });
    }
}

// A GPU kernel writes the state into pinned host memory, chunk by chunk, and the host copies it
// into the CPU lattice. Reading the GPU lattice from the host instead would move its managed
// pages to the host, and the next time steps would fault them back one by one (about 7 s for
// NAP180 at 5 M voxels, against 0.1 s for the copy itself).
void GpuSolver::copyTo(MultiBlockLattice3D<T, DESCRIPTOR>& cpu) const
{
    constexpr int numExt = D::ExternalField::numScalars;
    constexpr int perCell = D::q + numExt;
    constexpr plint chunk = plint(1) << 20;   // Cells per chunk: 88 MB of staging in single precision
    if (staging == nullptr && cudaMallocHost((void**)&staging, sizeof(T) * perCell * chunk) != cudaSuccess) {
        staging = nullptr;
        throw std::runtime_error("GPU lattice: no pinned host memory for the copy to the CPU.");
    }
    T* out = staging;
    for (auto const& entry : lattice->getBlockLattices()) {
        Block const* gpuBlock = entry.second;
        BlockLattice3D<T, DESCRIPTOR>& cpuBlock = cpu.getComponent(entry.first);
        plint ny = gpuBlock->getNy(), nz = gpuBlock->getNz(), n = gpuBlock->getN();
        for (plint begin = 0; begin < n; begin += chunk) {
            plint count = std::min(chunk, n - begin);
            auto indices = std::views::iota(plint(0), count);
            // Component-major within the chunk, so that neighbouring threads write neighbouring words
            std::for_each(std::execution::par_unseq, indices.begin(), indices.end(), [=](plint j) {
                Array<T, D::q> f;
                gpuBlock->pullPop(begin + j, f);
                for (int k = 0; k < D::q; ++k)
                    out[k * chunk + j] = f[k];
                for (int k = 0; k < numExt; ++k)
                    out[(D::q + k) * chunk + j] = gpuBlock->pullExt(begin + j, k);
            });
            for (plint j = 0; j < count; ++j) {
                plint i = begin + j, iX = i / (ny * nz), iY = (i / nz) % ny, iZ = i % nz;
                Cell<T, DESCRIPTOR>& cell = cpuBlock.get(iX, iY, iZ);
                for (int k = 0; k < D::q; ++k)
                    cell[k] = out[k * chunk + j];
                for (int k = 0; k < numExt; ++k)
                    *cell.getExternal(k) = out[(D::q + k) * chunk + j];
            }
        }
    }
}

T GpuSolver::averageEnergy() const
{
    double sum = 0.0;
    MultiBlockManagement3D const& management = lattice->getMultiBlockManagement();
    for (auto const& entry : lattice->getBlockLattices()) {
        Block* block = entry.second;
        Dot3D location = block->getLocation();
        Box3D bulk = management.getBulk(entry.first).shift(-location.x, -location.y, -location.z);
        plint ny = block->getNy(), nz = block->getNz();
        int const* models = block->getCollisionMatrix();
        auto indices = std::views::iota(plint(0), block->getN());
        sum += std::transform_reduce(
            std::execution::par_unseq, indices.begin(), indices.end(), 0.0, std::plus<double>(), [=](plint i) -> double {
                plint iX = i / (ny * nz), remainder = i % (ny * nz), iY = remainder / nz, iZ = remainder % nz;
                if (iX < bulk.x0 || iX > bulk.x1 || iY < bulk.y0 || iY > bulk.y1 || iZ < bulk.z0 || iZ > bulk.z1)
                    return 0.0;
                int model = models[i];
                Vec u((T)0, (T)0, (T)0);
                if (model == CollisionModel::BGK) {
                    u = sourceVelocity(block, nullptr, nullptr, Source {Source::FLUID, i, -1});
                }
                else if (isVelocityTag(model) || isPressureTag(model)) {
                    T* scalars = nullptr;
                    plint index = -1;
                    block->getDynamicScalar(i, scalars, index);
                    switch (model) {
                    case CollisionModel::Boundary_RegularizedDensity_0_M1__TRT: u = pressureNodeVelocity<0, -1>(block, i, scalars[index]); break;
                    case CollisionModel::Boundary_RegularizedDensity_0_1__TRT: u = pressureNodeVelocity<0, 1>(block, i, scalars[index]); break;
                    case CollisionModel::Boundary_RegularizedDensity_1_M1__TRT: u = pressureNodeVelocity<1, -1>(block, i, scalars[index]); break;
                    case CollisionModel::Boundary_RegularizedDensity_1_1__TRT: u = pressureNodeVelocity<1, 1>(block, i, scalars[index]); break;
                    case CollisionModel::Boundary_RegularizedDensity_2_M1__TRT: u = pressureNodeVelocity<2, -1>(block, i, scalars[index]); break;
                    case CollisionModel::Boundary_RegularizedDensity_2_1__TRT: u = pressureNodeVelocity<2, 1>(block, i, scalars[index]); break;
                    default: u = Vec(scalars[index], scalars[index + 1], scalars[index + 2]); break;   // Velocity node
                    }
                }
                else {
                    return 0.0;
                }
                return (double)VectorTemplate<T, DESCRIPTOR>::normSqr(u);
            });
    }
    global::mpi().reduceAndBcast(sum, MPI_SUM);   // Several GPUs
    return (T)(sum / (double)domain.nCells());
}
