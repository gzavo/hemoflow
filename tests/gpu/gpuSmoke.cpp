// GPU smoke test (doc/dev/TODO-GPU.md, Phases 0 and 1); built with -DHEMOFLOW_GPU=ON as gpuSmoke.
//
// Checks the toolchain and the pattern of the GPU port on a small case: a closed box with
// bounce-back walls, a bounce-back sphere, a NoDynamics region, a density bump and a constant
// body force.
//  - Reference: a CPU lattice with the real dynamics of the chosen model.
//  - GPU: a CPU lattice with placeholder dynamics (BGK background, which Palabos can map to its
//    GPU models) is copied to an AcceleratedLattice3D and stepped with a HemoFlow collision
//    functor; the force comes from the external field.
// Models:
//  - bgk:     Guo-forced BGK (GuoExternalForceBGKdynamics on the CPU).
//  - hemo:    HemoFlow's Newtonian collision, complete regularized BGK with Guo forcing, from
//             cellPhysics.h (GuoExternalForceCompleteRegularizedBGKdynamics on the CPU).
//  - carreau: the Carreau-Yasuda collision of cellPhysics.h (GuoForcedCarreauDynamics on the CPU).
// After the given number of steps the populations and forces of all bulk cells are compared;
// the exit code is 1 if they differ by more than round-off. The timings show the cost of the
// blocks: every block and every block overlap costs synchronous GPU kernels per step.
//
// Usage: [mpirun -np N] gpuSmoke [n=64] [blocksPerAxis=1] [steps=200] [bgk|hemo|carreau]

#include "globals.h"
#include "nonNewtonian.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <string>

typedef DESCRIPTOR<T> D;
typedef CellPhysics<T, DESCRIPTOR> Physics;

void bounceBack(Array<T, D::numPop> &f)
{
    for (int iPop = 1; iPop <= D::q / 2; ++iPop)
        std::swap(f[iPop], f[iPop + D::q / 2]);
}

// Guo-forced BGK with the arithmetic of GuoExternalForceBGKdynamics::collide. The model id of
// a cell is only a tag set by the placeholder dynamics.
struct GuoBGKCollision {
    T omega;
    void operator()(int model, Array<T, D::numPop> &f, Array<T, D::ExternalField::numScalars> &ext,
                    Array<T, GPUconst<T, DESCRIPTOR>::maxStaticScalars>, T *, plint) const
    {
        if (model == CollisionModel::BGK) {
            T rhoBar;
            Array<T, 3> j, u;
            momentTemplatesImpl<T, D::BaseDescriptor>::get_rhoBar_j(f, rhoBar, j);
            T invRho = D::invRho(rhoBar);
            for (int iD = 0; iD < 3; ++iD)
                u[iD] = j[iD] * invRho + ext[D::ExternalField::forceBeginsAt + iD] / (T)2;
            T rho = D::fullRho(rhoBar);
            for (int iD = 0; iD < 3; ++iD)
                j[iD] = rho * u[iD];
            dynamicsTemplatesImpl<T, D::BaseDescriptor>::bgk_ma2_collision(f, rhoBar, j, omega);
            externalForceTemplatesImpl<T, D>::addGuoForce(f, &ext[0], u, omega, (T)1);
        }
        else if (model == CollisionModel::BounceBack) {
            bounceBack(f);
        }
        // NoDynamics: nothing
    }
};

// HemoFlow's bulk collision (cellPhysics.h) with a Newtonian or Carreau relaxation rate
template <class OmegaFunction>
struct HemoCollision {
    T omegaRegularize;
    OmegaFunction omegaOf;
    void operator()(int model, Array<T, D::numPop> &f, Array<T, D::ExternalField::numScalars> &ext,
                    Array<T, GPUconst<T, DESCRIPTOR>::maxStaticScalars>, T *, plint) const
    {
        if (model == CollisionModel::BGK) {
            T rhoBar;
            Physics::collide(f, &ext[0], omegaRegularize, omegaOf, rhoBar);
        }
        else if (model == CollisionModel::BounceBack) {
            bounceBack(f);
        }
    }
};

struct InitialState {
    plint n;
    void operator()(plint iX, plint iY, plint iZ, T &rho, Array<T, 3> &u) const
    {
        T c = (T)0.5 * (T)(n - 1);
        T zc = (T)0.3 * (T)n;
        T r2 = ((iX - c) * (iX - c) + (iY - c) * (iY - c) + (iZ - zc) * (iZ - zc)) / (T)(n * n);
        rho = (T)1 + (T)0.02 * std::exp(-(T)40 * r2);
        u = Array<T, 3>((T)0.01 * std::sin((T)6.2831853 * (T)iY / (T)n), (T)0,
                        (T)0.01 * std::cos((T)6.2831853 * (T)iX / (T)n));
    }
};

// The box is closed by walls: cells beyond the domain (the envelope of the outer blocks) are
// collided and streamed on the GPU but not on the CPU, so an open box would differ there.
MultiBlockLattice3D<T, DESCRIPTOR> *makeLattice(plint n, plint nb, Dynamics<T, DESCRIPTOR> *background)
{
    MultiBlockManagement3D management = reparallelize(
        defaultMultiBlockPolicy3D().getMultiBlockManagement(n, n, n, 1), n / nb, n / nb, n / nb);
    auto *lattice = new MultiBlockLattice3D<T, DESCRIPTOR>(
        management, defaultMultiBlockPolicy3D().getBlockCommunicator(),
        defaultMultiBlockPolicy3D().getCombinedStatistics(),
        defaultMultiBlockPolicy3D().getMultiCellAccess<T, DESCRIPTOR>(), background);

    Box3D all = lattice->getBoundingBox();
    defineDynamics(*lattice, Box3D(0, 0, 0, n - 1, 0, n - 1), new BounceBack<T, DESCRIPTOR>(1.));
    defineDynamics(*lattice, Box3D(n - 1, n - 1, 0, n - 1, 0, n - 1), new BounceBack<T, DESCRIPTOR>(1.));
    defineDynamics(*lattice, Box3D(0, n - 1, 0, 0, 0, n - 1), new BounceBack<T, DESCRIPTOR>(1.));
    defineDynamics(*lattice, Box3D(0, n - 1, n - 1, n - 1, 0, n - 1), new BounceBack<T, DESCRIPTOR>(1.));
    defineDynamics(*lattice, Box3D(0, n - 1, 0, n - 1, 0, 0), new BounceBack<T, DESCRIPTOR>(1.));
    defineDynamics(*lattice, Box3D(0, n - 1, 0, n - 1, n - 1, n - 1), new BounceBack<T, DESCRIPTOR>(1.));
    plint c = n / 2, R = n / 8;
    for (plint iX = c - R; iX <= c + R; ++iX)
        for (plint iY = c - R; iY <= c + R; ++iY)
            for (plint iZ = c - R; iZ <= c + R; ++iZ)
                if ((iX - c) * (iX - c) + (iY - c) * (iY - c) + (iZ - c) * (iZ - c) <= R * R)
                    defineDynamics(*lattice, iX, iY, iZ, new BounceBack<T, DESCRIPTOR>(1.));
    defineDynamics(*lattice, Box3D(n / 8, n / 4, n - 1 - n / 4, n - 2 - n / 8, 1 + n / 16, n / 6),
                   new NoDynamics<T, DESCRIPTOR>(1.));

    initializeAtEquilibrium(*lattice, all, InitialState {n});
    setExternalVector(*lattice, all, D::ExternalField::forceBeginsAt, Array<T, 3>((T)2e-5, (T)-1e-5, (T)5e-6));
    lattice->initialize();
    return lattice;
}

int main(int argc, char *argv[])
{
    plbInit(&argc, &argv);
    // Needed for the GPU-side packing of MPI messages (see doc/dev/TODO-GPU.md)
    defaultMultiBlockPolicy3D().toggleBlockingCommunication(true);

    plint n = argc > 1 ? std::atoi(argv[1]) : 64;
    plint nb = argc > 2 ? std::atoi(argv[2]) : 1;
    plint steps = argc > 3 ? std::atoi(argv[3]) : 200;
    std::string model = argc > 4 ? argv[4] : "bgk";
    if (n < 16 || nb < 1 || steps < 2 || (model != "bgk" && model != "hemo" && model != "carreau")) {
        pcout << "Usage: gpuSmoke [n >= 16] [blocksPerAxis >= 1] [steps >= 2] [bgk|hemo|carreau]" << std::endl;
        return 1;
    }
    const T omega = 1.6;
    // Carreau-Yasuda in lattice units: omega 1.6 at infinite shear, 10 times the viscosity at zero shear
    CarreauModel carreau;
    carreau.nuInf = (1.0 / omega - 0.5) / 3.0;
    carreau.nu0 = 10.0 * carreau.nuInf;
    carreau.lambda = 2000.0;
    carreau.n = 0.3568;

    Dynamics<T, DESCRIPTOR> *referenceDynamics = nullptr;
    if (model == "bgk")
        referenceDynamics = new GuoExternalForceBGKdynamics<T, DESCRIPTOR>(omega);
    else if (model == "hemo")
        referenceDynamics = new GuoExternalForceCompleteRegularizedBGKdynamics<T, DESCRIPTOR>(omega);
    else
        referenceDynamics = new GuoForcedCarreauDynamics<T, DESCRIPTOR>(omega, carreau);
    MultiBlockLattice3D<T, DESCRIPTOR> *reference = makeLattice(n, nb, referenceDynamics);
    MultiBlockLattice3D<T, DESCRIPTOR> *placeholder = makeLattice(n, nb, new BGKdynamics<T, DESCRIPTOR>(omega));
    plint numBlocks = reference->getMultiBlockManagement().getSparseBlockStructure().getNumBlocks();
    pcout << "gpuSmoke: " << n << "^3 cells, " << numBlocks << " blocks, " << global::mpi().getSize() << " MPI ranks, "
          << steps << " steps, model " << model << ", " << (sizeof(T) == 4 ? "single" : "double") << " precision" << std::endl;

    auto *gpu = new AcceleratedLattice3D<T, DESCRIPTOR>(*placeholder);
    delete placeholder;

    global::mpi().barrier();
    auto t0 = std::chrono::steady_clock::now();
    for (plint i = 0; i < steps; ++i)
        reference->collideAndStream();
    global::mpi().barrier();
    auto t1 = std::chrono::steady_clock::now();
    // First step separately: the data migrates to the GPU
    auto t2 = t1, t3 = t1;
    auto runGpu = [&](auto const &collision) {
        gpu->collideAndStream(collision);
        global::mpi().barrier();
        t2 = std::chrono::steady_clock::now();
        for (plint i = 1; i < steps; ++i)
            gpu->collideAndStream(collision);
        global::mpi().barrier();
        t3 = std::chrono::steady_clock::now();
    };
    if (model == "bgk")
        runGpu(GuoBGKCollision {omega});
    else if (model == "hemo")
        runGpu(HemoCollision<Physics::ConstantOmega> {omega, Physics::ConstantOmega {omega}});
    else
        runGpu(HemoCollision<Physics::CarreauOmega> {omega, Physics::CarreauOmega {carreau}});

    double cpuSec = std::chrono::duration<double>(t1 - t0).count();
    double firstSec = std::chrono::duration<double>(t2 - t1).count();
    double gpuSec = std::chrono::duration<double>(t3 - t2).count();
    double cells = (double)n * n * n;
    pcout << "CPU: " << cpuSec / steps * 1e3 << " ms/step (" << cells * steps / cpuSec * 1e-6 << " MLUPS)" << std::endl;
    pcout << "GPU: first step " << firstSec * 1e3 << " ms, then " << gpuSec / (steps - 1) * 1e3 << " ms/step ("
          << cells * (steps - 1) / gpuSec * 1e-6 << " MLUPS)" << std::endl;

    // Compare the bulk cells of the local blocks
    double maxDiffF = 0.0, maxDiffExt = 0.0, maxF = 0.0, numNonFinite = 0.0;
    MultiBlockManagement3D const &management = reference->getMultiBlockManagement();
    for (auto const &entry : reference->getBlockLattices()) {
        plint id = entry.first;
        BlockLattice3D<T, DESCRIPTOR> &cpuBlock = *entry.second;
        AtomicAcceleratedLattice3D<T, DESCRIPTOR> &gpuBlock = gpu->getComponent(id);
        Box3D bulk = management.getBulk(id);
        Dot3D loc = cpuBlock.getLocation();
        plint ny = gpuBlock.getNy(), nz = gpuBlock.getNz();
        for (plint x = bulk.x0; x <= bulk.x1; ++x)
            for (plint y = bulk.y0; y <= bulk.y1; ++y)
                for (plint z = bulk.z0; z <= bulk.z1; ++z) {
                    plint iX = x - loc.x, iY = y - loc.y, iZ = z - loc.z;
                    Array<T, D::q> f;
                    gpuBlock.pullPop(iZ + nz * (iY + ny * iX), f);
                    Array<T, D::ExternalField::numScalars> ext;
                    gpuBlock.pullExt(iX, iY, iZ, ext);
                    Cell<T, DESCRIPTOR> const &cell = cpuBlock.get(iX, iY, iZ);
                    for (int k = 0; k < D::q; ++k) {
                        if (!std::isfinite((double)f[k]))
                            numNonFinite += 1.0;
                        maxDiffF = std::max(maxDiffF, std::fabs((double)f[k] - (double)cell[k]));
                        maxF = std::max(maxF, std::fabs((double)cell[k]));
                    }
                    for (int k = 0; k < D::ExternalField::numScalars; ++k)
                        maxDiffExt = std::max(maxDiffExt, std::fabs((double)ext[k] - (double)*cell.getExternal(k)));
                }
    }
    double localValues[4] = {maxDiffF, maxDiffExt, maxF, numNonFinite};
    double globalValues[4];
    for (int k = 0; k < 4; ++k)
        global::mpi().reduce(localValues[k], globalValues[k], k == 3 ? MPI_SUM : MPI_MAX);

    // Round-off after the given number of steps; FMA contraction on the GPU makes the results differ in the last bits
    double tolerance = (sizeof(T) == 4 ? 1e-5 : 1e-12) * globalValues[2];
    bool ok = globalValues[0] <= tolerance && globalValues[1] <= tolerance && globalValues[3] == 0.0;
    pcout << "max |f_GPU - f_CPU| = " << globalValues[0] << " (max |f| " << globalValues[2] << "), max |force difference| = "
          << globalValues[1] << ", non-finite GPU values: " << globalValues[3] << std::endl;
    pcout << (ok ? "PASSED" : "FAILED") << " (tolerance " << tolerance << ")" << std::endl;

    delete gpu;
    delete reference;
    int status = ok ? 0 : 1;
    global::mpi().bCast(&status, 1);
    return status;
}
