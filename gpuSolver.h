#ifndef GPU_SOLVER_H
#define GPU_SOLVER_H

// Time loop on an NVIDIA GPU (doc/dev/TODO-GPU.md). Compiled only in the GPU build
// (HEMOFLOW_GPU, nvc++ -stdpar).
//
// The GPU lattice is a Palabos AcceleratedLattice3D, copied from a CPU lattice with
// placeholder dynamics: every cell carries a dynamics whose name Palabos maps to one of its
// GPU collision models. The model id then only tags the cell; the collision functor of this
// module does HemoFlow's physics (cellPhysics.h) for each tag:
//
//   NoDynamics                                   outside voxels: nothing
//   BounceBack                                   walls: full-way bounce-back
//   BGK (background)                             fluid: forced complete regularized BGK,
//                                                Newtonian or Carreau-Yasuda
//   Boundary_RegularizedVelocity_<d>_<o> >> BGK  velocity and Murray openings: the fluid
//                                                collision; the imposed velocity is in 3 dynamic
//                                                scalars
//   Boundary_RegularizedDensity_<d>_<o> >> TRT   pressure openings: the fluid collision; the
//                                                imposed rhoBar is in 1 dynamic scalar
//
// A time step follows the CPU:
//   1. collision and streaming;
//   2. velocity nodes rebuilt from the imposed velocity (on the CPU WrappedLocalBoundaryFunctional3D
//      after streaming);
//   3. pressure nodes rebuilt by the finite-difference condition (on the CPU PlaneFdBoundaryFunctional3D
//      after streaming), in two passes: the velocities of all pressure nodes into a buffer, then the
//      rebuild from their neighbours' velocities. The CPU updates the nodes one after another in
//      place, so the two agree at round-off, not bit for bit;
//   4. porous force of a flow diverter (on the CPU PorousForceFunctional, after the opening
//      processors): from the velocity after streaming, written into the force field that the
//      next collision uses. Only porous fluid cells, from a list made at start-up;
//   5. envelope exchange (populations and external scalars).
//
// The CPU lattice with the real dynamics stays: copyTo() brings the GPU state into it for
// output, checkpoints and diagnostics, so that io.cpp and the checkpoint format do not change.

#include "globals.h"
#include "helper.h"
#include "opening.h"

#include <cstddef>
#include <string>
#include <vector>

// The GPU of this process. Returns false if there is none.
bool queryGpu(std::string& name, std::size_t& freeBytes, std::size_t& totalBytes);

// Release the CUDA context of this process (MPI ranks that wait while another rank uses the
// GPU). Afterwards this process must not touch memory it allocated before.
void releaseGpu();

// GPU memory per cell of the GPU lattice: two population arrays, the external scalars,
// the model id, the bounce-back links and the dynamic-scalar index.
std::size_t gpuBytesPerCell();

// GPU memory per porous cell of a flow diverter (cell index and two multipliers)
std::size_t gpuBytesPerPorousCell();

// Porous medium of a flow diverter (porous.h)
struct GpuPorousMedium {
    MultiTensorField3D<T, 2> const* multipliers = nullptr;   // (linear, quadratic) per voxel; nullptr: none
    Box3D box;                                             // Bounding box of the porous voxels
    T linCoeff = 0, quadCoeff = 0;                         // a and b in lattice units
};

class GpuSolver {
public:
    // placeholder: CPU lattice with placeholder dynamics and the state to start from; it is
    // copied and can be deleted afterwards. openings: their boundary values are written to the
    // GPU by imposeOpenings(), which the constructor calls once. porous: read once (the
    // multipliers of the porous fluid cells are copied).
    GpuSolver(MultiBlockLattice3D<T, DESCRIPTOR> const& placeholder, std::vector<OpeningHandler*> const& openings,
              SimPar const& sim, GpuPorousMedium const& porous);
    ~GpuSolver();

    // One time step (collision, streaming, opening nodes, porous force, envelope exchange)
    void collideAndStream();
    // Write the current boundary values of the openings (OpeningHandler::getVelocityScale,
    // getDensity) to the GPU
    void imposeOpenings(std::vector<OpeningHandler*> const& openings);
    // Populations and external scalars into a CPU lattice with the block structure of the placeholder
    void copyTo(MultiBlockLattice3D<T, DESCRIPTOR>& lattice) const;
    // Sum of u^2 over the domain divided by its number of cells, as Palabos' computeAverageEnergy
    T averageEnergy() const;

    // Velocity of a neighbour of a pressure node, by the type of the neighbour
    struct VelocitySource {
        enum Kind { ZERO, FLUID, VELOCITY_NODE, PRESSURE_NODE };
        int kind;
        plint cell;   // FLUID: linear index of the cell in its block
        plint aux;    // VELOCITY_NODE: dynamic-scalar index; PRESSURE_NODE: index in the velocity buffer
    };
    // Fluid cell in the porous medium
    struct PorousCell {
        plint cell;       // Linear index of the cell in its block
        T mLin, mQuad;    // Resistance multipliers
    };

private:
    struct VelocityNode {
        plint cell;             // Linear index of the cell in its block
        plint index;            // Of the node's dynamic scalars (imposed velocity)
        double ux, uy, uz;      // Unit-peak profile (OpeningHandler)
    };
    struct PressureNode {
        plint cell;
        plint index;            // Of the node's dynamic scalar (imposed rhoBar)
        plint buffer;           // Of the node's velocity in the velocity buffer
        VelocitySource inward;  // The neighbour inside the domain
        VelocitySource plus[3], minus[3];   // Neighbours at +-1 along the axes of the face
        int rebuild;            // 0 for copies in the outer envelope layer (neighbours outside the
                                // block): only their velocity is computed, the exchange sets them
    };
    // Nodes of one opening in one block (envelope copies included: dynamic scalars are not
    // exchanged between blocks)
    template <class Node>
    struct NodeGroup {
        int opening;                  // Index in the openings vector
        int direction, orientation;   // Face of the opening (Palabos convention: orientation = outward normal)
        AtomicAcceleratedLattice3D<T, DESCRIPTOR>* block;
        T* dynamicScalars;
        std::vector<Node> nodes;
    };

    struct PorousGroup {
        AtomicAcceleratedLattice3D<T, DESCRIPTOR>* block;
        std::vector<PorousCell> cells;   // In the bulk of the block; the exchange fills the envelope
    };

    void completeVelocityNodes();
    void completePressureNodes();
    void applyPorousForce();

    AcceleratedLattice3D<T, DESCRIPTOR>* lattice;
    std::vector<NodeGroup<VelocityNode> > velocityNodes;
    std::vector<NodeGroup<PressureNode> > pressureNodes;
    std::vector<Array<T, 3> > pressureVelocities;   // Velocity buffer of the pressure nodes
    std::vector<PorousGroup> porousCells;
    T linCoeff, quadCoeff;      // Porous coefficients a, b (lattice units)
    T omega;                    // Newtonian omega, or the one at infinite shear (Carreau)
    bool nonNewtonian;
    CarreauModel carreau;
    Box3D domain;
};

#endif
