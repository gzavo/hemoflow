# Running HemoFlow on NVIDIA GPUs (plan)

Status: investigated, not started (November 2026). The solver runs on CPUs only. This
note records what the vendored Palabos offers for GPUs, what HemoFlow needs beyond it, and
the steps of a port. Palabos stays an external dependency: the port should live in HemoFlow
and use the extension points that Palabos provides. Nothing below was compiled or run on a
GPU yet; it comes from reading the vendored Palabos (`palabos/src`) and the upstream GPU
examples.

## Short answer

- Palabos' GPU support alone does not run HemoFlow. It provides the GPU data layout,
  streaming, bounce-back, the block decomposition with MPI exchange and a fixed list of
  collision models. HemoFlow's bulk dynamics (Guo forcing), the porous force, Carreau-Yasuda,
  the pressure and free-outflow conditions and the per-node opening profiles are not in that
  list and have to be written as GPU code in HemoFlow.
- Several GPUs per node and several nodes work through MPI, one rank per GPU, with
  CUDA-aware MPI and the device chosen by the launcher.

## What the vendored Palabos provides

The GPU lattice is `AcceleratedLattice3D` (`multiBlock/acceleratedLattice3D.h`), made of
`AtomicAcceleratedLattice3D` blocks (`atomicBlock/atomicAcceleratedLattice3D.h`).

- It is a separate class from `MultiBlockLattice3D`. It is built by copying a CPU lattice
  (`AcceleratedLattice3D(MultiBlockLattice3D const&)`) and copied back with `writeBack()`.
  The copy keeps the block layout, also a sparse one.
- The kernels are C++ parallel algorithms (`std::for_each(std::execution::par_unseq, ...)`),
  which `nvc++ -stdpar` offloads to the GPU. Populations are stored as structure of arrays;
  collision and streaming are one kernel per block and time step.
- Collision models are a fixed list (`acceleratedLattice/acceleratedDefinitions.h`):
  - `NoDynamics`, `BounceBack`;
  - `BGK`, `BGK_ExternalMoment` (BGK with imposed density and momentum), `TRT`, and the
    moment-space models `RM`, `HM`, `CM`, `CHM`, `K`, `GH`, `RR`;
  - regularized velocity boundary on the six axis-aligned faces (with BGK or TRT) and
    regularized density boundary on the faces (TRT only);
  - half-way bounce-back combined with each bulk model.

  There is no body force, no non-Newtonian or LES dynamics, and no edge or corner boundary.
- When a CPU lattice is copied, the dynamics of every cell and the background dynamics are
  mapped to this list by name (`toCollisionModel`, `acceleratedLattice/acceleratedCollisions.h`).
  An unknown name fails a `PLB_ASSERT`; without assertions it reads past the end of the table.
- The copy also takes over the data processors stored on the CPU lattice (`MultiBlock3D`
  copy constructor), and they still point at the CPU lattice.
- Extension points:
  - `collideAndStream(collFun)`: a user collision function per cell, called with the
    collision-model id of the cell (usable as a tag), the populations, the cell's external
    fields (read and write), the static scalars of the background dynamics (ω) and the
    cell's dynamic scalars.
  - On a GPU block: `for_each(domain, f)`, `transform_reduce(domain, ...)`, `pullExt` /
    `pushExt`, `getDynamicScalar`.
  - Data processors on the GPU lattice derive from `BoxProcessingFunctional3D_A`, `_AA`, `_AS`
    or `_AT` (lattice with a scalar or tensor field) and are applied with the usual
    `applyProcessingFunctional` / `integrateProcessingFunctional`. A force computed on the GPU
    this way: `ShanChenMultiComponentAccelerated3D` (`acceleratedLattice/acceleratedProcessors3D.hh`).
  - GPU functionals for density, velocity, velocity norm, kinetic energy and averages.
- Communication: the envelope exchange is packed on the GPU and sent with MPI (`communicateRaw`
  in `parallelism/parallelBlockCommunicator3D.cpp`).
  - Only static variables (populations and external fields) can be exchanged. Processors that
    declare dynamic variables or the data structure as modified fail an assertion.
  - With `USE_CUDA_MALLOC` the MPI buffers are device memory (`allocateBytes` in
    `core/globalDefs.h`).
  - Sends are non-blocking, receives blocking, and the exchange does not overlap the computation.
- Palabos never selects a GPU (no `cudaSetDevice`).
- Upstream GPU examples (master, `examples/gpuExamples`: cavity3d, multiComponentPorous,
  rayleighTaylor3D, sandstone, tgv3d; not in the vendored copy):
  - they need the NVIDIA HPC SDK and are configured with `cmake .. -DCMAKE_CXX_COMPILER=nvc++`;
  - multiComponentPorous uses `-stdpar -std=c++20 -Msingle -Mfcon -fopenmp -DUSE_CUDA_MALLOC`
    and MPI;
  - upstream documents no device mapping, multi-GPU setup or list of supported features.

## HemoFlow against it

| HemoFlow uses | On the GPU in Palabos | Needed |
|---|---|---|
| Bulk: `GuoExternalForceCompleteRegularizedBGKdynamics` (`globals.h`) | No | Own collision function in HemoFlow |
| Porous force, a processor over the device box (`porous.h`) | No | Inside the collision function; multipliers as extra external fields |
| Carreau-Yasuda (`nonNewtonian.h`, `rheology.h`) | No | Same collision function |
| Smagorinsky LES (`LES` in `globals.h`, off) | No | Only if needed |
| Walls (bounce-back), outside cells | Yes | — |
| Velocity and Murray openings: local regularized condition on rows of face nodes, profile set every step (`ImposeOpeningValues3D`) | Partly: faces only, BGK/TRT only | Own boundary handling with the HemoFlow collision; GPU loop for the per-node values; openings that touch box edges or corners |
| Pressure openings: finite-difference condition, two cells deep | No (only a local regularized density condition) | Port it, or validate a local condition (Palabos' local one diverged on the CPU pipe case) |
| Free outflow (`VirtualOutlet`) | No | Port or keep CPU-only |
| Energy check every 200 steps, debug flow rates | Yes (reductions) | Small rewrite |
| Output (velocity, pressure, stress, viscosity) and checkpoints | Velocity and density only | `writeBack` to the CPU lattice at output steps, keep `io.cpp` |
| Sparse decomposition, blocks of 24³ | Layout kept | Much larger blocks (see below) |
| Geometry read on rank 0 and scattered, dynamics assignment | — | Stays on the CPU, with tag dynamics from the GPU list |

## Plan

### Build

- [ ] CMake switch (e.g. `HEMOFLOW_GPU`): HemoFlow, Palabos and cnpy compiled with `nvc++`
      and the upstream flags; the GCC build stays as it is. All code that allocates GPU data
      must be compiled by `nvc++`.
- [ ] Check that HemoFlow compiles with `nvc++` (EDG front end; the warning flags in
      `CMakeLists.txt` are GCC/Clang only). First with `-stdpar=multicore`, to separate
      compiler problems from GPU problems.
- [ ] Precision: HemoFlow uses double (`globals.h`). Keep double first; single precision
      (half the memory and bandwidth) as a later option, validated on the pipe and NAP180.

### Setup on the CPU, time loop on the GPU

- [ ] Read the geometry, decompose and assign dynamics as now, but with tag dynamics from the
      GPU list: `BGK` for fluid (also porous and Carreau cells), `BounceBack` for walls,
      `NoDynamics` outside, face boundary models (or `BGK` plus a node list) for openings.
- [ ] Do not integrate HemoFlow's processors on that CPU lattice (the copy would take them
      over); integrate GPU processors on the `AcceleratedLattice3D` after the copy.
- [ ] Keep the CPU lattice for output and checkpoints (`writeBack`).
- [ ] Descriptor: external fields for the force (3) and the two porous multipliers (2), so the
      multipliers live with the lattice (alternative: keep the `MultiTensorField3D` and an
      `_AT` processor).

### Collision function (HemoFlow code)

- [ ] Guo-forced complete regularized BGK as in `GuoExternalForceCompleteRegularizedBGKdynamics`,
      using Palabos' header-only template helpers (callable in GPU code).
- [ ] Porous force in the same function, -(m_lin a + m_quad b |u|) u / dx, replacing the
      separate pass over the device box. Keep the CPU order of operations (force from the
      post-streaming velocity, then the collision with the new force) so that GPU and CPU
      results can be compared closely.
- [ ] Carreau-Yasuda: local ω from the non-equilibrium stress with the Newton solve of
      `rheology.h`; the device path must not throw or use virtual calls.

### Openings

- [ ] Velocity and Murray openings: regularized velocity condition with HemoFlow's collision on
      the opening nodes. All nodes of an opening lie on one bounding-box face, so the normal is
      known. Per-node velocities in dynamic scalars or external fields, updated every step by a
      GPU loop over the node list (replaces `ImposeOpeningValues3D`).
- [ ] Opening nodes on an edge or corner of the bounding box (the preprocessor allows openings
      to touch them): decide how to treat the extra missing populations.
- [ ] Pressure openings: port the finite-difference condition as a GPU processor (reads two cells
      inward, envelope 2), or implement and validate a one-cell condition (Zou/He or
      anti-bounce-back, TODO.md section 4) on the pipe first.
- [ ] Free outflow (`VirtualOutlet`): port or document as CPU-only.
- [ ] Energy check and debug flow rates: `transform_reduce` / `AccKineticEnergyFunctional3D`.

### Output and checkpoints

- [ ] At output steps `writeBack` into the CPU lattice and reuse `io.cpp`; outputs are rare
      compared with time steps. GPU functionals for velocity and density exist if the copy
      turns out too slow; stress and viscosity would need custom ones.
- [ ] Checkpoints: `writeBack` and the current save; on restart load into the CPU lattice, then copy.

### Decomposition and several GPUs

- [ ] One MPI rank per GPU. The launcher binds the devices: `CUDA_VISIBLE_DEVICES` per local
      rank (e.g. from `OMPI_COMM_WORLD_LOCAL_RANK` in a wrapper script) or Slurm
      `--gpus-per-task=1`; otherwise all ranks of a node use GPU 0. Document it in the README;
      optionally check at start-up that the ranks of a node see different devices.
- [ ] CUDA-aware MPI (e.g. the HPC-X / OpenMPI of the NVIDIA HPC SDK), required with
      `USE_CUDA_MALLOC`. The conda OpenMPI used on the workstation would have to be checked.
- [ ] Block layout: the CPU sparse default (blockSize 24) gives hundreds of small blocks per
      GPU, each one kernel launch per step plus copies for every block face. GPUs need few
      large blocks: dense with one block per GPU, or sparse with large blocks, balanced by
      fluid count. Measure on NAP180.
- [ ] Several nodes: plain MPI (GPUDirect RDMA if the cluster has it). The exchange does not
      overlap the computation, so each GPU should hold several million cells; small cases
      scale over few GPUs only.

### Verification

- [ ] Pipe (`tests/verification/0-pipe`): GPU against CPU, pressure drop and profile; Carreau
      pipe (`carreau.xml`).
- [ ] NAP180 with the flow diverter: GPU against CPU (expect agreement to rounding where the
      order of operations is kept; report the differences).
- [ ] 1, 2 and 4 GPUs give the same results, as process counts do on the CPU.
- [ ] Performance: cell updates per second against block size and layout, compared with the
      CPU (NAP180 20 M sparse: 300 steps in 44.5 s on 20 cores, about 25 M updates/s of
      allocated cells).

## Hardware and expected speed

- LBM is limited by memory bandwidth. A cell update moves about 380 bytes in double (two
  population arrays, force, multipliers, model id, bounce-back links), so an A100 (1.5–2 TB/s)
  is bounded at about 4–5 G updates/s, probably half of that in practice, against about
  25 M/s measured on the workstation. To be measured.
- Double precision: data-centre GPUs (A100, H100) run it at full speed; GeForce and RTX
  workstation GPUs at 1/32 to 1/64 of single precision, so the collision becomes
  compute-bound there.
- Memory: NAP180 20 M voxels dense is about 7 GB in double; sparse less.

## Prerequisites on the workstation (November 2026)

- NVIDIA GPU (PCI ID 10de:25b0); driver packages 570 and 580 are installed, but the kernel
  module is not loaded (`nvidia-smi` fails, PRIME on-demand).
- CUDA 12.0 `nvcc` (Ubuntu package); the NVIDIA HPC SDK (`nvc++`) is not installed.

## Effort (rough)

- About 1 week: build and a first GPU pipe run (BGK, walls, velocity inlet, simple outlet).
- 2–3 weeks: full physics and boundary conditions, output, checkpoints.
- 1–2 weeks: verification and multi-GPU tuning.

## Open questions

1. Which GPUs (local or cluster, model, double-precision speed), and is CUDA-aware MPI available there?
2. Double or single precision?
3. Are free outflow and LES needed on the GPU?
4. One binary for CPU and GPU, or two builds?
