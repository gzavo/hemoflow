# Running HemoFlow on NVIDIA GPUs (plan)

Status: Phases 0-5 done (21-22 November 2026). The GPU build runs velocity, Murray and pressure
openings, flow diverters, Newtonian and Carreau fluid on the GPU, verified against the CPU on the
verification cases. Next: Phase 6 (output speed, several GPUs). First investigated on 10 November 2026 (that version is in git, commit bca490e);
revised on 21 November 2026 after a second reading of the vendored Palabos (`palabos/src`), the
upstream GPU examples and the workstation. Target precision: single (decided on 21 November
2026); double only if single turns out insufficient. Palabos stays an external dependency: the
port lives in HemoFlow and uses the extension points that Palabos provides.

## Short answer

- Palabos' GPU lattice (`AcceleratedLattice3D`) provides the GPU data layout, streaming,
  bounce-back, the block decomposition and the envelope exchange. HemoFlow's bulk dynamics
  (complete regularized BGK with Guo forcing), the porous force, Carreau-Yasuda, the opening
  boundary conditions and the output fields are not available in a form HemoFlow can use, so
  they have to be written as GPU code in HemoFlow: a collision functor plus a few GPU data
  processors.
- Toolchain: the NVIDIA HPC SDK (`nvc++ -stdpar`), installed on the workstation (26.11). The
  Ubuntu CUDA 12.0 `nvcc` cannot compile Palabos' GPU code, which is written with C++ parallel
  algorithms, not CUDA kernels.
- The workstation GPU (RTX A1000, 8 GB) is good for development. It runs double precision at
  1/64 of single precision. On a simple forced BGK test it reaches about 0.75 G cell updates/s
  in single and 0.15 G/s in double, against about 0.025 G/s for the 20-core CPU run of NAP180.
  HemoFlow's collision is heavier, and HemoFlow's layout will be less favourable (next point).
- Every block of the decomposition costs synchronous GPU kernels per step. The CPU default
  (sparse, blocks of 24) would make the GPU about as slow as the CPU, so the GPU build needs one
  or a few large blocks per GPU (measured, see "Phase 0 results").
- Several GPUs and several nodes work through MPI, one rank per GPU, with CUDA-aware MPI and
  the device chosen by the launcher.

## Changes from the first version (10 November)

The overall approach of the first version holds (setup on the CPU, time loop on the GPU, own
collision function, one MPI rank per GPU). These points were wrong, missing or out of date:

1. **Communicator.** The GPU-side packing of the MPI messages (`communicateRaw`) exists only in
   `BlockingCommunicator3D`. HemoFlow creates its lattice with the default non-blocking
   `ParallelBlockCommunicator3D`, which packs into host buffers. The GPU build must call
   `defaultMultiBlockPolicy3D().toggleBlockingCommunication(true)` before the lattice is
   created, as both upstream GPU examples do. With a single rank both paths do only local
   copies, so the difference shows only with several ranks.
2. **Compiler flags.** The upstream flags `-Msingle -Mfcon` are for the examples' single
   precision: `-Mfcon` turns every floating-point literal into a float. They must not be used
   in a double build.
3. **Output.** `writeBack()` also copies the dynamics objects of the GPU lattice into the CPU
   lattice. Those are the placeholder dynamics needed for the copy (see "Design"), with the
   opening values they had when the copy was made. `io.cpp` on such a lattice gives velocity
   and stress without the forcing correction in the porous voxels and stale values at the
   openings. Palabos' GPU velocity functional has no forcing correction either, and its density
   functional runs on the host. The plan now keeps a CPU lattice with the real dynamics and
   copies back only populations and external fields (see "Output and checkpoints").
4. **Palabos' GPU velocity boundary is not HemoFlow's.** `collideRegularizedBoundary`
   (`acceleratedCollisionsBoundary.h`) regularizes with the second-order equilibrium and
   collides with plain BGK/TRT. HemoFlow's CPU condition (`RegularizedVelocityBoundaryDynamics`
   around the bulk dynamics) uses the complete equilibrium of the bulk dynamics for both
   Π^neq and the regularization. HemoFlow can reuse the structure and
   `computeBoundaryRhoBar`, but not the rest.
5. **Per-cell data in the collision function.** The collision function receives the model id,
   the populations, the external fields, the static scalars and the cell's dynamic scalars, but
   not the cell index. Data that differs per cell must be an external field (read and written
   for every cell, every step) or dynamic scalars (only for cells whose placeholder model has
   them). Dynamic scalars are not exchanged between blocks, so per-step updates must also write
   the envelope copies of a node.
6. **Porous multipliers.** Making them external fields costs about 32 bytes per cell and step
   for every cell, which is roughly 9 % more memory traffic. The plan now keeps the 3-scalar
   `ForcedD3Q19Descriptor` and computes the porous force in a GPU processor over the device
   box with the existing multiplier field, which matches the CPU code step for step. Folding the
   force into the collision stays as an optimization. (Phase 4: a kernel over a list of the
   porous cells instead of a processor over the box.)
7. **Openings on box edges and corners cannot occur.** The preprocessor makes an opening voxel
   only where fluid lies directly behind the wall layer, and every fluid voxel has a wall voxel
   on each side inside the cropped domain (`detectOpenings.py`, `createFluidSolid.py`). So no
   opening node lies on an edge of the bounding box; only geometry files from the old
   preprocessor could have such nodes. The GPU build should check this and refuse such files;
   then no edge or corner treatment is needed.
8. **Block layout.** Each block costs one kernel for collision and streaming, and each block
   overlap costs one copy kernel per step. Every parallel algorithm call returns only when its
   kernel has finished, and the kernels sweep every allocated cell, including outside and wall
   cells. Measured in Phase 0: about 30 µs per block overlap and step. With the CPU default
   (blocks of 24, about 270 blocks for NAP180 20 M) this is roughly 100 ms per step, about as
   slow as the CPU. The GPU build therefore needs one or a few large blocks per GPU. Dense NAP180
   fits the 8 GB A1000 in single (3.8 GB) but not in double (7.5 GB).
9. **Pressure openings.** The finite-difference condition (`PlaneFdBoundaryFunctional3D`) reads
   the neighbours' velocities after streaming and before the exchange. On the GPU the outer
   envelope layer is incomplete after streaming, as on the CPU, so a faithful port also needs
   envelope width 2. On the GPU this costs little with large blocks.
10. **One physics source for CPU and GPU (new).** Write HemoFlow's own per-cell physics (bulk
    collision, Carreau ω, porous force) once as array-based inline functions, used by both the
    CPU dynamics and the GPU functor. The models still to come (coils, semi-implicit forcing,
    device-level flow diverter, Carreau-scaled Darcy term) then reach both builds without being
    written twice. The device ports of Palabos' CPU boundary conditions get CPU unit tests
    against those classes.
11. **Workstation.** The NVIDIA driver is now loaded (595.99, CUDA 13.2), and the HPC SDK is
    installed (26.11, with a CUDA-aware HPC-X Open MPI). The system OpenMPI 4.1.6 is not
    CUDA-aware. The proprietary kernel module gives no HMM, so the GPU can reach only heap
    memory allocated by `nvc++`-compiled code (CUDA managed memory).
12. **Testing.** Several ranks can share the one GPU, so decomposition independence can be
    tested on the workstation. A build without fused multiply-add may reproduce the CPU results
    bit for bit for Newtonian cases.

## Phase 0 results (21 November 2026)

**Toolchain.** NVIDIA HPC SDK 26.11 in `/opt/nvidia/hpc_sdk`:
- `nvc++` 26.11 with CUDA 13.4;
- HPC-X 2.50 with Open MPI 5.0.10, which has the CUDA accelerator component and UCX with CUDA.

The driver (595.99, CUDA 13.2) runs the code built with CUDA 13.4. `cudaMalloc` (Palabos'
`USE_CUDA_MALLOC`) compiles with `-stdpar=gpu` alone.

**Build.** New CMake options:
- `HEMOFLOW_GPU`, `HEMOFLOW_STDPAR` (`gpu` | `multicore`), `HEMOFLOW_GPU_ARCH` (default
  `ccnative`);
- `HEMOFLOW_PRECISION` (`single` | `double`; single by default with `HEMOFLOW_GPU`).
  `HEMOFLOW_SINGLE` sets `T` to float in `globals.h`.

Palabos compiles with `nvc++` unchanged. HemoFlow needed two changes:
- `using namespace std;` moved after the Palabos includes (`hemoFlow.cpp`, `helper.h`): before
  them it makes `__declval` in the SDK's `std::execution` implementation (stdexec) ambiguous;
- one `(T)` cast for single precision (`initializeAtEquilibrium`).

With these changes the GCC double build writes byte-identical checkpoints to the committed code
(pipe; pipe with a flow diverter of a = 100, b = 200; Carreau pipe; 1000 steps on 8 ranks).

**`nvc++` against GCC, CPU time loop, same pipe cases.**
- Double: velocity within 3·10⁻¹¹ relative; pressure identical in the (float) output.
- Single: within 2·10⁻⁶ relative, i.e. float rounding.
- Single against double (GCC) after 1000 steps: about 1·10⁻⁵ relative in velocity and pressure.
- Full pipe validation (`pipe.xml`, 1.6 s, 80 k steps, GCC, 16 ranks), single against double:
  - pressure drop between x = L/4 and 3L/4 88.457 against 88.479 Pa (0.025 %);
  - identical flow rates;
  - velocities within 4.4·10⁻⁵ m/s (5·10⁻⁵ relative), pressures within 0.045 Pa.

  So single precision is enough for this case. The evaluation (scratch script, mean over the
  fluid voxels of the planes) gives 88.5 Pa, not the 78.7 Pa of `pressure_drop.py`: that script
  averages cell data, in which the cells at the wall mix in the 0 Pa of the wall nodes.

**Speed of the CPU paths.**
- **GCC single.** In short runs (1000 steps from rest) the GCC single build was 2.7 times slower
  than double, because of denormal floats while the flow develops. Single-precision builds now
  flush denormals to zero on the CPU (GCC `-mdaz-ftz`, `nvc++` `-Mdaz -Mflushz`). With that,
  1000 steps take 5.0 s in single against 6.6 s in double, with results changed only at float
  rounding (3·10⁻⁶). Over the full pipe run, single was faster than double even without the
  flush (263 against 358 s).
- **`nvc++` CPU loop.** Its speed is fine: the `-stdpar=multicore` build is faster than GCC
  (5.7 against 6.6 s). The `-stdpar=gpu` build adds 2-3 s of CUDA start-up per run and about
  10 % on one rank (managed memory). With several ranks it did not scale at first: 8 ranks took
  26 s, as long as 1 rank. The cause is HPC-X's CUDA-aware UCX, which handles the MPI messages
  in managed host memory with its CUDA transports. With `UCX_TLS=^cuda` 8 ranks take 7.3 s,
  most of it CUDA start-up, against 3.5 s for 16 ranks of the GCC single build. Each rank also
  creates a CUDA context of about 120 MB on the GPU.

**GPU smoke test** (`tests/gpu/gpuSmoke.cpp`, target `gpuSmoke` of the GPU build). A closed box
with walls, a sphere, a NoDynamics region, a density bump and a body force. The CPU reference
uses `GuoExternalForceBGKdynamics`; the GPU runs a copy of a placeholder lattice (BGK background)
with a HemoFlow collision functor (Guo-forced BGK, force from the external field).
- **Accuracy.** GPU against CPU after 200 steps: 6.5·10⁻¹⁸ absolute in double (round-off) and
  3·10⁻⁹ in single (7·10⁻⁷ relative). Results are the same with 1, 8 and 64 blocks, and with 2
  MPI ranks sharing the GPU (HPC-X, buffers in managed memory). So placeholder dynamics, the
  functor, external fields and the envelope exchange on and between ranks all work.
- **Open box.** With an open box (a NoDynamics region replacing part of the outer wall) GPU and
  CPU differ. The GPU collides and streams the cells beyond the domain (the envelope of the
  outer blocks); the CPU does not. For HemoFlow this is harmless: the domain faces hold walls,
  which reflect, and openings, which rebuild the populations that come from outside.
- **Speed, one block.**
  - Single: 640-800 M cell updates/s at 64³-128³, i.e. 74-80 % of the 192 GB/s.
  - Double: 140-170 M/s, limited by double-precision arithmetic.
  - For comparison: 10-27 M/s for one CPU core in the same test, and about 25 M/s for NAP180
    on 20 cores.
- **Speed, several blocks.** 128³ in single costs 2.6 ms per step with 1 block, 4.6 ms with 8
  and 55.6 ms with 64; 64³ costs 0.35 ms with 1 block and 56 ms with 64. Nsight Systems on 64
  blocks shows about 1,900 small kernels per step, two per block overlap. Each takes about
  15 µs (launch, about 5 µs on the GPU, synchronization), so about 30 µs per overlap. NAP180
  20 M with the CPU default (about 270 blocks of 24, each with 10-20 neighbours) would spend
  roughly 100 ms per step on overlaps (estimate), close to the 148 ms of the 20-core CPU run.

**Consequences.**
- The GPU build needs one block per GPU (dense) or a few large blocks. NAP180 20 M dense in
  single needs about 190 bytes per cell, i.e. 3.8 GB, which fits the A1000. At the measured
  bandwidth a dense step would take roughly 25-35 ms, 4-6 times faster than the CPU, although
  80 % of the cells are outside the vessel.
- Sparse layouts with a handful of large blocks may do better and are measured in Phase 6.
- A fluid-only layout would cut this further, but Palabos' GPU lattice has no such layout;
  outside the scope of this plan.

## What the vendored Palabos provides

The GPU lattice is `AcceleratedLattice3D` (`multiBlock/acceleratedLattice3D.h`), made of
`AtomicAcceleratedLattice3D` blocks (`atomicBlock/atomicAcceleratedLattice3D.h`).

- **Copy in.** It is built by copying a CPU lattice (`AcceleratedLattice3D(MultiBlockLattice3D
  const&)`). The copy keeps the block layout, including a sparse one, and the envelope width.
  The processors themselves are not copied, so they would not run on the GPU. The copy does
  take over the CPU lattice's list of blocks that processors modify (`MultiBlock3D` copy
  constructor), and that list still points at the CPU lattice and its fields. Those would then
  be exchanged every step for nothing, or through dangling pointers once deleted. So the lattice
  that is copied must have no processors.
- **Kernels.** The kernels are C++ parallel algorithms (`std::for_each(std::execution::par_unseq,
  ...)`), which `nvc++ -stdpar` offloads to the GPU; every call is synchronous. Populations are
  stored as a structure of arrays (two copies), and collision plus push-streaming is one kernel
  per block over all its cells, envelope included. Memory comes from `new`, which `nvc++` turns
  into CUDA managed memory. If the host touches it during the run, the pages migrate.
- **Collision models.** There is a fixed list (`acceleratedLattice/acceleratedDefinitions.h`):
  `NoDynamics`, `BounceBack`, `BGK`, `BGK_ExternalMoment`, `TRT`, the moment-space models (`RM`,
  `HM`, `CM`, `CHM`, `K`, `GH`, `RR`), regularized velocity faces (BGK or TRT, 3 dynamic
  scalars: u), regularized density faces (TRT, 1 dynamic scalar: ρ), and half-way bounce-back
  with each bulk model. There is no body force, no non-Newtonian or LES model and no edge or
  corner boundary.
- **Model mapping.** When a CPU lattice is copied, every cell's dynamics and the background
  dynamics are mapped to this list by their name chain, e.g. `Boundary_RegularizedVelocity_2_-1
  >> BGK` (`toCollisionModel`, `acceleratedCollisions.h`). An unknown name fails a `PLB_ASSERT`;
  without assertions it reads past the end of the table. The static scalars passed to the
  collision come from the background dynamics (`staticScalars[0]` = ω). Dynamic scalars are
  filled once, at the copy, from the CPU dynamics (`getDynamicScalars`).
- **Extension points:**
  - `collideAndStream(collFun)`: a user functor called per cell as `collFun(model, f, ext,
    staticScalars, dynamicScalars, index)`. The model id can serve as a tag; `ext` is read and
    written. The functor is captured by value, so it must be trivially copyable.
  - On a GPU block: `for_each(domain, f)` (still one sweep over the whole block, filtered by
    `domain`), `transform_reduce(domain, ...)`, `pullPop`, `pullExt` / `pushExt` and
    `getDynamicScalar`.
  - Data processors on the GPU lattice derive from `BoxProcessingFunctional3D_A`, `_AA`, `_AS` or
    `_AT` (lattice with a scalar or tensor field), and reductions from
    `ReductiveBoxProcessingFunctional3D_A` / `_AS` / `_AT`. They are used with the usual
    `applyProcessingFunctional` / `integrateProcessingFunctional`. Example:
    `ShanChenMultiComponentAccelerated3D` (`acceleratedProcessors3D.hh`).
  - Array-based moment and collision templates (`momentTemplatesImpl`, `dynamicsTemplatesImpl`,
    `externalForceTemplatesImpl`, `offEquilibriumTemplates`) are header-only and usable on the
    device. The `Cell`-based versions that HemoFlow uses today only forward to them
    (`latticeBoltzmann/dynamicsTemplates.h`), so an array-based port can do exactly the same
    arithmetic.
  - GPU functionals: velocity, velocity norm and kinetic energy compute j/ρ from the
    populations only, without the forcing correction and whatever the cell type. Density uses
    `reconstructCell` on the host.
- **Communication.** In `BlockingCommunicator3D::communicateRaw`
  (`parallelism/parallelBlockCommunicator3D.cpp`) the envelope data is packed on the GPU and
  sent with MPI: sends non-blocking, receives blocking, no overlap with computation. Only
  static variables (populations and external fields) can go this way; modifying dynamic
  variables or the data structure fails an assertion. With `USE_CUDA_MALLOC` the message
  buffers are device memory (`allocateBytes`, `core/globalDefs.h`), which needs CUDA-aware MPI.
  Local copies between blocks of one rank are GPU kernels (`attribute_static`).
- **Checkpoints.** `saveBinaryBlock` / `loadBinaryBlock` serialize static variables through the
  block's data transfer (`multiBlockSerializer3D.cpp`), so they also accept the GPU lattice.
  Whether its files are byte-identical to the CPU lattice's is to be checked.
- Palabos never selects a GPU (no `cudaSetDevice`).
- **Upstream GPU examples** (master, `examples/gpuExamples`: cavity3d, multiComponentPorous,
  rayleighTaylor3D, sandstone, tgv3d; not in the vendored copy):
  - built with `cmake .. -DCMAKE_CXX_COMPILER=nvc++`, flags `-stdpar -std=c++20 -Msingle -Mfcon
    -fopenmp -DUSE_CUDA_MALLOC`;
  - `typedef float T`;
  - `toggleBlockingCommunication(true)`;
  - a CPU lattice set up with ordinary dynamics and copied once;
  - `collideAndStream(CollisionKernel<T, DESCRIPTOR, model, ...>())`;
  - a body force set as an external field before the copy;
  - output from the GPU functionals;
  - no device selection, multi-GPU notes or minimum SDK version.

## HemoFlow against it

| HemoFlow uses | On the GPU in Palabos | Needed in HemoFlow |
|---|---|---|
| Bulk: `GuoExternalForceCompleteRegularizedBGKdynamics` (`globals.h`) | No | Collision functor, array port of the same arithmetic |
| Carreau-Yasuda (`nonNewtonian.h`, `rheology.h`) | No | Same functor; local ω from Π^neq (device code: no throw, no virtual calls) |
| Porous force, a processor over the device box (`porous.h`) | No | GPU `_AT` processor with the multiplier field |
| Walls (full-way bounce-back), outside cells | Yes | — |
| Velocity and Murray openings: regularized condition on the opening nodes, profile set every step (`ImposeOpeningValues3D`) | Faces only, second-order regularization, BGK/TRT bulk | Own boundary collision with the complete regularization; GPU update of the per-node velocities |
| Pressure openings: finite-difference condition (`createInterpBoundaryCondition3D`) | No (only a local regularized density condition, the kind that diverged on the CPU pipe) | GPU port of `PlaneFdBoundaryFunctional3D`, envelope 2 |
| Free outflow (`VirtualOutlet`) | No | Refuse in the GPU build (port only if needed) |
| Smagorinsky LES (`LES` in `globals.h`, off) | No | Refuse in the GPU build |
| Energy check every 200 steps, debug flow rates | Reductions exist | Own reductions with HemoFlow's velocity definition |
| Output (velocity, pressure, stress, viscosity), checkpoints | Velocity and density without forcing correction | Copy populations and external fields back to a CPU lattice with the real dynamics; `io.cpp` unchanged |
| Sparse decomposition, blocks of 24³ | Layout kept | Larger blocks; block size tuned for the GPU |
| Geometry read on rank 0 and scattered, dynamics assignment | — | Stays on the CPU, with placeholder dynamics for the copy |

## Design

**Two builds.** The GCC build stays as it is, with bit-identical results and no new
dependency. A GPU build is made with `nvc++` (CMake option `HEMOFLOW_GPU`), and the time loop
is chosen at compile time. A single `nvc++` binary for both is possible, but `-stdpar=gpu`
places every heap allocation in managed memory, and the speed of the CPU path in such a binary
is untested.

**Three lattices, same block layout.**

1. The CPU lattice as today, with the real dynamics and boundary conditions but without the
   porous processor. It is used for output, checkpoints and debug measurements.
2. A CPU lattice with placeholder dynamics. It exists only to be copied and is deleted after
   the copy. It must have no data processors.
3. The GPU lattice, a copy of (2).

Placeholder dynamics (the model id is only a tag; HemoFlow's functor decides what happens):

| Cells | Placeholder on CPU lattice (2) | GPU model id | Per-cell data on the GPU |
|---|---|---|---|
| Outside (label 0) | `NoDynamics` | `NoDynamics` | — |
| Wall | `BounceBack` | `BounceBack` | — |
| Fluid, also porous and Carreau | background `BGKdynamics(ω)` | `BGK` | ext: force (3) |
| Velocity / Murray opening, face (d, o) | `RegularizedVelocityBoundaryDynamics<d,o>(new BGKdynamics(ω))` | `Boundary_RegularizedVelocity_d_o__BGK` | 3 dynamic scalars: u |
| Pressure opening, face (d, o) | `RegularizedDensityBoundaryDynamics<d,o>(new TRTdynamics(ω))` | `Boundary_RegularizedDensity_d_o__TRT` | 1 dynamic scalar: ρ |

**Time step on the GPU, in the CPU order:**

1. Write the opening values (velocity per node, density per pressure opening) into the dynamic
   scalars. This is one small kernel per block that holds opening nodes, its envelope copies
   included. The Murray split and the waveforms stay on the host, as now.
2. `collideAndStream(HemoCollision{...})`.
3. Level-0 GPU processors, after streaming: the pressure condition on the pressure nodes, then
   the porous force over the device box (bulk only).
4. Envelope exchange (Palabos, static variables).

The CPU does the same: `imposeOpenings`, collide, stream, then the finite-difference pressure
condition and the porous force (processors), then the exchange. This keeps GPU and CPU results
comparable to rounding.

**Device code rules.**

- The GPU sees only heap memory allocated by `nvc++`-compiled code. Lambdas capture by value,
  or by reference only to heap objects.
- No virtual calls, exceptions, `std::map` or `std::unordered_set`, and no `Dynamics` objects in
  kernels.
- No host loops over lattice data inside the time loop, because they migrate pages.
- `OpeningHandler`'s node sets and profile maps become flat per-block arrays (cell index,
  dynamic-scalar index, unit profile).
- Single precision: no double literals or double variables in device code (`0.5 * x` with a
  float `x` is computed in double, at 1/64 speed on the A1000). Write `(T)0.5` as Palabos does.
  GCC's `-Wdouble-promotion` on the single build finds such places (Phase 0):
  - Palabos' collision and Guo-forcing templates are clean;
  - `porous.h` and `nonNewtonian.h` have a few comparisons with `0.0`;
  - `rheology.h` is double throughout.

  Check the GPU code for double instructions with Nsight Compute or `cuobjdump`.

## Plan

### Phase 0: toolchain (done 21 November 2026, see "Phase 0 results")

- [x] NVIDIA HPC SDK 26.11 installed (user).
- [x] Palabos and HemoFlow build with `nvc++` (`-stdpar=gpu` and `-stdpar=multicore`); the CPU
      time loop of the `nvc++` build reproduces GCC on the pipe (double to round-off, single to
      float rounding).
- [x] Smoke test of the GPU lattice on the A1000 (`tests/gpu/gpuSmoke.cpp`): passes with 1-64
      blocks and 2 ranks.
- [x] CMake options `HEMOFLOW_GPU`, `HEMOFLOW_STDPAR`, `HEMOFLOW_GPU_ARCH`, `HEMOFLOW_PRECISION`;
      HDF5 is not found with HPC-X and is skipped. `-Mfcon` / `-Msingle` are not used. They
      would make every double literal a float, which suits a pure single build but also
      changes the host code (unit conversions, geometry).
- [ ] Later, when needed:
  - `-DUSE_CUDA_MALLOC` (only with CUDA-aware MPI, i.e. HPC-X);
  - `-DUSE_NVIDIA_PROFILER` for the NVTX ranges of Palabos (Nsight Systems already shows
    those of the stdpar runtime);
  - a validation variant without fused multiply-add (`-gpu=nofma`, `-Kieee`; check the flag
    names).
- [x] Speed of the CPU paths explained (see "Phase 0 results"): denormals in single
      (flush-to-zero added), and CUDA start-up plus managed memory in the `-stdpar=gpu` build.
- [x] Single precision validated on the CPU against double on the full pipe (0.025 % in the
      pressure drop).

### Phase 1: preparation on the CPU (done 21 November 2026, results bit-identical)

- [x] `main()` split into steps (`hemoFlow.cpp`):
  - `printBanner`, `parseCommandLine`;
  - `readSimulationSection` (`<simulation>`, output fields, `<fluid>`), `createOutputFolder`;
  - `setUpFromGeometry` (domain, output frame, flow diverter, LB parameters);
  - `readOpenings`, `createLattice`, `createFieldsAndOpenings` (labels, multipliers, boundary
    conditions), `setCyclePeriod`, `printSimulationParameters`;
  - `runSimulation` (initial state or checkpoint, time loop, output, checkpoints).

  The blocks moved unchanged and the globals stay (replacing them with a configuration struct
  is still TODO.md section 5). The "placeholder dynamics" switch of `createLattice` /
  `createFieldsAndOpenings` comes with Phase 2, where it is first used.
- [x] Shared per-cell physics in `cellPhysics.h` (`CellPhysics<T, Descriptor>`: no `Cell`, no
      virtual calls, no exceptions; array templates of Palabos in the same order as its `Cell`
      versions):
  - `collide` / `collideExternal`: forced complete-regularized BGK with an omega function
    (`ConstantOmega`, `CarreauOmega`);
  - `forcedVelocity`, `porousForce`;
  - `completeVelocityNode<dir, ori>` (velocity and Murray openings);
  - `pressureNodeVelocity`, `pressureNodeGradient`, `completePressureNode` (finite-difference
    pressure openings).

  `GuoForcedCarreauDynamics` and `PorousForceFunctional` now call it. Newtonian CPU runs keep
  Palabos' `GuoExternalForceCompleteRegularizedBGKdynamics` (less change to the CPU path).
  Instead, a unit test checks that `CellPhysics::collide` reproduces it bit for bit.
- [x] Unit tests (`tests/unit/cellPhysicsTest.cpp`, target `hemoUnitTests`, `ctest`): 72 checks
      on random cells, bit-identical in single and double:
  - forced velocity;
  - Newtonian `collide` and `collideExternal` against Palabos;
  - Carreau against the previous code;
  - porous force;
  - for all six faces: the velocity node (completion; with the Newtonian and the Carreau bulk
    collision), the pressure-node velocity, and the finite-difference pressure node against
    `PlaneFdBoundaryFunctional3D`.
- [x] Regression harness against the committed binary, byte-identical checkpoints after
      1000 / 400 steps:
  - pipe: Newtonian, flow diverter (a = 100, b = 200), Carreau, Carreau with flow diverter,
    100 Pa outlet with waveform;
  - NAP180 5 M with flow diverter and Murray outlets, Newtonian and Carreau.

  Logs, exit codes and error messages are unchanged, and a checkpoint restart is identical.
- [x] `OpeningHandler`: `progressTime` = `advanceTime` (waveform) + `imposeValues` (lattice);
      `getVelocityScale()` / `getDensity()` give the current boundary value for the GPU update.
- [ ] Refuse opening nodes on a bounding-box edge in GPU mode. This cannot happen with the
      current preprocessor; the check belongs to the GPU mode of Phase 2.

**Findings of Phase 1 for the GPU port.**
- **Pressure-node order.** Palabos' `PlaneFdBoundaryFunctional3D` rebuilds the pressure nodes in
  place, one after another. A node's in-plane neighbours may already be rebuilt when its
  gradient is taken, so the CPU result depends at round-off on the order of the nodes (also
  across blocks, which is probably why CPU process counts agree only to about 3·10⁻¹¹). On the
  GPU the nodes are processed in parallel. That needs two passes: first all node velocities
  into a buffer, then the rebuild from the buffer (in place it would be a race). GPU and CPU
  will therefore agree at round-off, not bitwise, on pressure openings.
- **Double literals.** The Palabos boundary code has `0.5 * force` (double literals). The port
  keeps them so that single-precision results stay bit-identical. They cost double-precision
  arithmetic only on opening nodes.
- **GPU check** (`gpuSmoke` models `hemo` and `carreau`, which use `cellPhysics.h` in the GPU
  functor). 128³, one block, against the CPU dynamics:

  | Model | Single | Double | GPU against CPU |
  |---|---|---|---|
  | HemoFlow Newtonian (`hemo`) | 805 M cell updates/s | 73 M/s | 1·10⁻⁸ (single) / 2·10⁻¹⁷ (double) |
  | Carreau (`carreau`) | 66 M/s | 36 M/s | 3·10⁻⁹ (single) / 4·10⁻¹⁸ (double) |

  - The Newtonian collision is bandwidth-bound in single, as fast as plain BGK.
  - Carreau in single is 12 times slower, because `rheology.h` computes in double. Templating
    it on `T` (Phase 4) is therefore the main speed item for Carreau runs.

### Phase 2: GPU skeleton

**Status: done (22 November 2026).**
- **Code.** `gpuSolver.h/.cpp`, the GPU mode in `hemoFlow.cpp`, `<simulation><device>auto|gpu|cpu`
  and `OpeningHandler::setPlaceholderBCType`.
- **Supported on the GPU.** Fluid (Newtonian and Carreau), walls, outside voxels, velocity and
  Murray openings. Other cases fall back to the CPU with the reason in the log: pressure openings
  (Phase 3), flow diverters with nonzero coefficients (Phase 4), free outflow, LES, openings on a
  box edge.
- **Accuracy.** On the velocity-velocity pipe, GPU double equals GCC double exactly in the output,
  and GPU single matches CPU single to 2·10⁻⁶.
- **Restart and checkpoints.** A GPU restart reproduces velocity, pressure and stress exactly.
  The checkpoint files differ only in non-fluid cells, which hold the populations that arrive
  from beyond the domain. CPU checkpoints continue on the GPU.
- **Several ranks.** With 8 ranks, rank 0 runs the GPU (identical output) and 7 ranks wait with
  their CUDA contexts released. An error after the decision releases them (exit 255).
- **CPU path.** Byte-identical to the committed code (harness, 7 cases).
- **Findings.**
  - Velocity openings: on the CPU, `createLocalBoundaryCondition3D` rebuilds the nodes after
    streaming (`WrappedLocalBoundaryFunctional3D`, `automaticPrepareCollision = false`). The GPU
    does the same: external collide-and-stream, then the node rebuild, then the exchange. A
    placeholder lattice must not carry that processor (its copy kept a dangling pointer).
  - HPC-X: `OMPI_MCA_coll_ucc_enable=0` is set before `MPI_Init`. UCC's CUDA transport crashed
    `MPI_Finalize` on the ranks that had released their context.
- **Memory fallback, tested.** A pipe of 45.5 M voxels (9.1 GB estimated in one single block) runs
  on the CPU with the sparse decomposition, on all ranks. With `<device>gpu</device>` it stops
  with the reason.


- [ ] `toggleBlockingCommunication(true)` before any lattice is created (GPU build only).
- [ ] Execution policy of the GPU build (decided 21 November 2026, README "GPU build"):
  - **Default.** The whole domain is one dense block on the GPU, in single precision.
  - **Fallback.** If that block does not fit the GPU memory, run on the CPU with the sparse
    multi-block decomposition of the CPU build, on all ranks of the launch.

  Implementation:
  - **Decision.** Rank 0 reads the geometry header (Nx, Ny, Nz), estimates the GPU memory of
    the block plus a margin for the CUDA context, and compares it with `cudaMemGetInfo`, then
    broadcasts the decision. GPU-resident per cell: about 180 bytes (two population arrays,
    force, model id, bounce-back links, dynamic-scalar index) plus the multiplier field.
  - **GPU mode.** Rank 0 owns the block and Palabos runs on a sub-communicator of the GPU ranks
    (`global::mpi().init(MPI_Comm)`). The other ranks release their CUDA context, so it does not
    take GPU memory (about 120 MB per rank, measured), and wait for the end without spinning on
    a core.
  - **CPU mode.** The current CPU code path, in the precision of the build.
  - **MPI with HPC-X.** On the CPU path, UCX's CUDA transports must be off (`UCX_TLS=^cuda`).
    With them, managed host memory is staged as if it were GPU memory: 8 ranks took 26 s
    instead of 7.3 s for 1000 pipe steps. The single-GPU mode has no MPI messages, so the
    solver can set `UCX_TLS=^cuda` before `MPI_Init` when the user has not set it. Multi-GPU
    runs (Phase 6) will need the CUDA transports.
  - **Log.** The chosen mode and the memory estimate.
  - **Override.** A config element or command-line option to force CPU or GPU mode.
- [ ] Layout: the GPU mode above, not the CPU's blocks of 24 (Phase 0). Large sparse blocks as
      an alternative are measured in Phase 6.
- [ ] Lattice (2) with the placeholder dynamics, equilibrium initialization and the initial
      opening values. On restart, the checkpoint is loaded into (2). Then copy it to the GPU and
      delete (2).
- [ ] Collision functor: `NoDynamics`, `BounceBack`, fluid (Newtonian first).
- [ ] Output and checkpoints through lattice (1) (next section).
- [ ] First GPU run: the pipe with a velocity inlet and a velocity outlet (no pressure opening),
      compared with the CPU.

### Phase 3: openings (done 22 November 2026)

- [x] **Velocity and Murray openings** (Phase 2).
  - Per-block node lists with the unit profile; before each step a kernel writes
    u = profile × (peak × waveform) into the dynamic scalars, rounded as on the CPU.
  - After streaming, `CellPhysics::completeVelocityNode` rebuilds the nodes, as
    `WrappedLocalBoundaryFunctional3D` does on the CPU.
- [x] **Pressure openings** (`gpuSolver.cpp`).
  - Placeholder `RegularizedDensityBoundaryDynamics >> TRT`, with the imposed rhoBar in its
    dynamic scalar.
  - After streaming, two passes: the velocities of all pressure nodes into a buffer, then the
    rebuild from the neighbour velocities.
  - The neighbour types (fluid, wall or outside, velocity node, pressure node) are resolved once
    on the host.
  - Envelope width 2 in the memory estimate. ω∞ in the condition, as on the CPU.
- [x] Free outflow, `LES 1` and openings on a box edge run on the CPU, with the reason in the log;
      with `<device>gpu</device>` they are an error.
- [x] **Energy check.** A GPU reduction with the velocity of each tag (`GpuSolver::averageEnergy`).
      The `-d` flow rates use the copy to the CPU lattice every 200 steps.
- **Verification** (GPU against CPU, same case):
  - **Pipe** with a 0 Pa outlet, with a 100 Pa outlet and waveform, and Carreau. Double: velocity
    within 3·10⁻¹¹, pressure identical in the output. Single: 3·10⁻⁶.
  - **Full pipe validation** (`pipe.xml`, 80 k steps, single): `pressure_drop.py` gives 78.70 Pa
    and flow rates 1.483 / 1.485·10⁻⁵ m³/s (CPU: 78.70 single, 78.72 double). It took 67 s,
    against 263 / 358 s for 16 CPU cores in single / double.
  - **NAP180 5 M, Murray mode** (velocity inlet, Murray outlets, 0 Pa outlet), flow diverter
    with zero coefficients, Newtonian and Carreau. Double: output identical. Single: 2-4·10⁻⁶.
- **Speed** (NAP180 5 M dense, 400 steps without output):
  - GPU single 2.8 s, i.e. 728 M updates/s of the box;
  - CPU single, 16 ranks, sparse: 9.8 s.

  With output every 200 steps the GPU took 18.6 s. Each output or checkpoint costs about 5 s:
  the copy from the GPU takes 0.1 s, the rest is `writeVTK` / `saveBinaryBlock` on one rank over
  the dense box (see Phase 6). Setup on one rank takes about 15 s at 5 M voxels.

### Phase 4: porous force and Carreau (done 22 November 2026)

- [x] **Porous force** (`gpuSolver.cpp`).
  - At start-up the host lists the porous fluid cells of each block, with their cell index and
    both multipliers (16 bytes per cell in single). This replaces the planned `_AT` processor
    over `porousBox`.
  - Each step a kernel over the list computes the force from the velocity after streaming
    (half-force correction with the force of the previous step) with
    `CellPhysics::porousForce`, as `PorousForceFunctional` does, and writes it into the force
    field for the next collision.
  - Order as on the CPU: collision and streaming, opening nodes, porous force, exchange. On the
    CPU the opening processors come first because they were added first. So on both, the
    finite-difference pressure condition reads its neighbours' velocities with the force of the
    previous step.
  - **Cells.**
    - Opening nodes in the porous medium get no force, as on the CPU, where
      `PorousForceFunctional` skips boundary dynamics.
    - Porous voxels are fluid or opening voxels (`GeometryFile::isPorous`), so any other tag
      stops the start-up with an error.
    - Only bulk cells are listed: the exchange copies the external scalars, so envelope copies
      get the force from their owner block.
  - **Zero coefficients.** The kernel also runs with zero coefficients, as the CPU, which
    writes a zero force every step.
  - **Memory.** The list counts in the GPU memory estimate; NAP180 5 M has 5205 porous voxels.
  - **Cost.** About 30 µs per step, one kernel launch: NAP180 5 M takes 2.82 s for 400 steps,
    against 2.80 s without the force (Phase 3). Folding the force into the collision is not
    worth it.
- [x] **Carreau in the precision of the lattice** (`rheology.h`).
  - `CarreauModelT<R>` is the solver in precision R. `CarreauModel` (`CarreauModelT<double>`)
    keeps the parameters for the configuration and the serialization.
    `CellPhysics::CarreauOmega` converts it to `T`, so single-precision builds solve in float,
    on the GPU and on the CPU. The Carreau output fields (`carreauCellOmega`) use the precision
    of the collision.
  - Newton tolerance: 1e-12 in double, unchanged, so double CPU results stay byte-identical. In
    float it is 64 ulp (7.6·10⁻⁶), as 1e-12 cannot be reached. ω in float is within 1.4·10⁻⁷
    of double from rest to g = 1/dt (unit test, Carreau and Yasuda a = 1.25).
  - **Speed.**
    - 128³ box (`gpuSmoke`): Carreau single runs at 702 M cell updates/s, against 66 before
      and 783 for Newtonian.
    - NAP180 5 M: Carreau takes 2.83 s for 400 steps, the same as Newtonian, so the Newton
      solve is hidden behind the memory traffic. A warm start is not needed.
- **Verification** (GPU against CPU in the same precision; flow diverter a = 100, b = 200, the
  pipe and NAP180 cases of Phase 3):
  - **Pipe** (`short_fd`, 370 porous voxels), Newtonian and Carreau.
    - Single: velocity within 2.8 / 3.5·10⁻⁶, pressure within 1.9 / 2.1·10⁻⁶.
    - Double: pressure identical in the output. The velocity differs only by one-ulp flips of
      the float32 output, at most 2.3·10⁻¹⁰ m/s, away from the device; nvc++ on the CPU against
      GCC: 7·10⁻¹².
  - **NAP180 5 M, Murray mode** (5205 porous voxels).
    - Newtonian double: output identical.
    - Carreau double: pressure identical; velocity one-ulp flips of the output, at most
      7.5·10⁻⁹ m/s, next to the device.
    - Single: velocity within 1.3 / 1.8·10⁻⁶, pressure within 3.7·10⁻⁶.
    - The flow diverter itself changes this flow by up to 0.039 m/s and 35 Pa (zero against
      these coefficients), four orders of magnitude above the GPU-CPU difference in single.
  - **Restart.** A GPU run restarted from a GPU checkpoint with the flow diverter gives the
    output of the uninterrupted run, exactly.
  - **CPU path.** Byte-identical to the committed code (harness, 7 cases). CPU single builds now
    also solve Carreau in float.
  - **Unit tests.** 74 checks: the 72 of Phase 1, plus the double solver bit for bit against the
    double-only one and the float solver against double.
    - GCC: all pass.
    - nvc++ (`-fast`): the double porous-force check differs in the last bit, as before
      Phase 4, because the compiler contracts the test's own reference expression into fused
      multiply-adds differently from the function.
- **Speed** (NAP180 5 M with the flow diverter, single, 400 steps without output):

  | Run | Newtonian | Carreau |
  |---|---|---|
  | GPU | 2.82 s | 2.83 s |
  | CPU, 16 ranks | 7.9 s | 8.1 s |

  With output every 200 steps the GPU took 19.3 s (18.6 s in Phase 3); output is the Phase 6
  item.

### Output and checkpoints

- [ ] At output steps:
  - copy populations and external fields from the GPU into lattice (1) with HemoFlow's own
    loop over `pullPop` / `pullExt` (not `writeBack`, which would replace the dynamics);
  - impose the current opening values on (1);
  - call `writeVTK` unchanged.

  Each output moves the lattice over PCIe and back, about 0.2-0.5 s for NAP180 sparse.
  Prefetching (`cudaMemPrefetchAsync`) can help.
- [ ] Checkpoints: the same copy, then `saveBinaryBlock(lattice (1))`, so the files are the same
      as from a CPU run and a run can be restarted on either build. Test CPU → GPU and
      GPU → CPU restarts.
- [ ] Later, only if output dominates: GPU functionals for the output fields, which then must
      reproduce the moments of each dynamics type.

### Phase 5: verification (done 22 November 2026)

All runs in single precision, GPU against the CPU build in the same precision (GCC, 16 ranks,
sparse decomposition) unless stated otherwise.

- [x] **Pipe** (`tests/verification/0-pipe`).
  - `pipe.xml`, 80 k steps: `pressure_drop.py` gives 78.70 Pa and flow rates
    1.483 / 1.485·10⁻⁵ m³/s, as in Phase 3 (CPU: 78.70 single, 78.72 double). The velocity
    profile is within 6.6·10⁻⁷ m/s of the CPU single profile and within 3.4·10⁻⁵ m/s (3.8·10⁻⁵
    relative) of CPU double. GPU 71 s.
  - Carreau pipe (`carreau.xml` with `--calibrate newtonian.xml`, 48 k steps each). The GPU
    gives the CPU's verification metrics:

    | Run | Pressure gradient | Peak velocity | Profile RMS |
    |---|---|---|---|
    | GPU single | 220.20 Pa/m (+0.13 %) | −0.79 % | 0.72 % |
    | CPU single | 220.19 Pa/m (+0.13 %) | −0.78 % | 0.72 % |
    | CPU double | 220.21 Pa/m (+0.12 %) | −0.79 % | 0.72 % |

    GPU 44 s for either fluid; CPU 150 s Newtonian, 225 s Carreau.
- [x] **Pressure waveform**: the pipe with the cerebral flow waveform at the inlet and a
      100 Pa × waveform outlet, one cycle (0.8 s, 40 k steps). Velocity within 2-5·10⁻⁶
      relative at every output. Pressure differs by a constant 0.0095 Pa: one float32 ulp of
      the density near 1, i.e. the pressure resolution of single precision at this dx and dt,
      the same on the CPU (see Phase 6, output).
- [x] **C0002, Murray mode** (`tests/verification/2-aneurisk_C0002`, voxelized at its
      configured resolution: 1.8 M voxels, 108 k fluid, four openings), one cycle (40 k steps,
      50 outputs).
  - Velocity within about 1·10⁻⁶ relative for most of the cycle. The largest difference,
    2.9·10⁻⁵, is at peak systole (1.0 m/s); it decays to 1·10⁻⁶ by the end of the cycle.
  - Pressure within 1-2 float ulps (0.004-0.008 Pa).
  - `Murray_check.py` gives the same split on both, to all printed digits. Murray outlets
    within 0.04 % of their d³ share at the end of the cycle (0.24-0.65 % at systole). The 0 Pa
    outlet deviates by 18-36 % on both (the smallest outlet, about 2.5 voxels in radius).
  - GPU 148 s, CPU 475 s.
- [x] **NAP180 5 M with the flow diverter** (a = 100, b = 200, Murray mode), one cycle of the
      cerebral waveform (80 k steps, outputs every 0.1 s). Velocity within
      1.2-5.5·10⁻⁶ relative at all eight outputs, pressure within one float ulp (0.018 Pa).
      GPU 569 s, CPU 1643 s.
- **No-FMA build**: not made. GPU double already gives the CPU double output to one ulp of the
  float32 output (Phase 4); bit-identical checkpoints would also need the two-pass pressure
  condition to run in the CPU's order.
- [x] **Decomposition**: several blocks on the one GPU, with the development switch
      `HEMOFLOW_GPU_BLOCK_SIZE=n` (sparse decomposition with blocks of n in GPU mode). Several
      ranks on the GPU are a Phase 6 item (one rank per GPU). This test found two bugs, both
      only with several blocks per process:
  - **Palabos** (`attribute_static`, the copy between GPU blocks on the same process) indexed
    the source block with the destination block's size. Blocks of different sizes therefore
    exchanged populations from the wrong cells; with 16-cell blocks on the pipe the velocity was
    off by up to 35 %, and NAP180 with 64-cell blocks hit an illegal memory access. Fixed in
    the vendored Palabos (local patch, TODO.md section 0). gpuSmoke did not catch it, as its
    blocks have equal sizes.
  - **HemoFlow**: Palabos' GPU exchange copies only the envelope populations that streaming
    within the receiving block cannot provide. HemoFlow rebuilds all populations of an opening
    node after streaming, so envelope copies of opening nodes kept a mix of their own and the
    owner's values. `GpuSolver` now copies all populations of these copies from their owner
    after the exchange (`copyNodeEnvelopes`). With MPI between GPUs (Phase 6) this copy needs
    communication.

  After the fixes, NAP180 (Newtonian and Carreau, 8 and 44 blocks) and three pipe cases
  (velocity only, 100 Pa outlet with waveform, flow diverter with Carreau; 27 blocks) give
  exactly the output of the single block. The single block is unchanged (byte-identical
  checkpoint).
- [x] **Restarts.**
  - GPU → GPU (Phase 4): identical output.
  - CPU → GPU (Phase 2).
  - GPU → CPU: a GPU double checkpoint of the pipe with flow diverter, continued on 16 CPU
    ranks with the sparse decomposition, gives the uninterrupted CPU run: pressure identical,
    velocity within 1.8·10⁻¹² m/s.
- **CPU path** byte-identical to the committed code (harness, 7 cases); unit tests 74/74 (GCC).

### Phase 6: performance and several GPUs

- [ ] Layout on NAP180 20 M: dense (one block per GPU) against sparse with a few large blocks
      (e.g. 96-160): time per step, memory, kernel count (Nsight Systems). Each block overlap
      costs about 30 µs per step (Phase 0), so the optimum is at few blocks. The GPU build
      should default to that layout and warn when a rank has many blocks.
      `HEMOFLOW_GPU_BLOCK_SIZE=n` (development switch, Phase 5) already runs the sparse
      decomposition with blocks of n on the GPU; several blocks give the single block's output.
- [ ] Single precision, the target: validate the GPU build against the CPU double results on
      the pipe, NAP180 with flow diverter (velocities near the device) and the Carreau pipe.
      Palabos stores f - t_i, which helps single precision. On the CPU, single precision
      already matches double on the pipe (Phase 0). If single is not enough somewhere, the
      fallbacks are double on a data-centre GPU, or double only where it matters (e.g.
      pressure output from ρ - 1).
- [ ] Several GPUs: one MPI rank per GPU; the launcher binds devices (`CUDA_VISIBLE_DEVICES`
      from `OMPI_COMM_WORLD_LOCAL_RANK` in a wrapper, or Slurm `--gpus-per-task=1`), otherwise
      all ranks of a node use GPU 0. Print the device of every rank at start-up. Use CUDA-aware
      MPI (the HPC SDK's HPC-X) with `USE_CUDA_MALLOC`; GPUDirect RDMA if the cluster has it.
      The exchange does not overlap the computation, so each GPU should hold several million
      cells; small cases scale over few GPUs only.
      `GpuSolver::copyNodeEnvelopes` (Phase 5) copies the envelope copies of opening nodes from
      their owner block on the same process only; across ranks it needs an MPI exchange of
      these cells after Palabos' exchange.
- [ ] Output in GPU mode: about 5 s per output at 5 M voxels, against about 1 s on 16 CPU ranks,
      because `writeVTK` runs on one rank over the dense box. Options:
  - compute the output fields on the GPU (velocity, pressure, stress, viscosity, with the
    moments of each tag);
  - write only the fluid voxels;
  - let the waiting ranks write.

  Pressure in single precision: the output takes it from ρ, whose float32 resolution near 1
  is one ulp, 0.004-0.01 Pa at the verification resolutions (Phase 5). Writing it from
  ρ̄ = ρ − 1, which Palabos stores, would remove this limit, also for CPU single builds.

  Setup on one rank (about 15 s at 5 M voxels) could use the same help.
- [x] GPU build and launch documented in the README (22 November 2026).

## Flow diverters and coils

The GPU port handles any porous region described by the multiplier field. Flow diverters
(shells) and coils (volumes, once the preprocessor voxelizes them as solids, TODO.md section 2)
therefore run the same way; for a coil the porous box covers the sac, which makes the porous
processor more expensive but still small next to the collision. The porous models still to
come have to be written in the shared porous-force function of Phase 1, so that both builds
get them:

- device-level, area-weighted coefficients ([TODO-FD-INHOMOGEN.md](TODO-FD-INHOMOGEN.md));
- Kozeny-Carman coils;
- semi-implicit forcing for high resistance;
- a Darcy term scaled with the Carreau viscosity.

## Hardware and expected speed

- **Memory traffic.** LBM is limited by memory bandwidth. A cell update moves about 370 bytes
  in double: populations read and written (304), force read and written (48), model id,
  bounce-back links and dynamic-scalar index (16). In single precision it is about 190 bytes.
  Every allocated cell counts (also wall and outside cells in a block), as in the CPU figure of
  about 25 M updates/s (NAP180 20 M sparse, 300 steps in 44.5 s on 20 cores).
- **RTX A1000** (workstation; GA107, 8 GB, 192 GB/s, 50 W).
  - Measured with a simple Guo-forced BGK on one block (Phase 0): 0.64-0.80 G updates/s in
    single (74-80 % of the bandwidth) and 0.14-0.17 G/s in double (arithmetic-bound).
  - HemoFlow's complete regularized collision has several times the arithmetic. In single it
    should stay close to bandwidth-bound; in double it would drop well below 0.1 G/s.
  - Carreau costs more than the +6 % on the CPU.
- **Data-centre GPUs.** A100 (1.5-2 TB/s, full-speed double) bound at about 4-5 G updates/s in
  double, probably half of that in practice; H100 about twice that.
- **Memory.** NAP180 20 M voxels dense needs about 7.5 GB in double (does not fit the A1000);
  sparse with blocks of 24 and envelope 2 needs about 2.3 GB; larger blocks need more.

## Prerequisites on the workstation (21 November 2026)

- GPU: NVIDIA RTX A1000, 8 GB, compute capability 8.6. Driver 595.99 (CUDA 13.2) loaded,
  `nvidia-smi` works.
- Proprietary kernel module, addressing mode "None": no HMM, so managed memory only.
- CUDA 12.0 `nvcc` (Ubuntu package): not used.
- NVIDIA HPC SDK 26.11 in `/opt/nvidia/hpc_sdk/Linux_x86_64/26.11`: `nvc++`, CUDA 13.4,
  HPC-X Open MPI 5.0.10 with CUDA support, Nsight Systems and Compute in `profilers/13.4`.
  - The paths are set at the end of `~/.bashrc`, after its "not interactive: return" guard, so
    scripts and non-interactive shells must set them themselves.
  - With the SDK first in `PATH`, `mpirun` is HPC-X's; GCC builds against the system Open MPI
    must be started with `/usr/bin/mpirun`.
- System OpenMPI 4.1.6 (`/usr/bin/mpirun`) built without CUDA support, used by the GCC build.
- GCC 13.3 (the `nvc++` host compiler and its C++ library), 24 cores, 125 GB RAM.

## Effort (rough)

- Phase 0: done (21 November 2026).
- Phase 1: done (21 November 2026).
- Phases 2-3: 2 weeks (first GPU pipe run after about one week).
- Phase 4 and output/checkpoints: 1 week.
- Phase 5: 1 week.
- Phase 6: 1 week on the workstation, plus cluster time for several GPUs.

## Open questions

1. Which GPUs for production (cluster, model, double-precision speed), and is CUDA-aware MPI
   available there? The A1000 is a development GPU.
2. ~~Double or single precision?~~ Single, decided 21 November 2026; double only if single is
   not enough (measured factor of about 5 in speed on the A1000).
3. Are free outflow and LES needed on the GPU?
4. Phase 1 changes the CPU code (bit-identical results, checked): acceptable as the first step?
