# Running HemoFlow on NVIDIA GPUs (plan)

Status: Phase 0 (toolchain) done on 21 November 2026, see "Phase 0 results"; Phase 1 not
started. First investigated on 10 November 2026 (that version is in git, commit bca490e);
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
   force into the collision stays as an optimization.
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
  10 % on one rank (managed memory). Several ranks sharing one GPU do not scale on the CPU side:
  8 ranks took 26 s, as long as 1 rank. This does not matter for production (one rank per GPU),
  only for tests of several ranks on the workstation.

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

### Phase 1: preparation on the CPU (results must stay bit-identical)

- [ ] Split `main()` (TODO.md section 5) at least into: configuration, geometry, lattice
      builder, boundary setup and time loop. The lattice builder and boundary setup take a
      "placeholder dynamics" switch, so that lattices (1) and (2) come from the same code.
- [ ] Shared per-cell physics as array-based inline functions in a new header (no `Cell`, no
      virtual calls, no exceptions), with the same order of operations as now:
  - forced complete-regularized BGK with a given ω;
  - the Carreau ω from Π^neq (`rheology.h` is already plain C++);
  - the porous force.

  A HemoFlow bulk dynamics class (Newtonian and Carreau; it replaces the direct use of
  `GuoExternalForceCompleteRegularizedBGKdynamics` and generalizes `GuoForcedCarreauDynamics`)
  and `PorousForceFunctional` call them. Check with the regression harness that the CPU stays
  bit-identical: pipe, NAP180 with flow diverter, C0002 Murray, pressure waveform, Carreau pipe.
- [ ] Array-based ports of the two Palabos boundary steps HemoFlow uses:
  - completion of a velocity-opening node (complete-equilibrium Π^neq and regularization, as in
    `RegularizedVelocityBoundaryDynamics` around the bulk dynamics);
  - the finite-difference pressure completion of one node from its neighbours' velocities
    (`PlaneFdBoundaryFunctional3D`).

  CPU unit tests compare them with the Palabos classes on random populations, for all six
  faces.
- [ ] `OpeningHandler`: split `progressTime` into advancing the waveform and imposing the values,
      so that the host can impose on the GPU lattice and on the output lattice separately.
- [ ] A check that refuses opening nodes on a bounding-box edge in the GPU build (this cannot
      happen with the current preprocessor).

### Phase 2: GPU skeleton

- [ ] `toggleBlockingCommunication(true)` before any lattice is created (GPU build only).
- [ ] Layout of the GPU build: dense (one block per rank) or large sparse blocks by default, not
      the CPU's blocks of 24 (Phase 0); tuned in Phase 6.
- [ ] Lattice (2) with the placeholder dynamics, equilibrium initialization and the initial
      opening values. On restart, the checkpoint is loaded into (2). Then copy it to the GPU and
      delete (2).
- [ ] Collision functor: `NoDynamics`, `BounceBack`, fluid (Newtonian first).
- [ ] Output and checkpoints through lattice (1) (next section).
- [ ] First GPU run: the pipe with a velocity inlet and a velocity outlet (no pressure opening),
      compared with the CPU.

### Phase 3: openings

- [ ] Velocity and Murray openings: per-block node lists with the unit profile; a per-step
      kernel writes u = profile × (peak × waveform) into the dynamic scalars of the bulk and
      envelope copies; the functor completes the node (Phase 1 function) and then applies the
      bulk collision.
- [ ] Pressure openings: GPU processor over a node list per block (not per row as on the CPU,
      which would be one full-block sweep per row).
  - Per node it uses the imposed ρ (dynamic scalar), the boundary velocity from the known
    populations, one-sided normal and central in-plane velocity differences, and
    complete-regularizes the populations.
  - The neighbour velocity depends on the tag: 0 for walls, j/ρ + F/2 for fluid, the boundary
    velocity for pressure nodes, the imposed u for velocity nodes.
  - Envelope width 2, as on the CPU.
  - The functor then applies only the bulk collision. Keep ω∞ in the condition, as on the CPU
    (TODO.md section 3), unless both builds change.
- [ ] Refuse free outflow (type 4) and `LES 1` in the GPU build with a clear message.
- [ ] Energy check and `-d` flow rates: GPU reductions (`transform_reduce`) with HemoFlow's
      velocity definition (j/ρ + F/2 for fluid, imposed u at velocity nodes). For the debug
      flow rates, a copy to lattice (1) every 200 steps is also acceptable.

### Phase 4: porous force and Carreau

- [ ] Porous force: `_AT` processor over `porousBox` with the existing
      `MultiTensorField3D<T,2>` of multipliers. It uses the Phase 1 function, applies to fluid
      tags only, in the bulk, and writes the force into `ext`. Same order as the CPU: force from
      the post-streaming velocity, then the next collision with that force.
- [ ] Alternative, only if profiling asks for it: fold the force into the collision with a
      5-scalar descriptor (force plus multipliers, about 9 % more traffic for every cell), or
      carry the multipliers in dynamic scalars of a separate placeholder model. With the
      in-collision variant the first step after a restart differs slightly from the CPU.
- [ ] Carreau: local ω in the functor. `rheology.h` computes in double with a 1e-12 tolerance.
      For the single build, template it on `T` with a tolerance near float resolution;
      otherwise every cell runs a double-precision Newton solve with `pow` at 1/64 speed. The
      warm start from the previous ω (TODO.md section 4) would need one more external field.

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

### Phase 5: verification

- [ ] Pipe (`tests/verification/0-pipe`): GPU against CPU, pressure drop (78.7 Pa) and profile;
      Carreau pipe (`carreau.xml`).
- [ ] NAP180 with the flow diverter, C0002 Murray, a pressure waveform: GPU against CPU.
      Newtonian runs of the validation build without FMA may be bit-identical; otherwise
      report the differences, which are expected near round-off. Carreau differs at round-off
      level because the device `pow` differs from glibc's.
- [ ] Decomposition: 1, 2 and 4 ranks on the one GPU (without `USE_CUDA_MALLOC`, as the system
      MPI is not CUDA-aware) and different block sizes must give bit-identical results, as on
      the CPU.
- [ ] Checkpoint restart reproduces an uninterrupted GPU run.

### Phase 6: performance and several GPUs

- [ ] Layout on NAP180 20 M: dense (one block per GPU) against sparse with a few large blocks
      (e.g. 96-160): time per step, memory, kernel count (Nsight Systems). Each block overlap
      costs about 30 µs per step (Phase 0), so the optimum is at few blocks. The GPU build
      should default to that layout and warn when a rank has many blocks.
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
- [ ] Document the GPU build and launch in the README.

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
- Phase 1: 3-5 days.
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
