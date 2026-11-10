# HemoFlow: open tasks

Status as of November 2026, after the review and fix round on `dev-refactor`
(see the README for what the code does now). Items are ordered by priority
within each section.

## How to resume

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build -j 8
conda env create -f preprocessor/environment.yml          # env "lbmpre", Python 3.12
conda create -n hemoflow-post -c conda-forge python=3.12 pyvista pandas matplotlib
```

Quick regression check: `tests/verification/0-pipe` (README, test case 1).
Reference: 80.6 Pa pressure drop, inlet/outlet flow 1.537/1.539·10⁻⁵ m³/s.

## 0. Palabos upgrade

The vendored Palabos was moved from the 2024 snapshot to upstream master (after
release 2.3.0; adds GPU-oriented accelerated lattices, needs C++20).

- [x] Replace `palabos/` with the new version
- [x] Review upstream changes that touch the models HemoFlow uses (collision templates
      only switched to an identical `constexpr` constant; the removed `addGuoForce` call
      fixes double forcing in the Smagorinsky variant, which HemoFlow does not use)
- [x] Adapt the build: C++20 default, skip the MPI C++ bindings
- [x] Local patch: `inline` added to `replaceAll()` in
      `palabos/src/acceleratedLattice/acceleratedProcessors3D.hh` (non-inline function in a
      header, so programs with several source files fail to link). Report upstream and drop
      the patch once fixed there.
- [x] Build with GCC 13 (system) and GCC 15 (conda)
- [x] Regression-test against the old Palabos: bit-identical results for the pipe,
      NAP180 with flow diverter, C0002 Murray split, 100 Pa outlet with waveform, free
      outflow, velocity outlet, checkpoint restart and sparse decomposition; pipe
      validation unchanged (80.6 Pa); same speed
- [x] Remove `palabos_old/` (the previous version is in git history, commit 9e06177)
- [x] Trim the vendored copy to what HemoFlow builds: `src/`, `externalLibraries/`,
      `dco/`, `COPYING`, `README.md`, `CHANGELOG.md` (examples, docker, utility,
      coupledSimulators, upstream CMake and CI files removed; 49 MB -> 25 MB). A pinned
      submodule remains an option (section 8).

## 1. Correctness and verification

- [ ] **Sparse decomposition differs from dense in one case.** `<sparse>true</sparse>`
      is bit-identical to dense for C0002 (block size 16) and NAP180 (block size 24),
      but for NAP180 (`testcase/`, 1 M elements, `test_Aneurysm.xml`) at block size 16
      a round-off-level difference appears at step 70 on fluid nodes next to the 0 Pa
      outlet (label 11, y = 0 face), and grows to 0.06 m/s by step 400. Dense runs on
      2 and 8 processes are identical, so the cause is the dropped blocks. Suspect the
      finite-difference pressure boundary (`createInterpBoundaryCondition3D`) reading
      envelope cells of a dropped neighbour block. Once fixed, make sparse the default
      (it saves 40–70 % of the cells in aneurysm geometries).
- [x] **Voxelizer dropped the partial top layer in y and z.** Fixed by the new VTK
      voxelizer, which covers the whole bounding box (phase 1, November 2026).
- [ ] Full-length verification runs not yet done with the current code:
      `1-flow_split` (2 cycles, ~1–2 h on 20 cores), `2-aneurisk_C0002`, and the
      `testcase/test_*.xml` variants on the 5 M geometry. Only short runs were checked.
- [ ] Re-run the V&V notebook (`doc/VV/hemoflow_VV.ipynb`) and the EasyVVUQ campaigns
      with the current code; its numbers come from older versions.
- [ ] Not tested: the HDF5 writer, the build with Ubuntu's `libhdf5-openmpi-dev`,
      checkpoint restart with a different number of MPI processes, and the ParaView WSS
      recipe in `doc/paraview_wss_calc.txt` (the stress output itself was checked
      against the pipe pressure gradient).

## 2. Flow diverter and coil modelling

- [ ] The porous coefficients act per voxel layer (force −(a u + b u|u|)/dx), so the
      resistance of a device depends on its voxelized thickness and therefore on dx and
      on the voxelization. Normalize by the physical device thickness, or specify
      permeability (Darcy K, Forchheimer C) instead. Note: the new ray-cast shell
      voxelization (phase 1) gives ~25 % fewer device voxels than the old one, whose three
      projections were shifted against each other; coefficients tuned on old geometry
      files need re-checking.
- [ ] Ray-crossing voxelization along the axes can let diagonal (D3Q19) lattice links
      cross the device without touching a porous voxel; consider a conservative
      (triangle–voxel overlap) voxelization.
- [ ] Make per-voxel coefficients the main path for deployment-simulation results: the
      preprocessor writes `linear`/`quadratic` (from `<stent_mesh_base>values.vtp` when
      `INHOMOGEN = True`), but the solver ignores them. Read them in `hemoFlow.cpp` and
      use them in `PorousForceFunctional`. The shipped `PED_5x20_values.vtp` carries the
      deployment descriptors `Density`, `Expansion` and `Force`, not `linearCoeff` /
      `quadraticCoeff`; a mapping from descriptors to coefficients is still needed (the
      interpolation itself was tested with test arrays).
- [ ] Coils: volumetric porous regions (e.g. Kozeny–Carman from the packing density)
      instead of shells; needs a solid (not shell) voxelization of the coil region.
- [ ] High-resistance devices will need a semi-implicit treatment of the porous force
      for stability.
- [ ] With `LES 1` the Smagorinsky dynamics has no force term, so devices have no effect.

## 3. Boundary conditions and openings

- [ ] Openings with overlapping bounding rectangles on the same face are refused.
      Supporting them needs per-node boundary instantiation (velocity BCs are
      dynamics-based and could be set per node; the interp pressure BC is a data
      processor per domain).
- [x] Openings on edges/corners of the bounding box: openings now come from the caps of
      the surface and only the wall in front of each cap is opened, so openings near
      edges and corners and vessels close to a face with an opening work (tested with
      synthetic cases, phase 1). Openings must still be flat cuts on bounding-box faces.
- [x] Velocity `parameter` stays the mean of an ideal paraboloid (peak = 2 × parameter),
      so existing configs keep their results; the log reports the actual mean over the
      discrete opening (decided November 2026).
- [ ] Optional warm-up ramp and a steady-run option (end time or convergence criterion)
      instead of a constant waveform file whose length defines a "cycle". Not needed
      for now (decided November 2026).
- [ ] Non-Newtonian rheology (Carreau parameters exist only as literature references).

## 4. Performance and scalability

- [ ] The whole geometry file (flags, stent, coefficients) is loaded on every MPI
      process. Read the `.npz` on one rank and distribute it (decided November 2026:
      keep the `.npz` format, no HDF5 input).
- [ ] `OpeningHandler` scans the full domain once per opening at start-up.
- [x] Voxelizer replaced by VTK (`vtkPolyDataToImageStencil` for the vessel, ray casting
      for flow diverters): NAP180 at 20 M voxels in ~4 s instead of ~53 s.
- [x] Recursive opening detection replaced: openings come from the surface caps; the
      connectivity checks use `scipy.ndimage.label`.

## 5. Code structure (C++)

- [ ] Split `main()` in `hemoFlow.cpp` (config, geometry, lattice, boundary setup, time
      loop) and replace the globals with a configuration/simulation struct.
- [ ] Remove `using namespace std` from headers; use `std::` consistently.
- [ ] Decide on the unused writers in `io.cpp` (`writeHDF5`, `writeNPZ`): wire them up as
      output options or remove them.
- [ ] Bump the version string (still `v0.31`).

## 6. Preprocessor (Python)

- [ ] Turn `main.py` into functions plus a CLI; move the constants (`SI_FACTOR`,
      `DEBUG_MODE`, `INHOMOGEN`, `DISTANCE_RADIUS_FACTOR`) into the JSON config or CLI
      options.
- [ ] Remove the ignored legacy keys (`cutWidth`, `distance`) from the shipped configs.
- [x] `slice.py`, `perimeter.py`, `util.py` removed with the old voxelizer (no more
      shadowing of the built-in `slice`, no more "watertight" messages).

## 7. Testing and CI

- [ ] CI (e.g. GitHub Actions): build with GCC 13 and GCC 15, with and without HDF5.
- [ ] Short solver regression test: coarse pipe, a few thousand steps, check the flow
      balance and the pressure drop against stored values.
- [ ] Preprocessor tests (pytest) on the pipe: opening count and labels, no fluid voxel
      next to an outside voxel, flow diverter voxels only in the lumen, origin aligned
      with the STL.
- [ ] Solver checks used during the review that are worth automating: non-zero pressure
      outlet (100 Pa imposed → 100 Pa measured), checkpoint restart reproduces an
      uninterrupted run, Murray outlets receive their d³ share (`-d` output).

## 8. Repository hygiene

- [ ] Track Palabos as a pinned git submodule or CMake `FetchContent` instead of a
      vendored snapshot (keeping the `replaceAll` patch, see section 0); record the
      HighFive and cnpy versions (cnpy has a local patch for ZIP64 `.npz` files written by
      Python ≥ 3.12).
- [ ] Remove or update stale scripts: `scripts/cartesius_job.sh` (Cartesius is
      decommissioned), `utils/update_times.sh`, the duplicate `compress_output.py`
      (root and `utils/`), the empty `preprocessor/run.sh`.
- [ ] Campaign scripts in `tests/sensitivity` and `tests/uncertainty` hard-code
      `SOFTWARE_PATH`; make it an argument or environment variable. The aneurysm
      campaigns depend on LBMpost, which is not public.
