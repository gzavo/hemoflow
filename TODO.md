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
Reference: 78.7 Pa pressure drop, inlet/outlet flow 1.483/1.485·10⁻⁵ m³/s (unchanged by phases 2–3).

## Progress

Phases agreed for sections 1, 3 and 4 (November 2026): 1 preprocessor rework (done),
2 solver geometry distribution and per-node opening boundaries (done), 3 sparse
decomposition bug (done), 4 non-Newtonian (Carreau–Yasuda) blood (done), 5 full verification.
After phases 3 and 4: spatially varying flow diverter resistance (section 2).

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

- [x] **Results depended on the block decomposition near pressure outlets** (phase 3,
      November 2026). Seen first as sparse ≠ dense for NAP180 at block size 16 (0.06 m/s
      after 400 steps next to the 0 Pa outlet), but not sparse-specific: a dense run with
      16-node blocks gave the same result as the sparse one, and a dense run of a small
      case on 1 vs 4 processes differed by 7·10⁻³ m/s (peak 0.66 m/s) once a block
      boundary cut the outlet. Cause: the finite-difference pressure condition
      (`createInterpBoundaryCondition3D`) reads two cells deep, but the lattice had a
      one-cell envelope. Fix: envelope width 2 when the case has a pressure opening. Now
      sparse and dense are bit-identical (15 configs: NAP180, branch and C0002 at block
      sizes 12–32), and process counts 1–24 agree to ≤ 3·10⁻¹¹ m/s. Sparse is now the
      default, block size 24 (NAP180 20 M: 18.6 % of the cells, 2.6× faster than dense
      on 20 processes). The pipe validation on 10 processes had its outlet cut at z = 12:
      the fix removes an artefact of up to 0.03 m/s and 27 Pa in the last three cells
      before the outlet; the pressure drop stays 78.7 Pa.
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
- [ ] **Spatially varying resistance (planned after phases 3 and 4).** Different parts of a
      deployed device have different braid densities, so the resistance must vary per voxel.
      Plan: the per-voxel `linear`/`quadratic` arrays of the geometry file become dimensionless
      multipliers of the `<flowdiverter>` `linCoeff`/`quadCoeff` (the uniform case stays 1.0);
      the main process reads them and distributes them like the other fields (geometry.cpp,
      `scatterToField`) into a two-component field that replaces the 0/1 porosity field, and
      `PorousForceFunctional` uses a(x), b(x). Check that the uniform case stays bit-identical.
      Needs from the physics side: the mapping from the deployment descriptors to the
      multipliers, normalized by the voxelized device thickness (see the first item above).
      The shipped `PED_5x20_values.vtp` carries the descriptors `Density`, `Expansion` and
      `Force`; the preprocessor's `INHOMOGEN` interpolation expects `linearCoeff` /
      `quadraticCoeff` (tested with test arrays).
- [ ] Coils: volumetric porous regions (e.g. Kozeny–Carman from the packing density)
      instead of shells; needs a solid (not shell) voxelization of the coil region.
- [ ] High-resistance devices will need a semi-implicit treatment of the porous force
      for stability.
- [ ] With `LES 1` the Smagorinsky dynamics has no force term, so devices have no effect.
- [ ] With the Carreau model the porous coefficients stay constant. A Darcy term that scales
      with the local viscosity (a = μ / permeability) would follow the fluid; decide together
      with the spatially varying resistance above.

## 3. Boundary conditions and openings

- [x] Openings with overlapping bounding rectangles on the same face: boundary
      conditions are now instantiated row by row on the opening voxels only (phase 2,
      November 2026; tested with a synthetic case whose outlet rectangles overlap).
      Results of existing cases are bit-identical.
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
- [x] **Non-Newtonian rheology** (phase 4, November 2026): `<fluid><model>carreau</model>`,
      Carreau–Yasuda with blood defaults after Cho and Kensey (1991) and an optional Yasuda
      exponent `a` (default 2 = Carreau). Every cell solves for its own shear rate and
      viscosity at every step (`rheology.h`: Newton with a bisection fallback, round-off
      accurate; `nonNewtonian.h`: the forced regularized BGK with the local ω). Palabos only
      has plain Carreau, and its standalone Carreau collisions have no force term, so the
      flow diverter would be ignored. Stress output uses the local viscosity; new output
      field `viscosity`. Verified on the pipe (`tests/verification/0-pipe/carreau.xml`,
      README test case 1); Newtonian runs are bit-identical to before.
- [ ] Carreau: the finite-difference pressure boundary rebuilds the stress of its cells with
      one ω (the one at μ∞), so the viscosity in the outlet layer is not the local one.
      Fixing it needs a pressure boundary condition on the HemoFlow side that takes ω per
      cell (Palabos' local pressure condition, which would, diverged in the pipe case).
- [ ] Carreau: velocity openings impose a parabolic profile. A developed Carreau–Yasuda
      profile (from the same semi-analytical solution as `carreau_profile.py`) would shorten
      the entrance region.

## 4. Performance and scalability

- [x] The geometry file is read on the main process only and each process receives
      its own blocks (phase 2, `geometry.h/.cpp`); the unused `linear`/`quadratic` arrays
      are not loaded. NAP180 20 M, sparse, 8 processes: 572 -> ~230 MB per process.
      Results bit-identical.
- [x] Opening voxels come from one scan of the label array on the main process
      (instead of one full-domain scan per opening on every process).
- [x] Voxelizer replaced by VTK (`vtkPolyDataToImageStencil` for the vessel, ray casting
      for flow diverters): NAP180 at 20 M voxels in ~4 s instead of ~53 s.
- [x] Recursive opening detection replaced: openings come from the surface caps; the
      connectivity checks use `scipy.ndimage.label`.
- [ ] Carreau: the shear-rate solve takes 3–4 Newton iterations per cell and step (+6 % per
      step for NAP180 at 20 M voxels, about +60 % in small cases such as the pipe). Starting
      from the previous step's ω (one stored scalar per cell, a change of the lattice
      descriptor and the checkpoint format) would need 1–2.
- [ ] The two-cell envelope needed by the finite-difference pressure condition costs
      20 % (dense) to 50 % (sparse) per time step (NAP180 20 M, 20 processes: sparse 29.8
      -> 44.5 s for 300 steps). Palabos' local pressure condition needs only one cell but
      diverged in the pipe case (NaN at step 1400) and changes NAP180 by 0.7 m/s, so it is
      not an option. Possible: a pressure condition that reads only one cell deep
      (e.g. Zou/He or an anti-bounce-back outlet) and re-verify against the pipe.

## 5. Code structure (C++)

- [ ] Split `main()` in `hemoFlow.cpp` (config, geometry, lattice, boundary setup, time
      loop) and replace the globals with a configuration/simulation struct.
- [ ] Remove `using namespace std` from headers; use `std::` consistently.
- [ ] Decide on the unused writers in `io.cpp` (`writeHDF5`, `writeNPZ`): wire them up as
      output options or remove them. `writeHDF5` converts the stress with the Newtonian ω,
      so it would need the Carreau path of `writeVTK`.
- [ ] Bump the version string (still `v0.31`).

## 6. Preprocessor (Python)

- [ ] Turn `main.py` into functions plus a CLI; move the constants (`SI_FACTOR`,
      `DEBUG_MODE`, `INHOMOGEN`, `DISTANCE_RADIUS_FACTOR`) into the JSON config or CLI
      options.
- [ ] Remove the ignored legacy keys (`cutWidth`, `distance`) from the shipped configs.
- [x] `slice.py`, `perimeter.py`, `util.py` removed with the old voxelizer (no more
      shadowing of the built-in `slice`, no more "watertight" messages).
- [x] Everything in the coordinates of the input STL (November 2026): the preprocessor
      writes the labels as `<base>geometry.vti` (and `.nrrd` with spacing and origin) and
      stores `lengthUnit` (= `SI_FACTOR`) in the `.npz`; the solver writes its grid in
      those units by default (`<outputCoordinates>stl|m</outputCoordinates>`). Note that
      the STL unit is still the `SI_FACTOR` constant (see the first item).

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
      `SOFTWARE_PATH`; make it an argument or environment variable. The pipe campaign's
      `pressure_drop.py` also assumes an output grid in meters starting at (0, 0, 0);
      take spacing and origin from the geometry file as `tests/verification/0-pipe` does. The aneurysm
      campaigns depend on LBMpost, which is not public.
