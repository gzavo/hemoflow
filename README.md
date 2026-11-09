# HemoFlow

HemoFlow is a lattice Boltzmann (LBM) blood-flow solver for patient-specific
cerebral vasculature, with a focus on intracranial aneurysms and their
endovascular treatment. It takes a segmented vessel surface (STL) and its
centerline, voxelizes it, imposes physiological inflow/outflow boundary
conditions and runs a transient, MPI-parallel CFD simulation. Deployed
devices such as flow diverters can be included as **porous volumes** that
add a Darcy–Forchheimer momentum sink to the flow, so that the
haemodynamic effect of a device can be compared against the untreated case.

The solver is built on the [Palabos](https://palabos.unige.ch) library
(D3Q19 lattice, regularized BGK collision with Guo forcing). Palabos is an
external dependency that ships with this repository in `palabos/`; it is not
developed here.

* Language: C++ (solver), Python (preprocessing and post-processing scripts)
* Parallelism: MPI (Palabos block decomposition)
* License: AGPL v3.0

---

## Contents

1. [Workflow overview](#workflow-overview)
2. [Repository layout](#repository-layout)
3. [Building the solver](#building-the-solver)
4. [Setting up the preprocessor](#setting-up-the-preprocessor)
5. [Preparing a geometry (preprocessor)](#preparing-a-geometry-preprocessor)
6. [Running a simulation (solver)](#running-a-simulation-solver)
7. [Test cases and verification](#test-cases-and-verification)
8. [Known issues and limitations](#known-issues-and-limitations)
9. [Troubleshooting](#troubleshooting)

---

## Workflow overview

```
 vessel surface (.stl, mm)  ─┐
 centerline (.vtp, VMTK)    ─┼─► preprocessor/main.py ──► <name>c.npz ──┐
 flow diverter shell (.stl) ─┘    (+ JSON config)         <name>geometry.nrrd (for inspection in 3D Slicer)
                                                                         │
 inflow waveform (.txt) ─────────────────────────────────┐               │
 simulation config (.xml) ───────────────────────────────┴─► hemoFlow ◄──┘
                                                               │  (mpirun)
                                                               ▼
                                                  <outputDir>/vtkNNNNNN.vti  ──► ParaView / pyvista
```

1. **Preprocessing** (`preprocessor/`, Python) voxelizes the vessel surface at
   the requested resolution, adds a one-voxel wall layer, opens the
   inlet/outlets on the faces of the bounding box, matches every opening to
   a centerline end point (radius and direction), labels the openings and,
   optionally, voxelizes a flow diverter into a porous mask. The result
   goes into a single compressed NumPy archive (`.npz`).
2. **Simulation** (`hemoFlow`, C++/MPI) reads an XML configuration and the
   `.npz` file, sets up the lattice, boundary conditions and porous forcing,
   and integrates a given number of cardiac cycles. VTK image files are
   written during the last cycle.
3. **Post-processing** is done in ParaView or with small `pyvista` scripts
   (see `tests/verification/*/*.py`).

## Repository layout

| Path | Content |
|---|---|
| `hemoFlow.cpp` | Main program: reads config and geometry, sets up the lattice, boundary conditions and the time loop |
| `opening.h/.cpp` | `OpeningHandler`: inlet/outlet boundary conditions, velocity/pressure profiles, time-scale (waveform) functions |
| `porous.h` | Porous (flow diverter) momentum sink as Palabos data processors |
| `io.h/.cpp` | VTK output (and unused HDF5/NPZ writers), directory helpers |
| `globals.h`, `helper.h` | Lattice/dynamics typedefs, compile-time switches, geometry labels, small vector types |
| `CMakeLists.txt` | Build script (builds Palabos as a static library and links `hemoFlow`) |
| `template_config.xml` | Annotated solver configuration template (also in `testcase/`) |
| `preprocessor/` | Python voxelization and opening-detection tools (`main.py` is the entry point) |
| `testcase/` | Example aneurysm geometries (NAP126, NAP180) with a PED 5×20 flow diverter and example configs |
| `tests/verification/` | Verification cases: straight pipe (Poiseuille), bifurcation flow split (Murray), Aneurisk C0002 |
| `tests/sensitivity/`, `tests/uncertainty/` | [EasyVVUQ](https://github.com/UCL-CCS/EasyVVUQ) sensitivity and uncertainty campaigns |
| `doc/VV/` | Verification and validation notebook (`hemoflow_VV.ipynb`) and its data |
| `doc/paraview_wss_calc.txt` | How to compute wall shear stress in ParaView from the stress output |
| `palabos/` | Vendored Palabos library (external dependency) |
| `external/` | Vendored [HighFive](https://github.com/BlueBrain/HighFive) (HDF5 C++ wrapper) and [cnpy](https://github.com/rogersce/cnpy) (`.npz` reader) |
| `scripts/`, `utils/` | Cluster job/environment scripts and small helpers (output compression) |

---

## Building the solver

### Requirements

| Dependency | Notes |
|---|---|
| C++ compiler | GCC or Clang with C++14 support |
| CMake | ≥ 3.13 |
| MPI | e.g. OpenMPI or MPICH |
| zlib | Needed to read the compressed geometry files |
| Parallel HDF5 | Optional. Enables the HDF5/XDMF writer in `io.cpp`, which the solver does not call at the moment. Without it the build simply skips that writer. |

Palabos, HighFive and cnpy come with the repository; nothing else needs
to be downloaded.

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j 8
```

The executable is `build/hemoFlow`. The first build compiles all of Palabos
and takes a few minutes; later builds only recompile the HemoFlow sources.

Tested toolchains (November 2026):

* **Ubuntu 24.04 system packages**: GCC 13.3, CMake 3.28, OpenMPI 4.1.6, no
  HDF5 (`sudo apt install build-essential cmake openmpi-bin libopenmpi-dev zlib1g-dev`).
  Add `libhdf5-openmpi-dev` for the HDF5 writer (not tested).
* **conda-forge**, no root access needed: GCC 15.3, CMake 4.4, OpenMPI 5.0, parallel HDF5 2.2:

  ```bash
  conda create -n hemoflow-build -c conda-forge --override-channels \
      cxx-compiler cmake make openmpi "hdf5=*=mpi_openmpi*" zlib
  conda activate hemoflow-build
  ```

On HPC clusters, load a compiler, CMake and MPI module and build as above.
`scripts/snellius_env.sh` lists the modules used on SURF Snellius in 2021;
current module versions will differ.

### Build options

| Option | Default | Meaning |
|---|---|---|
| `CMAKE_BUILD_TYPE` | `Release` | `Release` (`-O3`), `RelWithDebInfo`, or `Debug` (`-O0 -g -DPLB_DEBUG`, executable `hemoFlowd`) |
| `ENABLE_MPI` | `ON` | Build with MPI |
| `BUILD_HDF5` | `ON` | Use parallel HDF5 if it is found (a warning is printed otherwise) |
| `CMAKE_CXX_STANDARD` | `14` | A higher standard can be requested |

Compile-time switches in `globals.h`:

| Macro | Default | Meaning |
|---|---|---|
| `LES` | `0` | Smagorinsky LES. The Smagorinsky dynamics has no force term, so flow diverters have no effect with `LES 1`. |
| `DESCRIPTOR`, `BackgroundDynamics` | D3Q19, `GuoExternalForceCompleteRegularizedBGKdynamics` | Lattice and collision model |

---

## Setting up the preprocessor

The preprocessor needs Python 3.12 with NumPy, SciPy, VTK, `numpy-stl` and
`pynrrd`. The environment is in `preprocessor/environment.yml`:

```bash
conda env create -f preprocessor/environment.yml    # creates the env "lbmpre"
conda activate lbmpre
```

The verification scripts in `tests/` also need `pyvista`, `pandas` and
`matplotlib` (and `easyvvuq`, `dask` for the campaigns). Keep them in a
separate environment:

```bash
conda create -n hemoflow-post -c conda-forge python=3.12 pyvista pandas matplotlib
```

---

## Preparing a geometry (preprocessor)

```bash
conda run -n lbmpre python preprocessor/main.py path/to/case_vox.config
```

All file names in the config are relative to the config file's directory, and
all outputs are written there.

### Input requirements

* **Vessel surface** (`.stl`): watertight triangulated surface in
  **millimetres** (the factor `SI_FACTOR = 0.001` is hard-coded in `main.py`).
  The vessel must be clipped so that every inlet/outlet is a flat open cut
  that lies on a face of the axis-aligned bounding box, and **not on an edge
  or corner** of the bounding box. Triangles should not be much smaller than
  a voxel.
* **Centerline** (`.vtp`), e.g. from `vmtkcenterlines`: one polyline per
  outlet, each starting at the inlet and ending at that outlet, carrying the
  point array `MaximumInscribedSphereRadius`. The start point defines the
  inlet; every opening gets its radius and direction from the centerline.
* **Flow diverter** (optional): a surface mesh of the deployed device in the
  same coordinate frame as the vessel, as `<stent_mesh_base>mesh.stl`. It is
  voxelized from three projections; the voxels it touches inside the lumen
  become porous voxels.

### Configuration file (JSON)

Example (`tests/verification/0-pipe/input/input_pipe_vox.config`):

```json
{
"geometry_original_stl" : "pipe_Srf_prep.stl",
"centerline_vtp" : "pipe_Ctl_prep.vtp",
"stent_mesh_base" : "stent_",
"target_elements" : "1000000",
"target_dx" : "0.3",
"output_base_name" : "vox_pipe_"
}
```

| Key | Required | Meaning |
|---|---|---|
| `geometry_original_stl` | yes | Vessel surface STL (mm) |
| `centerline_vtp` | yes | Centerline VTP |
| `stent_mesh_base` | yes (may be `""`) | Prefix of the flow diverter files: `<prefix>mesh.stl` is the geometry. `""` = no device. A configured but missing file is an error. |
| `stent_folder` | no | Alternative to `stent_mesh_base`: uses `<stent_folder>/<stent_folder>_<config dir name>_stent_mesh.stl` |
| `target_elements` | yes | Target number of voxels in the bounding box. Sets the resolution when `target_dx` is absent. |
| `target_dx` | no | Target voxel size in **mm**. Overrides `target_elements`. The final `dx` is adjusted slightly so that the domain spans a whole number of voxels; the exact value is printed and stored in the output. |
| `output_base_name` | yes | Prefix for all output files |
| `cutWidth`, `distance` | – | Legacy keys, ignored (still present in several shipped configs) |

Further options are constants at the top of `preprocessor/main.py`:
`SI_FACTOR`, `DISTANCE_RADIUS_FACTOR` (opening matching tolerance relative to
the opening radius), `DEBUG_MODE` (writes intermediate `.nrrd` files) and
`INHOMOGEN` (reads per-point porous coefficients from `<stent_mesh_base>values.vtp`).

### Outputs

* `<output_base_name>c.npz`: simulation input for `hemoFlow`
* `<output_base_name>geometry.nrrd`: label volume, to check the voxelization
  (e.g. in 3D Slicer)

Contents of the `.npz` archive (all arrays C-ordered):

| Array | Type, shape | Meaning |
|---|---|---|
| `geometryFlag` | int16, (Nx, Ny, Nz) | Voxel labels (table below) |
| `dx` | float64, (1,) | Voxel size [m] |
| `origin` | float64, (3,) | Position of voxel (0, 0, 0) in the frame of the input STL [m] |
| `openingIndex` | int16, (N,) | Label of each opening |
| `openingRadius` | float64, (N,) | Centerline radius at each opening [m] |
| `openingNormal` | float64, (N, 3) | Unit direction of each opening, pointing **into** the domain |
| `openingCenter` | float64, (N, 3) | Centre of each opening [voxel indices] |
| `openingNormalizedQRatio` | float64, (N,) | r³ / Σ r³ over the outlets (diagnostics) |
| `stent` | int16, (Nx, Ny, Nz) | 1 in the porous (flow diverter) voxels, only if a device was given |
| `linear`, `quadratic` | float64, (Nx, Ny, Nz) | Per-voxel porous coefficients (`INHOMOGEN`); not read by the solver yet |

### Voxel labels

| Label | Meaning | Treatment in the solver |
|---|---|---|
| 0 | Outside the vessel | No dynamics |
| 1 | Wall | Bounce-back (no-slip) |
| 2 | Fluid | Fluid |
| 10 | **Inlet**: the opening at the start of the centerline | Boundary condition from config |
| 11 | The **smallest** outlet (by centerline radius) | Boundary condition from config |
| 12, 13, … | The other outlets, by decreasing radius | Boundary condition from config |
| 20 | Fluid voxels inside the flow diverter (for visualization) | Fluid; the porous forcing comes from the `stent` array |

`openingIndex` lists the labels in the order 10, 12, 13, …, 11. The solver's
start-up log prints the face, centre, direction and size of every opening,
which helps to identify the labels.

---

## Running a simulation (solver)

```bash
mpirun -np 8 ./build/hemoFlow path/to/case.xml [-d] [-c]
```

| Flag | Meaning |
|---|---|
| `-d`, `--debug` | Every 200 steps, print the imposed and the measured flow rate of every opening and the net mass balance |
| `-c`, `--checkpoint` | Restart from the checkpoint in `<outputDir>` (use the same number of MPI processes) |
| `-h`, `--help` | Usage |

Paths inside the XML file (geometry, waveforms, output directory) are
relative to the directory of the XML file, or absolute. Configuration errors
(missing or unknown tags, labels that are not in the geometry, unreadable
waveforms, …) stop the solver with a message before the simulation starts.

### Physical model and units

* All solver inputs are in SI units (m, s, Pa). (The preprocessor's
  `target_dx` and the STL are in mm.)
* The fluid is Newtonian. The defaults are blood, ν = 3.22·10⁻⁶ m²/s and
  ρ = 1055 kg/m³, and can be changed in the `<fluid>` section. ρ only scales
  the pressure and stress values.
* If `<dt>` is given, the relaxation rate follows from dx, dt and ν. If it is
  omitted, ω = 1 and dt = dx²/(6ν). The solver warns when ω > 1.9 or when
  the peak inlet velocity exceeds 0.1 in lattice units. Both warnings mean
  dt should be reduced.

### Configuration file (XML)

An annotated template is in `template_config.xml` (same as
`testcase/template_config.xml`). The top-level sections are `<simulation>`,
`<fluid>`, `<geometry>` and `<flowdiverter>`.

**`<simulation>`**

| Tag | Required | Meaning |
|---|---|---|
| `dt` | no | Time step [s] (see above) |
| `cycle` | yes | Number of cardiac cycles to simulate. Output is written only during the last cycle. |
| `saveFrequency` | yes | Output interval [s], ≥ `dt` |
| `outputDir` | yes | Output directory, created if missing |
| `outputFields` | no | Space-separated list of `velocity`, `pressure`, `stress`. Default: `velocity pressure`. |
| `checkpointFrequency` | no | Checkpoint interval [s]. Without it, checkpointing is disabled. |
| `mode` | no | Boundary-condition mode: none (default), `murray` or `pzero`, see below |
| `sparse` | no | `true`: only allocate the blocks that contain the vessel. Default `false`; see [Known issues](#known-issues-and-limitations). |
| `blockSize` | no | Block size of the sparse decomposition [lattice nodes], default 50 |

The **simulated time is `cycle` × period**, where the period is the last time
entry of the first waveform file found among the openings. At least one
opening must have a `timeScaleFunction`; for steady inflow use a two-point
constant file, whose last time sets the length of a "cycle".

**`<fluid>`** (optional): `kinematicViscosity` [m²/s], `density` [kg/m³].

**`<geometry>`**

* `file`: the `.npz` file from the preprocessor.
* `opening_0`, `opening_1`, …: one block per opening. The number in the tag
  is a running index from 0 to N−1 (N = number of openings in the geometry
  file) and is unrelated to the label:

| Tag | Meaning |
|---|---|
| `name` | Free text, used in the log |
| `type` | 1 = velocity, 2 = Murray outlet, 3 = pressure, 4 = free outflow |
| `label` | Voxel label of the opening (10, 11, …) |
| `timeScaleFunction` | Waveform file (optional) |
| `parameter` | Velocity openings (required): mean velocity [m/s]. Pressure openings: gauge pressure [Pa] (default 0). Ignored for Murray and free-outflow openings. |

**Opening types**

| Type | Boundary condition |
|---|---|
| 1 velocity | Parabolic (Poiseuille) profile on the opening, oriented along the centerline direction, zero at the fitted ellipse, with peak velocity 2 × `parameter`, multiplied by the waveform. A negative value gives an outflow. Because the opening is discrete, the actual mean velocity differs slightly from `parameter`; the log prints the actual flow rate and mean velocity. |
| 2 Murray | Velocity outlet whose flow rate is set every time step by Murray's law, Q_i = Q_in · d_i³ / Σ d_j³ (Q_in: net flow of all velocity openings; Σ over all other openings), [Chnafa et al., AJNR 2018](https://doi.org/10.3174/ajnr.A5484). The imposed flux of the discrete profile matches Q_i exactly. |
| 3 pressure | Gauge pressure `parameter` [Pa], multiplied by the waveform if one is given |
| 4 free outflow | Palabos `VirtualOutlet` |

**Boundary-condition modes** (`<simulation><mode>`)

| Mode | Behaviour |
|---|---|
| *(absent)* | Every opening must be described by an `opening_i` block |
| `murray` | Label 10 (inlet) becomes a velocity inlet with `parameter` and `timeScaleFunction` from `opening_0`. Label 11 (smallest outlet) becomes a 0 Pa pressure outlet, all other openings Murray outlets. Only `opening_0` is read. |
| `pzero` | Label 10 becomes a velocity inlet (from `opening_0`), all other openings 0 Pa pressure outlets |

**`<flowdiverter>`** (required if the geometry contains a `stent` array)

| Tag | Meaning |
|---|---|
| `linCoeff` | Linear (Darcy) coefficient *a* [Pa·s/m] |
| `quadCoeff` | Quadratic (Forchheimer) coefficient *b* [Pa·s²/m²] |

In every porous voxel the force density −(*a*·**u** + *b*·**u**|**u**|)/dx is
applied. The coefficients therefore describe the pressure drop across **one
voxel-thick layer**, Δp = *a*·u + *b*·u²; the total resistance of a device
scales with its thickness in voxels.

### Waveform (time-scale function) files

Plain text: the first line is the number of samples *N*, followed by *N*
lines `time[s] value`. Times start at 0 and increase; the last time is the
period. The value multiplies `parameter`, values are linearly interpolated,
and the curve repeats every period, so the first and last values should
match (the solver warns otherwise). Example
(`testcase/input/cerebral_flowrate.txt`, normalised ICA waveform with a
0.8 s period):

```
124
0	0.2161882
0.004713109	0.217917705
...
0.8	0.2161882
```

### Output

* `<outputDir>/vtkNNNNNN.vti` (NNNNNN = iteration): VTK image data, written
  every `saveFrequency` during the last cycle, with the fields selected in
  `outputFields`:
  * `velocity [m/s]`
  * `pressure [Pa]`: gauge pressure, 0 at the reference density
  * `stress [Pa]`: viscous stress tensor (components xx, xy, xz, yy, yz, zz),
    for wall shear stress; see `doc/paraview_wss_calc.txt`
* The grid has spacing dx and the `origin` of the geometry file, so it
  overlays the input STL once the STL is scaled from mm to m. (Geometry files
  from older preprocessor versions have no origin; the output then starts at
  (0, 0, 0).)
* Checkpoints (if enabled): `checkpoint_lattice.dat` and
  `checkpoint_parameters.dat`, with the previous checkpoint kept as `*_old.dat`.
* The log (stdout) contains the derived lattice parameters and, for every
  opening, its label, face, radius, centre, direction, size and imposed flow
  rate.

---

## Test cases and verification

All examples are run from inside their directory, so the relative paths in
the configs resolve. Voxelized geometries (`*.npz`, `*.nrrd`) and `output*/`
directories are git-ignored. Timings below are from a 24-core workstation.

### 1. Straight pipe, Poiseuille flow (`tests/verification/0-pipe`), quick check

220 k voxels (434 × 23 × 22, dx ≈ 0.3 mm). The pipe has a 6.5 mm diameter and
a 0.45 m/s mean inflow, with a 0 Pa outlet and a flow diverter shell with zero
resistance. The run covers 1.6 s (80 k steps).

```bash
cd tests/verification/0-pipe/input
conda run -n lbmpre python ../../../../preprocessor/main.py ./input_pipe_vox.config   # ~2 s
cd ..
mpirun -np 8 ../../../build/hemoFlow pipe.xml                                        # ~7 min on 10 cores
conda run -n hemoflow-post python pressure_drop.py .
```

`pressure_drop.py` compares the pressure drop between ¼ and ¾ of the pipe
length with Hagen–Poiseuille for the simulated lumen, and writes
`output.csv` and `profile_comp.png`.

**Reference result** (November 2026):

| Quantity | Value |
|---|---|
| Inlet / outlet flow rate | 1.537·10⁻⁵ / 1.539·10⁻⁵ m³/s (0.1 % apart) |
| Pressure drop over the central 65 mm | 80.6 Pa |
| Hagen–Poiseuille for the simulated lumen (6.50 mm) | 77.7 Pa (+3.8 %) |

At Re ≈ 900 the laminar entrance length (≈ 0.05·Re·D ≈ 290 mm) is longer
than the pipe (130 mm), so a pressure drop slightly above Poiseuille is
expected. The earlier V&V result (`doc/VV`, 90.7 Pa at 0.5 m/s and 2.5 s)
scales to 81.6 Pa at 0.45 m/s; the same 2.5 s run gives 81.4 Pa with the
current code.

### 2. Aneurysm with flow diverter (`testcase/`)

NAP180 sidewall aneurysm with a PED 5×20 flow diverter (NAP126 is also
included). Three voxelization configs are provided in `testcase/input/`:

| Config | Elements | Output | Used by |
|---|---|---|---|
| `voxel_input.config` | 1 M (80 × 69 × 181, dx ≈ 0.36 mm, ~15 s) | `vox_NAP180_ane_c.npz` | `test_Aneurysm.xml` |
| `voxel_input_5M.config` | 5 M | `vox_NAP180_ane_PED_5x20_5M_c.npz` | `test_Velocity/Pressure/Murray/Freeflow.xml` |
| `voxel_input_20M.config` | 20 M | `vox_NAP180_ane_PED_5x20_20M_c.npz` | – |

```bash
cd testcase/input
conda run -n lbmpre python ../../preprocessor/main.py ./voxel_input.config
cd ..
mpirun -np 8 ../build/hemoFlow test_Aneurysm.xml
```

`test_Aneurysm.xml` uses the `murray` mode with steady inflow (0.8 s,
dt = 10⁻⁵ s, 80 k steps on 1 M cells, about 25 min on 22 cores). The other
`test_*.xml` files show the different opening types with the ICA waveform.

### 3. Flow split / Murray's law (`tests/verification/1-flow_split`)

Bifurcation with three outlets, `murray` mode, two cycles of the ICA
waveform.

```bash
cd tests/verification/1-flow_split/input
conda run -n lbmpre python ../../../../preprocessor/main.py ./input_branch_vox.config   # 5 M elements, ~20 s
cd ..
mpirun -np 20 ../../../build/hemoFlow flow_split.xml      # 80 k steps on 5 M cells: ~1-2 h on 20 cores
conda run -n hemoflow-post python Murray_check.py .
```

`Murray_check.py` measures the flow through each branch 8 voxels inside the
domain and compares the outlets with their Murray share. The solver imposes
the Murray share exactly on the Murray outlets (check with `-d`). In a 0.1 s
run with steady inflow (November 2026), the interior measurement gave
+0.3 % and 0.0 % for the two Murray outlets, +9 % for the still-settling
0 Pa outlet, and inflow and outflow within 1 %.

### 4. Aneurisk C0002 (`tests/verification/2-aneurisk_C0002`)

Sidewall aneurysm from the
[Aneurisk](https://ecm2.mathcs.emory.edu/aneuriskweb/about) database with a
normalised ICA waveform and `murray` outlets. Voxelize with
`input/voxelize.sh` (~5 s), then run `aneurisk.xml`.

### 5. Sensitivity and uncertainty campaigns (`tests/sensitivity`, `tests/uncertainty`)

EasyVVUQ campaigns that voxelize, run and post-process parameter sweeps over
dx, dt, save interval and inflow velocity. Before running them, edit the
hard-coded `SOFTWARE_PATH` at the top of each `*_campaign.py` /
`*_test.py`. The aneurysm campaigns call LBMpost, an in-house tool that is not
public. Results are summarised in `doc/VV/hemoflow_VV.ipynb`.

---

## Known issues and limitations

The full list of open tasks is in [TODO.md](TODO.md).

* **Sparse decomposition is experimental** (`<sparse>true</sparse>`). It
  gave results bit-identical to the dense decomposition for Aneurisk C0002
  (block size 16, 29 % of the cells, 1.8× faster) and for NAP180 at block
  size 24. For NAP180 at block size 16, the two runs diverged from a
  round-off-level difference next to the 0 Pa outlet, reaching 0.06 m/s
  after 400 steps. Until this is understood, compare against a dense run
  before relying on it.
* Openings must lie on bounding-box faces, away from edges and corners.
  Boundary conditions are instantiated on the bounding rectangle of each
  opening, so the solver refuses two openings on the same face whose
  rectangles overlap.
* The velocity `parameter` is the mean velocity of an ideal paraboloid
  (peak = 2 × `parameter`). On the discrete opening the actual mean differs
  slightly (0.464 instead of 0.45 m/s in the pipe case); the log reports it.
* The porous coefficients act per voxel layer, so the resistance of a device
  depends on its voxelized thickness and therefore on dx. The per-voxel
  coefficients (`linear`, `quadratic`) are written but not used by the solver.
* The scan-line voxelizer drops up to one partial voxel layer at the maximum
  y and z sides of the vessel, is slow for large targets (tens of millions of
  elements) and prints "is the geometry watertight?" messages for scan lines
  through the open cuts.
* The whole geometry file is loaded on every MPI process, which limits the
  domain size per node.
* The fluid is Newtonian; non-Newtonian rheology is not implemented.

---

## Troubleshooting

**Parallel HDF5 from source.** Only needed for the (currently unused) HDF5
writer:

```bash
git clone https://github.com/HDFGroup/hdf5 && cd hdf5 && git checkout hdf5-1_14_0
cmake -S . -B build -DHDF5_ENABLE_PARALLEL=ON -DHDF5_ENABLE_Z_LIB_SUPPORT=ON -DCMAKE_INSTALL_PREFIX=$HOME/opt/hdf5
cmake --build build -j 8 --target install
```

Then configure HemoFlow with `-DHDF5_ROOT=$HOME/opt/hdf5`.

**`No opening has a <timeScaleFunction>`.** Give at least one opening
(normally the inlet) a waveform file; use a two-point constant file for
steady inflow.

**Walls appear inside the fluid.** The STL triangles are too small compared
to the voxel size, or the surface is not closed. Increase the resolution or
remesh the surface.

**Wrong number of openings ("number of outlets found ... differ").** An
opening lies on an edge or corner of the bounding box, or the centerline does
not end at every opening. Re-clip the geometry so each cut lies cleanly on one
face.

**The simulation diverges (NaN energy).** The solver writes the current state
and exits with code 1. Reduce `dt` (check the ω and velocity warnings at
start-up) or increase the resolution.

---

## License

HemoFlow is released under the GNU Affero General Public License v3.0.
Palabos is distributed under its own license (AGPL v3), see `palabos/COPYING`.
