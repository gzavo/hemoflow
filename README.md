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
| `io.h/.cpp` | VTK, HDF5 and NPZ output and directory helpers |
| `globals.h`, `helper.h` | Lattice/dynamics typedefs, compile-time switches, geometry labels, small vector types |
| `CMakeLists.txt` | Build script (builds Palabos as a static library and links `hemoFlow`) |
| `template_config.xml` | Annotated solver configuration template. Prefer `testcase/template_config.xml`: this copy skips `opening_2`, which the solver rejects for three openings. |
| `preprocessor/` | Python voxelization and opening-detection tools (`main.py` is the entry point) |
| `testcase/` | Example aneurysm geometries (NAP126, NAP180) with a PED 5×20 flow diverter and example configs |
| `tests/verification/` | Verification cases: straight pipe (Poiseuille), bifurcation flow split (Murray), Aneurisk C0002 |
| `tests/sensitivity/`, `tests/uncertainty/` | [EasyVVUQ](https://github.com/UCL-CCS/EasyVVUQ) sensitivity and uncertainty campaigns |
| `doc/VV/` | Verification and validation notebook (`hemoflow_VV.ipynb`) and its data |
| `doc/paraview_wss_calc.txt` | Notes on computing wall shear stress in ParaView |
| `palabos/` | Vendored Palabos library (external dependency) |
| `external/` | Vendored [HighFive](https://github.com/BlueBrain/HighFive) (HDF5 C++ wrapper) and [cnpy](https://github.com/rogersce/cnpy) (`.npz` reader) |
| `scripts/`, `utils/` | Cluster job/environment scripts and small helpers (output compression) |

---

## Building the solver

### Requirements

| Dependency | Notes |
|---|---|
| C++ compiler | GCC or Clang. GCC ≥ 15 needs a workaround flag, see below. |
| CMake | ≥ 3.5. CMake ≥ 4 needs `-DCMAKE_POLICY_VERSION_MINIMUM=3.5`, see below. |
| MPI | e.g. OpenMPI or MPICH. |
| **Parallel** HDF5 | Required at compile time: `io.h` includes HDF5/HighFive unconditionally, even though the HDF5 writer is currently not called. `-DBUILD_HDF5=OFF` does **not** produce a working build. |
| zlib | |

Palabos, HighFive and cnpy come with the repository; nothing else needs
to be downloaded.

### Option A: conda toolchain (tested, no root access needed)

This is the route used to verify the instructions in this README (conda-forge,
GCC 15.3, CMake 4.4, OpenMPI 5.0, HDF5 2.2):

```bash
conda create -n hemoflow-build -c conda-forge --override-channels \
    cxx-compiler cmake make openmpi "hdf5=*=mpi_openmpi*" zlib
conda activate hemoflow-build
```

From the repository root:

```bash
mkdir -p build
echo 'add_compile_options(-Wno-error=template-body)' > build/gcc15-workaround.cmake
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_POLICY_VERSION_MINIMUM=3.5 \
      -DCMAKE_PROJECT_INCLUDE=$PWD/build/gcc15-workaround.cmake
cmake --build build -j 8
```

The executable is `build/hemoFlow`. The first build compiles all of Palabos
and takes a few minutes; later builds only recompile the HemoFlow sources.

> **Why the two extra flags?**
> * `CMAKE_POLICY_VERSION_MINIMUM=3.5`: `CMakeLists.txt` declares
>   `cmake_minimum_required(VERSION 3.0)`, which CMake ≥ 4 rejects.
> * `-Wno-error=template-body`: GCC ≥ 15 rejects a template in the vendored
>   Palabos (`multiGrid/multiGridParameterManager.hh`) that is never used. `CMakeLists.txt`
>   overwrites `CMAKE_CXX_FLAGS`, so the flag cannot go through
>   `-DCMAKE_CXX_FLAGS=...`; the project-include file adds it instead.
>
> Both flags are harmless with older tools, so they can always be passed.

### Option B: system packages (Debian/Ubuntu)

```bash
sudo apt install build-essential cmake openmpi-bin libopenmpi-dev libhdf5-openmpi-dev zlib1g-dev
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j 8
```

With the Ubuntu 24.04 packages (GCC 13, CMake 3.28) neither workaround flag
is needed: Palabos compiles cleanly with GCC 13. The HDF5 1.10 package has not
been verified end-to-end. If CMake picks up a serial HDF5 or HighFive fails to
compile, build parallel HDF5 yourself (see [Troubleshooting](#troubleshooting)).

### Option C: HPC clusters

Load a compiler, CMake, MPI and a parallel HDF5 module, then build as in
Option B. `scripts/snellius_env.sh` lists the modules used on SURF Snellius
in 2021 (`source scripts/snellius_env.sh`); module versions on current systems
will differ.

### Build options

| Option | Default | Meaning |
|---|---|---|
| `CMAKE_BUILD_TYPE` | `Release` | `Release` (`-O3`) or `Debug` (`-g -O0 -DPLB_DEBUG`, executable `hemoFlowd`) |
| `ENABLE_MPI` | `ON` | Build with MPI |
| `BUILD_HDF5` | `ON` | Searches for parallel HDF5 (see the note in the requirements above) |

Some switches are compile-time macros in `globals.h`:

| Macro | Default | Meaning |
|---|---|---|
| `SPARSE` | `0` | `1` allocates only blocks that contain fluid (`<blockSize>` sets the block size). Strongly recommended for aneurysm geometries, where usually < 10 % of the bounding box is fluid. |
| `LES` | `0` | Smagorinsky LES. **Currently does not compile** (`cSmago` is not defined). |
| `DESCRIPTOR`, `BackgroundDynamics` | D3Q19, `GuoExternalForceCompleteRegularizedBGKdynamics` | Lattice and collision model |

---

## Setting up the preprocessor

The preprocessor needs Python 3.9 with VTK 9.1, NumPy, SciPy, `numpy-stl` and
`pynrrd`. The pinned environment is in `preprocessor/environment.yml`:

```bash
conda env create -f preprocessor/environment.yml    # creates the env "lbmpre"
conda activate lbmpre
```

The verification and post-processing scripts in `tests/` also need `pyvista`,
`pandas` and `matplotlib` (and `easyvvuq`, `dask` for the campaigns). Keep them
in a separate environment so the pinned preprocessor env stays unchanged:

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
  point array `MaximumInscribedSphereRadius`. Openings get their radius and
  direction from the centerline end points.
* **Flow diverter** (optional): a surface mesh of the deployed device in the
  same coordinate frame as the vessel, as `<stent_mesh_base>mesh.stl`. It is
  voxelized from three projections, and every voxel it touches becomes a porous voxel.

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
| `stent_mesh_base` | yes (may be `""`) | Prefix of the flow diverter files: `<prefix>mesh.stl` is the geometry. `""` = no device. |
| `stent_folder` | no | Alternative to `stent_mesh_base`: uses `<stent_folder>/<stent_folder>_<config dir name>_stent_mesh.stl` |
| `target_elements` | yes | Target number of voxels in the bounding box. Sets the resolution when `target_dx` is absent. |
| `target_dx` | no | Target voxel size in **mm**. Overrides `target_elements`. The final `dx` is adjusted slightly so that the domain spans a whole number of voxels; the exact value is printed and stored in the output. |
| `output_base_name` | yes | Prefix for all output files |
| `cutWidth`, `distance` | – | Legacy keys, **ignored** (still present in several shipped configs) |

Further options are hard-coded constants at the top of `preprocessor/main.py`:
`SI_FACTOR`, `DISTANCE_RADIUS_FACTOR` (opening detection tolerance relative to
the opening radius), `DEBUG_MODE` (writes intermediate `.nrrd` files) and
`INHOMOGEN` (reads per-point porous coefficients from `<stent_mesh_base>values.vtp`).

### Outputs

* `<output_base_name>c.npz`: simulation input for `hemoFlow`
* `<output_base_name>geometry.nrrd`: label volume, to check the voxelization
  (e.g. in 3D Slicer)

Contents of the `.npz` archive:

| Array | Type, shape | Meaning |
|---|---|---|
| `geometryFlag` | int16, (Nx, Ny, Nz) | Voxel labels (table below) |
| `dx` | float64, (1,) | Voxel size [m] |
| `openingIndex` | int16, (N,) | Label of each opening |
| `openingRadius` | float64, (N,) | Centerline radius at each opening [m] |
| `openingNormal` | float64, (N, 3) | Unit direction of each opening, pointing **into** the domain |
| `openingCenter`, `openingNormalizedQRatio` | float64 | Diagnostics only, not used by the solver |
| `stent` | int16, (Nx, Ny, Nz) | 1 inside the flow diverter (only if a device was given) |
| `linear`, `quadratic` | int32 | Per-voxel porous coefficients. Not used by the solver; see [Known issues](#known-issues-and-limitations). |

### Voxel labels

| Label | Meaning | Treatment in the solver |
|---|---|---|
| 0 | Outside the vessel (unused) | No dynamics |
| 1 | Wall | Bounce-back (no-slip) |
| 2 | Fluid | Fluid |
| 10 | Opening with the **largest** centerline radius (normally the inlet) | Boundary condition from config |
| 11 | Opening with the **smallest** centerline radius | Boundary condition from config |
| 12, 13, … | Remaining openings, by decreasing radius | Boundary condition from config |
| 20 | Flow diverter voxels (added for visualization) | Fluid. The porous forcing comes from the `stent` array. |

`openingIndex` lists the labels in the order 10, 12, 13, …, 11. The
`murray` and `pzero` modes of the solver rely on this order (see below).
Check the label of your inlet in the `.nrrd` file, or in the solver's
start-up log, which prints the centre and normal of every opening.

---

## Running a simulation (solver)

```bash
mpirun -np 8 ./build/hemoFlow path/to/case.xml [-d] [-c]
```

| Flag | Meaning |
|---|---|
| `-d`, `--debug` | Every 200 steps, also print the flow rate of every opening |
| `-c`, `--checkpoint` | Restart from the checkpoint in `<outputDir>` |
| `-h`, `--help` | Usage |

All paths inside the XML file (geometry, waveforms, output directory) are
**relative to the directory of the XML file**. Absolute paths are not
supported.

### Physical model and units

* All solver inputs are in SI units (m, s, Pa). (The preprocessor's
  `target_dx` and the STL are in mm.)
* The blood is Newtonian with ν = 3.22·10⁻⁶ m²/s, hard-coded in
  `hemoFlow.cpp`. ρ = 1055 kg/m³ (`globals.h`) is only used to convert
  pressures. Carreau parameters are set in the code but not used by the
  current collision model.
* If `<dt>` is given, the relaxation rate follows from dx, dt and ν. If it is
  omitted, ω = 1 and dt = dx²/(6ν). The solver warns when ω > 1.9 or when
  the peak inlet velocity exceeds 0.1 in lattice units. Both warnings mean
  dt should be reduced.

### Configuration file (XML)

An annotated template is in `testcase/template_config.xml`. The top-level
sections are `<simulation>`, `<geometry>` and `<flowdiverter>`.

**`<simulation>`**

| Tag | Required | Meaning |
|---|---|---|
| `dt` | no | Time step [s] (see above) |
| `cycle` | yes | Number of cardiac cycles to simulate. Output is written only during the last cycle. |
| `saveFrequency` | yes | Output interval [s]. Must be ≥ `dt`. |
| `outputDir` | yes | Output directory, created if missing |
| `blockSize` | yes | Atomic block size for the sparse decomposition (only used with `SPARSE 1`) |
| `checkpointFrequency` | no | Checkpoint interval [s]. Without it, checkpointing is disabled. |
| `mode` | no | Boundary-condition mode: none (default), `murray` or `pzero`, see below |

The **simulated time is `cycle` × period**, where the period is the last time
entry of the first waveform file found among the openings. At least one
opening must have a `timeScaleFunction`, otherwise the solver stops.

**`<geometry>`**

* `file`: the `.npz` file from the preprocessor.
* `opening_0`, `opening_1`, …: one block per opening. The number in the tag
  is just a running index from 0 to N−1 (N = number of openings in the
  geometry file) and is unrelated to the label:

| Tag | Meaning |
|---|---|
| `name` | Free text, used in the log |
| `type` | 1 = velocity, 2 = Murray outlet, 3 = pressure, 4 = free outflow |
| `label` | Voxel label of the opening (10, 11, …) |
| `timeScaleFunction` | Waveform file (may be empty) |
| `parameter` | Velocity openings: **mean** velocity [m/s]. Pressure openings: pressure [Pa]. Ignored for Murray and free-flow openings. |

**Opening types**

| Type | Boundary condition |
|---|---|
| 1 velocity | Parabolic (Poiseuille) profile, oriented along the centerline direction, with peak velocity 2 × `parameter`, multiplied by the waveform. A negative value gives an outflow. |
| 2 Murray | Velocity outlet whose flow rate is redistributed every step following Murray's law, Q_i = Q_in · d_i³ / Σ d_j³ (Σ over all non-velocity openings), [Chnafa et al., AJNR 2018](https://doi.org/10.3174/ajnr.A5484) |
| 3 pressure | Constant pressure (gauge, relative to the reference density). Use 0 Pa; see [Known issues](#known-issues-and-limitations) for non-zero or time-dependent values. |
| 4 free outflow | Palabos `VirtualOutlet` |

**Boundary-condition modes** (`<simulation><mode>`)

| Mode | Behaviour |
|---|---|
| *(absent)* | Every opening must be described by an `opening_i` block |
| `murray` | The first opening in `openingIndex` (label 10) becomes a velocity inlet using `parameter` and `timeScaleFunction` from `opening_0`. The last one (label 11, smallest) becomes a 0 Pa pressure outlet, and all others become Murray outlets. Only `opening_0` is needed. |
| `pzero` | Label 10 becomes a velocity inlet (parameters from `opening_0`), and all other openings become 0 Pa pressure outlets |

Any other value of `<mode>` is silently treated as "absent".

**`<flowdiverter>`** (required only if the geometry contains a `stent` array)

| Tag | Meaning |
|---|---|
| `linCoeff` | Linear (Darcy) coefficient *a* |
| `quadCoeff` | Quadratic (Forchheimer) coefficient *b* |

In every porous voxel the force −(*a*·**u** + *b*·**u**|**u**|)/dx is
applied. The coefficients therefore describe the pressure drop across **one
voxel-thick layer**, Δp = *a*·u + *b*·u² (*a* in Pa·s/m, *b* in
Pa·s²/m²). The total resistance of the device scales with its thickness in
voxels.

### Waveform (time-scale function) files

Plain text: the first line is the number of samples *N*, followed by *N*
lines `time[s] value`. The value multiplies `parameter`, values are linearly
interpolated, and the curve repeats with a period equal to the last time
entry, so the first and last values should match. Example
(`testcase/input/cerebral_flowrate.txt`, normalised ICA waveform with a
0.8 s period):

```
183
0	0.2161882
0.004713109	0.217917705
...
0.8	0.2161882
```

### Output

* `<outputDir>/vtkNNNNNN.vti` (NNNNNN = iteration): VTK image data with the
  field `velocity [m/s]`, written every `saveFrequency` during the last
  cycle. The grid spacing is dx and the origin is (0, 0, 0). This is the
  voxel frame of the preprocessed domain, **not** the coordinate frame of the
  input STL.
* Checkpoints (if enabled): `checkpoint_lattice.dat` and
  `checkpoint_parameters.dat`, with the previous checkpoint kept as `*_old.dat`.
* The log (stdout) contains the derived lattice parameters and, for every
  opening, its label, radius, centre, normal, area and imposed flow rate.

Pressure and stress fields are currently not written; see
[Known issues](#known-issues-and-limitations).

---

## Test cases and verification

All examples are run from inside their directory, so the relative paths in
the configs resolve. Voxelized geometries (`*.npz`, `*.nrrd`) and `output*/`
directories are git-ignored.

### 1. Straight pipe, Poiseuille flow (`tests/verification/0-pipe`), quick check

220 k voxels (434 × 23 × 22, dx ≈ 0.3 mm). The pipe has a 6.5 mm diameter and
a 0.45 m/s mean inflow, with a 0 Pa outlet and a flow diverter shell with zero
resistance.

```bash
cd tests/verification/0-pipe/input
conda run -n lbmpre python ../../../../preprocessor/main.py ./input_pipe_vox.config   # ~2 s
cd ..
mpirun -np 8 ../../../build/hemoFlow pipe.xml
conda run -n hemoflow-post python pressure_drop.py .
```

`pressure_drop.py` compares the pressure drop between ¼ and ¾ of the pipe
length with Hagen–Poiseuille and writes `output.csv` and `profile_comp.png`.

Before running it as shipped, note:

* `input/stationary.txt` has a 100 s period, so `<cycle>2</cycle>` means
  200 s of physical time (10 M steps). For a quick check, use a waveform
  file with a shorter period, e.g.
  `printf "2\n0.0 1.0\n1.25 1.0\n" > input/stationary_short.txt`, and point
  `timeScaleFunction` at it. That gives 2.5 s and 125 k steps, about 9 min
  on 11 cores.
* `pressure_drop.py` needs the `density [Pa]` field, which `writeVTK` in
  `io.cpp` currently does not write (the line is commented out).

**Reference result** (November 2026, the 2.5 s variant above, with the density
output line re-enabled):

| Quantity | Value |
|---|---|
| Inlet / outlet flow rate | 1.536·10⁻⁵ / 1.538·10⁻⁵ m³/s (0.1 % apart) |
| Pressure drop over the central 65 mm | 81.3 Pa |
| Hagen–Poiseuille, lumen diameter 6.48 mm | 78.5 Pa (+3.5 %) |
| Earlier V&V result (`doc/VV`, 90.7 Pa at 0.5 m/s), scaled to 0.45 m/s | 81.6 Pa |

At Re ≈ 900 the laminar entrance length (≈ 0.05·Re·D ≈ 290 mm) is longer
than the pipe (130 mm), so a pressure drop slightly above Poiseuille is
expected. `pressure_drop.py` itself reports a 12 % difference because it uses
the full domain width including the wall voxels (6.61 mm) as the pipe
diameter.

### 2. Aneurysm with flow diverter (`testcase/`)

NAP180 sidewall aneurysm with a PED 5×20 flow diverter (NAP126 is also
included). `input/voxel_input.config` voxelizes at 1 M bounding-box elements
(80 × 69 × 181, dx ≈ 0.36 mm, about 30 s). `input_P64_4_vox.config` is a 20 M
element variant.

```bash
cd testcase/input
conda run -n lbmpre python ../../preprocessor/main.py ./voxel_input.config   # -> vox_NAP180_ane_c.npz
cd ..
mpirun -np 8 ../build/hemoFlow test_Aneurysm.xml
```

`test_Aneurysm.xml` also uses the 100 s `stationary.txt` waveform with
dt = 10⁻⁵ s, i.e. 10 M steps. With a short waveform file (as in case 1), 5 000
steps take about 1.5 min on 22 cores (1 M lattice cells, dense
decomposition).

The `test_*.xml` files show the different opening types (`Velocity`,
`Pressure`, `Murray`, `Freeflow`). Apart from `test_Aneurysm.xml`, they
expect `input/vox_NAP180_ane_PED_5x20_5M_c.npz`. To create it, set
`target_elements` to `5000000` and `output_base_name` to
`vox_NAP180_ane_PED_5x20_5M_` in a copy of the voxelization config.

### 3. Flow split / Murray's law (`tests/verification/1-flow_split`)

Bifurcation with three outlets. Use `input/voxelize.sh` (5 M elements),
then run `flow_split.xml`. `Murray_check.py` compares the outlet flow
rates with Murray's law. It loads `input/vox_murray_5M_c.npz`, while the
voxelization config writes `input/vox_branch_5M_c.npz`, so adjust one of the
two names.

### 4. Aneurisk C0002 (`tests/verification/2-aneurisk_C0002`)

Sidewall aneurysm from the
[Aneurisk](https://ecm2.mathcs.emory.edu/aneuriskweb/about) database, with a
normalised ICA waveform. Voxelize with
`input/voxelize.sh`, run `aneurisk.xml`.

> In cases 2–4, `<mode>aneurysm</mode>` is an outdated name and is treated as
> "no mode". For cases 3 and 4, which only define `opening_0`, change it to
> `<mode>murray</mode>`; otherwise the solver stops with
> `Element opening_1 not found`. The voxelization configs of cases 3 and 4 set
> `"stent_mesh_base": "stent_"`, but their `input/` folders have no
> `stent_mesh.stl`, so the preprocessor stops. Set the key to `""` to
> voxelize without a device.

### 5. Sensitivity and uncertainty campaigns (`tests/sensitivity`, `tests/uncertainty`)

EasyVVUQ campaigns that voxelize, run and post-process parameter sweeps over
dx, dt, save interval and inflow velocity. Before running them, edit the
hard-coded `SOFTWARE_PATH` at the top of each `*_campaign.py` /
`*_test.py`. The aneurysm campaigns call LBMpost, an in-house tool that is not
public. Results are summarised in `doc/VV/hemoflow_VV.ipynb`.

---

## Known issues and limitations

These came up while reviewing the code (November 2026). They are worth
knowing before trusting results:

* **Flow diverter voxels break the wall.** The preprocessor writes label 20
  over every flow diverter voxel in `geometryFlag`, including wall (1) and
  outside (0) voxels. The solver treats label 20 as fluid, which opens holes
  in the vessel wall where the device touches it. There are 131 such voxels in
  the pipe case and 259 in the NAP180 testcase. In the pipe test with zero
  porous resistance this did not change the pressure drop measurably (81.256
  vs 81.255 Pa without the device), but it changes the geometry and applies
  porous forcing outside the lumen.
* **Small backflow at the rim of velocity openings.** The parabolic profile is
  not clamped at zero, so nodes outside the fitted ellipse get a small reverse
  velocity (7 of 367 inlet nodes in the pipe case, down to −0.035 m/s).
* Murray outlets are scaled with the nominal area π·r² (with π = 3.14) and
  an ideal paraboloid, not with the discrete flux actually imposed.
* **Non-zero pressure boundaries are wrong.** A pressure *p* is converted to
  the lattice density 1 + p/(3·C_p) instead of 1 + 3p/C_p (9× too small). A
  waveform on a pressure opening scales the absolute density, not the
  gauge pressure. 0 Pa outlets, as used in all shipped cases, are not
  affected.
* Only velocity is written to the output, so pressure and wall shear
  stress cannot be post-processed.
* After a checkpoint restart, the waveform restarts at t = 0 while the
  iteration counter continues, which shifts the phase.
* The `linear`/`quadratic` coefficient arrays (for `INHOMOGEN`) come out one
  voxel shorter than the geometry when the Xmax face is cut. The solver
  does not read them yet.
* Openings must lie on bounding-box faces, away from edges and corners.
  Openings are handled through their bounding rectangles, so two openings on
  the same face with overlapping rectangles can interfere.
* The pure-Python voxelizer is slow for large targets (tens of millions of
  elements), and the whole geometry is loaded on every MPI rank.

---

## Troubleshooting

**Parallel HDF5 from source.** If no suitable parallel HDF5 is available:

```bash
git clone https://github.com/HDFGroup/hdf5 && cd hdf5 && git checkout hdf5-1_14_0
cmake -S . -B build -DHDF5_ENABLE_PARALLEL=ON -DHDF5_ENABLE_Z_LIB_SUPPORT=ON -DCMAKE_INSTALL_PREFIX=$HOME/opt/hdf5
cmake --build build -j 8 --target install
```

Then configure HemoFlow with `-DHDF5_ROOT=$HOME/opt/hdf5`.

**`Compatibility with CMake < 3.5 has been removed`.** Add
`-DCMAKE_POLICY_VERSION_MINIMUM=3.5`.

**`error: ... [-Wtemplate-body]` in Palabos.** Use the project-include
workaround from [Option A](#option-a-conda-toolchain-tested-no-root-access-needed).

**`No opening with a timeScaleFunction was found`.** Give at least one
opening (normally the inlet) a waveform file; use a two-point constant file
for steady inflow.

**Walls appear inside the fluid, or "is the geometry watertight?"
messages.** The STL triangles are too small compared to the voxel size, or
the surface is not closed. Increase the resolution or remesh the surface.

**Wrong number of openings ("number of outlets found ... differ").** An
opening lies on an edge or corner of the bounding box, or the centerline does not
end at every opening. Re-clip the geometry so each cut lies cleanly on one face.

---

## License

HemoFlow is released under the GNU Affero General Public License v3.0.
Palabos is distributed under its own license (AGPL v3), see `palabos/COPYING`.
