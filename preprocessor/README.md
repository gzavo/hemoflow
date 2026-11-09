# HemoFlow preprocessor

Turns a vessel surface (STL), its centerline (VTP) and, optionally, a deployed flow
diverter (STL) into the voxelized geometry file (`.npz`) that the `hemoFlow` solver
reads. The main [README](../README.md) documents the configuration keys, the input
requirements and the output arrays in detail.

## Installation

```bash
conda env create -f environment.yml    # creates the env "lbmpre" (Python 3.12)
```

## Usage

```bash
conda run -n lbmpre python main.py path/to/case_vox.config
```

All file names in the JSON config are relative to the config file, and all outputs are
written next to it:

* `<output_base_name>c.npz`: the solver input
* `<output_base_name>geometry.nrrd`: the label volume, to inspect the voxelization
  (e.g. in 3D Slicer)

## What it does

1. `voxelizeStl.py`, `slice.py`, `perimeter.py`: scan-line voxelization of the STL
   (pure Python).
2. `createFluidSolid.py`: adds a one-voxel wall around the fluid, removes empty outer
   layers and cuts one layer from every bounding-box face that contains an opening.
3. `readCL.py`: reads the opening radii and directions from the centerline end points.
4. `detectOpenings.py`: finds the openings on the domain faces and labels them.
5. `main.py`: matches openings to centerline end points, voxelizes the flow diverter
   (union of three projections, restricted to the lumen), computes the physical origin
   of the grid and writes the outputs.

## Voxel labels

| Label | Meaning |
|---|---|
| 0 | Outside the vessel |
| 1 | Wall |
| 2 | Fluid |
| 10 | Inlet: the opening at the start of the centerline |
| 11 | The smallest outlet |
| 12, 13, ... | The other outlets, by decreasing radius |
| 20 | Fluid voxels inside the flow diverter (visualization; the solver treats them as fluid) |

## Known limitations

* Openings must be flat cuts that lie on a face of the axis-aligned bounding box,
  away from its edges and corners.
* Triangles should not be much smaller than a voxel. "is the geometry watertight?"
  messages are printed for scan lines that pass through the open cuts; they are only
  a problem if walls appear inside the fluid.
* The scan-line voxelizer allocates `int(extent / dx)` samples per axis, so up to
  about one voxel layer at the maximum y and z sides of the vessel is not voxelized.
* The voxelizer is slow for large targets (tens of millions of elements).
