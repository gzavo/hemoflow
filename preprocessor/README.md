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
* `<output_base_name>geometry.vti`: the voxel labels as VTK image data in the
  coordinates of the input STL; overlays the STL and the solver output in ParaView
* `<output_base_name>geometry.nrrd`: the same label volume for 3D Slicer (spacing and
  origin in STL coordinates)

## What it does

1. `voxelizeStl.py`: reads the STL with VTK, closes open ends with flat caps,
   voxelizes the closed surface (`vtkPolyDataToImageStencil`, voxel centres over
   the whole bounding box), finds the caps on the bounding-box faces and
   voxelizes flow diverter shells by ray casting along the three axes.
2. `createFluidSolid.py`: adds a one-voxel wall around the fluid and removes empty
   outer layers; also the leak and connectivity checks.
3. `readCL.py`: reads the opening radii and directions from the centerline end
   points (inlet first, outlets by decreasing radius).
4. `detectOpenings.py`: turns the wall voxels in front of each cap into opening
   voxels.
5. `main.py`: reads the config, matches the caps to the centerline end points,
   restricts the flow diverter to the lumen, computes the physical origin of the
   grid and writes the outputs.

A 20 M voxel geometry takes a few seconds.

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

* Openings must be flat cuts (open or capped) on faces of the axis-aligned
  bounding box. They may touch its edges and corners.
* Flat strips of the vessel wall that are flush with a face are told apart from
  caps by the angle at which the neighbouring wall meets them; a cap whose wall
  meets it almost tangentially (an extremely oblique cut) would be ignored.
