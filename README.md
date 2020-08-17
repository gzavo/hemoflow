# LBM_preprocess tools for HemoFlow

## Usage

See main.py for execution arguments.

Input: JSON config file
Output: voxelized geometry

Notes: 
- The nrrd output is saved since it can be opened in 3DSlicer to investigate the results of the voxelisation.
- The centerline must contain one line per outlet, with each line running from the inlet (source) to one outlet (target)
.
### Specify sides that should contain an opening.
This will practically cut away a layer of voxels from these sides.

Cutlist meaning -> cut one layer from the planes:
    # 0,1 => Xmin, Xmax
    # 2,3 => Ymin, Ymax
    # 4,5 => Zmin, Zmax

### Dependencies

Tested with Anaconda and Python 3.7

- pynrrd (PIP/conda-forge)
- numpy-stl (PIP/conda-forge)

### Known problems

- The triangle size cannot be smaller than the projected voxel size. I.e.: If you see walls appearing inside the fluid domain, or get complaints that the thing is not watertight (but you are sure it is), increase the resolution (also reduce the number of faces on the geometry?).

- The inlet and the pressure outlet are automatically detected as the largest and smallest are opening. Due to the various cut angles, this can be incorrect! Look for a better method based on centerline!

- Also automate the cutlist based on the centerline!

## Flag table

- 0 unused
- 1 wall
- 2 fluid
- 3 porous
- 10 inlet (velocity)
- 11 outlet 1 (pressure)
- 12... other outlets (velocity)

