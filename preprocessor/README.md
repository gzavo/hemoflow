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

## Installation
After cloning the directory use miniconda to create the required environment (run in repository directory):

    conda env create -f environment.yml
    conda activate lbmpre


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

## Coil porosity mode (optional)

With a coil/stent STL configured (`stent_mesh_base`), two extra config keys switch on porosity mode:

    "stent_porosity_subvoxels" : "5",   sub-voxels per voxel edge used to measure the coil
    "stent_porosity_seed" : "0",        random seed (default 0), so the draw is reproducible

The coil STL (must be watertight) is voxelized on a grid `stent_porosity_subvoxels` times finer than
the CFD grid; the share of sub-voxels inside the coil is the solid fraction of each voxel, and its
porosity is N = 1 - that. A fluid voxel then becomes a coil voxel (`stent` = 1, label 20) with
probability 1 - N. Wall and opening voxels are never changed. Use an odd number of sub-voxels: the
sub-voxels are then exactly centred on each voxel. The porosity is also written to
`<output_base_name>stent_porosity.nrrd`. See `stentPorosity.py`.

Three more keys (porosity mode only) choose other coil models from the same solid fraction:

    "stent_porosity_mode" : "draw",       draw (default): probability 1 - N
                                          threshold: solid fraction >= stent_porosity_threshold
                                          touched: solid fraction > 0, every voxel the wire reaches
    "stent_porosity_threshold" : "0.5",   used by "threshold" (default 0.5)
    "stent_as_wall" : "1",                coil voxels become wall (label 1, bounce-back) instead of
                                          porous voxels; the .npz then has no stent arrays

All modes only change fluid voxels.

Without these keys every voxel the coil surface touches is a coil voxel, as before.
