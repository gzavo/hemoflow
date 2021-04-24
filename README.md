# HemoFlow

Macroscopic flow simulation, aimed at vessel simulations. It relies on the 'preprocessorlbm' package to prepare the voxelized simulation domain.

The input parameters are read from an XML descriptor file. By default every path is relative to the XML file location.

The code supports MPI execution.

# Note
- All input parameters should be either SI or non-dimensional!

# Setup
The code builds on the Palabos open-source code. If not present, copy it to the 'palabos' directory, or use the 'setup.sh' script to clone it from the repository.
Afterwards use CMake to build the executable, e.g.:
> mkdir build
> cd build
> cmake ..
> make -j 4

## Shortcomings
- Openings must be on the axis aligned (AA) bounding box border for now to make geometry preparation automatic.
- An opening cannot fall to an edge or corner of the AA bounding box (or it can be detected on the wrong side).


## TODO
- [X] New inlet scale function (a more realistic one)
- [X] New outlet pressure distribution based on Murray-law (smallest outlet -> p=0)
- [X] Pass the angle of the openings based on centerline calculations (new ID / opening, every centerline goes from the inlet to an opening)
- [X] Calculate proper axis aligned Pouseuille profile even if the boundary is not perpendicular.
- [] Validate and verify

