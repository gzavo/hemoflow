# HemoFlow

Macroscopic flow simulation, aimed at vessel simulations. It relies on the 'preprocessorlbm' package to prepare the voxelized simulation domain.
The input parameters are read from an XML descriptor file. By default every path is relative to the XML file location.
The code supports MPI execution and sparse domain decomposition.

# Note
All input parameters should be either SI or non-dimensional!

# Setup
## Solver
The code builds on the Palabos open-source code. If not present, copy it to the 'palabos' directory, or use the 'setup.sh' script to clone it from the repository.
Afterwards use CMake to build the executable, e.g.:

```bash
mkdir build
cd build
cmake ..
make -j 4
```

## Voxelization
The voxelization code runs on python, the environment can be set up using conda:
```bash
conda env create -f preprocessor/environment.yml
conda activate lbmpre
```
## Troubleshooting
### HDF5 problems on Ubuntu
HighFive requires 1.13+ parallel HDF5 version which is not available from apt-get (yet).
Can be compiled manually using

```bash
git clone https://github.com/HDFGroup/hdf5
cd hdf5
git checkout hdf5-1_14_0
mkdir build
cd build
cmake -G "Unix Makefiles" -DHDF5_ENABLE_PARALLEL=ON -DHDF5_ENABLE_Z_LIB_SUPPORT=ON -DCMAKE_INSTALL_PREFIX=/opt/hdf5 ..
sudo make install -j 4
```

When compiling hemoflow cmake needs help to find our custom HDF5:
```bash
cmake -DHDF5_ROOT=/opt/hdf5 ..
make -j 4
```
### Compilation on newer systems

Newer compilers can throw an error, which can be disabled in the CmakeLists adding

    -Wno-error=template-body

to the flags


## Shortcomings
- Openings must be on the axis aligned (AA) bounding box border for now to make geometry preparation automatic.
- An opening cannot fall to an edge or corner of the AA bounding box (or it can be detected on the wrong side).

## TODO
- [X] New inlet scale function (a more realistic one)
- [X] New outlet pressure distribution based on Murray-law
- [X] Pass the angle of the openings based on centerline calculations (new ID / opening, every centerline goes from the inlet to an opening)
- [X] Calculate proper axis aligned Pouseuille profile even if the boundary is not perpendicular.
- [] Bump up preprocessor to python 3.12 and corresponding numpy.
- [] Validate and verify


## License: AGPL v3.0
