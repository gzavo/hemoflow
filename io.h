#ifndef IO_H
#define IO_H

#include "globals.h"
#include "helper.h"

// For directory manipulations
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>

#ifdef HEMOFLOW_WITH_HDF5
// HDF5 and HighFive includes
#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5File.hpp>
#include <highfive/H5PropertyList.hpp>
#include <hdf5.h>
#endif

#include "cnpy.h"
#include "io/xdmfDataOutput.h"

bool fileExists (const std::string& name);
int dirExists(const string& pathName);
int do_mkdir(const char *path, mode_t mode);
int mkpath(const char *path, mode_t mode);

// Fields written to the output files
struct OutputOptions {
    bool velocity = true;   // velocity [m/s]
    bool pressure = true;   // gauge pressure [Pa] (relative to the reference density)
    bool stress = false;    // viscous (deviatoric) stress tensor [Pa]: xx, xy, xz, yy, yz, zz
    bool viscosity = false; // dynamic viscosity [Pa s] (constant for a Newtonian fluid)
    Array<double, 3> origin = Array<double, 3>(0.0, 0.0, 0.0);  // Physical position of lattice node (0,0,0) [m]
    double lengthUnit = 1.0;    // Length of one unit of the output grid coordinates [m] (0.001: mm, as the input STL)
};

void writeVTK(MultiBlockLattice3D<T,DESCRIPTOR>& lattice, const SimPar &sim, plint iter, const OutputOptions &out);
void writeNPZ(MultiBlockLattice3D<T,DESCRIPTOR>& lattice, const SimPar &sim, plint iter);  // Velocity [m/s] as (Nx, Ny, Nz, 3)
#ifdef HEMOFLOW_WITH_HDF5
void writeHDF5(MultiBlockLattice3D<T,DESCRIPTOR>& lattice, const SimPar &sim, plint iter, string outDir,
               const Array<double, 3> &origin, MultiNTensorField3D<T> *field1 = nullptr);
#endif

#endif
