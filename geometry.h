#ifndef GEOMETRY_H
#define GEOMETRY_H

#include "globals.h"
#include "helper.h"
#include "cnpy.h"

#include <map>
#include <string>
#include <vector>

// The geometry file (.npz) written by the preprocessor.
//
// Only the main MPI process reads the voxel arrays. The other processes receive the
// small arrays (resolution, origin, openings) by broadcast and, through the fill*()
// functions, only the parts of the voxel fields that belong to their own blocks.
// All member functions are collective (every MPI process must call them).
class GeometryFile {
public:
    explicit GeometryFile(const std::string& fileName);

    int nx, ny, nz;              // Domain size [voxels]
    double dx;                   // Voxel size [m]
    bool hasOrigin;              // False for files from older preprocessor versions
    Array<double, 3> origin;     // Position of voxel (0,0,0) [m]
    bool hasLengthUnit;          // False for files from older preprocessor versions
    double lengthUnit;           // Length of one unit of the input STL [m] (0.001: mm)
    bool hasStent;               // The file contains a flow diverter

    std::vector<int> openingLabels;
    std::vector<double> openingRadius;      // [m]
    std::vector<vec3d> openingNormal;       // Unit vectors pointing into the domain
    std::map<int, std::vector<Index3D> > openingNodes;  // Voxels of each opening (all processes)

    // Fill a distributed field (any block structure) from the voxel arrays.
    void fillLabels(MultiScalarField3D<int>& field) const;      // Voxel labels
    void fillOccupancy(MultiScalarField3D<int>& field) const;   // 1 for wall, fluid or opening voxels
    void fillPorosity(MultiScalarField3D<T>& field) const;      // 1 in the fluid voxels of the flow diverter

    // Free the voxel arrays of the main process once all fields are filled.
    void release();

private:
    std::string fileName;
    cnpy::NpyArray labels;   // Main process only
    cnpy::NpyArray stent;    // Main process only (if hasStent)
};

#endif
