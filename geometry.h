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
    bool hasMultipliers;         // Per-voxel resistance multipliers (older files: 1 in every device voxel)
    long long numPorous;         // Porous voxels: flow diverter voxels in the fluid
    Box3D porousBox;             // Bounding box of the porous voxels (if numPorous > 0)
    Array<double, 2> multiplierMin, multiplierMax, multiplierMean;   // (linear, quadratic) over the porous voxels

    std::vector<int> openingLabels;
    std::vector<double> openingRadius;      // [m]
    std::vector<vec3d> openingNormal;       // Unit vectors pointing into the domain
    std::map<int, std::vector<Index3D> > openingNodes;  // Voxels of each opening (all processes)

    // Fill a distributed field (any block structure) from the voxel arrays.
    void fillLabels(MultiScalarField3D<int>& field) const;      // Voxel labels
    void fillOccupancy(MultiScalarField3D<int>& field) const;   // 1 for wall, fluid or opening voxels
    // (linear, quadratic) resistance multipliers in the porous voxels, 0 elsewhere
    void fillPorosity(MultiTensorField3D<T, 2>& field) const;

    // Free the voxel arrays of the main process once all fields are filled.
    void release();

private:
    bool isPorous(size_t id) const;                              // Main process only
    std::vector<double> summarizeDevice(bool withMultipliers);   // Main process only

    std::string fileName;
    cnpy::NpyArray labels;   // Main process only
    cnpy::NpyArray stent;    // Main process only (if hasStent)
    cnpy::NpyArray linearMultiplier, quadraticMultiplier;   // Main process only (if hasMultipliers)
};

#endif
