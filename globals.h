#ifndef GLOBALS_H
#define GLOBALS_H

#include "palabos3D.h"
#include "palabos3D.hh"

#include <cstdint>

using namespace plb;

typedef double T;

#define DESCRIPTOR descriptors::ForcedD3Q19Descriptor
//ForcedMRTD3Q19Descriptor
//ForcedD3Q19Descriptor

//typedef GuoExternalForceBGKdynamics<T,DESCRIPTOR> BackgroundDynamics;
typedef GuoExternalForceCompleteRegularizedBGKdynamics<T,DESCRIPTOR> BackgroundDynamics;
//GuoExternalForceCompleteRegularizedBGKdynamics
//GuoExternalForceMRTdynamics
//GuoExternalForceBGKdynamics
// Non-Newtonian (Carreau-Yasuda) bulk dynamics: CarreauBackgroundDynamics in nonNewtonian.h

// TODO: scalar field for rheModel implementation
#define CELLDESCRIPTOR descriptors::D3Q7Descriptor

// Enable Large Eddy simulation (constant Smagorinsky)?
// Note: the Smagorinsky dynamics has no external force term, so the porous
// (flow diverter) forcing has no effect when LES is enabled.
#define LES 0
#if LES
const T cSmago = 0.14;  // Smagorinsky constant
#endif

// Voxel labels in the geometry file (see preprocessor/constants.py)
enum GeometryLabel {
    UNUSED = 0,
    WALL = 1,
    FLUID = 2,
    FIRST_OPENING = 10,     // Openings are labelled 10, 11, 12, ... (10 is the centerline inlet)
    FLOW_DIVERTER = 20      // Fluid voxels inside a flow diverter (visualization label, treated as fluid)
};

enum OpeningType {
    OPENING_VELOCITY = 1,
    OPENING_MURRAY = 2,
    OPENING_PRESSURE = 3,
    OUTLET_FREEFLOW = 4
};

// Simulation domain size
extern int Nx;
extern int Ny;
extern int Nz;

// Default physical properties of blood (can be overridden in the <fluid> section of the config)
const double DEFAULT_BLOOD_DENSITY = 1055;                // [kg/m3]
const double DEFAULT_BLOOD_KINEMATIC_VISCOSITY = 3.22e-6; // [m2/s]
// Carreau model of blood (Cho and Kensey, Biorheology 28, 241, 1991), used with <model>carreau</model>
const double DEFAULT_CARREAU_MU0 = 0.056;       // Zero-shear viscosity [Pa s]
const double DEFAULT_CARREAU_MU_INF = 0.00345;  // Infinite-shear viscosity [Pa s]
const double DEFAULT_CARREAU_LAMBDA = 3.313;    // Relaxation time [s]
const double DEFAULT_CARREAU_N = 0.3568;        // Power-law index

#endif
