#ifndef __GLOBALS_H__
#define __GLOBALS_H__

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
//ForcedCarreauDynamics

#define CELLDESCRIPTOR descriptors::D3Q7Descriptor

// Enable Large Eddy simulation (constant Smagorinsky)?
#define LES 0

const T U_AVG_LB = 0.05;     // Re is computed in relation to this! This is the average velocity on the inlet, when the inlet flow function == 1.0
 
enum GeometryLabel {
    UNUSED = 0,
    WALL = 1,
    FLUID = 2,
    FIRST_OPENING = 10      // Opening IDs go up from 10. Usually 10 is an inlet, but it is not necessary anymore.
};

enum OpeningType {
    OPENING_VELOCITY,
    OPENING_PRESSURE,
    OUTLET_FREEFLOW 
};

// Simulation domain size
extern int Nx;
extern int Ny;
extern int Nz;

// Physical units
const double BLOOD_DENSITY = 1055;  // [kg/m3]

#endif