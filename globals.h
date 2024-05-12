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

// TODO: scalar field for rheModel implementation
#define CELLDESCRIPTOR descriptors::D3Q7Descriptor

// Enable Large Eddy simulation (constant Smagorinsky)?
#define LES 0
 
enum GeometryLabel {
    UNUSED = 0,
    WALL = 1,
    FLUID = 2,
    FIRST_OPENING = 10      // Opening IDs go up from 10. Usually 10 is an inlet, but it is not necessary anymore.
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

// Physical units
const double BLOOD_DENSITY = 1055;  // [kg/m3]

#endif