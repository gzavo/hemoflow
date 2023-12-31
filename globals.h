#ifndef __GLOBALS_H__
#define __GLOBALS_H__

#include "palabos3D.h"
#include "palabos3D.hh"

using namespace plb;

typedef double T;

#define DESCRIPTOR descriptors::ForcedD3Q19Descriptor
//ForcedMRTD3Q19Descriptor
//ForcedD3Q19Descriptor

typedef GuoExternalForceBGKdynamics<T,DESCRIPTOR> BackgroundDynamics;
//GuoExternalForceCompleteRegularizedBGKdynamics
//GuoExternalForceMRTdynamics
//GuoExternalForceBGKdynamics
//ForcedCarreauDynamics

#define CELLDESCRIPTOR descriptors::D3Q7Descriptor

const T U_AVG_LB = 0.05;     // Re is computed in relation to this! This is the average velocity on the inlet, when the inlet flow function == 1.0

// Geometry labels
const int UNUSED = 0;
const int WALL = 1;
const int FLUID = 2;
const int INLET = 10;			// This is always a velocity inlet
const int FIRST_OUTLET = 11;	// This is always the smallest outlet defined as pressure

// Simulation domain size
extern int Nx;
extern int Ny;
extern int Nz;

// Physical units
const double BLOOD_DENSITY = 1055;  // [kg/m3]

#endif