#ifndef __GLOBALS_H__
#define __GLOBALS_H__

#include "palabos3D.h"
#include "palabos3D.hh"

using namespace plb;

typedef double T;

#define DESCRIPTOR descriptors::ForcedD3Q19Descriptor

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

#endif