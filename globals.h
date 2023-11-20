#ifndef __GLOBALS_H__
#define __GLOBALS_H__

#include "palabos3D.h"
#include "palabos3D.hh"

#include <H5DataSet.hpp>
#include <H5DataSpace.hpp>
#include <H5File.hpp>
#include <H5PropertyList.hpp>
#include <hdf5.h>
#include <cstdint>

using namespace plb;

typedef double T;

#define DESCRIPTOR descriptors::ForcedD3Q19Descriptor
//ForcedMRTD3Q19Descriptor
//ForcedD3Q19Descriptor

typedef GuoExternalForceCompleteRegularizedBGKdynamics<T,DESCRIPTOR> BackgroundDynamics;
//GuoExternalForceCompleteRegularizedBGKdynamics
//GuoExternalForceMRTdynamics
//GuoExternalForceBGKdynamics
//ForcedCarreauDynamics

#define CELLDESCRIPTOR descriptors::D3Q7Descriptor

// const T U_AVG_LB = 0.1;     // Re is computed in relation to this! This is the average velocity on the inlet, when the inlet flow function == 1.0

// Geometry labels
const int UNUSED = 0;
const int WALL = 1;
const int FLUID = 2;
const int INLET_SVC = 10;
const int INLET_IVC = 11;
const int INLET_RHV = 12;
const int INLET_MHV = 13;
const int INLET_LHV = 14;
const int OUTLET_RPA = 15;
const int OUTLET_LPA = 16;
const int OUTLET_SIDE_V = 17;
const int OUTLET_SIDE_P = 18;	// Smallest outlet has zero pressure BC

// Simulation domain size
extern int Nx;
extern int Ny;
extern int Nz;

// Physical units
const double BLOOD_DENSITY = 1055;  // [kg/m3]

#endif