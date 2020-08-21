#ifndef __OPENING_H__
#define __OPENING_H__

#include "globals.h"
#include "helper.h"


// Define the opening handler

class OpeningHandler
{
public:
    OpeningHandler(unsigned short* flagAray, int flag_, T radius_lb, vec3d direction);
    ~OpeningHandler();

    void printOpeningDetails();
    void loadScaleFunction(string fileName);
    void setBC(MultiBlockLattice3D<T, DESCRIPTOR> *lattice, OnLatticeBoundaryCondition3D<T, DESCRIPTOR> *bc);
    void imposeBC(MultiBlockLattice3D<T, DESCRIPTOR> *lattice, T dt);
    void setVelocityProfile(MultiBlockLattice3D<T, DESCRIPTOR> *lattice, OnLatticeBoundaryCondition3D<T, DESCRIPTOR> *bc, field3D &velocityArr);
    void setPressureProfile(MultiBlockLattice3D<T, DESCRIPTOR> *lattice, OnLatticeBoundaryCondition3D<T, DESCRIPTOR> *bc, scalar3D &pressureArr);
    void createConstantPressureProfile(T density = 1.0);
    void createPoiseauilleProfile(T u_avg = U_AVG_LB);
    void createBluntVelocityProfile(T u_avg = U_AVG_LB);

    int getFlag() { return flag; }
    int getSurfaceSize() { return nodes.size(); }

    T getRadius() { return R; }

    vec3d getCenter() { return center; }

    Box3D *getBoundingBox() { return boundingBox; }

private:

    int flag;
    bool hasScaleFunction;
    
    vector<Index3D> nodes;
    vec3d center;       // LBM units
    vec3d direction;    // Normal vector
    T R;                // LBM units
    Box3D *boundingBox; 

    // Boundary condition values
    field3D velArr;
    scalar3D presArr;

    // Scale signal
    vector<T> scaleSignal;
    vector<T> scaleTime;

    int cTimePos;
    T cTimeVal;

};

// 3D functor to get velocity
template<typename T_, template<typename U> class Descriptor>
class VelocityProfile3D {
    public:
        VelocityProfile3D (field3D *velocityValues, T_ scaleVelocity) : velocity ( velocityValues ), scale ( scaleVelocity)
        { }

        void operator() (plint iX, plint iY, plint iZ, Array<T, 3>& u) const {
            u[0] = (*velocity)[iX][iY][iZ].x * scale;
            u[1] = (*velocity)[iX][iY][iZ].y * scale;
            u[2] = (*velocity)[iX][iY][iZ].z * scale;
        }

    private:
        field3D *velocity;
        T_ scale;
};

// 3D functor to extract pressure
template<typename T_, template<typename U> class Descriptor>
class PressureProfile3D {
    public:
        PressureProfile3D (scalar3D *pressureValues, T_ scalePressure) : pressure ( pressureValues ), scale (scalePressure)
        { }

        T_ operator() (plint iX, plint iY, plint iZ) const {       
            T_ prescPressure = (*pressure)[iX][iY][iZ];
            if(prescPressure <= 0)
                prescPressure = (T_)1.0;
            return prescPressure*scale;
        }

    private:
        scalar3D *pressure;
        T_ scale;
};

#endif