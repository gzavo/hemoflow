#ifndef __HELPER_H__
#define __HELPER_H__

#include <map>
#include <algorithm>
#include <cstdlib>
#include <vector>
#include <math.h>

using namespace std;

#include "globals.h"

// *** Function to convert array coordinates

// Fortran ordering - numpy
inline int gT(int x, int y, int z) {
    return x*Nz*Ny + y*Nz + z;
}

inline int gT2D(int size_y ,int x, int y){
    return x*size_y + y;
}

// C ordering
inline int gTC(int x, int y, int z) {
    return x + y*Nx + z*Nx*Ny;
}

// Fortran ordering - Array index conversion between 1D - 3+ D arrays
inline int gT(int d, int x, int y, int z) {
    return d*Nx*Ny*Nz + x*Nz*Ny + y*Nz + z;
}

// *** Simple linear interpolation
inline T interpolate(T x1, T x2, T xi, T y1, T y2)
{
    T r = (xi-x1)/(x2-x1);
    return y1+r*(y2-y1);
}


// *** Data structures

struct Index3D 
{ 
   int x, y, z; 
}; 

class vec3d 
{
public:
    double x, y, z;
    vec3d() {}
    vec3d(double x, double y, double z) {set(x,y,z);}
    
    void set(double X, double Y, double Z) {x=X; y=Y; z=Z;}
    void set(vec3d v) {x=v.x; y=v.y; z=v.z;}
    double norm() {return sqrt(x*x+y*y+z*z);}
    void normalize() {double n = norm(); if(n==0) return; x/=n; y/=n; z/=n;}
    void negate() {x=-x;y=-y;z=-z;}
    double dot(vec3d v) { return x*v.x + y*v.y + z*v.z; }
    vec3d cross(vec3d v) { return vec3d(y*v.z-z*v.y, z*v.x - x*v.z, x*v.y - y*v.x);}
    double operator[] ( int i ) { if ( i == 0 ) return x; else if ( i == 1 ) return y; else return z; }

    vec3d operator+ (vec3d v) {vec3d t; t.x=x+v.x; t.y=y+v.y; t.z=z+v.z; return t;}
    vec3d operator* (double s) {vec3d t; t.x=s*x; t.y=s*y;t.z=s*z; return t;}
    vec3d operator- (vec3d v) {vec3d t; t.x = x-v.x; t.y = y-v.y; t.z = z-v.z; return t;}
};

class symMtx3 
{
public:
    double m1,m2,m3,m4,m5,m6;
    symMtx3() {}
    symMtx3(double M1, double M2, double M3, double M4, double M5, double M6)
    {set(M1,M2,M3,M4,M5,M6);}

    double operator[] ( int i ) { 
        if ( i == 0 ) return m1;
        else if ( i == 1 ) return m2;
        else if ( i == 2 ) return m3;
        else if ( i == 3 ) return m4;
        else if ( i == 4 ) return m5;
        else return m6; }

    void set(double M1, double M2, double M3, double M4, double M5, double M6) {m1=M1;m2=M2;m3=M3;m4=M4;m5=M5;m6=M6;}

    vec3d vecMult(vec3d v) {
        vec3d r;
        r.x = m1*v.x + m2*v.y + m3*v.z;
        r.y = m2*v.x + m4*v.y + m5*v.z;
        r.z = m3*v.x + m5*v.y + m6*v.z;
        return r;
    }
};
typedef map<int, map<int, map<int, double > > > scalar3D;
typedef map<int, map<int, map<int, vec3d > > > field3D;
typedef map<int, map<int, map<int, symMtx3 > > > tensor3D;


// *** Functors to handle arrays in a parallel fashion

// 3D functor to select out all flags greater than a given flag
// This is useful to define the domain for decomposition
template<typename T_>
class FlagMaskDomain3D : public DomainFunctional3D {
    public:
        FlagMaskDomain3D(T_ *allGeometryFlags, int flagToMaskAbove = 0) : maskFlag(flagToMaskAbove), geometryFlags(allGeometryFlags)
        { }

        virtual bool operator () (plint iX, plint iY, plint iZ) const {
            if(geometryFlags[gT(iX,iY,iZ)] > maskFlag)
                return true;
            return false;
        }

        virtual FlagMaskDomain3D<T_>* clone() const {
            return new FlagMaskDomain3D<T_>(*this);
        }

    private:
        int maskFlag;
        T_ *geometryFlags;
};

// 3D functor to select out a single given flag
// This is useful to define the domain for decomposition
template<typename T_>
class FlagMaskSingleDomain3D : public DomainFunctional3D {
    public:
        FlagMaskSingleDomain3D(T_ *allGeometryFlags, int flagToMask = 0) : maskFlag(flagToMask), geometryFlags(allGeometryFlags)
        { }

        virtual bool operator () (plint iX, plint iY, plint iZ) const {
            if(geometryFlags[gT(iX,iY,iZ)] == maskFlag)
                return true;
            return false;
        }

        virtual FlagMaskSingleDomain3D<T_>* clone() const {
            return new FlagMaskSingleDomain3D<T_>(*this);
        }

    private:
        int maskFlag;
        T_ *geometryFlags;
};


#endif