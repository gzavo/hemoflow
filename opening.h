#ifndef OPENING_H
#define OPENING_H

#include "globals.h"
#include "helper.h"


// Define the opening handler

class OpeningHandler
{
public:
    OpeningHandler(const unsigned short* flagAray, GeometryLabel flag_, OpeningType type_, T radius_lb, vec3d direction);
    ~OpeningHandler();

    void printOpeningDetails(SimPar s);
    void loadScaleFunction(const string& fileName);    // Load in the scale function form a text file and MPI boradcast it
    void setBCParameter(T parameter_, SimPar s);  // Set Q or p or other BC parameter, convert it to LBM units.
    void setBCType(MultiBlockLattice3D<T, DESCRIPTOR> *lattice);    // Sets the BC type for LBM
    void progressTime(MultiBlockLattice3D<T, DESCRIPTOR> *lattice, T dt);     // Progress time by dt and impose time-dependent values on the opening
    void progressWarmup(MultiBlockLattice3D<T, DESCRIPTOR> *lattice, T dt, T rampingInterval);     // Progress warmup ramping
    void setExternalVelocityProfile(MultiBlockLattice3D<T, DESCRIPTOR> *lattice, field3D &velocityArr);     // Overwrites profile with external array
    void setExternalPressureProfile(MultiBlockLattice3D<T, DESCRIPTOR> *lattice, scalar3D &pressureArr);    // Overwrites profile with external array
    void setScaledBoundaryProfile(MultiBlockLattice3D<T, DESCRIPTOR> *lattice, T scale);
    void createConstantPressureProfile(T density = 1.0);
    void createPoiseauilleProfile();    // norm(v_max) = 1.0
    void createBluntVelocityProfile();  // norm(v) = 1.0
    void normalizeFlowRate();           // Set Q = 1.0
    void scaleFlowRate(T scale_);       // Scale velocity
    void scalePressure(T scale_);       // Scale pressure
    void setBCParameter (T bcParameter_) { bcParameter = bcParameter_; }

    int getGeometryLabel() { return flag; }
    int getOpeningType() { return type; }
    int getSurfaceSize() { return nodes.size(); }
    bool getIsMurrayOpening() { return type==OPENING_MURRAY; }
    T getProfileVelSum();
    T getFlowRate(SimPar s);
    T getScaledFlowRate();
    T getRadius() const { return R; }
    T getArea() const { return nodes.size(); }
    vec3d getCenter() { return center; }
    Box3D *getBoundingBox() { return boundingBox; }

    void setName(string name_) {name = name_;}
    string getName () {return name;}

private:
    string name;
    GeometryLabel flag;
    OpeningType type;
    
    vector<Index3D> nodes;  // List of LBM nodes on the opening
    vec3d center;           // LBM units
    vec3d direction;        // Normal vector
    T R;                    // LBM units (area-derived hydrodynamic radius)
    Box3D *boundingBox;     // BB for functionals

    // Boundary condition profiles
    field3D velArr;
    scalar3D presArr;
    T bcParameter;    // Normal velocity magnitude max (for Q) or absolute pressure (p). These multiply the scale functions.

    // Scale signal
    bool hasScaleFunction;
    vector<T> scaleSignal;
    vector<T> scaleTime;
    T cScale; // The current scale value of the function

    // Keeping track of the scale function (and looping it)
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