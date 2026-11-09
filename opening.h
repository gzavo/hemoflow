#ifndef OPENING_H
#define OPENING_H

#include "globals.h"
#include "helper.h"

#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

// Unit-peak velocity profile of an opening, keyed by the linear voxel index gT(x,y,z).
typedef std::unordered_map<long long, vec3d> VelocityProfileMap;

// Handles one inlet/outlet ("opening") of the domain: its lattice nodes, the
// boundary-condition type, the velocity/pressure profile and the optional
// time-dependent scale function (waveform).
//
// Conventions:
//  - All openings lie on a face of the domain bounding box.
//  - 'direction' (from the centerline) and 'faceNormal' both point INTO the domain.
//  - Positive flow rates are inflows, negative flow rates are outflows.
class OpeningHandler
{
public:
    OpeningHandler(const unsigned short* flagArray, int label, OpeningType type, T radius_lb, vec3d direction);

    void printOpeningDetails(const SimPar& s) const;
    void loadScaleFunction(const std::string& fileName);  // Load the waveform from a text file (read on the main process, broadcast)
    void setBCParameter(T parameter, const SimPar& s);    // SI value (velocity [m/s] or gauge pressure [Pa]) -> lattice units
    void setBCParameterLB(T parameter_lb) { bcParameter = parameter_lb; }
    void setBCType(MultiBlockLattice3D<T, DESCRIPTOR>& lattice);   // Instantiates the boundary condition on the lattice
    void progressTime(MultiBlockLattice3D<T, DESCRIPTOR>& lattice, T dt);  // Advance time by dt and impose the boundary values
    void createConstantPressureProfile();
    void createPoiseuilleProfile();     // Peak velocity 1 (lattice units), zero on the fitted ellipse
    void createBluntVelocityProfile();  // Velocity 1 (lattice units) on every node

    int getLabel() const { return label; }
    OpeningType getOpeningType() const { return type; }
    bool isVelocityType() const { return type == OPENING_VELOCITY || type == OPENING_MURRAY; }
    size_t getNumNodes() const { return nodes.size(); }
    T getRadius() const { return R; }
    vec3d getCenter() const { return center; }
    Box3D const& getBoundingBox() const { return boundingBox; }
    int getFaceAxis() const { return faceAxis; }
    int getFaceSign() const { return faceSign; }

    // Flow rate (lattice units, inflow positive) of the profile with unit peak velocity.
    T getProfileFlux() const { return profileFlux; }
    // Flow rate currently imposed on a velocity or Murray opening (lattice units, inflow positive).
    T getImposedFluxLB() const { return isVelocityType() ? profileFlux * bcParameter * cScale : 0.0; }
    T getImposedFlowRate(const SimPar& s) const;    // [m^3/s], inflow positive
    // Flow rate measured on the lattice nodes of the opening [m^3/s], inflow positive. Collective (all MPI ranks).
    T measureFlowRate(MultiBlockLattice3D<T, DESCRIPTOR>& lattice, const SimPar& s) const;

    bool getHasScaleFunction() const { return hasScaleFunction; }
    T getPeriod() const { return hasScaleFunction ? scaleTime.back() : -1.0; }

    // True if both openings lie on the same domain face and their bounding rectangles intersect.
    bool overlapsOnFace(const OpeningHandler& other) const;

    void setName(const std::string& name_) { name = name_; }
    std::string getName() const { return name; }

private:
    const unsigned short* flags;
    std::string name;
    int label;
    OpeningType type;

    std::vector<Index3D> nodes;  // LBM nodes of the opening
    vec3d center;                // Geometric centre [lattice units]
    vec3d direction;             // Unit vector along the centerline, pointing into the domain
    T R;                         // Radius from the centerline [lattice units]
    Box3D boundingBox;
    int faceAxis;                // 0, 1, 2: the domain face (x, y, z) the opening lies on
    int faceSign;                // +1 on the min face, -1 on the max face (direction of the inward face normal)

    // Boundary condition profile
    VelocityProfileMap profile;
    T profileFlux;   // Flux of the unit-peak profile through the face [lattice units]
    T bcParameter;   // Velocity openings: peak velocity [lb]. Pressure openings: density deviation rho-1 [lb]. Scaled by the waveform.

    // Scale function (waveform)
    bool hasScaleFunction;
    std::vector<T> scaleSignal;
    std::vector<T> scaleTime;
    T cScale;      // Current value of the scale function
    int cTimePos;  // Current segment of the waveform
    T cTimeVal;    // Current time within the period
};


// Imposes the boundary value on the cells of ONE opening (cells whose label matches),
// so that openings whose bounding boxes overlap cannot overwrite each other.
template<typename T_, template<typename U> class Descriptor>
class ImposeOpeningValues3D : public BoxProcessingFunctional3D_L<T_, Descriptor> {
public:
    // Velocity opening: u = profile * velocityScale
    ImposeOpeningValues3D(const unsigned short* flags_, int label_, const VelocityProfileMap* profile_, T_ velocityScale_)
        : flags(flags_), label(label_), profile(profile_), velocityScale(velocityScale_), isVelocity(true), density(1.0)
    { }
    // Pressure opening: rho = density
    ImposeOpeningValues3D(const unsigned short* flags_, int label_, T_ density_)
        : flags(flags_), label(label_), profile(nullptr), velocityScale(0.0), isVelocity(false), density(density_)
    { }

    virtual void process(Box3D domain, BlockLattice3D<T_, Descriptor>& lattice)
    {
        Dot3D offset = lattice.getLocation();
        for (plint iX = domain.x0; iX <= domain.x1; ++iX)
            for (plint iY = domain.y0; iY <= domain.y1; ++iY)
                for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                    plint x = iX + offset.x, y = iY + offset.y, z = iZ + offset.z;
                    if (x < 0 || y < 0 || z < 0 || x >= Nx || y >= Ny || z >= Nz)
                        continue;
                    long long id = gT(x, y, z);
                    if (flags[id] != label)
                        continue;
                    Cell<T_, Descriptor>& cell = lattice.get(iX, iY, iZ);
                    if (isVelocity) {
                        typename VelocityProfileMap::const_iterator it = profile->find(id);
                        if (it == profile->end())
                            continue;
                        Array<T_, 3> u(it->second.x * velocityScale, it->second.y * velocityScale, it->second.z * velocityScale);
                        cell.defineVelocity(u);
                    }
                    else {
                        cell.defineDensity(density);
                    }
                }
    }

    virtual ImposeOpeningValues3D<T_, Descriptor>* clone() const { return new ImposeOpeningValues3D<T_, Descriptor>(*this); }
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const { modified[0] = modif::staticVariables; }
    virtual BlockDomain::DomainT appliesTo() const { return BlockDomain::bulkAndEnvelope; }

private:
    const unsigned short* flags;
    int label;
    const VelocityProfileMap* profile;
    T_ velocityScale;
    bool isVelocity;
    T_ density;
};


// Sums the velocity component normal to the face (inflow positive) over the cells of one opening.
template<typename T_, template<typename U> class Descriptor>
class OpeningFluxFunctional3D : public ReductiveBoxProcessingFunctional3D_L<T_, Descriptor> {
public:
    OpeningFluxFunctional3D(const unsigned short* flags_, int label_, int axis_, int sign_)
        : flags(flags_), label(label_), axis(axis_), sign(sign_), sumId(this->getStatistics().subscribeSum())
    { }

    virtual void process(Box3D domain, BlockLattice3D<T_, Descriptor>& lattice)
    {
        BlockStatistics& statistics = this->getStatistics();
        Dot3D offset = lattice.getLocation();
        for (plint iX = domain.x0; iX <= domain.x1; ++iX)
            for (plint iY = domain.y0; iY <= domain.y1; ++iY)
                for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                    if (flags[gT(iX + offset.x, iY + offset.y, iZ + offset.z)] != label)
                        continue;
                    Array<T_, 3> u;
                    lattice.get(iX, iY, iZ).computeVelocity(u);
                    statistics.gatherSum(sumId, sign * u[axis]);
                }
    }

    virtual OpeningFluxFunctional3D<T_, Descriptor>* clone() const { return new OpeningFluxFunctional3D<T_, Descriptor>(*this); }
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const { modified[0] = modif::nothing; }
    virtual BlockDomain::DomainT appliesTo() const { return BlockDomain::bulk; }
    T_ getFlux() const { return this->getStatistics().getSum(sumId); }

private:
    const unsigned short* flags;
    int label;
    int axis;
    int sign;
    plint sumId;
};

#endif
