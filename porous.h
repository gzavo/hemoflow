#ifndef POROUS_H
#define POROUS_H

#include "globals.h"
#include "helper.h"


// Initializes the porosity field from the 'stent' array of the geometry file.
// Porosity is only assigned to fluid voxels (label >= FLUID); wall and outside
// voxels never receive porous forcing.
template<typename T, typename T_>
class InitializePorousField : public BoxProcessingFunctional3D_N<T> {
    public:
        InitializePorousField(T_ *porousGeometryData, const unsigned short *geometryFlags)
            : porousData(porousGeometryData), flags(geometryFlags)
        { }

    virtual void process(Box3D domain, NTensorField3D<T>& field)
    {
        Dot3D offset = field.getLocation();

        for (plint iX=domain.x0; iX<=domain.x1; ++iX)
            for(plint iY=domain.y0; iY<=domain.y1; ++iY) 
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    plint absX = offset.x + iX;
                    plint absY = offset.y + iY;
                    plint absZ = offset.z + iZ;

                    if (absX < 0 || absY < 0 || absZ < 0 || absX >= Nx || absY >= Ny || absZ >= Nz) {
                        *field.get(iX, iY, iZ) = 0;
                        continue;
                    }
                    long long id = gT(absX, absY, absZ);
                    T porosity = porousData[id];
                    if(porosity > 0.0 && flags[id] >= FLUID) {
                        *field.get(iX, iY, iZ) = porosity;
                    }
                    else {
                    	*field.get(iX, iY, iZ) = 0;
                    }
                }  
    }

    virtual InitializePorousField<T,T_>* clone() const
    {
        return new InitializePorousField<T,T_>(*this);
    }

    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        modified[0] = modif::staticVariables;
    }

    virtual BlockDomain::DomainT appliesTo() const {
        return BlockDomain::bulkAndEnvelope;
    }

    private:
        T_ *porousData;
        const unsigned short *flags;
};



/// Apply macroscopic porous material forcing
template<typename T, template<typename U> class Descriptor>
struct PorousForceFunctional : public BoxProcessingFunctional3D_LN<T, Descriptor, T> {

	public:
        PorousForceFunctional(T linCoeff_, T quadCoeff_) : linCoeff(linCoeff_), quadCoeff(quadCoeff_)
        { }

    virtual void process(Box3D domain, BlockLattice3D<T,Descriptor>& lattice,
                                       NTensorField3D<T>& field)
    {
        // The lattice and the porosity field blocks need not coincide (e.g. sparse decomposition)
        Dot3D offset = computeRelativeDisplacement(lattice, field);

    	for (plint iX=domain.x0; iX<=domain.x1; ++iX) 
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) 
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                
                T porosity = *field.get(iX + offset.x, iY + offset.y, iZ + offset.z);

                if(porosity == 0.0)
                	continue;
                	
                if(lattice.get(iX,iY,iZ).getDynamics().isBoundary())
                    continue;

                //If it has porosity
                Array<T,Descriptor<T>::d> vel;
            	lattice.get(iX,iY,iZ).computeVelocity(vel);

            	T *force = lattice.get(iX, iY, iZ).getExternal(Descriptor<T>::ExternalField::forceBeginsAt);

            	// Calculate the force (momentum loss)
            	for (pluint iD = 0; iD < Descriptor<T>::d; ++iD) 
            	{
                	force[iD] = -porosity * (linCoeff * vel[iD] + quadCoeff * copysign(vel[iD] * vel[iD], vel[iD]) );
            	}
            }
        
    }

    virtual PorousForceFunctional<T, Descriptor>* clone() const
    {
        return new PorousForceFunctional<T, Descriptor>(*this);
    }

    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        modified[0] = modif::staticVariables;
    }

    virtual BlockDomain::DomainT appliesTo() const {
        return BlockDomain::bulkAndEnvelope;
    }

	private:
		T linCoeff;
		T quadCoeff;
};


#endif
