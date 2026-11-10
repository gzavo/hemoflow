#ifndef POROUS_H
#define POROUS_H

#include "globals.h"
#include "helper.h"


// The porosity field (1 in porous voxels, 0 elsewhere) is filled from the geometry
// file by GeometryFile::fillPorosity().

/// Apply macroscopic porous material forcing
template<typename T, template<typename U> class Descriptor>
struct PorousForceFunctional : public BoxProcessingFunctional3D_LS<T, Descriptor, T> {

	public:
        PorousForceFunctional(T linCoeff_, T quadCoeff_) : linCoeff(linCoeff_), quadCoeff(quadCoeff_)
        { }

    virtual void process(Box3D domain, BlockLattice3D<T,Descriptor>& lattice,
                                       ScalarField3D<T>& field)
    {
        // The lattice and the porosity field blocks need not coincide (e.g. sparse decomposition)
        Dot3D offset = computeRelativeDisplacement(lattice, field);

    	for (plint iX=domain.x0; iX<=domain.x1; ++iX) 
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) 
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                
                T porosity = field.get(iX + offset.x, iY + offset.y, iZ + offset.z);

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
