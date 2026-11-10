#ifndef POROUS_H
#define POROUS_H

#include "globals.h"
#include "helper.h"


// The multiplier field holds (linear, quadratic) resistance multipliers in the porous
// voxels and 0 elsewhere. GeometryFile::fillPorosity() fills it from the geometry file
// (1 in every porous voxel for a uniform device or an older geometry file).

/// Apply macroscopic porous material forcing: in each porous cell the force density
/// -(mLin a + mQuad b |u|) u / dx, with a = linCoeff and b = quadCoeff (lattice units).
/// The force is antiparallel to u and does not depend on the flow direction relative to the grid.
template<typename T, template<typename U> class Descriptor>
struct PorousForceFunctional : public BoxProcessingFunctional3D_LT<T, Descriptor, T, 2> {

	public:
        PorousForceFunctional(T linCoeff_, T quadCoeff_) : linCoeff(linCoeff_), quadCoeff(quadCoeff_)
        { }

    virtual void process(Box3D domain, BlockLattice3D<T,Descriptor>& lattice,
                                       TensorField3D<T,2>& multipliers)
    {
        // The lattice and the multiplier field blocks need not coincide (e.g. sparse decomposition)
        Dot3D offset = computeRelativeDisplacement(lattice, multipliers);

    	for (plint iX=domain.x0; iX<=domain.x1; ++iX) 
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) 
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                
                Array<T,2> const& m = multipliers.get(iX + offset.x, iY + offset.y, iZ + offset.z);

                if(m[0] == 0.0 && m[1] == 0.0)
                	continue;
                	
                if(lattice.get(iX,iY,iZ).getDynamics().isBoundary())
                    continue;

                //If it has porosity
                Array<T,Descriptor<T>::d> vel;
            	lattice.get(iX,iY,iZ).computeVelocity(vel);

            	T *force = lattice.get(iX, iY, iZ).getExternal(Descriptor<T>::ExternalField::forceBeginsAt);

            	// Calculate the force (momentum loss)
            	T speed = std::sqrt(VectorTemplate<T, Descriptor>::normSqr(vel));
            	T resistance = m[0] * linCoeff + m[1] * quadCoeff * speed;
            	for (pluint iD = 0; iD < Descriptor<T>::d; ++iD) 
            	{
                	force[iD] = -resistance * vel[iD];
            	}
            }
        
    }

    virtual PorousForceFunctional<T, Descriptor>* clone() const
    {
        return new PorousForceFunctional<T, Descriptor>(*this);
    }

    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        modified[0] = modif::staticVariables;
        modified[1] = modif::nothing;
    }

    virtual BlockDomain::DomainT appliesTo() const {
        return BlockDomain::bulkAndEnvelope;
    }

	private:
		T linCoeff;
		T quadCoeff;
};


#endif
