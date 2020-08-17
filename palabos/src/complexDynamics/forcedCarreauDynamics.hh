
#ifndef EXTERNAL_FORCE_CARREAU_DYNAMICS_HH
#define EXTERNAL_FORCE_CARREAU_DYNAMICS_HH

#include "basicDynamics/isoThermalDynamics.h"
#include "core/cell.h"
#include "latticeBoltzmann/dynamicsTemplates.h"
#include "latticeBoltzmann/momentTemplates.h"
#include "latticeBoltzmann/offEquilibriumTemplates.h"
#include "latticeBoltzmann/d3q13Templates.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"
#include "latticeBoltzmann/externalForceTemplates.h"
#include "complexDynamics/forcedCarreauDynamics.h"
#include "core/latticeStatistics.h"
#include <algorithm>
#include "core/dynamicsIdentifiers.h"
#include <limits>

namespace plb {


template<typename T, template<typename U> class Descriptor>
int ForcedCarreauDynamics<T,Descriptor>::id =
    meta::registerGeneralDynamics<T,Descriptor,ForcedCarreauDynamics<T,Descriptor> > (std::string("Forced_CarreauDynamics_BGK_"));

/** \param omega_ relaxation parameter, related to the dynamic viscosity
 */
template<typename T, template<typename U> class Descriptor>
ForcedCarreauDynamics<T,Descriptor>::ForcedCarreauDynamics(T omega_ )
    : ExternalForceDynamics<T,Descriptor>(omega_)
{ }

template<typename T, template<typename U> class Descriptor>
ForcedCarreauDynamics<T,Descriptor>::ForcedCarreauDynamics(HierarchicUnserializer& unserializer)
    : ExternalForceDynamics<T,Descriptor>(T())
{
    this->unserialize(unserializer);
}

template<typename T, template<typename U> class Descriptor>
ForcedCarreauDynamics<T,Descriptor>* ForcedCarreauDynamics<T,Descriptor>::clone() const {
    return new ForcedCarreauDynamics<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
int ForcedCarreauDynamics<T,Descriptor>::getId() const {
    return id;
}

template<typename T, template<typename U> class Descriptor>
void ForcedCarreauDynamics<T,Descriptor>::collide (
        Cell<T,Descriptor>& cell,
        BlockStatistics& statistics )
{
    
    T nu0 = global::CarreauParameters().getNu0();
    T nuInf  = global::CarreauParameters().getNuInf();
    T nMinusOneOverTwo  = (global::CarreauParameters().getExponent() - (T)1)/(T)2;
    T lambda  = global::CarreauParameters().getLambda();

    T rhoBar;
    Array<T,Descriptor<T>::d> j;
    Array<T,SymmetricTensor<T,Descriptor>::n> PiNeq;
    momentTemplates<T,Descriptor>::compute_rhoBar_j_PiNeq(cell, rhoBar, j, PiNeq);
    
    T omega0 = this->getOmega();
    // T rateOfStrainNormSqr = SymmetricTensor<T,Descriptor>::tensorNormSqr(-1.5*omega0*PiNeq);
    T rateOfStrainNormSqr = SymmetricTensor<T,Descriptor>::tensorNormSqr(-omega0*Descriptor<T>::invCs2*Descriptor<T>::invRho(rhoBar)*0.5*PiNeq);
    T nuEff = nuInf * (nu0-nuInf)*pow(1 + lambda*lambda*2.0*rateOfStrainNormSqr, nMinusOneOverTwo);
    
    omega0 = 1.0/(3*nuEff + 0.5);
    
    this->setOmega(omega0);
    

    Array<T,Descriptor<T>::d> u;
    this->computeVelocity(cell, u);
    // T rho = Descriptor<T>::fullRho(rhoBar);

    // for (plint iD = 0; iD < Descriptor<T>::d; ++iD)
    // {
    //     j[iD] = rho * u[iD];
    // }
    
    T uSqr = dynamicsTemplates<T,Descriptor>::bgk_ma2_collision(cell, rhoBar, j, this->getOmega());
    externalForceTemplates<T,Descriptor>::addGuoForce(cell, u, this->getOmega(), (T)1);
    
    if (cell.takesStatistics()) {
        gatherStatistics(statistics, rhoBar, uSqr);
    }
}


template<typename T, template<typename U> class Descriptor>
T ForcedCarreauDynamics<T,Descriptor>::computeEquilibrium (
        plint iPop, T rhoBar, Array<T,Descriptor<T>::d> const& j,
        T jSqr, T thetaBar) const
{
    T invRho = Descriptor<T>::invRho(rhoBar);
    return dynamicsTemplates<T,Descriptor>::bgk_ma2_equilibrium(iPop, rhoBar, invRho, j, jSqr);
}

}  // namespace plb

#endif