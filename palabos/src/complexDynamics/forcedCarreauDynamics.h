#ifndef FORCED_CARREAU_DYNAMICS_H
#define FORCED_CARREAU_DYNAMICS_H

#include "core/globalDefs.h"
#include "core/dynamics.h"
#include "basicDynamics/externalForceDynamics.h"
#include "complexDynamics/carreauDynamics.h"

namespace plb {


/// Implementation of O(Ma^2) BGK dynamics with an external force (Guo approach)
template<typename T, template<typename U> class Descriptor>
class ForcedCarreauDynamics : public ExternalForceDynamics<T,Descriptor> {
public:
/* *************** Construction / Destruction ************************ */
    ForcedCarreauDynamics(T omega_=(T)1);
    ForcedCarreauDynamics(HierarchicUnserializer& unserializer);

    /// Clone the object on its dynamic type.
    virtual ForcedCarreauDynamics<T,Descriptor>* clone() const;

    /// Return a unique ID for this class.
    virtual int getId() const;

/* *************** Collision and Equilibrium ************************* */

    /// Implementation of the collision step
    virtual void collide(Cell<T,Descriptor>& cell,
                         BlockStatistics& statistics_);

    /// Compute equilibrium distribution function
    virtual T computeEquilibrium(plint iPop, T rhoBar, Array<T,Descriptor<T>::d> const& j,
                                 T jSqr, T thetaBar=T()) const;
private:
    static int id;
};


}

#endif