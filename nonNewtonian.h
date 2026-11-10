#ifndef NON_NEWTONIAN_H
#define NON_NEWTONIAN_H

#include "globals.h"
#include "rheology.h"

// Bulk dynamics of a Carreau-Yasuda fluid: the Newtonian bulk dynamics of HemoFlow
// (GuoExternalForceCompleteRegularizedBGKdynamics, regularized BGK with Guo forcing) with
// the relaxation rate of every cell computed from its own non-equilibrium stress
// (CarreauModel::omegaFromPi) at every collision.
//
// The stored omega (getOmega()) is the value at the infinite-shear viscosity. Palabos uses it
// only where a boundary condition needs one value for all its cells: the finite-difference
// pressure boundary builds the stress of the opening cells from the velocity gradient with it.
template <typename T_, template <typename U> class Descriptor>
class GuoForcedCarreauDynamics : public GuoExternalForceCompleteRegularizedBGKdynamics<T_, Descriptor> {
public:
    typedef GuoExternalForceCompleteRegularizedBGKdynamics<T_, Descriptor> Base;

    GuoForcedCarreauDynamics(T_ omegaInf, CarreauModel const& model_) : Base(omegaInf), model(model_) { }

    GuoForcedCarreauDynamics(HierarchicUnserializer& unserializer) : Base(T_())
    {
        this->unserialize(unserializer);
    }

    virtual GuoForcedCarreauDynamics<T_, Descriptor>* clone() const { return new GuoForcedCarreauDynamics<T_, Descriptor>(*this); }
    virtual int getId() const { return id; }

    // The model parameters travel with the dynamics when Palabos copies cells between blocks.
    virtual void serialize(HierarchicSerializer& serializer) const
    {
        Base::serialize(serializer);
        serializer.addValue(model.nu0);
        serializer.addValue(model.nuInf);
        serializer.addValue(model.lambda);
        serializer.addValue(model.n);
        serializer.addValue(model.a);
    }

    virtual void unserialize(HierarchicUnserializer& unserializer)
    {
        Base::unserialize(unserializer);
        unserializer.readValue(model.nu0);
        unserializer.readValue(model.nuInf);
        unserializer.readValue(model.lambda);
        unserializer.readValue(model.n);
        unserializer.readValue(model.a);
    }

    // Same steps as Base::collide, with the local relaxation rate
    virtual void collide(Cell<T_, Descriptor>& cell, BlockStatistics& statistics)
    {
        T_ rhoBar = this->computeRhoBar(cell);
        Array<T_, Descriptor<T_>::d> u, j;
        this->computeVelocity(cell, u);
        T_ rho = Descriptor<T_>::fullRho(rhoBar);
        for (plint iD = 0; iD < Descriptor<T_>::d; ++iD)
            j[iD] = rho * u[iD];
        collideWithOmega(cell, rhoBar, u, j, statistics);
    }

    // Same steps as Base::collideExternal, with the local relaxation rate
    virtual void collideExternal(Cell<T_, Descriptor>& cell, T_ rhoBar, Array<T_, Descriptor<T_>::d> const& j, T_,
                                 BlockStatistics& statistics)
    {
        Array<T_, Descriptor<T_>::d> u;
        this->computeVelocityExternal(cell, rhoBar, j, u);
        T_ rho = Descriptor<T_>::fullRho(rhoBar);
        Array<T_, Descriptor<T_>::d> newJ;
        for (plint iD = 0; iD < Descriptor<T_>::d; ++iD)
            newJ[iD] = rho * u[iD];
        collideWithOmega(cell, rhoBar, u, newJ, statistics);
    }

    CarreauModel const& getModel() const { return model; }

private:
    void collideWithOmega(Cell<T_, Descriptor>& cell, T_ rhoBar, Array<T_, Descriptor<T_>::d> const& u,
                          Array<T_, Descriptor<T_>::d> const& j, BlockStatistics& statistics)
    {
        T_ rhoBarLb;
        Array<T_, Descriptor<T_>::d> jLb;
        Array<T_, SymmetricTensor<T_, Descriptor>::n> piNeqLb;
        momentTemplates<T_, Descriptor>::compute_rhoBar_j_PiNeq(cell, rhoBarLb, jLb, piNeqLb);
        T_ jSqrLb = VectorTemplate<T_, Descriptor>::normSqr(jLb);

        T_ omega = model.omegaFromPi(SymmetricTensor<T_, Descriptor>::tensorNormSqr(piNeqLb),
                                     Descriptor<T_>::fullRho(rhoBarLb), Descriptor<T_>::cs2);

        // The regularization does not depend on omega (only on the ratio of its two omegas)
        this->regularize(cell, rhoBarLb, jLb, jSqrLb, piNeqLb);
        T_ uSqr = dynamicsTemplates<T_, Descriptor>::complete_bgk_ma2_collision(
            cell, rhoBar, Descriptor<T_>::invRho(rhoBar), j, omega);
        externalForceTemplates<T_, Descriptor>::addGuoForce(cell, u, omega, (T_)1);

        if (cell.takesStatistics())
            gatherStatistics(statistics, rhoBar, uSqr);
    }

    CarreauModel model;   // Lattice units
    static int id;
};

template <typename T_, template <typename U> class Descriptor>
int GuoForcedCarreauDynamics<T_, Descriptor>::id =
    meta::registerGeneralDynamics<T_, Descriptor, GuoForcedCarreauDynamics<T_, Descriptor> >("HemoFlow_GuoForcedCarreau");


// Local relaxation rate of a cell, as GuoForcedCarreauDynamics computes it. Applied to every
// cell (also walls and boundary cells), so that the output fields need no dynamics lookup.
template <typename T_, template <typename U> class Descriptor>
inline T_ carreauCellOmega(CarreauModel const& model, Cell<T_, Descriptor> const& cell, Array<T_, 6>& piNeq)
{
    cell.computePiNeq(piNeq);
    T_ rho = cell.computeDensity();
    if (!(rho > 0.0))
        rho = 1.0;
    return model.omegaFromPi(SymmetricTensor<T_, Descriptor>::tensorNormSqr(piNeq), rho, Descriptor<T_>::cs2);
}

// Viscous stress -(1 - omega/2) Pi^neq with the local relaxation rate (lattice units).
template <typename T_, template <typename U> class Descriptor>
class CarreauStressFunctional3D : public BoxProcessingFunctional3D_LT<T_, Descriptor, T_, 6> {
public:
    explicit CarreauStressFunctional3D(CarreauModel const& model_) : model(model_) { }

    virtual void process(Box3D domain, BlockLattice3D<T_, Descriptor>& lattice, TensorField3D<T_, 6>& stress)
    {
        static_assert(SymmetricTensor<T_, Descriptor>::n == 6, "3D lattice expected");
        Dot3D offset = computeRelativeDisplacement(lattice, stress);
        Array<T_, 6> piNeq;
        for (plint iX = domain.x0; iX <= domain.x1; ++iX)
            for (plint iY = domain.y0; iY <= domain.y1; ++iY)
                for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                    T_ omega = carreauCellOmega(model, lattice.get(iX, iY, iZ), piNeq);
                    stress.get(iX + offset.x, iY + offset.y, iZ + offset.z) = -((T_)1 - (T_)0.5 * omega) * piNeq;
                }
    }

    virtual CarreauStressFunctional3D<T_, Descriptor>* clone() const { return new CarreauStressFunctional3D<T_, Descriptor>(*this); }
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const
    {
        modified[0] = modif::nothing;
        modified[1] = modif::staticVariables;
    }
    virtual BlockDomain::DomainT appliesTo() const { return BlockDomain::bulk; }

private:
    CarreauModel model;
};

// Local kinematic viscosity cs^2 (1/omega - 1/2) (lattice units).
template <typename T_, template <typename U> class Descriptor>
class CarreauViscosityFunctional3D : public BoxProcessingFunctional3D_LS<T_, Descriptor, T_> {
public:
    explicit CarreauViscosityFunctional3D(CarreauModel const& model_) : model(model_) { }

    virtual void process(Box3D domain, BlockLattice3D<T_, Descriptor>& lattice, ScalarField3D<T_>& viscosity)
    {
        Dot3D offset = computeRelativeDisplacement(lattice, viscosity);
        Array<T_, 6> piNeq;
        for (plint iX = domain.x0; iX <= domain.x1; ++iX)
            for (plint iY = domain.y0; iY <= domain.y1; ++iY)
                for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                    T_ omega = carreauCellOmega(model, lattice.get(iX, iY, iZ), piNeq);
                    viscosity.get(iX + offset.x, iY + offset.y, iZ + offset.z) = Descriptor<T_>::cs2 * ((T_)1 / omega - (T_)0.5);
                }
    }

    virtual CarreauViscosityFunctional3D<T_, Descriptor>* clone() const { return new CarreauViscosityFunctional3D<T_, Descriptor>(*this); }
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const
    {
        modified[0] = modif::nothing;
        modified[1] = modif::staticVariables;
    }
    virtual BlockDomain::DomainT appliesTo() const { return BlockDomain::bulk; }

private:
    CarreauModel model;
};

typedef GuoForcedCarreauDynamics<T, DESCRIPTOR> CarreauBackgroundDynamics;

#endif
