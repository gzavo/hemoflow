/* This file is part of the Palabos library.
 *
 * The Palabos softare is developed since 2011 by FlowKit-Numeca Group Sarl
 * (Switzerland) and the University of Geneva (Switzerland), which jointly
 * own the IP rights for most of the code base. Since October 2019, the
 * Palabos project is maintained by the University of Geneva and accepts
 * source code contributions from the community.
 *
 * Contact:
 * Jonas Latt
 * Computer Science Department
 * University of Geneva
 * 7 Route de Drize
 * 1227 Carouge, Switzerland
 * jonas.latt@unige.ch
 *
 * The most recent release of Palabos can be downloaded at
 * <https://palabos.unige.ch/>
 *
 * The library Palabos is free software: you can redistribute it and/or
 * modify it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * The library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/** \file
 * Data processors for data analysis -- header file.
 */

#ifndef ACCELERATED_PROCESSORS_3D_H
#define ACCELERATED_PROCESSORS_3D_H

#include "core/globalDefs.h"
#include "core/array.h"
#include "atomicBlock/dataProcessingFunctional3D.h"
#include "atomicBlock/reductiveDataProcessingFunctional3D.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"
#include "atomicBlock/blockLattice3D.h"
#include "atomicBlock/atomicAcceleratedLattice3D.h"
#include "atomicBlock/dataField3D.h"
#include "multiBlock/multiBlockLattice3D.h"
#include "multiBlock/acceleratedLattice3D.h"
#include "multiBlock/multiDataField3D.h"
#include "multiPhysics/shanChenLattices3D.h"
#include "dataProcessors/dataAnalysisFunctional3D.h"

namespace plb {

/* *************** Density ******************************************* */

template <typename T, template <typename U> class Descriptor>
void computeDensity(
    AcceleratedLattice3D<T, Descriptor> &lattice, MultiScalarField3D<T> &density, Box3D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField3D<T> > computeDensity(
    AcceleratedLattice3D<T, Descriptor> &lattice, Box3D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField3D<T> > computeDensity(
    AcceleratedLattice3D<T, Descriptor> &lattice);

template <typename T, template <typename U> class Descriptor>
class AccDensityFunctional3D : public BoxProcessingFunctional3D_AS<T, Descriptor, T> {
public:
    virtual void process(
        Box3D domain, AtomicAcceleratedLattice3D<T, Descriptor> &lattice,
        ScalarField3D<T> &scalarField);
    virtual AccDensityFunctional3D<T, Descriptor> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
};

/* *************** Velocity Norm ************************************* */

template <typename T, template <typename U> class Descriptor>
void computeVelocityNorm(
    AcceleratedLattice3D<T, Descriptor> &lattice, MultiScalarField3D<T> &velocityNorm);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField3D<T> > computeVelocityNorm(
    AcceleratedLattice3D<T, Descriptor> &lattice);

template <typename T, template <typename U> class Descriptor>
class AccelVelocityNormFunctional3D : public BoxProcessingFunctional3D_AS<T, Descriptor, T> {
public:
    virtual void process(
        Box3D domain, AtomicAcceleratedLattice3D<T, Descriptor> &lattice,
        ScalarField3D<T> &scalarField);
    virtual AccelVelocityNormFunctional3D<T, Descriptor> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

/* *************** Velocity Component ************************************* */

template <typename T, template <typename U> class Descriptor>
void computeVelocityComponent(
    AcceleratedLattice3D<T, Descriptor> &lattice, MultiScalarField3D<T> &velocityComponent,
    int iComponent);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField3D<T> > computeVelocityComponent(
    AcceleratedLattice3D<T, Descriptor> &lattice, int iComponent);

template <typename T, template <typename U> class Descriptor>
class AccelVelocityComponentFunctional3D : public BoxProcessingFunctional3D_AS<T, Descriptor, T> {
public:
    AccelVelocityComponentFunctional3D(int iComponent_) : iComponent(iComponent_) { }
    virtual void process(
        Box3D domain, AtomicAcceleratedLattice3D<T, Descriptor> &lattice,
        ScalarField3D<T> &scalarField);
    virtual AccelVelocityComponentFunctional3D<T, Descriptor> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;

private:
    int iComponent;
};

/* *************** Velocity ****************************************** */

template <typename T, template <typename U> class Descriptor>
void computeVelocity(
    AcceleratedLattice3D<T, Descriptor> &lattice,
    MultiTensorField3D<T, Descriptor<T>::d> &velocity);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiTensorField3D<T, Descriptor<T>::d> > computeVelocity(
    AcceleratedLattice3D<T, Descriptor> &lattice);

template <typename T, template <typename U> class Descriptor>
class AccelVelocityFunctional3D :
    public BoxProcessingFunctional3D_AT<T, Descriptor, T, Descriptor<T>::d> {
public:
    virtual void process(
        Box3D domain, AtomicAcceleratedLattice3D<T, Descriptor> &lattice,
        TensorField3D<T, Descriptor<T>::d> &tensorField);
    virtual AccelVelocityFunctional3D<T, Descriptor> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

/* *************** Kinetic Energy ************************************ */

template <typename T, template <typename U> class Descriptor>
void computeKineticEnergy(
    AcceleratedLattice3D<T, Descriptor> &lattice, MultiScalarField3D<T> &energy, Box3D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField3D<T> > computeKineticEnergy(
    AcceleratedLattice3D<T, Descriptor> &lattice, Box3D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField3D<T> > computeKineticEnergy(
    AcceleratedLattice3D<T, Descriptor> &lattice);

template <typename T, template <typename U> class Descriptor>
class AccKineticEnergyFunctional3D : public BoxProcessingFunctional3D_AS<T, Descriptor, T> {
public:
    virtual void process(
        Box3D domain, AtomicAcceleratedLattice3D<T, Descriptor> &lattice,
        ScalarField3D<T> &scalarField);
    virtual AccKineticEnergyFunctional3D<T, Descriptor> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

/* *************** Initialization ************************************ */

template <typename T, template <class U> class Descriptor, class RhoUFunction>
void initializeAtEquilibrium(
    AcceleratedLattice3D<T, Descriptor> &lattice, Box3D domain, RhoUFunction f);

template <typename T, template <class U> class Descriptor, class RhoUFunction>
void initializeAtEquilibrium_o2(
    AcceleratedLattice3D<T, Descriptor> &lattice, Box3D domain, RhoUFunction f);

template <typename T, template <class U> class Descriptor, class RhoUFunction>
class AccCustomEquilibriumFunctional3D : public BoxProcessingFunctional3D_A<T, Descriptor> {
public:
    AccCustomEquilibriumFunctional3D(RhoUFunction const &f_) : f(f_) { }
    virtual void process(Box3D domain, AtomicAcceleratedLattice3D<T, Descriptor> &lattice);
    virtual AccCustomEquilibriumFunctional3D<T, Descriptor, RhoUFunction> *clone() const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    RhoUFunction f;
};

template <typename T, template <class U> class Descriptor, class RhoUFunction>
class AccCustomEquilibriumFunctional_o2_3D : public BoxProcessingFunctional3D_A<T, Descriptor> {
public:
    AccCustomEquilibriumFunctional_o2_3D(RhoUFunction const &f_) : f(f_) { }
    virtual void process(Box3D domain, AtomicAcceleratedLattice3D<T, Descriptor> &lattice);
    virtual AccCustomEquilibriumFunctional_o2_3D<T, Descriptor, RhoUFunction> *clone() const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    RhoUFunction f;
};

/* *************** Accelerated Scalar-Field Manipulation ************* */

template <typename T, int nDim>
class AccComputeNormSqrFunctional3D : public BoxProcessingFunctional3D_ST<T, T, nDim> {
public:
    virtual void process(
        Box3D domain, ScalarField3D<T> &scalarField, TensorField3D<T, nDim> &tensorField);
    virtual AccComputeNormSqrFunctional3D<T, nDim> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template <typename T, int nDim>
std::unique_ptr<MultiScalarField3D<T> > computeNormSqrAcc(
    MultiTensorField3D<T, nDim> &tensorField, Box3D domain);

template <typename T, int nDim>
std::unique_ptr<MultiScalarField3D<T> > computeNormSqrAcc(MultiTensorField3D<T, nDim> &tensorField);

template <typename T>
class AccScalarSumFunctional3D : public ReductiveBoxProcessingFunctional3D_S<T> {
public:
    AccScalarSumFunctional3D();
    virtual void process(Box3D domain, ScalarField3D<T> &scalarField);
    virtual AccScalarSumFunctional3D<T> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        modified[0] = modif::nothing;
    }
    T getSumScalar() const;

private:
    plint sumScalarId;
};

template <typename T>
T computeAverageAcc(MultiScalarField3D<T> &scalarField, Box3D domain);

template <typename T>
T computeAverageAcc(MultiScalarField3D<T> &scalarField);

/* *************** Accelerated Tensor-Field Manipulation ************* */

/** Attention: No matter what the type of T is (even if it is an integer type),
 *    the sum is computed in double-precision floating point numbers, and
 *    converted to T at the end (and rounded, if T is an integer).
 **/
template <typename T, int nDim>
class AccTensorSumFunctional3D : public ReductiveBoxProcessingFunctional3D_T<T, nDim> {
public:
    AccTensorSumFunctional3D();
    virtual void process(Box3D domain, TensorField3D<T, nDim> &tensorField);
    virtual AccTensorSumFunctional3D<T, nDim> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        modified[0] = modif::nothing;
    }
    Array<T, nDim> getSumTensor() const;

private:
    Array<plint, nDim> sumTensorId;
};

template <typename T, int nDim>
Array<T, nDim> computeAverageAcc(MultiTensorField3D<T, nDim> &tensorField, Box3D domain);

template <typename T, int nDim>
Array<T, nDim> computeAverageAcc(MultiTensorField3D<T, nDim> &tensorField);

/// Use of a 8 points stencil for computation of FD gradients for the vorticity.
template <typename T, int nDim>
class AccBulkVorticityOrderEightFunctional3D :
    public BoxProcessingFunctional3D_TT<T, nDim, T, nDim> {
public:
    virtual void process(
        Box3D domain, TensorField3D<T, nDim> &velocity, TensorField3D<T, nDim> &vorticity);
    virtual AccBulkVorticityOrderEightFunctional3D<T, nDim> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template <typename T>
std::unique_ptr<MultiTensorField3D<T, 3> > computeBulkVorticityOrderEightAcc(
    MultiTensorField3D<T, 3> &velocity);

/* *************** Multi-phase *************************************** */

/// Helper functions with full-lattice access
template <typename T, template <typename U> class Descriptor>
struct accelMultiPhaseTemplates3D {
    static void shanChenInteraction(
        AtomicAcceleratedLattice3D<T, Descriptor> &lattice,
        Array<T, Descriptor<T>::d> &rhoContribution, plint iX, plint iY, plint iZ)
    {
        enum { densityOffset = Descriptor<T>::ExternalField::densityBeginsAt };

        rhoContribution.resetToZero();
        for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
            plint nextX = iX + Descriptor<T>::c_gpu(iPop, 0);
            plint nextY = iY + Descriptor<T>::c_gpu(iPop, 1);
            plint nextZ = iZ + Descriptor<T>::c_gpu(iPop, 2);
            T rho = lattice.pullExt(nextX, nextY, nextZ, densityOffset);
            for (int iD = 0; iD < Descriptor<T>::d; ++iD) {
                rhoContribution[iD] +=
                    Descriptor<T>::t_gpu(iPop) * rho * Descriptor<T>::c_gpu(iPop, iD);
            }
        }
    }

    static void shanChenInteraction(
        AtomicAcceleratedLattice3D<T, Descriptor> &lattice, ScalarField3D<T> &rhoBar,
        Array<T, Descriptor<T>::d> &rhoContribution, plint iX, plint iY, plint iZ)
    {
        Dot3D ofs = computeRelativeDisplacement(lattice, rhoBar);
        rhoContribution.resetToZero();
        for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
            plint nextX = iX + Descriptor<T>::c_gpu(iPop, 0);
            plint nextY = iY + Descriptor<T>::c_gpu(iPop, 1);
            plint nextZ = iZ + Descriptor<T>::c_gpu(iPop, 2);
            T rho = Descriptor<T>::fullRho(rhoBar.get(nextX + ofs.x, nextY + ofs.y, nextZ + ofs.z));
            for (int iD = 0; iD < Descriptor<T>::d; ++iD) {
                rhoContribution[iD] +=
                    Descriptor<T>::t_gpu(iPop) * rho * Descriptor<T>::c_gpu(iPop, iD);
            }
        }
    }
};

template <typename T>
struct accelMultiPhaseTemplates3D<T, descriptors::ForcedShanChenD3Q19Descriptor> {
    typedef descriptors::ForcedShanChenD3Q19Descriptor<T> D;

    static void shanChenInteraction(
        AtomicAcceleratedLattice3D<T, descriptors::ForcedShanChenD3Q19Descriptor> &lattice,
        Array<T, D::d> &rhoContribution, plint iX, plint iY, plint iZ)
    {
        enum { densityOffset = D::ExternalField::densityBeginsAt };

        T rho;
        rho = lattice.pullExt(iX - 1, iY, iZ, densityOffset);
        rhoContribution[0] = -D::t_gpu(1) * rho;
        rho = lattice.pullExt(iX, iY - 1, iZ, densityOffset);
        rhoContribution[1] = -D::t_gpu(2) * rho;
        rho = lattice.pullExt(iX, iY, iZ - 1, densityOffset);
        rhoContribution[2] = -D::t_gpu(3) * rho;
        rho = lattice.pullExt(iX - 1, iY - 1, iZ, densityOffset);
        rhoContribution[0] -= D::t_gpu(4) * rho;
        rhoContribution[1] -= D::t_gpu(4) * rho;
        rho = lattice.pullExt(iX - 1, iY + 1, iZ, densityOffset);
        rhoContribution[0] -= D::t_gpu(5) * rho;
        rhoContribution[1] += D::t_gpu(5) * rho;
        rho = lattice.pullExt(iX - 1, iY, iZ - 1, densityOffset);
        rhoContribution[0] -= D::t_gpu(6) * rho;
        rhoContribution[2] -= D::t_gpu(6) * rho;
        rho = lattice.pullExt(iX - 1, iY, iZ + 1, densityOffset);
        rhoContribution[0] -= D::t_gpu(7) * rho;
        rhoContribution[2] += D::t_gpu(7) * rho;
        rho = lattice.pullExt(iX, iY - 1, iZ - 1, densityOffset);
        rhoContribution[1] -= D::t_gpu(8) * rho;
        rhoContribution[2] -= D::t_gpu(8) * rho;
        rho = lattice.pullExt(iX, iY - 1, iZ + 1, densityOffset);
        rhoContribution[1] -= D::t_gpu(9) * rho;
        rhoContribution[2] += D::t_gpu(9) * rho;

        rho = lattice.pullExt(iX + 1, iY, iZ, densityOffset);
        rhoContribution[0] += D::t_gpu(10) * rho;
        rho = lattice.pullExt(iX, iY + 1, iZ, densityOffset);
        rhoContribution[1] += D::t_gpu(11) * rho;
        rho = lattice.pullExt(iX, iY, iZ + 1, densityOffset);
        rhoContribution[2] += D::t_gpu(12) * rho;
        rho = lattice.pullExt(iX + 1, iY + 1, iZ, densityOffset);
        rhoContribution[0] += D::t_gpu(13) * rho;
        rhoContribution[1] += D::t_gpu(13) * rho;
        rho = lattice.pullExt(iX + 1, iY - 1, iZ, densityOffset);
        rhoContribution[0] += D::t_gpu(14) * rho;
        rhoContribution[1] -= D::t_gpu(14) * rho;
        rho = lattice.pullExt(iX + 1, iY, iZ + 1, densityOffset);
        rhoContribution[0] += D::t_gpu(15) * rho;
        rhoContribution[2] += D::t_gpu(15) * rho;
        rho = lattice.pullExt(iX + 1, iY, iZ - 1, densityOffset);
        rhoContribution[0] += D::t_gpu(16) * rho;
        rhoContribution[2] -= D::t_gpu(16) * rho;
        rho = lattice.pullExt(iX, iY + 1, iZ + 1, densityOffset);
        rhoContribution[1] += D::t_gpu(17) * rho;
        rhoContribution[2] += D::t_gpu(17) * rho;
        rho = lattice.pullExt(iX, iY + 1, iZ - 1, densityOffset);
        rhoContribution[1] += D::t_gpu(18) * rho;
        rhoContribution[2] -= D::t_gpu(18) * rho;
    }
};

/// Shan-Chen coupling for multi-component flow with or without external force
template <typename T, template <typename U> class Descriptor, int numSpecies>
class ShanChenMultiComponentAccelerated3D :
    public AcceleratedBoxProcessingFunctional3D<T, Descriptor> {
public:
    ShanChenMultiComponentAccelerated3D(std::vector<std::vector<T> > const &speciesG_);
    /// With these constructors, the values of the relaxation parameters omega are
    ///   taken to be species-dependent, but not space- or time-dependent. Their
    ///   value is imposed in the constructor.
    ShanChenMultiComponentAccelerated3D(T G_, std::vector<T> const &imposedOmega_);
    ShanChenMultiComponentAccelerated3D(
        std::vector<std::vector<T> > const &speciesG_, std::vector<T> const &imposedOmega_);
    virtual void process(
        Box3D domain, std::vector<AtomicAcceleratedLattice3D<T, Descriptor> *> lattices);
    virtual ShanChenMultiComponentAccelerated3D<T, Descriptor, numSpecies> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    T G;
    std::vector<T> speciesG;
    std::vector<T> imposedOmega;
};

template <typename T, template <typename U> class Descriptor>
void showTemplateArguments(MultiBlockLattice3D<T, Descriptor> &lattice);

}  // namespace plb

#endif  // ACCELERATED_PROCESSORS_3D_H
