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

#ifndef ACCELERATED_PROCESSORS_3D_HH
#define ACCELERATED_PROCESSORS_3D_HH

#include "core/globalDefs.h"
#include "acceleratedLattice/acceleratedProcessors3D.h"
#include "core/plbDebug.h"
#include "core/util.h"
#include "core/blockStatistics.h"
#include "latticeBoltzmann/momentTemplates.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"
#include "finiteDifference/fdStencils1D.h"
#include "atomicBlock/atomicBlock3D.h"
#include "atomicBlock/blockLattice3D.h"
#include "atomicBlock/atomicAcceleratedLattice3D.h"
#include "atomicBlock/dataField3D.h"
#include "atomicBlock/reductiveDataProcessorWrapper3D.h"
#include "atomicBlock/dataProcessorWrapper3D.h"
#include "multiBlock/reductiveMultiDataProcessorWrapper3D.h"
#include "multiBlock/multiDataProcessorWrapper3D.h"
#include "multiBlock/multiBlockGenerator3D.h"
#include <cmath>
#include <limits>
#include <map>
#include <string>

namespace plb {

/* *************** Density ******************************************* */

template <typename T, template <typename U> class Descriptor>
void computeDensity(
    AcceleratedLattice3D<T, Descriptor> &lattice, MultiScalarField3D<T> &density, Box3D domain)
{
    applyProcessingFunctional(new AccDensityFunctional3D<T, Descriptor>, domain, lattice, density);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField3D<T>> computeDensity(
    AcceleratedLattice3D<T, Descriptor> &lattice, Box3D domain)
{
    std::unique_ptr<MultiScalarField3D<T>> density = generateMultiScalarField<T>(lattice, domain);
    computeDensity(lattice, *density, domain);
    return density;
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField3D<T>> computeDensity(AcceleratedLattice3D<T, Descriptor> &lattice)
{
    return computeDensity(lattice, lattice.getBoundingBox());
}

template <typename T, template <typename U> class Descriptor>
void AccDensityFunctional3D<T, Descriptor>::process(
    Box3D domain, AtomicAcceleratedLattice3D<T, Descriptor> &lattice, ScalarField3D<T> &scalarField)
{
    Dot3D offset = computeRelativeDisplacement(lattice, scalarField);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                Cell<T, Descriptor> cell;
                lattice.reconstructCell(iX, iY, iZ, cell);
                scalarField.get(iX + offset.x, iY + offset.y, iZ + offset.z) =
                    cell.computeDensity();
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
AccDensityFunctional3D<T, Descriptor> *AccDensityFunctional3D<T, Descriptor>::clone() const
{
    return new AccDensityFunctional3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void AccDensityFunctional3D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

/* *************** Velocity Norm ************************************* */

template <typename T, template <typename U> class Descriptor>
void computeVelocityNorm(
    AcceleratedLattice3D<T, Descriptor> &lattice, MultiScalarField3D<T> &velocityNorm, Box3D domain)
{
    applyProcessingFunctional(
        new AccelVelocityNormFunctional3D<T, Descriptor>, domain, lattice, velocityNorm);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField3D<T>> computeVelocityNorm(
    AcceleratedLattice3D<T, Descriptor> &lattice, Box3D domain)
{
    std::unique_ptr<MultiScalarField3D<T>> velocityNorm =
        generateMultiScalarField<T>(lattice, domain);

    computeVelocityNorm(lattice, *velocityNorm, domain);
    return velocityNorm;
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField3D<T>> computeVelocityNorm(
    AcceleratedLattice3D<T, Descriptor> &lattice)
{
    return computeVelocityNorm(lattice, lattice.getBoundingBox());
}

template <typename T, template <typename U> class Descriptor>
void AccelVelocityNormFunctional3D<T, Descriptor>::process(
    Box3D domain, AtomicAcceleratedLattice3D<T, Descriptor> &lattice, ScalarField3D<T> &scalarField)
{
    Dot3D offset = computeRelativeDisplacement(lattice, scalarField);
    lattice.for_each(
        domain, [&lattice, &scalarField, offset](
                    plint i, plint iX, plint iY, plint iZ, int collisionModel) {
            Array<T, Descriptor<T>::q> f;
            lattice.pullPop(i, f);
            Array<T, Descriptor<T>::d> velocity;
            momentTemplatesImpl<T, typename Descriptor<T>::BaseDescriptor>::compute_uLb(
                f, velocity);
            // The type cast converts the result of normSqr to type U in case T is of type
            // Complex<U>. Otherwise, the call to std::sqrt would fail, because std::sqrt is
            // overloaded, but not for Palabos' Complex type.
            plint iX2 = iX + offset.x;
            plint iY2 = iY + offset.y;
            plint iZ2 = iZ + offset.z;
            scalarField.get(iX2, iY2, iZ2) = std::sqrt(
                (typename PlbTraits<T>::BaseType)VectorTemplate<T, Descriptor>::normSqr(velocity));
        });
}

template <typename T, template <typename U> class Descriptor>
AccelVelocityNormFunctional3D<T, Descriptor> *AccelVelocityNormFunctional3D<T, Descriptor>::clone()
    const
{
    return new AccelVelocityNormFunctional3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void AccelVelocityNormFunctional3D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T, template <typename U> class Descriptor>
BlockDomain::DomainT AccelVelocityNormFunctional3D<T, Descriptor>::appliesTo() const
{
    return BlockDomain::bulk;
}

/* *************** Velocity Component ************************************* */

template <typename T, template <typename U> class Descriptor>
void computeVelocityComponent(
    AcceleratedLattice3D<T, Descriptor> &lattice, MultiScalarField3D<T> &velocityComponent,
    int iComponent, Box3D domain)
{
    applyProcessingFunctional(
        new AccelVelocityComponentFunctional3D<T, Descriptor>(iComponent), domain, lattice,
        velocityComponent);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField3D<T>> computeVelocityComponent(
    AcceleratedLattice3D<T, Descriptor> &lattice, int iComponent, Box3D domain)
{
    std::unique_ptr<MultiScalarField3D<T>> velocityComponent =
        generateMultiScalarField<T>(lattice, domain);

    computeVelocityComponent(lattice, *velocityComponent, iComponent, domain);
    return velocityComponent;
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField3D<T>> computeVelocityComponent(
    AcceleratedLattice3D<T, Descriptor> &lattice, int iComponent)
{
    return computeVelocityComponent(lattice, iComponent, lattice.getBoundingBox());
}

template <typename T, template <typename U> class Descriptor>
void AccelVelocityComponentFunctional3D<T, Descriptor>::process(
    Box3D domain, AtomicAcceleratedLattice3D<T, Descriptor> &lattice, ScalarField3D<T> &scalarField)
{
    Dot3D offset = computeRelativeDisplacement(lattice, scalarField);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                Cell<T, Descriptor> cell;
                lattice.reconstructCell(iX, iY, iZ, cell);
                Array<T, Descriptor<T>::d> velocity;
                cell.computeVelocity(velocity);
                // The type cast converts the result of normSqr to type U in case T is of type
                // Complex<U>. Otherwise, the call to std::sqrt would fail, because std::sqrt is
                // overloaded, but not for Palabos' Complex type.
                scalarField.get(iX + offset.x, iY + offset.y, iZ + offset.z) = velocity[iComponent];
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
AccelVelocityComponentFunctional3D<T, Descriptor>
    *AccelVelocityComponentFunctional3D<T, Descriptor>::clone() const
{
    return new AccelVelocityComponentFunctional3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void AccelVelocityComponentFunctional3D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T, template <typename U> class Descriptor>
BlockDomain::DomainT AccelVelocityComponentFunctional3D<T, Descriptor>::appliesTo() const
{
    return BlockDomain::bulk;
}

/* *************** Velocity ****************************************** */

template <typename T, template <typename U> class Descriptor>
void computeVelocity(
    AcceleratedLattice3D<T, Descriptor> &lattice, MultiTensorField3D<T, Descriptor<T>::d> &velocity,
    Box3D domain)
{
    applyProcessingFunctional(
        new AccelVelocityFunctional3D<T, Descriptor>, domain, lattice, velocity);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiTensorField3D<T, Descriptor<T>::d>> computeVelocity(
    AcceleratedLattice3D<T, Descriptor> &lattice, Box3D domain)
{
    std::unique_ptr<MultiTensorField3D<T, Descriptor<T>::d>> velocity =
        generateMultiTensorField<T, Descriptor<T>::d>(lattice, domain);

    computeVelocity(lattice, *velocity, domain);
    return velocity;
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiTensorField3D<T, Descriptor<T>::d>> computeVelocity(
    AcceleratedLattice3D<T, Descriptor> &lattice)
{
    return computeVelocity(lattice, lattice.getBoundingBox());
}

template <typename T, template <typename U> class Descriptor>
void AccelVelocityFunctional3D<T, Descriptor>::process(
    Box3D domain, AtomicAcceleratedLattice3D<T, Descriptor> &lattice,
    TensorField3D<T, Descriptor<T>::d> &tensorField)
{
    Dot3D offset = computeRelativeDisplacement(lattice, tensorField);
    plint ny2 = tensorField.getNy();
    plint nz2 = tensorField.getNz();
    lattice.for_each(
        domain, [&lattice, &tensorField, ny2, nz2, offset](
                    plint i, plint iX, plint iY, plint iZ, int collisionModel) {
            plint iX2 = iX + offset.x;
            plint iY2 = iY + offset.y;
            plint iZ2 = iZ + offset.z;
            plint i2 = iZ2 + nz2 * (iY2 + ny2 * iX2);
            Array<T, Descriptor<T>::q> f;
            lattice.pullPop(i, f);
            Array<T, 3> u;
            momentTemplatesImpl<T, typename Descriptor<T>::BaseDescriptor>::compute_uLb(
                f, tensorField[i2]);
        });
}

template <typename T, template <typename U> class Descriptor>
AccelVelocityFunctional3D<T, Descriptor> *AccelVelocityFunctional3D<T, Descriptor>::clone() const
{
    return new AccelVelocityFunctional3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void AccelVelocityFunctional3D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T, template <typename U> class Descriptor>
BlockDomain::DomainT AccelVelocityFunctional3D<T, Descriptor>::appliesTo() const
{
    return BlockDomain::bulk;
}

/* *************** Kinetic Energy ************************************ */

template <typename T, template <typename U> class Descriptor>
void computeKineticEnergy(
    AcceleratedLattice3D<T, Descriptor> &lattice, MultiScalarField3D<T> &energy, Box3D domain)
{
    applyProcessingFunctional(
        new AccKineticEnergyFunctional3D<T, Descriptor>, domain, lattice, energy);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField3D<T>> computeKineticEnergy(
    AcceleratedLattice3D<T, Descriptor> &lattice, Box3D domain)
{
    std::unique_ptr<MultiScalarField3D<T>> energy = generateMultiScalarField<T>(lattice, domain);

    computeKineticEnergy(lattice, *energy, domain);

    return energy;
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField3D<T>> computeKineticEnergy(
    AcceleratedLattice3D<T, Descriptor> &lattice)
{
    return computeKineticEnergy(lattice, lattice.getBoundingBox());
}

template <typename T, template <typename U> class Descriptor>
void AccKineticEnergyFunctional3D<T, Descriptor>::process(
    Box3D domain, AtomicAcceleratedLattice3D<T, Descriptor> &lattice, ScalarField3D<T> &scalarField)
{
    Dot3D offset = computeRelativeDisplacement(lattice, scalarField);
    lattice.for_each(
        domain, [&lattice, &scalarField, offset](
                    plint i, plint iX, plint iY, plint iZ, int collisionModel) {
            T energy = T();
            if (collisionModel != CollisionModel::BounceBack) {
                Array<T, Descriptor<T>::d> u;
                Array<T, Descriptor<T>::q> f;
                lattice.pullPop(i, f);
                momentTemplatesImpl<T, typename Descriptor<T>::BaseDescriptor>::compute_uLb(f, u);
                energy = VectorTemplate<T, Descriptor>::normSqr(u) * (T)0.5;
            }
            scalarField.ScalarField3D<T>::get(iX + offset.x, iY + offset.y, iZ + offset.z) = energy;
        });
}

template <typename T, template <typename U> class Descriptor>
AccKineticEnergyFunctional3D<T, Descriptor> *AccKineticEnergyFunctional3D<T, Descriptor>::clone()
    const
{
    return new AccKineticEnergyFunctional3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void AccKineticEnergyFunctional3D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T, template <typename U> class Descriptor>
BlockDomain::DomainT AccKineticEnergyFunctional3D<T, Descriptor>::appliesTo() const
{
    return BlockDomain::bulk;
}

/* *************** InitializeAtEquilibrium ************************************ */

template <typename T, template <class U> class Descriptor, class RhoUFunction>
void initializeAtEquilibrium(
    AcceleratedLattice3D<T, Descriptor> &lattice, Box3D domain, RhoUFunction f)
{
    applyProcessingFunctional(
        new AccCustomEquilibriumFunctional3D<T, Descriptor, RhoUFunction>(f), domain, lattice);
}

template <typename T, template <class U> class Descriptor, class RhoUFunction>
void initializeAtEquilibrium_o2(
    AcceleratedLattice3D<T, Descriptor> &lattice, Box3D domain, RhoUFunction f)
{
    applyProcessingFunctional(
        new AccCustomEquilibriumFunctional_o2_3D<T, Descriptor, RhoUFunction>(f), domain, lattice);
}

template <typename T, template <class U> class Descriptor, class RhoUFunction>
void AccCustomEquilibriumFunctional3D<T, Descriptor, RhoUFunction>::process(
    Box3D domain, AtomicAcceleratedLattice3D<T, Descriptor> &lattice)
{
    Dot3D relativeOffset = lattice.getLocation();
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                T rho;
                Array<T, 3> j;
                f(iX + relativeOffset.x, iY + relativeOffset.y, iZ + relativeOffset.z, rho, j);
                j *= rho;
                T jSqr = VectorTemplate<T, Descriptor>::normSqr(j);
                Cell<T, Descriptor> cell(&lattice.getDynamics(iX, iY, iZ));
                for (int iPop = 0; iPop < Descriptor<T>::numPop; ++iPop) {
                    cell[iPop] = cell.computeEquilibrium(iPop, Descriptor<T>::rhoBar(rho), j, jSqr);
                }
                lattice.pushStatic(iX, iY, iZ, cell);
            }
        }
    }
}

template <typename T, template <class U> class Descriptor, class RhoUFunction>
AccCustomEquilibriumFunctional3D<T, Descriptor, RhoUFunction>
    *AccCustomEquilibriumFunctional3D<T, Descriptor, RhoUFunction>::clone() const
{
    return new AccCustomEquilibriumFunctional3D<T, Descriptor, RhoUFunction>(*this);
}

template <typename T, template <class U> class Descriptor, class RhoUFunction>
void AccCustomEquilibriumFunctional3D<T, Descriptor, RhoUFunction>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;
}

template <typename T, template <class U> class Descriptor, class RhoUFunction>
BlockDomain::DomainT AccCustomEquilibriumFunctional3D<T, Descriptor, RhoUFunction>::appliesTo()
    const
{
    return BlockDomain::bulk;
}

template <typename T, template <class U> class Descriptor, class RhoUFunction>
void AccCustomEquilibriumFunctional_o2_3D<T, Descriptor, RhoUFunction>::process(
    Box3D domain, AtomicAcceleratedLattice3D<T, Descriptor> &lattice)
{
    Dot3D relativeOffset = lattice.getLocation();
    lattice.for_each(
        domain, [relativeOffset, &lattice, f = f](
                    plint i, plint iX, plint iY, plint iZ, int collisionModel) {
            T rho;
            Array<T, 3> j;
            f(iX + relativeOffset.x, iY + relativeOffset.y, iZ + relativeOffset.z, rho, j);
            j *= rho;
            T jSqr = VectorTemplate<T, Descriptor>::normSqr(j);
            T rhoBar = Descriptor<T>::rhoBar(rho);
            T invRho = Descriptor<T>::invRho(rhoBar);
            Array<T, Descriptor<T>::q> fPop;
            for (int iPop = 0; iPop < Descriptor<T>::numPop; ++iPop) {
                fPop[iPop] = Descriptor<T>::bgk_ma2_equilibrium(iPop, rhoBar, invRho, j, jSqr);
            }
            lattice.pushPop(i, fPop);
        });
}

template <typename T, template <class U> class Descriptor, class RhoUFunction>
AccCustomEquilibriumFunctional_o2_3D<T, Descriptor, RhoUFunction>
    *AccCustomEquilibriumFunctional_o2_3D<T, Descriptor, RhoUFunction>::clone() const
{
    return new AccCustomEquilibriumFunctional_o2_3D<T, Descriptor, RhoUFunction>(*this);
}

template <typename T, template <class U> class Descriptor, class RhoUFunction>
void AccCustomEquilibriumFunctional_o2_3D<T, Descriptor, RhoUFunction>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;
}

template <typename T, template <class U> class Descriptor, class RhoUFunction>
BlockDomain::DomainT AccCustomEquilibriumFunctional_o2_3D<T, Descriptor, RhoUFunction>::appliesTo()
    const
{
    return BlockDomain::bulk;
}

template <typename T, int nDim>
void AccComputeNormSqrFunctional3D<T, nDim>::process(
    Box3D domain, ScalarField3D<T> &scalarField, TensorField3D<T, nDim> &tensorField)
{
    Dot3D offset = computeRelativeDisplacement(scalarField, tensorField);
    scalarField.for_each(
        domain, [&scalarField, &tensorField, offset](
                    plint i, plint iX, plint iY, plint iZ, T const &scalar) {
            scalarField.get(iX, iY, iZ) = VectorTemplateImpl<T, nDim>::normSqr(
                tensorField.get(iX + offset.x, iY + offset.y, iZ + offset.z));
        });
}

template <typename T, int nDim>
AccComputeNormSqrFunctional3D<T, nDim> *AccComputeNormSqrFunctional3D<T, nDim>::clone() const
{
    return new AccComputeNormSqrFunctional3D<T, nDim>(*this);
}

template <typename T, int nDim>
void AccComputeNormSqrFunctional3D<T, nDim>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
}

template <typename T, int nDim>
BlockDomain::DomainT AccComputeNormSqrFunctional3D<T, nDim>::appliesTo() const
{
    return BlockDomain::bulk;
}

template <typename T, int nDim>
void computeNormSqrAcc(MultiTensorField3D<T, nDim> &tensorField, MultiScalarField3D<T> &normSqr)
{
    applyProcessingFunctional(
        new AccComputeNormSqrFunctional3D<T, nDim>, tensorField.getBoundingBox(), normSqr,
        tensorField);
}

template <typename T, int nDim>
std::unique_ptr<MultiScalarField3D<T>> computeNormSqrAcc(
    MultiTensorField3D<T, nDim> &tensorField, Box3D domain)
{
    std::unique_ptr<MultiScalarField3D<T>> normSqr =
        generateMultiScalarField<T>(tensorField, domain);

    applyProcessingFunctional(
        new AccComputeNormSqrFunctional3D<T, nDim>, domain, *normSqr, tensorField);

    return normSqr;
}

template <typename T, int nDim>
std::unique_ptr<MultiScalarField3D<T>> computeNormSqrAcc(MultiTensorField3D<T, nDim> &tensorField)
{
    return computeNormSqrAcc(tensorField, tensorField.getBoundingBox());
}

template <typename T>
AccScalarSumFunctional3D<T>::AccScalarSumFunctional3D() :
    sumScalarId(this->getStatistics().subscribeSum())
{ }

template <typename T>
void AccScalarSumFunctional3D<T>::process(Box3D domain, ScalarField3D<T> &scalarField)
{
    T sumval = scalarField.transform_reduce(
        domain, (T)0, std::plus<T>(),
        [&scalarField](plint i, plint iX, plint iY, plint iZ, T value) {
            return value;
        });
    this->getStatistics().gatherSum(sumScalarId, (double)sumval);
}

template <typename T>
AccScalarSumFunctional3D<T> *AccScalarSumFunctional3D<T>::clone() const
{
    return new AccScalarSumFunctional3D<T>(*this);
}

template <typename T>
T AccScalarSumFunctional3D<T>::getSumScalar() const
{
    double doubleSum = this->getStatistics().getSum(sumScalarId);
    // The sum is internally computed on floating-point values. If T is
    //   integer, the value must be rounded at the end.
    if (std::numeric_limits<T>::is_integer) {
        return (T)util::roundToInt(doubleSum);
    }
    return (T)doubleSum;
}

template <typename T>
T computeAverageAcc(MultiScalarField3D<T> &scalarField, Box3D domain)
{
    AccScalarSumFunctional3D<T> functional;
    applyProcessingFunctional(functional, domain, scalarField);
    return functional.getSumScalar() / (T)domain.nCells();
}

template <typename T>
T computeAverageAcc(MultiScalarField3D<T> &scalarField)
{
    return computeAverageAcc(scalarField, scalarField.getBoundingBox());
}

template <typename T, int nDim>
AccTensorSumFunctional3D<T, nDim>::AccTensorSumFunctional3D()
{
    for (plint i = 0; i < nDim; i++) {
        sumTensorId[i] = this->getStatistics().subscribeSum();
    }
}

template <typename T, int nDim>
void AccTensorSumFunctional3D<T, nDim>::process(Box3D domain, TensorField3D<T, nDim> &tensorField)
{
    Array<T, nDim> sumval = tensorField.transform_reduce(
        domain, Array<T, nDim>::zero(), std::plus<Array<T, nDim>>(),
        [&tensorField](plint i, plint iX, plint iY, plint iZ, Array<T, nDim> value) {
            return value;
        });
    for (plint i = 0; i < nDim; i++) {
        this->getStatistics().gatherSum(sumTensorId[i], (double)sumval[i]);
    }
}

template <typename T, int nDim>
AccTensorSumFunctional3D<T, nDim> *AccTensorSumFunctional3D<T, nDim>::clone() const
{
    return new AccTensorSumFunctional3D<T, nDim>(*this);
}

template <typename T, int nDim>
Array<T, nDim> AccTensorSumFunctional3D<T, nDim>::getSumTensor() const
{
    Array<T, nDim> sum;
    for (plint i = 0; i < nDim; i++) {
        double doubleSum = this->getStatistics().getSum(sumTensorId[i]);
        // The sum is internally computed on floating-point values. If T is
        //   integer, the value must be rounded at the end.
        if (std::numeric_limits<T>::is_integer) {
            doubleSum = util::roundToInt(doubleSum);
        }
        sum[i] = doubleSum;
    }
    return sum;
}

template <typename T, int nDim>
Array<T, nDim> computeAverageAcc(MultiTensorField3D<T, nDim> &tensorField, Box3D domain)
{
    AccTensorSumFunctional3D<T, nDim> functional;
    applyProcessingFunctional(functional, domain, tensorField);
    auto tensor = functional.getSumTensor();
    for (int i = 0; i < nDim; ++i) {
        tensor[i] /= (T)domain.nCells();
    }
    return tensor;
}

template <typename T, int nDim>
Array<T, nDim> computeAverageAcc(MultiTensorField3D<T, nDim> &tensorField)
{
    return computeAverageAcc(tensorField, tensorField.getBoundingBox());
}

template <typename T, int nDim>
void AccBulkVorticityOrderEightFunctional3D<T, nDim>::process(
    Box3D domain, TensorField3D<T, nDim> &velocity, TensorField3D<T, nDim> &vorticity)
{
    Dot3D offset = computeRelativeDisplacement(velocity, vorticity);
    plint ny2 = vorticity.getNy();
    plint nz2 = vorticity.getNz();
    velocity.for_each(
        domain, [&velocity, &vorticity, offset, ny2, nz2](
                    plint i, plint iX, plint iY, plint iZ, Array<T, nDim> const &vel) {
            plint iX2 = iX + offset.x;
            plint iY2 = iY + offset.y;
            plint iZ2 = iZ + offset.z;
            vorticity.get(iX2, iY2, iZ2)[0] =
                fdDataField::bulkVorticityXOrderEight(velocity, iX, iY, iZ);
            vorticity.get(iX2, iY2, iZ2)[1] =
                fdDataField::bulkVorticityYOrderEight(velocity, iX, iY, iZ);
            vorticity.get(iX2, iY2, iZ2)[2] =
                fdDataField::bulkVorticityZOrderEight(velocity, iX, iY, iZ);
        });
}

template <typename T, int nDim>
AccBulkVorticityOrderEightFunctional3D<T, nDim>
    *AccBulkVorticityOrderEightFunctional3D<T, nDim>::clone() const
{
    return new AccBulkVorticityOrderEightFunctional3D<T, nDim>(*this);
}

template <typename T, int nDim>
void AccBulkVorticityOrderEightFunctional3D<T, nDim>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T, int nDim>
BlockDomain::DomainT AccBulkVorticityOrderEightFunctional3D<T, nDim>::appliesTo() const
{
    // Don't apply to envelope, because nearest neighbors need to be accessed.
    return BlockDomain::bulk;
}

template <typename T>
void computeBulkVorticityOrderEightAcc(
    MultiTensorField3D<T, 3> &velocity, MultiTensorField3D<T, 3> &vorticity)
{
    applyProcessingFunctional(
        new AccBulkVorticityOrderEightFunctional3D<T, 3>, velocity.getBoundingBox(), velocity,
        vorticity);
}

template <typename T>
std::unique_ptr<MultiTensorField3D<T, 3>> computeBulkVorticityOrderEightAcc(
    MultiTensorField3D<T, 3> &velocity)
{
    int envelopeWidth = 1;
    std::unique_ptr<MultiTensorField3D<T, 3>> vorticity =
        generateMultiTensorField<T, 3>((MultiBlock3D &)velocity, envelopeWidth);
    computeBulkVorticityOrderEightAcc(velocity, *vorticity);
    return vorticity;
}

/* *************** ShanChenMultiComponentAccelerated3D ***************** */

template <typename T, template <typename U> class Descriptor, int numSpecies>
ShanChenMultiComponentAccelerated3D<T, Descriptor, numSpecies>::ShanChenMultiComponentAccelerated3D(
    std::vector<std::vector<T>> const &speciesG_) :
    G((T)0)
{
    // Although speciesG_ has a 2D "matrix structure", speciesG has a 1D "array structure".
    speciesG.resize(numSpecies * numSpecies);
    for (pluint iSpecies = 0; iSpecies < numSpecies; iSpecies++) {
        PLB_ASSERT(speciesG_[iSpecies].size() == numSpecies);
        for (pluint jSpecies = 0; jSpecies < numSpecies; jSpecies++) {
            speciesG[iSpecies * numSpecies + jSpecies] = speciesG_[iSpecies][jSpecies];
        }
    }
}

template <typename T, template <typename U> class Descriptor, int numSpecies>
ShanChenMultiComponentAccelerated3D<T, Descriptor, numSpecies>::ShanChenMultiComponentAccelerated3D(
    T G_, std::vector<T> const &imposedOmega_) :
    G(G_), imposedOmega(imposedOmega_)
{
    speciesG.resize(numSpecies * numSpecies, G);
}

template <typename T, template <typename U> class Descriptor, int numSpecies>
ShanChenMultiComponentAccelerated3D<T, Descriptor, numSpecies>::ShanChenMultiComponentAccelerated3D(
    std::vector<std::vector<T>> const &speciesG_, std::vector<T> const &imposedOmega_) :
    G((T)0), imposedOmega(imposedOmega_)
{
    // Although speciesG_ has a 2D "matrix structure", speciesG has a 1D "array structure".
    speciesG.resize(numSpecies * numSpecies);
    for (pluint iSpecies = 0; iSpecies < numSpecies; iSpecies++) {
        PLB_ASSERT(speciesG_[iSpecies].size() == numSpecies);
        for (pluint jSpecies = 0; jSpecies < numSpecies; jSpecies++) {
            speciesG[iSpecies * numSpecies + jSpecies] = speciesG_[iSpecies][jSpecies];
        }
    }
}

template <typename T, template <typename U> class Descriptor, int numSpecies>
void ShanChenMultiComponentAccelerated3D<T, Descriptor, numSpecies>::process(
    Box3D domain, std::vector<AtomicAcceleratedLattice3D<T, Descriptor> *> lattices)
{
    // Short-hand notation for the lattice descriptor
    typedef Descriptor<T> D;
    // Handle to external scalars
    enum {
        densityOffset = D::ExternalField::densityBeginsAt,
        momentumOffset = D::ExternalField::momentumBeginsAt,
        forceOffset = D::ExternalField::forceBeginsAt
    };

    auto latticesPtr = &lattices[0];
    auto omegaPtr = &imposedOmega[0];
    auto speciesGptr = &speciesG[0];

    // Compute per-lattice density  and momentum on every site and on each
    //   lattice, and store result in external scalars;  envelope cells are included,
    //   because they are needed to compute the interaction potential in the following.
    //   Note that the per-lattice value of the momentum is stored temporarily only, as
    //   it is corrected later on, based on the common fluid velocity.
    for (plint iSpecies = 0; iSpecies < numSpecies; ++iSpecies) {
        lattices[0]->for_each([latticesPtr, iSpecies](plint i, int collisionModel) {
            T rho;
            Array<T, Descriptor<T>::d> j;
            Array<T, Descriptor<T>::q> f;
            latticesPtr[iSpecies]->pullPop(i, f);
            if (collisionModel == CollisionModel::BounceBack) {
                T *dynamicScalars;
                plint scalarIndex;
                latticesPtr[iSpecies]->getDynamicScalar(i, dynamicScalars, scalarIndex);
                PLB_ASSERT(scalarIndex != -1);
                rho = dynamicScalars[scalarIndex];
                momentTemplatesImpl<T, typename Descriptor<T>::BaseDescriptor>::get_j(f, j);
            } else {
                T rhoBar;
                momentTemplatesImpl<T, typename Descriptor<T>::BaseDescriptor>::get_rhoBar_j(
                    f, rhoBar, j);
                rho = Descriptor<T>::fullRho(rhoBar);
            }

            latticesPtr[iSpecies]->pushExt(i, densityOffset, rho);
            latticesPtr[iSpecies]->pushExt(i, momentumOffset + 0, j[0]);
            latticesPtr[iSpecies]->pushExt(i, momentumOffset + 1, j[1]);
            latticesPtr[iSpecies]->pushExt(i, momentumOffset + 2, j[2]);
        });
    }

    // Compute the interaction force between the species, and store it by
    //   means of a velocity correction in the external velocity field.
    lattices[0]->for_each(
        domain, [latticesPtr, omegaPtr, speciesGptr](
                    plint i, plint iX, plint iY, plint iZ, int collisionModel) {
            T weightedDensity = T();
            for (plint iSpecies = 0; iSpecies < numSpecies; ++iSpecies) {
                weightedDensity +=
                    omegaPtr[iSpecies] * latticesPtr[iSpecies]->pullExt(i, densityOffset);
            }
            // Computation of the common velocity, shared among all populations.
            Array<T, Descriptor<T>::d> uTot;
            for (int iD = 0; iD < Descriptor<T>::d; ++iD) {
                uTot[iD] = T();
                for (plint iSpecies = 0; iSpecies < numSpecies; ++iSpecies) {
                    Array<T, Descriptor<T>::ExternalField::numScalars> ext;
                    uTot[iD] +=
                        omegaPtr[iSpecies] * latticesPtr[iSpecies]->pullExt(i, momentumOffset + iD);
                }
                uTot[iD] /= weightedDensity;
            }

            // Computation of the interaction potential.
            Array<T, D::d * numSpecies> rhoContributions;
            for (plint iSpecies = 0; iSpecies < numSpecies; ++iSpecies) {
                Array<T, D::d> rhoContribution;
                accelMultiPhaseTemplates3D<T, Descriptor>::shanChenInteraction(
                    *latticesPtr[iSpecies], rhoContribution, iX, iY, iZ);
                rhoContributions[D::d * iSpecies + 0] = rhoContribution[0];
                rhoContributions[D::d * iSpecies + 1] = rhoContribution[1];
                rhoContributions[D::d * iSpecies + 2] = rhoContribution[2];
            }

            // Computation and storage of the final velocity, consisting
            //   of uTot plus the momentum difference due to interaction
            //   potential and external force
            for (plint iSpecies = 0; iSpecies < numSpecies; ++iSpecies) {
                for (int iD = 0; iD < D::d; ++iD) {
                    T momentumContribution = uTot[iD];
                    // Initialize force contribution with force from external fields if there
                    //   is any, or with zero otherwise.
                    T forceContribution = latticesPtr[iSpecies]->pullExt(i, forceOffset + iD);
                    // Then, add a contribution from the potential of all other species.
                    for (plint iPartnerSpecies = 0; iPartnerSpecies < numSpecies; ++iPartnerSpecies)
                    {
                        if (iPartnerSpecies != iSpecies) {
                            forceContribution -=
                                speciesGptr[iSpecies * numSpecies + iPartnerSpecies]
                                * rhoContributions[iPartnerSpecies * D::d + iD];
                            // rhoContribution[iPartnerSpecies][iD];
                        }
                    }
                    momentumContribution += 1. / omegaPtr[iSpecies] * forceContribution;
                    // Multiply by rho to convert from velocity to momentum.
                    momentumContribution *= latticesPtr[iSpecies]->pullExt(i, densityOffset);
                    latticesPtr[iSpecies]->pushExt(i, momentumOffset + iD, momentumContribution);
                }
            }
        });
}

template <typename T, template <typename U> class Descriptor, int numSpecies>
ShanChenMultiComponentAccelerated3D<T, Descriptor, numSpecies>
    *ShanChenMultiComponentAccelerated3D<T, Descriptor, numSpecies>::clone() const
{
    return new ShanChenMultiComponentAccelerated3D<T, Descriptor, numSpecies>(*this);
}

template <typename T, template <typename U> class Descriptor, int numSpecies>
void ShanChenMultiComponentAccelerated3D<T, Descriptor, numSpecies>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    // All blocks are modified by the Shan/Chen processor.
    for (pluint iBlock = 0; iBlock < modified.size(); ++iBlock) {
        modified[iBlock] = modif::staticVariables;
    }
}

// https://stackoverflow.com/questions/3418231/replace-part-of-a-string-with-another-string
// HemoFlow patch: 'inline' added, the header is included by several translation units.
inline void replaceAll(std::string &str, const std::string &from, const std::string &to)
{
    if (from.empty())
        return;
    size_t start_pos = 0;
    while ((start_pos = str.find(from, start_pos)) != std::string::npos) {
        str.replace(start_pos, from.length(), to);
        start_pos += to.length();  // In case 'to' contains 'from', like replacing 'x' with 'yx'
    }
}

template <typename T, template <typename U> class Descriptor>
void showTemplateArguments(MultiBlockLattice3D<T, Descriptor> &lattice)
{
    std::map<int, std::string> nameOfDynamics;
    auto dynField = extractDynamicsChain(lattice, nameOfDynamics);
    pcout << "collideAndStream(CollisionKernel<T, DESCRIPTOR";
    for (auto iter = nameOfDynamics.begin(); iter != nameOfDynamics.end(); ++iter) {
        pcout << ",\n                                 ";
        pcout << "CollisionModel::";
        std::string name = iter->second;
        replaceAll(name, " >> ", "__");
        replaceAll(name, "-", "M");
        pcout << name;
    }
    pcout << ">() );" << std::endl;
}

}  // namespace plb

#endif  // ACCELERATED_PROCESSORS_3D_HH
