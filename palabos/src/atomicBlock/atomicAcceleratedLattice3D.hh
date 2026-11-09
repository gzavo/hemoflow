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
 * The dynamics of a 3D block lattice -- generic implementation.
 */
#ifndef ATOMIC_ACCELERATED_LATTICE_3D_HH
#define ATOMIC_ACCELERATED_LATTICE_3D_HH

#include "atomicBlock/atomicAcceleratedLattice3D.h"
#include "acceleratedLattice/acceleratedDefinitions.h"
#include "core/dynamics.h"
#include "core/cell.h"
#include "core/plbTimer.h"
#include "latticeBoltzmann/latticeTemplates.h"
#include "latticeBoltzmann/indexTemplates.h"
#include "core/util.h"
#include "core/latticeStatistics.h"
#include "core/dynamicsIdentifiers.h"
#include "core/plbProfiler.h"
#include "acceleratedLattice/acceleratedCollisions.h"
#include <algorithm>
#include <numeric>
#include <execution>
#include <typeinfo>
#include <cmath>
#include <vector>
#include <functional>
#include <utility>
#include <ranges>

namespace plb {

// Class AtomicAcceleratedLattice3D /////////////////////////

template <typename T, template <typename U> class Descriptor>
template <class Functional>
void AtomicAcceleratedLattice3D<T, Descriptor>::for_each(Functional functional)
{
    std::for_each(
        std::execution::par_unseq, this->collisionMatrix, this->collisionMatrix + N,
        [this, functional](int const &collisionModel) {
            plint i = (plint)(&collisionModel - this->collisionMatrix);
            functional(i, collisionModel);
        });
}

template <typename T, template <typename U> class Descriptor>
template <class Functional>
void AtomicAcceleratedLattice3D<T, Descriptor>::for_each(Box3D const &domain, Functional functional)
{
    std::for_each(
        std::execution::par_unseq, this->collisionMatrix, this->collisionMatrix + N,
        [this, domain, functional](int const &collisionModel) {
            plint i = (plint)(&collisionModel - this->collisionMatrix);
            plint nynz = this->getNy() * this->getNz();
            plint iX = i / nynz;
            plint remainder = i % nynz;
            plint iY = remainder / this->getNz();
            plint iZ = remainder % this->getNz();
            if (iX >= domain.x0 && iX <= domain.x1 && iY >= domain.y0 && iY <= domain.y1
                && iZ >= domain.z0 && iZ <= domain.z1)
            {
                functional(i, iX, iY, iZ, this->collisionMatrix[i]);
            }
        });
}

template <typename T, template <typename U> class Descriptor>
template <class Functional, class BinaryReductionOp>
T AtomicAcceleratedLattice3D<T, Descriptor>::transform_reduce(
    Box3D const &domain, T neutral, BinaryReductionOp reduce, Functional functional)
{
    auto indices = std::views::iota(plint(0), N);
    return std::transform_reduce(
        std::execution::par_unseq, indices.begin(), indices.end(), neutral, reduce,
        [this, domain, functional, neutral](plint i) {
            plint nynz = this->getNy() * this->getNz();
            plint iX = i / nynz;
            plint remainder = i % nynz;
            plint iY = remainder / this->getNz();
            plint iZ = remainder % this->getNz();
            if (iX >= domain.x0 && iX <= domain.x1 && iY >= domain.y0 && iY <= domain.y1
                && iZ >= domain.z0 && iZ <= domain.z1)
            {
                return functional(i, iX, iY, iZ, this->collisionMatrix[i]);
            } else {
                return neutral;
            }
        });
}

/** \param nx_ lattice width (first index)
 *  \param ny_ lattice height (second index)
 *  \param nz_ lattice depth (third index)
 */
template <typename T, template <typename U> class Descriptor>
AtomicAcceleratedLattice3D<T, Descriptor>::AtomicAcceleratedLattice3D(
    plint nx_, plint ny_, plint nz_, Dynamics<T, Descriptor> *backgroundDynamics_) :
    AtomicBlock3D(nx_, ny_, nz_, new AcceleratedLatticeDataTransfer3D<T, Descriptor>()),
    backgroundDynamics(backgroundDynamics_),
    N(nx_ * ny_ * nz_),
    timeCounter(0)
{
    // Allocate memory, and initialize dynamics.
    allocateAndInitialize();
    std::for_each(
        std::execution::par_unseq, dynamicsArray, dynamicsArray + N, [this](auto &dynamics) {
            dynamics = backgroundDynamics;
        });
    int backgroundCollisionModel = toCollisionModel(*backgroundDynamics);
    std::for_each(
        std::execution::par_unseq, collisionMatrix, collisionMatrix + N,
        [this, backgroundCollisionModel](int &collisionModel) {
            collisionModel = backgroundCollisionModel;
            size_t i = &collisionModel - collisionMatrix;
            dynamicScalarIndex[i] = -1;
        });
    std::for_each(std::execution::par_unseq, hw_bb_links, hw_bb_links + N, [this](uint32_t &link) {
        link = uint32_t {};
    });
    int backgroundCollisionNumScalars = numDynamicScalars<T, Descriptor>(backgroundCollisionModel);
    if (backgroundCollisionNumScalars > 0) {
        std::vector<T> backgroundCollisionScalars =
            getDynamicScalars<T, Descriptor>(*backgroundDynamics, backgroundCollisionModel);
        for (plint i = 0; i < N; ++i) {
            dynamicScalarIndex[i] = (plint)dynamicScalars.size();
            dynamicScalars.insert(
                dynamicScalars.end(), backgroundCollisionScalars.begin(),
                backgroundCollisionScalars.end());
        }
    }

    // Attribute default value to the standard statistics (average uSqr,
    //   max uSqr, average rho). These have previously been subscribed
    //   in the constructor of BlockLatticeBase3D.
    std::vector<double> average, sum, max;
    std::vector<plint> intSum;
    average.push_back(Descriptor<double>::rhoBar(1.));
    // default average rho to 1, to avoid division by
    // zero in constRhoBGK and related models
    average.push_back(0.);  // default average uSqr to 0
    max.push_back(0.);      // default max uSqr to 0
    plint numCells = 1;     // pretend fictitious cell to evaluate statistics
    this->getInternalStatistics().evaluate(average, sum, max, intSum, numCells);
    global::plbCounter("MEMORY_LATTICE").increment(allocatedMemory());
}

/** During destruction, the memory for the lattice and the contained
 * cells is released. However, the dynamics objects pointed to by
 * the cells must be deleted manually by the user.
 */
template <typename T, template <typename U> class Descriptor>
AtomicAcceleratedLattice3D<T, Descriptor>::~AtomicAcceleratedLattice3D()
{
    global::plbCounter("MEMORY_LATTICE").increment(-allocatedMemory());
    releaseMemory();
}

/** The whole data of the lattice is duplicated. This includes
 * both particle distribution function and external fields.
 * \warning The dynamics objects and internalProcessors are not copied
 * \param rhs the lattice to be duplicated
 */
template <typename T, template <typename U> class Descriptor>
AtomicAcceleratedLattice3D<T, Descriptor>::AtomicAcceleratedLattice3D(
    AtomicAcceleratedLattice3D<T, Descriptor> const &rhs) :
    AtomicBlock3D(rhs),
    backgroundDynamics(rhs.backgroundDynamics->clone()),
    N(rhs.N),
    timeCounter(rhs.timeCounter)
{
    allocateAndInitialize();
    plint numPop = Descriptor<T>::numPop;
    std::for_each(
        std::execution::par_unseq, populations, populations + N * numPop, [this, &rhs](T &f) {
            size_t i = &f - populations;
            f = rhs.populations[i];
        });
    std::for_each(
        std::execution::par_unseq, collisionMatrix, collisionMatrix + N,
        [this, &rhs](int &collisionModel) {
            size_t i = &collisionModel - collisionMatrix;
            collisionModel = rhs.collisionMatrix[i];
        });
    std::for_each(
        std::execution::par_unseq, hw_bb_links, hw_bb_links + N, [this, &rhs](uint32_t &links) {
            size_t i = &links - hw_bb_links;
            links = rhs.hw_bb_links[i];
        });
    std::for_each(
        std::execution::par_unseq, dynamicScalarIndex, dynamicScalarIndex + N,
        [this, &rhs](plint &index) {
            size_t i = &index - dynamicScalarIndex;
            index = rhs.dynamicScalarIndex[i];
        });
    dynamicScalars = rhs.dynamicScalars;
    std::for_each(
        std::execution::seq, dynamicsArray, dynamicsArray + N, [this, &rhs](auto &dynamics) {
            size_t i = &dynamics - dynamicsArray;
            if (rhs.dynamicsArray[i] == rhs.backgroundDynamics) {
                dynamics = backgroundDynamics;
            } else {
                dynamics = rhs.dynamicsArray[i]->clone();
            }
        });
    global::plbCounter("MEMORY_LATTICE").increment(allocatedMemory());
}

template <typename T, template <typename U> class Descriptor>
void AtomicAcceleratedLattice3D<T, Descriptor>::writeBack(BlockLattice3D<T, Descriptor> &rhs)
{
    PLB_ASSERT(this->getNx() == rhs.getNx());
    PLB_ASSERT(this->getNy() == rhs.getNy());
    PLB_ASSERT(this->getNz() == rhs.getNz());

    plint numPop = Descriptor<T>::numPop;
    for (plint iX = 0; iX < this->getNx(); ++iX) {
        for (plint iY = 0; iY < this->getNy(); ++iY) {
            for (plint iZ = 0; iZ < this->getNz(); ++iZ) {
                Cell<T, Descriptor> &cell = rhs.get(iX, iY, iZ);
                for (plint iPop = 0; iPop < numPop; ++iPop) {
                    cell[iPop] = populationGrid[iPop][iX][iY][iZ];
                }
                for (int iExt = 0; iExt < Descriptor<T>::ExternalField::numScalars; ++iExt) {
                    *cell.getExternal(iExt) = *externalScalarGrid[iX][iY][iZ].get(iExt);
                }
                Dynamics<T, Descriptor> *previousDynamics = &cell.getDynamics();
                if (previousDynamics != rhs.backgroundDynamics) {
                    delete previousDynamics;
                }
                cell.attributeDynamics(dynamicsGrid[iX][iY][iZ]->clone());
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
AtomicAcceleratedLattice3D<T, Descriptor>::AtomicAcceleratedLattice3D(
    BlockLattice3D<T, Descriptor> const &rhs) :
    AtomicBlock3D(
        rhs.getNx(), rhs.getNy(), rhs.getNz(),
        new AcceleratedLatticeDataTransfer3D<T, Descriptor>()),
    backgroundDynamics(rhs.backgroundDynamics->clone()),
    N(rhs.getBoundingBox().nCells()),
    timeCounter(rhs.getTimeCounter().getTime())
{
    this->setLocation(rhs.getLocation());
    this->setFlag(rhs.getFlag());
    int backgroundCollisionModel = toCollisionModel(*backgroundDynamics);
    int backgroundCollisionNumScalars = numDynamicScalars<T, Descriptor>(backgroundCollisionModel);
    std::vector<T> backgroundCollisionScalars =
        getDynamicScalars<T, Descriptor>(*backgroundDynamics, backgroundCollisionModel);
    allocateAndInitialize();
    plint numPop = Descriptor<T>::numPop;
    plint i = 0;
    for (plint iX = 0; iX < this->getNx(); ++iX) {
        for (plint iY = 0; iY < this->getNy(); ++iY) {
            for (plint iZ = 0; iZ < this->getNz(); ++iZ) {
                Cell<T, Descriptor> const &cell = rhs.get(iX, iY, iZ);
                for (plint iPop = 0; iPop < numPop; ++iPop) {
                    populationGrid[iPop][iX][iY][iZ] = cell[iPop];
                }
                for (int iExt = 0; iExt < Descriptor<T>::ExternalField::numScalars; ++iExt) {
                    *externalScalarGrid[iX][iY][iZ].get(iExt) = *cell.getExternal(iExt);
                }
                if (&cell.getDynamics() == rhs.backgroundDynamics) {
                    dynamicsGrid[iX][iY][iZ] = backgroundDynamics;
                    collisionMatrix[i] = backgroundCollisionModel;
                    if (backgroundCollisionNumScalars > 0) {
                        dynamicScalarIndex[i] = (plint)dynamicScalars.size();
                        dynamicScalars.insert(
                            dynamicScalars.end(), backgroundCollisionScalars.begin(),
                            backgroundCollisionScalars.end());
                    } else {
                        dynamicScalarIndex[i] = -1;
                    }
                } else {
                    dynamicsGrid[iX][iY][iZ] = cell.getDynamics().clone();
                    int collisionModel = toCollisionModel(*dynamicsGrid[iX][iY][iZ]);
                    collisionMatrix[i] = collisionModel;
                    if (numDynamicScalars<T, Descriptor>(collisionModel) > 0) {
                        std::vector<T> newScalars = getDynamicScalars<T, Descriptor>(
                            *dynamicsGrid[iX][iY][iZ], collisionModel);
                        dynamicScalarIndex[i] = (plint)dynamicScalars.size();
                        dynamicScalars.insert(
                            dynamicScalars.end(), newScalars.begin(), newScalars.end());
                    } else {
                        dynamicScalarIndex[i] = -1;
                    }
                }
                ++i;
            }
        }
    }

    bool hasDynamicScalars = !dynamicScalars.empty();
    T *dynamicScalarsPtr = &dynamicScalars[0];

    std::for_each(
        std::execution::par_unseq, hw_bb_links, hw_bb_links + N,
        [this, dynamicScalarsPtr, hasDynamicScalars](uint32_t &links) {
            size_t i = &links - hw_bb_links;
            int collisionModel = collisionMatrix[i];
            links = uint32_t {};

            plint index = -1;
            if (hasDynamicScalars) {
                index = dynamicScalarIndex[i];
            }
            if (collisionModel >= CollisionModel::HalfwayBounceBack__TRT) {
                for (int iPop = 0; iPop < Descriptor<T>::numPop; ++iPop) {
                    if ((index >= 0) && (!std::isnan(dynamicScalarsPtr[index + iPop]))) {
                        setBit(links, iPop);
                    }
                }
            }
        });
}

/** The current lattice is deallocated, then the lattice from the rhs
 * is duplicated. This includes both particle distribution function
 * and external fields.
 * \warning The dynamics objects and internalProcessors are not copied
 * \param rhs the lattice to be duplicated
 */
template <typename T, template <typename U> class Descriptor>
AtomicAcceleratedLattice3D<T, Descriptor> &AtomicAcceleratedLattice3D<T, Descriptor>::operator=(
    AtomicAcceleratedLattice3D<T, Descriptor> const &rhs)
{
    AtomicAcceleratedLattice3D<T, Descriptor> tmp(rhs);
    swap(tmp);
    return *this;
}

/** The swap is efficient, in the sense that only pointers to the
 * lattice are copied, and not the lattice itself.
 */
template <typename T, template <typename U> class Descriptor>
void AtomicAcceleratedLattice3D<T, Descriptor>::swap(AtomicAcceleratedLattice3D &rhs)
{
    global::plbCounter("MEMORY_LATTICE").increment(-allocatedMemory());
    AtomicBlock3D::swap(rhs);
    std::swap(backgroundDynamics, rhs.backgroundDynamics);
    std::swap(populations, rhs.populations);
    std::swap(populationGrid, rhs.populationGrid);
    std::swap(externalScalars, rhs.externalScalars);
    std::swap(externalScalarGrid, rhs.externalScalarGrid);
    std::swap(dynamicsArray, rhs.dynamicsArray);
    std::swap(dynamicsGrid, rhs.dynamicsGrid);
    std::swap(collisionMatrix, rhs.collisionMatrix);
    std::swap(hw_bb_links, rhs.hw_bb_links);
    std::swap(dynamicScalarIndex, rhs.dynamicScalarIndex);
    dynamicScalars.swap(rhs.dynamicScalars);
    std::swap(N, rhs.N);
    std::swap(timeCounter, rhs.timeCounter);
    global::plbCounter("MEMORY_LATTICE").increment(allocatedMemory());
}

template <typename T, template <typename U> class Descriptor>
plint AtomicAcceleratedLattice3D<T, Descriptor>::getN() const
{
    return N;
}

template <typename T, template <typename U> class Descriptor>
int const *AtomicAcceleratedLattice3D<T, Descriptor>::getCollisionMatrix() const
{
    return collisionMatrix;
}

template <typename T, template <typename U> class Descriptor>
void AtomicAcceleratedLattice3D<T, Descriptor>::getDynamicScalar(
    plint i, T *&dynamicScalarsPtr, plint &index)
{
    PLB_ASSERT(i < N);
    index = dynamicScalarIndex[i];
    dynamicScalarsPtr = &dynamicScalars[0];
}

template <typename T, template <typename U> class Descriptor>
void AtomicAcceleratedLattice3D<T, Descriptor>::collideAndStream()
{
    nvtx a {"collideAndStream on domain"};
    plint nx = this->getNx();
    plint ny = this->getNy();
    plint nz = this->getNz();
    plint delta = Descriptor<T>::vicinity;
    Box3D bulk(delta, nx - 1 - delta, delta, ny - 1 - delta, delta, nz - 1 - delta);
    Box3D fullDomain(this->getBoundingBox());

    bool hasDynamicScalars = !dynamicScalars.empty();
    T *dynamicScalarsPtr = hasDynamicScalars ? &dynamicScalars[0] : nullptr;
    Array<T, GPUconst<T, Descriptor>::maxStaticScalars> staticScalars;
    int backgroundCollisionModel = toCollisionModel(this->getBackgroundDynamics());
    getStaticScalars(this->getBackgroundDynamics(), backgroundCollisionModel, staticScalars);
    std::for_each(std::execution::par_unseq, populations, populations + N, [=, this](T &f0) {
        size_t i = &f0 - populations;
        plint iX = i / (ny * nz);
        plint remainder = i % (ny * nz);
        plint iY = remainder / nz;
        plint iZ = remainder % nz;

        Array<T, Descriptor<T>::numPop> f;
        Array<T, Descriptor<T>::ExternalField::numScalars> ext;
        int collisionModel = collisionMatrix[i];
        // PULL
        for (int iPop = 0; iPop < Descriptor<T>::numPop; ++iPop) {
            f[iPop] = populations[iPop * N + i];
        }
        for (int iExt = 0; iExt < Descriptor<T>::ExternalField::numScalars; ++iExt) {
            ext[iExt] = *externalScalars[i].get(iExt);
        }

        plint index = -1;
        if (hasDynamicScalars) {
            index = dynamicScalarIndex[i];
        }

        collide<T, Descriptor>(collisionModel, f, ext, staticScalars, dynamicScalarsPtr, index);

        // PUSH
        uint32_t links = hw_bb_links[i];
        for (int iPop = 0; iPop < Descriptor<T>::numPop; ++iPop) {
            if (readBit(links, iPop)) {
                int iOpp = indexTemplates::opposite<Descriptor<T>>(iPop);
                tmpPopulations[iOpp * N + i] = f[iPop];
            } else {
                plint nextX = iX + Descriptor<T>::c_gpu(iPop, 0);
                plint nextY = iY + Descriptor<T>::c_gpu(iPop, 1);
                plint nextZ = iZ + Descriptor<T>::c_gpu(iPop, 2);
                if (contained(nextX, nextY, nextZ, fullDomain)) {
                    plint iNext = nextZ + nz * (nextY + ny * nextX);
                    tmpPopulations[iPop * N + iNext] = f[iPop];
                }
            }
        }

        for (int iExt = 0; iExt < Descriptor<T>::ExternalField::numScalars; ++iExt) {
            *externalScalars[i].get(iExt) = ext[iExt];
        }
    });

    std::swap(populations, tmpPopulations);
    std::swap(populationGrid, tmpPopulationGrid);
}

template <typename T, template <typename U> class Descriptor>
template <class CollFun>
void AtomicAcceleratedLattice3D<T, Descriptor>::collideAndStream(CollFun const &collFun)
{
    nvtx a {"collideAndStream on domain"};
    plint nx = this->getNx();
    plint ny = this->getNy();
    plint nz = this->getNz();
    plint delta = Descriptor<T>::vicinity;
    Box3D bulk(delta, nx - 1 - delta, delta, ny - 1 - delta, delta, nz - 1 - delta);
    Box3D fullDomain(this->getBoundingBox());
    bool hasDynamicScalars = !dynamicScalars.empty();
    T *dynamicScalarsPtr = &dynamicScalars[0];
    Array<T, GPUconst<T, Descriptor>::maxStaticScalars> staticScalars;
    int backgroundCollisionModel = toCollisionModel(this->getBackgroundDynamics());
    getStaticScalars(this->getBackgroundDynamics(), backgroundCollisionModel, staticScalars);

    std::for_each(
        std::execution::par_unseq, populations, populations + N,
        [this, nx, ny, nz, fullDomain, collFun, hasDynamicScalars, dynamicScalarsPtr, staticScalars,
         backgroundCollisionModel](T &f0) {
            size_t i = &f0 - populations;
            plint iX = i / (ny * nz);
            plint remainder = i % (ny * nz);
            plint iY = remainder / nz;
            plint iZ = remainder % nz;

            Array<T, Descriptor<T>::numPop> f;
            Array<T, Descriptor<T>::ExternalField::numScalars> ext;
            int collisionModel = collisionMatrix[i];
            // PULL
            for (int iPop = 0; iPop < Descriptor<T>::numPop; ++iPop) {
                f[iPop] = populations[iPop * N + i];
            }
            for (int iExt = 0; iExt < Descriptor<T>::ExternalField::numScalars; ++iExt) {
                ext[iExt] = *externalScalars[i].get(iExt);
            }

            plint index = -1;
            if (hasDynamicScalars) {
                index = dynamicScalarIndex[i];
            }

            collFun(collisionModel, f, ext, staticScalars, dynamicScalarsPtr, index);

            // PUSH
            // MODIF
            uint32_t links = hw_bb_links[i];
            for (int iPop = 0; iPop < Descriptor<T>::numPop; ++iPop) {
                // if ((collisionModel >= CollisionModel::HalfwayBounceBack__TRT) &&
                //(index >= 0) && (!std::isnan(dynamicScalarsPtr[index + iPop]))) {
                if (readBit(links, iPop)) {
                    int iOpp = indexTemplates::opposite<Descriptor<T>>(iPop);
                    tmpPopulations[iOpp * N + i] = f[iPop];
                } else {
                    plint nextX = iX + Descriptor<T>::c_gpu(iPop, 0);
                    plint nextY = iY + Descriptor<T>::c_gpu(iPop, 1);
                    plint nextZ = iZ + Descriptor<T>::c_gpu(iPop, 2);
                    if (contained(nextX, nextY, nextZ, fullDomain)) {
                        plint iNext = nextZ + nz * (nextY + ny * nextX);
                        tmpPopulations[iPop * N + iNext] = f[iPop];
                    }
                }
            }
            for (int iExt = 0; iExt < Descriptor<T>::ExternalField::numScalars; ++iExt) {
                *externalScalars[i].get(iExt) = ext[iExt];
            }
        });
    std::swap(populations, tmpPopulations);
    std::swap(populationGrid, tmpPopulationGrid);
}

template <typename T, template <typename U> class Descriptor>
void AtomicAcceleratedLattice3D<T, Descriptor>::incrementTime()
{
    timeCounter++;
}

template <typename T, template <typename U> class Descriptor>
void AtomicAcceleratedLattice3D<T, Descriptor>::resetTime(plint timeCounter_)
{
    timeCounter = timeCounter_;
}

template <typename T, template <typename U> class Descriptor>
void AtomicAcceleratedLattice3D<T, Descriptor>::reconstructCell(
    plint iX, plint iY, plint iZ, Cell<T, Descriptor> &cell) const
{
    cell.attributeDynamics(dynamicsGrid[iX][iY][iZ]);
    reconstructCellStatic(iX, iY, iZ, cell);
}

template <typename T, template <typename U> class Descriptor>
void AtomicAcceleratedLattice3D<T, Descriptor>::reconstructCellStatic(
    plint iX, plint iY, plint iZ, Cell<T, Descriptor> &cell) const
{
    for (int iPop = 0; iPop < Descriptor<T>::numPop; ++iPop) {
        cell[iPop] = populationGrid[iPop][iX][iY][iZ];
    }
    for (int iExt = 0; iExt < Descriptor<T>::ExternalField::numScalars; ++iExt) {
        *cell.getExternal(iExt) = *externalScalarGrid[iX][iY][iZ].get(iExt);
    }
}

template <typename T, template <typename U> class Descriptor>
void AtomicAcceleratedLattice3D<T, Descriptor>::pullPop(
    plint i, Array<T, Descriptor<T>::q> &f) const
{
    for (int iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
        f[iPop] = populations[iPop * N + i];
    }
}

template <typename T, template <typename U> class Descriptor>
void AtomicAcceleratedLattice3D<T, Descriptor>::pullExt(
    plint iX, plint iY, plint iZ, Array<T, Descriptor<T>::ExternalField::numScalars> &ext) const
{
    for (int iExt = 0; iExt < Descriptor<T>::ExternalField::numScalars; ++iExt) {
        ext[iExt] = *externalScalarGrid[iX][iY][iZ].get(iExt);
    }
}

template <typename T, template <typename U> class Descriptor>
T AtomicAcceleratedLattice3D<T, Descriptor>::pullExt(plint iX, plint iY, plint iZ, int offset) const
{
    return *externalScalarGrid[iX][iY][iZ].get(offset);
}

template <typename T, template <typename U> class Descriptor>
T AtomicAcceleratedLattice3D<T, Descriptor>::pullExt(plint i, int offset) const
{
    return *externalScalars[i].get(offset);
}

template <typename T, template <typename U> class Descriptor>
void AtomicAcceleratedLattice3D<T, Descriptor>::pushStatic(
    plint iX, plint iY, plint iZ, Cell<T, Descriptor> const &cell)
{
    for (int iPop = 0; iPop < Descriptor<T>::numPop; ++iPop) {
        populationGrid[iPop][iX][iY][iZ] = cell[iPop];
    }
    for (int iExt = 0; iExt < Descriptor<T>::ExternalField::numScalars; ++iExt) {
        *externalScalarGrid[iX][iY][iZ].get(iExt) = *cell.getExternal(iExt);
    }
}

template <typename T, template <typename U> class Descriptor>
void AtomicAcceleratedLattice3D<T, Descriptor>::pushPop(
    plint i, Array<T, Descriptor<T>::q> const &f)
{
    for (int iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
        populations[iPop * N + i] = f[iPop];
    }
}

template <typename T, template <typename U> class Descriptor>
void AtomicAcceleratedLattice3D<T, Descriptor>::pushPop(
    plint iX, plint iY, plint iZ, Array<T, Descriptor<T>::q> const &f)
{
    for (int iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
        populationGrid[iPop][iX][iY][iZ] = f[iPop];
    }
}

template <typename T, template <typename U> class Descriptor>
void AtomicAcceleratedLattice3D<T, Descriptor>::pushExt(
    plint iX, plint iY, plint iZ, int offset, T value)
{
    *externalScalarGrid[iX][iY][iZ].get(offset) = value;
}

template <typename T, template <typename U> class Descriptor>
void AtomicAcceleratedLattice3D<T, Descriptor>::pushExt(plint i, int offset, T value)
{
    *externalScalars[i].get(offset) = value;
}

template <typename T, template <typename U> class Descriptor>
void AtomicAcceleratedLattice3D<T, Descriptor>::allocateAndInitialize()
{
    this->getInternalStatistics().subscribeAverage();  // Subscribe average rho-bar
    this->getInternalStatistics().subscribeAverage();  // Subscribe average uSqr
    this->getInternalStatistics().subscribeMax();      // Subscribe max uSqr

    plint numPop = Descriptor<T>::numPop;

    plint nx = this->getNx();
    plint ny = this->getNy();
    plint nz = this->getNz();
    populations = new T[numPop * N];
    tmpPopulations = new T[numPop * N];
    externalScalars = new External[N];
    dynamicsArray = new Dynamics<T, Descriptor> *[N];
    collisionMatrix = new int[N];
    hw_bb_links = new uint32_t[N];
    dynamicScalarIndex = new plint[N];

    populationGrid = new T ***[numPop];
    tmpPopulationGrid = new T ***[numPop];
    for (plint iPop = 0; iPop < numPop; ++iPop) {
        populationGrid[iPop] = new T **[nx];
        tmpPopulationGrid[iPop] = new T **[nx];
        for (plint iX = 0; iX < nx; ++iX) {
            populationGrid[iPop][iX] = new T *[ny];
            tmpPopulationGrid[iPop][iX] = new T *[ny];
            for (plint iY = 0; iY < ny; ++iY) {
                populationGrid[iPop][iX][iY] = populations + nz * (iY + ny * (iX + nx * iPop));
                tmpPopulationGrid[iPop][iX][iY] =
                    tmpPopulations + nz * (iY + ny * (iX + nx * iPop));
            }
        }
    }

    externalScalarGrid = new External **[nx];
    dynamicsGrid = new Dynamics<T, Descriptor> ***[nx];
    for (plint iX = 0; iX < nx; ++iX) {
        externalScalarGrid[iX] = new External *[ny];
        dynamicsGrid[iX] = new Dynamics<T, Descriptor> **[ny];
        for (plint iY = 0; iY < ny; ++iY) {
            externalScalarGrid[iX][iY] = externalScalars + nz * (iY + ny * iX);
            dynamicsGrid[iX][iY] = dynamicsArray + nz * (iY + ny * iX);
        }
    }
}

template <typename T, template <typename U> class Descriptor>
void AtomicAcceleratedLattice3D<T, Descriptor>::releaseMemory()
{
    plint N = this->getBoundingBox().nCells();

    std::for_each(dynamicsArray, dynamicsArray + N, [this, N](auto &dynamics) {
        if (dynamics != backgroundDynamics)
            delete dynamics;
    });

    delete backgroundDynamics;
    delete[] populations;
    delete[] tmpPopulations;
    delete[] externalScalars;
    delete[] dynamicsArray;
    delete[] collisionMatrix;
    delete[] hw_bb_links;
    delete[] dynamicScalarIndex;

    plint numPop = Descriptor<T>::numPop;
    for (int iPop = 0; iPop < numPop; ++iPop) {
        for (plint iX = 0; iX < getNx(); ++iX) {
            delete[] populationGrid[iPop][iX];
            delete[] tmpPopulationGrid[iPop][iX];
        }
        delete[] populationGrid[iPop];
        delete[] tmpPopulationGrid[iPop];
    }
    delete[] populationGrid;
    delete[] tmpPopulationGrid;

    for (plint iX = 0; iX < getNx(); ++iX) {
        delete[] externalScalarGrid[iX];
        delete[] dynamicsGrid[iX];
    }
    delete[] externalScalarGrid;
    delete[] dynamicsGrid;
}

template <typename T, template <typename U> class Descriptor>
void AtomicAcceleratedLattice3D<T, Descriptor>::attributeDynamics(
    plint iX, plint iY, plint iZ, Dynamics<T, Descriptor> *dynamics)
{
    Dynamics<T, Descriptor> *previousDynamics = dynamicsGrid[iX][iY][iZ];
    if (previousDynamics != backgroundDynamics) {
        delete previousDynamics;
    }
    dynamicsGrid[iX][iY][iZ] = dynamics;
    plint nz = this->getNz();
    plint ny = this->getNy();
    int collisionModel = toCollisionModel(*dynamics);
    int collisionNumScalars = numDynamicScalars<T, Descriptor>(collisionModel);
    plint i = iZ + nz * (iY + ny * iX);
    if (collisionNumScalars > 0) {
        std::vector<T> collisionScalars =
            getDynamicScalars<T, Descriptor>(*dynamics, collisionModel);
        dynamicScalarIndex[i] = (plint)dynamicScalars.size();
        dynamicScalars.insert(
            dynamicScalars.end(), collisionScalars.begin(), collisionScalars.end());
    } else {
        dynamicScalarIndex[i] = -1;
    }
}

template <typename T, template <typename U> class Descriptor>
Dynamics<T, Descriptor> &AtomicAcceleratedLattice3D<T, Descriptor>::getBackgroundDynamics()
{
    return *backgroundDynamics;
}

template <typename T, template <typename U> class Descriptor>
Dynamics<T, Descriptor> const &AtomicAcceleratedLattice3D<T, Descriptor>::getBackgroundDynamics()
    const
{
    return *backgroundDynamics;
}

template <typename T, template <typename U> class Descriptor>
void AtomicAcceleratedLattice3D<T, Descriptor>::resetDynamics(
    Dynamics<T, Descriptor> const &dynamics)
{
    plint N = this->getBoundingBox().nCells();
    std::for_each(
        std::execution::seq, dynamicsArray, dynamicsArray + N, [this, &dynamics](auto &dyn) {
            if (dyn != backgroundDynamics) {
                delete dyn;
            }
            dyn = dynamics.clone();
        });

    std::vector<T>().swap(dynamicScalars);
    int backgroundCollisionModel = toCollisionModel(*backgroundDynamics);
    int backgroundCollisionNumScalars = numDynamicScalars<T, Descriptor>(backgroundCollisionModel);
    if (backgroundCollisionNumScalars > 0) {
        std::vector<T> backgroundCollisionScalars =
            getDynamicScalars<T, Descriptor>(*backgroundDynamics, backgroundCollisionModel);
        for (plint i = 0; i < N; ++i) {
            dynamicScalarIndex[i] = (plint)dynamicScalars.size();
            dynamicScalars.insert(
                dynamicScalars.end(), backgroundCollisionScalars.begin(),
                backgroundCollisionScalars.end());
        }
    }
}

template <typename T, template <typename U> class Descriptor>
Dynamics<T, Descriptor> const &AtomicAcceleratedLattice3D<T, Descriptor>::getDynamics(
    plint iX, plint iY, plint iZ) const
{
    return *dynamicsGrid[iX][iY][iZ];
}

template <typename T, template <typename U> class Descriptor>
Dynamics<T, Descriptor> &AtomicAcceleratedLattice3D<T, Descriptor>::getDynamics(
    plint iX, plint iY, plint iZ)
{
    return *dynamicsGrid[iX][iY][iZ];
}

template <typename T, template <typename U> class Descriptor>
plint AtomicAcceleratedLattice3D<T, Descriptor>::allocatedMemory() const
{
    return this->getBoundingBox().nCells() * sizeof(T)
           * (Descriptor<T>::numPop + Descriptor<T>::ExternalField::numScalars);
}

////////////////////// Class AcceleratedLatticeDataTransfer3D /////////////////////////

template <typename T, template <typename U> class Descriptor>
AcceleratedLatticeDataTransfer3D<T, Descriptor>::AcceleratedLatticeDataTransfer3D() :
    lattice(0), constLattice(0)
{ }

template <typename T, template <typename U> class Descriptor>
void AcceleratedLatticeDataTransfer3D<T, Descriptor>::setBlock(AtomicBlock3D &block)
{
    lattice = dynamic_cast<AtomicAcceleratedLattice3D<T, Descriptor> *>(&block);
    PLB_ASSERT(lattice);
    constLattice = lattice;
}

template <typename T, template <typename U> class Descriptor>
void AcceleratedLatticeDataTransfer3D<T, Descriptor>::setConstBlock(AtomicBlock3D const &block)
{
    constLattice = dynamic_cast<AtomicAcceleratedLattice3D<T, Descriptor> const *>(&block);
    PLB_ASSERT(constLattice);
}

template <typename T, template <typename U> class Descriptor>
AcceleratedLatticeDataTransfer3D<T, Descriptor>
    *AcceleratedLatticeDataTransfer3D<T, Descriptor>::clone() const
{
    return new AcceleratedLatticeDataTransfer3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
plint AcceleratedLatticeDataTransfer3D<T, Descriptor>::staticCellSize() const
{
    return sizeof(T) * (Descriptor<T>::numPop + Descriptor<T>::ExternalField::numScalars);
}

template <typename T, template <typename U> class Descriptor>
void AcceleratedLatticeDataTransfer3D<T, Descriptor>::send(
    Box3D domain, std::vector<char> &buffer, modif::ModifT kind) const
{
    PLB_PRECONDITION(constLattice);
    PLB_PRECONDITION(contained(domain, constLattice->getBoundingBox()));
    // It's the responsibility of the functions called below to allocate
    //   the right amount of memory for the buffer.
    buffer.clear();
    switch (kind) {
    case modif::staticVariables:
        send_static(domain, buffer);
        break;
    case modif::dynamicVariables:
        send_dynamic(domain, buffer);
        break;
    // Serialization is the same no matter if the dynamics object
    //   is being regenerated or not by the recipient.
    case modif::allVariables:
    case modif::dataStructure:
        send_all(domain, buffer);
        break;
    default:
        PLB_ASSERT(false);
    }
}

template <typename T, template <typename U> class Descriptor>
void AcceleratedLatticeDataTransfer3D<T, Descriptor>::send_raw(
    Box3D destDomain, Box3D destBoundingBox, plint deltaX, plint deltaY, plint deltaZ,
    char **buffer, plint &bufferSize, char **indices, modif::ModifT kind) const
{
    PLB_PRECONDITION(constLattice);
    // It's the responsibility of the functions called below to allocate
    //   the right amount of memory for the buffer.
    switch (kind) {
    case modif::staticVariables:
        send_static_raw(
            destDomain, destBoundingBox, deltaX, deltaY, deltaZ, buffer, bufferSize, indices);
        break;
    case modif::dynamicVariables:
        PLB_ASSERT(false);
        break;
    // Serialization is the same no matter if the dynamics object
    //   is being regenerated or not by the recipient.
    case modif::allVariables:
    case modif::dataStructure:
        PLB_ASSERT(false);
        break;
    default:
        PLB_ASSERT(false);
    }
}

template <typename T, template <typename U> class Descriptor>
void AcceleratedLatticeDataTransfer3D<T, Descriptor>::send_static(
    Box3D domain, std::vector<char> &buffer) const
{
    nvtx a {"packing on domain"};
    typedef ExternalFieldArray<T, typename Descriptor<T>::ExternalField> External;
    PLB_PRECONDITION(constLattice);
    plint cellSize = staticCellSize();
    pluint numBytes = domain.nCells() * cellSize;
    plint numExt = Descriptor<T>::ExternalField::numScalars;
    // Avoid dereferencing uninitialized pointer.
    if (numBytes == 0)
        return;

    if (buffer.size() < numBytes) {
        buffer.resize(numBytes);
    }

    plint Ndomain = domain.nCells();
    plint Nlattice = constLattice->getBoundingBox().nCells();
    T const *populations = constLattice->populations;
    External const *externalScalars = constLattice->externalScalars;
    plint domain_x0 = domain.x0;
    plint domain_y0 = domain.y0;
    plint domain_z0 = domain.z0;
    plint domain_ny = domain.getNy();
    plint domain_nz = domain.getNz();
    plint ny = constLattice->getNy();
    plint nz = constLattice->getNz();
    T *bufferPtr = (T *)&buffer[0];
    plint numScalarsInCell = Descriptor<T>::q + numExt;
    std::for_each(std::execution::par_unseq, bufferPtr, bufferPtr + Ndomain, [=](auto &value) {
        plint i = (plint)(&value - bufferPtr);
        plint iX = domain_x0 + i / (domain_ny * domain_nz);
        plint remainder = i % (domain_ny * domain_nz);
        plint iY = domain_y0 + remainder / domain_nz;
        plint iZ = domain_z0 + remainder % domain_nz;
        plint iAbsolute = iZ + nz * (iY + ny * iX);
        for (plint iPop = 0; iPop < Descriptor<T>::numPop; ++iPop) {
            bufferPtr[iPop + numScalarsInCell * i] = populations[iPop * Nlattice + iAbsolute];
        }
        for (plint iExt = 0; iExt < numExt; ++iExt) {
            bufferPtr[numScalarsInCell * i + Descriptor<T>::q + iExt] =
                *externalScalars[iAbsolute].get(iExt);
        }
    });
}

template <typename T, template <typename U> class Descriptor>
void AcceleratedLatticeDataTransfer3D<T, Descriptor>::send_static_raw(
    Box3D destDomain, Box3D destBoundingBox, plint deltaX, plint deltaY, plint deltaZ,
    char **buffer, plint &bufferSize, char **indices) const
{
    Box3D domain(destDomain.shift(deltaX, deltaY, deltaZ));
    nvtx a {"packing on domain"};
    typedef ExternalFieldArray<T, typename Descriptor<T>::ExternalField> External;
    PLB_PRECONDITION(constLattice);

    plint numExt = Descriptor<T>::ExternalField::numScalars;
    plint Ndomain = domain.nCells();

    // Avoid dereferencing uninitialized pointer.
    if (Ndomain == 0) {
        bufferSize = 0;
        *buffer = nullptr;
        return;
    }

    plint Nlattice = constLattice->getBoundingBox().nCells();
    T const *populations = constLattice->populations;
    External const *externalScalars = constLattice->externalScalars;

    plint to_x0 = destDomain.x0;
    plint to_y0 = destDomain.y0;
    plint to_z0 = destDomain.z0;

    plint domain_ny = domain.getNy();
    plint domain_nz = domain.getNz();
    plint from_ny = constLattice->getNy();
    plint from_nz = constLattice->getNz();

    plint *indexPtr = nullptr;
    if (*buffer == nullptr) {
        PLB_ASSERT(*indices == nullptr);
        plint *sizeOfCellPtr;
        allocateBytes((char **)&sizeOfCellPtr, Ndomain * sizeof(plint));

        // First, generate the sequence of i [0, Ndomain)
#ifdef __NVCOMPILER
        auto i = std::views::iota(plint(0), Ndomain);
#else
        std::vector<plint> i(Ndomain);
        std::iota(i.begin(), i.end(), plint(0));
#endif

        // Parallel transform step: fill sizeOfCellPtr based on the original lambda logic
        std::transform(
            std::execution::par_unseq, i.begin(), i.end(), sizeOfCellPtr,
            [=](plint iCell) -> plint {
                plint iX_to = to_x0 + iCell / (domain_ny * domain_nz);
                plint remainder = iCell % (domain_ny * domain_nz);
                plint iY_to = to_y0 + remainder / domain_nz;
                plint iZ_to = to_z0 + remainder % domain_nz;

                plint sizeOfCell = numExt;
                for (plint iPop = 0; iPop < Descriptor<T>::numPop; ++iPop) {
                    plint iX_prev = iX_to - Descriptor<T>::c_gpu(iPop, 0);
                    plint iY_prev = iY_to - Descriptor<T>::c_gpu(iPop, 1);
                    plint iZ_prev = iZ_to - Descriptor<T>::c_gpu(iPop, 2);
                    bool selectPopulation = !contained(iX_prev, iY_prev, iZ_prev, destBoundingBox);
                    if (selectPopulation) {
                        ++sizeOfCell;
                    }
                }
                return sizeOfCell;
            });

        // Parallel reduce step: sum all the elements of sizeOfCellPtr
        plint numFloats =
            std::reduce(std::execution::par, sizeOfCellPtr, sizeOfCellPtr + Ndomain, plint(0));

        allocateBytes((char **)indices, Ndomain * sizeof(plint));
        indexPtr = (plint *)*indices;
        std::exclusive_scan(
            std::execution::par_unseq, sizeOfCellPtr, sizeOfCellPtr + Ndomain, indexPtr, (plint)0);
        releaseBytes((char *)sizeOfCellPtr);

        bufferSize = numFloats * sizeof(T);
        allocateBytes(buffer, bufferSize);
    } else {
        PLB_ASSERT(*indices != nullptr);
        indexPtr = (plint *)*indices;
    }

    T *bufferPtr = (T *)*buffer;
    std::for_each(std::execution::par_unseq, indexPtr, indexPtr + Ndomain, [=](auto &bufferIndex) {
        plint i = (plint)(&bufferIndex - indexPtr);
        plint iX_to = to_x0 + i / (domain_ny * domain_nz);
        plint remainder = i % (domain_ny * domain_nz);
        plint iY_to = to_y0 + remainder / domain_nz;
        plint iZ_to = to_z0 + remainder % domain_nz;

        // Origin = Destination + delta
        plint iX_from = iX_to + deltaX;
        plint iY_from = iY_to + deltaY;
        plint iZ_from = iZ_to + deltaZ;
        plint iFrom = iZ_from + from_nz * (iY_from + from_ny * iX_from);

        plint nOffset = 0;
        for (plint iPop = 0; iPop < Descriptor<T>::numPop; ++iPop) {
            plint iX_prev = iX_to - Descriptor<T>::c_gpu(iPop, 0);
            plint iY_prev = iY_to - Descriptor<T>::c_gpu(iPop, 1);
            plint iZ_prev = iZ_to - Descriptor<T>::c_gpu(iPop, 2);
            bool selectPopulation = !contained(iX_prev, iY_prev, iZ_prev, destBoundingBox);
            if (selectPopulation) {
                bufferPtr[nOffset + bufferIndex] = populations[iPop * Nlattice + iFrom];
                ++nOffset;
            }
        }
        for (plint iExt = 0; iExt < numExt; ++iExt) {
            bufferPtr[nOffset + iExt + bufferIndex] = *externalScalars[iFrom].get(iExt);
        }
    });
}

template <typename T, template <typename U> class Descriptor>
void AcceleratedLatticeDataTransfer3D<T, Descriptor>::send_dynamic(
    Box3D domain, std::vector<char> &buffer) const
{
    PLB_PRECONDITION(constLattice);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                // The serialize function automatically reallocates memory for buffer.
                serialize(*constLattice->dynamicsGrid[iX][iY][iZ], buffer);
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
void AcceleratedLatticeDataTransfer3D<T, Descriptor>::send_all(
    Box3D domain, std::vector<char> &buffer) const
{
    PLB_PRECONDITION(constLattice);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                // 1. Send dynamic info (automaic allocation of buffer memory).
                serialize(*constLattice->dynamicsGrid[iX][iY][iZ], buffer);
                pluint pos = buffer.size();
                // 2. Send static info (needs manual allocation of buffer memory).
                if (staticCellSize() > 0) {
                    buffer.resize(pos + staticCellSize());
                    Cell<T, Descriptor> cell;
                    constLattice->reconstructCell(iX, iY, iZ, cell);
                    cell.serialize(&buffer[pos]);
                }
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
void AcceleratedLatticeDataTransfer3D<T, Descriptor>::receive(
    Box3D domain, std::vector<char> const &buffer, modif::ModifT kind,
    std::map<int, std::string> const &foreignIds)
{
    receive(domain, buffer, kind);
}

template <typename T, template <typename U> class Descriptor>
void AcceleratedLatticeDataTransfer3D<T, Descriptor>::receive(
    Box3D domain, std::vector<char> const &buffer, modif::ModifT kind)
{
    PLB_PRECONDITION(lattice);
    PLB_PRECONDITION(contained(domain, lattice->getBoundingBox()));
    switch (kind) {
    case modif::staticVariables:
        receive_static(domain, buffer);
        break;
    case modif::dynamicVariables:
        receive_dynamic(domain, buffer);
        break;
    case modif::allVariables:
        receive_all(domain, buffer);
        break;
    case modif::dataStructure:
        receive_regenerate(domain, buffer);
        break;
    default:
        PLB_ASSERT(false);
    }
}

template <typename T, template <typename U> class Descriptor>
void AcceleratedLatticeDataTransfer3D<T, Descriptor>::receive_raw(
    Box3D domain, char const *buffer, plint bufferSize, char **indices, modif::ModifT kind)
{
    PLB_PRECONDITION(lattice);
    PLB_PRECONDITION(contained(domain, lattice->getBoundingBox()));
    switch (kind) {
    case modif::staticVariables:
        receive_static_raw(domain, buffer, bufferSize, indices);
        break;
    case modif::dynamicVariables:
        PLB_ASSERT(false);
        break;
    case modif::allVariables:
        PLB_ASSERT(false);
        break;
    case modif::dataStructure:
        PLB_ASSERT(false);
        break;
    default:
        PLB_ASSERT(false);
    }
}

template <typename T, template <typename U> class Descriptor>
void AcceleratedLatticeDataTransfer3D<T, Descriptor>::receive_static_raw(
    Box3D domain, char const *buffer, plint bufferSize, char **indices)
{
    nvtx a {"un-packing on domain"};
    typedef ExternalFieldArray<T, typename Descriptor<T>::ExternalField> External;
    PLB_PRECONDITION(lattice);
    // PLB_PRECONDITION( (plint) buffer.size() == domain.nCells()*staticCellSize() );
    //  Avoid dereferencing uninitialized pointer.
    if (bufferSize == 0)
        return;
    plint numExt = Descriptor<T>::ExternalField::numScalars;
    plint Ndomain = domain.nCells();

    Box3D fullDomain = constLattice->getBoundingBox();
    plint Nlattice = fullDomain.nCells();
    plint domain_x0 = domain.x0;
    plint domain_y0 = domain.y0;
    plint domain_z0 = domain.z0;
    plint domain_ny = domain.getNy();
    plint domain_nz = domain.getNz();
    plint full_ny = fullDomain.getNy();
    plint full_nz = fullDomain.getNz();
    T *populations = lattice->populations;
    External *externalScalars = lattice->externalScalars;
    T const *bufferPtr = (T const *)buffer;

    plint *indexPtr = nullptr;
    if (*indices == nullptr) {
        plint *sizeOfCellPtr;
        allocateBytes((char **)&sizeOfCellPtr, Ndomain * sizeof(plint));

        // First, generate the sequence of i [0, Ndomain)
#ifdef __NVCOMPILER
        auto i = std::views::iota(plint(0), Ndomain);
#else
        std::vector<plint> i(Ndomain);
        std::iota(i.begin(), i.end(), plint(0));
#endif

        // Step 1: Parallel transform to compute sizeOfCell for each cell index
        std::transform(
            std::execution::par_unseq, i.begin(), i.end(), sizeOfCellPtr, [=](plint iCell) {
                plint iX = domain_x0 + iCell / (domain_ny * domain_nz);
                plint remainder = iCell % (domain_ny * domain_nz);
                plint iY = domain_y0 + remainder / domain_nz;
                plint iZ = domain_z0 + remainder % domain_nz;

                plint sizeOfCell = numExt;
                for (plint iPop = 0; iPop < Descriptor<T>::numPop; ++iPop) {
                    plint iX_prev = iX - Descriptor<T>::c_gpu(iPop, 0);
                    plint iY_prev = iY - Descriptor<T>::c_gpu(iPop, 1);
                    plint iZ_prev = iZ - Descriptor<T>::c_gpu(iPop, 2);
                    bool selectPopulation = !contained(iX_prev, iY_prev, iZ_prev, fullDomain);
                    if (selectPopulation) {
                        ++sizeOfCell;
                    }
                }
                return sizeOfCell;
            });

        allocateBytes((char **)indices, Ndomain * sizeof(plint));
        indexPtr = (plint *)*indices;
        std::exclusive_scan(
            std::execution::par_unseq, sizeOfCellPtr, sizeOfCellPtr + Ndomain, indexPtr, (plint)0);
        releaseBytes((char *)sizeOfCellPtr);
    } else {
        indexPtr = (plint *)*indices;
    }

    std::for_each(std::execution::par_unseq, indexPtr, indexPtr + Ndomain, [=](auto &bufferIndex) {
        plint i = (plint)(&bufferIndex - indexPtr);
        plint iX = domain_x0 + i / (domain_ny * domain_nz);
        plint remainder = i % (domain_ny * domain_nz);
        plint iY = domain_y0 + remainder / domain_nz;
        plint iZ = domain_z0 + remainder % domain_nz;
        plint iAbsolute = iZ + full_nz * (iY + full_ny * iX);
        plint nOffset = 0;
        for (plint iPop = 0; iPop < Descriptor<T>::numPop; ++iPop) {
            plint iX_prev = iX - Descriptor<T>::c_gpu(iPop, 0);
            plint iY_prev = iY - Descriptor<T>::c_gpu(iPop, 1);
            plint iZ_prev = iZ - Descriptor<T>::c_gpu(iPop, 2);
            bool selectPopulation = !contained(iX_prev, iY_prev, iZ_prev, fullDomain);
            if (selectPopulation) {
                populations[iPop * Nlattice + iAbsolute] = bufferPtr[nOffset + bufferIndex];
                ++nOffset;
            }
        }
        for (plint iExt = 0; iExt < numExt; ++iExt) {
            *externalScalars[iAbsolute].get(iExt) = bufferPtr[nOffset + iExt + bufferIndex];
        }
    });
}

template <typename T, template <typename U> class Descriptor>
void AcceleratedLatticeDataTransfer3D<T, Descriptor>::receive_static(
    Box3D domain, std::vector<char> const &buffer)
{
    nvtx a {"un-packing on domain"};
    typedef ExternalFieldArray<T, typename Descriptor<T>::ExternalField> External;
    PLB_PRECONDITION(lattice);
    // PLB_PRECONDITION( (plint) buffer.size() == domain.nCells()*staticCellSize() );
    //  Avoid dereferencing uninitialized pointer.
    if (buffer.empty())
        return;
    plint numExt = Descriptor<T>::ExternalField::numScalars;
    plint Ndomain = domain.nCells();

    plint Nlattice = constLattice->getBoundingBox().nCells();
    plint domain_x0 = domain.x0;
    plint domain_y0 = domain.y0;
    plint domain_z0 = domain.z0;
    plint domain_ny = domain.getNy();
    plint domain_nz = domain.getNz();
    plint ny = constLattice->getNy();
    plint nz = constLattice->getNz();
    T *populations = lattice->populations;
    External *externalScalars = lattice->externalScalars;
    T const *bufferPtr = (T const *)&buffer[0];
    plint numScalarsInCell = Descriptor<T>::q + numExt;
    std::for_each(std::execution::par_unseq, bufferPtr, bufferPtr + Ndomain, [=](auto &value) {
        plint i = (plint)(&value - bufferPtr);
        plint iX = domain_x0 + i / (domain_ny * domain_nz);
        plint remainder = i % (domain_ny * domain_nz);
        plint iY = domain_y0 + remainder / domain_nz;
        plint iZ = domain_z0 + remainder % domain_nz;
        plint iAbsolute = iZ + nz * (iY + ny * iX);
        for (plint iPop = 0; iPop < Descriptor<T>::numPop; ++iPop) {
            populations[iPop * Nlattice + iAbsolute] = bufferPtr[iPop + numScalarsInCell * i];
            // populations[iPop * Nlattice + iAbsolute] = tmpPop[iPop + Descriptor<T>::numPop * i];
        }
        for (plint iExt = 0; iExt < numExt; ++iExt) {
            *externalScalars[iAbsolute].get(iExt) =
                bufferPtr[numScalarsInCell * i + Descriptor<T>::q + iExt];
            //*externalScalars[iAbsolute].get(iExt) = *external.get(iExt);
        }
    });
}

template <typename T, template <typename U> class Descriptor>
void AcceleratedLatticeDataTransfer3D<T, Descriptor>::receive_dynamic(
    Box3D domain, std::vector<char> const &buffer)
{
    PLB_PRECONDITION(lattice);
    pluint serializerPos = 0;
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                // No assert is included here, because incompatible types of
                //   dynamics are detected by asserts inside HierarchicUnserializer.
                serializerPos =
                    unserialize(*lattice->dynamicsGrid[iX][iY][iZ], buffer, serializerPos);
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
void AcceleratedLatticeDataTransfer3D<T, Descriptor>::receive_all(
    Box3D domain, std::vector<char> const &buffer)
{
    PLB_PRECONDITION(lattice);
    pluint posInBuffer = 0;
    plint cellSize = staticCellSize();
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                // 1. Unserialize dynamic data.
                posInBuffer = unserialize(*lattice->dynamicsGrid[iX][iY][iZ], buffer, posInBuffer);
                // 2. Unserialize static data.
                if (staticCellSize() > 0) {
                    Cell<T, Descriptor> cell;
                    cell.unSerialize(&buffer[posInBuffer]);
                    lattice->pushStatic(iX, iY, iZ, cell);
                    posInBuffer += cellSize;
                }
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
void AcceleratedLatticeDataTransfer3D<T, Descriptor>::receive_regenerate(
    Box3D domain, std::vector<char> const &buffer, std::map<int, int> const &idIndirect)
{
    PLB_PRECONDITION(lattice);
    pluint posInBuffer = 0;
    plint cellSize = staticCellSize();
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                // 1. Generate dynamics object, and unserialize dynamic data.
                std::map<int, int> const *indirectPtr = idIndirect.empty() ? 0 : &idIndirect;
                HierarchicUnserializer unserializer(buffer, posInBuffer, indirectPtr);
                Dynamics<T, Descriptor> *newDynamics =
                    meta::dynamicsRegistration<T, Descriptor>().generate(unserializer);
                posInBuffer = unserializer.getCurrentPos();
                lattice->attributeDynamics(iX, iY, iZ, newDynamics);

                // 2. Unserialize static data.
                if (staticCellSize() > 0) {
                    PLB_ASSERT(!buffer.empty());
                    PLB_ASSERT(posInBuffer + cellSize <= buffer.size());
                    Cell<T, Descriptor> cell;
                    cell.unSerialize(&buffer[posInBuffer]);
                    lattice->pushStatic(iX, iY, iZ, cell);
                    posInBuffer += cellSize;
                }
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
void AcceleratedLatticeDataTransfer3D<T, Descriptor>::attribute(
    Box3D toDomain, plint deltaX, plint deltaY, plint deltaZ, AtomicBlock3D const &from,
    modif::ModifT kind)
{
    PLB_PRECONDITION(lattice);
    PLB_PRECONDITION(typeid(from) == typeid(AtomicAcceleratedLattice3D<T, Descriptor> const &));
    PLB_PRECONDITION(contained(toDomain, lattice->getBoundingBox()));
    AtomicAcceleratedLattice3D<T, Descriptor> const &fromLattice =
        (AtomicAcceleratedLattice3D<T, Descriptor> const &)from;
    switch (kind) {
    case modif::staticVariables:
        attribute_static(toDomain, deltaX, deltaY, deltaZ, fromLattice);
        break;
    case modif::dynamicVariables:
        attribute_dynamic(toDomain, deltaX, deltaY, deltaZ, fromLattice);
        break;
    case modif::allVariables:
        attribute_all(toDomain, deltaX, deltaY, deltaZ, fromLattice);
        break;
    case modif::dataStructure:
        attribute_regenerate(toDomain, deltaX, deltaY, deltaZ, fromLattice);
        break;
    default:
        PLB_ASSERT(false);
    }
}

template <typename T, template <typename U> class Descriptor>
void AcceleratedLatticeDataTransfer3D<T, Descriptor>::attribute_static(
    Box3D toDomain, plint deltaX, plint deltaY, plint deltaZ,
    AtomicAcceleratedLattice3D<T, Descriptor> const &from)
{
    PLB_PRECONDITION(lattice);
    plint numExt = Descriptor<T>::ExternalField::numScalars;

    plint to_x0 = toDomain.x0;
    plint to_y0 = toDomain.y0;
    plint to_z0 = toDomain.z0;
    plint domain_ny = toDomain.getNy();
    plint domain_nz = toDomain.getNz();
    plint lattice_ny = constLattice->getNy();
    plint lattice_nz = constLattice->getNz();
    plint Ndomain = toDomain.nCells();
    Box3D destBoundingBox = constLattice->getBoundingBox();
    plint Nlattice = destBoundingBox.nCells();
    T *to_populations = lattice->populations;
    T const *from_populations = from.populations;
    ExternalFieldArray<T, typename Descriptor<T>::ExternalField> *to_ext = lattice->externalScalars;
    ExternalFieldArray<T, typename Descriptor<T>::ExternalField> const *from_ext =
        from.externalScalars;
    std::for_each(
        std::execution::par_unseq, to_populations, to_populations + Ndomain, [=](auto &population) {
            plint i = (plint)(&population - to_populations);
            plint iX_to = to_x0 + i / (domain_ny * domain_nz);
            plint remainder_to = i % (domain_ny * domain_nz);
            plint iY_to = to_y0 + remainder_to / domain_nz;
            plint iZ_to = to_z0 + remainder_to % domain_nz;
            plint iTo = iZ_to + lattice_nz * (iY_to + lattice_ny * iX_to);

            plint iX_from = iX_to + deltaX;
            plint iY_from = iY_to + deltaY;
            plint iZ_from = iZ_to + deltaZ;
            plint iFrom = iZ_from + lattice_nz * (iY_from + lattice_ny * iX_from);

            for (plint iPop = 0; iPop < Descriptor<T>::numPop; ++iPop) {
                plint iX_prev = iX_to - Descriptor<T>::c_gpu(iPop, 0);
                plint iY_prev = iY_to - Descriptor<T>::c_gpu(iPop, 1);
                plint iZ_prev = iZ_to - Descriptor<T>::c_gpu(iPop, 2);
                bool selectPopulation = !contained(iX_prev, iY_prev, iZ_prev, destBoundingBox);
                if (selectPopulation) {
                    to_populations[iPop * Nlattice + iTo] =
                        from_populations[iPop * Nlattice + iFrom];
                }
            }
            for (plint iExt = 0; iExt < numExt; ++iExt) {
                *to_ext[iTo].get(iExt) = *from_ext[iFrom].get(iExt);
            }
        });
}

template <typename T, template <typename U> class Descriptor>
void AcceleratedLatticeDataTransfer3D<T, Descriptor>::attribute_dynamic(
    Box3D toDomain, plint deltaX, plint deltaY, plint deltaZ,
    AtomicAcceleratedLattice3D<T, Descriptor> const &from)
{
    PLB_PRECONDITION(lattice);
    std::vector<char> serializedData;
    for (plint iX = toDomain.x0; iX <= toDomain.x1; ++iX) {
        for (plint iY = toDomain.y0; iY <= toDomain.y1; ++iY) {
            for (plint iZ = toDomain.z0; iZ <= toDomain.z1; ++iZ) {
                serializedData.clear();
                serialize(
                    *from.dynamicsGrid[iX + deltaX][iY + deltaY][iZ + deltaZ], serializedData);
                unserialize(*lattice->dynamicsGrid[iX][iY][iZ], serializedData);
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
void AcceleratedLatticeDataTransfer3D<T, Descriptor>::attribute_all(
    Box3D toDomain, plint deltaX, plint deltaY, plint deltaZ,
    AtomicAcceleratedLattice3D<T, Descriptor> const &from)
{
    std::vector<char> serializedData;
    for (plint iX = toDomain.x0; iX <= toDomain.x1; ++iX) {
        for (plint iY = toDomain.y0; iY <= toDomain.y1; ++iY) {
            for (plint iZ = toDomain.z0; iZ <= toDomain.z1; ++iZ) {
                // 1. Attribute dynamic content.
                serializedData.clear();
                serialize(
                    *from.dynamicsGrid[iX + deltaX][iY + deltaY][iZ + deltaZ], serializedData);
                unserialize(*lattice->dynamicsGrid[iX][iY][iZ], serializedData);

                // 2. Attribute static content.
                for (plint iPop = 0; iPop < Descriptor<T>::numPop; ++iPop) {
                    lattice->populationGrid[iPop][iX][iY][iZ] =
                        from.populationGrid[iPop][iX + deltaX][iY + deltaY][iZ + deltaZ];
                }
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
void AcceleratedLatticeDataTransfer3D<T, Descriptor>::attribute_regenerate(
    Box3D toDomain, plint deltaX, plint deltaY, plint deltaZ,
    AtomicAcceleratedLattice3D<T, Descriptor> const &from)
{
    PLB_PRECONDITION(lattice);
    std::vector<char> serializedData;
    for (plint iX = toDomain.x0; iX <= toDomain.x1; ++iX) {
        for (plint iY = toDomain.y0; iY <= toDomain.y1; ++iY) {
            for (plint iZ = toDomain.z0; iZ <= toDomain.z1; ++iZ) {
                // 1. Generate new dynamics and attribute dynamic content.
                serializedData.clear();
                serialize(
                    *from.dynamicsGrid[iX + deltaX][iY + deltaY][iZ + deltaZ], serializedData);
                HierarchicUnserializer unserializer(serializedData, 0);
                Dynamics<T, Descriptor> *newDynamics =
                    meta::dynamicsRegistration<T, Descriptor>().generate(unserializer);
                lattice->attributeDynamics(iX, iY, iZ, newDynamics);

                // 2. Attribute static content.
                for (plint iPop = 0; iPop < Descriptor<T>::numPop; ++iPop) {
                    lattice->populationGrid[iPop][iX][iY][iZ] =
                        from.populationGrid[iPop][iX + deltaX][iY + deltaY][iZ + deltaZ];
                }
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
CachePolicy3D &AtomicAcceleratedLattice3D<T, Descriptor>::cachePolicy()
{
    static CachePolicy3D cachePolicySingleton(30);
    return cachePolicySingleton;
}

/////////// Free Functions //////////////////////////////

template <typename T, template <typename U> class Descriptor>
double getStoredAverageDensity(AtomicAcceleratedLattice3D<T, Descriptor> const &blockLattice)
{
    return Descriptor<T>::fullRho(
        blockLattice.getInternalStatistics().getAverage(LatticeStatistics::avRhoBar));
}

template <typename T, template <typename U> class Descriptor>
double getStoredAverageEnergy(AtomicAcceleratedLattice3D<T, Descriptor> const &blockLattice)
{
    return 0.5 * blockLattice.getInternalStatistics().getAverage(LatticeStatistics::avUSqr);
}

template <typename T, template <typename U> class Descriptor>
double getStoredMaxVelocity(AtomicAcceleratedLattice3D<T, Descriptor> const &blockLattice)
{
    return std::sqrt(blockLattice.getInternalStatistics().getMax(LatticeStatistics::maxUSqr));
}

}  // namespace plb

#endif  // ATOMIC_ACCELERATED_LATTICE_3D_HH
