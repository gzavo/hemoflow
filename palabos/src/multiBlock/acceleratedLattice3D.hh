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
 * A 3D multiblock lattice -- generic implementation.
 */
#ifndef ACCELERATED_LATTICE_3D_HH
#define ACCELERATED_LATTICE_3D_HH

#include "multiBlock/multiBlockLattice3D.h"
#include "atomicBlock/blockLattice3D.h"
#include "multiBlock/defaultMultiBlockPolicy3D.h"
#include "multiBlock/nonLocalTransfer3D.h"
#include "multiBlock/multiBlockGenerator3D.h"
#include "core/latticeStatistics.h"
#include "core/plbTypenames.h"
#include "core/multiBlockIdentifiers3D.h"
#include "core/plbProfiler.h"
#include "core/dynamicsIdentifiers.h"
#include "dataProcessors/metaStuffWrapper3D.h"
#include <algorithm>
#include <limits>
#include <cmath>

namespace plb {

////////////////////// Class AcceleratedLattice3D /////////////////////////

template <typename T, template <typename U> class Descriptor>
const int AcceleratedLattice3D<T, Descriptor>::staticId = meta::registerMultiBlock3D(
    AcceleratedLattice3D<T, Descriptor>::basicType(),
    AcceleratedLattice3D<T, Descriptor>::descriptorType(),
    AcceleratedLattice3D<T, Descriptor>::blockName(),
    defaultGenerateAcceleratedLattice3D<T, Descriptor>);

template <typename T, template <typename U> class Descriptor>
void AcceleratedLattice3D<T, Descriptor>::setExecutionMode(ExecutionMode executionMode_)
{
    executionMode = executionMode_;
    for (typename BlockMap::iterator it = atomicLattices.begin(); it != atomicLattices.end(); ++it)
    {
        it->second->setExecutionMode(executionMode);
    }
}

template <typename T, template <typename U> class Descriptor>
AcceleratedLattice3D<T, Descriptor>::AcceleratedLattice3D(
    MultiBlockManagement3D const &multiBlockManagement_, BlockCommunicator3D *blockCommunicator_,
    CombinedStatistics *combinedStatistics_, Dynamics<T, Descriptor> *backgroundDynamics_) :
    MultiBlock3D(multiBlockManagement_, blockCommunicator_, combinedStatistics_),
    backgroundDynamics(backgroundDynamics_),
    timeCounter(0)
{
    allocateAndInitialize();
}

template <typename T, template <typename U> class Descriptor>
AcceleratedLattice3D<T, Descriptor>::AcceleratedLattice3D(
    plint nx, plint ny, plint nz, Dynamics<T, Descriptor> *backgroundDynamics_) :
    MultiBlock3D(nx, ny, nz, Descriptor<T>::vicinity),
    backgroundDynamics(backgroundDynamics_),
    timeCounter(0)
{
    allocateAndInitialize();
}

template <typename T, template <typename U> class Descriptor>
AcceleratedLattice3D<T, Descriptor>::~AcceleratedLattice3D()
{
    for (typename BlockMap::iterator it = atomicLattices.begin(); it != atomicLattices.end(); ++it)
    {
        delete it->second;
    }
    delete backgroundDynamics;
}

template <typename T, template <typename U> class Descriptor>
AcceleratedLattice3D<T, Descriptor>::AcceleratedLattice3D(
    AcceleratedLattice3D<T, Descriptor> const &rhs) :
    MultiBlock3D(rhs),
    backgroundDynamics(rhs.backgroundDynamics->clone()),
    timeCounter(rhs.timeCounter)
{
    for (typename BlockMap::const_iterator it = rhs.atomicLattices.begin();
         it != rhs.atomicLattices.end(); ++it)
    {
        atomicLattices[it->first] = new AtomicAcceleratedLattice3D<T, Descriptor>(*it->second);
    }
}

template <typename T, template <typename U> class Descriptor>
AcceleratedLattice3D<T, Descriptor>::AcceleratedLattice3D(
    MultiBlockLattice3D<T, Descriptor> const &rhs) :
    MultiBlock3D(rhs),
    backgroundDynamics(rhs.getBackgroundDynamics().clone()),
    timeCounter(rhs.getTimeCounter().getTime())
{
    for (auto it = rhs.getBlockLattices().begin(); it != rhs.getBlockLattices().end(); ++it) {
        atomicLattices[it->first] = new AtomicAcceleratedLattice3D<T, Descriptor>(*it->second);
    }
}

template <typename T, template <typename U> class Descriptor>
void AcceleratedLattice3D<T, Descriptor>::writeBack(MultiBlockLattice3D<T, Descriptor> &rhs)
{
    for (auto it = rhs.getBlockLattices().begin(); it != rhs.getBlockLattices().end(); ++it) {
        atomicLattices[it->first]->writeBack(*it->second);
    }
}

template <typename T, template <typename U> class Descriptor>
AcceleratedLattice3D<T, Descriptor>::AcceleratedLattice3D(MultiBlock3D const &rhs)
    // Use MultiBlock's sub-domain constructor to avoid that the data-processors are copied
    :
    MultiBlock3D(rhs, rhs.getBoundingBox(), false),
    backgroundDynamics(new NoDynamics<T, Descriptor>),
    timeCounter(0)
{
    allocateAndInitialize();
}

template <typename T, template <typename U> class Descriptor>
AcceleratedLattice3D<T, Descriptor>::AcceleratedLattice3D(
    MultiBlock3D const &rhs, Box3D subDomain, bool crop) :
    MultiBlock3D(rhs, subDomain, crop),
    backgroundDynamics(new NoDynamics<T, Descriptor>),
    timeCounter(0)
{
    allocateAndInitialize();
}

template <typename T, template <typename U> class Descriptor>
void AcceleratedLattice3D<T, Descriptor>::swap(AcceleratedLattice3D<T, Descriptor> &rhs)
{
    MultiBlock3D::swap(rhs);
    std::swap(backgroundDynamics, rhs.backgroundDynamics);
    atomicLattices.swap(rhs.atomicLattices);
    std::swap(timeCounter, rhs.timeCounter);
}

template <typename T, template <typename U> class Descriptor>
AcceleratedLattice3D<T, Descriptor> &AcceleratedLattice3D<T, Descriptor>::operator=(
    AcceleratedLattice3D<T, Descriptor> const &rhs)
{
    AcceleratedLattice3D<T, Descriptor> tmp(rhs);
    swap(tmp);
    return *this;
}

template <typename T, template <typename U> class Descriptor>
AcceleratedLattice3D<T, Descriptor> *AcceleratedLattice3D<T, Descriptor>::clone() const
{
    return new AcceleratedLattice3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
AcceleratedLattice3D<T, Descriptor> *AcceleratedLattice3D<T, Descriptor>::clone(
    MultiBlockManagement3D const &newManagement) const
{
    AcceleratedLattice3D<T, Descriptor> *newLattice = new AcceleratedLattice3D<T, Descriptor>(
        newManagement, this->getBlockCommunicator().clone(), this->getCombinedStatistics().clone(),
        getBackgroundDynamics().clone());
    // Make sure background dynamics is not used in the cloned lattice, and instead every cell
    // has an individual dynamics instance. This is the behavior we expect most of the time
    // when cloning a lattice. As a matter of fact, the background dynamics is now mostly
    // considered to be dangereous, and is kept as a default behavior of the multi-block lattice
    // for legacy reasons only.
    newLattice->resetDynamics(getBackgroundDynamics());
    // Use the same domain in the "from" and "to" argument, so that the data is not shifted
    // in space during the creation of the new block.
    copy(
        *this, newLattice->getBoundingBox(), *newLattice, newLattice->getBoundingBox(),
        modif::dataStructure);
    return newLattice;
}

template <typename T, template <typename U> class Descriptor>
void AcceleratedLattice3D<T, Descriptor>::resetDynamics(Dynamics<T, Descriptor> const &dynamics)
{
    for (typename BlockMap::const_iterator it = atomicLattices.begin(); it != atomicLattices.end();
         ++it)
    {
        atomicLattices[it->first]->resetDynamics(dynamics);
    }
}

template <typename T, template <typename U> class Descriptor>
Dynamics<T, Descriptor> const &AcceleratedLattice3D<T, Descriptor>::getBackgroundDynamics() const
{
    return *backgroundDynamics;
}

template <typename T, template <typename U> class Descriptor>
template <class CollFun>
void AcceleratedLattice3D<T, Descriptor>::collideAndStream(CollFun const &collFun)
{
    global::profiler().start("cycle");
    {
        nvtx a {"COLLIDE-AND-STREAM"};
        collideAndStreamImplementation<CollFun>(collFun);
    }

    {
        nvtx b {"COMMUNICATE"};
        this->executeInternalProcessors();
    }
    this->incrementTime();
    global::profiler().stop("cycle");
    if (global::profiler().cyclingIsAutomatic()) {
        global::profiler().cycle();
    }
}

template <typename T, template <typename U> class Descriptor>
void AcceleratedLattice3D<T, Descriptor>::collideAndStream()
{
    global::profiler().start("cycle");
    collideAndStreamImplementation();
    this->executeInternalProcessors();
    this->incrementTime();
    global::profiler().stop("cycle");
    if (global::profiler().cyclingIsAutomatic()) {
        global::profiler().cycle();
    }
}

template <typename T, template <typename U> class Descriptor>
template <class CollFun>
void AcceleratedLattice3D<T, Descriptor>::externalCollideAndStream(CollFun const &collFun)
{
    global::profiler().start("cycle");
    collideAndStreamImplementation<CollFun>(collFun);
    if (global::profiler().cyclingIsAutomatic()) {
        global::profiler().cycle();
    }
    global::profiler().stop("cycle");
}

template <typename T, template <typename U> class Descriptor>
template <class CollFun>
void AcceleratedLattice3D<T, Descriptor>::collideAndStreamImplementation(CollFun const &collFun)
{
    for (typename BlockMap::iterator it = atomicLattices.begin(); it != atomicLattices.end(); ++it)
    {
        it->second->template collideAndStream<CollFun>(collFun);
    }
}

template <typename T, template <typename U> class Descriptor>
void AcceleratedLattice3D<T, Descriptor>::collideAndStreamImplementation()
{
    for (typename BlockMap::iterator it = atomicLattices.begin(); it != atomicLattices.end(); ++it)
    {
        it->second->collideAndStream();
    }
}

template <typename T, template <typename U> class Descriptor>
void AcceleratedLattice3D<T, Descriptor>::incrementTime()
{
    for (typename BlockMap::iterator it = atomicLattices.begin(); it != atomicLattices.end(); ++it)
    {
        it->second->incrementTime();
    }
    ++timeCounter;
}

template <typename T, template <typename U> class Descriptor>
void AcceleratedLattice3D<T, Descriptor>::resetTime(pluint value)
{
    for (typename BlockMap::iterator it = atomicLattices.begin(); it != atomicLattices.end(); ++it)
    {
        it->second->resetTime(value);
    }
    timeCounter = value;
}

template <typename T, template <typename U> class Descriptor>
void AcceleratedLattice3D<T, Descriptor>::allocateAndInitialize()
{
    for (pluint iBlock = 0; iBlock < this->getLocalInfo().getBlocks().size(); ++iBlock) {
        plint blockId = this->getLocalInfo().getBlocks()[iBlock];
        SmartBulk3D bulk(this->getMultiBlockManagement(), blockId);
        Box3D envelope = bulk.computeEnvelope();
        AtomicAcceleratedLattice3D<T, Descriptor> *newLattice =
            new AtomicAcceleratedLattice3D<T, Descriptor>(
                envelope.getNx(), envelope.getNy(), envelope.getNz(), backgroundDynamics->clone());
        newLattice->setLocation(Dot3D(envelope.x0, envelope.y0, envelope.z0));
        atomicLattices[blockId] = newLattice;
    }
}

template <typename T, template <typename U> class Descriptor>
std::map<plint, AtomicAcceleratedLattice3D<T, Descriptor> *>
    &AcceleratedLattice3D<T, Descriptor>::getBlockLattices()
{
    return atomicLattices;
}

template <typename T, template <typename U> class Descriptor>
std::map<plint, AtomicAcceleratedLattice3D<T, Descriptor> *> const &
    AcceleratedLattice3D<T, Descriptor>::getBlockLattices() const
{
    return atomicLattices;
}

template <typename T, template <typename U> class Descriptor>
void AcceleratedLattice3D<T, Descriptor>::getDynamicsDict(
    Box3D domain, std::map<std::string, int> &dict)
{
    std::vector<int> ids;
    uniqueDynamicsIds(*this, domain, ids);
    dict.clear();
    for (pluint i = 0; i < ids.size(); ++i) {
        int id = ids[i];
        std::string name = meta::dynamicsRegistration<T, Descriptor>().getName(id);
        dict.insert(std::pair<std::string, int>(name, id));
    }
}

template <typename T, template <typename U> class Descriptor>
std::string AcceleratedLattice3D<T, Descriptor>::getBlockName() const
{
    return std::string("AtomicAcceleratedLattice3D");
}

template <typename T, template <typename U> class Descriptor>
std::vector<std::string> AcceleratedLattice3D<T, Descriptor>::getTypeInfo() const
{
    std::vector<std::string> info;
    info.push_back(basicType());
    info.push_back(descriptorType());
    return info;
}

template <typename T, template <typename U> class Descriptor>
std::string AcceleratedLattice3D<T, Descriptor>::blockName()
{
    return std::string("AtomicAcceleratedLattice3D");
}

template <typename T, template <typename U> class Descriptor>
std::string AcceleratedLattice3D<T, Descriptor>::basicType()
{
    return std::string(NativeType<T>::getName());
}

template <typename T, template <typename U> class Descriptor>
std::string AcceleratedLattice3D<T, Descriptor>::descriptorType()
{
    return std::string(Descriptor<T>::name);
}

template <typename T, template <typename U> class Descriptor>
AtomicAcceleratedLattice3D<T, Descriptor> &AcceleratedLattice3D<T, Descriptor>::getComponent(
    plint blockId)
{
    typename BlockMap::iterator it = atomicLattices.find(blockId);
    PLB_ASSERT(it != atomicLattices.end());
    return *it->second;
}

template <typename T, template <typename U> class Descriptor>
AtomicAcceleratedLattice3D<T, Descriptor> const &AcceleratedLattice3D<T, Descriptor>::getComponent(
    plint blockId) const
{
    typename BlockMap::const_iterator it = atomicLattices.find(blockId);
    PLB_ASSERT(it != atomicLattices.end());
    return *it->second;
}

template <typename T, template <typename U> class Descriptor>
plint AcceleratedLattice3D<T, Descriptor>::sizeOfCell() const
{
    return sizeof(T) * (Descriptor<T>::numPop + Descriptor<T>::ExternalField::numScalars);
}

template <typename T, template <typename U> class Descriptor>
plint AcceleratedLattice3D<T, Descriptor>::getCellDim() const
{
    return Descriptor<T>::numPop + Descriptor<T>::ExternalField::numScalars;
}

template <typename T, template <typename U> class Descriptor>
int AcceleratedLattice3D<T, Descriptor>::getStaticId() const
{
    return staticId;
}

template <typename T, template <typename U> class Descriptor>
void AcceleratedLattice3D<T, Descriptor>::copyReceive(
    MultiBlock3D const &fromBlock, Box3D const &fromDomain, Box3D const &toDomain,
    modif::ModifT whichData)
{
    AcceleratedLattice3D<T, Descriptor> const *fromLattice =
        dynamic_cast<AcceleratedLattice3D<T, Descriptor> const *>(&fromBlock);
    PLB_ASSERT(fromLattice);
    copy(*fromLattice, fromDomain, *this, toDomain, whichData);
}

/////////// Free Functions //////////////////////////////

template <typename T, template <typename U> class Descriptor>
AcceleratedLattice3D<T, Descriptor> &findAcceleratedLattice3D(id_t id)
{
    MultiBlock3D *multiBlock = multiBlockRegistration3D().find(id);
    if (!multiBlock || multiBlock->getStaticId() != AcceleratedLattice3D<T, Descriptor>::staticId) {
        throw PlbLogicException("Trying to access a multi block lattice that is not registered.");
    }
    return (AcceleratedLattice3D<T, Descriptor> &)(*multiBlock);
}

}  // namespace plb

#endif  // ACCELERATED_LATTICE_3D_HH
