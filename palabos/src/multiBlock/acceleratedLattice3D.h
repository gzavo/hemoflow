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
 * A 3D accelerated lattice -- header file.
 */
#ifndef ACCELERATED_BLOCK_LATTICE_3D_H
#define ACCELERATED_BLOCK_LATTICE_3D_H

#include "core/globalDefs.h"
#include "multiBlock/multiBlock3D.h"
#include "core/blockLatticeBase3D.h"
#include "atomicBlock/blockLattice3D.h"
#include "acceleratedLattice/acceleratedDefinitions.h"
#include "core/blockStatistics.h"
#include "core/cell.h"
#include "core/dynamics.h"
#include <vector>

namespace plb {

template <typename T, template <typename U> class Descriptor>
class AtomicAcceleratedLattice3D;

/// A complex LatticeBase, itself decomposed into smaller components.
/** This extensible class can be used for example for cache-optimized
 * lattices, irregular domains (no memory allocation in areas exterior to
 * the domain) and parallel lattices. The actual behavior of the lattice
 * is parametrizable by a multiBlockHandler instance, which is given to
 * the constructor.
 *
 * The AcceleratedLattice does not itself possess LatticeProcessors. The Lattice-
 * Processors are delegated to the respective LatticeBases.
 */
template <typename T, template <typename U> class Descriptor>
class AcceleratedLattice3D : public MultiBlock3D {
public:
    typedef std::map<plint, AtomicAcceleratedLattice3D<T, Descriptor> *> BlockMap;

public:
    virtual bool hasRawCommunication() const
    {
        return true;
    }
    ExecutionMode getExecutionMode() const
    {
        return executionMode;
    }
    void setExecutionMode(ExecutionMode executionMode_);
    AcceleratedLattice3D(
        MultiBlockManagement3D const &multiBlockManagement, BlockCommunicator3D *blockCommunicator_,
        CombinedStatistics *combinedStatistics_, Dynamics<T, Descriptor> *backgroundDynamics_);
    AcceleratedLattice3D(
        plint nx, plint ny, plint nz, Dynamics<T, Descriptor> *backgroundDynamics_);
    ~AcceleratedLattice3D();
    AcceleratedLattice3D(AcceleratedLattice3D<T, Descriptor> const &rhs);
    AcceleratedLattice3D(MultiBlockLattice3D<T, Descriptor> const &rhs);
    void writeBack(MultiBlockLattice3D<T, Descriptor> &rhs);
    AcceleratedLattice3D(MultiBlock3D const &rhs);
    virtual AcceleratedLattice3D<T, Descriptor> *clone() const;
    virtual AcceleratedLattice3D<T, Descriptor> *clone(
        MultiBlockManagement3D const &newManagement) const;
    /// Extract sub-domain from rhs and construct a multi-block-lattice with the same
    ///  data distribution and policy-classes; but the data itself and the data-processors
    ///  are not copied. MultiCellAccess takes default value.
    AcceleratedLattice3D(MultiBlock3D const &rhs, Box3D subDomain, bool crop = true);
    /// Attention: data-processors of rhs, which were pointing at rhs, will continue pointing
    /// to rhs, and not to *this.
    void swap(AcceleratedLattice3D &rhs);
    /// Attention: data-processors of rhs, which were pointing at rhs, will continue pointing
    /// to rhs, and not to *this.
    AcceleratedLattice3D<T, Descriptor> &operator=(AcceleratedLattice3D<T, Descriptor> const &rhs);
    // Assign an individual clone of the new dynamics to every cell.
    void resetDynamics(Dynamics<T, Descriptor> const &dynamics);

    Dynamics<T, Descriptor> const &getBackgroundDynamics() const;
    template <class CollFun>
    void collideAndStream(CollFun const &collFun);
    void collideAndStream();
    template <class CollFun>
    void externalCollideAndStream(CollFun const &collFun);
    void incrementTime();
    void resetTime(pluint value);
    virtual AtomicAcceleratedLattice3D<T, Descriptor> &getComponent(plint blockId);
    virtual AtomicAcceleratedLattice3D<T, Descriptor> const &getComponent(plint blockId) const;
    virtual plint sizeOfCell() const;
    virtual plint getCellDim() const;
    virtual int getStaticId() const;
    virtual void copyReceive(
        MultiBlock3D const &fromBlock, Box3D const &fromDomain, Box3D const &toDomain,
        modif::ModifT whichData = modif::dataStructure);

public:
    BlockMap &getBlockLattices();
    BlockMap const &getBlockLattices() const;
    virtual void getDynamicsDict(Box3D domain, std::map<std::string, int> &dict);
    virtual std::string getBlockName() const;
    virtual std::vector<std::string> getTypeInfo() const;
    static std::string blockName();
    static std::string basicType();
    static std::string descriptorType();

private:
    template <class CollFun>
    void collideAndStreamImplementation(CollFun const &collFun);
    void collideAndStreamImplementation();
    void allocateAndInitialize();

private:
    Dynamics<T, Descriptor> *backgroundDynamics;
    BlockMap atomicLattices;
    plint timeCounter;
    ExecutionMode executionMode = ExecutionMode::stdpar;

public:
    static const int staticId;
};

template <typename T, template <typename U> class Descriptor>
AcceleratedLattice3D<T, Descriptor> &findAcceleratedLattice3D(id_t id);

}  // namespace plb

#endif  // ACCELERATED_BLOCK_LATTICE_3D_H
