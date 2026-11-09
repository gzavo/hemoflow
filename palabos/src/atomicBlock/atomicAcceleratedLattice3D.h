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
 * The dynamics of a 3D block lattice -- header file.
 */
#ifndef ATOMIC_ACCELERATED_LATTICE_3D_H
#define ATOMIC_ACCELERATED_LATTICE_3D_H

#include "core/globalDefs.h"
#include "acceleratedLattice/acceleratedDefinitions.h"
#include "core/plbDebug.h"
#include "core/cell.h"
#include "atomicBlock/dataField3D.h"
#include "core/blockLatticeBase3D.h"
#include "atomicBlock/atomicBlock3D.h"
#include "core/blockIdentifiers.h"
#include <vector>
#include <map>
#include <cstdint>

/// All Palabos code is contained in this namespace.
namespace plb {

template <typename T, template <typename U> class Descriptor>
struct Dynamics;
template <typename T, template <typename U> class Descriptor>
class AtomicAcceleratedLattice3D;

template <typename T, template <typename U> class Descriptor>
class AcceleratedLatticeDataTransfer3D : public BlockDataTransfer3D {
public:
    AcceleratedLatticeDataTransfer3D();
    virtual void setBlock(AtomicBlock3D &block);
    virtual void setConstBlock(AtomicBlock3D const &block);
    virtual AcceleratedLatticeDataTransfer3D<T, Descriptor> *clone() const;
    virtual plint staticCellSize() const;
    /// Send data from the lattice into a byte-stream.
    virtual void send(Box3D domain, std::vector<char> &buffer, modif::ModifT kind) const;
    virtual void send_raw(
        Box3D destDomain, Box3D destBoundingBox, plint deltaX, plint deltaY, plint deltaZ,
        char **buffer, plint &bufferSize, char **indices, modif::ModifT kind) const;
    /// Receive data from a byte-stream into the lattice.
    virtual void receive(Box3D domain, std::vector<char> const &buffer, modif::ModifT kind);
    virtual void receive(
        Box3D domain, std::vector<char> const &buffer, modif::ModifT kind, Dot3D absoluteOffset)
    {
        receive(domain, buffer, kind);
    }
    /// Receive data from a byte-stream into the block, and re-map IDs for dynamics if exist.
    virtual void receive(
        Box3D domain, std::vector<char> const &buffer, modif::ModifT kind,
        std::map<int, std::string> const &foreignIds);
    virtual void receive_raw(
        Box3D domain, char const *buffer, plint bufferSize, char **indices, modif::ModifT kind);
    virtual void receive_raw(
        Box3D domain, char const *buffer, plint bufferSize, char **indices, modif::ModifT kind,
        Dot3D absoluteOffset)
    {
        receive_raw(domain, buffer, bufferSize, indices, kind);
    }
    /// Attribute data between two lattices.
    virtual void attribute(
        Box3D toDomain, plint deltaX, plint deltaY, plint deltaZ, AtomicBlock3D const &from,
        modif::ModifT kind);
    virtual void attribute(
        Box3D toDomain, plint deltaX, plint deltaY, plint deltaZ, AtomicBlock3D const &from,
        modif::ModifT kind, Dot3D absoluteOffset)
    {
        attribute(toDomain, deltaX, deltaY, deltaZ, from, kind);
    }

private:
    void send_static(Box3D domain, std::vector<char> &buffer) const;
    void send_static_raw(
        Box3D destDomain, Box3D destBoundingBox, plint deltaX, plint deltaY, plint deltaZ,
        char **buffer, plint &bufferSize, char **indices) const;
    void send_dynamic(Box3D domain, std::vector<char> &buffer) const;
    void send_all(Box3D domain, std::vector<char> &buffer) const;

    void receive_static(Box3D domain, std::vector<char> const &buffer);
    void receive_static_raw(Box3D domain, char const *buffer, plint bufferSize, char **indices);
    void receive_dynamic(Box3D domain, std::vector<char> const &buffer);
    void receive_all(Box3D domain, std::vector<char> const &buffer);
    void receive_regenerate(
        Box3D domain, std::vector<char> const &buffer,
        std::map<int, int> const &idIndirect = (std::map<int, int>()));

    void attribute_static(
        Box3D toDomain, plint deltaX, plint deltaY, plint deltaZ,
        AtomicAcceleratedLattice3D<T, Descriptor> const &from);
    void attribute_dynamic(
        Box3D toDomain, plint deltaX, plint deltaY, plint deltaZ,
        AtomicAcceleratedLattice3D<T, Descriptor> const &from);
    void attribute_all(
        Box3D toDomain, plint deltaX, plint deltaY, plint deltaZ,
        AtomicAcceleratedLattice3D<T, Descriptor> const &from);
    void attribute_regenerate(
        Box3D toDomain, plint deltaX, plint deltaY, plint deltaZ,
        AtomicAcceleratedLattice3D<T, Descriptor> const &from);

private:
    AtomicAcceleratedLattice3D<T, Descriptor> *lattice;
    AtomicAcceleratedLattice3D<T, Descriptor> const *constLattice;
    template <typename T_, template <typename U_> class Descriptor_>
    friend class ExternalRhoJcollideAndStream3D;
    template <typename T_, template <typename U_> class Descriptor_>
    friend class ExternalCollideAndStream3D;
    template <typename T_, template <typename U_> class Descriptor_>
    friend class WaveAbsorptionExternalRhoJcollideAndStream3D;
};

/// A regular lattice for highly efficient 3D LB dynamics.
/** A block lattice contains a regular array of Cell objects and
 * some useful methods to execute the LB dynamics on the lattice.
 *
 * This class is not intended to be derived from.
 */
template <typename T, template <typename U> class Descriptor>
class AtomicAcceleratedLattice3D : public AtomicBlock3D {
public:
    template <class Functional>
    void for_each(Functional functional);
    template <class Functional>
    void for_each(Box3D const &domain, Functional functional);
    template <class Functional, class BinaryReductionOp>
    T transform_reduce(
        Box3D const &domain, T neutral, BinaryReductionOp reduce, Functional functional);

public:
    /// Construction of an nx_ by ny_ by nz_ lattice
    AtomicAcceleratedLattice3D(
        plint nx_, plint ny_, plint nz_, Dynamics<T, Descriptor> *backgroundDynamics_);
    /// Destruction of the lattice
    ~AtomicAcceleratedLattice3D();
    /// Copy construction
    AtomicAcceleratedLattice3D(AtomicAcceleratedLattice3D<T, Descriptor> const &rhs);
    AtomicAcceleratedLattice3D(BlockLattice3D<T, Descriptor> const &rhs);
    void writeBack(BlockLattice3D<T, Descriptor> &rhs);
    /// Copy assignment
    AtomicAcceleratedLattice3D &operator=(AtomicAcceleratedLattice3D<T, Descriptor> const &rhs);
    /// Swap the content of two BlockLattices
    void swap(AtomicAcceleratedLattice3D &rhs);

public:
    plint getN() const;
    int const *getCollisionMatrix() const;
    template <class CollFun>
    void collideAndStream(CollFun const &collFun);
    void collideAndStream();
    void incrementTime();
    void resetTime(plint timeCounter_);
    // Cell reconstruction can be used when using the accelerated lattice on CPU. On GPU,
    // this kind of pointer-to-dynamics operations are not supported.
    void reconstructCell(plint iX, plint iY, plint iZ, Cell<T, Descriptor> &cell) const;
    void reconstructCellStatic(plint iX, plint iY, plint iZ, Cell<T, Descriptor> &cell) const;
    void pullPop(plint i, Array<T, Descriptor<T>::q> &f) const;
    void pullExt(
        plint iX, plint iY, plint iZ,
        Array<T, Descriptor<T>::ExternalField::numScalars> &ext) const;
    T pullExt(plint iX, plint iY, plint iZ, int offset) const;
    T pullExt(plint i, int offset) const;
    void pushStatic(plint iX, plint iY, plint iZ, Cell<T, Descriptor> const &cell);
    void pushPop(plint i, Array<T, Descriptor<T>::q> const &f);
    void pushPop(plint iX, plint iY, plint iZ, Array<T, Descriptor<T>::q> const &f);
    void pushExt(plint iX, plint iY, plint iZ, int offset, T value);
    void pushExt(plint i, int offset, T value);
    void getDynamicScalar(plint i, T *&dynamicScalarsPtr, plint &index);

public:
    /// Attribute dynamics to a cell.
    void attributeDynamics(plint iX, plint iY, plint iZ, Dynamics<T, Descriptor> *dynamics);
    /// Get a reference to the background dynamics
    Dynamics<T, Descriptor> &getBackgroundDynamics();
    /// Get a const reference to the background dynamics
    Dynamics<T, Descriptor> const &getBackgroundDynamics() const;
    /// Assign an individual clone of the new dynamics to every cell.
    void resetDynamics(Dynamics<T, Descriptor> const &dynamics);
    Dynamics<T, Descriptor> const &getDynamics(plint iX, plint iY, plint iZ) const;
    Dynamics<T, Descriptor> &getDynamics(plint iX, plint iY, plint iZ);

private:
    /// Helper method for memory allocation
    void allocateAndInitialize();
    /// Helper method for memory de-allocation
    void releaseMemory();
    plint allocatedMemory() const;

    // Set the bit at the given position to 1
    inline void setBit(uint32_t &mask, uint8_t position)
    {
        mask |= (1ULL << position);
    }

    // Reset the bit at the given position to 0
    inline void resetBit(uint32_t &mask, uint8_t position)
    {
        mask &= ~(1ULL << position);
    }

    // Read the bit at the given position and return it as a bool
    inline bool readBit(const uint32_t &mask, uint8_t position)
    {
        return (mask >> position) & 1ULL;
    }

private:
    template <typename T_, template <typename U_> class Descriptor_>
    friend class AcceleratedLatticeDataTransfer3D;
    typedef ExternalFieldArray<T, typename Descriptor<T>::ExternalField> External;
    Dynamics<T, Descriptor> *backgroundDynamics;
    T *populations, *tmpPopulations;
    T ****populationGrid, ****tmpPopulationGrid;
    External *externalScalars;
    External ***externalScalarGrid;
    Dynamics<T, Descriptor> **dynamicsArray;
    Dynamics<T, Descriptor> ****dynamicsGrid;
    int *collisionMatrix;
    uint32_t *hw_bb_links;
    plint *dynamicScalarIndex;
    std::vector<T> dynamicScalars;
    plint N;
    plint timeCounter;

public:
    static CachePolicy3D &cachePolicy();
};

}  // namespace plb

#endif  // ATOMIC_ACCELERATED_LATTICE_3D_H
