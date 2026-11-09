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
#ifndef ACCELERATED_COLLISIONS_H
#define ACCELERATED_COLLISIONS_H

#include "core/globalDefs.h"
#include "core/dynamics.h"
#include "core/dynamics.hh"
#include "core/dynamicsIdentifiers.h"
#include "core/dynamicsIdentifiers.hh"
#include "acceleratedLattice/acceleratedCollisionModels.h"
#include "acceleratedLattice/acceleratedCollisionsBoundary.h"
#include <algorithm>

/// All Palabos code is contained in this namespace.
namespace plb {

inline int toCollisionModel(std::string modelName)
{
    auto result = CollisionModel::stringIDs().find(modelName);
    if (result == CollisionModel::stringIDs().end()) {
        pcout << "ERROR: Model \"" << modelName << "\" not implemented." << std::endl;
    }
    PLB_ASSERT(result != CollisionModel::stringIDs().end());
    return result->second;
}

inline std::string toCollisionString(int collisionModel)
{
    std::string result = "";
    for (auto it = CollisionModel::stringIDs().begin(); it != CollisionModel::stringIDs().end();
         ++it) {
        if (it->second == collisionModel) {
            result = it->first;
            break;
        }
    }
    return result;
}

template <typename T, template <typename U> class Descriptor>
int toCollisionModel(Dynamics<T, Descriptor> const &dynamics)
{
    std::vector<int> dynamicsChain;
    constructIdChain(dynamics, dynamicsChain);
    std::string dynamicsName = meta::constructIdNameChain<T, Descriptor>(dynamicsChain, " >> ");
    return toCollisionModel(dynamicsName);
}

// This version of the "collide" function contains a switch including all possible
// collision models. This can lead to a GPU kernel with quite a lot of code. As
// an alternative, you can call the "static switch" version below using variadic
// templates, in which the needed collision models are enumerated explicitly.
template <typename T, template <typename U> class Descriptor>
void collide(
    int collisionModel, Array<T, Descriptor<T>::numPop> &f,
    Array<T, Descriptor<T>::ExternalField::numScalars> &ext,
    Array<T, GPUconst<T, Descriptor>::maxStaticScalars> staticScalars, T *dynamicScalars,
    plint index)
{
    switch (collisionModel) {
    case CollisionModel::NoDynamics:
        Collision<T, Descriptor, CollisionModel::NoDynamics>::collide(
            f, ext, staticScalars, dynamicScalars, index);
        break;
    case CollisionModel::BounceBack:
        Collision<T, Descriptor, CollisionModel::BounceBack>::collide(
            f, ext, staticScalars, dynamicScalars, index);
        break;
    case CollisionModel::BGK:
        Collision<T, Descriptor, CollisionModel::BGK>::collide(
            f, ext, staticScalars, dynamicScalars, index);
        break;
    case CollisionModel::BGK_ExternalMoment:
        Collision<T, Descriptor, CollisionModel::BGK_ExternalMoment>::collide(
            f, ext, staticScalars, dynamicScalars, index);
        break;
    case CollisionModel::TRT:
        Collision<T, Descriptor, CollisionModel::TRT>::collide(
            f, ext, staticScalars, dynamicScalars, index);
        break;
    case CollisionModel::RM:
        Collision<T, Descriptor, CollisionModel::RM>::collide(
            f, ext, staticScalars, dynamicScalars, index);
        break;
    case CollisionModel::HM:
        Collision<T, Descriptor, CollisionModel::HM>::collide(
            f, ext, staticScalars, dynamicScalars, index);
        break;
    case CollisionModel::CM:
        Collision<T, Descriptor, CollisionModel::CM>::collide(
            f, ext, staticScalars, dynamicScalars, index);
        break;
    case CollisionModel::CHM:
        Collision<T, Descriptor, CollisionModel::CHM>::collide(
            f, ext, staticScalars, dynamicScalars, index);
        break;
    case CollisionModel::K:
        Collision<T, Descriptor, CollisionModel::K>::collide(
            f, ext, staticScalars, dynamicScalars, index);
        break;
    case CollisionModel::GH:
        Collision<T, Descriptor, CollisionModel::GH>::collide(
            f, ext, staticScalars, dynamicScalars, index);
        break;
    case CollisionModel::RR:
        Collision<T, Descriptor, CollisionModel::RR>::collide(
            f, ext, staticScalars, dynamicScalars, index);
        break;
    case CollisionModel::HalfwayBounceBack__TRT:
        Collision<T, Descriptor, CollisionModel::HalfwayBounceBack__TRT>::collide(
            f, ext, staticScalars, dynamicScalars, index);
        break;
    case CollisionModel::HalfwayBounceBack__BGK:
        Collision<T, Descriptor, CollisionModel::HalfwayBounceBack__BGK>::collide(
            f, ext, staticScalars, dynamicScalars, index);
        break;
    case CollisionModel::Boundary_RegularizedVelocity_0_1__TRT:
        Collision<T, Descriptor, CollisionModel::Boundary_RegularizedVelocity_0_1__TRT>::collide(
            f, ext, staticScalars, dynamicScalars, index);
        break;
    case CollisionModel::Boundary_RegularizedVelocity_0_M1__TRT:
        Collision<T, Descriptor, CollisionModel::Boundary_RegularizedVelocity_0_M1__TRT>::collide(
            f, ext, staticScalars, dynamicScalars, index);
        break;
    case CollisionModel::Boundary_RegularizedVelocity_1_1__TRT:
        Collision<T, Descriptor, CollisionModel::Boundary_RegularizedVelocity_1_1__TRT>::collide(
            f, ext, staticScalars, dynamicScalars, index);
        break;
    case CollisionModel::Boundary_RegularizedVelocity_1_M1__TRT:
        Collision<T, Descriptor, CollisionModel::Boundary_RegularizedVelocity_1_M1__TRT>::collide(
            f, ext, staticScalars, dynamicScalars, index);
        break;
    case CollisionModel::Boundary_RegularizedVelocity_2_1__TRT:
        Collision<T, Descriptor, CollisionModel::Boundary_RegularizedVelocity_2_1__TRT>::collide(
            f, ext, staticScalars, dynamicScalars, index);
        break;
    case CollisionModel::Boundary_RegularizedVelocity_2_M1__TRT:
        Collision<T, Descriptor, CollisionModel::Boundary_RegularizedVelocity_2_M1__TRT>::collide(
            f, ext, staticScalars, dynamicScalars, index);
        break;
    case CollisionModel::Boundary_RegularizedVelocity_0_1__BGK:
        Collision<T, Descriptor, CollisionModel::Boundary_RegularizedVelocity_0_1__BGK>::collide(
            f, ext, staticScalars, dynamicScalars, index);
        break;
    case CollisionModel::Boundary_RegularizedVelocity_0_M1__BGK:
        Collision<T, Descriptor, CollisionModel::Boundary_RegularizedVelocity_0_M1__BGK>::collide(
            f, ext, staticScalars, dynamicScalars, index);
        break;
    case CollisionModel::Boundary_RegularizedVelocity_1_1__BGK:
        Collision<T, Descriptor, CollisionModel::Boundary_RegularizedVelocity_1_1__BGK>::collide(
            f, ext, staticScalars, dynamicScalars, index);
        break;
    case CollisionModel::Boundary_RegularizedVelocity_1_M1__BGK:
        Collision<T, Descriptor, CollisionModel::Boundary_RegularizedVelocity_1_M1__BGK>::collide(
            f, ext, staticScalars, dynamicScalars, index);
        break;
    case CollisionModel::Boundary_RegularizedVelocity_2_1__BGK:
        Collision<T, Descriptor, CollisionModel::Boundary_RegularizedVelocity_2_1__BGK>::collide(
            f, ext, staticScalars, dynamicScalars, index);
        break;
    case CollisionModel::Boundary_RegularizedVelocity_2_M1__BGK:
        Collision<T, Descriptor, CollisionModel::Boundary_RegularizedVelocity_2_M1__BGK>::collide(
            f, ext, staticScalars, dynamicScalars, index);
        break;
    default:
        PLB_ASSERT(false);
    }
}

template <int...>
struct IntList { };

// Static switch technique from
// https://stackoverflow.com/questions/25202250/c-template-instantiation-avoiding-long-switches
template <typename T, template <typename U> class Descriptor>
struct StaticSwitch {
    static void static_switch(
        IntList<>, int collisionModel, Array<T, Descriptor<T>::numPop> &f,
        Array<T, Descriptor<T>::ExternalField::numScalars> &ext,
        Array<T, GPUconst<T, Descriptor>::maxStaticScalars> staticScalars, T *dynamicScalars,
        plint index)
    {
        printf("Collision model not implemented: %d\n", collisionModel);
    }

    template <int MODEL, int... N>
    static void static_switch(
        IntList<MODEL, N...>, int collisionModel, Array<T, Descriptor<T>::numPop> &f,
        Array<T, Descriptor<T>::ExternalField::numScalars> &ext,
        Array<T, GPUconst<T, Descriptor>::maxStaticScalars> staticScalars, T *dynamicScalars,
        plint index)
    {
        if (collisionModel == MODEL) {
            Collision<T, Descriptor, MODEL>::collide(f, ext, staticScalars, dynamicScalars, index);
        } else {
            static_switch(
                IntList<N...>(), collisionModel, f, ext, staticScalars, dynamicScalars, index);
        }
    }

    template <int MODEL>
    static void static_switch(
        IntList<MODEL>, int collisionModel, Array<T, Descriptor<T>::numPop> &f,
        Array<T, Descriptor<T>::ExternalField::numScalars> &ext,
        Array<T, GPUconst<T, Descriptor>::maxStaticScalars> staticScalars, T *dynamicScalars,
        plint index)
    {
        // PLB_ASSERT (collisionModel == MODEL);
        if (collisionModel == MODEL) {
            Collision<T, Descriptor, MODEL>::collide(f, ext, staticScalars, dynamicScalars, index);
        } else {
            printf("Collision model not implemented: %d\n", collisionModel);
        }
    }

    template <int... N>
    static void static_switch(
        int collisionModel, Array<T, Descriptor<T>::numPop> &f,
        Array<T, Descriptor<T>::ExternalField::numScalars> &ext,
        Array<T, GPUconst<T, Descriptor>::maxStaticScalars> staticScalars, T *dynamicScalars,
        plint index)
    {
        static_switch(
            IntList<N...>(), collisionModel, f, ext, staticScalars, dynamicScalars, index);
    }
};

template <typename T, template <typename U> class Descriptor, int... N>
struct CollisionKernel {
    void operator()(
        int collisionModel, Array<T, Descriptor<T>::numPop> &f,
        Array<T, Descriptor<T>::ExternalField::numScalars> &ext,
        Array<T, GPUconst<T, Descriptor>::maxStaticScalars> staticScalars, T *dynamicScalars,
        plint index) const
    {
        StaticSwitch<T, Descriptor>::template static_switch<N...>(
            collisionModel, f, ext, staticScalars, dynamicScalars, index);
    }
};

}  // namespace plb

#endif  // ACCELERATED_COLLISIONS_H
