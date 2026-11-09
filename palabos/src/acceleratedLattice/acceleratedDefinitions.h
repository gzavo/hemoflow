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

#ifndef ACCELERATED_DEFINITIONS_H
#define ACCELERATED_DEFINITIONS_H

#include "core/globalDefs.h"
#include "core/array.h"
#include "core/dynamics.h"
#include "basicDynamics/comprehensiveIsoThermalDynamics.h"
#include <string>
#include <map>
#include <algorithm>

namespace plb {

enum class ExecutionMode { seq, stdpar, openmp };

template <typename T, template <typename U> class Descriptor>
struct GPUconst {
    enum { maxStaticScalars = Descriptor<T>::numRelaxationTimes + 2 };
};

namespace CollisionModel {
enum {
    NoDynamics,
    BounceBack,
    BGK,
    BGK_ExternalMoment,
    TRT,
    RM,
    HM,
    CM,
    CHM,
    K,
    GH,
    RR,
    Boundary_RegularizedVelocity_0_1__TRT,
    Boundary_RegularizedVelocity_0_M1__TRT,
    Boundary_RegularizedVelocity_1_1__TRT,
    Boundary_RegularizedVelocity_1_M1__TRT,
    Boundary_RegularizedVelocity_2_1__TRT,
    Boundary_RegularizedVelocity_2_M1__TRT,
    Boundary_RegularizedVelocity_0_1__BGK,
    Boundary_RegularizedVelocity_0_M1__BGK,
    Boundary_RegularizedVelocity_1_1__BGK,
    Boundary_RegularizedVelocity_1_M1__BGK,
    Boundary_RegularizedVelocity_2_1__BGK,
    Boundary_RegularizedVelocity_2_M1__BGK,
    Boundary_RegularizedDensity_0_1__TRT,
    Boundary_RegularizedDensity_0_M1__TRT,
    Boundary_RegularizedDensity_1_1__TRT,
    Boundary_RegularizedDensity_1_M1__TRT,
    Boundary_RegularizedDensity_2_1__TRT,
    Boundary_RegularizedDensity_2_M1__TRT,
    // Important: all half-way bounce-back models must be listed at the
    // end, and start with HalfwayBounceBack__TRT: the streaming step,
    // which is modified in face of half-way bounce-back, depends on this.
    HalfwayBounceBack__TRT,
    HalfwayBounceBack__BGK,
    HalfwayBounceBack__RM,
    HalfwayBounceBack__HM,
    HalfwayBounceBack__CM,
    HalfwayBounceBack__CHM,
    HalfwayBounceBack__K,
    HalfwayBounceBack__GH,
    HalfwayBounceBack__RR

};

inline std::map<std::string, int> const &stringIDs()
{
    static const std::map<std::string, int> stringIDmap {
        {"NoDynamics", NoDynamics},
        {"BounceBack", BounceBack},
        {"BGK", BGK},
        {"BGK_ExternalMoment", BGK_ExternalMoment},
        {"TRT", TRT},
        {"RM", RM},
        {"HM", HM},
        {"CM", CM},
        {"CHM", CHM},
        {"K", K},
        {"GH", GH},
        {"RR", RR},
        {"Boundary_RegularizedVelocity_0_1 >> TRT", Boundary_RegularizedVelocity_0_1__TRT},
        {"Boundary_RegularizedVelocity_0_-1 >> TRT", Boundary_RegularizedVelocity_0_M1__TRT},
        {"Boundary_RegularizedVelocity_1_1 >> TRT", Boundary_RegularizedVelocity_1_1__TRT},
        {"Boundary_RegularizedVelocity_1_-1 >> TRT", Boundary_RegularizedVelocity_1_M1__TRT},
        {"Boundary_RegularizedVelocity_2_1 >> TRT", Boundary_RegularizedVelocity_2_1__TRT},
        {"Boundary_RegularizedVelocity_2_-1 >> TRT", Boundary_RegularizedVelocity_2_M1__TRT},
        {"Boundary_RegularizedVelocity_0_1 >> BGK", Boundary_RegularizedVelocity_0_1__BGK},
        {"Boundary_RegularizedVelocity_0_-1 >> BGK", Boundary_RegularizedVelocity_0_M1__BGK},
        {"Boundary_RegularizedVelocity_1_1 >> BGK", Boundary_RegularizedVelocity_1_1__BGK},
        {"Boundary_RegularizedVelocity_1_-1 >> BGK", Boundary_RegularizedVelocity_1_M1__BGK},
        {"Boundary_RegularizedVelocity_2_1 >> BGK", Boundary_RegularizedVelocity_2_1__BGK},
        {"Boundary_RegularizedVelocity_2_-1 >> BGK", Boundary_RegularizedVelocity_2_M1__BGK},
        {"Boundary_RegularizedDensity_0_1 >> TRT", Boundary_RegularizedDensity_0_1__TRT},
        {"Boundary_RegularizedDensity_0_-1 >> TRT", Boundary_RegularizedDensity_0_M1__TRT},
        {"Boundary_RegularizedDensity_1_1 >> TRT", Boundary_RegularizedDensity_1_1__TRT},
        {"Boundary_RegularizedDensity_1_-1 >> TRT", Boundary_RegularizedDensity_1_M1__TRT},
        {"Boundary_RegularizedDensity_2_1 >> TRT", Boundary_RegularizedDensity_2_1__TRT},
        {"Boundary_RegularizedDensity_2_-1 >> TRT", Boundary_RegularizedDensity_2_M1__TRT},
        {"HalfwayBounceBack >> TRT", HalfwayBounceBack__TRT},
        {"HalfwayBounceBack >> BGK", HalfwayBounceBack__BGK},
        {"HalfwayBounceBack >> RM", HalfwayBounceBack__RM},
        {"HalfwayBounceBack >> HM", HalfwayBounceBack__HM},
        {"HalfwayBounceBack >> CM", HalfwayBounceBack__CM},
        {"HalfwayBounceBack >> CHM", HalfwayBounceBack__CHM},
        {"HalfwayBounceBack >> K", HalfwayBounceBack__K},
        {"HalfwayBounceBack >> GH", HalfwayBounceBack__GH},
        {"HalfwayBounceBack >> RR", HalfwayBounceBack__RR}};

    return stringIDmap;
}

}  // namespace CollisionModel

template <typename T, template <typename U> class Descriptor>
int numDynamicScalars(int collisionModel)
{
    switch (collisionModel) {
    case CollisionModel::BounceBack:
        return 1;
    case CollisionModel::Boundary_RegularizedVelocity_0_1__TRT:
    case CollisionModel::Boundary_RegularizedVelocity_0_M1__TRT:
    case CollisionModel::Boundary_RegularizedVelocity_1_1__TRT:
    case CollisionModel::Boundary_RegularizedVelocity_1_M1__TRT:
    case CollisionModel::Boundary_RegularizedVelocity_2_1__TRT:
    case CollisionModel::Boundary_RegularizedVelocity_2_M1__TRT:
    case CollisionModel::Boundary_RegularizedVelocity_0_1__BGK:
    case CollisionModel::Boundary_RegularizedVelocity_0_M1__BGK:
    case CollisionModel::Boundary_RegularizedVelocity_1_1__BGK:
    case CollisionModel::Boundary_RegularizedVelocity_1_M1__BGK:
    case CollisionModel::Boundary_RegularizedVelocity_2_1__BGK:
    case CollisionModel::Boundary_RegularizedVelocity_2_M1__BGK:
        return 3;
    case CollisionModel::Boundary_RegularizedDensity_0_1__TRT:
    case CollisionModel::Boundary_RegularizedDensity_0_M1__TRT:
    case CollisionModel::Boundary_RegularizedDensity_1_1__TRT:
    case CollisionModel::Boundary_RegularizedDensity_1_M1__TRT:
    case CollisionModel::Boundary_RegularizedDensity_2_1__TRT:
    case CollisionModel::Boundary_RegularizedDensity_2_M1__TRT:
        return 1;
    case CollisionModel::HalfwayBounceBack__TRT:
    case CollisionModel::HalfwayBounceBack__BGK:
    case CollisionModel::HalfwayBounceBack__RM:
    case CollisionModel::HalfwayBounceBack__HM:
    case CollisionModel::HalfwayBounceBack__CM:
    case CollisionModel::HalfwayBounceBack__CHM:
    case CollisionModel::HalfwayBounceBack__K:
    case CollisionModel::HalfwayBounceBack__GH:
    case CollisionModel::HalfwayBounceBack__RR:
        return Descriptor<T>::q;
    default:
        return 0;
    }
}

template <typename T, template <typename U> class Descriptor>
std::vector<T> getDynamicScalars(Dynamics<T, Descriptor> const &dynamics, int collisionModel)
{
    static Cell<T, Descriptor> dummyCell;
    switch (collisionModel) {
    case CollisionModel::BounceBack:
        return std::vector<T> {dynamics.computeDensity(dummyCell)};
    case CollisionModel::Boundary_RegularizedVelocity_0_1__TRT:
    case CollisionModel::Boundary_RegularizedVelocity_0_M1__TRT:
    case CollisionModel::Boundary_RegularizedVelocity_1_1__TRT:
    case CollisionModel::Boundary_RegularizedVelocity_1_M1__TRT:
    case CollisionModel::Boundary_RegularizedVelocity_2_1__TRT:
    case CollisionModel::Boundary_RegularizedVelocity_2_M1__TRT:
    case CollisionModel::Boundary_RegularizedVelocity_0_1__BGK:
    case CollisionModel::Boundary_RegularizedVelocity_0_M1__BGK:
    case CollisionModel::Boundary_RegularizedVelocity_1_1__BGK:
    case CollisionModel::Boundary_RegularizedVelocity_1_M1__BGK:
    case CollisionModel::Boundary_RegularizedVelocity_2_1__BGK:
    case CollisionModel::Boundary_RegularizedVelocity_2_M1__BGK: {
        Array<T, Descriptor<T>::d> u;
        dynamics.computeVelocity(dummyCell, u);
        return std::vector<T> {u[0], u[1], u[2]};
    }
    case CollisionModel::Boundary_RegularizedDensity_0_1__TRT:
    case CollisionModel::Boundary_RegularizedDensity_0_M1__TRT:
    case CollisionModel::Boundary_RegularizedDensity_1_1__TRT:
    case CollisionModel::Boundary_RegularizedDensity_1_M1__TRT:
    case CollisionModel::Boundary_RegularizedDensity_2_1__TRT:
    case CollisionModel::Boundary_RegularizedDensity_2_M1__TRT: {
        return std::vector<T> {dynamics.computeDensity(dummyCell)};
    }
    case CollisionModel::HalfwayBounceBack__TRT:
    case CollisionModel::HalfwayBounceBack__BGK:
    case CollisionModel::HalfwayBounceBack__RM:
    case CollisionModel::HalfwayBounceBack__HM:
    case CollisionModel::HalfwayBounceBack__CM:
    case CollisionModel::HalfwayBounceBack__CHM:
    case CollisionModel::HalfwayBounceBack__K:
    case CollisionModel::HalfwayBounceBack__GH:
    case CollisionModel::HalfwayBounceBack__RR: {
        HalfwayBounceBack<T, Descriptor> const *hwBBdynamics {
            &dynamic_cast<HalfwayBounceBack<T, Descriptor> const &>(dynamics)};
        std::vector<T> bdData(Descriptor<T>::q);
        for (int iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
            bdData[iPop] = hwBBdynamics->getData(iPop);
        }
        return bdData;
    }
    // Uncomment this line if you want BGK to handle omega as a dynamic scalar
    // case CollisionModel::BGK: return std::vector<T>(1, dynamics.getOmega());
    default:
        return std::vector<T>();
    }
}

template <typename T, template <typename U> class Descriptor>
void getStaticScalars(
    Dynamics<T, Descriptor> const &backgroundDynamics, int collisionModel,
    Array<T, GPUconst<T, Descriptor>::maxStaticScalars> &staticScalars)
{
    std::fill(
        &staticScalars[0] + 1, &staticScalars[0] + GPUconst<T, Descriptor>::maxStaticScalars, T());
    staticScalars[0] = backgroundDynamics.getOmega();
    switch (collisionModel) {
    case CollisionModel::TRT:
    case CollisionModel::HalfwayBounceBack__TRT:
        staticScalars[1] = backgroundDynamics.getParameter(dynamicParams::omega_minus);
        break;
    case CollisionModel::RM:
        for (int i = 0; i < Descriptor<T>::numRelaxationTimes; ++i) {
            staticScalars[i] = RMdynamics<T, Descriptor>::allOmega[i];
        }
        break;
    case CollisionModel::HM:
        for (int i = 0; i < Descriptor<T>::numRelaxationTimes; ++i) {
            staticScalars[i] = HMdynamics<T, Descriptor>::allOmega[i];
        }
        break;
    case CollisionModel::CM:
        for (int i = 0; i < Descriptor<T>::numRelaxationTimes; ++i) {
            staticScalars[i] = CMdynamics<T, Descriptor>::allOmega[i];
        }
        break;
    case CollisionModel::CHM:
        for (int i = 0; i < Descriptor<T>::numRelaxationTimes; ++i) {
            staticScalars[i] = CHMdynamics<T, Descriptor>::allOmega[i];
        }
        break;
    case CollisionModel::K:
        for (int i = 0; i < Descriptor<T>::numRelaxationTimes; ++i) {
            staticScalars[i] = Kdynamics<T, Descriptor>::allOmega[i];
        }
        break;
    case CollisionModel::GH:
        for (int i = 0; i < Descriptor<T>::numRelaxationTimes; ++i) {
            staticScalars[i] = GHdynamics<T, Descriptor>::allOmega[i];
        }
        break;
    case CollisionModel::RR:
        for (int i = 0; i < Descriptor<T>::numRelaxationTimes; ++i) {
            staticScalars[i] = RRdynamics<T, Descriptor>::allOmega[i];
        }
        break;
    default:
        break;
    }
}

template <typename T, template <typename U> class Descriptor, int model>
struct Collision {
    static void collide(
        Array<T, Descriptor<T>::numPop> &f, Array<T, Descriptor<T>::ExternalField::numScalars> &ext,
        Array<T, GPUconst<T, Descriptor>::maxStaticScalars> staticScalars, T *dynamicScalars,
        plint index)
    {
        printf("Collision model is not implemented: %d\n", model);
        throw -1;
    }
};

}  // namespace plb

#endif
