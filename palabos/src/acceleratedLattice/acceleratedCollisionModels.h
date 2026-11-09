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
#ifndef ACCELERATED_COLLISION_MODELS_H
#define ACCELERATED_COLLISION_MODELS_H

#include "core/globalDefs.h"
#include "core/dynamics.h"
#include "core/dynamics.hh"
#include "core/dynamicsIdentifiers.h"
#include "core/dynamicsIdentifiers.hh"
#include "latticeBoltzmann/comprehensiveModelsTemplates.h"
#include <map>
#include <algorithm>
#include <cmath>

/// All Palabos code is contained in this namespace.
namespace plb {

template <typename T, template <typename U> class Descriptor>
struct Collision<T, Descriptor, CollisionModel::NoDynamics> {
    static void collide(
        Array<T, Descriptor<T>::numPop> &f, Array<T, Descriptor<T>::ExternalField::numScalars> &ext,
        Array<T, GPUconst<T, Descriptor>::maxStaticScalars> staticScalars, T *dynamicScalars,
        plint index)
    { }
};

template <typename T, template <typename U> class Descriptor>
struct Collision<T, Descriptor, CollisionModel::BounceBack> {
    static void collide(
        Array<T, Descriptor<T>::numPop> &f, Array<T, Descriptor<T>::ExternalField::numScalars> &ext,
        Array<T, GPUconst<T, Descriptor>::maxStaticScalars> staticScalars, T *dynamicScalars,
        plint index)
    {
        for (plint iPop = 1; iPop <= Descriptor<T>::q / 2; ++iPop) {
            std::swap(f[iPop], f[iPop + Descriptor<T>::q / 2]);
        }
    }
};

template <typename T, template <typename U> class Descriptor>
struct Collision<T, Descriptor, CollisionModel::BGK> {
    static void collide(
        Array<T, Descriptor<T>::numPop> &f, Array<T, Descriptor<T>::ExternalField::numScalars> &ext,
        Array<T, GPUconst<T, Descriptor>::maxStaticScalars> staticScalars, T *dynamicScalars,
        plint index)
    {
        T omega = staticScalars[0];
        // Uncomment this line if you want BGK to handle omega as a dynamic scalar
        // T omega = dynamicScalars[index + 1];
        T rhoBar;
        Array<T, Descriptor<T>::d> j;
        momentTemplatesImpl<T, typename Descriptor<T>::BaseDescriptor>::get_rhoBar_j(f, rhoBar, j);
        dynamicsTemplatesImpl<T, typename Descriptor<T>::BaseDescriptor>::bgk_ma2_collision(
            f, rhoBar, j, omega);
    }
};

template <int numScalars, typename T, template <typename U> class Descriptor>
struct CollideExternalMomentBGK {
    static void collide(
        Array<T, Descriptor<T>::numPop> &f, Array<T, Descriptor<T>::ExternalField::numScalars> &ext,
        T omega)
    {
        T rho = ext[Descriptor<T>::ExternalField::densityBeginsAt];
        T rhoBar = Descriptor<T>::rhoBar(rho);
        constexpr int m = Descriptor<T>::ExternalField::momentumBeginsAt;
        Array<T, Descriptor<T>::d> j {ext[m], ext[m + 1], ext[m + 2]};
        dynamicsTemplatesImpl<T, typename Descriptor<T>::BaseDescriptor>::bgk_ma2_collision(
            f, rhoBar, j, omega);
    }
};

template <typename T, template <typename U> class Descriptor>
struct CollideExternalMomentBGK<0, T, Descriptor> {
    static void collide(
        Array<T, Descriptor<T>::numPop> &f, Array<T, Descriptor<T>::ExternalField::numScalars> &ext,
        T omega)
    { }
};

template <typename T, template <typename U> class Descriptor>
struct CollideExternalMomentBGK<1, T, Descriptor> {
    static void collide(
        Array<T, Descriptor<T>::numPop> &f, Array<T, Descriptor<T>::ExternalField::numScalars> &ext,
        T omega)
    { }
};

template <typename T, template <typename U> class Descriptor>
struct CollideExternalMomentBGK<2, T, Descriptor> {
    static void collide(
        Array<T, Descriptor<T>::numPop> &f, Array<T, Descriptor<T>::ExternalField::numScalars> &ext,
        T omega)
    { }
};

template <typename T, template <typename U> class Descriptor>
struct CollideExternalMomentBGK<3, T, Descriptor> {
    static void collide(
        Array<T, Descriptor<T>::numPop> &f, Array<T, Descriptor<T>::ExternalField::numScalars> &ext,
        T omega)
    { }
};

template <typename T, template <typename U> class Descriptor>
struct Collision<T, Descriptor, CollisionModel::BGK_ExternalMoment> {
    static void collide(
        Array<T, Descriptor<T>::numPop> &f, Array<T, Descriptor<T>::ExternalField::numScalars> &ext,
        Array<T, GPUconst<T, Descriptor>::maxStaticScalars> staticScalars, T *dynamicScalars,
        plint index)
    {
        T omega = staticScalars[0];
        CollideExternalMomentBGK<Descriptor<T>::ExternalField::numScalars, T, Descriptor>::collide(
            f, ext, omega);
    }
};

template <typename T, template <typename U> class Descriptor>
struct Collision<T, Descriptor, CollisionModel::TRT> {
    static void collide(
        Array<T, Descriptor<T>::numPop> &f, Array<T, Descriptor<T>::ExternalField::numScalars> &ext,
        Array<T, GPUconst<T, Descriptor>::maxStaticScalars> staticScalars, T *dynamicScalars,
        plint index)
    {
        const T omegaPlus = staticScalars[0];
        const T omegaMinus = staticScalars[1];

        Array<T, Descriptor<T>::d> j;
        T rhoBar;
        momentTemplatesImpl<T, typename Descriptor<T>::BaseDescriptor>::get_rhoBar_j(f, rhoBar, j);
        dynamicsTemplatesImpl<T, typename Descriptor<T>::BaseDescriptor>::trt_ma2_collision(
            f, rhoBar, j, omegaPlus, omegaMinus);
    }
};

template <typename T, template <typename U> class Descriptor>
struct Collision<T, Descriptor, CollisionModel::RM> {
    static void collide(
        Array<T, Descriptor<T>::numPop> &f, Array<T, Descriptor<T>::ExternalField::numScalars> &ext,
        Array<T, GPUconst<T, Descriptor>::maxStaticScalars> staticScalars, T *dynamicScalars,
        plint index)
    {
        T rho;
        Array<T, Descriptor<T>::q> RM, RMeq;
        comprehensiveDynamicsTemplatesImpl<
            T, typename Descriptor<T>::BaseDescriptor>::RMcomputeMoments(f, RM, rho);
        Array<T, Descriptor<T>::d> u;
        for (int i = 0; i < Descriptor<T>::d; ++i) {
            u[i] = RM[i + 1];
        }
        comprehensiveDynamicsTemplates<T, Descriptor>::RMcomputeEquilibriumMoments(u, RMeq);
        // All relaxation times are static class members, except for the first one
        // (the viscous relaxation time), which is dynamic.
        Array<T, Descriptor<T>::numRelaxationTimes> allOmegaTmp;
        for (int i = 0; i < Descriptor<T>::numRelaxationTimes; ++i) {
            allOmegaTmp[i] = staticScalars[i];
        }
        comprehensiveDynamicsTemplatesImpl<T, typename Descriptor<T>::BaseDescriptor>::RMcollide(
            f, rho, u, RM, RMeq, allOmegaTmp);
    }
};

template <typename T, template <typename U> class Descriptor>
struct Collision<T, Descriptor, CollisionModel::HM> {
    static void collide(
        Array<T, Descriptor<T>::numPop> &f, Array<T, Descriptor<T>::ExternalField::numScalars> &ext,
        Array<T, GPUconst<T, Descriptor>::maxStaticScalars> staticScalars, T *dynamicScalars,
        plint index)
    {
        T rho;
        Array<T, Descriptor<T>::q> HM, HMeq;
        comprehensiveDynamicsTemplatesImpl<
            T, typename Descriptor<T>::BaseDescriptor>::HMcomputeMoments(f, HM, rho);
        Array<T, Descriptor<T>::d> u;
        for (int i = 0; i < Descriptor<T>::d; ++i) {
            u[i] = HM[i + 1];
        }
        comprehensiveDynamicsTemplates<T, Descriptor>::HMcomputeEquilibriumMoments(u, HMeq);
        // All relaxation times are static class members, except for the first one
        // (the viscous relaxation time), which is dynamic.
        Array<T, Descriptor<T>::numRelaxationTimes> allOmegaTmp;
        for (int i = 0; i < Descriptor<T>::numRelaxationTimes; ++i) {
            allOmegaTmp[i] = staticScalars[i];
        }
        comprehensiveDynamicsTemplatesImpl<T, typename Descriptor<T>::BaseDescriptor>::HMcollide(
            f, rho, u, HM, HMeq, allOmegaTmp);
    }
};

template <typename T, template <typename U> class Descriptor>
struct Collision<T, Descriptor, CollisionModel::CM> {
    static void collide(
        Array<T, Descriptor<T>::numPop> &f, Array<T, Descriptor<T>::ExternalField::numScalars> &ext,
        Array<T, GPUconst<T, Descriptor>::maxStaticScalars> staticScalars, T *dynamicScalars,
        plint index)
    {
        T rho;
        Array<T, Descriptor<T>::q> CM, CMeq;
        // Here we cannot use CMs to compute u... That is why it is done during CMcomputeMoments
        Array<T, Descriptor<T>::d> u;
        comprehensiveDynamicsTemplatesImpl<
            T, typename Descriptor<T>::BaseDescriptor>::CMcomputeMoments(f, CM, rho, u);
        comprehensiveDynamicsTemplates<T, Descriptor>::CMcomputeEquilibriumMoments(CMeq);
        // All relaxation times are static class members, except for the first one
        // (the viscous relaxation time), which is dynamic.
        Array<T, Descriptor<T>::numRelaxationTimes> allOmegaTmp;
        for (int i = 0; i < Descriptor<T>::numRelaxationTimes; ++i) {
            allOmegaTmp[i] = staticScalars[i];
        }
        comprehensiveDynamicsTemplatesImpl<T, typename Descriptor<T>::BaseDescriptor>::CMcollide(
            f, rho, u, CM, CMeq, allOmegaTmp);
    }
};

template <typename T, template <typename U> class Descriptor>
struct Collision<T, Descriptor, CollisionModel::CHM> {
    static void collide(
        Array<T, Descriptor<T>::numPop> &f, Array<T, Descriptor<T>::ExternalField::numScalars> &ext,
        Array<T, GPUconst<T, Descriptor>::maxStaticScalars> staticScalars, T *dynamicScalars,
        plint index)
    {
        T rho;
        Array<T, Descriptor<T>::q> CHM, CHMeq;
        // Here we cannot use CHMs to compute u... That is why it is done during CHMcomputeMoments
        Array<T, Descriptor<T>::d> u;
        comprehensiveDynamicsTemplatesImpl<
            T, typename Descriptor<T>::BaseDescriptor>::CHMcomputeMoments(f, CHM, rho, u);
        comprehensiveDynamicsTemplates<T, Descriptor>::CHMcomputeEquilibriumMoments(CHMeq);
        // All relaxation times are static class members, except for the first one
        // (the viscous relaxation time), which is dynamic.
        Array<T, Descriptor<T>::numRelaxationTimes> allOmegaTmp;
        for (int i = 0; i < Descriptor<T>::numRelaxationTimes; ++i) {
            allOmegaTmp[i] = staticScalars[i];
        }
        comprehensiveDynamicsTemplatesImpl<T, typename Descriptor<T>::BaseDescriptor>::CHMcollide(
            f, rho, u, CHM, CHMeq, allOmegaTmp);
    }
};

template <typename T, template <typename U> class Descriptor>
struct Collision<T, Descriptor, CollisionModel::K> {
    static void collide(
        Array<T, Descriptor<T>::numPop> &f, Array<T, Descriptor<T>::ExternalField::numScalars> &ext,
        Array<T, GPUconst<T, Descriptor>::maxStaticScalars> staticScalars, T *dynamicScalars,
        plint index)
    {
        T rho;
        Array<T, Descriptor<T>::q> K, Keq;
        // Here we cannot use Ks to compute u... That is why it is done during KcomputeMoments
        Array<T, Descriptor<T>::d> u;
        comprehensiveDynamicsTemplatesImpl<
            T, typename Descriptor<T>::BaseDescriptor>::KcomputeMoments(f, K, rho, u);
        comprehensiveDynamicsTemplates<T, Descriptor>::KcomputeEquilibriumMoments(u, Keq);
        // All relaxation times are static class members, except for the first one
        // (the viscous relaxation time), which is dynamic.
        Array<T, Descriptor<T>::numRelaxationTimes> allOmegaTmp;
        for (int i = 0; i < Descriptor<T>::numRelaxationTimes; ++i) {
            allOmegaTmp[i] = staticScalars[i];
        }
        comprehensiveDynamicsTemplatesImpl<T, typename Descriptor<T>::BaseDescriptor>::Kcollide(
            f, rho, u, K, Keq, allOmegaTmp);
    }
};

template <typename T, template <typename U> class Descriptor>
struct Collision<T, Descriptor, CollisionModel::GH> {
    static void collide(
        Array<T, Descriptor<T>::numPop> &f, Array<T, Descriptor<T>::ExternalField::numScalars> &ext,
        Array<T, GPUconst<T, Descriptor>::maxStaticScalars> staticScalars, T *dynamicScalars,
        plint index)
    {
        T rho;
        Array<T, Descriptor<T>::q> GH, GHeq;
        comprehensiveDynamicsTemplatesImpl<
            T, typename Descriptor<T>::BaseDescriptor>::GHcomputeMoments(f, GH, rho);
        Array<T, Descriptor<T>::d> u;
        for (int i = 0; i < Descriptor<T>::d; ++i) {
            u[i] = GH[i + 1];
        }
        comprehensiveDynamicsTemplates<T, Descriptor>::GHcomputeEquilibriumMoments(u, GHeq);
        // All relaxation times are static class members, except for the first one
        // (the viscous relaxation time), which is dynamic.
        Array<T, Descriptor<T>::numRelaxationTimes> allOmegaTmp;
        for (int i = 0; i < Descriptor<T>::numRelaxationTimes; ++i) {
            allOmegaTmp[i] = staticScalars[i];
        }
        comprehensiveDynamicsTemplatesImpl<T, typename Descriptor<T>::BaseDescriptor>::GHcollide(
            f, rho, u, GH, GHeq, allOmegaTmp);
    }
};

template <typename T, template <typename U> class Descriptor>
struct Collision<T, Descriptor, CollisionModel::RR> {
    static void collide(
        Array<T, Descriptor<T>::numPop> &f, Array<T, Descriptor<T>::ExternalField::numScalars> &ext,
        Array<T, GPUconst<T, Descriptor>::maxStaticScalars> staticScalars, T *dynamicScalars,
        plint index)
    {
        T rho;
        Array<T, Descriptor<T>::q> RR, RReq;
        comprehensiveDynamicsTemplatesImpl<
            T, typename Descriptor<T>::BaseDescriptor>::RRcomputeMoments(f, RR, rho);
        Array<T, Descriptor<T>::d> u;
        for (int i = 0; i < Descriptor<T>::d; ++i) {
            u[i] = RR[i + 1];
        }
        comprehensiveDynamicsTemplates<T, Descriptor>::RRcomputeEquilibriumMoments(u, RReq);
        // All relaxation times are static class members, except for the first one
        // (the viscous relaxation time), which is dynamic.
        Array<T, Descriptor<T>::numRelaxationTimes> allOmegaTmp;
        for (int i = 0; i < Descriptor<T>::numRelaxationTimes; ++i) {
            allOmegaTmp[i] = staticScalars[i];
        }
        comprehensiveDynamicsTemplatesImpl<T, typename Descriptor<T>::BaseDescriptor>::RRcollide(
            f, rho, u, RR, RReq, allOmegaTmp);
    }
};

template <typename T, template <typename U> class Descriptor>
struct Collision<T, Descriptor, CollisionModel::HalfwayBounceBack__TRT> {
    static void collide(
        Array<T, Descriptor<T>::numPop> &f, Array<T, Descriptor<T>::ExternalField::numScalars> &ext,
        Array<T, GPUconst<T, Descriptor>::maxStaticScalars> staticScalars, T *dynamicScalars,
        plint index)
    {
        Collision<T, Descriptor, CollisionModel::TRT>::collide(
            f, ext, staticScalars, dynamicScalars, index);
        for (int iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
            if (!std::isnan(dynamicScalars[index + iPop])) {
                f[iPop] += dynamicScalars[index + iPop];
            }
        }
    }
};

template <typename T, template <typename U> class Descriptor>
struct Collision<T, Descriptor, CollisionModel::HalfwayBounceBack__BGK> {
    static void collide(
        Array<T, Descriptor<T>::numPop> &f, Array<T, Descriptor<T>::ExternalField::numScalars> &ext,
        Array<T, GPUconst<T, Descriptor>::maxStaticScalars> staticScalars, T *dynamicScalars,
        plint index)
    {
        Collision<T, Descriptor, CollisionModel::BGK>::collide(
            f, ext, staticScalars, dynamicScalars, index);
        for (int iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
            if (!std::isnan(dynamicScalars[index + iPop])) {
                f[iPop] += dynamicScalars[index + iPop];
            }
        }
    }
};

template <typename T, template <typename U> class Descriptor>
struct Collision<T, Descriptor, CollisionModel::HalfwayBounceBack__RM> {
    static void collide(
        Array<T, Descriptor<T>::numPop> &f, Array<T, Descriptor<T>::ExternalField::numScalars> &ext,
        Array<T, GPUconst<T, Descriptor>::maxStaticScalars> staticScalars, T *dynamicScalars,
        plint index)
    {
        Collision<T, Descriptor, CollisionModel::RM>::collide(
            f, ext, staticScalars, dynamicScalars, index);
        for (int iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
            if (!std::isnan(dynamicScalars[index + iPop])) {
                f[iPop] += dynamicScalars[index + iPop];
            }
        }
    }
};

template <typename T, template <typename U> class Descriptor>
struct Collision<T, Descriptor, CollisionModel::HalfwayBounceBack__HM> {
    static void collide(
        Array<T, Descriptor<T>::numPop> &f, Array<T, Descriptor<T>::ExternalField::numScalars> &ext,
        Array<T, GPUconst<T, Descriptor>::maxStaticScalars> staticScalars, T *dynamicScalars,
        plint index)
    {
        Collision<T, Descriptor, CollisionModel::HM>::collide(
            f, ext, staticScalars, dynamicScalars, index);
        for (int iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
            if (!std::isnan(dynamicScalars[index + iPop])) {
                f[iPop] += dynamicScalars[index + iPop];
            }
        }
    }
};

template <typename T, template <typename U> class Descriptor>
struct Collision<T, Descriptor, CollisionModel::HalfwayBounceBack__CM> {
    static void collide(
        Array<T, Descriptor<T>::numPop> &f, Array<T, Descriptor<T>::ExternalField::numScalars> &ext,
        Array<T, GPUconst<T, Descriptor>::maxStaticScalars> staticScalars, T *dynamicScalars,
        plint index)
    {
        Collision<T, Descriptor, CollisionModel::CM>::collide(
            f, ext, staticScalars, dynamicScalars, index);
        for (int iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
            if (!std::isnan(dynamicScalars[index + iPop])) {
                f[iPop] += dynamicScalars[index + iPop];
            }
        }
    }
};

template <typename T, template <typename U> class Descriptor>
struct Collision<T, Descriptor, CollisionModel::HalfwayBounceBack__CHM> {
    static void collide(
        Array<T, Descriptor<T>::numPop> &f, Array<T, Descriptor<T>::ExternalField::numScalars> &ext,
        Array<T, GPUconst<T, Descriptor>::maxStaticScalars> staticScalars, T *dynamicScalars,
        plint index)
    {
        Collision<T, Descriptor, CollisionModel::CHM>::collide(
            f, ext, staticScalars, dynamicScalars, index);
        for (int iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
            if (!std::isnan(dynamicScalars[index + iPop])) {
                f[iPop] += dynamicScalars[index + iPop];
            }
        }
    }
};

template <typename T, template <typename U> class Descriptor>
struct Collision<T, Descriptor, CollisionModel::HalfwayBounceBack__K> {
    static void collide(
        Array<T, Descriptor<T>::numPop> &f, Array<T, Descriptor<T>::ExternalField::numScalars> &ext,
        Array<T, GPUconst<T, Descriptor>::maxStaticScalars> staticScalars, T *dynamicScalars,
        plint index)
    {
        Collision<T, Descriptor, CollisionModel::K>::collide(
            f, ext, staticScalars, dynamicScalars, index);
        for (int iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
            if (!std::isnan(dynamicScalars[index + iPop])) {
                f[iPop] += dynamicScalars[index + iPop];
            }
        }
    }
};

template <typename T, template <typename U> class Descriptor>
struct Collision<T, Descriptor, CollisionModel::HalfwayBounceBack__GH> {
    static void collide(
        Array<T, Descriptor<T>::numPop> &f, Array<T, Descriptor<T>::ExternalField::numScalars> &ext,
        Array<T, GPUconst<T, Descriptor>::maxStaticScalars> staticScalars, T *dynamicScalars,
        plint index)
    {
        Collision<T, Descriptor, CollisionModel::GH>::collide(
            f, ext, staticScalars, dynamicScalars, index);
        for (int iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
            if (!std::isnan(dynamicScalars[index + iPop])) {
                f[iPop] += dynamicScalars[index + iPop];
            }
        }
    }
};

template <typename T, template <typename U> class Descriptor>
struct Collision<T, Descriptor, CollisionModel::HalfwayBounceBack__RR> {
    static void collide(
        Array<T, Descriptor<T>::numPop> &f, Array<T, Descriptor<T>::ExternalField::numScalars> &ext,
        Array<T, GPUconst<T, Descriptor>::maxStaticScalars> staticScalars, T *dynamicScalars,
        plint index)
    {
        Collision<T, Descriptor, CollisionModel::RR>::collide(
            f, ext, staticScalars, dynamicScalars, index);
        for (int iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
            if (!std::isnan(dynamicScalars[index + iPop])) {
                f[iPop] += dynamicScalars[index + iPop];
            }
        }
    }
};

}  // namespace plb

#endif  // ACCELERATED_COLLISION_MODELS_H
