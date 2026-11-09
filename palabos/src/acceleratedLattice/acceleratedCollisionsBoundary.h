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
#ifndef ACCELERATED_COLLISIONS_BOUNDARY_H
#define ACCELERATED_COLLISIONS_BOUNDARY_H

#include "core/globalDefs.h"
#include "core/dynamics.h"
#include "core/dynamics.hh"
#include "core/dynamicsIdentifiers.h"
#include "core/dynamicsIdentifiers.hh"
#include "boundaryCondition/regularizedBoundaryDynamics.h"
#include "boundaryCondition/regularizedBoundaryDynamics.hh"
#include "boundaryCondition/boundaryTemplates.h"

/// All Palabos code is contained in this namespace.
namespace plb {

template <typename T, template <typename U> class Descriptor>
Array<T, Descriptor<T>::d> computeBoundaryJ(
    Array<T, Descriptor<T>::numPop> &f, T rhoBar, int direction, int orientation)
{
    typedef Descriptor<T> L;
    T rhoOnWall = T();
    for (plint iPop = 0; iPop < L::q; ++iPop) {
        if (L::c_gpu(iPop, direction) == 0) {
            rhoOnWall += f[iPop];
        }
    }

    T rhoNormal = T();
    for (plint iPop = 0; iPop < L::q; ++iPop) {
        if (L::c_gpu(iPop, direction) == orientation) {
            rhoNormal += f[iPop];
        }
    }

    // All velocity components parallel to the wall are zero by definition.
    Array<T, Descriptor<T>::d> j_;
    for (int iD = 0; iD < Descriptor<T>::d; ++iD) {
        j_[iD] = T();
    }

    j_[direction] = (T)orientation * ((T)2 * rhoNormal + rhoOnWall - rhoBar);
    return j_;
}

template <typename T, template <typename U> class Descriptor>
T computeBoundaryRhoBar(
    Array<T, Descriptor<T>::numPop> &f, Array<T, 3> const &velocity, int direction, int orientation)
{
    typedef Descriptor<T> L;
    T rhoOnWall = T();
    for (plint iPop = 0; iPop < L::q; ++iPop) {
        if (L::c_gpu(iPop, direction) == 0) {
            rhoOnWall += f[iPop];
        }
    }

    T rhoNormal = T();
    for (plint iPop = 0; iPop < L::q; ++iPop) {
        if (L::c_gpu(iPop, direction) == orientation) {
            rhoNormal += f[iPop];
        }
    }

    T velNormal = (T)orientation * (velocity[direction]);
    T rhoBar = ((T)2 * rhoNormal + rhoOnWall - L::SkordosFactor() * velNormal) / ((T)1 + velNormal);
    return rhoBar;
}

template <typename T, class Descriptor, int direction, int orientation>
struct BoundaryTemplate {
    static void computeBoundaryPiNeqMa2(
        Array<T, Descriptor::q> const &f, T rhoBar, T invRho, Array<T, Descriptor::d> const &j,
        T jSqr, Array<T, SymmetricTensorImpl<T, Descriptor::d>::n> &PiNeq)
    {
        typedef Descriptor L;

        // Compute off-equilibrium for known particle populations.
        Array<T, L::q> fNeq;
        for (plint iPop = 0; iPop < L::q; ++iPop) {
            int c_i_alpha = L::c_gpu(iPop, direction);
            if (c_i_alpha == 0 || c_i_alpha == orientation) {
                fNeq[iPop] =
                    f[iPop]
                    - dynamicsTemplatesImpl<T, typename L::BaseDescriptor>::bgk_ma2_equilibrium(
                        iPop, rhoBar, invRho, j, jSqr);
            }
        }

        int iPi = 0;
        for (int iAlpha = 0; iAlpha < L::d; ++iAlpha) {
            for (int iBeta = iAlpha; iBeta < L::d; ++iBeta) {
                PiNeq[iPi] = T();
                for (plint iPop = 0; iPop < L::q; ++iPop) {
                    if (L::c_gpu(iPop, direction) == 0) {
                        PiNeq[iPi] += L::c_gpu(iPop, iAlpha) * L::c_gpu(iPop, iBeta) * fNeq[iPop];
                    } else if (L::c_gpu(iPop, direction) == orientation) {
                        PiNeq[iPi] +=
                            (T)2 * L::c_gpu(iPop, iAlpha) * L::c_gpu(iPop, iBeta) * fNeq[iPop];
                    }
                }
                ++iPi;
            }
        }
    }
};

template <typename T>
struct BoundaryTemplate<T, descriptors::D3Q19Descriptor<T>, 0, 1> {
    enum { direction = 0, orientation = 1 };

    static void computeBoundaryPiNeqMa2(
        Array<T, 19> const &f, T rhoBar, T invRho, Array<T, 3> const &j, T jSqr, Array<T, 6> &PiNeq)
    {
        typedef descriptors::D3Q19Descriptor<T> L;

        // Compute off-equilibrium for known particle populations.
        Array<T, L::q> fNeq;
        for (plint iPop : {2, 3, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18}) {
            fNeq[iPop] =
                f[iPop]
                - dynamicsTemplatesImpl<T, typename L::BaseDescriptor>::bgk_ma2_equilibrium(
                    iPop, rhoBar, invRho, j, jSqr);
        }

        // Compute PiNeq from fNeq, by using "bounce-back of off-equilibrium part" rule.
        PiNeq[0] = 2. * (fNeq[10] + fNeq[13] + fNeq[14] + fNeq[15] + fNeq[16]);
        PiNeq[1] = 2. * (fNeq[13] - fNeq[14]);
        PiNeq[2] = 2. * (fNeq[15] - fNeq[16]);
        PiNeq[3] = fNeq[2] + fNeq[8] + fNeq[9] + fNeq[11] + fNeq[17] + fNeq[18]
                   + 2. * (fNeq[13] + fNeq[14]);
        PiNeq[4] = fNeq[8] - fNeq[9] + fNeq[17] - fNeq[18];
        PiNeq[5] = fNeq[3] + fNeq[8] + fNeq[9] + fNeq[12] + fNeq[17] + fNeq[18]
                   + 2. * (fNeq[15] + fNeq[16]);
    }
};

template <typename T>
struct BoundaryTemplate<T, descriptors::D3Q19Descriptor<T>, 0, -1> {
    enum { direction = 0, orientation = -1 };

    static void computeBoundaryPiNeqMa2(
        Array<T, 19> const &f, T rhoBar, T invRho, Array<T, 3> const &j, T jSqr, Array<T, 6> &PiNeq)
    {
        typedef descriptors::D3Q19Descriptor<T> L;

        // Compute off-equilibrium for known particle populations.
        Array<T, L::q> fNeq;
        for (plint iPop : {1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12, 17, 18}) {
            fNeq[iPop] =
                f[iPop]
                - dynamicsTemplatesImpl<T, typename L::BaseDescriptor>::bgk_ma2_equilibrium(
                    iPop, rhoBar, invRho, j, jSqr);
        }

        // Compute PiNeq from fNeq, by using "bounce-back of off-equilibrium part" rule.
        PiNeq[0] = 2. * (fNeq[1] + fNeq[4] + fNeq[5] + fNeq[6] + fNeq[7]);
        PiNeq[1] = 2. * (fNeq[4] - fNeq[5]);
        PiNeq[2] = 2. * (fNeq[6] - fNeq[7]);
        PiNeq[3] =
            fNeq[2] + fNeq[8] + fNeq[9] + fNeq[11] + fNeq[17] + fNeq[18] + 2. * (fNeq[4] + fNeq[5]);
        PiNeq[4] = fNeq[8] - fNeq[9] + fNeq[17] - fNeq[18];
        PiNeq[5] =
            fNeq[3] + fNeq[8] + fNeq[9] + fNeq[12] + fNeq[17] + fNeq[18] + 2. * (fNeq[6] + fNeq[7]);
    }
};

template <typename T, template <typename U> class Descriptor>
void regularizeCell(
    Array<T, Descriptor<T>::numPop> &f, T rhoBar, T invRho, Array<T, Descriptor<T>::d> const &j,
    T jSqr, Array<T, SymmetricTensor<T, Descriptor>::n> const &PiNeq)
{
    typedef Descriptor<T> L;
    f[0] = dynamicsTemplatesImpl<T, typename L::BaseDescriptor>::bgk_ma2_equilibrium(
               0, rhoBar, invRho, j, jSqr)
           + offEquilibriumTemplates<T, Descriptor>::fromPiToFneq(0, PiNeq);
    for (plint iPop = 1; iPop <= L::q / 2; ++iPop) {
        f[iPop] =
            dynamicsTemplatesImpl<T, typename Descriptor<T>::BaseDescriptor>::bgk_ma2_equilibrium(
                iPop, rhoBar, invRho, j, jSqr);
        f[iPop + L::q / 2] =
            dynamicsTemplatesImpl<T, typename Descriptor<T>::BaseDescriptor>::bgk_ma2_equilibrium(
                iPop + L::q / 2, rhoBar, invRho, j, jSqr);
        T fNeq = offEquilibriumTemplates<T, Descriptor>::fromPiToFneq(iPop, PiNeq);
        f[iPop] += fNeq;
        f[iPop + L::q / 2] += fNeq;
    }
}

template <
    typename T, template <typename U> class Descriptor, int BaseCollision, int direction,
    int orientation>
void collideRegularizedBoundary(
    Array<T, Descriptor<T>::numPop> &f, Array<T, Descriptor<T>::ExternalField::numScalars> &ext,
    Array<T, GPUconst<T, Descriptor>::maxStaticScalars> staticScalars, T *dynamicScalars,
    plint index)
{
    Array<T, 3> velocity(
        dynamicScalars[index], dynamicScalars[index + 1], dynamicScalars[index + 2]);
    T rhoBar = computeBoundaryRhoBar<T, Descriptor>(f, velocity, direction, orientation);
    T rho = Descriptor<T>::fullRho(rhoBar);
    Array<T, Descriptor<T>::d> j(velocity[0] * rho, velocity[1] * rho, velocity[2] * rho);
    T invRho = Descriptor<T>::invRho(rhoBar);
    T jSqr = VectorTemplate<T, Descriptor>::normSqr(j);

    Array<T, SymmetricTensor<T, Descriptor>::n> PiNeq;
    BoundaryTemplate<T, Descriptor<T>, direction, orientation>::computeBoundaryPiNeqMa2(
        f, rhoBar, invRho, j, jSqr, PiNeq);
    regularizeCell<T, Descriptor>(f, rhoBar, invRho, j, jSqr, PiNeq);

    Collision<T, Descriptor, BaseCollision>::collide(f, ext, staticScalars, dynamicScalars, index);
}

template <
    typename T, template <typename U> class Descriptor, int BaseCollision, int direction,
    int orientation>
void collideRegularizedDensityBoundary(
    Array<T, Descriptor<T>::numPop> &f, Array<T, Descriptor<T>::ExternalField::numScalars> &ext,
    Array<T, GPUconst<T, Descriptor>::maxStaticScalars> staticScalars, T *dynamicScalars,
    plint index)
{
    T rho = dynamicScalars[index];
    T rhoBar = Descriptor<T>::rhoBar(rho);
    auto j = computeBoundaryJ<T, Descriptor>(f, rhoBar, direction, orientation);
    T invRho = Descriptor<T>::invRho(rhoBar);
    T jSqr = VectorTemplate<T, Descriptor>::normSqr(j);

    Array<T, SymmetricTensor<T, Descriptor>::n> PiNeq;
    BoundaryTemplate<T, Descriptor<T>, direction, orientation>::computeBoundaryPiNeqMa2(
        f, rhoBar, invRho, j, jSqr, PiNeq);
    regularizeCell<T, Descriptor>(f, rhoBar, invRho, j, jSqr, PiNeq);

    Collision<T, Descriptor, BaseCollision>::collide(f, ext, staticScalars, dynamicScalars, index);
}

template <typename T, template <typename U> class Descriptor>
struct Collision<T, Descriptor, CollisionModel::Boundary_RegularizedVelocity_0_1__TRT> {
    static void collide(
        Array<T, Descriptor<T>::numPop> &f, Array<T, Descriptor<T>::ExternalField::numScalars> &ext,
        Array<T, GPUconst<T, Descriptor>::maxStaticScalars> staticScalars, T *dynamicScalars,
        plint index)
    {
        collideRegularizedBoundary<T, Descriptor, CollisionModel::TRT, 0, 1>(
            f, ext, staticScalars, dynamicScalars, index);
    }
};

template <typename T, template <typename U> class Descriptor>
struct Collision<T, Descriptor, CollisionModel::Boundary_RegularizedVelocity_0_M1__TRT> {
    static void collide(
        Array<T, Descriptor<T>::numPop> &f, Array<T, Descriptor<T>::ExternalField::numScalars> &ext,
        Array<T, GPUconst<T, Descriptor>::maxStaticScalars> staticScalars, T *dynamicScalars,
        plint index)
    {
        collideRegularizedBoundary<T, Descriptor, CollisionModel::TRT, 0, -1>(
            f, ext, staticScalars, dynamicScalars, index);
    }
};

template <typename T, template <typename U> class Descriptor>
struct Collision<T, Descriptor, CollisionModel::Boundary_RegularizedVelocity_1_1__TRT> {
    static void collide(
        Array<T, Descriptor<T>::numPop> &f, Array<T, Descriptor<T>::ExternalField::numScalars> &ext,
        Array<T, GPUconst<T, Descriptor>::maxStaticScalars> staticScalars, T *dynamicScalars,
        plint index)
    {
        collideRegularizedBoundary<T, Descriptor, CollisionModel::TRT, 1, 1>(
            f, ext, staticScalars, dynamicScalars, index);
    }
};

template <typename T, template <typename U> class Descriptor>
struct Collision<T, Descriptor, CollisionModel::Boundary_RegularizedVelocity_1_M1__TRT> {
    static void collide(
        Array<T, Descriptor<T>::numPop> &f, Array<T, Descriptor<T>::ExternalField::numScalars> &ext,
        Array<T, GPUconst<T, Descriptor>::maxStaticScalars> staticScalars, T *dynamicScalars,
        plint index)
    {
        collideRegularizedBoundary<T, Descriptor, CollisionModel::TRT, 1, -1>(
            f, ext, staticScalars, dynamicScalars, index);
    }
};

template <typename T, template <typename U> class Descriptor>
struct Collision<T, Descriptor, CollisionModel::Boundary_RegularizedVelocity_2_1__TRT> {
    static void collide(
        Array<T, Descriptor<T>::numPop> &f, Array<T, Descriptor<T>::ExternalField::numScalars> &ext,
        Array<T, GPUconst<T, Descriptor>::maxStaticScalars> staticScalars, T *dynamicScalars,
        plint index)
    {
        collideRegularizedBoundary<T, Descriptor, CollisionModel::TRT, 2, 1>(
            f, ext, staticScalars, dynamicScalars, index);
    }
};

template <typename T, template <typename U> class Descriptor>
struct Collision<T, Descriptor, CollisionModel::Boundary_RegularizedVelocity_2_M1__TRT> {
    static void collide(
        Array<T, Descriptor<T>::numPop> &f, Array<T, Descriptor<T>::ExternalField::numScalars> &ext,
        Array<T, GPUconst<T, Descriptor>::maxStaticScalars> staticScalars, T *dynamicScalars,
        plint index)
    {
        collideRegularizedBoundary<T, Descriptor, CollisionModel::TRT, 2, -1>(
            f, ext, staticScalars, dynamicScalars, index);
    }
};

template <typename T, template <typename U> class Descriptor>
struct Collision<T, Descriptor, CollisionModel::Boundary_RegularizedVelocity_0_1__BGK> {
    static void collide(
        Array<T, Descriptor<T>::numPop> &f, Array<T, Descriptor<T>::ExternalField::numScalars> &ext,
        Array<T, GPUconst<T, Descriptor>::maxStaticScalars> staticScalars, T *dynamicScalars,
        plint index)
    {
        collideRegularizedBoundary<T, Descriptor, CollisionModel::BGK, 0, 1>(
            f, ext, staticScalars, dynamicScalars, index);
    }
};

template <typename T, template <typename U> class Descriptor>
struct Collision<T, Descriptor, CollisionModel::Boundary_RegularizedVelocity_0_M1__BGK> {
    static void collide(
        Array<T, Descriptor<T>::numPop> &f, Array<T, Descriptor<T>::ExternalField::numScalars> &ext,
        Array<T, GPUconst<T, Descriptor>::maxStaticScalars> staticScalars, T *dynamicScalars,
        plint index)
    {
        collideRegularizedBoundary<T, Descriptor, CollisionModel::BGK, 0, -1>(
            f, ext, staticScalars, dynamicScalars, index);
    }
};

template <typename T, template <typename U> class Descriptor>
struct Collision<T, Descriptor, CollisionModel::Boundary_RegularizedVelocity_1_1__BGK> {
    static void collide(
        Array<T, Descriptor<T>::numPop> &f, Array<T, Descriptor<T>::ExternalField::numScalars> &ext,
        Array<T, GPUconst<T, Descriptor>::maxStaticScalars> staticScalars, T *dynamicScalars,
        plint index)
    {
        collideRegularizedBoundary<T, Descriptor, CollisionModel::BGK, 1, 1>(
            f, ext, staticScalars, dynamicScalars, index);
    }
};

template <typename T, template <typename U> class Descriptor>
struct Collision<T, Descriptor, CollisionModel::Boundary_RegularizedVelocity_1_M1__BGK> {
    static void collide(
        Array<T, Descriptor<T>::numPop> &f, Array<T, Descriptor<T>::ExternalField::numScalars> &ext,
        Array<T, GPUconst<T, Descriptor>::maxStaticScalars> staticScalars, T *dynamicScalars,
        plint index)
    {
        collideRegularizedBoundary<T, Descriptor, CollisionModel::BGK, 1, -1>(
            f, ext, staticScalars, dynamicScalars, index);
    }
};

template <typename T, template <typename U> class Descriptor>
struct Collision<T, Descriptor, CollisionModel::Boundary_RegularizedVelocity_2_1__BGK> {
    static void collide(
        Array<T, Descriptor<T>::numPop> &f, Array<T, Descriptor<T>::ExternalField::numScalars> &ext,
        Array<T, GPUconst<T, Descriptor>::maxStaticScalars> staticScalars, T *dynamicScalars,
        plint index)
    {
        collideRegularizedBoundary<T, Descriptor, CollisionModel::BGK, 2, 1>(
            f, ext, staticScalars, dynamicScalars, index);
    }
};

template <typename T, template <typename U> class Descriptor>
struct Collision<T, Descriptor, CollisionModel::Boundary_RegularizedVelocity_2_M1__BGK> {
    static void collide(
        Array<T, Descriptor<T>::numPop> &f, Array<T, Descriptor<T>::ExternalField::numScalars> &ext,
        Array<T, GPUconst<T, Descriptor>::maxStaticScalars> staticScalars, T *dynamicScalars,
        plint index)
    {
        collideRegularizedBoundary<T, Descriptor, CollisionModel::BGK, 2, -1>(
            f, ext, staticScalars, dynamicScalars, index);
    }
};

template <typename T, template <typename U> class Descriptor>
struct Collision<T, Descriptor, CollisionModel::Boundary_RegularizedDensity_0_1__TRT> {
    static void collide(
        Array<T, Descriptor<T>::numPop> &f, Array<T, Descriptor<T>::ExternalField::numScalars> &ext,
        Array<T, GPUconst<T, Descriptor>::maxStaticScalars> staticScalars, T *dynamicScalars,
        plint index)
    {
        collideRegularizedDensityBoundary<T, Descriptor, CollisionModel::TRT, 0, 1>(
            f, ext, staticScalars, dynamicScalars, index);
    }
};

template <typename T, template <typename U> class Descriptor>
struct Collision<T, Descriptor, CollisionModel::Boundary_RegularizedDensity_0_M1__TRT> {
    static void collide(
        Array<T, Descriptor<T>::numPop> &f, Array<T, Descriptor<T>::ExternalField::numScalars> &ext,
        Array<T, GPUconst<T, Descriptor>::maxStaticScalars> staticScalars, T *dynamicScalars,
        plint index)
    {
        collideRegularizedDensityBoundary<T, Descriptor, CollisionModel::TRT, 0, -1>(
            f, ext, staticScalars, dynamicScalars, index);
    }
};

template <typename T, template <typename U> class Descriptor>
struct Collision<T, Descriptor, CollisionModel::Boundary_RegularizedDensity_1_1__TRT> {
    static void collide(
        Array<T, Descriptor<T>::numPop> &f, Array<T, Descriptor<T>::ExternalField::numScalars> &ext,
        Array<T, GPUconst<T, Descriptor>::maxStaticScalars> staticScalars, T *dynamicScalars,
        plint index)
    {
        collideRegularizedDensityBoundary<T, Descriptor, CollisionModel::TRT, 1, 1>(
            f, ext, staticScalars, dynamicScalars, index);
    }
};

template <typename T, template <typename U> class Descriptor>
struct Collision<T, Descriptor, CollisionModel::Boundary_RegularizedDensity_1_M1__TRT> {
    static void collide(
        Array<T, Descriptor<T>::numPop> &f, Array<T, Descriptor<T>::ExternalField::numScalars> &ext,
        Array<T, GPUconst<T, Descriptor>::maxStaticScalars> staticScalars, T *dynamicScalars,
        plint index)
    {
        collideRegularizedDensityBoundary<T, Descriptor, CollisionModel::TRT, 1, -1>(
            f, ext, staticScalars, dynamicScalars, index);
    }
};

template <typename T, template <typename U> class Descriptor>
struct Collision<T, Descriptor, CollisionModel::Boundary_RegularizedDensity_2_1__TRT> {
    static void collide(
        Array<T, Descriptor<T>::numPop> &f, Array<T, Descriptor<T>::ExternalField::numScalars> &ext,
        Array<T, GPUconst<T, Descriptor>::maxStaticScalars> staticScalars, T *dynamicScalars,
        plint index)
    {
        collideRegularizedDensityBoundary<T, Descriptor, CollisionModel::TRT, 2, 1>(
            f, ext, staticScalars, dynamicScalars, index);
    }
};

template <typename T, template <typename U> class Descriptor>
struct Collision<T, Descriptor, CollisionModel::Boundary_RegularizedDensity_2_M1__TRT> {
    static void collide(
        Array<T, Descriptor<T>::numPop> &f, Array<T, Descriptor<T>::ExternalField::numScalars> &ext,
        Array<T, GPUconst<T, Descriptor>::maxStaticScalars> staticScalars, T *dynamicScalars,
        plint index)
    {
        collideRegularizedDensityBoundary<T, Descriptor, CollisionModel::TRT, 2, -1>(
            f, ext, staticScalars, dynamicScalars, index);
    }
};

}  // namespace plb

#endif  // ACCELERATED_COLLISIONS_BOUNDARY_H
