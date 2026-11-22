#ifndef CELL_PHYSICS_H
#define CELL_PHYSICS_H

#include "palabos3D.h"
#include "palabos3D.hh"
#include "rheology.h"

#include <cmath>

using namespace plb;

// HemoFlow's per-cell physics on plain arrays: the populations f and the external scalars ext
// of one cell (ext points at the first external scalar; the force starts at forceBeginsAt).
// No Cell, no virtual calls, no exceptions, so that the CPU dynamics and the GPU collision
// functor (doc/dev/TODO-GPU.md) run the same code. Every function does the arithmetic of the
// Palabos code it reproduces in the same order, through the same array templates that the Cell
// versions of Palabos forward to, so CPU results stay bit-identical.
//
// The relaxation rate comes from an omega function, called as omegaOf(piNeq, rhoBar) with the
// non-equilibrium momentum flux and rhoBar of the cell before the collision:
//   ConstantOmega   Newtonian fluid
//   CarreauOmega    Carreau-Yasuda fluid (rheology.h)
template <typename T_, template <typename U> class Descriptor>
struct CellPhysics {
    typedef Descriptor<T_> D;
    typedef typename D::BaseDescriptor BaseD;
    typedef Array<T_, D::q> Populations;
    typedef Array<T_, D::d> Vector;
    typedef Array<T_, SymmetricTensor<T_, Descriptor>::n> Tensor;

    struct ConstantOmega {
        T_ omega;
        T_ operator()(Tensor const &, T_) const { return omega; }
    };

    // Solved in the precision of the lattice: float in single-precision builds (rheology.h)
    struct CarreauOmega {
        CarreauModelT<T_> model;   // Lattice units
        T_ operator()(Tensor const &piNeq, T_ rhoBar) const
        {
            return model.omegaFromPi(SymmetricTensor<T_, Descriptor>::tensorNormSqr(piNeq), D::fullRho(rhoBar), D::cs2);
        }
    };

    // Velocity with the half-force correction of the Guo scheme, u = j / rho + F / 2
    // (ExternalForceDynamics::computeVelocity).
    static void forcedVelocity(Populations const &f, T_ const *ext, Vector &u)
    {
        T_ rhoBar;
        Vector j;
        momentTemplatesImpl<T_, BaseD>::get_rhoBar_j(f, rhoBar, j);
        forcedVelocityExternal(ext, rhoBar, j, u);
    }

    // The same with given rhoBar and j (ExternalForceDynamics::computeVelocityExternal).
    static void forcedVelocityExternal(T_ const *ext, T_ rhoBar, Vector const &j, Vector &u)
    {
        T_ const *force = ext + D::ExternalField::forceBeginsAt;
        T_ invRho = D::invRho(rhoBar);
        for (int iD = 0; iD < D::d; ++iD)
            u[iD] = j[iD] * invRho + force[iD] / (T_)2;
    }

    // Complete regularized BGK with Guo forcing (GuoExternalForceCompleteRegularizedBGKdynamics::collide),
    // relaxation rate from omegaOf. omegaRegularize is the stored omega of the dynamics; the
    // regularization uses it twice, and only their ratio (1) enters. Returns u^2 (for the block
    // statistics) and the rhoBar of the collision.
    template <class OmegaFunction>
    static T_ collide(Populations &f, T_ *ext, T_ omegaRegularize, OmegaFunction const &omegaOf, T_ &rhoBar)
    {
        rhoBar = momentTemplatesImpl<T_, BaseD>::get_rhoBar(f);
        Vector u, j;
        forcedVelocity(f, ext, u);
        T_ rho = D::fullRho(rhoBar);
        for (int iD = 0; iD < D::d; ++iD)
            j[iD] = rho * u[iD];
        return collideWithMoments(f, ext, rhoBar, u, j, omegaRegularize, omegaOf);
    }

    // The same with given rhoBar and j (GuoExternalForceCompleteRegularizedBGKdynamics::collideExternal).
    template <class OmegaFunction>
    static T_ collideExternal(Populations &f, T_ *ext, T_ rhoBar, Vector const &jIn, T_ omegaRegularize,
                              OmegaFunction const &omegaOf)
    {
        Vector u, j;
        forcedVelocityExternal(ext, rhoBar, jIn, u);
        T_ rho = D::fullRho(rhoBar);
        for (int iD = 0; iD < D::d; ++iD)
            j[iD] = rho * u[iD];
        return collideWithMoments(f, ext, rhoBar, u, j, omegaRegularize, omegaOf);
    }

    // Porous resistance force -(mLin a + mQuad b |u|) u, antiparallel to the velocity u
    // (PorousForceFunctional). a = linCoeff and b = quadCoeff in lattice units, mLin and mQuad
    // the multipliers of the voxel.
    static void porousForce(Vector const &u, T_ mLin, T_ mQuad, T_ linCoeff, T_ quadCoeff, T_ *ext)
    {
        T_ *force = ext + D::ExternalField::forceBeginsAt;
        T_ speed = std::sqrt(VectorTemplate<T_, Descriptor>::normSqr(u));
        T_ resistance = mLin * linCoeff + mQuad * quadCoeff * speed;
        for (int iD = 0; iD < D::d; ++iD)
            force[iD] = -resistance * u[iD];
    }

    // ---- Opening nodes on a face of the domain ----
    // (direction, orientation) as in Palabos: the face normal is along `direction`, and
    // `orientation` is the sign of the outward normal (-1 on the min face, +1 on the max face).
    // The known populations are those with c[direction] == 0 or == orientation (they come from
    // inside the domain). The double literals 0.5 are kept from Palabos, so that single-precision
    // results stay bit-identical; they occur only on opening nodes.

    // Velocity node: rebuilds all populations from the imposed velocity uBoundary
    // (RegularizedVelocityBoundaryDynamics::completePopulations around the bulk dynamics:
    // VelocityDirichletBoundaryDynamics::computeRhoBarJ, boundaryTemplates::compute_PiNeq with the
    // complete equilibrium, regularization). The bulk collision follows (collide).
    template <int direction, int orientation>
    static void completeVelocityNode(Populations &f, T_ const *ext, Vector const &uBoundary, T_ omegaRegularize)
    {
        T_ const *force = ext + D::ExternalField::forceBeginsAt;
        T_ rhoOnWall, rhoNormal;
        wallSums<direction, orientation>(f, rhoOnWall, rhoNormal);
        T_ velNormal = (T_)orientation * (uBoundary[direction] - 0.5 * force[direction]);
        T_ rhoBar = ((T_)2 * rhoNormal + rhoOnWall - D::SkordosFactor() * velNormal) / ((T_)1 + velNormal);
        T_ rho = D::fullRho(rhoBar);
        Vector j;
        for (int iD = 0; iD < D::d; ++iD)
            j[iD] = rho * (uBoundary[iD] - 0.5 * force[iD]);
        T_ jSqr = VectorTemplate<T_, Descriptor>::normSqr(j);

        // Pi^neq from the known populations, "bounce-back of the off-equilibrium part"
        T_ invRho = D::invRho(rhoBar);
        Populations fNeq;
        for (int iPop = 0; iPop < D::q; ++iPop) {
            int c = D::c_gpu(iPop, direction);
            if (c == 0 || c == orientation)
                fNeq[iPop] = f[iPop] - dynamicsTemplatesImpl<T_, BaseD>::complete_bgk_ma2_equilibrium(iPop, rhoBar, invRho, j, jSqr);
        }
        Tensor piNeq;
        int iPi = 0;
        for (int iAlpha = 0; iAlpha < D::d; ++iAlpha) {
            for (int iBeta = iAlpha; iBeta < D::d; ++iBeta) {
                piNeq[iPi] = T_();
                for (int iPop = 0; iPop < D::q; ++iPop)
                    if (D::c_gpu(iPop, direction) == 0)
                        piNeq[iPi] += D::c_gpu(iPop, iAlpha) * D::c_gpu(iPop, iBeta) * fNeq[iPop];
                for (int iPop = 0; iPop < D::q; ++iPop)
                    if (D::c_gpu(iPop, direction) == orientation)
                        piNeq[iPi] += (T_)2 * D::c_gpu(iPop, iAlpha) * D::c_gpu(iPop, iBeta) * fNeq[iPop];
                ++iPi;
            }
        }
        regularize(f, rhoBar, j, jSqr, piNeq, omegaRegularize);
    }

    // Pressure node, finite-difference ("interp") condition: velocity of the node from its known
    // populations and the imposed rhoBar (DensityDirichletBoundaryDynamics::computeVelocity).
    template <int direction, int orientation>
    static void pressureNodeVelocity(Populations const &f, T_ const *ext, T_ rhoBar, Vector &u)
    {
        T_ const *force = ext + D::ExternalField::forceBeginsAt;
        for (int iD = 0; iD < D::d; ++iD)
            u[iD] = T_();
        T_ rhoOnWall, rhoNormal;
        wallSums<direction, orientation>(f, rhoOnWall, rhoNormal);
        u[direction] = (T_)orientation * ((T_)2 * rhoNormal + rhoOnWall - rhoBar);
        T_ invRho = D::invRho(rhoBar);
        for (int iD = 0; iD < D::d; ++iD) {
            u[iD] *= invRho;
            u[iD] += 0.5 * force[iD];
        }
    }

    // Velocity gradient of a pressure node (fd::DirectedGradients3D::o1_velocityDerivative):
    // grad[a][b] = d u_b / d x_a. One-sided first order along the normal (u0: the node,
    // uInward: the neighbour inside the domain), central differences along the face (uPlus[a]
    // and uMinus[a]: the neighbours at +1 and -1 along axis a; unused for a == direction).
    template <int direction, int orientation>
    static void pressureNodeGradient(Vector const &u0, Vector const &uInward, Vector const *uPlus, Vector const *uMinus,
                                     Vector *grad)
    {
        for (int a = 0; a < D::d; ++a) {
            for (int iD = 0; iD < D::d; ++iD) {
                if (a == direction)
                    grad[a][iD] = -orientation * (-u0[iD] + uInward[iD]);
                else
                    grad[a][iD] = (uPlus[a][iD] - uMinus[a][iD]) / (T_)2;
            }
        }
    }

    // Pressure node: rebuilds all populations from rho, the node velocity and its gradient
    // (PlaneFdBoundaryFunctional3D::processCell). omega is the stored omega of the bulk dynamics;
    // it converts the strain rate to Pi^neq, as in Palabos (for Carreau the one at muInf).
    // Afterwards the node gets only the bulk collision.
    static void completePressureNode(Populations &f, T_ rho, Vector const &vel, Vector const *grad, T_ omega)
    {
        typedef SymmetricTensorImpl<T_, D::d> S;
        T_ sToPi = -rho / D::invCs2 / omega;
        Tensor pi;
        pi[S::xx] = (T_)2 * grad[0][0] * sToPi;
        pi[S::yy] = (T_)2 * grad[1][1] * sToPi;
        pi[S::zz] = (T_)2 * grad[2][2] * sToPi;
        pi[S::xy] = (grad[0][1] + grad[1][0]) * sToPi;
        pi[S::xz] = (grad[0][2] + grad[2][0]) * sToPi;
        pi[S::yz] = (grad[1][2] + grad[2][1]) * sToPi;
        Vector j;
        for (int iD = 0; iD < D::d; ++iD)
            j[iD] = rho * vel[iD];
        T_ jSqr = VectorTemplate<T_, Descriptor>::normSqr(j);
        regularize(f, D::rhoBar(rho), j, jSqr, pi, omega);
    }

private:
    // Sums of the populations parallel to the face and of those leaving through it (ascending
    // population index, as indexTemplates::subIndex)
    template <int direction, int orientation>
    static void wallSums(Populations const &f, T_ &rhoOnWall, T_ &rhoNormal)
    {
        rhoOnWall = T_();
        for (int iPop = 0; iPop < D::q; ++iPop)
            if (D::c_gpu(iPop, direction) == 0)
                rhoOnWall += f[iPop];
        rhoNormal = T_();
        for (int iPop = 0; iPop < D::q; ++iPop)
            if (D::c_gpu(iPop, direction) == orientation)
                rhoNormal += f[iPop];
    }

    // GuoExternalForceCompleteRegularizedBGKdynamics::regularize
    static void regularize(Populations &f, T_ rhoBar, Vector const &j, T_ jSqr, Tensor const &piNeq, T_ omega)
    {
        dynamicsTemplatesImpl<T_, BaseD>::complete_bgk_ma2_regularize(
            f, rhoBar, D::invRho(rhoBar), j, jSqr, piNeq, omega, omega, D::d + SymmetricTensor<T_, Descriptor>::n);
    }

    template <class OmegaFunction>
    static T_ collideWithMoments(Populations &f, T_ *ext, T_ rhoBar, Vector const &u, Vector const &j,
                                 T_ omegaRegularize, OmegaFunction const &omegaOf)
    {
        T_ rhoBarLb;
        Vector jLb;
        Tensor piNeqLb;
        momentTemplatesImpl<T_, BaseD>::compute_rhoBar_j_PiNeq(f, rhoBarLb, jLb, piNeqLb);
        T_ jSqrLb = VectorTemplate<T_, Descriptor>::normSqr(jLb);
        T_ omega = omegaOf(piNeqLb, rhoBarLb);

        regularize(f, rhoBarLb, jLb, jSqrLb, piNeqLb, omegaRegularize);
        T_ uSqr = dynamicsTemplatesImpl<T_, BaseD>::complete_bgk_ma2_collision(f, rhoBar, D::invRho(rhoBar), j, omega);
        externalForceTemplatesImpl<T_, D>::addGuoForce(f, ext, u, omega, (T_)1);
        return uSqr;
    }
};

#endif
