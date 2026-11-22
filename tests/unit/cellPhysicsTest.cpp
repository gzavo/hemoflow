// Unit tests of cellPhysics.h (doc/dev/TODO-GPU.md, Phase 1): the array-based per-cell physics
// must reproduce the Palabos and HemoFlow CPU code bit for bit, in single and double precision.
// Random cells near equilibrium with a random force; every result is compared with memcmp.
// The float Carreau solver is compared with the double one to a tolerance.
// Built as hemoUnitTests; run directly or with ctest.

#include "cellPhysics.h"
#include "nonNewtonian.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <iostream>
#include <random>
#include <string>

namespace {

int failures = 0;
int checks = 0;

template <typename V>
bool sameBits(V const &a, V const &b)
{
    return std::memcmp(&a, &b, sizeof(V)) == 0;
}

void report(std::string const &name, bool ok)
{
    ++checks;
    if (!ok) {
        ++failures;
        std::cout << "FAILED: " << name << std::endl;
    }
}

// Populations and external scalars of src into dst (dst keeps its own dynamics)
template <typename Real>
void copyState(Cell<Real, descriptors::ForcedD3Q19Descriptor> &dst, Cell<Real, descriptors::ForcedD3Q19Descriptor> const &src)
{
    dst.getRawPopulations() = src.getRawPopulations();
    for (int k = 0; k < descriptors::ForcedD3Q19Descriptor<Real>::ExternalField::numScalars; ++k)
        *dst.getExternal(k) = *src.getExternal(k);
}

template <typename Real>
struct RandomCells {
    typedef descriptors::ForcedD3Q19Descriptor<Real> D;
    std::mt19937 gen;
    std::uniform_real_distribution<double> uniform;
    RandomCells() : gen(20261121), uniform(-1.0, 1.0) { }
    double r() { return uniform(gen); }

    // Equilibrium with density 1 +- 5 %, velocity up to 0.1, off-equilibrium noise and a force
    void fill(Cell<Real, descriptors::ForcedD3Q19Descriptor> &cell)
    {
        Array<Real, 3> u((Real)(0.1 * r()), (Real)(0.1 * r()), (Real)(0.1 * r()));
        iniCellAtEquilibrium(cell, (Real)(1.0 + 0.05 * r()), u);
        for (int i = 0; i < D::q; ++i)
            cell[i] += (Real)(1e-3 * r());
        for (int k = 0; k < 3; ++k)
            *cell.getExternal(D::ExternalField::forceBeginsAt + k) = (Real)(1e-4 * r());
    }
};

// Copy of GuoForcedCarreauDynamics::collideWithOmega before it moved to cellPhysics.h (Cell-based
// Palabos templates), to check that the move did not change a bit. Since Phase 4 the Carreau
// solver runs in the precision of the lattice (rheology.h), here as well.
template <typename Real>
void oldCarreauCollide(Cell<Real, descriptors::ForcedD3Q19Descriptor> &cell, CarreauModel const &model, Real omegaInf)
{
    typedef descriptors::ForcedD3Q19Descriptor<Real> D;
    GuoExternalForceCompleteRegularizedBGKdynamics<Real, descriptors::ForcedD3Q19Descriptor> base(omegaInf);
    Real rhoBar = base.computeRhoBar(cell);
    Array<Real, 3> u, j;
    base.computeVelocity(cell, u);
    Real rho = D::fullRho(rhoBar);
    for (int iD = 0; iD < 3; ++iD)
        j[iD] = rho * u[iD];
    Real rhoBarLb;
    Array<Real, 3> jLb;
    Array<Real, 6> piNeqLb;
    momentTemplates<Real, descriptors::ForcedD3Q19Descriptor>::compute_rhoBar_j_PiNeq(cell, rhoBarLb, jLb, piNeqLb);
    Real jSqrLb = VectorTemplate<Real, descriptors::ForcedD3Q19Descriptor>::normSqr(jLb);
    Real omega = CarreauModelT<Real>(model).omegaFromPi(
        SymmetricTensor<Real, descriptors::ForcedD3Q19Descriptor>::tensorNormSqr(piNeqLb), D::fullRho(rhoBarLb), D::cs2);
    base.regularize(cell, rhoBarLb, jLb, jSqrLb, piNeqLb);
    dynamicsTemplates<Real, descriptors::ForcedD3Q19Descriptor>::complete_bgk_ma2_collision(cell, rhoBar, D::invRho(rhoBar), j, omega);
    externalForceTemplates<Real, descriptors::ForcedD3Q19Descriptor>::addGuoForce(cell, u, omega, (Real)1);
}

template <typename Real>
void testPrecision(std::string const &precision)
{
    typedef descriptors::ForcedD3Q19Descriptor<Real> D;
    typedef CellPhysics<Real, descriptors::ForcedD3Q19Descriptor> Physics;
    typedef Cell<Real, descriptors::ForcedD3Q19Descriptor> CellT;
    const int numCells = 2000;
    const Real omega = (Real)1.7;
    RandomCells<Real> random;

    CarreauModel carreau;   // Blood after Cho and Kensey in lattice units of the pipe case (dx 0.3 mm, dt 1e-4 s)
    carreau.nu0 = 0.056 / 1055 * 1e-4 / 9e-8;
    carreau.nuInf = 0.00345 / 1055 * 1e-4 / 9e-8;
    carreau.lambda = 3.313 / 1e-4;
    carreau.n = 0.3568;
    const Real omegaInf = (Real)(1.0 / (3.0 * carreau.nuInf + 0.5));

    GuoExternalForceCompleteRegularizedBGKdynamics<Real, descriptors::ForcedD3Q19Descriptor> newtonian(omega);
    GuoForcedCarreauDynamics<Real, descriptors::ForcedD3Q19Descriptor> carreauDynamics(omegaInf, carreau);
    BlockStatistics statistics;

    bool velocityOk = true, collideOk = true, externalOk = true, carreauOk = true, carreauExternalOk = true, porousOk = true;
    for (int n = 0; n < numCells; ++n) {
        CellT cell(&newtonian);
        cell.specifyStatisticsStatus(false);
        random.fill(cell);

        // Velocity with the half-force correction
        Array<Real, 3> uPalabos, uArray;
        newtonian.computeVelocity(cell, uPalabos);
        Physics::forcedVelocity(cell.getRawPopulations(), cell.getExternal(0), uArray);
        velocityOk &= sameBits(uPalabos, uArray);

        // Newtonian collision against Palabos' GuoExternalForceCompleteRegularizedBGKdynamics
        CellT reference(cell), shared(cell);
        newtonian.collide(reference, statistics);
        Real rhoBar;
        Physics::collide(shared.getRawPopulations(), shared.getExternal(0), omega, typename Physics::ConstantOmega {omega}, rhoBar);
        collideOk &= sameBits(reference.getRawPopulations(), shared.getRawPopulations());

        // collideExternal with given rhoBar and j
        CellT referenceExt(cell), sharedExt(cell);
        Real rhoBarIn = (Real)(0.02 * random.r());
        Array<Real, 3> jIn((Real)(0.05 * random.r()), (Real)(0.05 * random.r()), (Real)(0.05 * random.r()));
        newtonian.collideExternal(referenceExt, rhoBarIn, jIn, (Real)0, statistics);
        Physics::collideExternal(sharedExt.getRawPopulations(), sharedExt.getExternal(0), rhoBarIn, jIn, omega,
                                 typename Physics::ConstantOmega {omega});
        externalOk &= sameBits(referenceExt.getRawPopulations(), sharedExt.getRawPopulations());

        // Carreau: GuoForcedCarreauDynamics (now on cellPhysics.h) against the previous Cell-based code
        CellT oldCarreau(cell), newCarreau(&carreauDynamics);
        copyState(newCarreau, cell);
        newCarreau.specifyStatisticsStatus(false);
        oldCarreauCollide(oldCarreau, carreau, omegaInf);
        carreauDynamics.collide(newCarreau, statistics);
        carreauOk &= sameBits(oldCarreau.getRawPopulations(), newCarreau.getRawPopulations());

        // Carreau collideExternal against the Carreau collision of the shared code with the same moments
        CellT carreauExt(&carreauDynamics), carreauExtShared(cell);
        copyState(carreauExt, cell);
        carreauExt.specifyStatisticsStatus(false);
        carreauDynamics.collideExternal(carreauExt, rhoBarIn, jIn, (Real)0, statistics);
        Physics::collideExternal(carreauExtShared.getRawPopulations(), carreauExtShared.getExternal(0), rhoBarIn, jIn,
                                 omegaInf, typename Physics::CarreauOmega {carreau});
        carreauExternalOk &= sameBits(carreauExt.getRawPopulations(), carreauExtShared.getRawPopulations());

        // Porous force against the formula of the previous PorousForceFunctional
        Real mLin = (Real)(1.0 + random.r()), mQuad = (Real)(1.0 + random.r());
        Real a = (Real)(0.01 * (1.0 + random.r())), b = (Real)(0.1 * (1.0 + random.r()));
        Array<Real, D::ExternalField::numScalars> ext;
        Physics::porousForce(uPalabos, mLin, mQuad, a, b, &ext[0]);
        Real speed = std::sqrt(VectorTemplate<Real, descriptors::ForcedD3Q19Descriptor>::normSqr(uPalabos));
        Real resistance = mLin * a + mQuad * b * speed;
        for (int iD = 0; iD < 3; ++iD) {
            Real expected = -resistance * uPalabos[iD];
            porousOk &= sameBits(expected, ext[D::ExternalField::forceBeginsAt + iD]);
        }
    }
    report(precision + ": forced velocity", velocityOk);
    report(precision + ": Newtonian collision", collideOk);
    report(precision + ": Newtonian collideExternal", externalOk);
    report(precision + ": Carreau collision", carreauOk);
    report(precision + ": Carreau collideExternal", carreauExternalOk);
    report(precision + ": porous force", porousOk);
}

// Opening nodes on the face (dir, ori) against the Palabos boundary classes HemoFlow uses:
// RegularizedVelocityBoundaryDynamics around the bulk dynamics (velocity and Murray openings) and
// DensityDirichletBoundaryDynamics with PlaneFdBoundaryFunctional3D (pressure openings).
template <typename Real, int dir, int ori>
void testFace(std::string const &precision)
{
    typedef descriptors::ForcedD3Q19Descriptor<Real> D;
    typedef CellPhysics<Real, descriptors::ForcedD3Q19Descriptor> Physics;
    typedef Cell<Real, descriptors::ForcedD3Q19Descriptor> CellT;
    typedef GuoExternalForceCompleteRegularizedBGKdynamics<Real, descriptors::ForcedD3Q19Descriptor> Bulk;
    typedef Array<Real, 3> Vec;
    const Real omega = (Real)1.7;
    RandomCells<Real> random;
    std::string face = precision + ", face (" + std::to_string(dir) + ", " + std::to_string(ori) + ")";

    CarreauModel carreau;
    carreau.nu0 = 0.056 / 1055 * 1e-4 / 9e-8;
    carreau.nuInf = 0.00345 / 1055 * 1e-4 / 9e-8;
    carreau.lambda = 3.313 / 1e-4;
    carreau.n = 0.3568;
    const Real omegaInf = (Real)(1.0 / (3.0 * carreau.nuInf + 0.5));

    RegularizedVelocityBoundaryDynamics<Real, descriptors::ForcedD3Q19Descriptor, dir, ori> velocityBc(new Bulk(omega));
    RegularizedVelocityBoundaryDynamics<Real, descriptors::ForcedD3Q19Descriptor, dir, ori> velocityBcCarreau(
        new GuoForcedCarreauDynamics<Real, descriptors::ForcedD3Q19Descriptor>(omegaInf, carreau));
    DensityDirichletBoundaryDynamics<Real, descriptors::ForcedD3Q19Descriptor, dir, ori> pressureBc(new Bulk(omega));
    BlockStatistics statistics;

    bool velocityNodeOk = true, velocityCollideOk = true, carreauNodeOk = true, pressureVelocityOk = true, fdOk = true;
    for (int n = 0; n < 500; ++n) {
        Vec uImposed((Real)(0.05 * random.r()), (Real)(0.05 * random.r()), (Real)(0.05 * random.r()));

        // Velocity node: completion, then completion + bulk collision
        CellT cell(&velocityBc);
        cell.specifyStatisticsStatus(false);
        random.fill(cell);
        velocityBc.defineVelocity(cell, uImposed);
        CellT completed(cell), collided(cell), shared(cell);
        velocityBc.completePopulations(completed);
        velocityBc.collide(collided, statistics);
        Physics::template completeVelocityNode<dir, ori>(shared.getRawPopulations(), shared.getExternal(0), uImposed, omega);
        velocityNodeOk &= sameBits(completed.getRawPopulations(), shared.getRawPopulations());
        Real rhoBar;
        Physics::collide(shared.getRawPopulations(), shared.getExternal(0), omega, typename Physics::ConstantOmega {omega}, rhoBar);
        velocityCollideOk &= sameBits(collided.getRawPopulations(), shared.getRawPopulations());

        // The same around the Carreau dynamics
        CellT cellC(&velocityBcCarreau);
        cellC.specifyStatisticsStatus(false);
        copyState(cellC, cell);
        velocityBcCarreau.defineVelocity(cellC, uImposed);
        CellT collidedC(cellC), sharedC(cellC);
        velocityBcCarreau.collide(collidedC, statistics);
        Physics::template completeVelocityNode<dir, ori>(sharedC.getRawPopulations(), sharedC.getExternal(0), uImposed, omegaInf);
        Physics::collide(sharedC.getRawPopulations(), sharedC.getExternal(0), omegaInf, typename Physics::CarreauOmega {carreau}, rhoBar);
        carreauNodeOk &= sameBits(collidedC.getRawPopulations(), sharedC.getRawPopulations());

        // Pressure node velocity from the known populations
        Real rhoImposed = (Real)(1.0 + 0.02 * random.r());
        CellT cellP(&pressureBc);
        cellP.specifyStatisticsStatus(false);
        copyState(cellP, cell);
        pressureBc.defineDensity(cellP, rhoImposed);
        Vec uPalabos, uShared;
        pressureBc.computeVelocity(cellP, uPalabos);
        Physics::template pressureNodeVelocity<dir, ori>(cellP.getRawPopulations(), cellP.getExternal(0), D::rhoBar(rhoImposed), uShared);
        pressureVelocityOk &= sameBits(uPalabos, uShared);
    }

    // Finite-difference pressure condition on single nodes of a small lattice (one node at a time:
    // Palabos updates nodes in place, so with several nodes the result depends on their order)
    const plint N = 5;
    for (int n = 0; n < 40; ++n) {
        BlockLattice3D<Real, descriptors::ForcedD3Q19Descriptor> lattice(N, N, N, new Bulk(omega));
        for (plint x = 0; x < N; ++x)
            for (plint y = 0; y < N; ++y)
                for (plint z = 0; z < N; ++z)
                    random.fill(lattice.get(x, y, z));
        plint plane = (ori < 0) ? 0 : N - 1;
        Box3D faceBox(0, N - 1, 0, N - 1, 0, N - 1);
        (dir == 0 ? faceBox.x0 : dir == 1 ? faceBox.y0 : faceBox.z0) = plane;
        (dir == 0 ? faceBox.x1 : dir == 1 ? faceBox.y1 : faceBox.z1) = plane;
        defineDynamics(lattice, faceBox, pressureBc.clone());
        for (plint x = faceBox.x0; x <= faceBox.x1; ++x)
            for (plint y = faceBox.y0; y <= faceBox.y1; ++y)
                for (plint z = faceBox.z0; z <= faceBox.z1; ++z)
                    lattice.get(x, y, z).defineDensity((Real)(1.0 + 0.02 * random.r()));

        plint node[3] = {2, 2, 2};
        node[dir] = plane;
        auto velocityAt = [&](int axis, int shift) {
            plint p[3] = {node[0], node[1], node[2]};
            p[axis] += shift;
            Vec u;
            lattice.get(p[0], p[1], p[2]).computeVelocity(u);
            return u;
        };
        CellT &target = lattice.get(node[0], node[1], node[2]);
        Array<Real, D::q> f = target.getRawPopulations();
        Real rho = target.computeDensity();
        Vec u0, uPlus[3], uMinus[3];
        Physics::template pressureNodeVelocity<dir, ori>(f, target.getExternal(0), D::rhoBar(rho), u0);
        Vec uInward = velocityAt(dir, -ori);
        for (int a = 0; a < 3; ++a)
            if (a != dir) {
                uPlus[a] = velocityAt(a, 1);
                uMinus[a] = velocityAt(a, -1);
            }
        Vec grad[3];
        Physics::template pressureNodeGradient<dir, ori>(u0, uInward, uPlus, uMinus, grad);
        Physics::completePressureNode(f, rho, u0, grad, omega);

        PlaneFdBoundaryFunctional3D<Real, descriptors::ForcedD3Q19Descriptor, dir, ori> fd;
        fd.process(Box3D(node[0], node[0], node[1], node[1], node[2], node[2]), lattice);
        fdOk &= sameBits(f, target.getRawPopulations());
    }

    report(face + ": velocity node completion", velocityNodeOk);
    report(face + ": velocity node collision", velocityCollideOk);
    report(face + ": velocity node, Carreau", carreauNodeOk);
    report(face + ": pressure node velocity", pressureVelocityOk);
    report(face + ": finite-difference pressure node", fdOk);
}

// The Carreau-Yasuda solver before it was templated on the precision (double only)
double oldShearRateFromK(CarreauModel const &m, double K, double invCs2)
{
    if (!(K > 0.0))
        return 0.0;
    if (m.isNewtonian())
        return K / (m.nu0 * invCs2 + 0.5);
    double lo = K / (m.nu0 * invCs2 + 0.5);
    double hi = K / (m.nuInf * invCs2 + 0.5);
    double exponent = (m.n - 1.0) / m.a;
    double g = hi;
    for (int it = 0; it < 100; ++it) {
        double x = m.lambda * g;
        x = m.a == 2.0 ? x * x : std::pow(x, m.a);
        double base = 1.0 + x;
        double powTerm = std::pow(base, exponent);
        double nu = m.nuInf + (m.nu0 - m.nuInf) * powTerm;
        double f = g * (nu * invCs2 + 0.5) - K;
        if (f == 0.0)
            return g;
        if (f < 0.0)
            lo = g;
        else
            hi = g;
        double dStress = m.nuInf + (m.nu0 - m.nuInf) * powTerm / base * (1.0 + m.n * x);
        double next = g - f / (dStress * invCs2 + 0.5);
        if (!(next >= lo && next <= hi))
            next = 0.5 * (lo + hi);
        if (std::fabs(next - g) <= 1e-12 * g)
            return next;
        g = next;
    }
    return g;
}

double oldOmegaFromPi(CarreauModel const &m, double piNeqNormSqr, double rho, double cs2)
{
    double invCs2 = 1.0 / cs2;
    double K = std::sqrt(0.5 * piNeqNormSqr) / (rho * cs2);
    if (!(K > 0.0))
        return 1.0 / (m.nu0 * invCs2 + 0.5);
    return oldShearRateFromK(m, K, invCs2) / K;
}

// Carreau-Yasuda relaxation rate: in double bit-identical to the double-only solver; in float
// within 1e-6 of double, from rest to shear rates far beyond the shear-thinning range, for the
// Carreau (a = 2) and a Yasuda (a = 1.25) model
void testCarreauSolver()
{
    CarreauModel models[2];
    models[0].nu0 = 0.056 / 1055 * 1e-4 / 9e-8;
    models[0].nuInf = 0.00345 / 1055 * 1e-4 / 9e-8;
    models[0].lambda = 3.313 / 1e-4;
    models[0].n = 0.3568;
    models[1] = models[0];
    models[1].a = 1.25;
    bool identical = true, close = true;
    double worst = 0.0;
    for (CarreauModel const &model : models) {
        CarreauModelT<float> single(model);
        for (int k = 0; k <= 400; ++k) {
            double g = 1e-10 * std::pow(10.0, k / 40.0);   // Shear rate 1e-10 .. 1 [1/dt]
            double rho = 1.0 + 0.05 * std::sin(0.37 * k);
            double cs2 = 1.0 / 3.0;
            double K = g * (model.viscosity(g) / cs2 + 0.5);
            double piSqr = 2.0 * (K * rho * cs2) * (K * rho * cs2);
            double omega = model.omegaFromPi(piSqr, rho, cs2);
            identical &= sameBits(omega, oldOmegaFromPi(model, piSqr, rho, cs2));
            float omegaSingle = single.omegaFromPi((float)piSqr, (float)rho, (float)cs2);
            double error = std::fabs(omegaSingle - omega) / omega;
            worst = std::max(worst, error);
            close &= error < 1e-6;
        }
    }
    report("double: Carreau solver", identical);
    report("single: Carreau solver within 1e-6 of double", close);
    std::cout << "Carreau solver, single against double: largest relative difference of omega " << worst << std::endl;
}

template <typename Real>
void testFaces(std::string const &precision)
{
    testFace<Real, 0, -1>(precision);
    testFace<Real, 0, 1>(precision);
    testFace<Real, 1, -1>(precision);
    testFace<Real, 1, 1>(precision);
    testFace<Real, 2, -1>(precision);
    testFace<Real, 2, 1>(precision);
}

}  // namespace

int main(int argc, char *argv[])
{
    plbInit(&argc, &argv);
    testPrecision<double>("double");
    testPrecision<float>("single");
    testFaces<double>("double");
    testFaces<float>("single");
    testCarreauSolver();
    std::cout << (failures == 0 ? "PASSED" : "FAILED") << ": " << checks - failures << " of " << checks
              << " checks passed" << std::endl;
    return failures == 0 ? 0 : 1;
}
