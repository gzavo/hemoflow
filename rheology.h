#ifndef RHEOLOGY_H
#define RHEOLOGY_H

#include <cmath>

// Carreau-Yasuda model of a shear-thinning fluid:
//
//     nu(g) = nuInf + (nu0 - nuInf) (1 + (lambda g)^a)^((n - 1) / a),
//
// with the shear rate g = sqrt(2 S:S). a = 2 is the Carreau model (e.g. Cho and Kensey,
// Biorheology 28, 241, 1991). Plain C++ (no Palabos), so that the solver can be tested on
// its own. The solver uses lattice units (nu [dx^2/dt], lambda [dt], g [1/dt]).
struct CarreauModel {
    double nu0 = 0.0;      // Zero-shear kinematic viscosity
    double nuInf = 0.0;    // Infinite-shear kinematic viscosity (nuInf <= nu0)
    double lambda = 0.0;   // Relaxation time
    double n = 1.0;        // Power-law index (0 < n <= 1)
    double a = 2.0;        // Yasuda exponent (a > 0; 2 = Carreau)

    // Constant viscosity nu0 (no shear dependence)
    bool isNewtonian() const { return lambda == 0.0 || n == 1.0 || nu0 == nuInf; }

    // (lambda g)^a
    double yasudaTerm(double shearRate) const
    {
        double x = lambda * shearRate;
        return a == 2.0 ? x * x : std::pow(x, a);
    }

    double viscosity(double shearRate) const
    {
        return nuInf + (nu0 - nuInf) * std::pow(1.0 + yasudaTerm(shearRate), (n - 1.0) / a);
    }

    // In the lattice Boltzmann method the non-equilibrium momentum flux is
    // Pi = -2 rho cs^2 S / omega with 1 / omega = nu / cs^2 + 1/2, so the shear rate is
    // g = omega K with K = |Pi| / (sqrt(2) rho cs^2). As omega depends on g, this solves
    //
    //     g (nu(g) / cs^2 + 1/2) = K
    //
    // for g. With x = (lambda g)^a, d/dg [g nu(g)] = nuInf + (nu0 - nuInf) (1 + x)^((n-1)/a - 1)
    // (1 + n x) > 0 for 0 < n <= 1 and any a > 0, so the root is unique and lies between the
    // Newtonian solutions for nu0 and nuInf. Newton's method, with bisection whenever a step
    // leaves the bracket.
    double shearRateFromK(double K, double invCs2) const
    {
        if (!(K > 0.0))
            return 0.0;
        if (isNewtonian())
            return K / (nu0 * invCs2 + 0.5);
        double lo = K / (nu0 * invCs2 + 0.5);
        double hi = K / (nuInf * invCs2 + 0.5);
        double exponent = (n - 1.0) / a;
        double g = hi;
        for (int it = 0; it < 100; ++it) {
            double x = yasudaTerm(g);
            double base = 1.0 + x;
            double powTerm = std::pow(base, exponent);
            double nu = nuInf + (nu0 - nuInf) * powTerm;
            double f = g * (nu * invCs2 + 0.5) - K;
            if (f == 0.0)
                return g;
            if (f < 0.0)
                lo = g;
            else
                hi = g;
            double dStress = nuInf + (nu0 - nuInf) * powTerm / base * (1.0 + n * x);
            double next = g - f / (dStress * invCs2 + 0.5);
            if (!(next >= lo && next <= hi))
                next = 0.5 * (lo + hi);
            if (std::fabs(next - g) <= 1e-12 * g)
                return next;
            g = next;
        }
        return g;
    }

    // Relaxation rate of a cell from the full contraction Pi:Pi of its non-equilibrium
    // momentum flux and its density. At the root g = omega K, so omega = g / K without
    // evaluating the viscosity again.
    double omegaFromPi(double piNeqNormSqr, double rho, double cs2) const
    {
        double invCs2 = 1.0 / cs2;
        double K = std::sqrt(0.5 * piNeqNormSqr) / (rho * cs2);
        if (!(K > 0.0))
            return 1.0 / (nu0 * invCs2 + 0.5);
        return shearRateFromK(K, invCs2) / K;
    }
};

#endif
