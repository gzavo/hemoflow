#ifndef RHEOLOGY_H
#define RHEOLOGY_H

#include <cmath>
#include <limits>

// Carreau-Yasuda model of a shear-thinning fluid:
//
//     nu(g) = nuInf + (nu0 - nuInf) (1 + (lambda g)^a)^((n - 1) / a),
//
// with the shear rate g = sqrt(2 S:S). a = 2 is the Carreau model (e.g. Cho and Kensey,
// Biorheology 28, 241, 1991). Plain C++ (no Palabos), so that the solver can be tested on
// its own. The solver uses lattice units (nu [dx^2/dt], lambda [dt], g [1/dt]).
//
// R is the precision of the computation. CarreauModel (double) holds the parameters; the
// collision converts it to the precision of the lattice, so that single-precision builds do
// not solve in double (on the GPU double runs at 1/64 of the float speed). In double the
// arithmetic is that of the earlier double-only version, bit for bit.
template <typename R>
struct CarreauModelT {
    R nu0 = R(0);      // Zero-shear kinematic viscosity
    R nuInf = R(0);    // Infinite-shear kinematic viscosity (nuInf <= nu0)
    R lambda = R(0);   // Relaxation time
    R n = R(1);        // Power-law index (0 < n <= 1)
    R a = R(2);        // Yasuda exponent (a > 0; 2 = Carreau)

    CarreauModelT() = default;
    // The same parameters in another precision
    template <typename S>
    CarreauModelT(CarreauModelT<S> const& m) : nu0((R)m.nu0), nuInf((R)m.nuInf), lambda((R)m.lambda), n((R)m.n), a((R)m.a)
    { }

    // Relative step at which Newton's method stops: 1e-12 in double, 64 ulp (7.6e-6) in float,
    // where 1e-12 cannot be reached. Newton's method converges quadratically, so the root is
    // then accurate to round-off.
    static constexpr R tolerance()
    {
        return R(64) * std::numeric_limits<R>::epsilon() > R(1e-12) ? R(64) * std::numeric_limits<R>::epsilon() : R(1e-12);
    }

    // Constant viscosity nu0 (no shear dependence)
    bool isNewtonian() const { return lambda == R(0) || n == R(1) || nu0 == nuInf; }

    // (lambda g)^a
    R yasudaTerm(R shearRate) const
    {
        R x = lambda * shearRate;
        return a == R(2) ? x * x : std::pow(x, a);
    }

    R viscosity(R shearRate) const
    {
        return nuInf + (nu0 - nuInf) * std::pow(R(1) + yasudaTerm(shearRate), (n - R(1)) / a);
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
    R shearRateFromK(R K, R invCs2) const
    {
        if (!(K > R(0)))
            return R(0);
        if (isNewtonian())
            return K / (nu0 * invCs2 + R(0.5));
        R lo = K / (nu0 * invCs2 + R(0.5));
        R hi = K / (nuInf * invCs2 + R(0.5));
        R exponent = (n - R(1)) / a;
        R g = hi;
        for (int it = 0; it < 100; ++it) {
            R x = yasudaTerm(g);
            R base = R(1) + x;
            R powTerm = std::pow(base, exponent);
            R nu = nuInf + (nu0 - nuInf) * powTerm;
            R f = g * (nu * invCs2 + R(0.5)) - K;
            if (f == R(0))
                return g;
            if (f < R(0))
                lo = g;
            else
                hi = g;
            R dStress = nuInf + (nu0 - nuInf) * powTerm / base * (R(1) + n * x);
            R next = g - f / (dStress * invCs2 + R(0.5));
            if (!(next >= lo && next <= hi))
                next = R(0.5) * (lo + hi);
            if (std::fabs(next - g) <= tolerance() * g)
                return next;
            g = next;
        }
        return g;
    }

    // Relaxation rate of a cell from the full contraction Pi:Pi of its non-equilibrium
    // momentum flux and its density. At the root g = omega K, so omega = g / K without
    // evaluating the viscosity again.
    R omegaFromPi(R piNeqNormSqr, R rho, R cs2) const
    {
        R invCs2 = R(1) / cs2;
        R K = std::sqrt(R(0.5) * piNeqNormSqr) / (rho * cs2);
        if (!(K > R(0)))
            return R(1) / (nu0 * invCs2 + R(0.5));
        return shearRateFromK(K, invCs2) / K;
    }
};

// The parameters as the solver reads them from the configuration (lattice units)
typedef CarreauModelT<double> CarreauModel;

#endif
