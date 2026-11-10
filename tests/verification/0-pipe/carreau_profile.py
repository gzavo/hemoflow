"""Fully developed pipe flow of a Carreau-Yasuda (or Newtonian) fluid.

Compares the simulated flow in the middle of the pipe with the semi-analytical
solution for the same flow rate: in fully developed flow the shear stress is
tau(r) = G r / 2 (G = -dp/dx), the shear rate follows from mu(g) g = tau, and
u(r) is the integral of the shear rate from the wall inwards.

The grid makes the pipe act slightly narrower than the radius of its lumen area.
With --calibrate newtonian.xml the radius is instead the one for which
Hagen-Poiseuille reproduces the Newtonian twin run, which removes this grid
effect from the comparison of the Carreau run.

Usage: python carreau_profile.py carreau.xml [--calibrate newtonian.xml]
       python carreau_profile.py newtonian.xml
"""

import argparse
import glob
import os
import re

import numpy as np
import pyvista as pv

# Defaults of the solver (globals.h)
DEFAULTS = {"density": 1055.0, "kinematicViscosity": 3.22e-6, "mu0": 0.056, "muInf": 0.00345,
            "lambda": 3.313, "n": 0.3568, "a": 2.0}


def read_config(xml_file):
    text = re.sub(r"<!--.*?-->", "", open(xml_file).read(), flags=re.S)

    def tag(name, section=text, default=None):
        m = re.search(rf"<{name}>\s*([^<]*?)\s*</{name}>", section)
        return m.group(1) if m else default

    m = re.search(r"<fluid>(.*?)</fluid>", text, flags=re.S)
    fluid = m.group(1) if m else ""
    cfg = {"outputDir": tag("outputDir"), "geometry": tag("file"), "model": tag("model", fluid, "newtonian")}
    for key, value in DEFAULTS.items():
        cfg[key] = float(tag(key, fluid, value))
    base = os.path.dirname(os.path.abspath(xml_file))
    cfg["outputDir"] = os.path.join(base, cfg["outputDir"])
    cfg["geometry"] = os.path.join(base, cfg["geometry"])
    return cfg


def viscosity_function(cfg):
    """Dynamic viscosity mu(g) [Pa s] and its bounds."""
    if cfg["model"] == "carreau":
        mu0, mu_inf, lam, n, a = (cfg[k] for k in ("mu0", "muInf", "lambda", "n", "a"))
        return (lambda g: mu_inf + (mu0 - mu_inf) * (1.0 + (lam * g) ** a) ** ((n - 1.0) / a)), mu_inf, mu0
    mu = cfg["kinematicViscosity"] * cfg["density"]
    return (lambda g: np.full_like(g, mu)), mu, mu


def shear_rate(tau, mu, mu_min, mu_max):
    """Solve mu(g) g = tau (stress increases with g) by bisection in log space."""
    tau = np.asarray(tau, dtype=float)
    g = np.zeros_like(tau)
    pos = tau > 0
    lo, hi = tau[pos] / mu_max, tau[pos] / mu_min
    for _ in range(200):
        mid = np.sqrt(lo * hi)
        low = mu(mid) * mid < tau[pos]
        lo, hi = np.where(low, mid, lo), np.where(low, hi, mid)
    g[pos] = np.sqrt(lo * hi)
    return g


def analytical(G, R, mu, mu_min, mu_max, n=4001):
    """Profile u(r) and flow rate Q of fully developed flow for the pressure gradient G."""
    r = np.linspace(0.0, R, n)
    g = shear_rate(G * r / 2.0, mu, mu_min, mu_max)
    seg = 0.5 * (g[1:] + g[:-1]) * np.diff(r)              # integral of g over each segment
    u = np.concatenate([np.cumsum(seg[::-1])[::-1], [0.0]])  # u(r) = int_r^R g ds
    Q = np.trapezoid(2.0 * np.pi * r * u, r)
    return r, u, g, Q


def gradient_for_flow(Q, R, mu, mu_min, mu_max):
    lo, hi = 8.0 * mu_min * Q / (np.pi * R**4), 8.0 * mu_max * Q / (np.pi * R**4)
    for _ in range(100):
        G = np.sqrt(lo * hi)
        if analytical(G, R, mu, mu_min, mu_max, 801)[3] < Q:
            lo = G
        else:
            hi = G
    return np.sqrt(lo * hi)


def load(vti, shape):
    """Point data of the output as (nx, ny, nz[, 3]) arrays."""
    grid = pv.read(vti)
    nx, ny, nz = shape

    def field(name):
        a = np.asarray(grid.point_data[name])
        a = a.reshape((nz, ny, nx) + a.shape[1:])
        return np.moveaxis(a, (0, 2), (2, 0))

    out = {"velocity": field("velocity [m/s]"), "pressure": field("pressure [Pa]")}
    if "viscosity [Pa s]" in grid.point_data:
        out["viscosity"] = field("viscosity [Pa s]")
    if "stress [Pa]" in grid.point_data:
        out["stress"] = field("stress [Pa]")   # xx, xy, xz, yy, yz, zz
    return out


def measure(cfg):
    """Pressure gradient, flow rate and lumen-area radius in the fully developed part."""
    geom = np.load(cfg["geometry"])
    flags, dx = geom["geometryFlag"], float(geom["dx"][0])
    fluid = flags >= 2
    nx = flags.shape[0]
    files = sorted(glob.glob(os.path.join(cfg["outputDir"], "*.vti")))
    last, prev = load(files[-1], flags.shape), load(files[-2], flags.shape)

    # Fully developed part: 45-80 % of the length (entrance length ~0.3 L at muInf)
    xs = np.arange(int(0.45 * nx), int(0.80 * nx))
    p = np.array([last["pressure"][i][fluid[i]].mean() for i in xs])
    G_sim = -np.polyfit(xs * dx, p, 1)[0]
    Q_sim = np.mean([last["velocity"][i][..., 0][fluid[i]].sum() * dx**2 for i in xs])
    R_area = np.sqrt(np.mean([fluid[i].sum() for i in xs]) * dx**2 / np.pi)
    return dict(fluid=fluid, dx=dx, nx=nx, files=files, last=last, prev=prev, G=G_sim, Q=Q_sim, R=R_area)


def main(xml_file, calibrate=None):
    cfg = read_config(xml_file)
    m = measure(cfg)
    fluid, dx, nx, files, last, prev = (m[k] for k in ("fluid", "dx", "nx", "files", "last", "prev"))
    G_sim, Q_sim, R = m["G"], m["Q"], m["R"]
    radius_note = "radius of the lumen area"
    if calibrate:
        ref_cfg = read_config(calibrate)
        assert ref_cfg["model"] == "newtonian", "--calibrate needs a Newtonian run"
        ref = measure(ref_cfg)
        mu_ref = ref_cfg["kinematicViscosity"] * ref_cfg["density"]
        R = (8.0 * mu_ref * ref["Q"] / (np.pi * ref["G"])) ** 0.25
        radius_note = f"hydraulic radius from {os.path.basename(calibrate)}; lumen-area radius {m['R'] * 1e3:.4f} mm"

    mu, mu_min, mu_max = viscosity_function(cfg)
    G_an = gradient_for_flow(Q_sim, R, mu, mu_min, mu_max)
    r_an, u_an, g_an, _ = analytical(G_an, R, mu, mu_min, mu_max)
    G_newt = 8.0 * mu_min * Q_sim / (np.pi * R**4)       # Hagen-Poiseuille with the lowest viscosity

    # Profile on the middle cross-section, against the distance from the lumen centroid
    i = int(0.625 * nx)
    yy, zz = np.nonzero(fluid[i])
    r = np.hypot(yy - yy.mean(), zz - zz.mean()) * dx
    u_sim = last["velocity"][i][..., 0][fluid[i]]
    u_ref = np.interp(r, r_an, u_an)
    peak = u_an[0]
    interior = r < R - 1.5 * dx                          # exclude the wall-adjacent ring
    steady = np.abs(last["velocity"] - prev["velocity"]).max() / peak

    print(f"Model: {cfg['model']}, outputs {os.path.basename(files[-2])} and {os.path.basename(files[-1])}")
    print(f"Radius {R * 1e3:.4f} mm ({radius_note})")
    print(f"Flow rate {Q_sim:.4e} m^3/s, mean velocity over the lumen {Q_sim / (np.pi * m['R']**2):.4f} m/s")
    print(f"Change between the last two outputs: {steady:.1e} of the peak velocity")
    print(f"Pressure gradient: simulated {G_sim:.2f} Pa/m, analytical {G_an:.2f} Pa/m "
          f"({(G_sim / G_an - 1) * 100:+.2f} %); Newtonian at mu = {mu_min} Pa s: {G_newt:.2f} Pa/m")
    print(f"Peak velocity: simulated {u_sim.max():.5f} m/s, analytical {peak:.5f} m/s ({(u_sim.max() / peak - 1) * 100:+.2f} %)")
    err = (u_sim - u_ref)[interior] / peak
    print(f"Profile error (not counting the wall-adjacent ring): RMS {np.sqrt(np.mean(err**2)) * 100:.2f} %, "
          f"max {np.abs(err).max() * 100:.2f} % of the peak velocity")
    if "viscosity" in last and cfg["model"] == "carreau":
        mu_sim = last["viscosity"][i][fluid[i]]
        mu_ref = np.interp(r, r_an, mu(g_an))
        band = (r > 0.2 * R) & interior
        rel = mu_sim[band] / mu_ref[band] - 1
        print(f"Viscosity at 0.2 R < r < R - 1.5 dx: median {np.median(np.abs(rel)) * 100:.2f} %, "
              f"max {np.abs(rel).max() * 100:.2f} % from mu(g_analytical(r)); "
              f"range {mu_sim[band].min():.5f}-{mu_sim[band].max():.5f} Pa s")

    if "stress" in last:
        # Fully developed flow: the shear stress balances the pressure gradient, |sigma_rx| = G r / 2
        sigma = last["stress"][i][fluid[i]]
        shear = np.hypot(sigma[:, 1], sigma[:, 2])
        band = (r > 0.2 * R) & interior
        rel = shear[band] / (G_sim * r[band] / 2) - 1
        print(f"Shear stress at 0.2 R < r < R - 1.5 dx against G r / 2: median {np.median(np.abs(rel)) * 100:.2f} %, "
              f"max {np.abs(rel).max() * 100:.2f} %")

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots()
        ax.plot(r * 1e3, u_sim, ".", ms=3, label="simulation")
        ax.plot(r_an * 1e3, u_an, label=f"analytical ({cfg['model']})")
        if cfg["model"] == "carreau":
            ax.plot(r_an * 1e3, 2 * Q_sim / (np.pi * R**2) * (1 - (r_an / R) ** 2), "--", label="Poiseuille, same flow")
        ax.set(xlabel="r [mm]", ylabel="axial velocity [m/s]", title=os.path.basename(xml_file))
        ax.legend()
        ax.grid()
        name = f"profile_{os.path.splitext(os.path.basename(xml_file))[0]}.png"
        fig.savefig(os.path.join(os.path.dirname(os.path.abspath(xml_file)), name), dpi=120)
    except ImportError:
        pass


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("config", help="Solver configuration (carreau.xml or newtonian.xml)")
    parser.add_argument("--calibrate", metavar="NEWTONIAN_XML", help="Take the effective radius from this Newtonian run")
    args = parser.parse_args()
    main(args.config, args.calibrate)
