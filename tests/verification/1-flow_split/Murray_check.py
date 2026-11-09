"""Check the outlet flow split against Murray's law.

For every output time step, the flow rate through each opening is measured on a
plane parallel to the opening, a few voxels inside the domain (default 8), over the
cross-section of the branch that ends at the opening. Murray's law distributes the inlet flow over the
outlets proportionally to d^3 (Chnafa et al., AJNR 2018); Murray outlets should
match their share exactly once the flow is established, the pressure outlet
approximately.

Usage: python Murray_check.py <case directory> [--npz input/vox_branch_5M_c.npz] [--depth 8]
"""
import argparse
import glob
import os

import numpy as np
import pyvista as pv


def load_velocity(vti_file):
    """Velocity [m/s] as an (Nx, Ny, Nz, 3) array."""
    grid = pv.read(vti_file)
    nx, ny, nz = grid.dimensions
    u = np.asarray(grid.point_data["velocity [m/s]"]).reshape((nz, ny, nx, 3))
    return np.moveaxis(u, (0, 1, 2), (2, 1, 0))


def connected_region(section, seeds):
    """All fluid pixels of a 2D plane that are 4-connected to the seed pixels."""
    region = np.zeros_like(section)
    stack = [tuple(p) for p in seeds if section[tuple(p)]]
    while stack:
        i, j = stack.pop()
        if region[i, j]:
            continue
        region[i, j] = True
        for di, dj in ((1, 0), (-1, 0), (0, 1), (0, -1)):
            a, b = i + di, j + dj
            if 0 <= a < section.shape[0] and 0 <= b < section.shape[1] and section[a, b] and not region[a, b]:
                stack.append((a, b))
    return region


def opening_planes(flags, labels, normals, depth):
    """For each opening: (axis, inward sign, plane index, mask of the vessel cross-section on that plane).

    The plane lies 'depth' voxels inside the domain, parallel to the opening. The opening
    nodes are moved along the opening direction (which can be oblique to the face) onto the
    plane, and the connected fluid region they hit is the cross-section of that branch."""
    planes = {}
    for label, normal in zip(labels, normals):
        nodes = np.argwhere(flags == label)
        axis = [a for a in range(3)
                if nodes[:, a].min() == nodes[:, a].max() and nodes[0, a] in (0, flags.shape[a] - 1)][0]
        sign = 1 if nodes[0, axis] == 0 else -1
        plane = nodes[0, axis] + sign * depth
        others = [a for a in range(3) if a != axis]
        shift = depth * np.array([normal[a] for a in others]) / abs(normal[axis])
        seeds = np.rint(nodes[:, others] + shift).astype(int)
        inside = np.all((seeds >= 0) & (seeds < [flags.shape[a] for a in others]), axis=1)
        index = [slice(None)] * 3
        index[axis] = plane
        section = flags[tuple(index)] >= 2   # fluid voxels on the plane
        planes[label] = (axis, sign, plane, connected_region(section, seeds[inside]))
    return planes


def main(path, npz, depth):
    geom = np.load(os.path.join(path, npz))
    flags, dx = geom["geometryFlag"], float(geom["dx"][0])
    labels, radii = [int(l) for l in geom["openingIndex"]], list(geom["openingRadius"])
    normals = geom["openingNormal"]
    inlet, outlets = labels[0], labels[1:]
    pressure_outlet = outlets[-1]   # label 11, the smallest outlet, is the 0 Pa outlet in murray mode
    d3 = {l: (2 * r) ** 3 for l, r in zip(labels, radii)}
    sum_d3 = sum(d3[l] for l in outlets)
    planes = opening_planes(flags, labels, normals, depth)

    timesteps = sorted(glob.glob(os.path.join(path, "output", "**", "*.vti"), recursive=True))
    if not timesteps:
        raise SystemExit(f"No .vti files found in {os.path.join(path, 'output')}")

    for timestep in timesteps:
        u = load_velocity(timestep)
        flux = {}
        for label, (axis, sign, plane, mask) in planes.items():
            index = [slice(None)] * 3
            index[axis] = plane
            flux[label] = sign * u[tuple(index)][..., axis][mask].sum() * dx * dx   # inflow positive
        q_in = flux[inlet]
        q_out = -sum(flux[l] for l in outlets)
        print(f"{os.path.basename(timestep)}: inlet VFR {q_in:.4e} m^3/s, sum of outlet VFR {q_out:.4e} m^3/s, "
              f"difference {(q_in - q_out) / q_in:.2%}")
        for l in outlets:
            share = q_in * d3[l] / sum_d3
            kind = "pressure outlet" if l == pressure_outlet else "Murray outlet"
            print(f"    {kind} {l}: VFR {-flux[l]:.4e} m^3/s, Murray share {share:.4e} m^3/s, "
                  f"difference {(-flux[l] - share) / share:.2%}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("d", nargs="?", default="./", help="Case directory (containing output/ and input/)")
    parser.add_argument("--npz", default=os.path.join("input", "vox_branch_5M_c.npz"),
                        help="Geometry file, relative to the case directory")
    parser.add_argument("--depth", type=int, default=8, help="Distance of the measurement planes from the openings [voxels]")
    args = parser.parse_args()
    main(os.path.abspath(args.d), args.npz, args.depth)
