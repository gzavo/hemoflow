"""HemoFlow preprocessor.

Turns a vessel surface (STL), its centerline (VTP) and optionally a flow diverter
(STL) into the voxelized geometry file (.npz) read by the hemoFlow solver.

Usage: python main.py case_vox.config
"""
import json
import os
import sys
import time

import nrrd
import numpy as np
import scipy.interpolate as scpinter
from scipy.optimize import linear_sum_assignment
from vtk.numpy_interface import dataset_adapter as dsa
from vtkmodules.vtkIOXML import vtkXMLPolyDataReader

import voxelizeStl as vox
from constants import const
from createFluidSolid import add_walls, crop_like, fluid_components, leak_check
from detectOpenings import paint_openings
from readCL import getOpeningsFromCenterline

#############################
# Parameters to check before execution:
SI_FACTOR = 0.001   # Ratio of the STL/centerline units to [m]. Most STL is in [mm]
DEBUG_MODE = False  # Writes intermediate nrrd output to check with e.g. 3D Slicer
INHOMOGEN = False   # Per-point porous coefficients from <stent_mesh_base>values.vtp
#############################

CONSTANTS = const()


def fail(message):
    print("!!! ERROR:", message)
    sys.exit(-1)


def read_config(confFile):
    """Paths and parameters from the JSON config (paths are relative to the config file)."""
    workDir = os.path.dirname(os.path.abspath(confFile))
    with open(confFile) as f:
        conf = json.load(f)

    stentFile = ""
    if conf.get("stent_folder"):
        dirName = os.path.split(workDir)[1]
        stentFile = os.path.join(workDir, conf["stent_folder"], conf["stent_folder"] + "_" + dirName + "_stent_mesh.stl")
    elif conf.get("stent_mesh_base"):
        stentFile = os.path.join(workDir, conf["stent_mesh_base"] + "mesh.stl")
    # A configured but missing flow diverter file is an error, not a silent "no device".
    if stentFile and not os.path.isfile(stentFile):
        fail(f"stent geometry file specified in config but not found: {stentFile}")

    targetDx = float(conf["target_dx"]) if conf.get("target_dx") not in (None, "") else None
    return {
        "vessel": os.path.join(workDir, conf["geometry_original_stl"]),
        "centerline": os.path.join(workDir, conf["centerline_vtp"]),
        "stent": stentFile,
        "output": os.path.join(workDir, conf["output_base_name"]),
        "targetElements": int(conf["target_elements"]),
        "targetDx": targetDx,
    }


def opening_labels(count):
    """Labels in centerline order (inlet first, outlets by decreasing radius):
    10 for the inlet, 11 for the smallest outlet, 12, 13, ... for the others."""
    labels = [CONSTANTS.INLET_VOXEL]
    middle = count - 2
    labels += [CONSTANTS.OUTLET_REST_VOXEL + i for i in range(middle)]
    if count > 1:
        labels.append(CONSTANTS.OUTLET_VOXEL)
    return labels


def match_caps(caps, clOpenings, dx):
    """Assign a cap to every centerline end point (optimal assignment by distance)."""
    if len(caps) != len(clOpenings):
        fail(f"the surface has {len(caps)} flat caps on its bounding-box faces, but the centerline has "
             f"{len(clOpenings)} end points. Every opening must be a flat cut on a bounding-box face.")
    cost = np.array([[np.linalg.norm(cap.centroid - pos) for cap in caps] for (_, pos, _) in clOpenings])
    rows, cols = linear_sum_assignment(cost)
    order = [None] * len(clOpenings)
    for r, c in zip(rows, cols):
        radius = clOpenings[r][0]
        capRadius = np.sqrt(caps[c].area / np.pi)
        tolerance = 1.5 * max(radius, capRadius) + 2 * dx
        if cost[r, c] > tolerance:
            fail(f"centerline end point {np.round(clOpenings[r][1], 3).tolist()} is {cost[r, c]:.3f} mm from the "
                 f"nearest cap (on face {caps[c].face}); tolerance {tolerance:.3f} mm. Does the centerline end at every opening?")
        order[r] = c
    return order


def interpolate_coefficients(stentBase, grid, shell):
    """Per-voxel porous coefficients from <stentBase>values.vtp (INHOMOGEN mode)."""
    valuesFile = stentBase + "values.vtp"
    if not os.path.isfile(valuesFile):
        print("INHOMOGEN resistance values not found:", valuesFile)
        return shell.astype(np.double), shell.astype(np.double)
    reader = vtkXMLPolyDataReader()
    reader.SetFileName(valuesFile)
    reader.Update()
    data = dsa.WrapDataObject(reader.GetOutput())
    available = list(data.PointData.keys())
    for name in ("linearCoeff", "quadraticCoeff"):
        if name not in available:
            fail(f"{valuesFile} has no point array '{name}' (available: {', '.join(available)})")
    points = (np.asarray(data.GetPoints()) - grid.corner) / grid.dx - 0.5   # voxel index coordinates
    xg, yg, zg = np.mgrid[0:grid.shape[0], 0:grid.shape[1], 0:grid.shape[2]]
    linear = scpinter.griddata(points, np.asarray(data.PointData["linearCoeff"]), (xg, yg, zg), method="nearest")
    quadratic = scpinter.griddata(points, np.asarray(data.PointData["quadraticCoeff"]), (xg, yg, zg), method="nearest")
    print("Linear and quadratic coefficients interpolated from", valuesFile)
    return linear * shell, quadratic * shell


def main(confFile):
    config = read_config(confFile)
    out = config["output"]
    startTime = time.time()

    print("\n### Reading the vessel surface ###")
    surface = vox.read_surface(config["vessel"])
    box = vox.bounds(surface)
    print("Surface:", surface.GetNumberOfCells(), "triangles, bounding box [mm]:", np.round(box, 4).tolist())

    # Open ends (holes) are closed with flat caps; they must lie on bounding-box faces.
    loops = vox.boundary_loops(surface)
    extent = (box[:, 1] - box[:, 0]).max()
    for loop in loops:
        onFace = any(np.abs(loop[:, a] - box[a, s]).max() <= 1e-6 * extent for a in range(3) for s in (0, 1))
        if not onFace:
            print("WARNING: closing a hole that is not on a bounding-box face (mesh defect?) at",
                  np.round(loop.mean(axis=0), 3).tolist())
    if loops:
        print(f"Closing {len(loops)} open end(s) of the surface")
        surface = vox.close_holes(surface)

    print("\n### Voxelizing the vessel ###")
    grid = vox.make_grid(surface, config["targetElements"], config["targetDx"])
    dx = grid.dx * SI_FACTOR
    print(f"Grid: {grid.shape[0]} x {grid.shape[1]} x {grid.shape[2]} voxels, dx = {dx} m")
    fluid = vox.voxelize_closed(surface, grid)
    if DEBUG_MODE:
        nrrd.write(out + "fluid_only.nrrd", fluid.astype(np.uint8))

    caps, rejected = vox.find_caps(surface, tolerance=max(1e-6 * extent, 1e-3 * grid.dx))
    print(f"Found {len(caps)} cap(s) on the bounding-box faces")
    for patch in rejected:
        print(f"Note: ignoring a flat patch of the vessel wall flush with face {patch.face} "
              f"({patch.area:.3g} mm^2); it is not an opening")

    print("\n### Creating walls ###")
    labels, offset = add_walls(fluid)
    originStl = grid.corner + (offset + 0.5) * grid.dx   # Centre of voxel (0, 0, 0), STL units
    origin = originStl * SI_FACTOR
    print("Domain size:", labels.shape, " centre of voxel (0, 0, 0) [STL units]:", np.round(originStl, 6).tolist())
    if DEBUG_MODE:
        nrrd.write(out + "wall_fluid.nrrd", labels)

    print("\n### Matching openings to the centerline ###")
    clOpenings = getOpeningsFromCenterline(config["centerline"])
    capOrder = match_caps(caps, clOpenings, grid.dx)
    openingIndex = opening_labels(len(clOpenings))
    orderedCaps = [caps[c] for c in capOrder]
    painted, info = paint_openings(labels, offset, grid, orderedCaps, openingIndex)

    rOutlets3 = sum(r ** 3 for (r, _, _) in clOpenings[1:])
    print(f"{'label':>5} {'face':>4} {'voxels':>7} {'cap area [mm2]':>15} {'cap r [mm]':>10} {'centerline r [mm]':>17}")
    for (r, _, _), cap, label, inf in zip(clOpenings, orderedCaps, openingIndex, info):
        print(f"{label:>5} {cap.face:>4} {inf['voxels']:>7} {cap.area:>15.3f} {np.sqrt(cap.area / np.pi):>10.3f} {r:>17.3f}")
        if inf["pieces"] > 1:
            print(f"WARNING: opening {label} consists of {inf['pieces']} separate pieces of voxels")

    haveStent = bool(config["stent"])
    if haveStent:
        print("\n### Voxelizing the flow diverter ###")
        stentSurface = vox.read_surface(config["stent"])
        shell = vox.voxelize_shell(stentSurface, grid)
        if INHOMOGEN:
            linear, quadratic = interpolate_coefficients(config["stent"].replace("mesh.stl", ""), grid, shell)
        else:
            linear, quadratic = shell.astype(np.double), shell.astype(np.double)
        shellL = crop_like(shell, offset, painted.shape)
        linearL = crop_like(linear, offset, painted.shape)
        quadraticL = crop_like(quadratic, offset, painted.shape)

        # Porous voxels only inside the lumen; walls, outside and openings stay as they are.
        inLumen = shellL & (painted == CONSTANTS.FLUID_VOXEL)
        print("Flow diverter voxels:", int(shell.sum()), "of which in the lumen:", int(inLumen.sum()),
              "(discarded: wall", int((shellL & (painted == CONSTANTS.WALL_VOXEL)).sum()),
              ", outside", int((shell.sum() - shellL.sum()) + (shellL & (painted == CONSTANTS.UNUSED_VOXEL)).sum()),
              ", openings", int((shellL & (painted >= CONSTANTS.INLET_VOXEL)).sum()), ")")
        linearL[~inLumen] = 0.0
        quadraticL[~inLumen] = 0.0
        print("Min and max of linear coefficients:", np.nanmin(linearL), np.nanmax(linearL))
        print("Min and max of quadratic coefficients:", np.nanmin(quadraticL), np.nanmax(quadraticL))
        if DEBUG_MODE:
            nrrd.write(out + "stent_final.nrrd", inLumen.astype(np.short))
            nrrd.write(out + "stent_linear.nrrd", linearL.astype(np.single))
            nrrd.write(out + "stent_quadratic.nrrd", quadraticL.astype(np.single))

        if max(openingIndex) >= CONSTANTS.FLOW_DIVERTER_VOXEL:
            fail(f"opening labels reach the flow diverter label {CONSTANTS.FLOW_DIVERTER_VOXEL} (too many openings)")
        print("Marking flow diverter voxels in geometry flag with label", CONSTANTS.FLOW_DIVERTER_VOXEL)
        painted[inLumen] = CONSTANTS.FLOW_DIVERTER_VOXEL

    print("\n### Checks ###")
    volume = painted.size
    fluidCount = int(np.count_nonzero(painted >= CONSTANTS.FLUID_VOXEL))
    print(f"Voxels: {volume}, fluid (incl. openings): {fluidCount} ({fluidCount / volume:.1%}), "
          f"walls: {int(np.count_nonzero(painted == CONSTANTS.WALL_VOXEL))}")
    leaks = leak_check(painted)
    if leaks:
        fail(f"{leaks} fluid voxel(s) touch the outside of the vessel or the domain border")
    print("Leak check passed: fluid only meets the outside through the openings")
    components = fluid_components(painted)
    if len(components) > 1:
        print(f"WARNING: the fluid has {len(components)} disconnected parts (sizes {components[:5]}...); "
              "isolated parts carry no flow")
    else:
        print("The fluid is one connected region")

    print("\n### Saving ###")
    # Label volumes in the coordinates of the input STL (points at the voxel centres), so
    # that they overlay the STL and the solver output in ParaView (.vti) or 3D Slicer (.nrrd).
    vox.write_label_image(out + "geometry.vti", painted.astype(np.short, copy=False), originStl, grid.dx)
    nrrd.write(out + "geometry.nrrd", painted, header={
        "space": "left-posterior-superior",   # Slicer reads STL files as LPS by default
        "space directions": np.diag([grid.dx] * 3),
        "space origin": originStl,
        "kinds": ["domain"] * 3,
    })
    arrays = dict(
        geometryFlag=painted.astype(np.short, copy=False),
        dx=np.array([dx], dtype=np.double),
        origin=origin.astype(np.double),
        lengthUnit=np.array([SI_FACTOR], dtype=np.double),
        openingIndex=np.array(openingIndex, dtype=np.short),
        openingRadius=np.array([r * SI_FACTOR for (r, _, _) in clOpenings], dtype=np.double),
        openingNormalizedQRatio=np.array([r ** 3 / rOutlets3 for (r, _, _) in clOpenings], dtype=np.double),
        openingCenter=np.array([inf["center"] for inf in info], dtype=np.double),
        openingNormal=np.array([d for (_, _, d) in clOpenings], dtype=np.double),
    )
    if haveStent:
        arrays.update(stent=inLumen.astype(np.short), linear=linearL, quadratic=quadraticL)
    np.savez_compressed(out + "c.npz", **arrays)
    print("Files:", out + "c.npz", "(solver input),", out + "geometry.vti", "and", out + "geometry.nrrd", "(labels)")
    print(f"\n### DONE :) -> Elapsed total time: {time.time() - startTime:.1f} [s]")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage:", sys.argv[0], "input.config")
        sys.exit(-1)
    main(sys.argv[1])
