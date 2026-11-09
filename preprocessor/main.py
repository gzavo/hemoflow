import sys
import time
import numpy as np
import json
import os
import scipy.interpolate as scpinter
from scipy.optimize import linear_sum_assignment
import vtk

from readCL import getOpeningsFromCenterline, convertToVoxelspace
from voxelizeStl import voxelize
from createFluidSolid import createWalls
from detectOpenings import detectOpenings, paint_inlets_outlets
from constants import const
from vtk.numpy_interface import dataset_adapter as dsa
from vtkmodules.vtkIOXML import vtkXMLPolyDataReader
from vtk.util import numpy_support

#############################
# Parameters to check before execution:
SI_FACTOR = 0.001 # Ratio to [m]. Most STL is in [mm]

# 'distance' tolerance (opening-detection voxel slack), in units of the opening's own
# local radius. Always used (no longer configurable) instead of one flat voxel count,
# since a flat count doesn't scale with vessel size/resolution: a value tuned for one
# opening can be wrong by an order of magnitude for a smaller branch in the same geometry.
# Keep this well under 1.0: at 1.0x, a large vessel's own tolerance can exceed its
# distance to an unrelated domain face it just happens to sit near coordinate-wise,
# spuriously flagging a face that isn't a real opening at all.
DISTANCE_RADIUS_FACTOR = 0.9

# Position of the voxel centre relative to its index, in voxel units, per axis. The scan-line
# voxelizer samples x and z at integer coordinates; along y (the scan direction) voxel i
# covers [i, i+1), so its centre is at i + 0.5.
VOXEL_SAMPLE_OFFSET = (0.0, 0.5, 0.0)

DEBUG_MODE = False # This will enable additional intermediate nrrd output to check with e.g. 3DSlicer
INHOMOGEN = False
#############################
import nrrd

CONSTANTS = const()


def inRange(value, rangeValue, distance):
    if np.abs(rangeValue-value) < distance:
        return True
    return False

def inRange3D(value3D, rangeValue3D, distance):
    isInRange = True
    for i in range(3):
        isInRange = (isInRange and inRange(value3D[i], rangeValue3D[i], distance) )
    
    return isInRange

def generateCutList(voxelDomainSize, radiusTangentVoxelList, perOpeningDistance):
    sidesToCut = np.zeros(6)

    # perOpeningDistance: one tolerance per centerline point (voxels). If a point is
    # within its own tolerance of a boundary, that side is considered an opening.
    # Note: cutting away unused layers might influence this!

    if DEBUG_MODE:
            print("-> (DEBUG) generatin cutlist -> voxelDomainSize:", voxelDomainSize)

    for o, distance in zip(radiusTangentVoxelList, perOpeningDistance):
        pos = o[1]

        if DEBUG_MODE:
            print("-> (DEBUG) generatin cutlist -> centerline point:", pos, "distance tolerance:", distance)

        for j in range(3):
            if inRange(pos[j], 0, distance):
                sidesToCut[j*2]=1
            if inRange(pos[j], voxelDomainSize[j], distance):
                sidesToCut[j*2+1]=1

    return np.where(sidesToCut == 1)[0]

def cropLikeGeometry(array, sliced, cutList, cutWidth):
    """Apply the cropping of createWalls() (removed empty layers and cut opening layers)
    to an array defined on the uncropped voxel grid."""
    array = array[sliced[0]:sliced[1], sliced[2]:sliced[3], sliced[4]:sliced[5]]
    for axis in range(3):
        index = [slice(None)] * 3
        if 2 * axis in cutList:
            index[axis] = slice(cutWidth, None)
            array = array[tuple(index)]
            index[axis] = slice(None)
        if 2 * axis + 1 in cutList:
            index[axis] = slice(None, -cutWidth)
            array = array[tuple(index)]
    return array

def computeGridOffset(sliced, cutList, cutWidth):
    """Offset (in voxels, per axis) from an index of the final geometry to the corresponding
    index of the voxelization grid: the one-voxel padding of the voxelizer, the removed
    empty layers and the cut opening layers."""
    return np.array([sliced[2 * axis] + (cutWidth if 2 * axis in cutList else 0) - 1 for axis in range(3)], dtype=float)

def computeOrigin(domainData, sliced, cutList, cutWidth):
    """Physical position [m] of voxel (0,0,0) of the final geometry.

    Voxel i of the voxelization samples the scaled and shifted STL at coordinate i
    (x_mm = i / scale - shift). The final array is offset from that grid by the one-voxel
    padding of the voxelizer, the removed empty layers and the cut opening layers."""
    scale, shift = domainData[0], domainData[1]
    gridOffset = computeGridOffset(sliced, cutList, cutWidth)
    origin = np.zeros(3)
    for axis in range(3):
        origin[axis] = ((gridOffset[axis] + VOXEL_SAMPLE_OFFSET[axis]) / scale[axis] - shift[axis]) * SI_FACTOR
    return origin

def scaleAndShiftData(points, scale, shift):
    for i in range(len(points)):
        pts = points[i]
        for j in range(3):
            pts[j] = (pts[j] + shift[j]) * scale[j]
        points[i] = pts
    return points

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage:", sys.argv[0], "input.config")
        sys.exit(-1)

    cutWidth = 1  # wall layers to cut open at each boundary. No longer config-driven.

    confFile = sys.argv[1]
    workDir = os.path.dirname(os.path.abspath(confFile))

    with open(confFile) as json_file:
        confData = json.load(json_file)

    # 'distance' (voxel opening-detection tolerance) and 'cutWidth' are no longer read from
    # the config. 'distance' is always resolved per-opening below (once each opening's
    # physical radius is known) as DISTANCE_RADIUS_FACTOR times that opening's own radius in
    # voxels, instead of one flat default shared by every opening regardless of size.
    print("Using radius-relative 'distance' default: {}x each opening's local radius".format(DISTANCE_RADIUS_FACTOR))
    print("Using fixed 'cutWidth':", cutWidth)

    # Cutlist meaning -> cut one layer from the planes:
    # 0,1 => Xmin, Xmax
    # 2,3 => Ymin, Ymax
    # 4,5 => Zmin, Zmax
    # cutList =  [int(x) for x in confData["cut_list"].split()]

    vesselGeomFile = workDir + "/" + confData["geometry_original_stl"]

    haveStent = False
    stentGeomFile = ""
    if "stent_folder" in confData.keys() and len(confData["stent_folder"]) > 0:
        dirName = os.path.split(workDir)[1]
        stentFileName = confData["stent_folder"] + "_" + dirName + "_stent_mesh.stl"

        stentGeomFile = os.path.join(workDir,confData["stent_folder"],stentFileName)
        haveStent = True
    elif (len(confData["stent_mesh_base"]) > 0):
        stentGeomFile = workDir + "/" + confData["stent_mesh_base"] + "mesh.stl"
        haveStent = True

    # A non-empty 'stent_folder'/'stent_mesh_base' means the config explicitly asks for a
    # stent. Silently falling back to "no stent" when that file is missing (e.g. a typo or
    # a path-resolution bug) would produce a result without a flow diverter and no warning
    # that anything was wrong, so treat a configured-but-missing stent file as fatal.
    if haveStent and not os.path.isfile(stentGeomFile):
        print("!!! ERROR: stent geometry file specified in config but not found:", stentGeomFile)
        sys.exit(-1)

    stentGeomBase = stentGeomFile.replace('mesh.stl','')

    centerLineFile = workDir + "/" + confData["centerline_vtp"]

    outputBaseName =  workDir + "/" + confData["output_base_name"]

    targetElem = int(confData["target_elements"])

    target_dx = None
    try:
        target_dx = float(confData["target_dx"])  # milimeter
    except:
        print("No dx data available")

    voxel_stent_final = np.zeros(0)

    startTime = time.time()

    print("\n### Voxelizing vessel geometry ###")
    voxelVol, domainData = voxelize(vesselGeomFile, targetElem, target_dx=target_dx)

    if DEBUG_MODE:
        print("-> (DEBUG) Saving voxelization result")
        nrrd.write(outputBaseName+"fluid_only.nrrd", voxelVol)

    sx,sy,sz = domainData[0]
    tx,ty,tz = domainData[1]

    print("Voxelized domain size:", voxelVol.shape)
    print("Domain:", domainData[2])
    print("Bounding box:", domainData[3])
    print("Transformation matrix from voxelization:")
    print("Translate: ",tz, tx, ty)
    print("Rotate: ", 90, 0, 90)
    print("Scale: ", sz, sx, sy)

    DX = (1.0/sx)*SI_FACTOR
    print("dx [m]:", DX)

    print("\n### Extracting information on openings from centerline ###")
    radiusTangentList = getOpeningsFromCenterline(centerLineFile)
    print("scale", domainData[0])
    print("translate", domainData[1])
    radiusTangentVoxelList = convertToVoxelspace(radiusTangentList, domainData[0], domainData[1])

    voxScale = domainData[0][0]  # isotropic: same scale on all 3 axes
    perOpeningDistance = [DISTANCE_RADIUS_FACTOR * rtv[0] * voxScale for rtv in radiusTangentVoxelList]

    print("Per-opening distance tolerance [voxels]:", perOpeningDistance)

    cutList = generateCutList(domainData[2], radiusTangentVoxelList, perOpeningDistance)

    print("Computed list of sides to cut away for openings:", cutList)

    print("\n### Creating walls and opening in/outlets ###")
    volWithWalls, sliced = createWalls(voxelVol, cutList, cutWidth)

    if DEBUG_MODE:
        print("\n-> (DEBUG) Saving nrrd wall geometry ###")
        nrrd.write(outputBaseName+"wall_fluid.nrrd", volWithWalls)

    gridOffset = computeGridOffset(sliced, cutList, cutWidth)
    origin = computeOrigin(domainData, sliced, cutList, cutWidth)
    print("Origin of the voxel grid [m]:", origin)

    print("Size after cutting layers for openings:", volWithWalls.shape)
    volume = np.prod(volWithWalls.shape)
    fluids = np.count_nonzero(volWithWalls == 2)
    print("Volume:", volume)
    print("Fluid nodes:", fluids)
    print("Fluid ratio:", fluids / volume)
    print("Walls:", np.count_nonzero(volWithWalls == 1))

    print("\n### Detecting and assigning voxel openings ###")
    inlet_outlets, data = detectOpenings(volWithWalls)

    openingCenters = []
    for io in inlet_outlets:
        oC = np.zeros(3)
        for (x,y,z) in io:
            oC += np.array((z,x,y))
        openingCenters.append(oC/len(io) - 1)  # Remove the one-voxel padding added by detectOpenings()

    # TODO: Assign tangents and radii to voxelized openings

    if len(openingCenters) != len(radiusTangentVoxelList):
        print("!!! ERROR: the number of outlets found on the voxelized domain sides differ from the number found along the centerline! :", len(radiusTangentVoxelList), len(openingCenters))
        sys.exit(-1)

    # The combined information about openings in the correct order (Inlet, Pressure outlet, Other velocity outlets)
    openingIndex = []
    openingRadius = []
    openingNormalizedQratio = []
    openingCenter = []
    inlets_outlets_sorted = []
    openingNormal = []

    # Radial ratio of outlets, note: Qinlet = 1, so it is not included
    r3Tot = np.sum([x[0]**3 for x in radiusTangentVoxelList[1:]])

    # Match each centerline point to its closest detected opening via a global optimal
    # (Hungarian) assignment, instead of a distance-threshold scan. A threshold scan can
    # silently drop a point that has no opening within 'distance' (leaving that opening
    # unpainted, with no BC at all) while still passing the earlier length check, since
    # that check only compares counts, not whether every point actually got matched.
    # The Hungarian assignment always produces a full one-to-one pairing between the two
    # equal-length lists, so no opening can be silently skipped.
    costMatrix = np.zeros((len(radiusTangentVoxelList), len(openingCenters)))
    for ccCL in range(len(radiusTangentVoxelList)):
        cCL = radiusTangentVoxelList[ccCL][1]
        for ccVox in range(len(openingCenters)):
            cVox = openingCenters[ccVox]
            # Compare in the voxelization frame (the centerline points are not cropped)
            costMatrix[ccCL, ccVox] = np.linalg.norm(np.array(cVox) + gridOffset - np.array((cCL[0], cCL[1], cCL[2])))

    clMatchIdx, voxMatchIdx = linear_sum_assignment(costMatrix)

    for ccCL, ccVox in zip(clMatchIdx, voxMatchIdx):
        matchDistance = costMatrix[ccCL, ccVox]
        matchTolerance = perOpeningDistance[ccCL]
        if matchDistance > matchTolerance:
            print("!!! ERROR: closest opening for centerline point is farther than the 'distance' tolerance!")
            print("    Centerline point", ccCL, "at", radiusTangentVoxelList[ccCL][1], "radius", radiusTangentVoxelList[ccCL][0])
            print("    Closest opening", ccVox, "at", openingCenters[ccVox], "distance", matchDistance, "> tolerance", matchTolerance)
            sys.exit(-1)

        rCL = radiusTangentVoxelList[ccCL]
        cVox = openingCenters[ccVox]

        openingRadius.append(rCL[0]*SI_FACTOR)
        openingNormalizedQratio.append(rCL[0]**3/r3Tot)  # TODO: it also assigns a number to the inlet, disredards that
        openingCenter.append(cVox)
        openingNormal.append( np.array((rCL[2][0], rCL[2][1], rCL[2][2])) )

        inlets_outlets_sorted.append(inlet_outlets[ccVox])

    openingIndex, openingCenters, paintedOpenings = paint_inlets_outlets(inlets_outlets_sorted, data, findBoundaryByArea=False)

    if len(openingIndex) != len(openingCenters):
        print("-> (DEBUG) Number of matched openings is incorrect:", len(openingIndex), "insead of", len(openingCenters))

    if haveStent:
        print("\n### Voxelizing flow diverter geometry from 3 projections ###")
        print("-> Voxelizing flow diverter geometry projection #1")
        voxelStent, domainData_stent = voxelize(stentGeomFile, targetElem, True, domainData)

        print("-> Voxelizing flow diverter geometry projection #2")
        voxelStent2, domainData_stent = voxelize(stentGeomFile, targetElem, True, domainData, 0)

        print("-> Voxelizing flow diverter geometry projection #3")
        voxelStent3, domainData_stent = voxelize(stentGeomFile, targetElem, True, domainData, 1)

        print("Starting stent interpolation")
        xg, yg, zg = np.mgrid[0:voxelStent3.shape[0], 0:voxelStent3.shape[1], 0:voxelStent3.shape[2]]

        stentVoxelLinearInterpolate=1
        stentVoxelQuadraticInterpolate=1

        if (INHOMOGEN and (os.path.isfile(stentGeomBase + "values.vtp"))):
            print("INHOMOGEN resistance values found")
            reader = vtkXMLPolyDataReader()
            reader.SetFileName(stentGeomBase + "values.vtp")
            reader.Update()
            vtpdata = dsa.WrapDataObject(reader.GetOutput())
            stentPoints = vtpdata.GetPoints()

            stentLinear = vtpdata.PointData['linearCoeff']
            stentQuadratic = vtpdata.PointData['quadraticCoeff']

            stentPoints = scaleAndShiftData(stentPoints, domainData[0], domainData[1])

            stentVoxelLinearInterpolate = scpinter.griddata(stentPoints, stentLinear, (xg, yg, zg), method="nearest")
            stentVoxelQuadraticInterpolate = scpinter.griddata(stentPoints, stentQuadratic, (xg, yg, zg), method="nearest")
            print("Linear and quadratic coefficients interpolated")
        elif (INHOMOGEN):
            print("INHOMOGEN resistance values not found!")

        print("-> Merging projections")
        sdomain_full = np.logical_or(np.logical_or(voxelStent, voxelStent2), voxelStent3)
        linear = stentVoxelLinearInterpolate*sdomain_full
        quadratic = stentVoxelQuadraticInterpolate*sdomain_full

        # Apply the same cropping as createWalls() applied to the vessel geometry
        voxel_stent_final = cropLikeGeometry(sdomain_full, sliced, cutList, cutWidth)
        voxel_linear_final = cropLikeGeometry(linear, sliced, cutList, cutWidth).astype(np.double)
        voxel_quadratic_final = cropLikeGeometry(quadratic, sliced, cutList, cutWidth).astype(np.double)

        if voxel_stent_final.shape != paintedOpenings.shape:
            print("!!! ERROR: flow diverter grid", voxel_stent_final.shape, "does not match the geometry", paintedOpenings.shape)
            sys.exit(-1)

        # Porous voxels are only allowed inside the lumen. The device voxelization also touches
        # wall and outside voxels; those must stay walls / outside.
        stentVoxels = voxel_stent_final.astype(bool)
        inLumen = stentVoxels & (paintedOpenings == CONSTANTS.FLUID_VOXEL)
        print("Flow diverter voxels:", int(stentVoxels.sum()), "of which in the lumen:", int(inLumen.sum()),
              "(discarded: wall", int((stentVoxels & (paintedOpenings == CONSTANTS.WALL_VOXEL)).sum()),
              ", outside", int((stentVoxels & (paintedOpenings == CONSTANTS.UNUSED_VOXEL)).sum()),
              ", openings", int((stentVoxels & (paintedOpenings >= CONSTANTS.INLET_VOXEL)).sum()), ")")
        voxel_stent_final = inLumen
        voxel_linear_final[~inLumen] = 0.0
        voxel_quadratic_final[~inLumen] = 0.0

        print("Min and max of linear coefficients:", np.nanmin(voxel_linear_final), np.nanmax(voxel_linear_final))
        print("Min and max of quadratic coefficients:", np.nanmin(voxel_quadratic_final), np.nanmax(voxel_quadratic_final))

        print("Flow diverter domain size after cutting layers for openings:", voxel_stent_final.shape)

        if DEBUG_MODE:
            print("-> (DEBUG) Saving voxelized flow diverter")
            nrrd.write(outputBaseName + "stent_final.nrrd", voxel_stent_final.astype(np.short, copy=False))
            nrrd.write(outputBaseName + "stent_linear.nrrd", voxel_linear_final.astype(np.single, copy=False))
            nrrd.write(outputBaseName + "stent_quadratic.nrrd", voxel_quadratic_final.astype(np.single, copy=False))

        # Label the flow diverter voxels for visualization. They are fluid voxels in the lumen,
        # so the solver treats them as fluid; the porous forcing comes from the 'stent' array.
        if max(openingIndex) >= CONSTANTS.FLOW_DIVERTER_VOXEL:
            print("!!! ERROR: opening labels reach the flow diverter label", CONSTANTS.FLOW_DIVERTER_VOXEL,
                  "(too many openings):", openingIndex)
            sys.exit(-1)
        print("Marking flow diverter voxels in geometry flag with label", CONSTANTS.FLOW_DIVERTER_VOXEL)
        paintedOpenings[voxel_stent_final] = CONSTANTS.FLOW_DIVERTER_VOXEL

    print("Saving nrrd geometry flag")
    nrrd.write(outputBaseName+"geometry.nrrd", paintedOpenings)

    print("\n### Saving final output ###")
    print("File:", outputBaseName + "c.npz")

    arrays = dict(
        geometryFlag=paintedOpenings.astype(np.short, copy=False),
        dx=np.array([DX]).astype(np.double, copy=False),
        origin=origin.astype(np.double, copy=False),
        openingIndex=np.array(openingIndex).astype(np.short, copy=False),
        openingRadius=np.array(openingRadius).astype(np.double, copy=False),
        openingNormalizedQRatio=np.array(openingNormalizedQratio).astype(np.double, copy=False),
        openingCenter=np.array(openingCenter).astype(np.double, copy=False),
        openingNormal=np.array(openingNormal).astype(np.double, copy=False),
    )
    if haveStent:
        arrays.update(
            stent=voxel_stent_final.astype(np.short, copy=False),
            linear=voxel_linear_final,
            quadratic=voxel_quadratic_final,
        )
    np.savez_compressed(outputBaseName + "c.npz", **arrays)

    endTime = time.time()
    timeElapsed = int(round((endTime - startTime)))
    print("\n### DONE :) -> Elapsed total time:", timeElapsed, "[s]")
