import sys
import time
import numpy as np
import json
import os
from readCL import getOpeningsFromCenterline, convertToVoxelspace
from voxelizeStl import voxelize
from createFluidSolid import createWalls
from detectOpenings import detectOpenings
import h5py


#############################
# Parameters to check before execution:
SI_FACTOR = 0.001 # Ratio to [m]. Most STL is in [mm]

DEBUG_MODE = False # This will enable additional intermediate nrrd output to check with e.g. 3DSlicer
#############################

if DEBUG_MODE:
    import nrrd

def inRange(value, rangeValue, distance):
    if np.abs(rangeValue-value) < distance:
        return True
    return False

def inRange3D(value3D, rangeValue3D, distance):
    isInRange = True
    for i in range(3):
        isInRange = (isInRange and inRange(value3D[i], rangeValue3D[i], distance) )
    
    return isInRange

def generateCutList(voxelDomainSize, radiusTangentVoxelList):
    sidesToCut = np.zeros(6)
    
    # If centerline point is within this distance of the boundary it is considered an opening
    distance = 4  # 4 voxel distance: note, cutting away unused layers might influence this!

    if DEBUG_MODE:
            print("-> (DEBUG) generatin cutlist -> voxelDomainSize:", voxelDomainSize) 
    
    for o in radiusTangentVoxelList:
        pos = o[1]

        if DEBUG_MODE:
            print("-> (DEBUG) generatin cutlist -> centerline point:", pos)    

        for j in range(3):
            if inRange(pos[j], 0, distance):
                sidesToCut[j*2]=1
            if inRange(pos[j], voxelDomainSize[j], distance):
                sidesToCut[j*2+1]=1
            
    return np.where(sidesToCut == 1)[0]

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage:", sys.argv[0], "input.config")
        sys.exit(-1) 

    cutWidth = 1 # Might need to set this to 2 if there is more than 1 padding layer for some reason
    distance = 4

    confFile = sys.argv[1]
    saveDir = os.path.dirname(confFile)

    with open(confFile) as json_file:
        confData = json.load(json_file)
        
    workDir = confData["geometry_dir"]

    # Cutlist meaning -> cut one layer from the planes:
    # 0,1 => Xmin, Xmax
    # 2,3 => Ymin, Ymax
    # 4,5 => Zmin, Zmax
    # cutList =  [int(x) for x in confData["cut_list"].split()]

    vesselGeomFile = workDir + "/" + confData["geometry_original_stl"]
    
    haveStent = False
    if(len(confData["stent_mesh_base"]) > 0):
        haveStent = True
        stentGeomFile = workDir + "/" + confData["stent_mesh_base"] + "mesh.stl"
    
    centerLineFile = workDir + "/" + confData["centerline_vtp"]

    outputBaseName =  saveDir + "/" + confData["output_base_name"]

    targetElem = int(confData["target_elements"])

    voxel_stent_final = np.zeros(0)

    startTime = time.time()

    print("\n### Voxelizing vessel geometry ###")
    voxelVol, domainData = voxelize(vesselGeomFile, targetElem)
    # domainData = voxelize(vesselGeomFile, targetElem)

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
    
    cutList = generateCutList(domainData[2], radiusTangentVoxelList)
    
    print("Computed list of sides to cut away for openings:", cutList)
    
    print("\n### Creating walls and opening in/outlets ###")
    start = time.time()
    volWithWalls, sliced = createWalls(voxelVol, cutList, cutWidth)
    end = time.time()
    print(f"Creating wall time: {end - start}")

    if DEBUG_MODE:
        print("\n-> (DEBUG) Saving nrrd wall geometry ###")    
        nrrd.write(outputBaseName+"wall_fluid.nrrd", volWithWalls)

    print("Size after cutting layers for openings:", volWithWalls.shape)
    volume = np.product(volWithWalls.shape)
    fluids = np.count_nonzero(volWithWalls == 2)
    print("Volume:", volume)
    print("Fluid nodes:", fluids)
    print("Fluid ratio:", fluids / volume)
    print("Walls:", np.count_nonzero(volWithWalls == 1))

    print("\n### Detecting and assigning voxel openings ###")
    start = time.time()
    opening_ids, openingCenters, Inlets_outlets_voxles, paintedOpenings = detectOpenings(volWithWalls, outputBaseName+"wall_fluid.nrrd")
    end = time.time()
    print(f"Detecting opening time: {end - start}")

    # TODO: Assign tangents and radii to voxelized openings
    
    if len(openingCenters) != len(radiusTangentVoxelList):
        print("!!! ERROR: the number of outlets found on the voxelized domain sides differ from the number found along the centerline! :", len(radiusTangentVoxelList), len(openingCenters))
        sys.exit(-1)
    
    # The combined information about openings in the correct order (Inlet, Pressure outlet, Other velocity outlets)
    openingIndex = []
    openingRadius = []
    openingNormalizedQratio = []
    openingCenter = []
    openingTangent = []
    
    # Radial ratio of outlets, note: Qinlet = 1, so it is not included
    r3Tot = np.sum([x[0]**3 for x in radiusTangentVoxelList[1:]])
    
    for ccVox in range(len(openingCenters)):
        for ccCL in range(len(radiusTangentVoxelList)):
            cVox = openingCenters[ccVox]
            rCL = radiusTangentVoxelList[ccCL]
            cCL = rCL[1]

            if inRange3D(cVox, (cCL[0], cCL[1], cCL[2]), distance) is True:
                openingIndex.append(opening_ids[ccVox])
                openingRadius.append(rCL[0]*SI_FACTOR)
                openingNormalizedQratio.append(rCL[0]**3/r3Tot)  # TODO: it also assigns a number to the inlet, disredards that
                openingCenter.append(cVox)
                openingTangent.append( np.array((rCL[2][0], rCL[2][1], rCL[2][2])) )

    if DEBUG_MODE:
        print("-> (DEBUG) Saving nrrd geometry flag")    
        nrrd.write(outputBaseName+"geometry.nrrd", paintedOpenings)
        if len(openingIndex) != len(openingCenters):
            print("-> (DEBUG) Number of matched openings is incorrect:", len(openingIndex), "insead of", len(openingCenters))

    print("\n### Saving final output ###")
    print("File:", outputBaseName+"c.npz")

    #np.savez(sys.argv[2]+".npz", geometryFlag=paintedOpenings, openingDescr=openingDescr, stent=voxel_stent_final.astype(np.short, copy=False))
    np.savez_compressed(outputBaseName+"ops.npz",
                        dx=np.array([DX]).astype(np.double, copy=False),
                        openingIndex=np.array(openingIndex).astype(np.int8, copy=False), 
                        openingRadius=np.array(openingRadius).astype(np.double, copy=False), 
                        openingNormalizedQRatio=np.array(openingNormalizedQratio).astype(np.double, copy=False), 
                        openingCenter=np.array(openingCenter).astype(np.double, copy=False), 
                        openingTangent=np.array(openingTangent).astype(np.double, copy=False), 
                        stent=voxel_stent_final.astype(np.short, copy=False),
                        Nx=paintedOpenings.shape[0], Ny=paintedOpenings.shape[1], Nz=paintedOpenings.shape[2])
    
    print(paintedOpenings.max())
    
    # Create an HDF5 file to save the geometry flag seperately
    with h5py.File(outputBaseName+"gf.h5", "w") as f:
        # Create a dataset in the file and store the array
        dset = f.create_dataset("geometryFlag", shape=paintedOpenings.flatten().shape, 
                                data=paintedOpenings.flatten().astype(np.int8), dtype=np.int8,
                                compression="gzip", compression_opts=9)

    endTime = time.time()
    timeElapsed = int(round((endTime - startTime)))
    print("\n### DONE :) -> Elapsed total time:", timeElapsed, "[s]")