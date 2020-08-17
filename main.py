import sys
import time
import numpy as np
import json
import os
from readCL import getOpeningsFromCenterline, convertToVoxelspace
from voxelizeStl import voxelize
from createFluidSolid import createWalls
from detectOpenings import detectOpenings

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
    
    distance = 4  # 4 voxel distance: note, cutting away unused layers might influence this!
    
    for o in radiusTangentVoxelList:
        pos = o[1]
        for j in range(3):
            if inRange(pos[j], 0, distance):
                sidesToCut[j*2+1]=1
            if inRange(pos[j], voxelDomainSize[j], distance):
                sidesToCut[j*2]=1
            
    return np.where(sidesToCut == 1)[0]

if __name__ == "__main__":
    if len(sys.argv) != 2:
        # print("Usage:", sys.argv[0], "input.stl output targetElements \"<optional list of sides to cut>\"  <optional stent_mesh>")
        # print("    e.g.:", sys.argv[0], "input.stl output 6000000 \"0 5\" ")
        print("Usage:", sys.argv[0], "input.config")
        sys.exit(-1) 

    cutWidth = 1 # Might need to set this to 2 if there is more than 1 padding layer for some reason
    distance = 4

    confFile = sys.argv[1]
    workDir = os.path.dirname(confFile)

    with open(confFile) as json_file:
        confData = json.load(json_file)

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

    outputBaseName =  workDir + "/" + confData["output_base_name"]

    targetElem = int(confData["target_elements"])

    voxel_stent_final = np.zeros(0)

    startTime = time.time()

    print("\n### Voxelizing vessel geometry ###")
    voxelVol, domainData = voxelize(vesselGeomFile, targetElem)
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
    radiusTangentVoxelList = convertToVoxelspace(radiusTangentList, domainData[0], domainData[1])
    
    cutList = generateCutList(domainData[2], radiusTangentVoxelList)
    
    print("Computed list of sides to cut away for openings:", cutList)
    
    print("\n### Creating walls and opening in/outlets ###")
    volWithWalls, sliced = createWalls(voxelVol, cutList, cutWidth)

    if DEBUG_MODE:
    	print("\n### Saving nrrd wall geometry ###")    
    	nrrd.write(outputBaseName+"wall_fluid_only.nrrd", volWithWalls)

    print("Size after cutting layers for openings:", volWithWalls.shape)
    volume = np.product(volWithWalls.shape)
    fluids = np.count_nonzero(volWithWalls == 2)
    print("Volume:", volume)
    print("Fluid nodes:", fluids)
    print("Fluid ratio:", fluids / volume)
    print("Walls:", np.count_nonzero(volWithWalls == 1))

    print("\n### Detecting and assigning voxel openings ###")
    openingIdxs, openingCenters, paintedOpenings = detectOpenings(volWithWalls)

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

            # Values from the centerline must be rotated to the same coordinate system as the voxelized geometry
            # Z, X, Y  -> TODO: I think this was for efficiency reason, should be checked at some point
            if inRange3D(cVox, (cCL[2], cCL[0], cCL[1]), distance) is True:
                openingIndex.append(openingIdxs[ccVox])
                openingRadius.append(rCL[0]*SI_FACTOR)
                openingNormalizedQratio.append(rCL[0]**3/r3Tot)  # TODO: it also assigns a number to the inlet, disredards that
                openingCenter.append(cVox)
                openingTangent.append( np.array((rCL[2][2], rCL[2][0], rCL[2][1])) )

    if DEBUG_MODE:
    	print("\n### Saving nrrd geometry flag ###")    
    	nrrd.write(outputBaseName+"geometry.nrrd", paintedOpenings)


    if haveStent:
        print("\n### Voxelizing flow diverter geometry from 3 projections ###")
        print("-> Voxelizing flow diverter geometry projection #1")
        voxelStent, domainData_stent = voxelize(stentGeomFile, targetElem, True, domainData) 
    
        print("-> Voxelizing flow diverter geometry projection #2")
        voxelStent2, domainData_stent = voxelize(stentGeomFile, targetElem, True, domainData, 0) 
    
        print("-> Voxelizing flow diverter geometry projection #3")
        voxelStent3, domainData_stent = voxelize(stentGeomFile, targetElem, True, domainData, 1) 
    
        print("-> Merging projections")
        sdomain_full = np.logical_or(np.logical_or(voxelStent, voxelStent2), voxelStent3)

        sdomain_full = sdomain_full[sliced[0]:sliced[1], sliced[2]:sliced[3], sliced[4]:sliced[5]]

        # Cut outer layer (z,y,x)
        if 0 in cutList:
            sdomain_full = sdomain_full[:,:,cutWidth:]
        if 1 in cutList:
            sdomain_full = sdomain_full[:,:,:-cutWidth]
        if 2 in cutList:
            sdomain_full = sdomain_full[:,cutWidth:,:]
        if 3 in cutList:
            sdomain_full = sdomain_full[:,:-cutWidth,:]
        if 4 in cutList:
            sdomain_full = sdomain_full[cutWidth:,:,:]
        if 5 in cutList:
            sdomain_full = sdomain_full[:-cutWidth,:,:]

        # TODO: why would this be needed, some old bug?
        # voxel_stent_final = sdomain_full[1:-1, 1:-1, 1:-1]
        voxel_stent_final = sdomain_full

        print("Flow diverter domain size after cutting layers for openings:", voxel_stent_final.shape)

        if DEBUG_MODE:
        	print("\n### Saving voxelized flow diverter ###") 
        	nrrd.write(outputBaseName + "stent_final.nrrd", voxel_stent_final.astype(np.short, copy=False))  

    print("\n### Saving final output ###")
    print("File:", outputBaseName+"c.npz")

    # TODO - change into SI + save dx

    #np.savez(sys.argv[2]+".npz", geometryFlag=paintedOpenings, openingDescr=openingDescr, stent=voxel_stent_final.astype(np.short, copy=False))
    np.savez_compressed(outputBaseName+"c.npz", geometryFlag=paintedOpenings, 
                        dx=np.array([DX]).astype(np.double, copy=False),
                        openingIndex=np.array(openingIndex).astype(np.short, copy=False), 
                        openingRadius=np.array(openingRadius).astype(np.double, copy=False), 
                        openingNormalizedQRatio=np.array(openingNormalizedQratio).astype(np.double, copy=False), 
                        openingCenter=np.array(openingCenter).astype(np.double, copy=False), 
                        openingTangent=np.array(openingTangent).astype(np.double, copy=False), 
                        stent=voxel_stent_final.astype(np.short, copy=False))

    endTime = time.time()
    timeElapsed = int(round((endTime - startTime)))
    print("\n### DONE :) -> Elapsed total time:", timeElapsed, "[s]")