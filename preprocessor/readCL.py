# -*- coding: utf-8 -*-
"""
Created on Wed Aug  5 13:50:33 2020

@author: gzavo
"""

# import sys
import numpy as np
import vtk

# vesselGeomFile = "../working_data/Process/NAP287_vessel3.stl"
# targetElem = 10000000
# voxelVol, domainData = voxelize(vesselGeomFile, targetElem)

# clFileName = "../working_data/Process/NAP287_centerline.vtp"

def getIdList(vtkIdList):
    l = []
    for i in range(vtkIdList.GetNumberOfIds()):
        l.append(vtkIdList.GetId(i))

    return l

def normalize(v):
    norm = np.linalg.norm(v)
    if norm == 0: 
       return v
    return v / norm

def getOpeningsFromCenterline(fileName):
    """
    Parameters
    ----------
    fileName : String
        Name of the centerline file (VTP).

    Returns
    -------
    List of (radius, position, direction) tuples, one per opening. All centerlines must
    run from the source (inlet) to a target (outlet). The first entry is the inlet (the
    start point of the first centerline); the outlets follow, sorted by decreasing radius.
    All directions point into the domain.

    """
    
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(fileName)
    reader.Update()

    data = reader.GetOutput()
    points = data.GetPoints()
    pdata = data.GetPointData()
    rdata = pdata.GetArray("MaximumInscribedSphereRadius") 

    # Read the number of lines in the dataset
    nLines = data.GetNumberOfLines()
    
    rTanData = []
    
    for ll in range(nLines):
        line = data.GetCell(ll)
        # numPoints = line.GetNumberOfPoints()
        
        # Get the IDs of the datapoints on the selected line
        idList = getIdList(line.GetPointIds())
        
        # Get radii on the selected line
        rArray = []
        for i in idList:
            rArray.append(rdata.GetTuple(i)[0])

        rArray = np.array(rArray)
        
        # Get points on the selected line
        pArray = []
        for i in idList:
            pArray.append(np.array(points.GetPoint(i)))
        
        
        r1 = rArray[0]
        r2 = rArray[-1]
        
        if ll == 0:
            v1 = pArray[0]-pArray[1]
            rTanData.append((r1, pArray[0], -v1))  # The inlet should point inwards
        
        # Outlets: the tangent points out of the domain, so it is negated to point inwards
        v2 = pArray[-1]-pArray[-3] # Note, the last 2 points are sometimes the same for some reason
        rTanData.append((r2, pArray[-1], -v2))

    # Keep the inlet first, sort the outlets by decreasing radius
    inlet, outlets = rTanData[0], rTanData[1:]
    outlets.sort(key=lambda opening: opening[0], reverse=True)

    return [inlet] + outlets


def convertToVoxelspace(r_tan_data, scale, shift):
    
    voxelList = []
    for opening in r_tan_data:
        r, pos, tan = opening
        newpos = (pos + shift) * scale
        voxelList.append((r, newpos, normalize(tan)))
    
    return voxelList
        