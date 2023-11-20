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
    We'll see. But we assume all centerline runs from source (inlet) to target.

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
    
    # inlet line ids. Note, only is correct with the current centerline file !!!!!!!!
    # How did I know? Hardcoded...... dumb but useful, need to find a portable way
    # VMTK has some problems with multiple inlets geometry, therefore need to hardcode
    inlet_line_ids = [4, 5, 6, 7] # and unlabelled inlet
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
        
        if ll == 0: # The first inlet!
            v1 = pArray[0]-pArray[1]
            rTanData.append((r1, pArray[0], -v1))  # The inlet should point inwards
        
        v2 = pArray[-1]-pArray[-3] # Note, the last 2 points are sometimes the same for some reason
        # Inlet! point inwards!
        if ll in(inlet_line_ids):
            rTanData.append((r2, pArray[-1], -v2)) 
        # All other tangents are pointing outwards
        else: 
            rTanData.append((r2, pArray[-1], v2))
    
        # print("Line: ", ll)
        # print(r1, pArray[0], v1)
        # print(r2, pArray[-1], v2)
        
    return rTanData


def convertToVoxelspace(r_tan_data, scale, shift):
    
    voxelList = []
    for opening in r_tan_data:
        r, pos, tan = opening
        newpos = (pos + shift) * scale
        voxelList.append((r, newpos, normalize(tan)))
    
    return voxelList
        