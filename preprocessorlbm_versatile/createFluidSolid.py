# -*- coding: utf-8 -*-
"""
Created on Wed Sep 07 16:32:07 2011

@author: Gábor Závodszky

Input labeling:
    0 - unused
    1 - fluidx

Output labeling:
    0 - unused
    1 - wall
    2 - fluidx

	
Short description:
	Takes a boolean voxel description as input and outputs
    wall and fluidx cells to nrrd.
	
"""

import numpy as np
import multiprocessing
import time

def roll_zeropad(a, shift, axis=None):
    """
    Roll array elements along a given axis.

    Elements off the end of the array are treated as zeros.

    Parameters
    ----------
    a : array_like
        Input array.
    shift : int
        The number of places by which elements are shifted.
    axis : int, optional
        The axis along which elements are shifted.  By default, the array
        is flattened before shifting, after which the original
        shape is restored.

    Returns
    -------
    res : ndarray
        Output array, with the same shape as `a`.

    See Also
    --------
    roll     : Elements that roll off one end come back on the other.
    rollaxis : Roll the specified axis backwards, until it lies in a
               given position.

    """
    
    a = np.asanyarray(a)
    if shift == 0: return a
    if axis is None:
        n = a.size
        reshape = True
    else:
        n = a.shape[axis]
        reshape = False
    if np.abs(shift) > n:
        res = np.zeros_like(a)
    elif shift < 0:
        shift += n
        zeros = np.zeros_like(a.take(np.arange(n-shift), axis))
        res = np.concatenate((a.take(np.arange(n-shift,n), axis), zeros), axis)
    else:
        zeros = np.zeros_like(a.take(np.arange(n-shift,n), axis))
        res = np.concatenate((zeros, a.take(np.arange(n-shift), axis)), axis)
    if reshape:
        return res.reshape(a.shape)
    else:
        return res
    
def removeUnusedOuterLayers(data):

    print("Removing unused outer layers...")

    x,y,z = data.shape
    sliced = np.array([0,x,0,y,0,z])

    print("Original bounds:", sliced)
    
    while(np.max(data[sliced[0],:,:]) == 0):
        sliced[0] = sliced[0]+1
        #data = data[1:,:,:]
    while(np.max(data[sliced[1]-1,:,:]) == 0):
        sliced[1] = sliced[1]-1
        #data = data[:-1,:,:]
    while(np.max(data[:,sliced[2],:]) == 0):
        sliced[2] = sliced[2]+1
        #data = data[:,1:,:]
    while(np.max(data[:,sliced[3]-1,:]) == 0):
        sliced[3] = sliced[3]-1
        #data = data[:,:-1,:]
    while(np.max(data[:,:,sliced[4]]) == 0):
        sliced[4] = sliced[4]+1
        #data = data[:,:,1:]
    while(np.max(data[:,:,sliced[5]-1]) == 0):
        sliced[5] = sliced[5]-1
        #data = data[:,:,:-1]

    data = data[sliced[0]:sliced[1],sliced[2]:sliced[3],sliced[4]:sliced[5]]

    print("Bounds after removing unused layers:", sliced)

    return data, sliced

def createWalls(inputArray, cutList, cutWidxth = 1):
    
    data1 = inputArray
    # Add extra layer of 0 (widxth 1) to allow for fluidx next to the domain boundary.
    # data1 = np.pad(inputArray, 1, 'constant')

    #Create a cpoy
    data2 = np.copy(data1)

    print("Calculating boundary...")
    
    start = time.time()

    # Serial
    #Circular bool shift to paint walls
    for l in range(-1,2):
        a = roll_zeropad(data1,l,0)
        for m in range(-1,2):
            b = roll_zeropad(a,m,1)
            for n in range(-1,2):
                c = roll_zeropad(b,n,2)
                data2 = data2 | c
    
    end = time.time()
    print(f"roll_zeropad time: {end - start}")
                
    #cut out fluidx part, so wall remains
    data3 = data2-data1
    
    #Convert datatype to int to allow third label
    data3 = data3.astype('int8')
    
    #Paint fluidx cells back
    data3[data1==1]=2

    # Remove layers with no information
    start = time.time()
    data3, sliced = removeUnusedOuterLayers(data3)
    end = time.time()
    print(f"removeUnusedOuterLayers time: {end - start}")

    # Cutlist meaning:
    # 0,1 => Xmin, Xmax
    # 2,3 => Ymin, Ymax
    # 4,5 => Zmin, Zmax

    # Cut outer layer (x,y,z)
    if 0 in cutList:
        data3 = data3[cutWidxth:,:,:]
    if 1 in cutList:
        data3 = data3[:-cutWidxth,:,:]
    if 2 in cutList:
        data3 = data3[:,cutWidxth:,:]
    if 3 in cutList:
        data3 = data3[:,:-cutWidxth,:]
    if 4 in cutList:
        data3 = data3[:,:,cutWidxth:]
    if 5 in cutList:
        data3 = data3[:,:,:-cutWidxth]

    return data3.astype('int8', copy=False), sliced


if __name__ == "__main__":
    import sys
    if len(sys.argv) < 3:
        print("Usage:", sys.argv[0], "input.nrrd output.nrrd \"<list of sidxes to cut>\" ")
        print("    e.g.:", sys.argv[0], "input.nrrd output.nrrd \"0 5\" ")
        sys.exit(-1) 
    
    import nrrd

    print("Loading data files...")
    
    #Lodaing input array
    data1, data_header = nrrd.read(sys.argv[1])

    if len(sys.argv) > 3:
        cutList = [int(x) for x in sys.argv[3].split()]
    else:
        cutList = []

    data2 = createWalls(data1, cutList)

    print("Saving...")    
    
    nrrd.write(sys.argv[2], data2)
        
    print("Done.")
