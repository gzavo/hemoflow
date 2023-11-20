# import os.path
# import sys
import numpy as np
import time

import slice
import perimeter
import multiprocessing
from multiprocessing import set_start_method
from util import padVoxelArray

from stl import mesh

def voxelize(inputFile, targetElements, isMeshAShell = False, domainData = None, rotation =-1):
    mesh = list(import_stl_file(inputFile))
    
    if domainData is None:
        scale, shift, domain, bounding_box = slice.calculateScaleAndShift(mesh, targetElements)
    else:
        scale, shift, domain, bounding_box = domainData
    
    mesh = list(slice.scaleAndShiftMesh(mesh, scale, shift))

    # vol = np.zeros((domain[2],domain[0],domain[1]), dtype=bool)
    #Note: vol should be addressed with vol[z][x][y]
    
    # If it is a shell mesh, rotate here for different projections!
    if rotation == 0:
        mesh = list(slice.swapAxis(mesh, 0, 1))
        domain = [domain[1], domain[0], domain[2]]
    elif rotation == 1:
        mesh = list(slice.swapAxis(mesh, 1, 2))
        domain = [domain[0], domain[2], domain[1]]
    
    start = time.time()
    
    # Parallel
    multiprocessing.set_start_method('fork')
    pool = multiprocessing.Pool(16)
    import functools
    vol = np.array(pool.map(functools.partial(render_slice, mesh=mesh, domain=domain, isMeshAShell=isMeshAShell), range(int(domain[2]))))
    
    pool.close()
    
    # # Serial
    # vol_list = []
    # for height in range(int(domain[2])):
    #     vol_list.append(render_slice(height, mesh, domain, isMeshAShell))
    # vol = np.stack(np.array(vol_list), dtype=bool, axis=0)
    # del vol_list
    
    print(vol.shape)

    vol, domain = padVoxelArray(vol)
    
    end = time.time()
    print(f"Voxelization time: {end - start}")

    # If it is a shell mesh, rotate back here the different projections!
    if rotation == 0:
        vol = np.swapaxes(vol, 1, 2)
    elif rotation == 1:
        #mesh = slice.swapAxisWithZ(mesh, 1)
        vol = np.swapaxes(vol, 0, 2)


    # At the very end here fix the [z][x][y] back to [x][y][z]
    vol = np.swapaxes(vol, 0, 2)
    vol = np.swapaxes(vol, 0, 1)

    return (vol, (scale, shift, domain, bounding_box))
    # return (scale, shift, domain, bounding_box)

# def render_slice(mesh, height, domain, return_dict, isMeshAShell):
#     lines = slice.toIntersectingLines(mesh, height)
#     prepixel = np.zeros((domain[0], domain[1]), dtype=bool)
#     perimeter.linesToVoxels(lines, prepixel, isMeshAShell)
#     return_dict[height] = prepixel

def render_slice(height, mesh, domain, isMeshAShell):
    lines = slice.toIntersectingLines(mesh, height)
    prepixel = np.zeros((domain[0], domain[1]), dtype=bool)
    perimeter.linesToVoxels(lines, prepixel, isMeshAShell)
    return prepixel

def import_stl_file(inputFile):
    imported_mesh = mesh.Mesh.from_file(inputFile)
    return get_stl_faces(imported_mesh)

def get_stl_faces(mesh):
    for i, j, k in zip(mesh.v0, mesh.v1, mesh.v2):
        yield (tuple(i), tuple(j), tuple(k))


if __name__ == '__main__':
    #voxelize('Files/Mesh_30000_faces/Mesh_21000.stl', 'Files/meshmesh.nrrd', pow(189, 3))
    #voxelize('test2.stl', 'test2.nrrd', pow(200,3))

    import sys
    if len(sys.argv) < 4:
        print("Usage:", sys.argv[0], "input.nrrd output.nrrd targetElementNum")        
        sys.exit(-1)

    from nrrd import write

    outputVolume, domainData = voxelize(sys.argv[1], int(sys.argv[3])) 

    write(sys.argv[2], outputVolume)