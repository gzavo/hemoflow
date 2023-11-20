import sys
import numpy as np
from constants import const
import pickle
import time
import multiprocessing
# import pyvista as pv
# from operator import itemgetter

CONSTANTS = const()

def has_neighboring_unused_voxel(x, y, z, data):
    if data[z - 1][x][y] == CONSTANTS.UNUSED_VOXEL:
        return True
    elif data[z + 1][x][y] == CONSTANTS.UNUSED_VOXEL:
        return True
    elif data[z][x - 1][y] == CONSTANTS.UNUSED_VOXEL:
        return True
    elif data[z][x + 1][y] == CONSTANTS.UNUSED_VOXEL:
        return True
    elif data[z][x][y - 1] == CONSTANTS.UNUSED_VOXEL:
        return True
    elif data[z][x][y + 1] == CONSTANTS.UNUSED_VOXEL:
        return True
    return False

def get_all_inlet_outlet_fluid_voxels(walls, data):
    for z, x, y in walls: 
        if data[z - 1][x][y] == CONSTANTS.FLUID_VOXEL:
            if has_neighboring_unused_voxel(x, y, z - 1, data):
                yield (z - 1, x, y)
        elif data[z + 1][x][y] == CONSTANTS.FLUID_VOXEL:
            if has_neighboring_unused_voxel(x, y, z + 1, data):
                yield (z + 1, x, y)
        elif data[z][x - 1][y] == CONSTANTS.FLUID_VOXEL:
            if has_neighboring_unused_voxel(x - 1, y, z, data):
                yield (z, x - 1, y)
        elif data[z][x + 1][y] == CONSTANTS.FLUID_VOXEL:
            if has_neighboring_unused_voxel(x + 1, y, z, data):
                yield (z, x + 1, y)
        elif data[z][x][y - 1] == CONSTANTS.FLUID_VOXEL:
            if has_neighboring_unused_voxel(x, y - 1, z, data):
                yield (z, x, y - 1)
        elif data[z][x][y + 1] == CONSTANTS.FLUID_VOXEL:
            if has_neighboring_unused_voxel(x, y + 1, z, data):
                yield (z, x, y + 1)

def get_neighbouring_unused_voxel_plane(x, y, z, data):
    if data[z - 1][x][y] == CONSTANTS.UNUSED_VOXEL or data[z + 1][x][y] == CONSTANTS.UNUSED_VOXEL:
        return "z"
    elif data[z][x - 1][y] == CONSTANTS.UNUSED_VOXEL or data[z][x + 1][y] == CONSTANTS.UNUSED_VOXEL:
        return "x"
    elif data[z][x][y - 1] == CONSTANTS.UNUSED_VOXEL or data[z][x][y + 1] == CONSTANTS.UNUSED_VOXEL:
        return "y"

def find_inlet_outlet(x, y, z, data, plane):
    """ Returns a list of voxels (x, y, z) that make up the inlet/outlet for a voxel """
    inlet_outlet = []
    inlet_outlet += recur_find_inlet_outlet(x, y, z, None, None, None, [], data, plane)
    
    if len(inlet_outlet) > 1:
        print(" -> Found opening of size:", len(inlet_outlet))
        return inlet_outlet

    return None

def show_progress(progress, total):
    sys.stdout.flush()
    sys.stdout.write('\r[{0}] {1:.2f}%'.format('#' * int((progress * 10) / total), (progress * 100.) / total))
    sys.stdout.flush()

def voxel_is_not_yet_processed(x, y, z, inlets_outlets):
    for l in inlets_outlets:
        if (x, y, z) in l:
            return False
    return True

def recur_find_inlet_outlet(x, y, z, prev_x, prev_y, prev_z, result, data, plane):
    """ Recursive function for finding the inlet/outlet for a voxel """
    result.append((x, y, z))

    dz_list = [-1, 0, 1]
    if plane == "z":
        dz_list = [0]

    dx_list = [-1, 0, 1]
    if plane == "x":
        dx_list = [0]

    dy_list = [-1, 0, 1]
    if plane == "y":
        dy_list = [0]

    for dz in dz_list:
        for dx in dx_list:
            for dy in dy_list:
                if (dx, dy, dz) != (0, 0, 0):
                    x_n = x + dx
                    y_n = y + dy
                    z_n = z + dz
                    if ((x_n, y_n, z_n) != (prev_x, prev_y, prev_z) and
                            data[z_n][x_n][y_n] == CONSTANTS.FLUID_VOXEL and
                            has_neighboring_unused_voxel(x_n, y_n, z_n, data) and
                            (x_n, y_n, z_n) not in result):
                        result = recur_find_inlet_outlet(x_n, y_n, z_n, x, y, z, result, data, plane)
    return result

def find_openings(fluid_voxel, data):
    z, x, y = fluid_voxel
    plane = get_neighbouring_unused_voxel_plane(x, y, z, data)
    inlet_outlet = find_inlet_outlet(x, y, z, data, plane)
    # if(inlet_outlet is not None):
    return inlet_outlet

def detect_inlets_outlets(data):
    """ Detects all inlets and outlets of the vessel model. Assumes the following:
        - The vessel has only one inlet.
        - All outlets are smaller in diameter than the inlets.
        - The inlets and outlets are cut parallel to the x, y or z plane.
        - The openings are located on the boundaries
    """
    print("Detecting all fluid voxels surrounded by an unused voxel...")
    #print(len(data) * len(data[0]) * len(data[0][0]))
    
    start = time.time()
    
    
    ############################ Method 1: Connected component analysis (very fast) ############################
    # Note: Only works if all openings lying on the boundary surfces
    # Drawbacks: Problems when detecting edge openings
    
    # boundaries = []
    # # Extract the boundary 2d matrix
    # # Because we know the openings can only locate on the boundaries
    # # Note data is (z,x,y), but inlet_outlet should be (x,y,z)
    # boundaries.append(data[0,:,:])
    # boundaries.append(data[-1,:,:])
    # boundaries.append(data[:,0,:])
    # boundaries.append(data[:,-1,:])
    # boundaries.append(data[:,:,0])
    # boundaries.append(data[:,:,-1])
    # # Validation
    # print(len(boundaries[0]), len(boundaries[2]), len(boundaries[4]))
    
    # # Optmized Connected-component labeling
    # inlets_outlets = []
    # from scipy.ndimage import label
    
    # for i, boundary in enumerate(boundaries):
    #     labeled_array, num_features = label(boundary)
    #     print(boundary.shape, num_features)
    #     for j in range(1, num_features+1):
    #         # A 2xN indexing array, first row is x indices, second row is y indices
    #         inlet_outlet = np.where(labeled_array == j)
    #         if i == 0:
    #             inlet_outlet = np.vstack((inlet_outlet, np.repeat(0, inlet_outlet[0].shape))) # (x,y) -> (x,y,z)
    #         elif i == 1:
    #             inlet_outlet = np.vstack((inlet_outlet, np.repeat(data.shape[0]-1, inlet_outlet[0].shape))) # (x,y) -> (x,y,z)
    #         elif i == 2:
    #             inlet_outlet = np.vstack((np.repeat(0, inlet_outlet[0].shape), inlet_outlet[1:], inlet_outlet[:1])) # (z,y) -> (x,y,z)
    #         elif i == 3:
    #             inlet_outlet = np.vstack((np.repeat(data.shape[1]-1, inlet_outlet[0].shape), inlet_outlet[1:], inlet_outlet[:1])) # (z,y) -> (x,y,z)
    #         elif i == 4:
    #             inlet_outlet = np.vstack((inlet_outlet[1:], np.repeat(0, inlet_outlet[0].shape), inlet_outlet[:1])) # (z,x) -> (x,y,z)
    #         else:
    #             inlet_outlet = np.vstack((inlet_outlet[1:], np.repeat(data.shape[2]-1, inlet_outlet[0].shape), inlet_outlet[:1])) # (z,x) -> (x,y,z)
                
    #         inlets_outlets.append(np.transpose(inlet_outlet))
    
    
    ############################ Method 2: Iterate through all surface voxels (slow) ############################
    # Note: Only work if all openings are parallel to normal planes.
    # Drawbacks: Very slow, large memory consumption
    
    wall_voxels = zip(*np.where(data == CONSTANTS.WALL_VOXEL))
    found_fluid_voxels = list(get_all_inlet_outlet_fluid_voxels(wall_voxels, data))

    inlets_outlets = []
    z = 0
    for z, x, y in found_fluid_voxels: 
        if (voxel_is_not_yet_processed(x, y, z, inlets_outlets)):
            plane = get_neighbouring_unused_voxel_plane(x, y, z, data)
            inlet_outlet = find_inlet_outlet(x, y, z, data, plane)
            # if(inlet_outlet is not None):
            inlets_outlets.append(inlet_outlet)
    
    end = time.time()
    print(f"Find opening time: {end - start} sec")
            
    print(len(inlets_outlets))

    return paint_inlets_outlets(inlets_outlets, data)

def paint_inlets_outlets(inlets_outlets, data):
    """ Determines the inlet by picking the inlet/outlet with the greatest area
        Then paints the inlet with the value 6 and the outlets with 7
    """

    openingIdx = [] # Label identifying this opening
    openingCenter = []
    # openingQ = [] # Volume ratio assigned to this opening based on area (TODO: angle or R needed for this one)
    
    data_result = np.copy(data)
    
    ''' Label the openings based on their area size '''
    areas = np.array([len(x) for x in inlets_outlets])
    # Sort the inlets_outlets based on their area size
    inlets_outlets = [x for _,x in sorted(zip(list(areas), inlets_outlets))]
    
    # Only uncomment when using Method 1: Connected component analysis.
    # # The rationale is, using connected component analysis, the algorithm will detect two clusters
    # # for the opening at the domain edge, so we need to remove it.
    # # Note: flawed, use carefully
    # inlets_outlets.pop(0)
    
    
    # Hardcoded
    '''
    inlet 6: SVC
    inlet 7: IVC
    inlet 4: RHV
    inlet 3: MHV
    inlet 2: LHV
    '''
    
    '''
    outlet 5: LPA
    outlet 8: RPA
    outlet 0/1: RPA side branches
    '''
    
    paint_dic = {6: CONSTANTS.INLET_SVC, 7: CONSTANTS.INLET_IVC, 4: CONSTANTS.INLET_RHV, 3: CONSTANTS.INLET_MHV,
                 2: CONSTANTS.INLET_LHV, 5: CONSTANTS.OUTLET_LPA, 8: CONSTANTS.OUTLET_RPA, 0: CONSTANTS.OUTLET_SIDE_P,
                 1: CONSTANTS.OUTLET_SIDE_V}

    print("Number of inlets: ", 5)
    print("Number of pressure outlets: 0")
    print("Number of velocity outlet(s): ", 4)
    
    for i, opening in enumerate(inlets_outlets):
        openingC = np.zeros(3)
        # Find the corresponding paint id
        paint_id = paint_dic[i]
        openingIdx.append(paint_id)
        
        for (x, y, z) in opening:
            data_result[z][x][y] = paint_id
            openingC += np.array((z,x,y))
            
        openingCenter.append(openingC / len(opening))

    
    return (openingIdx, openingCenter, inlets_outlets, data_result)


def detectOpenings(inputArray, wall_path):
    sys.setrecursionlimit(10000000)
    
    global wall_nrrd_path
    wall_nrrd_path = wall_path

    # # Append surrounding layer to avoid boundary checking (THIS SHOULD NOT BE NECESSARY! - > FLUID ON THE EDGE OF THE DOMAIN!)
    # # Note that if you use Method 1 to detect opening voxels, you can not append surrounding layers, just comment the following line out
    data = np.pad(inputArray, 1, 'constant')
    data = inputArray

    # Detect openings
    opening_ids, openingCenter, inlets_outlets, result = detect_inlets_outlets(data)

    # Return result without outer layer (For Method 2)
    return (opening_ids, openingCenter, inlets_outlets, result[1:-1, 1:-1, 1:-1].astype('int8', copy=False) )
    # (For Method 1)
    # return (opening_ids, openingCenter, inlets_outlets, result.astype('int8', copy=False) )


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage:", sys.argv[0], "input.nrrd output.nrrd")
        sys.exit(-1) 

    import time
    import nrrd
    
    startTime = time.time()
    #data, data_header = nrrd.read('Files/reference_geom_capped.nrrd')
    data, data_header = nrrd.read(sys.argv[1])
    
    openingDescr, result = detectOpenings(data)
    
    # Write result without the surrounding layer
    nrrd.write(sys.argv[2], result)

    endTime = time.time()
    timeElapsed = int(round((endTime - startTime) * 1000))
    print("Elapsed time:", timeElapsed, "[ms]")