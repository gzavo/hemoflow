import sys
import numpy as np
from constants import const
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

def detect_inlets_outlets(data):
    """ Detects all inlets and outlets of the vessel model. Assumes the following:
        - The vessel has only one inlet.
        - All outlets are smaller in diameter than the inlets.
        - The inlets and outlets are cut parallel to the x, y or z plane.
    """
    print("Detecting all fluid voxels surround by an unused voxel...")
    #print(len(data) * len(data[0]) * len(data[0][0]))
    wall_voxels = zip(*np.where(data == CONSTANTS.WALL_VOXEL))
    found_fluid_voxels = list(get_all_inlet_outlet_fluid_voxels(wall_voxels, data))
    # voxel_count = len(found_fluid_voxels)
    inlets_outlets = []
    count = 1
    z = 0
    for z, x, y in found_fluid_voxels: 
        #show_progress(count, voxel_count)
        if (voxel_is_not_yet_processed(x, y, z, inlets_outlets)):
            plane = get_neighbouring_unused_voxel_plane(x, y, z, data)
            inlet_outlet = find_inlet_outlet(x, y, z, data, plane)
            inlets_outlets.append(inlet_outlet)
        count += 1

    return paint_inlets_outlets(inlets_outlets, data)

def paint_inlets_outlets(inlets_outlets, data):
    """ Determines the inlet by picking the inlet/outlet with the greatest area
        Then paints the inlet with the value 6 and the outlets with 7
    """

    openingIdx = [] # Label identifying this opening
    openingCenter = []
    # openingQ = [] # Volume ratio assigned to this opening based on area (TODO: angle or R needed for this one)

    data_result = np.copy(data)
    
    inlet = max(inlets_outlets, key=len)
    pressureOutlet = min(inlets_outlets, key=len)

    velocityOutlets = [o for o in inlets_outlets if o != inlet and o != pressureOutlet]

    print("Number of inlets: 1")
    print("Number of pressure outlets: 1")
    print("Number of velocity outlet(s): ", len(velocityOutlets))

    openingC = np.zeros(3)
    openingIdx.append(CONSTANTS.INLET_VOXEL)
    for (x, y, z) in inlet:
        data_result[z][x][y] = CONSTANTS.INLET_VOXEL
        openingC += np.array((z,x,y))
    openingCenter.append(openingC / len(inlet))

    openingC = np.zeros(3)
    openingIdx.append(CONSTANTS.OUTLET_VOXEL)
    for (x, y, z) in pressureOutlet:
        data_result[z][x][y] = CONSTANTS.OUTLET_VOXEL
        openingC += np.array((z,x,y))
    openingCenter.append(openingC / len(pressureOutlet))

    outletCount = 0
    for outlet in velocityOutlets:
        openingC = np.zeros(3)
        openingIdx.append(CONSTANTS.OUTLET_REST_VOXEL + outletCount)
      
        for (x, y, z) in outlet:
            data_result[z][x][y] = CONSTANTS.OUTLET_REST_VOXEL + outletCount
            openingC += np.array((z,x,y))

        openingCenter.append(openingC / len(outlet))
        outletCount += 1
    
    return (openingIdx, openingCenter, data_result)

def detectOpenings(inputArray):
    sys.setrecursionlimit(100000)

    # Append surrounding layer to avoid boundary checking (THIS SHOULD NOT BE NECESSARY! - > FLUID ON THE EDGE OF THE DOMAIN!)
    data = np.pad(inputArray, 1, 'constant')

    # Detect openings
    openingIdx, openingCenter, result = detect_inlets_outlets(data)

    # Return result without outer layer
    return (openingIdx, openingCenter, result[1:-1, 1:-1, 1:-1].astype(np.short, copy=False) )


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