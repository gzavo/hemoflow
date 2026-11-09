def constant(f):
    def fset(self, value):
        raise TypeError
    def fget(self):
        return f()
    return property(fget, fset)

# Voxel labels. Openings: 10 = inlet (start of the centerline), 11 = smallest outlet,
# 12, 13, ... = remaining outlets by decreasing radius. 20 = fluid voxels inside a flow diverter.
class const(object):
    @constant
    def UNUSED_VOXEL():
        return 0
    @constant
    def WALL_VOXEL():
        return 1
    @constant
    def FLUID_VOXEL():
        return 2
    @constant
    def INLET_VOXEL():
        return 10
    @constant
    def OUTLET_VOXEL():
        return 11
    @constant
    def OUTLET_REST_VOXEL():
        return 12
    @constant
    def FLOW_DIVERTER_VOXEL():
        return 20
