def constant(f):
    def fset(self, value):
        raise TypeError
    def fget(self):
        return f()
    return property(fget, fset)

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
