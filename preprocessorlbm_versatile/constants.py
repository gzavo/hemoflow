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
    def INLET_SVC():
        return 10
    @constant
    def INLET_IVC():
        return 11
    @constant
    def INLET_RHV():
        return 12
    @constant
    def INLET_MHV():
        return 13
    @constant
    def INLET_LHV():
        return 14
    @constant
    def OUTLET_RPA():
        return 15
    @constant
    def OUTLET_LPA():
        return 16
    @constant
    def OUTLET_SIDE_V():
        return 17
    @constant
    def OUTLET_SIDE_P():
        return 18
