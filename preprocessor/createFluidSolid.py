"""Wall layer around the voxelized fluid.

Labels: 0 unused (outside), 1 wall, 2 fluid (see constants.py).
"""
import numpy as np
from scipy import ndimage

from constants import const

CONSTANTS = const()


def add_walls(fluid):
    """Surround the fluid with a one-voxel wall and crop the empty outer layers.

    The fluid grid is padded by one voxel on every side so that the wall fits, the
    wall is every non-fluid voxel that touches a fluid voxel (26-neighbourhood), and
    layers that contain neither wall nor fluid are removed.

    Returns (labels, offset): labels is an int16 array; voxel 'index' of labels
    corresponds to voxel 'index + offset' of the fluid grid (offset can be -1)."""
    padded = np.pad(fluid, 1, mode="constant", constant_values=False)
    wall = ndimage.binary_dilation(padded, structure=np.ones((3, 3, 3), dtype=bool)) & ~padded

    labels = np.full(padded.shape, CONSTANTS.UNUSED_VOXEL, dtype=np.int16)
    labels[wall] = CONSTANTS.WALL_VOXEL
    labels[padded] = CONSTANTS.FLUID_VOXEL

    used = labels != CONSTANTS.UNUSED_VOXEL
    lo, hi = [], []
    for axis in range(3):
        other = tuple(a for a in range(3) if a != axis)
        layers = np.where(used.any(axis=other))[0]
        lo.append(layers[0])
        hi.append(layers[-1] + 1)
    labels = labels[lo[0]:hi[0], lo[1]:hi[1], lo[2]:hi[2]]
    offset = np.array(lo) - 1  # padded index p corresponds to fluid-grid index p - 1
    return np.ascontiguousarray(labels), offset


def crop_like(array, offset, shape):
    """Bring an array defined on the fluid grid onto the labels grid of add_walls()."""
    out = np.zeros(shape, dtype=array.dtype)
    src, dst = [], []
    for axis in range(3):
        start = offset[axis]                       # fluid-grid index of labels index 0
        s0, s1 = max(start, 0), min(start + shape[axis], array.shape[axis])
        src.append(slice(s0, s1))
        dst.append(slice(s0 - start, s1 - start))
    out[tuple(dst)] = array[tuple(src)]
    return out


def _touches(mask, borderValue):
    """Voxels with a face neighbour in 'mask'; positions outside the array count as borderValue."""
    padded = np.pad(mask, 1, mode="constant", constant_values=borderValue)
    out = np.zeros(mask.shape, dtype=bool)
    for axis in range(3):
        for shift in (-1, 1):
            out |= np.roll(padded, shift, axis=axis)[1:-1, 1:-1, 1:-1]
    return out


def leak_check(labels):
    """Number of leaking voxels: fluid voxels next to an outside voxel or the domain
    border, and opening voxels next to an outside voxel. Should be 0."""
    unused = labels == CONSTANTS.UNUSED_VOXEL
    fluid = labels >= CONSTANTS.FLUID_VOXEL
    opening = (labels >= CONSTANTS.INLET_VOXEL) & (labels != CONSTANTS.FLOW_DIVERTER_VOXEL)
    leaks = (fluid & ~opening & _touches(unused, True)) | (opening & _touches(unused, False))
    return int(leaks.sum())


def fluid_components(labels):
    """Sizes of the face-connected components of the fluid (incl. openings), largest first."""
    components, n = ndimage.label(labels >= CONSTANTS.FLUID_VOXEL)
    sizes = np.bincount(components.ravel())[1:]
    return sorted(sizes.tolist(), reverse=True)
