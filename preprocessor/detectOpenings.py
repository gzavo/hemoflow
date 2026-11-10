"""Opening voxels from the caps of the vessel surface.

The vessel surface is closed with flat caps on faces of its bounding box. After the
wall has been added (createFluidSolid.add_walls), every cap is covered by a layer of
wall voxels on the outermost layer of the domain. Only the wall voxels in front of the
cap (centre inside the cap, fluid behind it) are turned into opening voxels; the rest
of that layer stays wall, so other vessels close to the same face keep their wall.
"""
import numpy as np
from scipy import ndimage

from constants import const

CONSTANTS = const()


def opening_voxels(labels, offset, grid, cap):
    """Indices (k, 3) of the opening voxels of one cap on the labels grid."""
    axis, side = cap.axis, cap.side
    a1, a2 = [a for a in range(3) if a != axis]
    layer = 0 if side == 0 else labels.shape[axis] - 1
    inward = layer + 1 if side == 0 else layer - 1

    # The outermost layer must lie just outside the cap plane
    layerCenter = grid.centers(axis, layer + offset[axis])
    capPlane = cap.triangles[0, 0, axis]
    if abs(layerCenter - capPlane) > grid.dx:
        raise RuntimeError(f"Cap on face {cap.face} is not on the outermost voxel layer "
                           f"(layer centre {layerCenter:.4g}, cap plane {capPlane:.4g}).")

    i1 = np.arange(labels.shape[a1])
    i2 = np.arange(labels.shape[a2])
    u, v = np.meshgrid(grid.centers(a1, i1 + offset[a1]), grid.centers(a2, i2 + offset[a2]), indexing="ij")
    inside = cap.contains(u, v)

    def plane(index):
        s = [slice(None)] * 3
        s[axis] = index
        return labels[tuple(s)]

    candidate = inside & (plane(layer) == CONSTANTS.WALL_VOXEL) & (plane(inward) == CONSTANTS.FLUID_VOXEL)
    j1, j2 = np.nonzero(candidate)
    voxels = np.zeros((len(j1), 3), dtype=int)
    voxels[:, axis] = layer
    voxels[:, a1] = j1
    voxels[:, a2] = j2
    return voxels, candidate


def paint_openings(labels, offset, grid, caps, capLabels):
    """Paint the opening voxels of every cap with its label.

    Returns, per cap, the number of voxels, the number of connected pieces and the
    centre of the opening (voxel indices of the labels grid)."""
    painted = labels.copy()
    info = []
    for cap, label in zip(caps, capLabels):
        voxels, mask2d = opening_voxels(labels, offset, grid, cap)
        if len(voxels) == 0:
            raise RuntimeError(f"Cap on face {cap.face} at {np.round(cap.centroid, 3).tolist()} has no opening voxels; "
                               "the vessel does not reach the bounding-box face there.")
        current = painted[tuple(voxels.T)]
        if np.any(current != CONSTANTS.WALL_VOXEL):
            raise RuntimeError(f"Openings overlap at the cap on face {cap.face}.")
        painted[tuple(voxels.T)] = label
        pieces = ndimage.label(mask2d, structure=np.ones((3, 3), dtype=bool))[1]
        info.append({"voxels": len(voxels), "pieces": pieces, "center": voxels.mean(axis=0)})
    return painted, info
