"""Coil porosity per voxel, and the random choice of which voxels become coil voxels.

Used by main.py when the config has "stent_porosity_subvoxels".

1. coilSolidFraction: voxelize the coil STL on a grid that is `subvoxels` times finer
   than the CFD grid (only around the coil, to keep it fast), then average the fine
   voxels inside each CFD voxel. The result is the solid volume fraction (svf) of each
   CFD voxel; its porosity is N = 1 - svf.

2. drawStentVoxels: a fluid voxel becomes a coil voxel with probability 1 - N (= svf).
   Wall, outside and opening voxels are never chosen, so the coil can not replace the
   one-voxel-thick vessel wall and open a hole to the outside of the vessel.

The idea of a per-voxel solid fraction / porosity comes from CoilPathPlanning
scripts/analysis/coil_porosity_statistics.py (commit fbe7110). That script estimates svf
from the coil centerline and wire radius; here it is measured from the deployed coil STL
instead, so no wire radius is needed.
"""

import numpy as np

from voxelizeStl import voxelize, import_stl_file


def coilSolidFraction(stentGeomFile, domainData, volShape, subvoxels):
    """Solid volume fraction of the coil in every CFD voxel, shape `volShape`."""
    scale = domainData[0]
    k = int(subvoxels)
    if k < 1:
        raise ValueError("stent_porosity_subvoxels must be >= 1, got %s" % subvoxels)

    fineSolid, lo, hi = fineCoilVoxels(stentGeomFile, domainData, k)
    block = averageFineBlocks(fineSolid, hi - lo, k, windowShift(k))

    svf = np.zeros(volShape, dtype=np.float64)
    # CFD padded index P <-> CFD coordinate P - 1
    svf[lo[0] + 1:hi[0] + 1, lo[1] + 1:hi[1] + 1, lo[2] + 1:hi[2] + 1] = block

    fineVolume = fineSolid.sum() / (scale[0] * k) ** 3
    print("-> coil volume from fine voxels: %.4f mm^3 (sum of svf: %.4f mm^3)"
          % (fineVolume, svf.sum() / scale[0] ** 3))
    return svf


def fineCoilVoxels(stentGeomFile, domainData, k):
    """Filled coil voxels on a grid k times finer than the CFD grid, around the coil only.

    Returns the fine boolean array and the coil box [lo, hi) in CFD voxel coordinates.
    """
    scale, shift, domain, boundingBox = domainData

    # The coil's bounding box in CFD voxel coordinates, plus one voxel of margin.
    points = np.array(list(import_stl_file(stentGeomFile))).reshape(-1, 3)
    coarse = (points + np.array(shift)) * np.array(scale)
    lo = np.maximum(np.floor(coarse.min(axis=0)).astype(int) - 1, 0)
    hi = np.minimum(np.ceil(coarse.max(axis=0)).astype(int) + 1, np.array(domain, dtype=int))
    print("-> coil spans CFD voxels", lo, "to", hi, "; sub-voxels per voxel edge:", k)

    # The fine grid starts one CFD voxel before `lo`, so the centred window never
    # reaches below index 0.
    fineScale = [s * k for s in scale]
    fineShift = [shift[i] - (lo[i] - 1) / scale[i] for i in range(3)]
    fineDomain = [int(n) for n in (hi - lo + 2) * k]
    fineDomainData = (fineScale, fineShift, fineDomain, boundingBox)

    # Solid (filled) voxelization along three projections; keep a fine voxel when at
    # least two of the three agree, so a single wrong ray can not add a streak.
    votes = None
    for projection, rotation in ((1, -1), (2, 0), (3, 1)):
        print("-> voxelizing coil on the fine grid, projection #%d" % projection)
        fine, _ = voxelize(stentGeomFile, None, False, fineDomainData, rotation)
        fine = fine.astype(np.int8)
        votes = fine if votes is None else votes + fine
    fineSolid = votes >= 2
    return fineSolid, lo, hi


def windowShift(k):
    """How many fine voxels the averaging window is moved back from the CFD voxel's start.

    The voxelizer samples each voxel at a point, so the k fine samples of a CFD voxel are
    centred on its sample point when the window is moved back by (k - 1) // 2. Exactly
    centred for odd k, half a fine voxel off for even k. Found by comparing with the coil
    voxelized directly on the CFD grid (best match at shift 1 for k = 4, 3 for k = 8).
    """
    return -((k - 1) // 2)


def averageFineBlocks(fineSolid, n, k, shift):
    """Mean of the k^3 fine voxels belonging to each of the n[0] x n[1] x n[2] CFD voxels.

    Fine padded index F <-> fine coordinate F - 1, and the fine grid starts one CFD voxel
    before the coil box, so CFD voxel number i of the box starts at fine index
    1 + (i + 1) * k, moved by `shift` fine voxels.
    """
    start = 1 + k + shift
    block = fineSolid[start:start + n[0] * k,
                      start:start + n[1] * k,
                      start:start + n[2] * k]
    return block.reshape(n[0], k, n[1], k, n[2], k).mean(axis=(1, 3, 5))


SELECT_MODES = ("draw", "threshold", "touched")


def selectStentVoxels(svf, geometryFlag, mode, seed=0, threshold=0.5):
    """Pick the coil voxels from the solid fraction, fluid voxels only.

    draw      : picked with probability svf (= 1 - N), see drawStentVoxels
    threshold : picked when svf >= threshold (0.5: the voxel is at least half wire)
    touched   : picked when svf > 0 (every voxel the wire reaches)
    """
    if mode == "draw":
        return drawStentVoxels(svf, geometryFlag, seed)
    if mode == "threshold":
        if not 0.0 < threshold <= 1.0:
            raise ValueError("stent_porosity_threshold must be in (0, 1], got %s" % threshold)
        chosen = svf >= threshold
        rule = "svf >= %g" % threshold
    elif mode == "touched":
        chosen = svf > 0.0
        rule = "svf > 0"
    else:
        raise ValueError("stent_porosity_mode must be one of %s, got %r" % (SELECT_MODES, mode))

    fluid = geometryFlag == 2
    stent = fluid & chosen
    print("-> rule: %s" % rule)
    print("-> voxels touched by the coil (svf > 0):    %8d" % np.count_nonzero(svf > 0.0))
    print("-> voxels meeting the rule:                 %8d" % np.count_nonzero(chosen))
    print("   of which not fluid (wall/opening, kept): %8d" % np.count_nonzero(chosen & ~fluid))
    print("-> coil voxels:                             %8d" % np.count_nonzero(stent))
    return stent


def drawStentVoxels(svf, geometryFlag, seed):
    """Pick coil voxels: a fluid voxel with porosity N is picked with probability 1 - N."""
    fluid = geometryFlag == 2
    candidates = fluid & (svf > 0.0)

    rng = np.random.default_rng(seed)
    draw = rng.random(svf.shape)
    stent = candidates & (draw < svf)

    print("-> random seed:", seed)
    print("-> voxels touched by the coil (svf > 0):    %8d" % np.count_nonzero(svf > 0.0))
    print("   of which not fluid (wall/opening, kept): %8d" % np.count_nonzero((svf > 0.0) & ~fluid))
    print("   fluid candidates:                        %8d" % np.count_nonzero(candidates))
    print("-> expected coil voxels, sum(svf) on fluid: %10.1f" % svf[candidates].sum())
    print("-> drawn coil voxels:                       %8d" % np.count_nonzero(stent))
    return stent
