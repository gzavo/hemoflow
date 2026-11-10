"""Voxelization of the vessel surface and of flow diverter shells (VTK based).

All lengths are in the units of the STL files (millimetres). Voxel (i, j, k) of a
VoxelGrid is the cell [corner + (i, j, k) * dx, corner + (i + 1, j + 1, k + 1) * dx);
its centre is at corner + (i + 0.5, j + 0.5, k + 0.5) * dx.
"""
import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy


class VoxelGrid:
    """Uniform grid of cubic voxels covering a bounding box."""

    def __init__(self, corner, dx, shape):
        self.corner = np.asarray(corner, dtype=float)
        self.dx = float(dx)
        self.shape = tuple(int(n) for n in shape)

    def centers(self, axis, indices):
        return self.corner[axis] + (np.asarray(indices) + 0.5) * self.dx

    def __repr__(self):
        return f"VoxelGrid(shape={self.shape}, dx={self.dx:.6g}, corner={self.corner.tolist()})"


def read_surface(fileName):
    """Read an STL file as a clean triangle mesh."""
    reader = vtk.vtkSTLReader()
    reader.SetFileName(fileName)
    reader.Update()
    if reader.GetOutput().GetNumberOfCells() == 0:
        raise RuntimeError(f"Cannot read any triangles from {fileName}")
    clean = vtk.vtkCleanPolyData()
    clean.SetInputConnection(reader.GetOutputPort())
    triangles = vtk.vtkTriangleFilter()
    triangles.SetInputConnection(clean.GetOutputPort())
    triangles.Update()
    out = vtk.vtkPolyData()
    out.DeepCopy(triangles.GetOutput())
    return out


def bounds(polyData):
    """Bounding box as a (3, 2) array of [min, max] per axis."""
    return np.array(polyData.GetBounds()).reshape(3, 2)


def boundary_loops(polyData):
    """Open boundary loops (holes) of a surface, as a list of (n, 3) point arrays."""
    edges = vtk.vtkFeatureEdges()
    edges.SetInputData(polyData)
    edges.BoundaryEdgesOn()
    edges.FeatureEdgesOff()
    edges.NonManifoldEdgesOff()
    edges.ManifoldEdgesOff()
    stripper = vtk.vtkStripper()
    stripper.SetInputConnection(edges.GetOutputPort())
    stripper.JoinContiguousSegmentsOn()
    stripper.Update()
    out = stripper.GetOutput()
    if out.GetNumberOfLines() == 0:
        return []
    points = vtk_to_numpy(out.GetPoints().GetData()).astype(float)
    loops = []
    for i in range(out.GetNumberOfLines()):
        ids = out.GetCell(i).GetPointIds()
        loops.append(points[[ids.GetId(k) for k in range(ids.GetNumberOfIds())]])
    return loops


def close_holes(polyData):
    """Close open boundary loops (open vessel ends) with planar caps."""
    fill = vtk.vtkFillHolesFilter()
    fill.SetInputData(polyData)
    fill.SetHoleSize(1e30)
    triangles = vtk.vtkTriangleFilter()
    triangles.SetInputConnection(fill.GetOutputPort())
    triangles.Update()
    out = vtk.vtkPolyData()
    out.DeepCopy(triangles.GetOutput())
    return out


def triangles_of(polyData):
    """Triangle vertex coordinates as an (n, 3, 3) array."""
    points = vtk_to_numpy(polyData.GetPoints().GetData()).astype(float)
    connectivity = vtk_to_numpy(polyData.GetPolys().GetConnectivityArray()).reshape(-1, 3)
    return points, connectivity


def _edge_connected_groups(cells):
    """Groups of triangles (given as vertex index triples) connected through shared edges."""
    parent = list(range(len(cells)))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    owner = {}
    for t, (a, b, c) in enumerate(cells):
        for edge in ((a, b), (b, c), (c, a)):
            edge = (min(edge), max(edge))
            if edge in owner:
                parent[find(t)] = find(owner[edge])
            else:
                owner[edge] = t
    groups = {}
    for t in range(len(cells)):
        groups.setdefault(find(t), []).append(t)
    return list(groups.values())


class Cap:
    """A flat cap of the closed vessel surface on a face of its bounding box (an opening)."""

    def __init__(self, axis, side, triangles):
        self.axis = axis          # 0, 1, 2
        self.side = side          # 0: min face, 1: max face
        self.triangles = triangles  # (n, 3, 3) vertex coordinates
        cross = np.cross(triangles[:, 1] - triangles[:, 0], triangles[:, 2] - triangles[:, 0])
        areas = 0.5 * np.linalg.norm(cross, axis=1)
        self.area = areas.sum()
        self.centroid = (triangles.mean(axis=1) * areas[:, None]).sum(axis=0) / self.area

    @property
    def face(self):
        return "xyz"[self.axis] + ("-" if self.side == 0 else "+")

    def contains(self, u, v):
        """Which 2D points (coordinates of the two other axes) lie inside the cap."""
        a1, a2 = [a for a in range(3) if a != self.axis]
        inside = np.zeros(np.shape(u), dtype=bool)
        for tri in self.triangles:
            (x0, y0), (x1, y1), (x2, y2) = tri[:, [a1, a2]]
            det = (y1 - y2) * (x0 - x2) + (x2 - x1) * (y0 - y2)
            if det == 0.0:
                continue
            l0 = ((y1 - y2) * (u - x2) + (x2 - x1) * (v - y2)) / det
            l1 = ((y2 - y0) * (u - x2) + (x0 - x2) * (v - y2)) / det
            l2 = 1.0 - l0 - l1
            eps = -1e-9
            inside |= (l0 >= eps) & (l1 >= eps) & (l2 >= eps)
        return inside


def _triangle_neighbours(cells, numPoints):
    """Pairs (a, b) of triangles that share an edge."""
    edges = np.sort(cells[:, [[0, 1], [1, 2], [2, 0]]].reshape(-1, 2), axis=1)
    owner = np.repeat(np.arange(len(cells)), 3)
    key = edges[:, 0].astype(np.int64) * numPoints + edges[:, 1]
    order = np.argsort(key, kind="stable")
    same = key[order][1:] == key[order][:-1]
    return owner[order[:-1][same]], owner[order[1:][same]]


def find_caps(polyData, tolerance, maxWallAlignment=0.95):
    """Flat caps of a closed surface that lie on faces of its bounding box.

    A cap candidate is a group of edge-connected triangles whose vertices all lie on the
    same bounding-box face (within 'tolerance'). It is a cap (an opening) if the wall
    around it meets it at a clear angle: the median |normal . face normal| of the
    neighbouring triangles is below maxWallAlignment. Flat strips of a vessel wall that
    happen to be flush with the face continue smoothly into the wall and are rejected.

    Returns (caps, rejected) where rejected lists the faces of the rejected patches."""
    points, cells = triangles_of(polyData)
    box = bounds(polyData)
    tri = points[cells]
    normals = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
    lengths = np.linalg.norm(normals, axis=1)
    normals = normals / np.where(lengths > 0, lengths, 1.0)[:, None]
    pairA, pairB = _triangle_neighbours(cells, len(points))

    caps, rejected = [], []
    for axis in range(3):
        for side in (0, 1):
            onFace = np.all(np.abs(tri[:, :, axis] - box[axis, side]) <= tolerance, axis=1)
            if not onFace.any():
                continue
            faceCells = np.where(onFace)[0]
            for group in _edge_connected_groups(cells[faceCells]):
                members = np.zeros(len(cells), dtype=bool)
                members[faceCells[group]] = True
                border = np.concatenate([pairB[members[pairA] & ~members[pairB]],
                                         pairA[members[pairB] & ~members[pairA]]])
                alignment = np.median(np.abs(normals[border, axis])) if len(border) else 1.0
                cap = Cap(axis, side, tri[faceCells[group]])
                if alignment < maxWallAlignment:
                    caps.append(cap)
                else:
                    rejected.append(cap)
    return caps, rejected


def make_grid(polyData, targetElements, target_dx=None):
    """Voxel grid over the bounding box of the surface.

    The resolution comes from target_dx [mm] if given, otherwise from the number of
    voxels targetElements in the bounding box. dx is then adjusted so that the x extent
    is a whole number of voxels; the grid covers the full extent along y and z."""
    box = bounds(polyData)
    extent = box[:, 1] - box[:, 0]
    if target_dx is not None:
        voxelsPerUnit = 1.0 / target_dx
    else:
        voxelsPerUnit = (targetElements / np.prod(extent)) ** (1.0 / 3.0)
    nx = max(int(voxelsPerUnit * extent[0]), 1)
    dx = extent[0] / nx
    shape = [nx] + [max(int(np.ceil(extent[a] / dx - 1e-9)), 1) for a in (1, 2)]
    return VoxelGrid(box[:, 0], dx, shape)


def voxelize_closed(polyData, grid):
    """Voxels whose centres lie inside a closed surface, as a boolean (nx, ny, nz) array."""
    nx, ny, nz = grid.shape
    stencil = vtk.vtkPolyDataToImageStencil()
    stencil.SetInputData(polyData)
    stencil.SetOutputOrigin(*(grid.corner + 0.5 * grid.dx))
    stencil.SetOutputSpacing(grid.dx, grid.dx, grid.dx)
    stencil.SetOutputWholeExtent(0, nx - 1, 0, ny - 1, 0, nz - 1)
    image = vtk.vtkImageStencilToImage()
    image.SetInputConnection(stencil.GetOutputPort())
    image.SetInsideValue(1)
    image.SetOutsideValue(0)
    image.SetOutputScalarTypeToUnsignedChar()
    image.Update()
    values = vtk_to_numpy(image.GetOutput().GetPointData().GetScalars())
    return values.reshape(nz, ny, nx).transpose(2, 1, 0).astype(bool)


def voxelize_shell(polyData, grid):
    """Voxels crossed by a (possibly open) surface, as a boolean (nx, ny, nz) array.

    Rays along each axis through the voxel centres are intersected with the surface;
    the voxel containing each intersection is marked. Used for flow diverter meshes."""
    tree = vtk.vtkModifiedBSPTree()
    tree.SetDataSet(polyData)
    tree.BuildLocator()
    box = bounds(polyData)
    dx, corner, shape = grid.dx, grid.corner, grid.shape
    mask = np.zeros(shape, dtype=bool)
    points, cellIds = vtk.vtkPoints(), vtk.vtkIdList()
    for axis in range(3):
        a1, a2 = [a for a in range(3) if a != axis]
        idx1 = [i for i in range(shape[a1]) if box[a1, 0] - dx <= grid.centers(a1, i) <= box[a1, 1] + dx]
        idx2 = [j for j in range(shape[a2]) if box[a2, 0] - dx <= grid.centers(a2, j) <= box[a2, 1] + dx]
        for i in idx1:
            for j in idx2:
                p1 = np.zeros(3)
                p1[a1], p1[a2] = grid.centers(a1, i), grid.centers(a2, j)
                p2 = p1.copy()
                p1[axis], p2[axis] = box[axis, 0] - dx, box[axis, 1] + dx
                points.Reset()
                cellIds.Reset()
                tree.IntersectWithLine(p1, p2, 1e-9, points, cellIds)
                for k in range(points.GetNumberOfPoints()):
                    c = int(np.floor((points.GetPoint(k)[axis] - corner[axis]) / dx))
                    if 0 <= c < shape[axis]:
                        index = [0, 0, 0]
                        index[axis], index[a1], index[a2] = c, i, j
                        mask[tuple(index)] = True
    return mask
