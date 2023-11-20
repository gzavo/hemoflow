import math
import itertools

import numpy as np

import perimeter
from util import manhattanDistance, removeDupsFromPointList

def toIntersectingLines(mesh, height):
    relevantTriangles = list(filter(lambda tri: isAboveAndBelow(tri, height), mesh))
    notSameTriangles = filter(lambda tri: not isIntersectingTriangle(tri, height), relevantTriangles)
    lines = list(map(lambda tri: triangleToIntersectingLines(tri, height), notSameTriangles))
    return lines

def drawLineOnPixels(p1, p2, pixels):
    lineSteps = math.ceil(manhattanDistance(p1, p2))
    if lineSteps == 0:
        pixels[int(p1[0]), int(p2[1])] = True
        return
    for j in range(lineSteps + 1):
        point = linearInterpolation(p1, p2, j / lineSteps)
        pixels[int(point[0]), int(point[1])] = True

def linearInterpolation(p1, p2, distance):
    '''
    :param p1: Point 1
    :param p2: Point 2
    :param distance: Between 0 and 1, Lower numbers return points closer to p1.
    :return: A point on the line between p1 and p2
    '''
    slopex = (p1[0] - p2[0])
    slopey = (p1[1] - p2[1])
    slopez = p1[2] - p2[2]
    return (
        p1[0] - distance * slopex,
        p1[1] - distance * slopey,
        p1[2] - distance * slopez
    )


def isAboveAndBelow(pointList, height):
    '''

    :param pointList: Can be line or triangle
    :param height:
    :return: true if any line from the triangle crosses or is on the height line,
    '''
    above = list(filter(lambda pt: pt[2] > height, pointList))
    below = list(filter(lambda pt: pt[2] < height, pointList))
    same = list(filter(lambda pt: pt[2] == height, pointList))
    if len(same) == 3 or len(same) == 2:
        return True
    elif (above and below):
        return True
    else:
        return False

def isIntersectingTriangle(triangle, height):
    assert (len(triangle) == 3)
    same = list(filter(lambda pt: pt[2] == height, triangle))
    return len(same) == 3


def triangleToIntersectingLines(triangle, height):
    assert (len(triangle) == 3)
    above = list(filter(lambda pt: pt[2] > height, triangle))
    below = list(filter(lambda pt: pt[2] < height, triangle))
    same = list(filter(lambda pt: pt[2] == height, triangle))
    assert len(same) != 3
    if len(same) == 2:
        return same[0], same[1]
    elif len(same) == 1:
        side1 = whereLineCrossesZ(above[0], below[0], height)
        return side1, same[0]
    else:
        lines = []
        for a in above:
            for b in below:
                lines.append((b, a))
        side1 = whereLineCrossesZ(lines[0][0], lines[0][1], height)
        side2 = whereLineCrossesZ(lines[1][0], lines[1][1], height)
        return side1, side2


def whereLineCrossesZ(p1, p2, z):
    if (p1[2] > p2[2]):
        t = p1
        p1 = p2
        p2 = t
    # now p1 is below p2 in z
    if p2[2] == p1[2]:
        distance = 0
    else:
        distance = (z - p1[2]) / (p2[2] - p1[2])
    return linearInterpolation(p1, p2, distance)


def calculateScaleAndShift(mesh, targetElements):
    allPoints = [item for sublist in mesh for item in sublist]
    mins = [0, 0, 0]
    maxs = [0, 0, 0]
    ds = [0, 0, 0]
    for i in range(3):
        mins[i] = min(allPoints, key=lambda tri: tri[i])[i]
        maxs[i] = max(allPoints, key=lambda tri: tri[i])[i]
        ds[i] = maxs[i] - mins[i]
    
    bounding_box = [mins, maxs]
    ds3 = ds[0]*ds[1]*ds[2]
    vox_scale = (targetElements / ds3 ) ** (1. / 3)
    
    #domain = map(int, [vox_scale * ds[0], vox_scale * ds[1], vox_scale * ds[2]]) 
    domain = [int(x) for x in [vox_scale * ds[0], vox_scale * ds[1], vox_scale * ds[2]]] 
    shift = [-minimum for minimum in mins]
    
    #xyscale = (domain[0] - 1.0) / ds[0]
    xyscale = domain[0] / ds[0]
    scale = [xyscale, xyscale, xyscale] #TODO Something is fishy here, what is this xyscale???

    return (scale, shift, domain, bounding_box)


def scaleAndShiftMesh(mesh, scale, shift):
    for tri in mesh:
        newTri = []
        for pt in tri:
            newpt = [0, 0, 0]
            for i in range(3):
                newpt[i] = (pt[i] + shift[i]) * scale[i]
            newTri.append(tuple(newpt))
        if len(removeDupsFromPointList(newTri)) == 3:
            yield newTri
        else:
            pass


def swapAxisWithZ(mesh, axis):
    # axis = 0 -> swap x - z
    # axis = 1 -> swap y - z
    for tri in mesh:
        newTri = []
        for pt in tri:
            newpt = [0, 0, 0]
            if axis == 0:
                newpt[0] = pt[2]
                newpt[1] = pt[1]
                newpt[2] = pt[0]
            elif axis == 1:
                newpt[0] = pt[0]
                newpt[1] = pt[2]
                newpt[2] = pt[1]
            else:
                newpt = pt
            newTri.append(tuple(newpt))
        if len(removeDupsFromPointList(newTri)) == 3:
            yield newTri
        else:
            pass

def swapAxis(mesh, axis1, axis2):
    for tri in mesh:
        newTri = []
        for pt in tri:
            newpt = [pt[0], pt[1], pt[2]]
            newpt[axis1] = pt[axis2]
            newpt[axis2] = pt[axis1]
            newTri.append(tuple(newpt))
        if len(removeDupsFromPointList(newTri)) == 3:
            yield newTri
        else:
            pass