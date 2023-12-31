import math
from collections import defaultdict
import numpy as np

def linesToVoxels(lineList, pixels, isShell):
    if isShell:
        for x in range(len(pixels)):
            #isBlack = False
            lines = list(findRelevantLines(lineList, x))
            targetYs = list(map(lambda line:int(generateY(line,x)),lines))
            for y in range(len(pixels[x])):
                #if isBlack:
                #    pixels[x][y] = True
                if y in targetYs:
                    for line in lines:
                        if onLine(line, x, y):
                            #isBlack = not isBlack
                            pixels[x][y] = True
    else:
        for x in range(len(pixels)):
            isBlack = False
            lines = list(findRelevantLines(lineList, x))
            targetYs = list(map(lambda line:int(generateY(line,x)),lines))
            for y in range(len(pixels[x])):
                if isBlack:
                    pixels[x][y] = True
                if y in targetYs:
                    for line in lines:
                        if onLine(line, x, y):
                            isBlack = not isBlack
                            pixels[x][y] = True

            if isBlack:
                print("An error has occured at x%sz%s - is the geometry watertight?"%(x,lineList[0][0][2]))


# Voxelize solid, watertight body
def linesToVoxelsSolid(lineList, pixels):
    for x in range(len(pixels)):
        isBlack = False
        lines = list(findRelevantLines(lineList, x))
        targetYs = list(map(lambda line:int(generateY(line,x)),lines))
        for y in range(len(pixels[x])):
            if isBlack:
                pixels[x][y] = True
            if y in targetYs:
                for line in lines:
                    if onLine(line, x, y):
                        isBlack = not isBlack
                        pixels[x][y] = True

        if isBlack:
            print("An error has occured at x%sz%s - is the geometry watertight?"%(x,lineList[0][0][2]))


def linesToVoxelsShell(lineList, pixels):
    for x in range(len(pixels)):
        #isBlack = False
        lines = list(findRelevantLines(lineList, x))
        targetYs = list(map(lambda line:int(generateY(line,x)),lines))
        for y in range(len(pixels[x])):
            #if isBlack:
            #    pixels[x][y] = True
            if y in targetYs:
                for line in lines:
                    if onLine(line, x, y):
                        #isBlack = not isBlack
                        pixels[x][y] = True

# This only voxelises the outline. However, it is not watertight!
def linesToVoxelsShell2(lineList, pixels):
    for line in lineList:
        x1, y1 = line[0][0], line[0][1]
        x2, y2 = line[1][0], line[1][1]
        pixels[int(x1)][int(y1)] = True
        pixels[int(x2)][int(y2)] = True


def findRelevantLines(lineList, x, ind=0):
    for line in lineList:
        same = False
        above = False
        below = False
        for pt in line:
            if pt[ind] > x:
                above = True
            elif pt[ind] == x:
                same = True
            else:
                below = True
        if above and below:
            yield line
        elif same and above:
            yield line


def generateY(line, x):
    if line[1][0] == line[0][0]:
        return -1
    ratio = (x - line[0][0]) / (line[1][0] - line[0][0])
    ydist = line[1][1] - line[0][1]
    newy = line[0][1] + ratio * ydist
    return newy


def onLine(line, x, y):
    newy = generateY(line, x)
    if int(newy) != y:
        return False
    if int(line[0][0]) != x and int(line[1][0]) != x and (max(line[0][0], line[1][0]) < x or min(line[0][0], line[1][0]) > x):
        return False
    if int(line[0][1]) != y and int(line[1][1]) != y and (max(line[0][1], line[1][1]) < y or min(line[0][1], line[1][1]) > y):
        return False
    return True

