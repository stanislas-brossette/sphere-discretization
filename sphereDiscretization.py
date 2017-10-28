#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt
import argparse

#Default values
default_height = 0.5
default_nIter = 3
default_radius = 1.0
default_center = [0.0, 0.0, 0.0]

#Argument parsing
parser = argparse.ArgumentParser()
parser.add_argument("--height", help="Height of the half-space",
        type=float, default=default_height)
parser.add_argument("--nIter", help="Number of iteration of the spheres discretization",
        type=int, default=default_nIter)
parser.add_argument("--radius", help="Radius of the sphere",
        type=float, default=default_radius)
parser.add_argument("--center", help="Center of the sphere", nargs=3,
        type=float, default=default_center)
args = parser.parse_args()

def main(args):
    center = np.array(args.center)
    radius = args.radius
    nIter = args.nIter
    faces = discretizeSphere(center, radius, nIter)
    faces = cutSphere(faces, args.height)
    plotFaces(faces)
    return

def discretizeSphere(center, radius, nIter):
    pXp = center + radius * np.array([1.0, 0.0, 0.0])
    pXm = center + radius * np.array([-1.0, 0.0, 0.0])
    pYp = center + radius * np.array([0.0, 1.0, 0.0])
    pYm = center + radius * np.array([0.0, -1.0, 0.0])
    pZp = center + radius * np.array([0.0, 0.0, 1.0])
    pZm = center + radius * np.array([0.0, 0.0, -1.0])
    faces = [
            [pXp, pYp, pZp],
            [pXp, pYm, pZp],
            [pXp, pYp, pZm],
            [pXp, pYm, pZm],
            [pXm, pYp, pZp],
            [pXm, pYm, pZp],
            [pXm, pYp, pZm],
            [pXm, pYm, pZm]
            ]
    iter = 0
    while iter < nIter:
        newFaces = []
        for face in faces:
            m01 = (face[0] + face[1])/2.0
            m12 = (face[1] + face[2])/2.0
            m20 = (face[2] + face[0])/2.0
            m01 = center + radius * (m01-center)/np.linalg.norm(m01-center)
            m12 = center + radius * (m12-center)/np.linalg.norm(m12-center)
            m20 = center + radius * (m20-center)/np.linalg.norm(m20-center)
            newFaces.append([m01, m20, face[0]])
            newFaces.append([m12, m01, face[1]])
            newFaces.append([m20, m12, face[2]])
            newFaces.append([m01, m12, m20])
        faces = newFaces
        iter = iter + 1
    return faces

def cutSphere(faces, height):
    newFaces = []
    for face in faces:
        if (face[0][2]>height and face[1][2]>height and face[2][2]>height):
            newFaces.append(face)
        elif (face[0][2]>height or face[1][2]>height or face[2][2]>height):
            newFaces.append(computeCutFace(face, height))
    return newFaces

def cutSegment(p0, p1, height):
    if(p0[2]>height and p1[2]>height):
        return [p0, p1]
    elif(p0[2]>height or p1[2]>height):
        pUp = []
        pDown = []
        if(p0[2]>height):
            pUp = p0
            pDown = p1
        if(p1[2]>height):
            pUp = p1
            pDown = p0
        alpha = (height - pDown[2])/(pUp[2] - pDown[2])
        pMid = pDown + alpha*(pUp - pDown)
        return [pMid, pUp]

def computeCutFace(face, height):
    newFace = []
    nPointsAbove = 0
    pointsAbove = []
    pointsBelow = []
    if(face[0][2]>height):
        nPointsAbove = nPointsAbove + 1
        pointsAbove.append(face[0])
    else:
        pointsBelow.append(face[0])
    if(face[1][2]>height):
        nPointsAbove = nPointsAbove + 1
        pointsAbove.append(face[1])
    else:
        pointsBelow.append(face[1])
    if(face[2][2]>height):
        nPointsAbove = nPointsAbove + 1
        pointsAbove.append(face[2])
    else:
        pointsBelow.append(face[2])

    if(nPointsAbove == 3):
        return face
    elif(nPointsAbove == 1):
        seg0 = cutSegment(pointsAbove[0],pointsBelow[0],height)
        seg1 = cutSegment(pointsAbove[0],pointsBelow[1],height)
        newFace.append(pointsAbove[0])
        newFace.append(seg0[0])
        newFace.append(seg1[0])
    elif(nPointsAbove == 2):
        seg0 = cutSegment(pointsAbove[0],pointsBelow[0],height)
        seg1 = cutSegment(pointsAbove[1],pointsBelow[0],height)
        newFace.append(pointsAbove[0])
        newFace.append(pointsAbove[1])
        newFace.append(seg1[0])
        newFace.append(seg0[0])
    return newFace

def plotFaces(faces):
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    for face in faces:
        lineX = np.array([])
        lineY = np.array([])
        lineZ = np.array([])
        for vertex in face:
            lineX = np.append(lineX, vertex[0])
            lineY = np.append(lineY, vertex[1])
            lineZ = np.append(lineZ, vertex[2])
        lineX = np.append(lineX, face[0][0])
        lineY = np.append(lineY, face[0][1])
        lineZ = np.append(lineZ, face[0][2])
        ax.plot(lineX, lineY, lineZ)
    ax.set_xlim(-1.0, 1.0)
    ax.set_ylim(-1.0, 1.0)
    ax.set_zlim(-1.0, 1.0)
    plt.show()


main(args)
