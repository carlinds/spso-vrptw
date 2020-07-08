#!/usr/bin/env python3
import math
import numpy as np

def evaluateParticle(particlePosition, distances, depotIndex):

    numberOfVehicles = 0
    totalDistance = 0

    for arc in particlePosition[depotIndex]:
        if arc[0] == depotIndex:
            numberOfVehicles += 1
            currentNode = arc[0]
            nextNode = arc[1]

            routeDistance = 0
            while nextNode != depotIndex:
                routeDistance += distances[currentNode][nextNode]
                currentNode = nextNode
                nextNode = particlePosition[currentNode][1][1]
            routeDistance += distances[currentNode][nextNode]
            totalDistance += routeDistance


    return numberOfVehicles, totalDistance


def objective(numberOfVehicles, totalDistance, mode='minNV'):

    if mode == 'minNV':
        cost = numberOfVehicles + np.arctan(totalDistance)/(math.pi/2)
    elif mode == 'minDist':
        cost = totalDistance

    return cost
    

def updateBest(positions, bestPositions, bestCosts, iterationsWithoutRefreshing, distances, depotIndex, mode='minNV'):

    for i, position in enumerate(positions):
        numberOfVehicles, totalDistance = evaluateParticle(position, distances, depotIndex)
        cost = objective(numberOfVehicles, totalDistance, mode)

        if cost < bestCosts[i]:
            bestPositions[i] = [position[a][:] for a in range(len(position))]
            bestCosts[i] = cost
            iterationsWithoutRefreshing[i] = 0
        else:
            iterationsWithoutRefreshing[i] += 1

    return bestPositions, bestCosts, iterationsWithoutRefreshing