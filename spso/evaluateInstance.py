#!/usr/bin/env python3
import decodeInstance
import initializePopulation
import constraintBasedDecoder
import plotTools
import numpy as np
import matplotlib.pyplot as plt
from random import randint
import update
import evaluate
import time
import math
import copy
import csv



def run(filename, numberOfCustomers, stoppingGap, useLocalSearch, plotBest, saveGlobalUpdates):
    t0 = time.time()
    
    # Parameters
    inertiaWeight = 0.9
    inertiaWeightLowerBound = 0.4
    linearDecreaseStep = 0.000225
    c = 2.0
    phi = 0.3
    numberOfParticles = 20
    depotIndex = 0
    pTournament = 1
    tournamentSize = 2
    pC = [0.05 + 0.45 * (math.exp(10*(i-1)/(numberOfParticles-1))-1) / (math.exp(10)-1) for i in range(1,numberOfParticles+1)]
    iterationsWithoutImprovement = 0
    iterationsWithoutRefreshing = 8
    refreshingGap = 7
    numberOfGenerations = 0
    oldGlobalBestCost = np.Inf
    timeOut = 10000
    newBestNV = []
    newBestTD = []
    genIndex = []

    globalUpdateCost = []
    globalUpdateNV = []
    globalUpdateTD = []
    globalUpdateTime = []

    # Algorithm mode
    mode = 'parallel' # 'parallel' or 'sequential'
    selectionMode = 'random' # 'heuristic' or 'random'
    objectiveFunction = 'minNV'

    # Read and decode instance
    instanceName, vehicleNumber, capacity, customerParams = decodeInstance.readSolomonInstance(filename)
    xCoord, yCoord, demand, readyTime, dueDate, service = decodeInstance.extractParams(customerParams, numberOfCustomers)

    # Init prints
    print('----------------------')
    print('Population size: ', numberOfParticles)
    print('Solomon instance: ', instanceName)
    print('Number of customers (including depot): ', numberOfCustomers)


    # Initialize population
    distances = initializePopulation.getDistances(xCoord, yCoord)
    positions = initializePopulation.initPositions(numberOfParticles, numberOfCustomers, phi, distances, demand, readyTime, dueDate, service, depotIndex)
    velocities = initializePopulation.initVelocities(numberOfParticles, numberOfCustomers)
    particleToFollow = [np.zeros(numberOfCustomers).tolist() for x in range(numberOfParticles)]

    # Constraint-based decoder
    for i, particle in enumerate(positions):
        positions[i] = constraintBasedDecoder.decodeBasedOnConstraints(particle, distances, capacity, demand, readyTime, dueDate, service, mode)

    # Best positions
    bestPositions = [[positions[x][y][:] for y in range(numberOfCustomers)] for x in range(numberOfParticles)]
    bestCosts = [np.Inf for x in range(numberOfParticles)]

    arrivalTimeList = [[] for x in range(numberOfParticles)]
    capacityList = [[] for x in range(numberOfParticles)]
    iterationsWithoutRefreshing = [refreshingGap+1 for i in range(numberOfParticles)]
    while iterationsWithoutImprovement < stoppingGap:

        # Refresh particle to follow
        for i in range(numberOfParticles):
            if iterationsWithoutRefreshing[i] >= refreshingGap:
                for d in range(numberOfCustomers):
                    particleToFollow[i][d] = update.getIndexToLearnFrom(i, bestCosts, pTournament, tournamentSize, pC[i])
                iterationsWithoutRefreshing[i] = 0


        # Update velocity
        for i,iVelocity in enumerate(velocities):
            for d,idVelocity in enumerate(iVelocity):
                velocities[i][d] = update.velocityUpdate(idVelocity, inertiaWeight, c, positions[i][d], bestPositions, particleToFollow[i][d], d)
          
        # Update position
        for i, iPosition in enumerate(positions):
            positions[i], arrivalTimeList[i], capacityList[i] = update.positionUpdate(iPosition, velocities[i], distances, capacity, demand, 
                readyTime, dueDate, service, mode, selectionMode)
            # Local search
            positions[i], capacityList[i], arrivalTimeList[i] = update.localSearch(positions[i], depotIndex, arrivalTimeList[i], capacityList[i], 
                readyTime, service, distances, demand, capacity, dueDate)

        # Evaluate feasible particles
        bestPositions, bestCosts, iterationsWithoutRefreshing = evaluate.updateBest(positions, bestPositions, bestCosts, iterationsWithoutRefreshing, distances, depotIndex, objectiveFunction)


        # Get global best
        globalBestCost = min(bestCosts)
        if globalBestCost < oldGlobalBestCost:
            oldGlobalBestCost = globalBestCost
            iterationsWithoutImprovement = 0

            bestIndex = np.argmin(bestCosts)
            numberOfVehicles, totalDistance = evaluate.evaluateParticle(bestPositions[bestIndex], distances, depotIndex)
            newBestNV.append(numberOfVehicles)
            newBestTD.append(totalDistance)
            genIndex.append(numberOfGenerations)

            if saveGlobalUpdates:
                time_stamp = time.time()
                globalUpdateCost.append(globalBestCost)
                globalUpdateNV.append(numberOfVehicles)
                globalUpdateTD.append(totalDistance)
                globalUpdateTime.append(time_stamp - t0)

        else:
            iterationsWithoutImprovement += 1

        # Linearly decrease inertia weight
        inertiaWeight -= linearDecreaseStep
        #inertiaWeight = beta * inertiaWeight
        
        if inertiaWeight < inertiaWeightLowerBound:
            inertiaWeight = inertiaWeightLowerBound

        numberOfGenerations += 1

        if numberOfGenerations > timeOut:
            break

    bestIndex = np.argmin(bestCosts)
    numberOfVehicles, totalDistance = evaluate.evaluateParticle(bestPositions[bestIndex], distances, depotIndex)

    if plotBest:
        # Plot
        plotTools.plotCustomers(xCoord, yCoord, demand, readyTime, dueDate, service, arrivalTimeList[bestIndex], capacityList[bestIndex], numberOfCustomers)
        plotTools.plotParticle(xCoord, yCoord, bestPositions[bestIndex], color='#%06X' % randint(0, 0xFFFFFF))
        plt.show(block=True)

    t1 = time.time()
    totalRunningTime = t1-t0

    if saveGlobalUpdates:

        savePath = filename[20:-4] + '_' +str(numberOfCustomers-1) + '.csv'
        with open(savePath, mode='w') as csvFile:
                filewriter = csv.writer(csvFile, delimiter=',',
                                quotechar='|', quoting=csv.QUOTE_MINIMAL)
                
                filewriter.writerow(['Global cost', 'Global NV', 'Global TD', 'Global time'])
        csvFile.close()


        for iUpdate in range(len(globalUpdateCost)):
            with open(savePath, 'a') as csvFile:
                filewriter = csv.writer(csvFile, delimiter=',',
                                quotechar='|', quoting=csv.QUOTE_MINIMAL)
                
                filewriter.writerow([globalUpdateCost[iUpdate], globalUpdateNV[iUpdate], globalUpdateTD[iUpdate], globalUpdateTime[iUpdate]])
            csvFile.close()

    print("Elapsed time: ", totalRunningTime)
    print('Number of vehicles: ', numberOfVehicles)
    print('Distance traveled: ', totalDistance)
    print('Number of generations: ', numberOfGenerations)

    return numberOfVehicles, totalDistance, bestIndex, numberOfGenerations, totalRunningTime


if __name__ == "__main__":
    filename = './solomon-instances/r101.txt'
    numberOfCustomers = 12
    stoppingGap = 1000
    useLocalSearch = True
    plotBest = True
    saveGlobalUpdates = True
    run(filename, numberOfCustomers, stoppingGap, useLocalSearch, plotBest, saveGlobalUpdates)