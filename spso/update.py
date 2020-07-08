#!/usr/bin/env python3
import numpy as np
from initializePopulation import velocity


def coeffTimesVelocity(coefficient, idVelocity):   
    for i in range(len(idVelocity)):
        possibility = idVelocity[i].possibility
        candidatePossibility = coefficient*possibility

        if candidatePossibility > 1:
            updatedPossibility = 1
        else:
            updatedPossibility = candidatePossibility

        idVelocity[i].possibility = updatedPossibility

    return idVelocity


def addVelocities(idVelocity, jdVelocity):

    # Add all velocities
    velocitySum = idVelocity + jdVelocity

    # Loop over both velocities
    for i, iObject in enumerate(idVelocity):
        for jObject in jdVelocity:
            
            # Set possibility to max and remove duplicate
            if iObject.arc == jObject.arc:
                pArc = max(iObject.possibility, jObject.possibility)
                velocitySum[i].possibility = pArc
                velocitySum.remove(jObject)

    return velocitySum


def subtractPositions(idPosition, jdPosition):
    differencePosition = idPosition[:] 
    for iArc in idPosition:
        for jArc in jdPosition:
            if iArc == jArc:
                differencePosition.remove(iArc)
                break

    return differencePosition


def coeffTimesPosition(coefficient, dPosition):
    particleVelocity = []

    assert(coefficient >= 0)

    for arc in dPosition:
        if coefficient > 1:
            velocityObject = velocity(arc, 1)
            particleVelocity.append(velocityObject)        
        else:
            velocityObject = velocity(arc, coefficient)
            particleVelocity.append(velocityObject)
    return particleVelocity


def tournamentSelection(particleIndexBeingUpdated, bestCosts, pTournament, tournamentSize):

    numberOfParticles = len(bestCosts)
    possiblePositions = [i for i in range(numberOfParticles)]
    possiblePositions.remove(particleIndexBeingUpdated)
    candidatePositions = np.random.choice(possiblePositions, tournamentSize, False).tolist()

    individualSelected = False

    while not individualSelected:
        bestIndividualIndex = np.argmin([bestCosts[x] for x in candidatePositions])

        r = np.random.rand()

        if len(candidatePositions) > 1:
            if (r < pTournament):
                iSelected = candidatePositions[bestIndividualIndex]
                individualSelected = True
            else:
                candidatePositions.pop(bestIndividualIndex)
        else:
            iSelected = candidatePositions[0]
            individualSelected = True

    return iSelected

def tryInsertionNewCustomer(newCustomer, route, newArcs, oldArcs, vehicleRoutes, endTimeList, capacityList, arrivalTimeList, distance, readyTime, service, demand, capacity, dueDate):
    for j, vehicleRoute in enumerate(vehicleRoutes):
            for i, currentCustomer in enumerate(vehicleRoute[:-1]):
                nextCustomer = vehicleRoute[i+1]
                
                arrivalTimeNewCustomer = endTimeList[currentCustomer] + distance[currentCustomer][newCustomer]
                waitTimeNewCustomer = max(readyTime[newCustomer] - arrivalTimeNewCustomer, 0)
                
                arrivalTimeAtNextCustomer = (arrivalTimeNewCustomer + waitTimeNewCustomer + service[newCustomer] +
                                    distance[newCustomer][nextCustomer])
                waitTimeNextCustomer = max(readyTime[nextCustomer] - arrivalTimeAtNextCustomer, 0)
                newEndTimeNextCustomer = arrivalTimeAtNextCustomer + waitTimeNextCustomer + service[nextCustomer]

                newCapacitySum = capacityList[vehicleRoute[-2]] + demand[newCustomer]

                # Constraints
                arrivalTimeFeasibleNewCustomer = (arrivalTimeNewCustomer < dueDate[newCustomer])
                endTimeNotChangedNextCustomer = (abs(newEndTimeNextCustomer - endTimeList[nextCustomer]) < 1e-5)
                capacityFeasible = (newCapacitySum <= capacity)



                if endTimeNotChangedNextCustomer and capacityFeasible and arrivalTimeFeasibleNewCustomer:
                    newArcs.append((currentCustomer, newCustomer))
                    newArcs.append((newCustomer, nextCustomer))
                    oldArcs.append((currentCustomer, nextCustomer))

                    # Remove the arcs that previously went to and from the inserted customer
                    routeIndex = route.index(newCustomer)
                    oldArcs.append((route[routeIndex-1], newCustomer))
                    oldArcs.append((newCustomer, route[routeIndex+1]))

                    # Update vehicle routes and vehicle capacity at each customer
                    capacityList[newCustomer] = capacityList[currentCustomer] + demand[newCustomer]
                    k = i + 1
                    while vehicleRoute[k]!=0:
                        capacityList[k] += demand[newCustomer]
                        k += 1
                    vehicleRoutes[j] = vehicleRoute[:i+1] + [newCustomer] + vehicleRoute[i+1:]

                    # Update end time for new customer
                    endTimeList[newCustomer] = arrivalTimeNewCustomer + waitTimeNewCustomer + service[newCustomer]
                    arrivalTimeList[newCustomer] = arrivalTimeNewCustomer
                    arrivalTimeList[nextCustomer] = arrivalTimeAtNextCustomer

                    route.remove(newCustomer)
                    return newArcs, oldArcs, route, vehicleRoutes, capacityList, endTimeList, arrivalTimeList

    return newArcs, oldArcs ,route, vehicleRoutes, capacityList, endTimeList, arrivalTimeList

def localSearch(particlePosition, depotIndex, arrivalTimeList, capacityList, readyTime, service, distance, demand, capacity, dueDate):

    # Extract all vehicle routes from particle
    vehicleRoutes = []
    routeIndex = -1
    leastNumberOfCustomers = np.inf

    for arc in particlePosition[depotIndex]:
        if arc[0] == depotIndex:
            routeIndex += 1
            vehicleRoutes.append([])
            currentNode = arc[0]
            nextNode = arc[1]

            vehicleRoutes[routeIndex].append(depotIndex)
            vehicleRoutes[routeIndex].append(nextNode)
            numberOfCustomers = 0

        while nextNode != depotIndex:
            currentNode = nextNode
            nextNode = particlePosition[currentNode][1][1]
            vehicleRoutes[routeIndex].append(nextNode)
            numberOfCustomers += 1

        if numberOfCustomers < leastNumberOfCustomers:
            leastNumberOfCustomers = numberOfCustomers
            leastRouteIndex = routeIndex

    # Choose the vehicle route with fewest customers
    route = vehicleRoutes[leastRouteIndex]
    vehicleRoutes.remove(route)

    endTimeList = [0 for i in range(len(arrivalTimeList))]
    # Compute end time for all customers
    for k in range(len(arrivalTimeList)):
        waitTime = max(readyTime[k] - arrivalTimeList[k], 0)
        endTimeList[k] = arrivalTimeList[k] + waitTime + service[k]

    counter = 0
    newArcs = []
    oldArcs = []

    # Try to insert customers into other routes
    for newCustomer in route[1:-1]:
        counter += 1
        newArcs, oldArcs, route, vehicleRoutes, capacityList, endTimeList, arrivalTimeList = tryInsertionNewCustomer(
            newCustomer, route, newArcs, oldArcs, vehicleRoutes, endTimeList, capacityList, arrivalTimeList,
            distance, readyTime, service, demand, capacity, dueDate)

        # If the current customer can't be inserted into any other routes, break the loop
        if len(route[1:-1]) != leastNumberOfCustomers - counter:
            break        

    # If the insertions of all customers were successful then construct the new particlePosition
    if route[1:-1] == []:
        newParticlePosition = np.empty(len(arrivalTimeList), dtype=object)
        newParticlePosition[depotIndex] = []
        for vechicleRoute in vehicleRoutes:
            for i, currentCustomer in enumerate(vechicleRoute[:-1]):

                if currentCustomer != depotIndex:
                    nextCustomer = vechicleRoute[i+1]
                    previousCustomer = vechicleRoute[i-1]

                    newArcFrom = (previousCustomer, currentCustomer)
                    newArcTo = (currentCustomer, nextCustomer)

                    newParticlePosition[currentCustomer] = [newArcFrom, newArcTo]
                
            # Each vehicle route only has one arc going from the depotIndex and one arc going to the depotIndex
            arcFromDepot = (depotIndex, vechicleRoute[1])
            arcToDepot = (vechicleRoute[-2], depotIndex)
            newParticlePosition[depotIndex].append(arcToDepot)
            newParticlePosition[depotIndex].append(arcFromDepot)

        return newParticlePosition, capacityList, arrivalTimeList

    return particlePosition, capacityList, arrivalTimeList


def getIndexToLearnFrom(particleIndexBeingUpdated, bestCosts, pTournament, tournamentSize, pC):

    assert(0 <= pC <= 1)

    r = np.random.rand()
    if r > pC:
        indexToFollow = particleIndexBeingUpdated
    else:
        indexToFollow = tournamentSelection(particleIndexBeingUpdated, bestCosts, pTournament, tournamentSize)

    return indexToFollow

def velocityUpdate(idVelocity, inertiaWeight, c, idPosition, bestPositions, particleToFollow, d):
    """
        CLPSO implementation
    """
    inertiaTimesidVelocity = coeffTimesVelocity(inertiaWeight, idVelocity)
    pBest = bestPositions[particleToFollow][d]
    pBestCopy = pBest[:]
    differencePosition = subtractPositions(pBestCopy, idPosition)
    
    # Draw random number
    r = np.random.rand()
    cRandTimesPositionDiff = coeffTimesPosition(c*r, differencePosition)
    idVelocityUpdated = addVelocities(inertiaTimesidVelocity, cRandTimesPositionDiff)

    return idVelocityUpdated


def cut(vij, alpha):
    cutVij = []
    for e in vij:
        if e.possibility >= alpha:
            cutVij.append(e.arc)

    return cutVij

def createFeasibleSet(setOfArcs, currentCustomer, arrivalTime, capacitySum, distances, capacity, demand, readyTime, dueDate, service, unvisitedNodes):
    
    feasibleArcs = []    

    for arc in setOfArcs:
        nextCustomer = arc[1]
        feasible = True

        if (arc[0] != currentCustomer) or (nextCustomer == arc[0]) or (nextCustomer not in unvisitedNodes):
            feasible = False
        else:
            capacitySum += demand[nextCustomer]
            waitTime = max(readyTime[currentCustomer] - arrivalTime, 0)
            arrivalTimeNextCustomer = arrivalTime + waitTime + service[currentCustomer] + distances[currentCustomer][nextCustomer]
            
            # Demand constraint
            if capacitySum > capacity:
                feasible = False
            
            # Time window constraint
            if arrivalTimeNextCustomer > dueDate[nextCustomer]: # or arrivalTimeNextCustomer < readyTime[nextCustomer]
                feasible = False
        
        if feasible:
            feasibleArcs.append(arc)

    return feasibleArcs
    
def updateFeasibleSets(cutVelocity, particlePosition, k, arrivalTimeCurrentCustomer, capacitySum, distances, capacity, demand, readyTime, dueDate, service, unvisitedNodes):
    numberOfCustomers = len(demand)
    feasibleSetV = createFeasibleSet(cutVelocity[k], k, arrivalTimeCurrentCustomer, capacitySum, distances, capacity, demand, readyTime, dueDate, service, unvisitedNodes)
    feasibleSetX = createFeasibleSet(particlePosition[k], k, arrivalTimeCurrentCustomer, capacitySum, distances, capacity, demand, readyTime, dueDate, service, unvisitedNodes)
    setE = [(k, x) for x in range(numberOfCustomers)]
    setE.remove((k,k))
    feasibleSetE = createFeasibleSet(setE, k, arrivalTimeCurrentCustomer, capacitySum, distances, capacity, demand, readyTime, dueDate, service, unvisitedNodes)

    return feasibleSetV, feasibleSetX, feasibleSetE

def select(arcSet, k, arrivalTimeAtCurrentCustomer, distances, readyTime, dueDate, service, selectionMode):
    numberOfArcs = len(arcSet)
    maxVal = np.inf

    assert (selectionMode == 'random' or selectionMode == 'heuristic'), "Unavailable selection mode"
 
    if selectionMode == 'random':
        randomIndex = np.random.randint(0, numberOfArcs)
        m = arcSet[randomIndex][1]

    if selectionMode == 'heuristic':
        heuristic = np.ones(len(dueDate)) * maxVal
        waitTimeAtCurrentCustomer = max(readyTime[k] - arrivalTimeAtCurrentCustomer, 0)
        for arc in arcSet:
            nextNode = arc[1]
            travelTime = distances[k][nextNode]
            arrivalTimeAtNextCustomer = arrivalTimeAtCurrentCustomer + waitTimeAtCurrentCustomer + service[k] + travelTime
            waitTime = max(readyTime[nextNode] - arrivalTimeAtNextCustomer, 0)

            heuristic[nextNode] = (travelTime + waitTime) + (dueDate[nextNode] - arrivalTimeAtNextCustomer)
        m = np.argmin(heuristic)

    return m


def positionUpdate(particlePosition, particleVelocity, distances, capacity, demand, readyTime, dueDate, service, mode, selectionMode):
    
    numberOfCustomers = len(particleVelocity)
    arrivalTimeCurrentCustomer = 0
    arrivalTimeList = [0 for x in range(numberOfCustomers)]
    capacityList = [0 for x in range(numberOfCustomers)]
    depotIndex = 0
    capacitySum = demand[depotIndex]
    unvisitedNodes = [node for node in range(numberOfCustomers)]
    unvisitedNodes.remove(depotIndex)
    alpha = np.random.rand()

    cutVelocity = np.empty(numberOfCustomers, dtype=object)
    for j, vij in enumerate(particleVelocity):
        cutVelocity[j] = cut(vij, alpha)
    
    newXi = np.empty(numberOfCustomers, dtype=object)
    for x in range(numberOfCustomers):
        newXi[x] = []

    k = depotIndex
    t = depotIndex
    m = -1

    feasibleSetV, feasibleSetX, feasibleSetE = updateFeasibleSets(cutVelocity, particlePosition, k, arrivalTimeCurrentCustomer, capacitySum, distances, capacity, demand, readyTime, dueDate, service, unvisitedNodes)

    while not((k == depotIndex) and (unvisitedNodes == [])):
        if feasibleSetV != []:
            m = select(feasibleSetV, k, arrivalTimeCurrentCustomer, distances, readyTime, dueDate, service, selectionMode)
            if m != depotIndex:
                unvisitedNodes.remove(m)
            newArc = (t, m)
            newXi[t].append(newArc)
            newXi[m].append(newArc)
            waitTime = max(readyTime[k] - arrivalTimeCurrentCustomer, 0)
            arrivalTimeCurrentCustomer += waitTime + service[k] + distances[k][m]
            arrivalTimeList[m] = arrivalTimeCurrentCustomer
            capacitySum += demand[m]
            capacityList[m] = capacitySum
            k = m
            t = m
            feasibleSetV, feasibleSetX, feasibleSetE = updateFeasibleSets(cutVelocity, particlePosition, k, arrivalTimeCurrentCustomer, capacitySum, distances, capacity, demand, readyTime, dueDate, service, unvisitedNodes)

        elif feasibleSetX != []:
            m = select(feasibleSetX, k, arrivalTimeCurrentCustomer, distances, readyTime, dueDate, service, selectionMode)
            if m != depotIndex:
                unvisitedNodes.remove(m)
            newArc = (t, m)
            newXi[t].append(newArc)
            newXi[m].append(newArc)
            waitTime = max(readyTime[k] - arrivalTimeCurrentCustomer, 0)
            arrivalTimeCurrentCustomer += waitTime + service[k] + distances[k][m]
            arrivalTimeList[m] = arrivalTimeCurrentCustomer
            capacitySum += demand[m]
            capacityList[m] = capacitySum
            k = m
            t = m
            feasibleSetV, feasibleSetX, feasibleSetE = updateFeasibleSets(cutVelocity, particlePosition, k, arrivalTimeCurrentCustomer, capacitySum, distances, capacity, demand, readyTime, dueDate, service, unvisitedNodes)
        
        elif feasibleSetE != []:
            m = select(feasibleSetE, k, arrivalTimeCurrentCustomer, distances, readyTime, dueDate, service, selectionMode)
            if m != depotIndex:
                unvisitedNodes.remove(m)
            newArc = (t, m)
            newXi[t].append(newArc)
            newXi[m].append(newArc)
            waitTime = max(readyTime[k] - arrivalTimeCurrentCustomer, 0)
            arrivalTimeCurrentCustomer += waitTime + service[k] + distances[k][m]
            arrivalTimeList[m] = arrivalTimeCurrentCustomer
            capacitySum += demand[m]
            capacityList[m] = capacitySum
            k = m
            t = m
            feasibleSetV, feasibleSetX, feasibleSetE = updateFeasibleSets(cutVelocity, particlePosition, k, arrivalTimeCurrentCustomer, capacitySum, distances, capacity, demand, readyTime, dueDate, service, unvisitedNodes)

        else:
            newXi[t].append((t, depotIndex))
            newXi[depotIndex].append((t, depotIndex))

            if mode == 'sequential':
                waitTime = max(readyTime[k] - arrivalTimeCurrentCustomer, 0)
                arrivalTimeCurrentCustomer += waitTime + service[k] + distances[k][depotIndex]
            elif mode == 'parallel':
                arrivalTimeCurrentCustomer = 0

            capacitySum = demand[depotIndex]
            k = depotIndex
            t = depotIndex
            feasibleSetV, feasibleSetX, feasibleSetE = updateFeasibleSets(cutVelocity, particlePosition, k, arrivalTimeCurrentCustomer, capacitySum, distances, capacity, demand, readyTime, dueDate, service, unvisitedNodes)

    return newXi, arrivalTimeList, capacityList