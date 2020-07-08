#!/usr/bin/env python3
import numpy as np
import itertools
import random

class velocity:
    def __init__(self, arc, possibility):
        self.arc = arc
        self.possibility = possibility

    def __repr__(self):
        print_str = str(self.arc) + '/' + str(round(self.possibility, 2))
        return print_str

    def __str__(self):
        return 'Test velocity'

def initPositions(numberOfParticles, numberOfCustomers, phi, distances, demand, readyTime, dueDate, service, depotIndex):
    # Check that phi is within allowed boundaries
    assert(0 <= phi <= 1)

    particles = np.empty(numberOfParticles, dtype=object)

    # Loop over the number of particles and intialize positions. A particle position is either intialized randomly or according
    # to a nearest neighbour heuristic depending on if phi is smaller or larger than a uniformly random number in the range
    # [0, 1]
    for i in range(numberOfParticles):
        r = np.random.rand()
        if r < phi:
            particlePosition = nearestNeighbourHeuristic(distances, demand, readyTime, dueDate, service, depotIndex)
            particles[i] = particlePosition
        else:
            particlePosition = randomInitialization(numberOfCustomers)
            particles[i] = particlePosition

    return particles


def initVelocities(numberOfParticles, numberOfCustomers):
    customerList = [x for x in range(numberOfCustomers)]
    universalSet = list(itertools.permutations(customerList,2))
    velocities = np.empty(numberOfParticles, dtype=object)
    
    # For each particle, sample random arcs and possibilities.
    for i in range(numberOfParticles):
        particleVelocity = [[] for x in range(numberOfCustomers)]
        randomArcs = random.sample(universalSet, numberOfCustomers)
        for arc in randomArcs:
            possibility = np.random.rand()
            velocityObject = velocity(arc, possibility)

            particleVelocity[arc[0]].append(velocityObject)
            particleVelocity[arc[1]].append(velocityObject)

            

        velocities[i] = particleVelocity
 
    return velocities


def randomInitialization(numberOfCustomers):
    # Starting at customer 0 i.e. the depot
    startCustomer = 0
    
    path = [startCustomer] + np.random.permutation(range(1, numberOfCustomers)).tolist()

    # Representation of the position according to the paper
    xi = np.empty(numberOfCustomers, dtype=object)
    for d in range(numberOfCustomers):
        customerIndex = path.index(d)
        if customerIndex == numberOfCustomers-1:
            xid = [(path[customerIndex-1], d), (d, startCustomer)]
        else:
            xid = [(path[customerIndex-1], d), (d, path[customerIndex+1])]
        xi[d] = xid

    particlePosition = xi

    return particlePosition
    

def getDistances(xCoord, yCoord):

    numberOfCustomers = len(xCoord)
    distances = np.zeros((numberOfCustomers, numberOfCustomers))

    # Calculates the Euclidian distance between all nodes
    for i in range(numberOfCustomers):
        for j in range(numberOfCustomers):
            if i == j:
                distances[i][j] = 0
            else:
                iCity = np.array((xCoord[i], yCoord[i]))
                jCity = np.array((xCoord[j], yCoord[j]))
                distances[i][j] = np.linalg.norm(iCity - jCity)

    return distances

def nearestNeighbourHeuristic(distances, demand, readyTime, dueDate, service, depotIndex):
    startCustomer = depotIndex
    maxVal = np.inf
    numberOfCustomers = len(demand)
    # arrivalTime and waitTime at depot is 0
    arrivalTimePreviousCustomer = 0
    waitTimePreviousCustomer = 0
    servicePreviousCustomer = 0
    
    customersLeft = [i for i in range(numberOfCustomers)]
    heuristic = np.zeros(numberOfCustomers)
    
    nhhPath = []
    nhhPath.append(startCustomer)
    heuristic[startCustomer] = maxVal
    currentCustomer = startCustomer
    customersLeft.remove(startCustomer)

    for k in range(numberOfCustomers-1):
        for city in customersLeft:
            travelTime = distances[currentCustomer][city]
            arrivalTime = arrivalTimePreviousCustomer + waitTimePreviousCustomer + servicePreviousCustomer + travelTime
            waitTime = max(readyTime[city] - arrivalTime, 0)

            heuristic[city] = (travelTime + waitTime) + (dueDate[city] - arrivalTime)

        nearestCustomer = np.argmin(heuristic)
        arrivalTimePreviousCustomer = arrivalTimePreviousCustomer + waitTimePreviousCustomer + servicePreviousCustomer + distances[currentCustomer][nearestCustomer]
        nhhPath.append(nearestCustomer)
        currentCustomer = nearestCustomer
        # waitTime and service for the current customer (which is the previous customer in the next iteration)
        waitTimePreviousCustomer = max(readyTime[nearestCustomer] - arrivalTimePreviousCustomer, 0)
        servicePreviousCustomer  = service[nearestCustomer]
        heuristic[nearestCustomer] = maxVal
        customersLeft.remove(nearestCustomer)

    # Representation of the position according to the paper
    xi = np.empty(numberOfCustomers, dtype=object)
    for d in range(numberOfCustomers):
        pathIndex = nhhPath.index(d)

        if pathIndex == numberOfCustomers-1:
            xid = [(nhhPath[pathIndex-1], d), (d, startCustomer)]
        else:
            xid = [(nhhPath[pathIndex-1], d), (d, nhhPath[pathIndex+1])]
            
        xi[d] = xid

    particlePosition = xi

    return particlePosition


    
