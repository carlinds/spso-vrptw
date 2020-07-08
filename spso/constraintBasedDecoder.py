#!/usr/bin/env python3
import numpy as np


def decodeBasedOnConstraints(particle, distances, capacity, demand, readyTime, dueDate, service, mode):

    arrivalTimeCurrentCustomer = 0
    depotIndex = 0
    capacitySum = demand[depotIndex]

    currentCustomer = particle[depotIndex][1][0]
    nextCustomer = particle[depotIndex][1][1]


    while nextCustomer != depotIndex:

        constraintViolated = False

        capacitySum += demand[nextCustomer]
        waitTime = max(readyTime[currentCustomer] - arrivalTimeCurrentCustomer, 0)
        arrivalTimeNextCustomer = arrivalTimeCurrentCustomer + waitTime + service[currentCustomer] + distances[currentCustomer][nextCustomer]


        # Demand constraint
        if capacitySum > capacity:
            constraintViolated = True

        if arrivalTimeNextCustomer > dueDate[nextCustomer]:
            constraintViolated = True

        # Add new vehicle route
        if constraintViolated:
            arcToDepot = (currentCustomer, depotIndex)
            newArc = (depotIndex, nextCustomer)

            particle[depotIndex].append(arcToDepot)
            particle[depotIndex].append(newArc)
            particle[currentCustomer][1] = arcToDepot
            particle[nextCustomer][0] = newArc

            capacitySum = demand[nextCustomer]

            if mode == 'sequential':
                arrivalTimeCurrentCustomer = arrivalTimeNextCustomer
            elif mode == 'parallel':
                arrivalTimeCurrentCustomer = distances[depotIndex][nextCustomer]

        else:
            arrivalTimeCurrentCustomer = arrivalTimeNextCustomer
            
        currentCustomer = nextCustomer
        nextCustomer = particle[currentCustomer][1][1]

    return particle