#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import os.path
import csv

def readSolomonInstance(filename):

    distances = []
    f = open(filename, "r")

    if f.mode == 'r':
        lines = f.readlines()

        instanceName = lines[0]
        settings = lines[4].split()
        vehicleNumber = int(settings[0])
        capacity = int(settings[1])
        customerParams = [[int(x) for x in line.split()] for line in lines[9:]] # List of customer lists
    

    return instanceName, vehicleNumber, capacity, customerParams

def extractParams(customerParams, numberOfCustomers):

    xCoord = []
    yCoord = []
    demand = []
    readyTime = []
    dueDate = []
    service = []

    for customer in customerParams[0:numberOfCustomers]:
        xCoord.append(customer[1])
        yCoord.append(customer[2])
        demand.append(customer[3])
        readyTime.append(customer[4])
        dueDate.append(customer[5])
        service.append(customer[6])

    return xCoord, yCoord, demand, readyTime, dueDate, service

if __name__ == "__main__":
    filename = './solomon-instances/r101.txt'
    numberOfCustomers = 7
    instanceName, vehicleNumber, capacity, customerParams = readSolomonInstance(filename)
    xCoord, yCoord, demand, readyTime, dueDate, service = extractParams(customerParams, numberOfCustomers)