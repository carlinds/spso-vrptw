#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt

def plotCustomers(xCoord, yCoord, demand, readyTime, dueDate, service, arrivalTimeListBest, capacityListBest, numberOfCustomers):
    fig,ax = plt.subplots()
    sc = plt.scatter(xCoord, yCoord)
    annot = ax.annotate("", xy=(0,0), xytext=(20,20),textcoords="offset points", bbox=dict(boxstyle="round", fc="w"))
    annot.set_visible(False)


    def update_annot(ind):
        pos = sc.get_offsets()[ind["ind"][0]]
        annot.xy = pos
        customerIndex = ind["ind"][0]
        text = ('Customer '+ str(customerIndex) + '\nDemand: ' + str(demand[customerIndex]) 
            + '\nReady time: ' + str(readyTime[customerIndex]) + '\nDue date: ' + str(dueDate[customerIndex]) 
            + '\nService: ' + str(service[customerIndex]) + '\nArrival time: ' + str(arrivalTimeListBest[customerIndex])
            + '\nVehicle capacity: ' + str(capacityListBest[customerIndex]))
        annot.set_text(text)
        annot.get_bbox_patch().set_alpha(0.4)

    def hover(event):
        vis = annot.get_visible()
        if event.inaxes == ax:
            cont, ind = sc.contains(event)
            if cont:
                update_annot(ind)
                annot.set_visible(True)
                fig.canvas.draw_idle()
            else:
                if vis:
                    annot.set_visible(False)
                    fig.canvas.draw_idle()

    fig.canvas.mpl_connect("motion_notify_event", hover)

def plotParticle(xCoord, yCoord, particlePosition, color = 'b', linestyle = '-'):
    for dimension in particlePosition:
        for arc in dimension:
            x1 = xCoord[arc[0]]
            y1 = yCoord[arc[0]]
            x2 = xCoord[arc[1]]
            y2 = yCoord[arc[1]]
            plt.arrow(x1, y1, x2-x1, y2-y1, color=color, linestyle = linestyle, 
            head_width = 0.7, length_includes_head = True, alpha = 0.5, shape = 'full')
