from __future__ import division, print_function
import numpy as np
from glob import glob
import os
from AIA_tools import *

pathToData ="/home/john/gitRepos/REU/jwaczak/data/correctedSyntheticObservations/"
pathToObs = "/home/john/gitRepos/REU/jwaczak/data/observationData/obs_data.txt"

dataFiles = os.listdir(pathToData)

obs = np.loadtxt(pathToObs, delimiter=',')

data = []

def chiSquared(obs, sim):
    X_sq = 0

    # NOTE for some reason we have an extra data point in the sims than in the obs data
    if len(sim) == len(obs)+1:
        sim = sim[:-1]
    assert(len(obs) == len(sim))
    for i in range(int(len(obs)/2)):
        dev = obs[i]-sim[i]
        X_sq += (dev)**2/obs[i]  # sqrt(obs[i]) measures the error in the observational data

    return X_sq



for file_ in dataFiles:
    freeParams = 2
    d = np.loadtxt(pathToData+file_, delimiter=',')

    sum = 0
    for i in range(1,6):
        X = chiSquared(obs[:,i], d[:,i])
        sum += X

    sum = sum/(5*(len(d[:,0])/2)-freeParams)
    data.append(sum)

data = np.array(data)
try:
    min_index = analysis.getNearestValue(data, 1.0)
    print("SUM:  ", dataFiles[min_index], "\tX2:", data[min_index])


except:
    print("Could not find minimum value(s)")

