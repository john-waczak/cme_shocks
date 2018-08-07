from __future__ import division 
from __future__ import print_function 
import matplotlib.pyplot as plt 
import numpy as np 
import os 
from glob import glob 
from sunpy.net import Fido, attrs
from datetime import datetime, timedelta
import astropy.units as u 
import sunpy.map as smap 
from astropy.coordinates import SkyCoord 
from scipy.interpolate import interp1d
from astropy.convolution import convolve, Box1DKernel
from datetime import datetime, timedelta
import matplotlib.dates as mdates 
import pickle
import AIA_tools as aia
import multiprocessing as mp 


obs_data = np.loadtxt('/data/khnum/REU2018/jwaczak/data/observationData/shockData.txt', delimiter=',')

time_vals = obs_data[:,0]
print(len(time_vals))
Dt = np.mean(time_vals[1:]-time_vals[0:-1])


# trying out temperature from Ma et al. paper 
# te_sta = np.array([1.7*1e6, 1.8*1e6, 1.9*1e6])
# te_end = aia.simulation.getTempVals()
# low_t_index = aia.analysis.getNearestValue(te_end, te_sta.max())
# te_end = te_end[low_t_index:]
# N = [7]
te_sta = [1.8e6, 1.5e6, 1.3e6]
te_end = [2.8e6]
N = [8, 7.5, 7]


pathToChiantiEmiss = '/data/khnum/REU2018/jwaczak/data/chiantiEmissData'

r0 = np.loadtxt('/data/khnum/REU2018/jwaczak/data/initialRadii.txt', delimiter=',')
dx = r0[0]/20

for t0 in te_sta:
    for t1 in te_end:
        for n in N:
            print('t0: {}  t1: {}  n: {}'.format(t0, t1, n))
            fname = 't0--{:.2E}__t1--{:.2E}__n--{:.2E}'.format(t0, t1, 10**n)
            element_list = '2, 6, 7, 8, 10, 12, 13, 14, 16, 18, 20, 26, 28'  # He, C, N, O, Ne, Mg, Al, Si, S, Ar, Ca, Fe, Ni
            simDataFile = aia.simulation.run(t0, t1, 10**n, num=13,
                                             indices=element_list, ntime=len(time_vals), dt= Dt, filename=fname+'.dat')

            obs, simParams = aia.simulation.getSyntheticObservation_II(t0, t1, n, simDataFile) 

            data = [[obs['time'][i], obs['171'][i], obs['193'][i], obs['211'][i], obs['304'][i], obs['335'][i]] for i in range(len(obs['time']))]
            data = np.asarray(data)
            np.savetxt('/data/khnum/REU2018/jwaczak/data/simOutput/'+fname+'.txt', data, delimiter=',', header='t0: {} t1: {} n: {} \n times, 171, 193, 211, 304, 335'.format(simParams[0], simParams[1], simParams[2]))


print('All done!')
