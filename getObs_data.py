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






# get the observation data

pickles = glob('./*.pickle')
print(pickles) 
pickles = {'boxInfo':pickles[0], 'intensities':pickles[1], 'times':pickles[2]}
data = {'boxInfo':None, 'intensities':None, 'times':None}
for key in data:
    file_ = open(pickles[key], 'r')
    data[key] = pickle.load(file_ )

# get the times
times = {'171':[], '193':[], '211':[], '304':[], '335':[]}

for key in times:
    # start the times at zero seconds and add that to the list 
    tot = 0.0
    times[key].append(tot)

    # get difference and add to tot
    for i in range(1, len(data['times'][key])): 
        t_diffs =(data['times'][key][i]-data['times'][key][i-1]).total_seconds()
        tot += t_diffs
        times[key].append(tot) 

intensities = data['intensities']

# times go from 540 to 1200 [s] 
# use the 171 A times as that is how we defined the start of the cme 

# create interpolation objects 
interp_171 = interp1d(times['171'], intensities['171'])
interp_193 = interp1d(times['193'], intensities['193']) 
interp_211 = interp1d(times['211'], intensities['211'])
interp_304 = interp1d(times['304'], intensities['304']) 
interp_335 = interp1d(times['335'], intensities['335']) 


# interpolate the times to the 171 times 
timeIndex_low = aia.analysis.getNearestValue(np.asarray(times['171']), 540)
timeIndex_high = aia.analysis.getNearestValue(np.asarray(times['171']), 1100) 
time_vals =  np.asarray(times['171'][timeIndex_low:timeIndex_high])
intensities['171'] = intensities['171'][timeIndex_low:timeIndex_high]
intensities['193'] = interp_193(time_vals) 
intensities['211'] = interp_211(time_vals) 
intensities['304'] = interp_304(time_vals) 
intensities['335'] = interp_335(time_vals) 

Dt = np.abs(time_vals[1:]-time_vals[:-1]).mean()
time_vals = np.arange(0, len(time_vals)*Dt, Dt)


obs_out = np.asarray([time_vals, intensities['171'],
            intensities['193'], intensities['211'],
            intensities['304'], intensities['335']])

obs_out = np.transpose(obs_out) 
np.savetxt('/data/khnum/REU2018/jwaczak/data/obs_data.txt', obs_out, delimiter=',')

