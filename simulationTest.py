
# coding: utf-8

# In[1]:


#get_ipython().magic(u'matplotlib notebook')
#%matplotlib inline
# standard imports 
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


# In[2]:


import AIA_tools as at


# In[3]:


#help(at.simulation)


# In[4]:


pickles = glob('./*.pickle')
print(pickles) 
pickles = {'boxInfo':pickles[0], 'intensities':pickles[1], 'times':pickles[2]}
#print(pickles)


# In[5]:


data = {'boxInfo':None, 'intensities':None, 'times':None}    

for key in data: 
    #print(key) 
    file_ = open(pickles[key], 'r') 
    data[key] = pickle.load(file_ ) 

    


# ** get array of time steps and then create a time axis in seconds ** 

# In[6]:


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


# In[7]:


#print(times['171'])


# ** Plot the observational data ** 

# In[8]:


intensities = data['intensities']

fig, ax = plt.subplots() 

ax.plot(times['171'], intensities['171'], label='$171$ $\AA$') 
ax.plot(times['193'], intensities['193'], label='$193$ $\AA$')
ax.plot(times['211'], 5*np.asarray(intensities['211']), label='$211 \cdot 5$ $\AA$')
ax.plot(times['304'], 50*np.asarray(intensities['304']), label='$304 \cdot 50$ $\AA$')
ax.plot(times['335'], 125*np.asarray(intensities['335']), label='$335 \cdot 125$ $\AA$') 
ax.set_ylabel("DN $\mathrm{s}^{-1} \mathrm{pix}^{-1}$", fontsize=18)
ax.set_xlabel("Time [s]", fontsize=18)
ax.set_xlim(541, 1200) 
ax.legend() 


# ** Now take ratios and plot -- FIXME -- ask how to deal with times for ratio plots (interpolation?) ** 

# In[9]:


intensities = data['intensities']

fig, ax = plt.subplots(nrows=2, ncols=2) 

ax[0,0].plot(np.mean([times['193'], times['171']], axis=0),
             np.divide(intensities['193'],intensities['171']), label='$193/171$') 
ax[0,1].plot(np.mean([times['211'], times['193']], axis=0),
             np.divide(intensities['211'],intensities['193']), label='$211/193$')
ax[1,0].plot(np.mean([times['304'], times['193']], axis=0),
             np.divide(intensities['304'],intensities['193']), label='$304/193$')
ax[1,1].plot(np.mean([times['335'], times['304']], axis=0),
             np.divide(intensities['335'],intensities['304']), label='$335/304$') 

ax[0,0].set_ylabel("Intensity Ratio")
ax[0,0].set_xlabel("Time [s]")
ax[0,0].set_xlim(0, 1200)
ax[0,0].set_title('$193/171$')

ax[0,1].set_ylabel("Intensity Ratio")
ax[0,1].set_xlabel("Time [s]")
ax[0,1].set_xlim(0, 1200)
ax[0,1].set_title('$211/193$')

ax[1,0].set_ylabel("Intensity Ratio")
ax[1,0].set_xlabel("Time [s]")
ax[1,0].set_xlim(0, 1200)
ax[1,0].set_title('$304/193$')

ax[1,1].set_ylabel("Intensity Ratio")
ax[1,1].set_xlabel("Tijme [s]")
ax[1,1].set_xlim(0, 1200)
ax[1,1].set_title('$335/304$')

plt.tight_layout()
#plt.show()


# ** Now let's try to run a simulation ** 

# In[10]:


# trying out temperature from Ma et al. paper 
te_sta = 1e6 
te_end = 2.8e6 
n = 1e7

element_list = '2, 6, 7, 8, 10, 12, 13, 14, 16, 18, 20, 26, 28'  # He, C, N, O, Ne, Mg, Al, Si, S, Ar, Ca, Fe, Ni
simDataFile = at.simulation.run(te_sta, te_end, n, num=13,
               indices=element_list, ntime=10, dt=0.75, filename='te_end-1e6__te_sta-2.8e6__n-1e7.dat')


# ** Now load in the data ** 

# In[11]:


simData = at.simulation.getSimulationData(simDataFile)
#print(np.shape(simData['fractions']))
#print(np.shape(simData['times']))


# ** So it looks like indexing is [time, element, ion]. Now lets grab the emissivities.** 

# In[12]:


pathToChiantiEmiss = '/data/khnum/REU2018/jwaczak/data/chiantiEmissData'
emissFiles = at.simulation.getEmissList(pathToChiantiEmiss)


# ** Figure out how to loop through our triple dictionary ** 

# for key in emissFiles.keys():
#     for key2 in emissFiles[key].keys(): 
#         for key3 in emissFiles[key][key2].keys(): 
#             #print(key, key2, key3, emissFiles[key][key2][key3])

# ** Yay the file lookup dictionary worked! Test plotting for simulation output for one element at one time.** 

# In[13]:


obs_times = simData['times']
#print(obs_times) 


# fig, ax = plt.subplots() 
# plt.yscale("log")
# plt.ylim(1.0e-5, 1.0)
# 
# 

# ** load in the AIA response sav file and abundance files ** 

# In[14]:


from scipy.io import readsav 

eff_area = readsav('/data/khnum/REU2018/jwaczak/data/June_2010_CME_burst/aia_response/eff_area.sav')
aia_channels = ['A94', 'A131', 'A171', 'A193', 'A211', 'A304', 'A335']

log10_abund = {} 
pathToAbund = '/data/khnum/REU2018/jwaczak/data/CHIANTI_8.0.7_database/abundance/sun_coronal_1992_feldman.abund'
with open(pathToAbund, 'r') as f: 
    for line in f: 
        splitline = line.split(' ')
        lineData = [] 
        for thing in splitline: 
            if not (thing is ''): 
                lineData.append(thing) 
        if lineData[0] == '-1\n': 
            break
        else:
            log10_abund.update({lineData[0]: float(lineData[1])})
                    
#print(log10_abund)


# In[15]:


# now get actual abundance 
abund = {} 
for key in log10_abund.keys(): 
    abund.update({key: np.power(10.0,(log10_abund[key]-12))})
    
#print(abund)
#rint(len(abund))


# ** Let's open one of these emissivity files to see what's going on ** 

# In[16]:


testEmiss = at.simulation.getEmissData(emissFiles[2][2][7])
em_2d = testEmiss['em_2d']
wavelengths = testEmiss['lambda_1d']
log_temps = testEmiss['logte_1d']
temps = np.power(10.0, log_temps) 
fig, ax = plt.subplots(nrows=1, ncols=2) 
ax[0].plot(wavelengths, 'k.')
ax[0].set_title('emiss wavelengths [$\AA$]')
ax[1].plot(temps) 
ax[1].set_title('emiss temp [K]')


print("Nearest temperature to te_end is:{}".format(temps[at.analysis.getNearestValue(temps, te_end)]))
print(wavelengths.max(), wavelengths.min())


# ** Now lets loop through and create the synthetic observation ** 

# In[17]:


import multiprocessing as mp 
from functools import partial 
nproc = int(mp.cpu_count()/4.0)
print(nproc) 


# In[18]:


obs = [] 

def mp_wrapper(inp): 
    return at.simulation.getSyntheticObservation(*inp)

args = [(timeIndex, te_sta, te_end, 7, eff_area  )for timeIndex in range(len(obs_times))]

#args = [(timeIndex, te_sta, te_end, 7,  eff_area, emissFiles, abund)for timeIndex in range(len(obs_times))]

for args_ in args[:1]: 
    obs.append(mp_wrapper(args_))






# ** It looks like the wrapper is working ** 

# In[ ]:


pool = mp.Pool(processes = 4) 
out = pool.map(mp_wrapper, args) 
pool.close() 
pool.join() 


# In[ ]:


print(obs)

