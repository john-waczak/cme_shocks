import numpy as np
#import sunpy.map as smap
import os, sys
#from astropy.coordinates import SkyCoord
#from astropy.visualization.mpl_normalize import ImageNormalize  
#import astropy.units as u
from glob import glob
#import matplotlib.pyplot as plt
#import matplotlib.dates as mdates
from scipy.interpolate import interp1d
from astropy.convolution import convolve, Box1DKernel
from datetime import datetime, timedelta 
from scipy.io.idl import readsav
from scipy.io import FortranFile 
import re
import analysis 


def getEmissList(pathToSavFiles):
    fileList = glob(pathToSavFiles+'/*.sav')
    savFiles = {} 
 
    # loop through .save files and categorize into a big dictionary
    # format: savFiles[iz][ion][n]
    for file_ in fileList:
        splitLine = re.split('_|\.', file_)
        iz = int(splitLine[1][2:])
        ion = int(splitLine[2][3:])
        n = int(splitLine[3][1:])

        if not (iz in savFiles.keys()):
            savFiles.update({iz:{ion:{n:file_}}})
        else:
            if not (ion in savFiles[iz].keys()):
                savFiles[iz].update({ion:{n:file_}})
            else:
                if not (n in savFiles[iz][ion].keys()):
                    savFiles[iz][ion].update({n:file_})


    return savFiles


def getEmissData(pathToSavFile):
    data = readsav(pathToSavFile)
    return data


def getSimulationData(dotDatFile):
    try:
        f = FortranFile(dotDatFile, 'r')
    except:
        raise SystemExit

    # get simulation parameters
    [te_sta, te_end, ne] = f.read_reals(dtype=np.float64)

    # get initial and final fractions for equilibrium ionization calculation 
    fraction_initial_ei  = f.read_reals(dtype=np.float64).reshape(30,30) # reshape for 30 elements by 30 ions grid
    fraction_ei_final = f.read_reals(dtype=np.float64).reshape(30,30)

    n_timeSteps = f.read_ints()
    times = []
    NEI_fractions = []

    # add initial fractions to lists
    times.append(0.0)
    NEI_fractions.append(fraction_initial_ei)

    # loop through the data and append to our lists
    for i in range(n_timeSteps[0]): # unfortunately n_timeSteps in an array with a single element
        time = f.read_reals(dtype=np.float64)
        current_nei_fractions = f.read_reals(dtype=np.float64).reshape(30,30)
        times.append(time[0])
        NEI_fractions.append(current_nei_fractions)

    # add the final ei fractions to the list
    NEI_fractions.append(fraction_ei_final)

    #*** note that the NEI_fractions list now has one extra element as it isn't reasonable to attach a time to the
    #*** final equilibrium ionization fractions

    simData = {'times':times, 'fractions': NEI_fractions}
    return simData


def run(te_sta = 1e6, te_end = 2.8*1e6, n = 1e7, num = 2, indices = '2, 26', ntime = 10, dt = 0.75,
        filename = 'nei_onestep.dat',
        outputPath = '/data/khnum/REU2018/jwaczak/data/simRuns',
        pathToSimCode = '/data/khnum/REU2018/jwaczak/time_dependent_fortran/Applications/NEI_Onestep',
        pathToChianti = "\'/data/khnum/REU2018/jwaczak/time_dependent_fortran/eigendata/chianti_8_07/\'"): 

    # open the main.f90 file 
    mainFile = open('{}/main.f90'.format(pathToSimCode), 'r')
    file_lines = mainFile.read().splitlines()

    # find parameters and change them 
    for i in range(len(file_lines)):
        if ("te_start = ") in file_lines[i]:
            file_lines[i] = "  te_start = {}\t! (K) used to set the initial charge state".format(te_sta)
        elif ("te_end = ") in file_lines[i]:
            file_lines[i] = "  te_end = {}\t! (K)".format(te_end)
        elif ("rhone = ") in file_lines[i]:
            file_lines[i] = "  rhone = {}\t! (cm^-3)".format(n)
        elif ("open(12, file=") in file_lines[i]:
            file_lines[i] = "open(12, file='{}', form=\'unformatted\')".format(filename)
        elif ("integer, parameter:: ntime=") in file_lines[i]:
            file_lines[i] = "integer, parameter:: ntime={}".format(ntime)
        elif ("real*8, parameter:: dt0 =") in file_lines[i]:
            file_lines[i] = "real*8, parameter:: dt0 = {}".format(dt) 

    mainFile.close() 

    # re open the file in write mode and update the file 
    mainFile = open('{}/main.f90'.format(pathToSimCode), 'w')
    mainFile.write('\n'.join(file_lines))
    mainFile.close()

    # open the init.txt file and change the input parameters.
    inputFile = open('{}/input.txt'.format(pathToSimCode), 'r')
    file_lines = inputFile.read().splitlines()

    # the lines we need to change are at index 1 and 2.
    file_lines[1] = str(num)
    file_lines[2] = indices
    file_lines[4] = pathToChianti 
    inputFile.close()

    # re open file in write mode and update the file
    inputFile = open('{}/input.txt'.format(pathToSimCode), 'w')
    inputFile.write('\n'.join(file_lines))
    inputFile.close()

    # -----------------------------------------------------------------------------------
    # Now we want to recompile using 'make all'
    # -----------------------------------------------------------------------------------
    os.chdir(pathToSimCode) 
    os.system("make all")
    os.system("./main")
    os.chdir("../..")
    os.system("mv {}/{} {}".format(pathToSimCode, filename, outputPath))

    return "{}/{}".format(outputPath, filename) 



def getCoronalAbundances(pathToAbund = '/data/khnum/REU2018/jwaczak/data/CHIANTI_8.0.7_database/abundance/sun_coronal_1992_feldman.abund'):
    log10_abund = {}
    with open(pathToAbund, 'r') as f:
        for line in f:
            splitline = line.split(' ')
            lineData = []
            for thing in splitline:
                if not (thing is ''):
                    lineData.append(thing)
            if len(lineData) == 1:
                break
            else:
                log10_abund.update({lineData[0]: float(lineData[1])})
    abund = {}
    for key in log10_abund.keys():
        abund.update({key: np.power(10.0, log10_abund[key]-12.0)})

    return abund 


def getSyntheticObservation(timeIndex, te_sta, te_end, n, simDataFile): 
    print(timeIndex)


    # get the simulation data
    simData = getSimulationData(simDataFile) 

    # get emissivity files
    pathToChiantiEmiss = '/data/khnum/REU2018/jwaczak/data/chiantiEmissData'
    emissFiles = getEmissList(pathToChiantiEmiss) 

    # get wavelength response from AIA
    eff_area = readsav('/data/khnum/REU2018/jwaczak/data/June_2010_CME_burst/aia_response/eff_area.sav')

    # get the coronal abundance data 
    abundances = getCoronalAbundances()

    # create observation data structure for output
    observation = {'time':simData['times'][timeIndex], '171':0.0, '193':0.0, '211':0.0, '304':0.0, '335':0.0}

    aia_channels = ['A94', 'A131', 'A171', 'A193', 'A211', 'A304', 'A335']
    elem_list = [2, 6, 7, 8, 10, 12, 13, 14, 16, 18, 20, 26, 28]


    # loop through the elements in the list
    for elem in elem_list:
        print('\t{}'.format(elem))

        ions = [ion for ion in range(30) if ion in emissFiles[elem].keys()]
        for ion in ions: 

            # get emissivity data
            emissData = getEmissData(emissFiles[elem][ion][n])
            emiss_wavelengths = emissData['lambda_1d']
            emiss_log_temps = emissData['logte_1d'] 

            # get nearest temperature index
            temp_index = analysis.getNearestValue(np.power(10.0, emiss_log_temps), te_end)

            # truncate wavlengths to allow for interpolation (scipy doesn't automatically extrapolate)
            wav_indices = np.asarray(range(len(emiss_wavelengths)))
            wav_indices = wav_indices[(emiss_wavelengths[wav_indices]<=900) & (emiss_wavelengths[wav_indices]>=25)]

            # get correct emissivity data
            emiss_wavelengths = emiss_wavelengths[wav_indices]
            emiss0 = emissData['em_2d'][wav_indices, temp_index]

            # calculate new emissivity using ion fraction and abundances
            emiss0 = emiss0*simData['fractions'][timeIndex][elem-1, ion]*abundances[str(elem)]

            # loop through channels to get synthetic counts
            #for channel in aia_channels:

            for channel in observation.keys():
                if not channel == 'time': 
                    ea_data = eff_area['effarea']['A'+channel][0]
                    ea_wavelengths = ea_data['wave'][0]
                    ea_values = ea_data['ea'][0]

                    # interpolate to match emissivity wavelengths
                    interp_obj = interp1d(ea_wavelengths, ea_values, kind='quadratic')
                    interp_ea = interp_obj(emiss_wavelengths)

                    # get real emissivity using aia response fucntion
                    em_new = interp_ea * emiss0

                    # get total counts by summing up whole band
                    em_tot = np.sum(em_new)

                    # add this to appropriate spot in observation dictionary
                    observation[channel] += em_tot 

    return [observation['time'], observation['171'], observation['193'], observation['211'], observation['304'], observation['335'] ]




if __name__ == '__main__':
    pathToSavFiles = "/data/khnum/REU2018/jwaczak/data/chiantiEmissData"
    savFiles = getIonList(pathToSavFiles) 
    print(savFiles[28][20][9])
    data = getIonData(savFiles[28][20][9])
    wavelength = data.lambda_1d
    log_temps = data.logte_1d
    emissivities = data.em_2d
    File = run() 
    print(File)




















