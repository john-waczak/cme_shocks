from __future__ import division, print_function
import numpy as np
import matplotlib.pyplot as plt
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
import multiprocessing as mp


# global constants
# see Ma_et_al.ipynb for calculation from fits file meta data

R_w = 898058.574287  # [km]
r_w = 219174.480844  # [km]
R_sun = 696000.0     # [km] i.e. km per solar Radii
arcsec_per_pixel = 0.600000023842  # [arcsec / pixel]
arcsec_per_solarRadii = 944.760129  # [arcsec / R_sun]
c = (1/3600)**2 * (np.pi/180)**2 * 0.6**2  # (square deg per square arcsec)(steradian per square deg)(square arcsec per pixel) i.e. (steradian per pixel)
X = 1.56  # shock compression ratio 

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
    [Te_sta, Te_end, ne] = f.read_reals(dtype=np.float64)
    simParams = [Te_sta, Te_end, ne]

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

    simData = {'times':times, 'fractions': NEI_fractions, 'simParams':simParams}
    return simData


def run(Te_sta = 1.8*1e6, Te_end = 2.8*1e6, n = 1e7, num = 2, indices = '2, 26', ntime = 10, dt = 0.75,
        filename = 'nei_onestep.dat',
        outputPath = '/home/john/gitRepos/REU/jwaczak/data/simRuns',
        pathToSimCode = '/home/john/gitRepos/REU/jwaczak/time_dependent_fortran/Applications/NEI_Onestep',
        pathToChianti = "\'/home/john/gitRepos/REU/jwaczak/time_dependent_fortran/eigendata/chianti_8_07/\'"):

    # open the main.f90 file
    mainFile = open('{}/main.f90'.format(pathToSimCode), 'r')
    file_lines = mainFile.read().splitlines()

    # find parameters and change them
    for i in range(len(file_lines)):
        if ("te_start = ") in file_lines[i]:
            file_lines[i] = "  te_start = {}\t! (K) used to set the initial charge state".format(Te_sta)
        elif ("te_end = ") in file_lines[i]:
            file_lines[i] = "  te_end = {}\t! (K)".format(Te_end)
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


def getCoronalAbundances(pathToAbund = '/home/john/gitRepos/REU/jwaczak/data/abundance/sun_coronal_1992_feldman.abund'):
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


def getSyntheticObservation_II(Te_sta, Te_end, n, simDataFile):
    # get the simulation data
    simData = getSimulationData(simDataFile)

    n_s = [7, 8, 9]  # densities in emissivity tables
    n_i = analysis.getNearestValue(np.asarray(n_s), n)  # index of closest density
    N = n_s[n_i]


    simParams = simData['simParams']

    # get the time values from the data
    time_vals = simData['times']

    # create emissivity file list
    pathToChiantiEmiss = '/home/john/gitRepos/REU/jwaczak/data/chiantiEmissData'
    emissFiles = getEmissList(pathToChiantiEmiss)

    # get wavelength response from AIA
    eff_area = readsav('/home/john/gitRepos/REU/jwaczak/data/June_2010_CME_burst/aia_response/eff_area.sav')

    # get the coronal abundance data
    abundances = getCoronalAbundances()

    elem_list = [2, 6, 7, 8, 10, 12, 13, 14, 16, 18, 20, 26, 28]

    # set up output dictionary
    obs = {'time':simData['times'], '171':np.zeros(len(time_vals)), '193':np.zeros(len(time_vals)),
           '211':np.zeros(len(time_vals)), '304':np.zeros(len(time_vals)), '335':np.zeros(len(time_vals))}

    # loop through all of the elements
    for elem in elem_list:
        print(elem)
        # make sure we are actually simulating the ion
        ions = [ion for ion in range(30) if ion in emissFiles[elem].keys()]

        # loop through each ion
        for ion in ions:
            # get emiss data
            emissData = getEmissData(emissFiles[elem][ion][N])
            emiss_wavelengths = emissData['lambda_1d']
            emiss_log_temps = emissData['logte_1d']

            # get index to temp nearest simulation Te_end
            temp_index = analysis.getNearestValue(np.power(10.0, emiss_log_temps), Te_end)

            # truncate the wavelengths to allow for interpolation
            wav_indices = np.asarray(range(len(emiss_wavelengths)))
            wav_indices = wav_indices[(emiss_wavelengths[wav_indices]<=400) & (emiss_wavelengths[wav_indices]>=90)]

           # get the correct emissivity data
            emiss_wavelengths = emiss_wavelengths[wav_indices]
            emiss0 = emissData['em_2d'][wav_indices, temp_index]  # =A_ji*(n_j/n_ion)(1/(n_e*4*pi))  [photons s^-1 cm^3 sr^-1]
            ea_interpolated = {}

            #create multidimensional array with emiss0*eff_area for each AIA band
            for channel in obs.keys():
               if channel is not 'time':
                   ea_data = eff_area['effarea']['A'+channel][0]
                   ea_wavelengths = ea_data['wave'][0]
                   ea_values = ea_data['ea'][0]

                   # interpolate to match emissivity wavelengths
                   interp_obj = interp1d(ea_wavelengths, ea_values, kind='quadratic')
                   interp_ea = interp_obj(emiss_wavelengths)  # I think this has units of [cm^2]

                   # get real emissivity using aia response fucntion
                   ea_interpolated.update({channel:interp_ea})

            for t in range(len(time_vals)):
                for channel in obs.keys():
                    if channel is not 'time':
                        em0_new = emiss0*ea_interpolated[channel]  # [photons s^-1 cm^5 sr^-1]
                        em_sum = np.sum(np.asarray(em0_new))
                        emout = em_sum*simData['fractions'][t][elem-1, ion-1]*abundances[str(elem)]  # [photons s^-1 cm^5 sr^-1 (n_ion/n_elem)*(n_elem)]
                        obs[channel][t] += emout
    return obs, simParams


def getTempVals(pathToTemps = '/home/john/gitRepos/REU/jwaczak/data/tempVals.txt'):
    temps = np.loadtxt(pathToTemps, delimiter=',')
    return temps


def getRadiusToShock(timeIndex, dt, r0, v_shock):
    return r0+v_shock*timeIndex*dt


def getShockSpeed(timeIndex):
    return 600.0 #  km/s


def getDensity(h, R_w, R_sun):  # k is a fitting constant
    """ given a height along line of sight, the distance to the observation window
    and the solar radius, return the coronal density based off of streamer data fromLimb

    Gibson, S.E., Fludra, A., et all 1999: Solar minimum streamer densities and temperatures...
    """

    R = np.sqrt(R_w**2+h**2)
    r = R/R_sun  # units need to be in solar radii 
    a = 77.1
    b = 31.4
    c = 0.954
    d = 8.30
    e = 0.550
    f = 4.63

    return (a*r**(-b) + c*r**(-d) + e*r**(-f))*10**8  #[cm-3]


# def fitBackgroundEmissivity(pathToBackground, R_w, R_sun, N_w, n_sim, R_max = 1.5, tol=0.05): #e.g. tolerance is 5%  k is a constant for adjusting density
#     background = np.loadtxt(pathToBackground, delimiter=',')
#     print('Background: {}'.format(background))
#     background_new = []
#     for channel in background:
#         dh = (R_max*R_sun) / (N_w/2)  # cell size [km]
#         densities = []
#         for j in range(int(N_w/2)):
#             densities.append(getDensity(j*dh, R_w, R_sun))  # density at each cell

#         # #--- adjust the density to match the simulation ---#
#         # densities = np.asarray(densities)
#         # adjust = (float(n_sim))/np.amax(densities)
#         # densities = [adjust*d for d in densities]


#         def get_error(em):
#             bckgnd_tot = 0
#             for j in range(int(N_w/2)):
#                 bckgnd_tot += X**2 * densities[j]**2 * dh * em * c * (10.0**5)  # include c and 10^5 to convert km to cm
#             return channel - bckgnd_tot

#         # establish lower and upper bound
#         em_max = 1
#         em_min = em_max

#         while(get_error(em_max) > 0):
#             em_max += 10*em_max

#         while(get_error(em_min) < 0):
#             em_min -= 10*em_min

#         # reduce until within tolerance
#         em = em_max
#         tolerance = tol*channel  # define the tolerance in terms of the background level
#         print(get_error(em), tolerance)
#         while(np.abs(get_error(em))>tolerance):
#             em = (em_max+em_min)/2
#             error = get_error(em)
#             if (error<0):
#                 em_max = em
#             elif (error>0):
#                 em_min = em
#             else:
#                 print("something went wrong! em: {}  error: {}".format(em, error))

#         background_new.append(em)
#     return background_new



def applySphericalCorrection2(syntheticObservationFile, N_w, T0, T1, n, R_max = 1.5,  X=1.56):
    pathToBackground = '/home/john/gitRepos/REU/jwaczak/data/background.txt'
    outputPath = "/home/john/gitRepos/REU/jwaczak/data/correctedSyntheticObservations/"
    bckgnd = np.loadtxt(pathToBackground, delimiter=',')
    syntheticObservation_Old = np.loadtxt(syntheticObservationFile, delimiter=',')

    h_max = np.sqrt((R_max*R_sun)**2-R_w**2)  # [km] maximum height of line of sight -- total length is 2*h_max

    cells = []  #  define the cells of plasma along the line of sight
    densities = []  # collect density at each cell of plasma
    dh = (2*h_max) / (N_w)
    for j in range(int(N_w/2)):
        densities.append(getDensity(j*dh, R_w, R_sun))  # density at each cell
        cells.append(j*dh)


    # #--- adjust the density to match the simulation ---#
    # we are going to fit the streamer density to the initial values of the simulation.
    ratio = []
    for chan in range(1,6):
        I_0 = np.sum([ X**2 * densities[j]**2 * dh * syntheticObservation_Old[0,chan] * c * (10.0**5) for j in range(len(cells))])
        ratio.append(bckgnd[chan-1]/I_0)
    ratio = np.asarray(ratio)
    d_factor = np.mean(ratio)
    print(ratio, d_factor)
    densities = [d*np.sqrt(d_factor) for d in densities] # square root factor due to density squared depend.

    dt = 12.0  # [s]  AIA time cadence

    #--------------------------------------------------------------------------------------------#
    # loop through each time step
    outData = []
    for i in range(len(syntheticObservation_Old[:,0])):
        v_s = getShockSpeed(i)  # [km/s] shock speed
        v_p = v_s*(1-(1/X))  # [km/s] speed of plasma behind shock
        r_s = getRadiusToShock(i, dt, r_w, getShockSpeed(i))  # distance from window to shock  [km]

        # NOTE I am only working with Half of the line of sight due to symmetry.

        l = np.sqrt(r_s**2-r_w**2)  # [km] length of line of sight inside of shock

        times_corrected = []  # correct times are independent of channel
        for j in range(int(N_w/2)):
            d = r_s - np.sqrt(r_w**2 + cells[j]**2)  # [km]
            t_obs = (X*d)/v_s  # [s]
            times_corrected.append(t_obs)

        # loop through each channel and determine new emissivity
        cell_emissivities = []
        cell_emissivities.append(syntheticObservation_Old[i,0])  # add the time to the data

        for chan in range(1,6):  # time i channel j
            em = []
            # loop through the cells and grab the  appropriate emissivity.

            # at our disposal are cells, times, densities, and shock distances
            for j in range(len(cells)):
                if cells[j] <= l:
                    t  = analysis.getNearestValue(syntheticObservation_Old[:,0], times_corrected[j])
                    em_new = syntheticObservation_Old[t, chan]
                    em.append(em_new)

                elif cells[j] > l:
                    em_new = bckgnd[chan-1]
                    em.append(bckgnd[chan-1])

                else:
                    print("Something went wrong! Could not assign an emissivity")
            Intensity_real = np.sum([ X**2 * densities[j]**2 * dh * em[j] * c * (10.0**5) for j in range(len(cells))])
            cell_emissivities.append(Intensity_real)

        outData.append(cell_emissivities)

    outData = np.asarray(outData)
    fileName = outputPath+"T0--{0:.2E}__T1--{1:.2E}__n--{2:.2E}.txt".format(float(T0), float(T1), float(n))
    return outData, fileName





if __name__ == '__main__':
    print("Testing...")






















