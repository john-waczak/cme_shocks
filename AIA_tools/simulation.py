from __future__ import division
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
    simParams = [te_sta, te_end, ne]

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


def run(te_sta = 1.8*1e6, te_end = 2.8*1e6, n = 1e7, num = 2, indices = '2, 26', ntime = 10, dt = 0.75,
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


def getCoronalAbundances(pathToAbund = '/data/khnum/REU2018/jwaczak/data/abundance/sun_coronal_1992_feldman.abund'):
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
            wav_indices = wav_indices[(emiss_wavelengths[wav_indices]<=400) & (emiss_wavelengths[wav_indices]>=90)]

            # get correct emissivity data
            emiss_wavelengths = emiss_wavelengths[wav_indices]
            emiss0 = emissData['em_2d'][wav_indices, temp_index]

            # calculate new emissivity using ion fraction and abundances
            emiss0 = emiss0*simData['fractions'][timeIndex][elem-1, ion]*abundances[str(elem)]

            # loop through channels to get synthetic counts
            #for channel in aia_channels:

            for channel in observation.keys():
                if channel is not 'time': 
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



def mp_wrapper(inp):
        return getSyntheticObservation(*inp)


def runParallel(simDataFile, nproc, obs_times, te_sta, te_end, n):

    Args = [(timeIndex, te_sta, te_end, n, simDataFile) for timeIndex in range(len(obs_times))]
    pool = mp.Pool(processes = nproc)
    obs = pool.map(mp_wrapper, Args) 
    pool.close()
    pool.join()

    data = np.asarray(obs)

    fname= '/data/khnum/REU2018/jwaczak/data/simOutput/t0--{:.2E}__t1--{:.2E}__n--{:.2E}'.format(te_sta, te_end, 10**n) 
    np.savetxt(fname+'.txt', data, delimiter=',', header='time, 171, 193, 211, 304, 335')

    return data


def getSyntheticObservation_II(te_sta, te_end, n, simDataFile):
    # get the simulation data
    simData = getSimulationData(simDataFile)

    n_s = [7, 8, 9]
    n_i = analysis.getNearestValue(np.asarray(n_s), n)
    N = n_s[n_i]


    simParams = simData['simParams']

    # get the time values from the data
    time_vals = simData['times']

    # create emissivity file list
    pathToChiantiEmiss = '/data/khnum/REU2018/jwaczak/data/chiantiEmissData'
    emissFiles = getEmissList(pathToChiantiEmiss) 

    # get wavelength response from AIA
    eff_area = readsav('/data/khnum/REU2018/jwaczak/data/June_2010_CME_burst/aia_response/eff_area.sav')

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

            # get index to temp nearest simulation te_end
            temp_index = analysis.getNearestValue(np.power(10.0, emiss_log_temps), te_end)

            # truncate the wavelengths to allow for interpolation
            wav_indices = np.asarray(range(len(emiss_wavelengths)))
            wav_indices = wav_indices[(emiss_wavelengths[wav_indices]<=400) & (emiss_wavelengths[wav_indices]>=90)]

           # get the correct emissivity data
            emiss_wavelengths = emiss_wavelengths[wav_indices]
            emiss0 = emissData['em_2d'][wav_indices, temp_index]  # =A_ji*(n_j/n_ion)(1/(n_e*4*pi))  [photons s^-1 cm^3 sr^-1]
            ea_interpolated = {}


# --------------------------------------------------#

            testDict = {26:{16:[], 8:[], 21:[], 22:[], 12:[]}, 12:{8:[]}, 8:{6:[]}} 

# --------------------------------------------------#

 
            testOut = [] 
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



                        # if elem in testDict.keys():
                        #     if ion in testDict[elem].keys():
                        #         if channel is '335' and (t % 5 is 0):
                        #             print("\t{} {} {} {}".format(elem, ion, channel, t))
                        #             testDict[elem][ion].append('channel: {}  t: {:.2f}  em0_new max: {} [photons s^-1 cm^5 sr^-1]  em_sum: {} [photons s^-1 cm^5 sr^-1]  emout: {} [photons s^-1 cm^5 sr^-1 (n_ion/n_elem)*(n_elem)]'.format(channel, time_vals[t], np.asarray(em0_new).max(), em_sum, emout))


                        obs[channel][t] += emout

            # if elem in testDict.keys():
            #     if ion in testDict[elem].keys():
            #         with open ('/data/khnum/REU2018/jwaczak/data/emissivityTests/{}_{}.txt'.format(elem, ion), 'w') as f:
            #             for line in testDict[elem][ion]:
            #                 f.write(line+'\n')

    return obs, simParams


def getTempVals(pathToTemps = '/data/khnum/REU2018/jwaczak/data/tempVals.txt'):
    temps = np.loadtxt(pathToTemps, delimiter=',')
    return temps


def getRadiusToShock(timeIndex, dt, r0, v_shock):
    return r0+v_shock*timeIndex*dt


def getShockSpeed(timeIndex):
    return 600.0 #  km/s


def applySphericalCorrection(syntheticObservationFile, dx, t0, t1, n, limb=False, X=1.56):
    pathToBackground = '/data/khnum/REU2018/jwaczak/data/background.txt'
    pathToInitialRadii = '/data/khnum/REU2018/jwaczak/data/initialRadii.txt'
    outputPath = '/data/khnum/REU2018/jwaczak/data/correctedSyntheticObservations/'

    syntheticObservation_Old = np.loadtxt(syntheticObservationFile, delimiter=',')
    background = np.loadtxt(pathToBackground, delimiter=',')
    R0 = np.loadtxt(pathToInitialRadii, delimiter=',')

    print(np.shape(syntheticObservation_Old), np.shape(background), np.shape(R0))


    if limb is False:
        r0 = R0[0]
    else:
        r0 = R0[1]

    dt = 12.0  # AIA time cadence


    #---- make a vector of shock distance at each time step -------#
    shockDistances = [r0]
    for i in range(1,len(syntheticObservation_Old[:,0])):
        #r_new = shockDistances[i-1]+dt*getShockSpeed(i)/X  # divide by compression ratio -- see notes
        shockDistances.append(getRadiusToShock(i, dt, r0, getShockSpeed(i)))

    print(shockDistances[0], shockDistances[-1])

    #---- loop through data and apply correction ---------#
    outData = []
    for i in range(len(syntheticObservation_Old[:,0])):
        print(i) 

        v_shock = getShockSpeed(i)
        r = getRadiusToShock(i, dt, r0, v_shock) 
        l = 2*np.sqrt(r**2-r0**2) # NOTE: x_i = r0

        I_back = background*(1-(l/(2*r)))  # note the background has 5 values, 1 for each channel 
        for I in I_back:
            if I < 0:
                I = 0 

        N_w = int(l/dx) # number of cells

        if N_w is 0: N_w = 1  # insure we always have at least 1 box
        if N_w > 1 and N_w % 2 ==0:  # insure we always have odd number of cells 
            N_w = N_w -1

        print('  {}'.format(N_w)) 


        cells = [[r0, 0.0]]
        if N_w > 1:
            for j in range(int((N_w-1)/2)):
                cells.append([r0, j*dx])
                cells.append([r0, -1*j*dx])

        d_i = [(r-np.sqrt(c[0]**2+c[1]**2)) for c in cells]  # distance of each shell to the shock
        d_i = np.asarray(d_i) 


        # now we need to get total emissivity for all celss in each channel 
        cell_emissivities = []
        for j in range(1,6):
            em = []
            for d in d_i:
                index = analysis.getNearestValue(shockDistances-r0, d)
                em_adjusted = syntheticObservation_Old[index,j]
                em.append(em_adjusted)  # modulate by compression ratio, original density,
            c = (1/3600)**2*(np.pi/180)**2*(1.5)
            em_new = np.sum(np.asarray(em)*c*(X**2)*(n**2)*dx*(10**3)) #+ I_back[j-1]
            em_w_back = em_new #+ I_back[j-1]
            print('\t{} {}'.format(em_new, em_w_back))
            cell_emissivities.append(em_w_back)
        # add new values to output data
        o_dat = np.array([syntheticObservation_Old[i,0], cell_emissivities[0], cell_emissivities[1],
                        cell_emissivities[2], cell_emissivities[3], cell_emissivities[4]])
        outData.append(o_dat)

    outData = np.asarray(outData)
    fileName = outputPath+'t0--{:.2E}__t1--{:.2E}__n--{:.2E}.txt'.format(t0, t1, n)
    np.savetxt(fileName, outData, delimiter=',')
    return outData, fileName


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





















