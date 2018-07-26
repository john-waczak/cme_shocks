import numpy as np
import matplotlib.pyplot as plt
from glob import glob
import os
import AIA_tools as aia


pathToData = '/data/khnum/REU2018/jwaczak/data/simOutput/'
pathToFigFolder = '/data/khnum/REU2018/jwaczak/data/figs/simulation/'

os.chdir(pathToData)
dataFiles = glob('*.*')

obs_data = np.loadtxt('/data/khnum/REU2018/jwaczak/data/obs_data.txt', delimiter=',')
Dt = np.abs(obs_data[0,1:]-obs_data[0,:-1]).mean()
times = np.arange(0, len(obs_data[0,:])*Dt, Dt)


for file_ in dataFiles:
    data = np.loadtxt(file_, delimiter=',')
    print(np.shape(data), np.shape(obs_data), np.shape(times))
    fig, ax = plt.subplots(nrows=2, ncols=2)
    ax[0,0].plot(data[:,0], np.divide(data[:,1], data[:,2]),'k', label='simulation')
    ax[0,0].plot(times, np.divide(obs_data[1,:], obs_data[2,:]), 'b', label='observation')
    ax[0,1].plot(data[:,0], np.divide(data[:,3], data[:,2]),'k', label='simulation')
    ax[0,1].plot(times, np.divide(obs_data[3,:], obs_data[2,:]), 'b', label='observation')
    ax[1,0].plot(data[:,0], np.divide(data[:,4], data[:,2]), 'k', label='simulation')
    ax[1,0].plot(times, np.divide(obs_data[4,:], obs_data[2,:]), 'b', label='observation')
    ax[1,1].plot(data[:,0], np.divide(data[:,5], data[:,2]), 'k', label='simulation')
    ax[1,1].plot(times, np.divide(obs_data[5,:], obs_data[2,:]), 'b', label='observation')


    ax[0,0].legend()
    ax[0,1].legend()
    ax[1,0].legend()
    ax[1,1].legend()

    ax[0,0].set_ylabel("Intensity Ratio")
    ax[0,0].set_xlabel("Time [s]")
    ax[0,0].set_title('$171/193$')

    ax[0,1].set_ylabel("Intensity Ratio")
    ax[0,1].set_xlabel("Time [s]")
    ax[0,1].set_title('$211/193$')

    ax[1,0].set_ylabel("Intensity Ratio")
    ax[1,0].set_xlabel("Time [s]")
    ax[1,0].set_title('$304/193$')

    ax[1,1].set_ylabel("Intensity Ratio")
    ax[1,1].set_xlabel("Time [s]")
    ax[1,1].set_title('$335/193$')

    fig.suptitle(file_[:-4])
    plt.tight_layout()
    plt.subplots_adjust(left=None, bottom=None, right=None, top=0.85,
                wspace=None, hspace=None)
    fileName = pathToFigFolder+'/intensityRatios/'+file_[:-4]
    plt.savefig(fileName+'.png')
    plt.close()

    fig_tot, ax_tot = plt.subplots()
    ax_tot.plot(data[:,0], data[:,1]/data[:,1].max(), 'k', label='$171$ $\AA$') 
    ax_tot.plot(times, obs_data[1,:]/obs_data[1,:].max(), 'k--')

    ax_tot.plot(data[:,0], data[:,2]/data[:,2].max(), 'b', label='$193$ $\AA$')
    ax_tot.plot(times, obs_data[2,:]/obs_data[2,:].max(), 'b--')

    ax_tot.plot(data[:,0], data[:,3]/data[:,3].max(), 'r',label='$211$ $\AA$')
    ax_tot.plot(times, obs_data[3,:]/obs_data[3,:].max(), 'r--')

    ax_tot.plot(data[:,0], data[:,4]/data[:,4].max(), 'g', label='$304$ $\AA$')
    ax_tot.plot(times, obs_data[4,:]/obs_data[4,:].max(), 'g--')

    ax_tot.plot(data[:,0], data[:,5]/data[:,5].max(), 'm', label='$335$ $\AA$')
    ax_tot.plot(times, obs_data[5,:]/obs_data[5,:].max(), 'm--')

    ax_tot.set_ylabel("Tot Intensity (normalized to 1) ")
    ax_tot.set_xlabel("Time [s]", fontsize=18)

    ax_tot.legend()
    ax_tot.set_title(file_[:-4])

    fileName = pathToFigFolder+'/totalIntensity/'+file_[:-4]
    plt.savefig(fileName+'.png')
    plt.close()

    simDataFile = '/data/khnum/REU2018/jwaczak/data/simRuns/'+ file_[:-4]+'.dat'
    print(os.path.isfile(simDataFile), simDataFile)
    simData = aia.simulation.getSimulationData(simDataFile)
    times_frac = simData['times']
    fractions = np.asarray(simData['fractions'])
    fractions = fractions[:-1, :, :]
    fe_9 = fractions[:, 25, 9]
    fe_10 = fractions[:, 25, 10] 
    fe_12 = fractions[:, 25, 12]
    fe_14 = fractions[:, 25, 14]
    fe_16 = fractions[:, 25, 16]
    si_12 = fractions[:, 13, 12]
    o_6 = fractions[:, 7, 6]

    plt.figure()
    plt.plot(times_frac, fe_9, label='Fe IX')
    plt.plot(times_frac, fe_10, label='Fe X')
    plt.plot(times_frac, fe_12, label='Fe XII')
    plt.plot(times_frac, fe_14, label='Fe XIV')
    plt.plot(times_frac, fe_16, label='Fe XVI')
    plt.plot(times_frac, si_12, label='Si XII')
    plt.plot(times_frac, o_6, label='O VI')
    plt.legend()
    plt.title(file_[:-4])
    plt.xlabel('Time [s]')
    plt.ylabel('Ion fraction')
    plt.savefig('/data/khnum/REU2018/jwaczak/data/figs/simulation/simIonFracs/'+file_[:-4]+'.png') 
    plt.close()
