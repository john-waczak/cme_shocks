from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from glob import glob
import os
import AIA_tools as aia
import matplotlib

pathToData = '/data/khnum/REU2018/jwaczak/data/simOutput/'
pathToFigFolder = '/data/khnum/REU2018/jwaczak/data/figs/simulation/data/'
os.chdir(pathToData)
dataFiles = glob('*.txt')

print(dataFiles)

obs_data = np.loadtxt('/data/khnum/REU2018/jwaczak/data/observationData/shockData.txt', delimiter=',')

for file_ in dataFiles:
    data = np.loadtxt(file_, delimiter=',')

    fig, ax = plt.subplots(nrows=2, ncols=2)

    data[:,0] = data[:,0]
    d_171 = data[:,1]
    d_193 = data[:,2]
    d_211 = data[:,3]
    d_304 = data[:,4]
    d_335 = data[:,5]

    r_171_193 = np.divide(d_171, d_193)
    r_211_193 = np.divide(d_211, d_193)
    r_304_193 = np.divide(d_304, d_193)
    r_335_193 = np.divide(d_335, d_193) 

    with matplotlib.pyplot.style.context(("dark_background")):
        for i in range(len(d_171)):
            print(r_171_193[i], d_171[i]/d_193[i])

            ax[0,0].plot(data[:,0], r_171_193, 'b', label='simulation')
            ax[0,0].plot(obs_data[:,0], np.divide(obs_data[:,1], obs_data[:,2]), 'r', label='observation')

            ax[0,1].plot(data[:,0], r_211_193, 'b', label='simulation')
            ax[0,1].plot(obs_data[:,0], np.divide(obs_data[:,3], obs_data[:,2]), 'r', label='observation')

            ax[1,0].plot(data[:,0], r_304_193, 'b', label='simulation')
            ax[1,0].plot(obs_data[:,0], np.divide(obs_data[:,4], obs_data[:,2]), 'r', label='observation')

            ax[1,1].plot(data[:,0], r_335_193, 'b', label='simulation')
            ax[1,1].plot(obs_data[:,0], np.divide(obs_data[:,5], obs_data[:,2]), 'r', label='observation')

            ax[0,0].legend(frameon=False) 
            ax[0,1].legend(frameon=False)
            ax[1,0].legend(frameon=False)
            ax[1,1].legend(frameon=False)

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
     #   plt.tight_layout()
        plt.subplots_adjust(left=None, bottom=None, right=None, top=0.85,
                            wspace=None, hspace=None)
        fileName = pathToFigFolder+'/intensityRatios/'+file_[:-4]
        plt.savefig(fileName+'.png', transparent=True)
        plt.close()

        fig_tot, ax_tot = plt.subplots()
        ax_tot.plot(data[:,0], data[:,1]/data[:,1].max(), 'c', label='$171$ $\AA$') 
        ax_tot.plot(obs_data[:,0], obs_data[:,1]/obs_data[:,1].max(), 'c--')

        ax_tot.plot(data[:,0], data[:,2]/data[:,2].max(), 'b', label='$193$ $\AA$')
        ax_tot.plot(obs_data[:,0], obs_data[:,2]/obs_data[:,2].max(), 'b--')

        ax_tot.plot(data[:,0], data[:,3]/data[:,3].max(), 'r',label='$211$ $\AA$')
        ax_tot.plot(obs_data[:,0], obs_data[:,3]/obs_data[:,3].max(), 'r--') 

        ax_tot.plot(data[:,0], data[:,4]/data[:,4].max(), 'g', label='$304$ $\AA$')
        ax_tot.plot(obs_data[:,0], obs_data[:,4]/obs_data[:,4].max(), 'g--')

        ax_tot.plot(data[:,0], data[:,5]/data[:,5].max(), 'm', label='$335$ $\AA$')
        ax_tot.plot(obs_data[:,0], obs_data[:,5]/obs_data[:,5].max(), 'm--')

        ax_tot.set_ylabel("Tot Intensity (normalized to 1) ")
        ax_tot.set_xlabel("Time [s]", fontsize=18)

        ax_tot.legend(frameon=False)
        ax_tot.set_title(file_[:-4])

        fileName = pathToFigFolder+'/totalIntensity/'+file_[:-4]
        plt.savefig(fileName+'.png', transparent=True)
        plt.close()








pathToData = '/data/khnum/REU2018/jwaczak/data/correctedSyntheticObservations'
pathToFigFolder = '/data/khnum/REU2018/jwaczak/data/figs/simulation/correctedData/'
os.chdir(pathToData)
dataFiles = glob('*.txt')

print(dataFiles)


for file_ in dataFiles:
    data = np.loadtxt(file_, delimiter=',')

    fig, ax = plt.subplots(nrows=2, ncols=2)

    data[:,0] = data[:,0]
    d_171 = data[:,1]
    d_193 = data[:,2]
    d_211 = data[:,3]
    d_304 = data[:,4]
    d_335 = data[:,5]

    r_171_193 = np.divide(d_171, d_193)
    r_211_193 = np.divide(d_211, d_193)
    r_304_193 = np.divide(d_304, d_193)
    r_335_193 = np.divide(d_335, d_193) 

    for i in range(len(d_171)):
        print(r_171_193[i], d_171[i]/d_193[i])

    ax[0,0].plot(data[:,0], r_171_193, 'k', label='simulation')
    ax[0,0].plot(obs_data[:,0], np.divide(obs_data[:,1], obs_data[:,2]), 'b', label='observation')

    ax[0,1].plot(data[:,0], r_211_193, 'k', label='simulation')
    ax[0,1].plot(obs_data[:,0], np.divide(obs_data[:,3], obs_data[:,2]), 'b', label='observation')

    ax[1,0].plot(data[:,0], r_304_193, 'k', label='simulation')
    ax[1,0].plot(obs_data[:,0], np.divide(obs_data[:,4], obs_data[:,2]), 'b', label='observation')

    ax[1,1].plot(data[:,0], r_335_193, 'k', label='simulation')
    ax[1,1].plot(obs_data[:,0], np.divide(obs_data[:,5], obs_data[:,2]), 'b', label='observation')

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
    ax_tot.plot(obs_data[:,0], obs_data[:,1]/obs_data[:,1].max(), 'k--')

    ax_tot.plot(data[:,0], data[:,2]/data[:,2].max(), 'b', label='$193$ $\AA$')
    ax_tot.plot(obs_data[:,0], obs_data[:,2]/obs_data[:,2].max(), 'b--')

    ax_tot.plot(data[:,0], data[:,3]/data[:,3].max(), 'r',label='$211$ $\AA$')
    ax_tot.plot(obs_data[:,0], obs_data[:,3]/obs_data[:,3].max(), 'r--') 

    ax_tot.plot(data[:,0], data[:,4]/data[:,4].max(), 'g', label='$304$ $\AA$')
    ax_tot.plot(obs_data[:,0], obs_data[:,4]/obs_data[:,4].max(), 'g--')

    ax_tot.plot(data[:,0], data[:,5]/data[:,5].max(), 'm', label='$335$ $\AA$')
    ax_tot.plot(obs_data[:,0], obs_data[:,5]/obs_data[:,5].max(), 'm--')

    ax_tot.set_ylabel("Tot Intensity (normalized to 1) ")
    ax_tot.set_xlabel("Time [s]", fontsize=18)

    ax_tot.legend()
    ax_tot.set_title(file_[:-4])

    fileName = pathToFigFolder+'/totalIntensity/'+file_[:-4]
    plt.savefig(fileName+'.png')
    plt.close()


