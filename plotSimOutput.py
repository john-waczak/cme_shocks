import numpy as np
import matplotlib.pyplot as plt
from glob import glob
import os 

pathToData = '/data/khnum/REU2018/jwaczak/data/simOutput/'
pathToFigFolder = '/data/khnum/REU2018/jwaczak/data/figs/simulation/'

os.chdir(pathToData) 
dataFiles = glob('*.*')
print(dataFiles)



for file_ in dataFiles:
    data = np.loadtxt(file_, delimiter=',')

    fig, ax = plt.subplots(nrows=2, ncols=2)
    ax[0,0].plot(data[:,0], np.divide(data[:,1], data[:,2]), label='$171/193$') 
    ax[0,1].plot(data[:,0], np.divide(data[:,3], data[:,2]), label='$211/193$')
    ax[1,0].plot(data[:,0], np.divide(data[:,4], data[:,2]), label='$304/193$')
    ax[1,1].plot(data[:,0], np.divide(data[:,5], data[:,2]), label='$335/193$') 

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
    print(fileName)
    plt.savefig(fileName+'.png')

    fig_tot, ax_tot = plt.subplots()
    ax_tot.plot(data[:,0], data[:,1]/data[:,1].max(), label='$171$ $\AA$') 
    ax_tot.plot(data[:,0], data[:,2]/data[:,2].max(), label='$193$ $\AA$')
    ax_tot.plot(data[:,0], data[:,3]/data[:,3].max(), label='$211$ $\AA$')
    ax_tot.plot(data[:,0], data[:,4]/data[:,4].max(), label='$304$ $\AA$')
    ax_tot.plot(data[:,0], data[:,5]/data[:,5].max(), label='$335$ $\AA$') 
    #ax_tot.set_ylabel("DN $\mathrm{s}^{-1} \mathrm{pix}^{-1}$", fontsize=18)
    ax_tot.set_ylabel("Tot Intensity (normalized to 1) [FIXME: Units]")
    ax_tot.set_xlabel("Time [s]", fontsize=18)
    ax_tot.legend()
    ax_tot.set_title(file_[:-4])

    fileName = pathToFigFolder+'/totalIntensity/'+file_[:-4]
    print(fileName)
    plt.savefig(fileName+'.png')






