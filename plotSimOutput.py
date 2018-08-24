from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from glob import glob
import os
import AIA_tools as aia
import matplotlib

pathToData = '/home/john/gitRepos/REU/jwaczak/data/simOutput/'
pathToCorrectedData = '/home/john/gitRepos/REU/jwaczak/data/correctedSyntheticObservations/'

pathToFigFolder = '/home/john/gitRepos/REU/jwaczak/data/figs/simulation/'

dataFiles = glob(pathToData+'*.txt')
correctedDataFiles = glob(pathToCorrectedData+'*.txt')


obs_data = np.loadtxt('/home/john/gitRepos/REU/jwaczak/data/observationData/shockData.txt', delimiter=',')

for dataFile in dataFiles:
    for correctedDataFile in correctedDataFiles:
        if dataFile[-43:] == correctedDataFile[-43:]:
            print(dataFile)

            # grab the data

            data = np.loadtxt(dataFile, delimiter=',')
            correctedData = np.loadtxt(correctedDataFile, delimiter=',')



            # define the data 

            o_times = obs_data[:,0]
            o_171 = obs_data[:,1]
            o_193 = obs_data[:,2]
            o_211 = obs_data[:,3]
            o_304 = obs_data[:,4]
            o_335 = obs_data[:,5]

            d_times = data[:,0]
            d_171 = data[:,1]
            d_193 = data[:,2]
            d_211 = data[:,3]
            d_304 = data[:,4]
            d_335 = data[:,5]

            c_times = correctedData[:,0]
            c_171 = correctedData[:,1]
            c_193 = correctedData[:,2]
            c_211 = correctedData[:,3]
            c_304 = correctedData[:,4]
            c_335 = correctedData[:,5]



            # calculate ratios
            o_171_193 = np.divide(o_171, o_193)
            o_211_193 = np.divide(o_211, o_193)
            o_304_193 = np.divide(o_304, o_193)
            o_335_193 = np.divide(o_335, o_193) 

            d_171_193 = np.divide(d_171, d_193)
            d_211_193 = np.divide(d_211, d_193)
            d_304_193 = np.divide(d_304, d_193)
            d_335_193 = np.divide(d_335, d_193) 

            c_171_193 = np.divide(c_171, c_193)
            c_211_193 = np.divide(c_211, c_193)
            c_304_193 = np.divide(c_304, c_193)
            c_335_193 = np.divide(c_335, c_193)

            with matplotlib.pyplot.style.context(("dark_background")):

                #---------INTENSITY RATIOS--------------#

                fig, ax = plt.subplots(nrows=2, ncols=2)

                ax[0,0].plot(o_times, o_171_193, 'b--', label='observation')
                ax[0,0].plot(d_times, d_171_193, 'b', label='simulation')
                ax[0,0].plot(c_times, c_171_193, 'r', label='correction')

                ax[0,1].plot(o_times, o_211_193, 'b--', label='observation')
                ax[0,1].plot(d_times, d_211_193, 'b', label='simulation')
                ax[0,1].plot(c_times, c_211_193, 'r', label='correction')

                ax[1,0].plot(o_times, o_304_193, 'b--', label='observation')
                ax[1,0].plot(d_times, d_304_193, 'b', label='simulation')
                ax[1,0].plot(c_times, c_304_193, 'r', label='correction')

                ax[1,1].plot(o_times, o_335_193, 'b--', label='observation')
                ax[1,1].plot(d_times, d_335_193, 'b', label='simulation')
                ax[1,1].plot(c_times, c_335_193, 'r', label='correction')


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

                fig.suptitle(dataFile[-43:-4])
                plt.tight_layout()
                plt.subplots_adjust(left=None, bottom=None, right=None, top=0.85,
                                    wspace=None, hspace=None)
                fileName = pathToFigFolder+'intensityRatios/'+dataFile[-43:-4]
                plt.savefig(fileName+'.png', transparent=True)
                plt.close()



                #---------TOTAL INTENSITY-----------#

                fig_tot, ax_tot = plt.subplots(nrows=2, ncols=3)

                # ax_tot[0,0].plot(o_times, o_171, 'r--', label='observation')
                # ax_tot[0,0].plot(d_times, d_171, 'b', label='simulation')
                ax_tot[0,0].plot(c_times, c_171, 'r', label='correction')

                # ax_tot[0,1].plot(o_times, o_193, 'r--', label='observation')
                # ax_tot[0,1].plot(d_times, d_193, 'b', label='simulation')
                ax_tot[0,1].plot(c_times, c_193, 'r', label='correction')

                # ax_tot[0,2].plot(o_times, o_211, 'r--', label='observation')
                # ax_tot[0,2].plot(d_times, d_211, 'b', label='simulation')
                ax_tot[0,2].plot(c_times, c_211, 'r', label='correction')

                # ax_tot[1,0].plot(o_times, o_304, 'r--', label='observation')
                # ax_tot[1,0].plot(d_times, d_304, 'b', label='simulation')
                ax_tot[1,0].plot(c_times, c_304, 'r', label='correction')

                # ax_tot[1,1].plot(o_times, o_335, 'r--', label='observation')
                # ax_tot[1,1].plot(d_times, d_335, 'b', label='simulation')
                ax_tot[1,1].plot(c_times, c_335, 'r', label='correction')


                ax_tot[0,0].legend(frameon=False)
                ax_tot[0,0].set_ylabel("Total intensity [Dn s^-1 px^-1]")
                ax_tot[0,0].set_xlabel("Time [s]")
                ax_tot[0,0].set_title('171 $\mathrm{\AA}$')

                ax_tot[0,1].legend(frameon=False)
                ax_tot[0,1].set_ylabel("Total intensity [Dn s^-1 px^-1]")
                ax_tot[0,1].set_xlabel("Time [s]")
                ax_tot[0,1].set_title('193 $\mathrm{\AA}$')


                ax_tot[0,2].legend(frameon=False)
                ax_tot[0,2].set_ylabel("Total intensity [Dn s^-1 px^-1]")
                ax_tot[0,2].set_xlabel("Time [s]")
                ax_tot[0,2].set_title('211 $\mathrm{\AA}$')


                ax_tot[1,0].legend(frameon=False)
                ax_tot[1,0].set_ylabel("Total intensity [Dn s^-1 px^-1]")
                ax_tot[1,0].set_xlabel("Time [s]")
                ax_tot[1,0].set_title('304 $\mathrm{\AA}$')

                ax_tot[1,1].legend(frameon=False)
                ax_tot[1,1].set_ylabel("Total intensity [Dn s^-1 px^-1]")
                ax_tot[1,1].set_xlabel("Time [s]")
                ax_tot[1,1].set_title('335 $\mathrm{\AA}$')

                fileName = pathToFigFolder+'totalIntensity/'+dataFile[-43:-4]
                plt.tight_layout()
                plt.savefig(fileName+'.png', transparent=True)
                plt.close()

