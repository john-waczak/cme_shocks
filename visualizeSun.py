from __future__ import division, print_function

import numpy as np
import sunpy.map as smap
import os
from astropy.coordinates import SkyCoord
from glob import glob
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.visualization.mpl_normalize import ImageNormalize # for the radial filter technique
from AIA_tools import *


# get the data

pathTo193 = "../data/June_2010_CME_burst/193_fits/"
fits_193 = glob(pathTo193 + "*.fits")

pathTo211 = "../data/June_2010_CME_burst/211_fits/"
fits_211 = glob(pathTo211 + "*.fits")

print("Gathering fits files") 

# turn them into sunpy map objects
data_193 = smap.Map(fits_193)
data_211 = smap.Map(fits_211)



#  coordinates for submape box 
bl_x = 600
bl_y = -650
tr_x = 1250
tr_y = 0


#  coordinates for radial slit 
x_i = 990*u.arcsec 
y_i = -450*u.arcsec 
x_f = 1220*u.arcsec 
y_f = -245*u.arcsec 


# create image directories
print("Creating figure directories")
if not os.path.isdir("../data/June_2010_CME_burst/193_figs/"):
    os.mkdir("../data/June_2010_CME_burst/193_figs/")
else:
    os.system("rm -r ../data/June_2010_CME_burst/193_figs/*.*")

if not os.path.isdir("../data/June_2010_CME_burst/211_figs/"):
    os.mkdir("../data/June_2010_CME_burst/211_figs") 
else:
    os.system("rm -r ../data/June_2010_CME_burst/211_figs/*.*")


scale_factor_193 = analysis.getScaleFactor(data_193[0]) 
scale_factor_211 = analysis.getScaleFactor(data_211[0])

# create and save figures 
for i in range(len(data_193)):
    print("{:2.2f}%".format((i/len(data_193))*100))

    print("\tApplying radial filter")
    filtered_193 = analysis.radialFilter(data_193[i], scale_factor_193)
    filtered_211 = analysis.radialFilter(data_211[i], scale_factor_211) 

    print("\tCreating submaps")
    submap_193 = analysis.makeSubmap(filtered_193, bl_x, bl_y, tr_x, tr_y)
    submap_211 = analysis.makeSubmap(filtered_211, bl_x, bl_y, tr_x, tr_y)

    print("\tGenerating radial slits") 
    slit_193 = analysis.getRadialSlit(x_i, x_f, y_i, y_f, 200, submap_193)
    slit_211 = analysis.getRadialSlit(x_i, x_f, y_i, y_f, 200, submap_211)

    print("\tGetting slit intensities and distances")
    slit_intensity_193, distances_193 = analysis.getRadialSlitIntensity(slit_193, submap_193)
    slit_intensity_211, distances_211 = analysis.getRadialSlitIntensity(slit_211, submap_211)

    print("\tSmoothing and normalizing")
    smoothed_193 = analysis.smoothAndNormalizeSlit(slit_intensity_193)
    smoothed_211 = analysis.smoothAndNormalizeSlit(slit_intensity_211)

    print("\tCreating figures")

    # 193 figure
    fig1 = plt.figure(figsize=(15,10)) 
    ax11 = plt.subplot(121, projection=submap_193) 
    submap_193.plot(axes=ax11, cmap=plt.cm.Greys)
    ax11.plot_coord(slit_193, color='c', linewidth=0.5)

    ax12 = plt.subplot(122)  
    ax12.plot(distances_193, smoothed_193, 'k')
    ax12.set_xlabel("Distance [$R_\odot$]")
    ax12.set_ylabel("Normalized Intenisty")
    ax12.set_ylim(0, 1.1)
    ax12.set_xlim(distances_193[0], distances_193[-1])
    ax12.set_aspect(np.abs(distances_193[0]-distances_193[-1])/1.1)
    ax12.set_title("Intensity along slit (w/ radial filter)")
    plt.tight_layout()
    plt.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.9,
                wspace=None, hspace=None)
    plt.savefig("../data/June_2010_CME_burst/193_figs/image_{0:04d}".format(i))
    plt.close(fig1) 

    # 211 figure
    fig2 = plt.figure(figsize=(15,10))
    ax21 = plt.subplot(121, projection=submap_211) 
    submap_211.plot(axes=ax21, cmap=plt.cm.Greys)
    ax21.plot_coord(slit_211, color='c', linewidth=0.5)
    
    ax22 = plt.subplot(122) 
    ax22.plot(distances_211, smoothed_211, 'k')
    ax22.set_xlabel("Distance [$R_\odot$]")
    ax22.set_ylabel("Normalized Intensity")
    ax22.set_ylim (0, 1.1)
    ax22.set_xlim(distances_211[0], distances_211[-1])
    ax22.set_aspect(np.abs(distances_211[0]-distances_211[-1])/1.1)
    ax22.set_title("Intensity along slit (w/ radial fitler)")
    plt.tight_layout()
    plt.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.9,
                wspace=None, hspace=None)
    plt.savefig("../data/June_2010_CME_burst/211_figs/image_{0:04d}".format(i)) 
    plt.close(fig2) 

os.system("ffmpeg -framerate 20 -qscale 1 -i ../data/June_2010_CME_burst/193_figs/image_%04d.png ../data/June_2010_CME_burst/radialSlit_193.mp4 ")
os.system("ffmpeg -framerate 20 -qscale 1 -i ../data/June_2010_CME_burst/211_figs/image_%04d.png ../data/June_2010_CME_burst/radialSlit_211.mp4" )
















