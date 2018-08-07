import numpy as np
import sunpy.map as smap
import os
from astropy.coordinates import SkyCoord
from astropy.visualization.mpl_normalize import ImageNormalize  
import astropy.units as u
from glob import glob
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from scipy.interpolate import interp1d
from astropy.convolution import convolve, Box1DKernel
from datetime import datetime, timedelta 



def getScaleFactor(data_map):
    x, y = np.meshgrid(*[np.arange(v.value) for v in data_map.dimensions]) * u.pix

    hpc_coords = data_map.pixel_to_world(x, y)
    r = np.sqrt(hpc_coords.Tx ** 2 + hpc_coords.Ty ** 2) / data_map.rsun_obs

    rsun_step_size = 0.01
    rsun_array = np.arange(1, r.max(), rsun_step_size)
    y = np.array([data_map.data[(r > this_r) * (r < this_r + rsun_step_size)].mean()
              for this_r in rsun_array])

    params = np.polyfit(rsun_array[rsun_array < 1.5],
                    np.log(y[rsun_array < 1.5]), 1)


    scale_factor = np.exp((r-1)*-params[0])
    scale_factor[r < 1] = 1
    return scale_factor


def radialFilter(data_map, scale_factor):
    """
    apply lograithmic radial filter to data in a fits file in order to enhance
    off limb features 
    """
    scaled_map = smap.Map(data_map.data * scale_factor, data_map.meta)
    scaled_map.plot_settings['norm'] = ImageNormalize(stretch=data_map.plot_settings['norm'].stretch, vmin=data_map.data.min(), vmax=data_map.data.max())
    return scaled_map


def makeSubmap(data_map, bl_x, bl_y, tr_x, tr_y):
    """
    Make a submap from a larger map object via specified top right and bottom left positions
    """

    bl = SkyCoord(bl_x*u.arcsec, bl_y*u.arcsec, frame = data_map.coordinate_frame)
    tr = SkyCoord(tr_x*u.arcsec, tr_y*u.arcsec, frame = data_map.coordinate_frame)
    return data_map.submap(bl, tr)


def animate(maps, xlabel, ylabel, path, blackWhite = False, movieName = 'movie'):
    """
    Create series of figures and turn them into a movie using ffmpeg 
    """

    for i in range(len(maps)):
        fig, ax = plt.subplots(figsize=(10,10))

        if blackWhite:
            maps[i].plot(cmap=plt.cm.Greys)
        else:
            maps[i].plot()

        plt.savefig(path+"/{0:04d}".format(i) )
        plt.close()

    os.system("ffmpeg -framerate 20 -qscale 1 -i {}/%04d.png {}.mp4".format(path, movieName))



def getRadialSlit(x_i, x_f, y_i, y_f, resolution, data):
    """
    description: Take coordinates for two points and return a SkyCoord object containing the coordinates for
                 a line that can be ploted on top of a smap image 
    """
    t = np.linspace(0, 1, resolution)
    radialSlit = SkyCoord(x_i+t*(x_f-x_i), y_i+t*(y_f-y_i), frame=data.coordinate_frame)
    return radialSlit


def getRadialSlitIntensity(radialSlit, data):
    """
    description: Take line from getRadialSlit and recover the intensity data for analysis
    """
    pixels = np.asarray(np.rint(data.world_to_pixel(radialSlit)), dtype=int)
    x = pixels[0, :]
    y = pixels[1, :]
    intensity_along_slit = data.data[y,x]

    #length_in_pixels = len(radialSlit)
    lengths_in_pixels=[0.0]
    for i in range(1,len(x)):
        lengths_in_pixels.append(np.sqrt((x[i]-x[0])**2+(y[i]-y[0])**2))
    lengths_in_pixels=np.asarray(lengths_in_pixels)
    
    arcsec_per_pixel = data.meta['cdelt1']
    arcsec_per_solarRadii = data.meta['rsun_obs']
    solarRadii_per_pixel = arcsec_per_pixel / arcsec_per_solarRadii

    distances = solarRadii_per_pixel*lengths_in_pixels

    return intensity_along_slit, distances


def smoothAndNormalizeSlit(raw_slit_intensity):
    smoothed =  np.asarray(convolve(raw_slit_intensity, Box1DKernel(15)))
    #normalized = [elem/smoothed.max() for elem in smoothed]
    #return normalized
    return smoothed

def getTimes(data_maps):
    times = []
    for map_ in data_maps:
        time = map_.meta['date_obs']
        time = datetime.strptime(time, '%Y-%m-%dT%H:%M:%S.%f')
        times.append(time)

    return times 


def plotRadialStack(times, distances, intensities, fig, ax):
    assert(np.shape(distances) == np.shape(intensities))
    shape = np.shape(distances)
    ts = [[times[i] for i in range(shape[0])] for j in range(shape[1])]

    ax.pcolormesh(ts, np.transpose(distances), np.transpose(intensities)) 
    fig.autofmt_xdate()
    myFmt = mdates.DateFormatter('%H:%M:%S')
    ax.xaxis.set_major_formatter(myFmt)



def getNearestValue(array, value):
    """ Return the index of the element in the array closest to value"""
    return np.abs(array-value).argmin() 
















