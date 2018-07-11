from __future__ import division, print_function
from sunpy.net import Fido, attrs
from datetime import datetime, timedelta
import astropy.units as u
import os
from glob import glob


def download(dtime, dt, dir_path, *args):
    """
    Description:
              download fits files for specified time window and wavelength channels
              from the VSO


    dtime:    should be a datetime object. EX: datetime(2010, 6, 13, 5, 25, 0)
    dt:       should be a time delta object in correct units. EX: timedelta(seconds=6)
    dir_path: path to directory for storing fits files
    *argv:    wavelength channels you want fits files for. EX: 193, 211, etc... 
    """


    sform = '%Y/%m/%d %H:%M:%S'  # format datetime objects correctly

    start = datetime.strftime(dtime-dt, sform)
    end = datetime.strftime(dtime+dt, sform)

    time = attrs.Time(start, end)
    ins = attrs.Instrument('aia') 


    for thing in args:
        print("Downloading {} angstrom data\n".format(thing)) 
        wave = attrs.Wavelength(int(thing)*u.AA)
        searchResults = Fido.search(time, ins, wave)

        for file_ in searchResults:
            print("{}".format(file_))

        dl_fil = Fido.fetch(searchResults, path=dir_path, wait=True) 



