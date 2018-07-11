#!usr/bin/python3


'''
Name: John Waczak
Date: 2018/06/14
Purpose: Simple text parser to make it easy to run application one
         of the time dependent fortran NEI code. Type --help for usage
         details 
'''



import os, sys
import argparse
import numpy as np



# -----------------------------------------------------------------------------------
# set up parser for getting command line arguments
# -----------------------------------------------------------------------------------
parser = argparse.ArgumentParser()

# starting temperature
parser.add_argument("--tempStart", type = str, action = "store", default = "10.0**6", help = "Starting temperature for simulation")

# final temperature
parser.add_argument("--tempEnd", type = np.float64, action = "store", default = 10.0**6.8, help = "ending temperature for simulation")

# density
parser.add_argument("--rho", type = str, action = "store", default = "10.0**7.0", help= "starting density for simulation") 

# num elements
parser.add_argument("--num", type = int, action = "store", default = 2, help = "number of total elements")

# atomic indices
parser.add_argument("--indices", type = str, action = "store", default = "2, 26", help = "atomic indices for desired elements. Separate by single comma. Ex: 2, 26 ")

# output file name
parser.add_argument("--filename", type=str, action = "store", default = "test_onestep_ei.dat", help="Specify the filename for data file output")

# number of time steps
parser.add_argument("--ntime", type=int, action = "store", default = 10, help="Specify the number of time steps")

# dt
parser.add_argument("--dt", type=np.float64, action="store", default=0.75, help="Specify the size of each time step")

args = parser.parse_args()


pathToFiles = "/data/khnum/REU2018/jwaczak/time_dependent_fortran/Applications/NEI_Onestep"


# -----------------------------------------------------------------------------------
# Now we want to open the files and change the lines
# -----------------------------------------------------------------------------------


# open the main.f90 file 
mainFile = open('{}/main.f90'.format(pathToFiles), 'r')
file_lines = mainFile.read().splitlines()

# find parameters and change them 
for i in range(len(file_lines)):
    if ("te_start = ") in file_lines[i]:
        file_lines[i] = "  te_start = {}\t! (K) used to set the initial charge state".format(args.tempStart)
    elif ("te_end = ") in file_lines[i]:
        file_lines[i] = "  te_end = {}\t! (K)".format(args.tempEnd)
    elif ("rhone = ") in file_lines[i]:
        file_lines[i] = "  rhone = {}\t! (cm^-3)".format(args.rho)
    elif ("open(12, file=") in file_lines[i]:
        file_lines[i] = "open(12, file='{}', form=\'unformatted\')".format(args.filename)
    elif ("integer, parameter:: ntime=") in file_lines[i]:
        file_lines[i] = "integer, parameter:: ntime={}".format(args.ntime)
    elif ("real*8, parameter:: dt0 =") in file_lines[i]:
        file_lines[i] = "real*8, parameter:: dt0 = {}".format(str(args.dt))

mainFile.close() 

# re open the file in write mode and update the file 
mainFile = open('{}/main.f90'.format(pathToFiles), 'w')
mainFile.write('\n'.join(file_lines))
mainFile.close()

# open the init.txt file and change the input parameters.
inputFile = open('{}/input.txt'.format(pathToFiles), 'r')
file_lines = inputFile.read().splitlines()



# the lines we need to change are at index 1 and 2.
file_lines[1] = str(args.num)
file_lines[2] = args.indices
file_lines[4] = "\'/data/khnum/REU2018/jwaczak/time_dependent_fortran/eigendata/chianti_8_07/\'"
inputFile.close()


# re open file in write mode and update the file
inputFile = open('{}/input.txt'.format(pathToFiles), 'w')
inputFile.write('\n'.join(file_lines))
inputFile.close()


# -----------------------------------------------------------------------------------
# Now we want to recompile using 'make all'
# -----------------------------------------------------------------------------------
os.chdir(pathToFiles)
os.system("make all")
os.system("./main")
os.chdir("../..")
os.system("mv {}/{} /data/khnum/REU2018/jwaczak/time_dependent_fortran/runs/ ".format(pathToFiles, args.filename))
























