import AIA_tools as aia
import numpy as np
from glob import glob
import argparse


parser = argparse.ArgumentParser(description="calculate 3D model AIA intensities for all present simulation files")
parser.add_argument('-n', dest='N_w', type=int, default=30, help="The number of plasma windows along the line of sight")

args = parser.parse_args()


simFiles = glob('/home/john/gitRepos/REU/jwaczak/data/simOutput/*.txt')

for sim in simFiles:
    minusExtension = sim[:-4]
    split = minusExtension.split('__')
    T0 = split[0][-8:]
    T1 = split[1][4:]
    n = split[2][3:]

    print('T0: {}, T1: {}, n: {}'.format(T0, T1, n))

    corrected_data, fileName = aia.simulation.applySphericalCorrection2(sim, args.N_w, T0, T1, n)
    np.savetxt(fileName , corrected_data, delimiter=',')
