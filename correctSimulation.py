import AIA_tools as aia
import numpy as np
from glob import glob


simFiles = glob('/home/john/gitRepos/REU/jwaczak/data/simOutput/*.txt')

for sim in simFiles:
    minusExtension = sim[:-4]
    split = minusExtension.split('__')
    print(split)
    T0 = split[0][-8:]
    T1 = split[1][4:]
    n = split[2][3:]

    print('T0: {}, T1: {}, n: {}'.format(T0, T1, n))

    corrected_data, fileName = aia.simulation.applySphericalCorrection2(sim, 30, T0, T1, n)
    np.savetxt(fileName , corrected_data, delimiter=',')
