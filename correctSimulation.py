import AIA_tools as aia
import numpy as np
from glob import glob

r0 = np.loadtxt('/data/khnum/REU2018/jwaczak/data/initialRadii.txt', delimiter=',')
dx = r0[0]/20


simFiles = glob('/data/khnum/REU2018/jwaczak/data/simOutput/*.txt')

for sim in simFiles:
    minusExtension = sim[:-4]
    split = minusExtension.split('__')
    print(split)
    t0 = split[0][-8:]
    t1 = split[1][4:]
    n = split[2][3:]

    print('t0: {}, t1: {}, n: {}'.format(t0, t1, n))

    corrected_data, fileName = aia.simulation.applySphericalCorrection(sim, dx, float(t0), float(t1), float(n))
    
