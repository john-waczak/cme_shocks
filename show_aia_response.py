import numpy as np
import matplotlib.pyplot as plt
import scipy.io as sio

filename = 'aia_response.sav'
res = sio.readsav(filename)

logte = res['logte']
a94 = res['a94']
a131 = res['a131']
a171 = res['a171']
a193 = res['a193']
a211 = res['a211']
a304 = res['a304']
a335 = res['a335']

fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(3, 2, 1)
plt.plot(logte, a94, label='94', c='red')
plt.plot(logte, a131, label='131', c='blue')
plt.plot(logte, a171, label='171', c='magenta')
plt.plot(logte, a193, label='193', c='green')
plt.plot(logte, a211, label='211', c='black')
plt.plot(logte, a304, label='304')
plt.plot(logte, a335, label='335')
plt.xlabel('Log$_{10}$ T (K)')
plt.ylabel('Response (DN cm$^5$ s$^{-1}$ pix$^{-1}$)')
plt.yscale("log")
plt.ylim((1.0e-28, 1.0e-23))
plt.xlim((5.5, 7.0))
plt.legend()

ax2 = fig.add_subplot(3, 2, 3)
plt.plot(logte, a171/a193)
plt.xlabel('Log$_{10}$ T (K)')
plt.ylabel('171/193')
plt.xlim((5.5, 7.0))

ax3 = fig.add_subplot(3, 2, 4)
plt.plot(logte, a211/a193)
plt.xlabel('Log$_{10}$ T (K)')
plt.ylabel('211/193')
plt.xlim((5.5, 7.0))

ax4 = fig.add_subplot(3, 2, 5)
plt.plot(logte, a304/a193)
plt.xlabel('Log$_{10}$ T (K)')
plt.ylabel('304/193')
plt.ylim((0, 0.5))
plt.xlim((5.5, 7.0))

ax5 = fig.add_subplot(3, 2, 6)
plt.plot(logte, a335/a193)
plt.xlabel('Log$_{10}$ T (K)')
plt.ylabel('335/193')
plt.xlim((5.5, 7.0))


plt.tight_layout()
fig.savefig('aia_response_and_ratio.png', dpi=150)
