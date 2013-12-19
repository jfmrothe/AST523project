#!/usr/bin/env python 

import subprocess
import numpy as np
import matplotlib.pyplot as plt
from pylab import meshgrid
from mpl_toolkits.mplot3d import Axes3D

subprocess.call('./drive')

print 'plotting joint posterior...'

datafile = 'posterior_pdfs.dat'

x, y, post = np.loadtxt(datafile, unpack=True)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(x, y, post, color='blue', s=1.5)

ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('posterior')

ax.set_xlim3d(-2, 2)
ax.set_ylim3d(0, 2)

try: import RaiseWindow
except: pass
plt.show()
