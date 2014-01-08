#!/usr/bin/env python 

import sys
import subprocess
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlb
from pylab import meshgrid
from mpl_toolkits.mplot3d import Axes3D

if len(sys.argv) != 2:
    print "usage: "+sys.argv[0]+" NumPoints"
    exit(1)

filename = sys.argv[1]

subprocess.call(['./multinest', filename])

print 'plotting likelihood...'

datafile = 'posterior_pdfs.dat'

x, y, logL, post = np.loadtxt(datafile, unpack=True)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(x, y, logL, color='blue', s=1.5)

ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('L/Z')

#ax.set_xlim3d(-2, 2)
#ax.set_ylim3d(0, 2)

try: import RaiseWindow
except: pass
plt.show()
