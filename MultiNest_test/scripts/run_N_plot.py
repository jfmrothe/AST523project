#!/usr/bin/env python 

import sys, os
import subprocess
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlb
from pylab import meshgrid
from mpl_toolkits.mplot3d import Axes3D

if (len(sys.argv) == 1 or len(sys.argv) > 4):
    print "usage: "+sys.argv[0]+" runtime_filename logL|post [run code? (0|1)]"
    exit(1)

os.chdir('..')
z_axis = sys.argv[2]
run = 1

if len(sys.argv) == 4:
    run = int(sys.argv[3])

filename = sys.argv[1]

if run:
    subprocess.call(['./multinest', filename])

filename = 'RunTimeFiles/'+filename

print 'plotting likelihood...'

datafile = 'posterior_pdfs.dat'


theta_min, theta_max = np.loadtxt(filename, skiprows=6, usecols=(3,4), unpack=True)
x, y, logL, post = np.loadtxt(datafile, unpack=True)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
plt.title("Gaussian shells: 2000 active points")

if z_axis == 'logL':
    ax.scatter(x, y, logL, color='blue', s=1.5)
    ax.set_zlabel('log(L)')
elif z_axis == 'post':
    ax.scatter(x, y, post, color='blue', s=1.5)
    ax.set_zlabel('L')
else:
    "nothing to plot"
    exit(1)

ax.set_xlabel('x')
ax.set_ylabel('y')

ax.set_xlim3d(theta_min[0], theta_max[0])
ax.set_ylim3d(theta_min[1], theta_max[1])

try: import RaiseWindow
except: pass
plt.show()
