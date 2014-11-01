#!/usr/bin/env python
import matplotlib 
from matplotlib import pyplot as plt
import matplotlib.mlab as mlb
from pylab import meshgrid
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
def plotlighthouse():
    infile = 'lighthouse.tab'
    #infile = 'temp.txt'
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    a = 305.0
    e = 65.
    ax.view_init(azim=a, elev=e)
     #x,y,prob = np.loadtxt(infile,usecols=[0,1,2],unpack=True)
    ##ax.scatter(x,y,np.exp(prob+160.455184),color='blue',s=1.5)
    #ax.scatter(x,y,np.exp(prob+235.83),color='blue',s=1.5)
    #infile = 'data/gaussianshell_norescale.tab'
    #infile = 'eggbox.tab'
    #infile = 'gaussianshell.tab'
    #x,y,prob = np.loadtxt(infile,usecols=[0,1,2],unpack=True)
    #ax.scatter(x,y,np.exp(prob+160.455184),color='red',s=1.5)
    #ax.scatter(x,y,np.exp(prob+235.83),color='red',s=1.5)
    #infile = 'data/gaussianshell_rescale1.tab'
    #infile = 'toygauss.tab'
    x,y,prob = np.loadtxt(infile,usecols=[0,1,2],unpack=True)
    pmin = -200.
    pmax = +300.
#    ax.scatter(x,y,np.exp(prob+160.455184),color='blue',s=1.5)
    #ax.scatter(x,y,[np.exp(min(pmax,max(p,pmin))) for p in prob],color='blue',s=1.5)
    ax.scatter(x,y,[min(pmax,max(p,pmin)) for p in prob],color='blue',s=1.5)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("L")
    plt.show()
    return

if __name__=='__main__':
    plotlighthouse()
