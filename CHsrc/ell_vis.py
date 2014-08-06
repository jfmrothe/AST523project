import matplotlib
from matplotlib import pyplot
import numpy as np

class clustering(object):
    def __init__(self,points,centers,covMats,enlFacs):
        self.points = points
        self.centers = centers
        self.covMats = covMats
        self.enlFacs = enlFacs


def lineEllipse(center,covMat,enlFac):
    U, s , Vh = np.linalg.svd(covMat)
    x0,y0 = center
    phi = np.math.atan2(U[1,0],U[0,0])
    a=np.math.sqrt(s[0])*np.sqrt(enlFac)
    b=np.math.sqrt(s[1])*np.sqrt(enlFac)
    Nsteps = 25
    line = [(x0+a*np.cos(alpha)*np.cos(phi)-b*np.sin(alpha)*np.sin(phi),y0+a*np.cos(alpha)*np.sin(phi)+b*np.sin(alpha)*np.cos(phi)) for alpha in np.arange(0,2*np.pi*(1+1./Nsteps),2*np.pi/Nsteps)]
    pyplot.plot(*zip(*line),color='b',linestyle="-")
    return

def readClustering(infile,Npts):
    #expected format: for each iteration, one empty line, Npts points, then an arbitrary no of ellipsoids (enlfac,center,covmat), finally an empty line
    data = []
    f = file(infile,'r')
    for line in f:
        line = line.strip()
        if line.startswith("#"):
            continue
        data.append(line)

    starts = [i for i, j in enumerate(data) if j == ""]
    clusterings = []
    for i,start in enumerate(starts[:-1]):
        points = [[float(el) for el in x.split()] for x in data[start+1:start+Npts+1]]
        enlFacs = [float(x) for x in data[start+Npts+1:starts[i+1]-1:4]]
        centers = [[float(el) for el in x.split()] for x in data[start+Npts+2:starts[i+1]:4]]
        covMats = [[[float(el) for el in data[i].split()],[float(el) for el in data[i+1].split()]] for i in range(start+Npts+3,starts[i+1],4)]
        clusterings.append(clustering(points,centers,covMats,enlFacs))     
    return clusterings

def plotCluster(cluster):
    for i in range(len(cluster.enlFacs)):
        lineEllipse(cluster.centers[i],cluster.covMats[i],cluster.enlFacs[i])
    return

clusters = readClustering("test/egg2000.txt",2000)
print len(clusters)
no = 305

fig = pyplot.figure()
plotCluster(clusters[no])
pyplot.scatter(*zip(*clusters[no].points),marker='+',
     color='black', s=5)
pyplot.show()
