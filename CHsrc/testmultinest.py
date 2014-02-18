#!/usr/bin/env python
#os.chdir("oblateTransit")
#os.system("make clean")
#os.system("make all")
#os.chdir("..")
import sys
sys.path.append("oblateTransit")
#import transitmodel as model #module that call occultquad, and also hold the data
from lighthouse import lighthousemodel as Model #module that call occultquad, and also hold the data
#from eggbox import eggboxmodel as Model #module that call occultquad, and also hold the data
#from gaussianshell import guaussianshellmodel as Model #module that call occultquad, and also hold the data
#from occultquad import occultquad
from Point import Point
from Ellipsoid import Ellipsoid
from Samplers import Samplers as MNest
import numpy as np
def main():
    model = Model('example.cfg')
    minvals,maxvals,e,Np = model.Getinitial()
    #print minvals,maxvals,e
    sampler = MNest(minvals,maxvals,e,Np, "uniform")
    #what I intended to use, for sampler does not need to know the names
    #sampler = MNest(priortypes,minvals,maxvals,guessvals) 
    print "running MultiNest algorithm... this may take a few minutes\n"
    nest = 0
    Flag = True
    NumRecluster = 0
    Alltheta = np.zeros([model.D,int(Np)])
    theta = np.zeros(model.D)
    logL = np.zeros(Np)
    logLmax = 0.
    logLmin = 0.
    THRESH  = 1.0e-5
    sampler.getAlltheta(Alltheta)
    model.Get_L(Alltheta.ravel(),logL,Np)
    sampler.SetAllPoint(logL) 
    #print Alltheta[0,0],Alltheta[1,0]
    #Flag = False
    while Flag: 
        X_i = np.exp(-1.0*nest/(1.0*Np))
        #print Alltheta.shape
        #print 'after Alltheta'
        #print Alltheta
        #discard and resample, get logLmax for convergence check as byproduct
        sampler.DisgardWorstPoint(int(nest)) 
        logLmax = sampler.getlogLmax() 
        logLmin = sampler.getlogLmin()
        #print logLmax,logLmin
        FlagSample = True
        templogL = np.array([0.])
        #print 'before reset'
        while FlagSample:
            sampler.ResetWorstPoint(theta) 
            model.Get_L(theta,templogL,1)
            #print theta, templogL,logLmin
            FlagSample=(logLmin > templogL[0])
        temp = templogL[0]
        sampler.ResetWorstPointLogL(temp) 
        #print 'after reset'
        
        #ellipsoidal partitioning 
        # I will check the criterion inside the function Recluster. I feel the nest==0 exception is unnecessary, reclustering should only happen if quality is bad. talk to chelsea
        #if(nest==0 or sampler.ClusteringQuality(X_i) > model.repartition): 
        # recluster?
        #print 'before recluster'
        NumRecluster += sampler.Recluster(X_i,model.repartition)
        #print 'after recluster'
        #NumRecluster+=1
        nest+=1 
        #print nest,X_i
        #Flag = False
        Flag = THRESH < abs(X_i*logLmax)
    #print NumRecluster
    #print 'before output'
    #output
    ntotal = sampler.countTotal()
    posterior = np.zeros([model.D,ntotal])
    prob = np.zeros(ntotal)
    sampler.getPosterior(posterior,prob)
    print posterior.shape
    logzinfo = np.zeros(3)
    sampler.getlogZ(logzinfo)
    model.Output(posterior.ravel(),prob)
    print "information: H=%f bits" % logzinfo[0]
    print "global evidence: logZ = %f +/- %f" % (logzinfo[1],logzinfo[2])
    return

def testL():
    model = Model('example.cfg')
    minvals,maxvals,e,Np = model.Getinitial()
    #x = 4.0*np.random.randn(10000)-(minvals[0]+maxvals[0])/2.
    #x = 16.0*np.random.randn(10000)+(minvals[0]+maxvals[0])/2.
    #y = 16.0*np.random.randn(10000)+(minvals[0]+maxvals[0])/2.
    x = 6.0*np.random.randn(10000)+(minvals[0]+maxvals[0])/2.
    y = 6.0*np.random.randn(10000)+(minvals[0]+maxvals[0])/2.
    #y = np.random.randn(10000)*4.0
    #lighthouse
    if(0):
        indexa = y>0.5
        indexb = y<1.5
        indexc = x>-2
        indexd = x<2
    #eggbox
    if(0):
        indexa = y>0
        indexb = y<31.415
        indexc = x>0
        indexd = x<31.415
    if(1):
        indexa = y>-6
        indexb = y<6
        indexc = x>-6
        indexd = x<6

    y = y[indexa*indexb*indexc*indexd]
    x = x [indexa*indexb*indexc*indexd]
    logL = np.zeros(len(y))
    theta = np.vstack((x,y)).ravel(order='F')
    #print theta.shape
    #.ravel()
    #print theta[0],theta[1],x[0],y[0]
    #return
    #print theta.shape
    #print x
    #return
    model.Get_L(theta,logL,len(y))
    for i in xrange(len(y)):
        print x[i],y[i],logL[i]


if __name__=='__main__':
    main()
    #testL()
    
