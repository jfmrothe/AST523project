#!/usr/bin/env python
import os
#os.chdir("oblateTransit")
#os.system("make clean")
#os.system("make all")
#os.chdir("..")
import sys
sys.path.append("oblateTransit")
#from transitmodel import transitmodel as Model #module that call occultquad, and also hold the data
#from lighthouse import lighthousemodel as Model #module that call occultquad, and also hold the data
#from eggbox import eggboxmodel as Model #module that call occultquad, and also hold the data
#from toygauss import toygaussmodel as Model
from poly import polymodel as Model #module that call occultquad, and also hold
#from gaussianshell import guaussianshellmodel as Model #module that call occultquad, and also hold the data
#from occultquad import occultquad
from Point import Point
from Ellipsoid import Ellipsoid
from Samplers import Samplers as MNest
import numpy as np
def main():
    #model = Model('toygauss.cfg')
    model = Model('example.cfg')
    #model = Model('oblateTransit/sorted_example.cfg')
    #model = Model('oblateTransit/testmarching_example.cfg')
    minvals,maxvals,e,Np = model.Getinitial()
    #print minvals,maxvals,e
    sampler = MNest(minvals,maxvals,e,Np, "uniform")
    #what I intended to use, for sampler does not need to know the names
    #sampler = MNest(priortypes,minvals,maxvals,guessvals) 
    print "running MultiNest algorithm... this may take a few minutes\n"
    nest = 0
    Flag = True
    NumRecluster = 0
    NumLeval = Np
    Alltheta = np.zeros([model.D,int(Np)])
    theta = np.zeros(model.D)
    logL = np.zeros(Np)*1.0
    logLmax = 0.
    logLmin = 0.
    THRESH  = 1.0e-7
    FullReclusterPeriod = 15000
    sampler.getAlltheta(Alltheta)
    model.Get_L(Alltheta.ravel(),logL,Np)
    sampler.SetAllPoint(logL) 
    #print Alltheta[0,0],Alltheta[1,0]
    #Flag = False
    Debug = False
    count=0
    logzinfo = np.zeros(3)
    count = 0
    while Flag:

        X_i = np.exp(-1.0*nest/(1.0*Np))
        #discard and resample, get logLmax for convergence check as byproduct
        if(Debug):
            print '----------------------'
            ntotal = sampler.countTotal()
            posterior = np.zeros([model.D,ntotal])
            prob = np.zeros(ntotal)
            sampler.getPosterior(posterior,prob)
            print posterior
            print prob
        
        sampler.DisgardWorstPoint(int(nest)) 
        if(Debug): 
            print '----------------------'
            ntotal = sampler.countTotal()
            posterior = np.zeros([model.D,ntotal])
            prob = np.zeros(ntotal)
            sampler.getPosterior(posterior,prob)
            print posterior
            print prob

        logLmax = sampler.getlogLmax() 
        logLmin = sampler.getlogLmin()
        #print logLmax,logLmin
        FlagSample = True
        templogL = np.array([0.])
        #print 'before reset'
        while FlagSample:
            sampler.ResetWorstPoint(theta) 
            model.Get_L(theta,templogL,1)
            NumLeval += 1
            #print theta, templogL,logLmin
            FlagSample=(logLmin > templogL[0])
            #FlagSample=False
        
        temp = templogL[0]
        #print templogL[0]
        sampler.ResetWorstPointLogL(temp) 
        #print 'after reset'
        if (Debug): 
            print '----------------------'
            ntotal = sampler.countTotal()
            posterior = np.zeros([model.D,ntotal])
            prob = np.zeros(ntotal)
            sampler.getPosterior(posterior,prob)
            print posterior
            print prob
        if nest % FullReclusterPeriod == 0:
            NumRecluster += 1
            sampler.FullRecluster(X_i)
        else:
            NumRecluster += sampler.Recluster(X_i,model.repartition)
        #ellipsoidal partitioning 
        #if(nest==0 or sampler.ClusteringQuality(X_i) > model.repartition):
        #    count+=1
        #    # recluster?
        #    if count==5:
        #        print 'before recluster'
        #        sampler.Recluster(X_i)
        #        print 'after recluster'
        #        NumRecluster+=1
        #        count=0
        #        #sampler.EllipsoidalRescaling(X_i);

        #    Flag1 = True
        #else:
        #    print sampler.ClusteringQuality(X_i)
        #    sampler.EllipsoidalRescaling(X_i);
        #    print sampler.ClusteringQuality(X_i)
         #   Flag1 = False
        nest+=1 
        zold = logzinfo[1]
        sampler.getlogZ(logzinfo)
        if(nest%1000==0):
            print "logZ after ",nest+Np," iterations: %f" % logzinfo[1]
        #print nest,X_i
        #Flag = False
        Flag = THRESH < abs(zold-logzinfo[1])
        count +=1
        if(count>2000):
            for i in xrange(1000):
                sampler.ResetWorstPoint(theta) 
                model.Get_L(theta,templogL,1)
                print theta[0],theta[1],templogL
            break
        #Flag = THRESH < abs(X_i*logLmax)
    #print 'before output'
    #output
    return


if __name__=='__main__':
    main()
    
