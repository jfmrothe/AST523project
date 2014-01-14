#!/usr/bin/env python
from multinest import MNest
os.chdir("oblateTransit")
#os.system("make clean")
os.system("make all")
os.chdir("..")
sys.path.append("oblateTransit")
import transitmodel as model #module that call occultquad, and also hold the data
#from occultquad import occultquad
from Point import Point
import Ellipsoid
from Samplers import Samplers as MNest

def main():
    minvals,maxvals,guessvals = model.Getinitial()
    sampler = MNest(priortypes,minvals,maxvals,guessvals,Np)
    #what I intended to use, for sampler does not need to know the names
    #sampler = MNest(priortypes,minvals,maxvals,guessvals) 
    print "running MultiNest algorithm... this may take a few minutes\n"
    nest = 0
    Flag = True
    NumRecluster = 0
    while Flag: 
        X_i = exp(-nest/N)
        sampler.getAlltheta(Alltheta, nx, ny)
        model.Get_L(Alltheta,nx,ny,logL,nl)
        #discard and resample, get logLmax for convergence check as byproduct
        sampler.DisgardWorstPoint(logL, nl,nest)
        #CH: need to take the data_obj out, and think about how to do that
        sampler.ResetWorstPoint(data_obj)
        #ellipsoidal partitioning 
        if(nest==0 or sampler.ClusteringQuality(X_i) > RepartitionFactor): 
            # recluster?
            sampler.Recluster(X_i)
            NumRecluster+=1
       nest+=1 
       Flag = THRESH < abs(X_i*logLmax)
    
    #output
    return

if __name__=='__main__'
    main()
