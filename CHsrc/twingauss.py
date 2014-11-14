#!/usr/bin/env python
import numpy as np
import scipy as sp
import cfg_parse as cfgp
from dataio import readcolumn
class twingaussmodel():
    def __init__(self,cfgfile):
        #self.infile_ = cfgp.File_parse(cfgfile,'infile')
        self.outfile_ = cfgp.File_parse(cfgfile,'outfile')
        self.inpath_ = cfgp.File_parse(cfgfile,'inpath')
        #self.data_ = []; readcolumn(self.data_,1,self.infile_); self.data_ = np.array(self.data_)
        self.D=2
        self.Np_=1000
        self.var0_=[1.0/2. for i in range(self.D)]
        self.varerr_=self.var0_
        self.repartition = 1.2
	self.NL_=0
        self.thresh = 1.0*10**-8.0
        return
    
    def Getinitial(self):
        minvals = np.array(self.var0_)-np.array(self.varerr_)
        maxvals = np.array(self.var0_)+np.array(self.varerr_)
        guessvals = np.array(self.var0_)
        eff = 0.3
        return [minvals,maxvals,eff,self.Np_]


    def Get_L(self,model_params, logL, nl):
	self.NL_+=1
        for i in xrange(nl):
            #print model_params[i*self.D],model_params[i*self.D+1]
            coor = model_params[i*self.D:(i+1)*self.D]
            dist1 = (coor[0]-0.25)**2
            dist2 = (coor[0]-0.75)**2
            for c in coor[1:]:
                dist1 += (c-0.5)**2
                dist2 += (c-0.5)**2
            L = np.exp(-100.*np.sum(dist1))+np.exp(-100.*np.sum(dist2))
            #x = model_params[i*self.D]
            #y = model_params[i*self.D+1]
            #L=-100*((x-0.5)**2+(y-0.5)**2)
            ##L = np.log(y/np.pi/((self.data_-x)**2.+y**2.))
            logL[i] = np.log(L) 
            #print x,y,logL[i]
        return
    def Output(self,posterior, prob, Wts):
        fout = open(self.outfile_,mode='w')
        print "output to %s" % self.outfile_
        #logLtemp = 0. 
        for i in xrange(posterior.shape[0]/self.D):
            for j in xrange(int(self.D)):
                fout.write('%f ' % posterior[i*self.D+j])
            fout.write('%f %f %9.6e\n' % (prob[i],-1.0*i/self.Np_,Wts[i]))
            # format: theta_j,i logL_i logX_i DeltaX*L_i
            #logLtemp=np.array([0])
            #self.Get_L(posterior[i*self.D:(i+1)*self.D],logLtemp,1)
            #fout.write('%f %f\n' % (prob[i],logLtemp))
        fout.close()
