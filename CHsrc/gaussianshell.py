#!/usr/bin/env python
import numpy as np
import scipy as sp
import cfg_parse as cfgp
from dataio import readcolumn
class gaussianshellmodel():
    def __init__(self,cfgfile):
        #self.infile_ = cfgp.File_parse(cfgfile,'infile')
        self.outfile_ = cfgp.File_parse(cfgfile,'outfile')
        self.inpath_ = cfgp.File_parse(cfgfile,'inpath')
        #self.data_ = []; readcolumn(self.data_,1,self.infile_); self.data_ = np.array(self.data_)
        self.D=3
        self.Np_=1000
        self.var0_=[0. for i in range(self.D)]
        self.varerr_=[6.  for i in range(self.D)]
        self.repartition = 1.2
        self.thresh = 1.0*10**-2.5
        self.fixparams_=[2.0,0.1,-3.5,0,3.5,0.0]; #r,ww,c1[0],c1[1],c2[0],c2[1], all other center coors are 0
        self.neginf = -1.e7
	self.NL_=0
        self.c1 = np.zeros(self.D)
        self.c1[:2] = self.fixparams_[2:4]
        self.c2 = np.zeros(self.D)
        self.c2[:2] = self.fixparams_[4:6]
        return
    
    def Getinitial(self):
        minvals = np.array(self.var0_)-np.array(self.varerr_)
        maxvals = np.array(self.var0_)+np.array(self.varerr_)
        guessvals = np.array(self.var0_)
        eff = 0.3
        return [minvals,maxvals,eff,self.Np_]

    def Get_dist(self,r1,r2):
        #print r1,r2
        if len(r1)!=len(r2):
            print "Gaussianshell Get_dist dimension mismatch."
            return
        return np.sum([(r1[i]-r2[i])**2 for i in range(len(r1))])**0.5
    def Get_L(self,model_params, logL, nl):
	self.NL_+=1
        norm = 1.0/np.sqrt(2*np.pi*self.fixparams_[1]**2.)
        
        for i in xrange(nl):
            #print model_params[i*self.D],model_params[i*self.D+1]
            coor = model_params[i*self.D:(i+1)*self.D]
            index1 = (self.Get_dist(coor,self.c1)-self.fixparams_[0])**2.
            index2 = (self.Get_dist(coor,self.c2)-self.fixparams_[0])**2.
            L = norm*(np.exp(-index1/2./self.fixparams_[1]**2.)+np.exp(-index2/2./self.fixparams_[1]**2.))
            #L = np.log(y/np.pi/((self.data_-x)**2.+y**2.))
            if(L ==0):
                logL[i] = self.neginf
            else:
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
            logLtemp=np.array([0])
            self.Get_L(posterior[i*self.D:(i+1)*self.D],logLtemp,1)
            fout.write('%f %f %9.6e\n' % (prob[i],-1.0*i/self.Np_,Wts[i]))
            # format: theta_j,i logL_i logX_i DeltaX*L_i

        fout.close()
