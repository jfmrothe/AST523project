cd#!/usr/bin/env python 

import os
import sys
import subprocess
import numpy as np
import time


SampleEffs = [2.0,1.5,1.25,1.1,1.0,0.9,0.75,0.5]
N = 100
Nrep = 100

for SampleEff in SampleEffs:
    for i in range(Nrep):
        time.sleep(1)
        os.system('echo "./multinest '+str(N)+' '+str(SampleEff)+' " >> desired_sampl_eff_auto.txt') 
        os.system('(time ./multinest '+str(N)+' '+str(SampleEff)+'  >> desired_sampl_eff_auto.txt) &>> desired_sampl_eff_auto.txt ')
