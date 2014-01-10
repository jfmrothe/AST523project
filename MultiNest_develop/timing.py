#!/usr/bin/env python 

import os
import sys
import subprocess
import numpy as np
import time


Ns = [3,5,10,50,100,500,1000]

for N in Ns:
    time.sleep(1)
    os.system('echo "time ./multinest '+str(N)+' " >> timing.txt')
    os.system('(time ./multinest '+str(N)+' >> timing.txt) &>> timing.txt')
    os.system('echo "\n" >> timing.txt' )
