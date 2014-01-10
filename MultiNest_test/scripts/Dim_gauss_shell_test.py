#!/usr/bin/env python

import os
import subprocess
import numpy as np

os.chdir('..')

for D in ['2']:

    runtime = open("RunTimeFiles/gauss_shell.txt", 'w')
    runtime.write("problem     : gauss_shell\n")
    runtime.write("n_column    : 0\n")
    runtime.write("Dimension   : "+D+"\n")
    runtime.write("Npoints     : 50\n")
    runtime.write("e           : 1.0\n")
    runtime.write("repartition : 1.2\n")
    for d in range(int(D)):
        runtime.write("param_"+str(d+1)+"     : uniform     -6.0    6.0\n")

    runtime.close()

    subprocess.call(["./multinest", "gauss_shell.txt"])
