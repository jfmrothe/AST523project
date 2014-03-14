import matplotlib
from matplotlib import pyplot as plt
import numpy as np


logfile = open("log",'r')
logdata = logfile.read()
logfile.close()
logdata = logdata.split("\n")
lines=[]
for line in logdata:
    lines.append(line.split(" "))

print lines[1:20]
logZ=[[]]
Niter=[[]]

for line in lines:  
    if len(line) > 1:
        if line[0] == "logZ":
            Niter[-1].append(line[3])
            logZ[-1].append(line[6])
        if line[2] == "iterations":
            Niter[-1].append(line[5])
        if line[0] == "global":
            logZ[-1].append(line[4])
            Niter.append([])
            logZ.append([])

Niter.pop()
logZ.pop()

#logZ = np.array(logZ,float)
#Niter = np.array(Niter,int)
print logZ
print Niter

for i in range(len(logZ)-1):
    logZ[i] =  np.array(logZ[i],float)
for i in range(len(Niter)-1):
    Niter[i] =  np.array(Niter[i],int)

fig = plt.figure()

plt.ylabel("logZ")
plt.xlabel("number of iterations")
for i in range(len(Niter)):
    plt.plot(Niter[i],logZ[i])
for i in range(len(Niter)):
    plt.scatter((int) ( Niter[i][-1]),(float) (logZ[i][-1]),s=12)
#plt.plot(Niter[0],logZ[0])
#plt.plot(Niter[1],logZ[1])
#plt.plot(Niter[2],logZ[2])


plt.show()
