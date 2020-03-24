import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.integrate import odeint

# Calculate the predictions

mt_length = 10
kon=9
koff=0.2
t = np.arange(0,20,0.2)
model=[None,None,None]
lattice_length = 0.008

lattice_sites = mt_length/lattice_length


# Case1 : Continuous
steady_bound=mt_length*kon/koff
model[0] = steady_bound*(1-np.exp(-t*koff))
# Case2 : Digits that do not occupy space in the lattice
model[1] =model[0]
# Case3: Digits that exclude each other

kon=np.array([4,9])
def myODE(bound,t):
    free_sites = lattice_sites-np.sum(bound)
    speed = -bound*koff+kon*free_sites*lattice_length
    return speed

model[2]=odeint(myODE,[0,0],t)

color = ["blue","orange"]
for j,f in enumerate(['runs1','runs2','runs3']):
    scan_folder = os.path.join(f,'scan')
    run_folders = os.listdir(scan_folder)
    run_folders.sort()
    plt.figure()
    for run_folder in  run_folders:
        report_file = os.path.join(scan_folder,run_folder,'results_single.txt')
        arr=np.genfromtxt(report_file,delimiter=" ")
        for i in range(2):
            plt.plot(arr[:,i],c=color[i],alpha=0.4)
        if model[j] is not None:
            plt.plot(model[j],'black')
    plt.savefig('plot%u.svg'%j)



plt.show()