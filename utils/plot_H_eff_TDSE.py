import sys
import os.path
import numpy as np
import scipy.integrate as si
import matplotlib.pyplot as plt
import matplotlib.cm
from matplotlib.colors import Normalize
from itertools import cycle

#sys.path.append('./format-plot/')
#from plot_utils import format_plot

def parse_H_eff(path):
    with open(path,'r') as file:
        lines = file.readlines()

    line = lines[2]
    line = line.split()
    n_quasi = int(line[0])
    n_ess = int(line[1])
    n_calc = int(line[2])

    projections = np.zeros((n_calc,n_quasi,n_ess),dtype=np.complex128)

    i = 0
    j = 0
    for line in lines[4:]:
        line = line.split()
        if len(line) > 1:
            for k in range(len(line)):
                projections[i,j,k] = complex(line[k])
            j += 1
        else:
            i += 1
            j = 0

    return projections

def RHS(t,y,H_eff):
    return -1j*H_eff@y

def solve_TDSE(H_eff,T,max_dt):
    y_0 = np.zeros(H_eff.shape[0],dtype=np.complex128)
    y_0[0] = 1
    t_lim = (0,T)
    sol  = si.solve_ivp(RHS,t_lim,y_0,method='DOP853',args=(H_eff,),max_step=max_dt)
    return sol

def plot_TDSE(sol,title):
    lines = ["-","--","-.",":"]
    linecycler = cycle(lines)
    fontsize = 18
    fig,ax = plt.subplots()

    for i in range(sol.y.shape[0]):
        ax.plot(sol.t,np.abs(sol.y[i,:])**2,linestyle = next(linecycler),label = str(i))
    ax.plot(sol.t,np.sum(np.abs(sol.y)**2,axis=0),linestyle = next(linecycler),label = 'Total')

    ax.set_title(title,fontsize=fontsize)
    ax.set_xlabel('Time [a.u.]', fontsize = fontsize)
    ax.set_ylabel('Population',fontsize = fontsize)
    ax.tick_params(labelsize = fontsize)
    ax.legend(framealpha = 0.0,fontsize = 16)
    fig.tight_layout()


T = float(sys.argv[1])
max_dt = float(sys.argv[2])
folder = sys.argv[3]

os.path.isfile(f'{folder}/H_eff.out')
H_eff = parse_H_eff(f'{folder}/H_eff.out')

for i in range(H_eff.shape[0]):
    print(i)
    sol = solve_TDSE(H_eff[i,:,:],T,max_dt)
    plot_TDSE(sol,i)


plt.show()
