import sys
import os.path
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm
from matplotlib.colors import Normalize
from itertools import cycle
from plot_utils import format_plot

def parse_projections(path):
    with open(path,'r') as file:
        lines = file.readlines()

    line = lines[2]
    line = line.split()
    n_quasi = int(line[0])
    n_ess = int(line[1])
    n_calc = int(line[2])

    projections = np.zeros((n_quasi,n_ess,n_calc))

    i = 0
    j = 0
    for line in lines[3:]:
        line = line.split()
        if len(line) > 1:
            for k in range(len(line)):
                projections[k,j,i] = float(line[k])
            j += 1
        else:
            i += 1
            j = 0

    return projections

lines = ["-","--","-.",":"]
linecycler = cycle(lines)


fig_re,ax_re = plt.subplots()
fig_im,ax_im = plt.subplots()
fig_trajs,ax_trajs = plt.subplots()
ax_im.set_yscale('log')

create_projs_plot = True

for folder in sys.argv[1:]:
    energies = np.loadtxt(f'{folder}/energies.out',dtype = np.complex128)
    energies_real = np.real(energies)
    rates = -2*np.imag(energies)

    #cmap = matplotlib.cm.get_cmap('viridis')
    #max_block = np.loadtxt(f'{folder}/max_block.out')
    #colors = cmap(max_block)

    sim_type = None
    if os.path.isfile(f'{folder}/omega.out'):
        sim_type = 'omega'
        x_variable = np.loadtxt(f'{folder}/omega.out')

    elif os.path.isfile(f'{folder}/intensity.out'):
        sim_type = 'intensity'
        x_variable = np.loadtxt(f'{folder}/intensity.out')

    elif os.path.isfile(f'{folder}/static_field.out'):
        sim_type = 'static'
        x_variable = np.loadtxt(f'{folder}/static_field.out')

    else:
        raise RuntimeError('Could not find appropriate data file for x values')

    projs = False
    if os.path.isfile(f'{folder}/projections.out'):
        projs = True
        projections = parse_projections(f'{folder}/projections.out')
        n_ess = projections.shape[1]
        if create_projs_plot:
            fig_projs,ax_projs = plt.subplots(1,n_ess,sharey='row',figsize=(4.0*n_ess,4.8))
            for ax in ax_projs:
                ax.set_yscale('log')
            create_projs_plot = False


    linestyle = next(linecycler)
    for i in range(energies.shape[1]):
    #for j in range(x_variable.shape[0]):
        #ax_re.plot(x_variable[j],energies_real[j,i],'s',color= colors[j,i])
        #ax_im.plot(x_variable[j],rates[j,i],'s',color = colors[j,i])
        #ax_re.plot(x_variable,energies_real[:,i],'+')
        ax_re.plot(x_variable,energies_real[:,i],linestyle=linestyle)
        ax_im.plot(x_variable,rates[:,i],linestyle=linestyle)
        ax_trajs.plot(energies_real[:,i],rates[:,i],linestyle=linestyle)
        if projs:
            for j in range(n_ess):
                ax_projs[j].plot(x_variable,projections[i,j,:],linestyle=linestyle)
                ax_projs[j].set_title(str(j+1),fontsize = 18)

if sim_type == 'omega':
    print('')
    format_plot(fig_re,ax_re,'omega [a.u.]','Re $E$ [a.u.]')
    format_plot(fig_im,ax_im,'omega [a.u.]','Rate [a.u.]')
    if projs:
        format_plot(fig_projs,ax_projs[0],'$\omega$ [a.u.]','Projections')
        for j in range(1,n_ess):
            format_plot(fig_projs,ax_projs[j],'$\omega$ [a.u.]','')

elif sim_type == 'intensity':
    ax_re.set_xscale('log')
    ax_im.set_xscale('log')
    format_plot(fig_re,ax_re,'Intensity [W/cm$^2$]','Re $E$ [a.u.]')
    format_plot(fig_im,ax_im,'Intensity [W/cm$^2$]','Rate [a.u.]')
    if projs:
        format_plot(fig_projs,ax_projs[0],'Intensity [W/cm$^2$]','Projections')
        for j in range(1,n_ess):
            format_plot(fig_projs,ax_projs[j],'Intensity [W/cm$^2$]','')

elif sim_type == 'static':
    ax_re.set_xscale('log')
    ax_im.set_xscale('log')

#format_plot(fig_trajs,ax_trajs,'Re $E$ [a.u.]','Rate [a.u.]')

plt.show()
