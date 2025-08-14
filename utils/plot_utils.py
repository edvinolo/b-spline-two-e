#!/usr/bin/python3
#Collection of functions that are useful for plotting in python
#Mainly related to stuff that useful for formatting the plot, and does not need to be repated in every program
import matplotlib.pyplot as plt


def format_plot(fig,ax,xlabel,ylabel,**kwargs):
    #Pass fig and ax objects to format the plot

    #Set optional arguments
    title = kwargs.get('title',False)
    grid = kwargs.get('grid',True)
    legend = kwargs.get('legend',False)
    fontsize = kwargs.get('fontsize',18)
    labelsize = kwargs.get('labelsize',18)

    #Format the plot
    if len(xlabel) > 0:
        ax.set_xlabel(xlabel,fontsize = fontsize)
    if len(ylabel) > 0:
        ax.set_ylabel(ylabel,fontsize = fontsize)
    if title:
        ax.set_title(title,fontsize = fontsize)
    ax.tick_params(axis = 'both', labelsize = labelsize)
    if grid:
        ax.grid()
    if legend:
        ax.legend(fontsize=fontsize*0.8)
    fig.tight_layout()
    return 