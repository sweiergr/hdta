"""
	Plot the heatmaps from our simulations.

"""
import numpy as np
from bld.project_paths import project_paths_join
#print("BLD IMPORT SUCCESSFUL!")
from matplotlib import pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.dates as mdates
import pandas as pd
import scipy.io as sio

## Create heatmap for Figure E.4.
# Load data from MATLAB export.
l3_fig_dict = sio.loadmat(project_paths_join('OUT_ANALYSIS', 'L3naive_data.mat'))
l3_data = l3_fig_dict['L3naive_data']
bb_data = l3_data[:,:,0]
dd_data = l3_data[:,:,1]
objfun_data = l3_data[:,:,2]
contour_levels_small = np.linspace(0,0.175,15)
contour_levels_large = np.linspace(0.2,1,30)
contour_levels = np.concatenate( (contour_levels_small, contour_levels_large),axis=0)
plt.clf()
fig, ax = plt.subplots(figsize=(8,6))
cax = ax.pcolormesh(bb_data,dd_data,objfun_data,shading='auto',linewidth=0,rasterized=True)
cax.set_edgecolor('face')
plt.contour(bb_data, dd_data, objfun_data,colors='k',linewidths=0.4,levels = contour_levels)
fig.colorbar(cax)
ax.set_ylabel(r'$\delta$')
ax.set_xlabel(r'$\beta$')
plt.savefig(project_paths_join('OUT_FIGURES', 'FigureE4.pdf'), bbox_inches ='tight')

## Create heatmap for Figure E.3
# Load data from MATLAB export.
l3_fig_dict = sio.loadmat(project_paths_join('OUT_ANALYSIS', 'P4_data.mat'))
# Just for convenience I extract the different objects from the dictionary.
l3_data = l3_fig_dict['P4_data']
bb_data = l3_data[:,:,0]
dd_data = l3_data[:,:,1]
objfun_data = l3_data[:,:,2]
contour_levels_small = np.linspace(0,0.175,15)
contour_levels_large = np.linspace(0.2,1,30)
contour_levels = np.concatenate( (contour_levels_small, contour_levels_large),axis=0)

plt.clf()
fig, ax = plt.subplots(figsize=(8,6))
cax = ax.pcolormesh(bb_data,dd_data,objfun_data,shading='auto',linewidth=0,rasterized=True)
cax.set_edgecolor('face')
plt.contour(bb_data, dd_data, objfun_data,colors='k',linewidths=0.4,levels = contour_levels)
fig.colorbar(cax)
ax.set_ylabel(r'$\delta$')
ax.set_xlabel(r'$\beta$')
plt.savefig(project_paths_join('OUT_FIGURES', 'FigureE3.pdf'), bbox_inches ='tight')
####### END PLOTTING CODE FOR HEATMAPS #########################################

################################################################################
# BEGIN LINE PLOTS.
# Load data from MATLAB export.
fig1_dict = sio.loadmat(project_paths_join('OUT_ANALYSIS', 'fig1Data.mat'))
fig1_data = fig1_dict['CCP_data']
# Slice CPP data.
CCP_total = fig1_data[:,0:2]
CCP_total_sophi = fig1_data[:,2:4]
CCP_total_naive = fig1_data[:,4:]
supportX = ["$x_1 = x_{1L}, x_2 = x_{2L}$","$x_1 = x_{1L}, x_2 = x_{2H}$","$x_1 = x_{1H}, x_2 = x_{2L}$","$x_1 = x_{1H}, x_2 = x_{2H}$"]
nPeriods = 4
periods = np.linspace(1,nPeriods,nPeriods)
nperiod_axis = [];
for s in periods:
	nperiod_axis.append('T-'+str(nPeriods - int(s)))
nperiod_axis[-1] = 'T'
nSuppX = len(supportX)
nXT = len(fig1_data)
Ttemp = np.reshape(np.linspace(1,nXT,nXT)-1,[nSuppX,nPeriods]).astype(int)
action = 1 - 1
legend_fig1 = ['time-consistent','sophisticated','naive']
lwidth=0.5

## Create Figure D.1
# Loop over subplots/state values
fig, axs = plt.subplots(2, 2)
for i, val in enumerate(supportX):
	row_idx = np.floor(i / 2).astype(int)
	col_idx = (i%2)
	axs[row_idx, col_idx].plot(periods,CCP_total[Ttemp[:,i],action],'g',linewidth=lwidth,label='time-consistent')
	axs[row_idx, col_idx].set_title(str(val))
	axs[row_idx, col_idx].plot(periods,CCP_total_sophi[Ttemp[:,i],action],'b',linewidth=lwidth,label='sophisticated')
	axs[row_idx, col_idx].set_title(str(val))
	axs[row_idx, col_idx].plot(periods,CCP_total_naive[Ttemp[:,i],action],'r',linewidth=lwidth,label='naive')
	axs[row_idx, col_idx].set_title(str(val))
	axs[row_idx, col_idx].set_ylim([0,1])
	axs[row_idx, col_idx].set_xticks(periods, minor=False)
	axs[row_idx, col_idx].set_xticklabels(nperiod_axis, fontdict=None, minor=False)
	handles, labels = axs[row_idx,col_idx].get_legend_handles_labels()
fig.legend(handles, labels, loc='lower center',ncol=3,bbox_to_anchor=(0.5, -0.075) )
fig.tight_layout()
plt.savefig(project_paths_join('OUT_FIGURES', 'FigureD1.pdf'), bbox_inches ='tight')

## Create Figure D.2
# Load data from MATLAB export.
fig1_dict = sio.loadmat(project_paths_join('OUT_ANALYSIS', 'CCPCounterdiffData.mat'))
fig1_data = fig1_dict['CCP_counter_diff_data']
# Slice CPP data.
CCP_true_counter_diff = fig1_data[:,0:2]
CCP_counter_diff = fig1_data[:,2:4]
CCP_counter_diff_normTerm = fig1_data[:,4:]
supportX = ["$x_1 = x_{1L}, x_2 = x_{2L}$","$x_1 = x_{1L}, x_2 = x_{2H}$","$x_1 = x_{1H}, x_2 = x_{2L}$","$x_1 = x_{1H}, x_2 = x_{2H}$"]
nPeriods = 10
periods = np.linspace(1,nPeriods,nPeriods)
nperiod_axis = [];
for s in periods:
	nperiod_axis.append('T-'+str(nPeriods - int(s)))
nperiod_axis[-1] = 'T'
nSuppX = len(supportX)
nXT = len(fig1_data)
Ttemp = np.reshape(np.linspace(1,nXT,nXT)-1,[nPeriods,nSuppX]).astype(int)
action = 1 - 1
lwidth=0.5

# Loop over subplots/state values
fig, axs = plt.subplots(2, 2)
for i, val in enumerate(supportX):
	row_idx = np.floor(i / 2).astype(int)
	col_idx = (i%2)
	axs[row_idx, col_idx].plot(periods,CCP_true_counter_diff[Ttemp[:,i],action],'g',linewidth=lwidth,label='cf w/ true model')
	axs[row_idx, col_idx].set_title(str(val))
	axs[row_idx, col_idx].plot(periods,CCP_counter_diff[Ttemp[:,i],action],'b',linewidth=lwidth,label='cf w/o norm')
	axs[row_idx, col_idx].set_title(str(val))
	axs[row_idx, col_idx].plot(periods,CCP_counter_diff_normTerm[Ttemp[:,i],action],'r',linewidth=lwidth,label='cf w/ norm')
	axs[row_idx, col_idx].set_title(str(val))
	axs[row_idx, col_idx].set_ylim([0,1])
	axs[row_idx, col_idx].set_xticks(periods, minor=False)
	axs[row_idx, col_idx].set_xticklabels(nperiod_axis, fontdict=None, minor=False)
	axs[row_idx, col_idx].xaxis.set_major_locator(ticker.MultipleLocator(2))
	handles, labels = axs[row_idx,col_idx].get_legend_handles_labels()
fig.legend(handles, labels, loc='lower center',ncol=3,bbox_to_anchor=(0.5, -0.075) )
fig.tight_layout()
plt.savefig(project_paths_join('OUT_FIGURES', 'FigureD2.pdf'), bbox_inches ='tight')


