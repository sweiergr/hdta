"""

	Plot the heatmaps for our simulations.

	Input files from estimation:
	- fig1Data.mat
	- CCPCounterData.mat
	- CCPCounterNormTermData.mat
	- L3_data.mat
	- L3naive_data.mat
	- P4_data.mat   
	
	Output files:
	- UdistanceL3.pdf
    - UdistanceNaiveL3.pdf
    - UdistanceP4.pdf
    - CCPAdoption_counter.pdf
    - CCPAdoption_counter_NormTerm.pdf
    - CCPAdoption.pdf

"""

import numpy as np
from bld.project_paths import project_paths_join
from matplotlib import pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.dates as mdates
import pandas as pd
import scipy.io as sio

## Create first heatmap: UdistanceL3.

# Load data from MATLAB export.
l3_fig_dict = sio.loadmat(project_paths_join('OUT_ANALYSIS', 'L3_data.mat'))
# Just for convenience I extract the different objects from the dictionary.
l3_data = l3_fig_dict['L3_data']
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
# gaussian_filter(Z, 5.), 4, colors='k',interpolation='none'
fig.colorbar(cax)
ax.set_ylabel(r'$\delta$')
ax.set_xlabel(r'$\beta$')
# Save the figure.
plt.savefig(project_paths_join('OUT_FIGURES', 'UdistanceL3.pdf'), bbox_inches ='tight')


## Create second heatmap: UdistanceL3naive.
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
# gaussian_filter(Z, 5.), 4, colors='k',interpolation='none'
fig.colorbar(cax)
ax.set_ylabel(r'$\delta$')
ax.set_xlabel(r'$\beta$')
# Save the figure.
plt.savefig(project_paths_join('OUT_FIGURES', 'UdistanceNaiveL3.pdf'), bbox_inches ='tight')

## Create third heatmap: UdistanceP4.
# Load data from MATLAB export.
l3_fig_dict = sio.loadmat(project_paths_join('OUT_ANALYSIS', 'P4_data.mat'))
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
# gaussian_filter(Z, 5.), 4, colors='k',interpolation='none'
fig.colorbar(cax)
ax.set_ylabel(r'$\delta$')
ax.set_xlabel(r'$\beta$')
# Save the figure.
plt.savefig(project_paths_join('OUT_FIGURES', 'UdistanceP4.pdf'), bbox_inches ='tight')
####### END PLOTTING CODE FOR HEATMAPS #########################################

################################################################################
# LINE PLOTS.
# Load data from MATLAB export.
fig1_dict = sio.loadmat(project_paths_join('OUT_ANALYSIS', 'fig1Data.mat'))
# Just for convenience I extract the different objects from the dictionary.
fig1_data = fig1_dict['CCP_data']
# Slice CPP data.
CCP_total = fig1_data[:,0:2]
CCP_total_sophi = fig1_data[:,2:4]
CCP_total_naive = fig1_data[:,4:]
# This should eventually be automated, but to get an example graph quickly, I hardcode it for now.
supportX = np.array([2,3,7,9])
nPeriods = 4
periods = np.linspace(1,nPeriods,nPeriods)
nperiod_axis = [];
for s in periods:
	nperiod_axis.append('T-'+str(nPeriods - int(s)))
nperiod_axis[-1] = 'T'
nSuppX = len(supportX)
nXT = len(fig1_data)
Ttemp = np.reshape(np.linspace(1,nXT,nXT)-1,[nSuppX,nPeriods]).astype(int)
# Recall that Python starts at 0, so if we want action 1, we need to use column 0
action = 1 - 1
legend_fig1 = ['time-consistent','sophisticated','naive']
lwidth=0.5

## Create Figure 1.
# Loop over subplots/state values
fig, axs = plt.subplots(2, 2)
for i, val in enumerate(supportX):
	# print (i, ",",val)
	row_idx = np.floor(i / 2).astype(int)
	col_idx = (i%2)
	axs[row_idx, col_idx].plot(periods,CCP_total[Ttemp[:,i],action],'g',linewidth=lwidth,label='time-consistent')
	axs[row_idx, col_idx].set_title('X='+str(val))
	axs[row_idx, col_idx].plot(periods,CCP_total_sophi[Ttemp[:,i],action],'b',linewidth=lwidth,label='sophisticated')
	axs[row_idx, col_idx].set_title('X='+str(val))
	axs[row_idx, col_idx].plot(periods,CCP_total_naive[Ttemp[:,i],action],'r',linewidth=lwidth,label='naive')
	axs[row_idx, col_idx].set_title('X='+str(val))
	axs[row_idx, col_idx].set_ylim([0,1])
	axs[row_idx, col_idx].set_xticks(periods, minor=False)
	axs[row_idx, col_idx].set_xticklabels(nperiod_axis, fontdict=None, minor=False)
	handles, labels = axs[row_idx,col_idx].get_legend_handles_labels()
# box = ax.get_position()
fig.legend(handles, labels, loc='lower center',ncol=3,bbox_to_anchor=(0.5, -0.075) )
fig.tight_layout()
# Save the figure.
plt.savefig(project_paths_join('OUT_FIGURES', 'CCPAdoption.pdf'), bbox_inches ='tight')



## Create Figure CCPs normalization on non-adoption.
# Load data from MATLAB export.
fig1_dict = sio.loadmat(project_paths_join('OUT_ANALYSIS', 'CCPCounterData.mat'))
# Just for convenience I extract the different objects from the dictionary.
fig1_data = fig1_dict['CCP_counter_data']
# Slice CPP data.
CCP_sophi_longer = fig1_data[:,0:2]
CCP_counter = fig1_data[:,2:4]
CCP_counter_norm = fig1_data[:,4:]
# This should eventually be automated, but to get an example graph quickly, I hardcode it for now.
supportX = np.array([2,3,7,9])
nPeriods = 10
periods = np.linspace(1,nPeriods,nPeriods)
nperiod_axis = [];
for s in periods:
	nperiod_axis.append('T-'+str(nPeriods - int(s)))
nperiod_axis[-1] = 'T'
nSuppX = len(supportX)
nXT = len(fig1_data)
Ttemp = np.reshape(np.linspace(1,nXT,nXT)-1,[nPeriods,nSuppX]).astype(int)
# Recall that Python starts at 0, so if we want action 1, we need to use column 0
action = 1 - 1
lwidth=0.5

# Loop over subplots/state values
fig, axs = plt.subplots(2, 2)
for i, val in enumerate(supportX):
	# print (i, ",",val)
	row_idx = np.floor(i / 2).astype(int)
	col_idx = (i%2)
	axs[row_idx, col_idx].plot(periods,CCP_sophi_longer[Ttemp[:,i],action],'g',linewidth=lwidth,label='true model')
	axs[row_idx, col_idx].set_title('X='+str(val))
	axs[row_idx, col_idx].plot(periods,CCP_counter[Ttemp[:,i],action],'b',linewidth=lwidth,label='cf w/o norm.')
	axs[row_idx, col_idx].set_title('X='+str(val))
	axs[row_idx, col_idx].plot(periods,CCP_counter_norm[Ttemp[:,i],action],'r',linewidth=lwidth,label='cf w/ norm.')
	axs[row_idx, col_idx].set_title('X='+str(val))
	axs[row_idx, col_idx].set_ylim([0,1])
	axs[row_idx, col_idx].set_xticks(periods, minor=False)
	axs[row_idx, col_idx].set_xticklabels(nperiod_axis, fontdict=None, minor=False)
	axs[row_idx, col_idx].xaxis.set_major_locator(ticker.MultipleLocator(2))
	handles, labels = axs[row_idx,col_idx].get_legend_handles_labels()
# box = ax.get_position()
fig.legend(handles, labels, loc='lower center',ncol=3,bbox_to_anchor=(0.5, -0.075) )
fig.tight_layout()
# Save the figure.
plt.savefig(project_paths_join('OUT_FIGURES', 'CCPAdoption_counter.pdf'), bbox_inches ='tight')


## Create Figure CCPs normalization on adoption.
# Load data from MATLAB export.
fig1_dict = sio.loadmat(project_paths_join('OUT_ANALYSIS', 'CCPCounterNormTermData.mat'))
# Just for convenience I extract the different objects from the dictionary.
fig1_data = fig1_dict['CCP_counter_normterm_data']
# Slice CPP data.
CCP_sophi_longer = fig1_data[:,0:2]
CCP_counter = fig1_data[:,2:4]
CCP_counter_normTerm = fig1_data[:,4:]
# This should eventually be automated, but to get an example graph quickly, I hardcode it for now.
supportX = np.array([2,3,7,9])
nPeriods = 10
periods = np.linspace(1,nPeriods,nPeriods)
nperiod_axis = [];
for s in periods:
	nperiod_axis.append('T-'+str(nPeriods - int(s)))
nperiod_axis[-1] = 'T'
nSuppX = len(supportX)
nXT = len(fig1_data)
Ttemp = np.reshape(np.linspace(1,nXT,nXT)-1,[nPeriods,nSuppX]).astype(int)
# Recall that Python starts at 0, so if we want action 1, we need to use column 0
action = 1 - 1
lwidth=0.5

# Loop over subplots/state values
fig, axs = plt.subplots(2, 2)
for i, val in enumerate(supportX):
	# print (i, ",",val)
	row_idx = np.floor(i / 2).astype(int)
	col_idx = (i%2)
	axs[row_idx, col_idx].plot(periods,CCP_sophi_longer[Ttemp[:,i],action],'g',linewidth=lwidth,label='true model')
	axs[row_idx, col_idx].set_title('X='+str(val))
	axs[row_idx, col_idx].plot(periods,CCP_counter[Ttemp[:,i],action],'b',linewidth=lwidth,label='cf w/o norm.')
	axs[row_idx, col_idx].set_title('X='+str(val))
	axs[row_idx, col_idx].plot(periods,CCP_counter_normTerm[Ttemp[:,i],action],'r',linewidth=lwidth,label='cf w/ norm.')
	axs[row_idx, col_idx].set_title('X='+str(val))
	axs[row_idx, col_idx].set_ylim([0,1])
	axs[row_idx, col_idx].set_xticks(periods, minor=False)
	axs[row_idx, col_idx].set_xticklabels(nperiod_axis, fontdict=None, minor=False)
	axs[row_idx, col_idx].xaxis.set_major_locator(ticker.MultipleLocator(2))
	handles, labels = axs[row_idx,col_idx].get_legend_handles_labels()
# box = ax.get_position()
fig.legend(handles, labels, loc='lower center',ncol=3,bbox_to_anchor=(0.5, -0.075) )
fig.tight_layout()
# Save the figure.
plt.savefig(project_paths_join('OUT_FIGURES', 'CCPAdoption_counter_NormTerm.pdf'), bbox_inches ='tight')


## Create Figure CCPs normalization on adoption with confidence band.
# Load data from MATLAB export.
fig1_dict = sio.loadmat(project_paths_join('OUT_ANALYSIS', 'CCPCounterNormTermBandData.mat'))
# Just for convenience I extract the different objects from the dictionary.
fig1_data = fig1_dict['CCP_counter_normTerm_banddata']
# Slice CPP data.
CCP_sophi_longer = fig1_data[:,0:1]
CCP_counter = fig1_data[:,1:4]
CCP_counter_normTerm = fig1_data[:,4:]
# This should eventually be automated, but to get an example graph quickly, I hardcode it for now.
supportX = np.array([2,3,7,9])
nPeriods = 10
periods = np.linspace(1,nPeriods,nPeriods)
nperiod_axis = [];
for s in periods:
	nperiod_axis.append('T-'+str(nPeriods - int(s)))
nperiod_axis[-1] = 'T'
nSuppX = len(supportX)
nXT = len(fig1_data)
Ttemp = np.reshape(np.linspace(1,nXT,nXT)-1,[nPeriods,nSuppX]).astype(int)
# Recall that Python starts at 0, so if we want action 1, we need to use column 0
# action = 1 - 1
lwidth=0.5

# Loop over subplots/state values
fig, axs = plt.subplots(2, 2)
for i, val in enumerate(supportX):
	# print (i, ",",val)
	row_idx = np.floor(i / 2).astype(int)
	col_idx = (i%2)
	axs[row_idx, col_idx].plot(periods,CCP_sophi_longer[Ttemp[:,i]],'g',linewidth=lwidth,label='true model')
	axs[row_idx, col_idx].set_title('X='+str(val))
	axs[row_idx, col_idx].plot(periods,CCP_counter[Ttemp[:,i],0],'b',linewidth=lwidth,label='cf w/o norm.')
	axs[row_idx, col_idx].fill_between(periods, (CCP_counter[Ttemp[:,i],1]), (CCP_counter[Ttemp[:,i],2]), color='b', alpha=.1)
	axs[row_idx, col_idx].set_title('X='+str(val))
	axs[row_idx, col_idx].plot(periods,CCP_counter_normTerm[Ttemp[:,i],0],'r',linewidth=lwidth,label='cf w/ norm.')
	axs[row_idx, col_idx].fill_between(periods, (CCP_counter_normTerm[Ttemp[:,i],1]), (CCP_counter_normTerm[Ttemp[:,i],2]), color='r', alpha=.1)
	axs[row_idx, col_idx].set_title('X='+str(val))
	axs[row_idx, col_idx].set_ylim([0,1])
	axs[row_idx, col_idx].set_xticks(periods, minor=False)
	axs[row_idx, col_idx].set_xticklabels(nperiod_axis, fontdict=None, minor=False)
	axs[row_idx, col_idx].xaxis.set_major_locator(ticker.MultipleLocator(2))
	handles, labels = axs[row_idx,col_idx].get_legend_handles_labels()
# box = ax.get_position()
fig.legend(handles, labels, loc='lower center',ncol=3,bbox_to_anchor=(0.5, -0.075) )
fig.tight_layout()
# Save the figure.
plt.savefig(project_paths_join('OUT_FIGURES', 'CCPAdoption_counter_NormTermBand.pdf'), bbox_inches ='tight')



