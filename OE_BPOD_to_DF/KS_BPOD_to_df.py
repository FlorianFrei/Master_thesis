# -*- coding: utf-8 -*-
"""
Created on Mon Aug  5 15:15:16 2024

@author: FlorianFreitag
"""
import os

abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)

import numpy as np
import pandas as pd
from helpers import*
basefolder ="E:/Florian/Data/Opto_1/O2_1"
mat_name = next(file for file in os.listdir(basefolder) if file.endswith('.mat'))

#%% load data


cluster_info= pd.read_csv(basefolder + str('/sorted/phy/cluster_info.tsv'),sep = '\t')
clust = np.array (np.load(basefolder + str('/sorted/phy/spike_clusters.npy')))
times = np.array(np.load(basefolder + str('/sorted/phy/spike_times.npy')))[:,0]
ITI = pd.read_csv(basefolder + str('/Meta/ttl_2.csv'))
raw_BPOD = load_mat(basefolder + '/' + mat_name)
surface_estimate = pd.read_csv(basefolder + str('/Meta/surface_estimate.csv'),header=None).iloc[0, 0]

sample_freq = 30000

#%% format ITI


ITI = ITI.iloc[:,0]
#%%
raw_BPOD = is_M9_1(raw_BPOD,mat_name)

#%%
#ITI =ITI[1:].reset_index(drop=True)
proceed = check_trialnumber_matches(ITI, raw_BPOD)

#%% wrangle BPOD
BPOD = BPOD_wrangle(raw_BPOD,ITI,proceed)

#%%
BPOD = adjust_BPOD_with_dead_time(BPOD, ITI)

#%% wrangle EPhys
Ephys_good = Ephys_wrangle(cluster_info,clust,times,sample_freq,surface_estimate)

#%% bin EPhys

Ephys_binned = bin_Ephys(Ephys_good,bin_size=0.01)
#%% add BEHV
EPHYS_trimmed = Add_Behv(Ephys_binned, BPOD, raw_BPOD, ITI)


#%% to csv

EPHYS_trimmed.to_csv(basefolder + str('/O2_1.csv'))


#%%

ITIdiff = [ITI[i + 1] - ITI[i] for i in range(len(ITI) - 1)]
Bpoddiff = [raw_BPOD['SessionData']['TrialStartTimestamp'][i + 1] - raw_BPOD['SessionData']['TrialStartTimestamp'][i] for i in range(len(raw_BPOD['SessionData']['TrialStartTimestamp']) - 1)]