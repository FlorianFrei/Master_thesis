# -*- coding: utf-8 -*-
"""
Created on Wed Dec 11 16:33:52 2024

@author: deepl
"""

import numpy as np
import pandas as pd
from helpers2 import*
basefolder ="E:/Florian/Data/Opto_1/O2_T"
mat_name1 = 'MLA46_OptoAir_20241121_155023'
mat_name2 = 'MLA46_OptoAir_20241121_155900'
#%% load data


cluster_info= pd.read_csv(basefolder + str('/sorted/phy/cluster_info.tsv'),sep = '\t')
clust = np.array (np.load(basefolder + str('/sorted/phy/spike_clusters.npy')))
times = np.array(np.load(basefolder + str('/sorted/phy/spike_times.npy')))[:,0]
ITI = pd.read_csv(basefolder + str('/Meta/ttl_2.csv'))
raw_BPOD1 = load_mat(basefolder + '/' + mat_name1)
raw_BPOD2 = load_mat(basefolder + '/' + mat_name2)
surface_estimate = pd.read_csv(basefolder + str('/Meta/surface_estimate.csv'),header=None).iloc[0, 0]

sample_freq = 30000

#%% format ITI

ITI = ITI.iloc[:,0]
time_differences = [ITI[i + 1] - ITI[i] for i in range(len(ITI) - 1)]

# Find the index of the maximum time difference
max_diff_index = time_differences.index(max(time_differences))

# Split the timestamps into two lists
ITI1 = ITI[:max_diff_index + 1]  # Up to the max difference
ITI2 = ITI[max_diff_index + 1:].reset_index(drop=True) # After the max difference


#%%
ITI1 = ITI1[1:].reset_index(drop=True) 

proceed1 = check_trialnumber_matches(ITI1, raw_BPOD1)
proceed2 = check_trialnumber_matches(ITI2, raw_BPOD2)
#%%
BPOD1 = BPOD_wrangle(raw_BPOD1,ITI1, proceed1)
BPOD2 = BPOD_wrangle(raw_BPOD2,ITI2, proceed2)

#%%
BPOD1 = adjust_BPOD_with_dead_time(BPOD1, ITI1)
BPOD2 = adjust_BPOD_with_dead_time(BPOD2, ITI2)

#%%
Ephys_good = Ephys_wrangle(cluster_info,clust,times,sample_freq,surface_estimate)



#%%
Ephys_binned = bin_Ephys(Ephys_good,bin_size=0.01)

#%%
EPHYS_trimmed1 = Add_Behv(Ephys_binned, BPOD1, raw_BPOD1, ITI1)
EPHYS_trimmed2 = Add_Behv(Ephys_binned, BPOD2, raw_BPOD2, ITI2)

#%%
EPHYS_trimmed1['phase'] = 'cortex'
EPHYS_trimmed2['phase'] = 'Thalamus'

EPhys_combined =  pd.concat([EPHYS_trimmed1, EPHYS_trimmed2],ignore_index=True)

#%%
EPhys_combined.to_csv(basefolder + str('/O1_T.csv'))