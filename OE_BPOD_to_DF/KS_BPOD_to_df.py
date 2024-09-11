# -*- coding: utf-8 -*-
"""
Created on Mon Aug  5 15:15:16 2024

@author: FlorianFreitag
"""

import numpy as np
import pandas as pd
from helpers import*

#%%
cluster_group= pd.read_csv('E:/Florian/SNA-127316_2023-12-16_14-41-01/sorted/phy/cluster_group.tsv',sep = '\t')
clust = np.array(np.load('E:/Florian/SNA-127316_2023-12-16_14-41-01/sorted/phy/spike_clusters.npy'))
times = np.array(np.load('E:/Florian/SNA-127316_2023-12-16_14-41-01/sorted//phy/spike_times.npy'))[:,0]
ITI = np.load('E:/Florian/SNA-127316_2023-12-16_14-41-01/Record Node 171/experiment1/recording1/events/Intan_RHD_USB-170.Rhythm Data/TTL/timestamps.npy')
TTL_states = np.load('E:/Florian/SNA-127316_2023-12-16_14-41-01/Record Node 171/experiment1/recording1//events/Intan_RHD_USB-170.Rhythm Data/TTL/states.npy')
raw_BPOD = load_mat('E:/Florian/SNA-127316_2023-12-16_14-41-01/m3_MA_20231216_144158')
sample_nums = np.load('E:/Florian/SNA-127316_2023-12-16_14-41-01/Record Node 171/experiment1/recording1/continuous/Intan_RHD_USB-170.Rhythm Data/sample_numbers.npy')
TTL_sample_times = np.load('E:/Florian/SNA-127316_2023-12-16_14-41-01/Record Node 171/experiment1/recording1/events/Intan_RHD_USB-170.Rhythm Data/TTL/sample_numbers.npy')
sample_freq = 30000

#%%

ITI = get_ITI_starts(TTL_states, TTL_sample_times, sample_nums,sample_freq)


#%%
BPOD = BPOD_wrangle(raw_BPOD,ITI)
BPOD = adjust_BPOD_with_dead_time(BPOD, ITI)

#%%

Ephys_good = Ephys_wrangle(cluster_group,clust,times,ITI,sample_freq)
#%%

Ephys_binned = bin_Ephys(Ephys_good,bin_size=0.01)
#%%
EPHYS_trimmed = Add_Behv(Ephys_binned, BPOD, raw_BPOD, ITI)

#%%
EPHYS_trimmed2 = Add_trial(EPHYS_trimmed)

#%%

EPHYS_trimmed = Raw_to_df_all(cluster_group, clust, times, TTL_states, TTL_sample_times, sample_nums, sample_freq, raw_BPOD)