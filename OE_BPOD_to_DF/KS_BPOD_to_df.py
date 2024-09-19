# -*- coding: utf-8 -*-
"""
Created on Mon Aug  5 15:15:16 2024

@author: FlorianFreitag
"""

import numpy as np
import pandas as pd
from helpers import*

#%%
cluster_group= pd.read_csv("D:/M8_SNA-135381_17092024_1_g0_imec1/sorted/phy/cluster_group.tsv",sep = '\t')
clust = np.array(np.load("D://M8_SNA-135381_17092024_1_g0_imec1/sorted/phy/spike_clusters.npy"))
times = np.array(np.load("D:/M8_SNA-135381_17092024_1_g0_imec1/sorted/phy/spike_times.npy"))[:,0]
ITI = pd.read_csv("D:/M8_SNA-135381_17092024_1_g0_imec1/Meta/ttl_2.csv")
raw_BPOD = load_mat("D:/M8_SNA-135381_17092024_1_g0_imec1/M8_MA_20240917_181329.mat")

sample_freq = 30000

#%%

ITI = ITI.iloc[:,1]


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

EPHYS_trimmed2.to_csv("D:/M8_SNA-135381_17092024_1_g0_imec1" + str('/M8_1.csv'))