# -*- coding: utf-8 -*-
"""
Created on Mon Aug  5 15:15:16 2024

@author: FlorianFreitag
"""

import numpy as np
import pandas as pd
from helpers import*
basefolder ="F:/copydaya\M9_1/M9_SNA-135383_19092024_1_g0_imec1"
mat_name = 'maybeM9_MA_20240919_180430'
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


#%% wrangle BPOD
BPOD = BPOD_wrangle(raw_BPOD,ITI)
BPOD = adjust_BPOD_with_dead_time(BPOD, ITI)

#%% wrangle EPhys

Ephys_good = Ephys_wrangle(cluster_info,clust,times,sample_freq,surface_estimate)

#%% bin EPhys

Ephys_binned = bin_Ephys(Ephys_good,bin_size=0.01)
#%% add BEHV
EPHYS_trimmed = Add_Behv(Ephys_binned, BPOD, raw_BPOD, ITI)

#%% re-add trials
EPHYS_trimmed2 = Add_trial(EPHYS_trimmed)

#%% to csv

EPHYS_trimmed2.to_csv(basefolder + str('/M9_1.csv'))
globals().clear()

#%%

import numpy as np
import pandas as pd
from helpers import*
basefolder = "F:/copydaya/M9_2"
mat_name = 'M9_MA_20240920_151618'


cluster_info= pd.read_csv(basefolder + str('/sorted/phy/cluster_info.tsv'),sep = '\t')
clust = np.array (np.load(basefolder + str('/sorted/phy/spike_clusters.npy')))
times = np.array(np.load(basefolder + str('/sorted/phy/spike_times.npy')))[:,0]
ITI = pd.read_csv(basefolder + str('/Meta/ttl_2.csv'))
raw_BPOD = load_mat(basefolder + '/' + mat_name)
surface_estimate = pd.read_csv(basefolder + str('/Meta/surface_estimate.csv'),header=None).iloc[0, 0]

sample_freq = 30000


ITI = ITI.iloc[:,0]


BPOD = BPOD_wrangle(raw_BPOD,ITI)
BPOD = adjust_BPOD_with_dead_time(BPOD, ITI)



Ephys_good = Ephys_wrangle(cluster_info,clust,times,sample_freq,surface_estimate)


Ephys_binned = bin_Ephys(Ephys_good,bin_size=0.01)

EPHYS_trimmed = Add_Behv(Ephys_binned, BPOD, raw_BPOD, ITI)


EPHYS_trimmed2 = Add_trial(EPHYS_trimmed)


EPHYS_trimmed2.to_csv(basefolder + str('/M9_2.csv'))