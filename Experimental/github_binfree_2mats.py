# -*- coding: utf-8 -*-
"""
Created on Wed Dec 18 17:12:33 2024

@author: deepl
"""

import os



import numpy as np
import pandas as pd
from helpers import*
basefolder ="D:/Florian_paper/Florian/Data/Opto_2/O5_2"
mat_name1 = 'MLA5686_OptoAir_20241219_130645'
mat_name2 = 'MLA5686_OptoAir_20241219_131846'

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
proceed1 = check_trialnumber_matches(ITI1, raw_BPOD1)
proceed2 = check_trialnumber_matches(ITI2, raw_BPOD2)
#%%
BPOD1 = BPOD_wrangle(raw_BPOD1,ITI1, proceed1)
BPOD2 = BPOD_wrangle(raw_BPOD2,ITI2, proceed2)

#%%
BPOD1 = adjust_BPOD_with_dead_time(BPOD1, ITI1)
BPOD2 = adjust_BPOD_with_dead_time(BPOD2, ITI2)

#%% wrangle EPhys
Ephys_good = Ephys_wrangle(cluster_info,clust,times,sample_freq,surface_estimate)

#%%

filtered_df = BPOD[BPOD['Trialtype'].isin([1, 4])]
filtered_df = filtered_df[filtered_df['type'].isin(['AirTop_noOpto', 'AirTop_Opto'])]
filtered_df['event_label'] = filtered_df['Trialtype'].replace({1: 'Airpuff', 4: 'Opto_air'})

#%%

result_rows = []
for _, row in filtered_df.iterrows():
    trial = row['Trial']
    cont_start = row['Cont_start']
    event_label = row['event_label']

    # Get rows in Ephys_good where seconds are within ±1 of Cont_start
    matches = Ephys_good[(Ephys_good['seconds'] >= cont_start - 1) & (Ephys_good['seconds'] <= cont_start + 1)].copy()

    # Add the Trial, Cont_start, and event_label columns
    matches['Trial'] = trial
    matches['time_from_event'] = matches['seconds'] - cont_start
    matches['event_time'] = cont_start
    matches['event_label'] = event_label

    result_rows.append(matches)

# Concatenate all matched rows into a new DataFrame
result_df = pd.concat(result_rows, ignore_index=True)


#%% save

result_csv_path = "D:/Florian_paper/Florian/Aligned_Data/otpo/wrangled_Data/Binfree/" +str(os.path.basename(basefolder)) +'.csv'
result_df.to_csv(result_csv_path, index=False)

