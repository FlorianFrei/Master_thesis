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
basefolder ="D:/Florian_paper/Florian/Data/Opto_2/M5686_g0"
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

EPHYS_trimmed.to_csv(basefolder + str('/O5_1.csv'))



#%%

ITI_idx = 0
# List to store new rows
new_rows = []

# Iterate through BPOD and identify ITI mismatches
for i in range(len(BPOD)):  # Start from 1 because we need to check the previous row
    row = BPOD.iloc[i]
    
    if row['type'] == 'ITI':
        # Find the true ITI start time
        true_iti = ITI.iloc[ITI_idx]
        ITI_idx += 1
        
        # If BPOD cont_start is earlier than true ITI, we have dead time
        if row['Cont_start'] < true_iti and BPOD.iloc[i - 1]['type'] == 'dead_time':
            dead_time2 = true_iti - row['Cont_start']
            
            # Add a 'dead_time2' row with cont_start from previous row's cont_time
            prev_row = BPOD.iloc[i - 1]
            new_row = {
                'type': 'dead_time2',
                'Cont_start': prev_row['Cont_time'],
                'Cont_time': true_iti,
                'Trial': prev_row['Trial'],   # Carry over trial from the previous row
                0: prev_row[1],           # '0' is '1' from the previous row
                1: prev_row[1] + (true_iti - prev_row['Cont_time']),  # Calculate '1' using state_len
                'state_len': true_iti - prev_row['Cont_time']  # Calculate the new state_len
            }
            new_rows.append((i, new_row))
            
            # Adjust ITI row's Cont_start and Cont_time
            BPOD.at[i, 'Cont_start'] = true_iti
            BPOD.at[i, 'Cont_time'] += dead_time2  # Adjust the Cont_time based on lost time
            BPOD.at[i, 'state_len'] = BPOD.at[i, 'Cont_time'] - BPOD.at[i, 'Cont_start']  # Update state_len
   
            # Now adjust the following rows
            for j in range(i + 1, len(BPOD)):
                BPOD.at[j, 'Cont_start'] += dead_time2
                BPOD.at[j, 'Cont_time'] += dead_time2
                BPOD.at[j, 'state_len'] = BPOD.at[j, 'Cont_time'] - BPOD.at[j, 'Cont_start']  # Update state_len
    elif row['Cont_start'] > true_iti and BPOD.iloc[i - 1]['type'] == 'dead_time':
            print(f"Warning: BPOD Cont_start ({row['Cont_start']}) is later than true ITI ({true_iti}).")
   
            # Add a 'dead_time2' row with negative state_len
            prev_row = BPOD.iloc[i - 1]
            dead_time2 = row['Cont_start'] - true_iti
            new_row = {
                'type': 'dead_time2',
                'Cont_start': true_iti,
                'Cont_time': row['Cont_start'],
                'Trial': prev_row['Trial'],
                0: prev_row[1],
                1: prev_row[1] - dead_time2,
                'state_len': -dead_time2
            }
            new_rows.append((i, new_row))
   
            # Adjust the ITI row's Cont_start to align with true ITI
            row_shift = row['Cont_start'] - true_iti
            BPOD.at[i, 'Cont_start'] = true_iti
            BPOD.at[i, 'Cont_time'] -= row_shift
            BPOD.at[i, 'state_len'] = BPOD.at[i, 'Cont_time'] - BPOD.at[i, 'Cont_start']
   
            # Update subsequent rows
            for j in range(i + 1, len(BPOD)):
                BPOD.at[j, 'Cont_start'] -= row_shift
                BPOD.at[j, 'Cont_time'] -= row_shift
                BPOD.at[j, 'state_len'] = BPOD.at[j, 'Cont_time'] - BPOD.at[j, 'Cont_start']

# Insert new 'dead_time2' rows into the BPOD dataframe
shift = 0  # Track how much we are shifting the index as rows are inserted
for idx, new_row in new_rows:
    idx += shift  # Adjust index by the shift since we're adding rows
    BPOD = pd.concat([BPOD.iloc[:idx], pd.DataFrame([new_row]), BPOD.iloc[idx:]]).reset_index(drop=True)
    shift += 1  # Increment shift since we inserted a row
print(f"Mean value of dead_time2: {BPOD[BPOD['type'] == 'dead_time2']['state_len'].mean()}")
print(f"Max absolute value of dead_time2: {BPOD[BPOD['type'] == 'dead_time2']['state_len'].abs().max()}")
filtered_data = BPOD[BPOD['type'] == 'dead_time2']['state_len'].abs()
max_excluded = filtered_data.nlargest(2).iloc[-1]  # Get the second largest value

print(f"Max absolute value of dead_time2 (excluding the highest): {max_excluded}")
