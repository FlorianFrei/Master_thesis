''' Mouse_Data.py

    Contains the Mouse_Data class that is used for analysing the MatLab data .txt files.
    @Mik Schutte
'''
import numpy as np
import pandas as pd
import os, re, datetime, scipy.io
from scipy.io import matlab
import zetapy as zp


def load_mat(filename): #TODO NOT MY CODE NORA GAVE TO ME
    """
    This function should be called instead of direct scipy.io.loadmat
    as it cures the problem of not properly recovering python dictionaries
    from mat files. It calls the function check keys to cure all entries
    which are still mat-objects
    """

    def _check_vars(d):
        """
        Checks if entries in dictionary are mat-objects. If yes
        todict is called to change them to nested dictionaries
        """
        for key in d:
            if isinstance(d[key], matlab.mio5_params.mat_struct):
                d[key] = _todict(d[key])
            elif isinstance(d[key], np.ndarray):
                d[key] = _toarray(d[key])
        return d
    
    def _todict(matobj):
        """
        A recursive function which constructs from matobjects nested dictionaries
        """
        d = {}
        for strg in matobj._fieldnames:
            elem = matobj.__dict__[strg]
            if isinstance(elem, matlab.mio5_params.mat_struct):
                d[strg] = _todict(elem)
            elif isinstance(elem, np.ndarray):
                d[strg] = _toarray(elem)
            else:
                d[strg] = elem
        return d

    def _toarray(ndarray):
        """
        A recursive function which constructs ndarray from cellarrays
        (which are loaded as numpy ndarrays), recursing into the elements
        if they contain matobjects.
        """
        if ndarray.dtype != 'float64':
            elem_list = []
            for sub_elem in ndarray:
                if isinstance(sub_elem, matlab.mio5_params.mat_struct):
                    elem_list.append(_todict(sub_elem))
                elif isinstance(sub_elem, np.ndarray):
                    elem_list.append(_toarray(sub_elem))
                else:
                    elem_list.append(sub_elem)
            return np.array(elem_list, dtype='object')
        else:
            return ndarray

    data = scipy.io.loadmat(filename, struct_as_record=False, squeeze_me=True)
    return _check_vars(data)

def get_ITI_starts(TTL_states,TTL_sample_times, sample_nums,sample_freq):
    t1 = pd.DataFrame()
    t1['state'] = TTL_states
    t1['time'] = TTL_sample_times
    t2 = t1[t1['state'].isin([1])]
    ITI =(t2['time'] - sample_nums[0]) /  sample_freq
    return ITI


def BPOD_wrangle(raw_BPOD,ITI):
    #takes raw MAtlab list of lists and transforms it into a Dataframe of all trials 
    
    # start and end time in same time as OE data
    start = raw_BPOD['SessionData']['TrialStartTimestamp'] +ITI[0] - raw_BPOD['SessionData']['TrialStartTimestamp'][0] 
    end = raw_BPOD['SessionData']['TrialEndTimestamp'] + ITI[0] - raw_BPOD['SessionData']['TrialStartTimestamp'][0] 

    # calculates the time between Trial end and next trial start
    dead_time=[]
    i=0
    for i in range(len(start) - 1):  # Iterate up to the second-to-last element
        diff = start[i + 1] - end[i]  # Compute the difference between the start of the next trial and the end of the current trial
        dead_time.append(diff)

    # Append 0 for the last trial since there's no "dead time" after the last trial
    dead_time.append(0)
            
    # list to Dataframe
    llist=[]
    i=0
    for i in range(len(raw_BPOD['SessionData']['RawEvents']['Trial'])):
        states = raw_BPOD['SessionData']['RawEvents']['Trial'][i]['States']
        hgh = pd.DataFrame.from_dict(states).transpose()
        hgh['type'] = hgh.index
        hgh['Trial'] = i
        hgh.loc[len(hgh.index)] = [hgh[1].max(),hgh[1].max() + dead_time[i],'dead_time',i]
        llist.append(hgh)
    BPOD = pd.concat(llist)
    BPOD.dropna(inplace=True)
    BPOD.reset_index(inplace=True,drop=True)


    # adds continous time for later steps
    cum =[]
    i= 0
    for i in range(len(BPOD)):
        if i == 0:
            cum.append(BPOD[1][i])
        else:
            temp = BPOD[1][i] - BPOD[0][i]
            diff = temp + cum[i-1]
            cum.append(diff) 
    BPOD['Cont_time'] = cum+ ITI[0]
    BPOD['state_len'] = BPOD[1] - BPOD[0]
    BPOD['Cont_start'] = BPOD['Cont_time']-BPOD['state_len']
    return BPOD

def Ephys_wrangle(cluster_group,clust,times,ITI,sample_freq):
    #takes raw OE vectos and turns them into a Dataframe
    #selects only clusters that have the PHY good label
    good_cluster = cluster_group.query('group=="good"').loc[:,['cluster_id']]
    Ephys_raw = pd.DataFrame({'cluster_id': clust, 'times':times}, columns=['cluster_id', 'times'])
    Ephys_raw = Ephys_raw.assign(seconds = Ephys_raw['times']/sample_freq)
    Ephys_good = (
        Ephys_raw.merge(good_cluster, 
                  on=['cluster_id'],
                  how='left', 
                  indicator=True)
        .query('_merge == "both"')
        .drop(columns='_merge')
    )
    return Ephys_good

def bin_Ephys(Ephys_good,bin_size=0.01):

    max_seconds = Ephys_good['seconds'].max()
    bin_edges = np.arange(0, max_seconds + bin_size, bin_size)

    # Create a new column with bin labels representing the end of each bin
    Ephys_good['time_bin'] = pd.cut(Ephys_good['seconds'], bins=bin_edges, right=False, labels=bin_edges[1:])

    # Group by cluster_id and time_bin, then count the number of events
    result = Ephys_good.groupby(['cluster_id', 'time_bin']).size().reset_index(name='event_count')

    # If you want to fill in missing bins with zeros
    all_combinations = pd.MultiIndex.from_product([Ephys_good['cluster_id'].unique(), bin_edges[1:]], names=['cluster_id', 'time_bin'])
    result = result.set_index(['cluster_id', 'time_bin']).reindex(all_combinations, fill_value=0).reset_index()

    return(result)


def Add_Behv(Ephys_binned,BPOD,raw_BPOD,ITI):
    
    start = raw_BPOD['SessionData']['TrialStartTimestamp'] +ITI[0] - raw_BPOD['SessionData']['TrialStartTimestamp'][0] 
    end = raw_BPOD['SessionData']['TrialEndTimestamp'] + ITI[0] - raw_BPOD['SessionData']['TrialStartTimestamp'][0] 
    
    EPHYS_trimmed = Ephys_binned.loc[(Ephys_binned['time_bin'] <= end[-1]) & (Ephys_binned['time_bin'] >= start[0])]
    EPHYS_trimmed = EPHYS_trimmed.sort_values(by='time_bin').reset_index(drop=True)
    i = 0
    j = 0
    trialtype = []
    
    while i < len(EPHYS_trimmed) and j < len(BPOD):
        if BPOD['Cont_start'][j] <= EPHYS_trimmed.iloc[i, 1] <= BPOD['Cont_time'][j]:
            trialtype.append(BPOD['type'][j])
            i += 1
        else:
            j += 1
            print(j)
    EPHYS_trimmed['Trialtype'] = trialtype 
    return(EPHYS_trimmed)

def Add_trial(EPHYS_trimmed):
    k=0
    trial =[]
    for i in range(len(EPHYS_trimmed)):
        # Check the condition for ITI
        if EPHYS_trimmed.iloc[i, 3] != 'ITI' and (i + 1 < len(EPHYS_trimmed) and EPHYS_trimmed.iloc[i + 1, 3] == 'ITI'):
            trial.append(k)
            print(k)
            k += 1
        else:
            trial.append(k)

    # Ensure the lengths match before assigning the trial list to the DataFrame
    if len(trial) == len(EPHYS_trimmed):
        EPHYS_trimmed['trial'] = trial
    else:
        # Handle the case where lengths do not match
        raise ValueError("Length of `trial` does not match the length of `EPHYS_trimmed`")
    return(EPHYS_trimmed)


def Raw_to_df_all(cluster_group,clust,times,TTL_states,TTL_sample_times, sample_nums,sample_freq,raw_BPOD,bin_size=0.01):
    ITI = get_ITI_starts(TTL_states, TTL_sample_times, sample_nums,sample_freq)
    BPOD = BPOD_wrangle(raw_BPOD,ITI)
    Ephys_good = Ephys_wrangle(cluster_group,clust,times,ITI,sample_freq)
    Ephys_binned = bin_Ephys(Ephys_good,bin_size=0.01)
    EPHYS_trimmed = Add_Behv(Ephys_binned, BPOD, raw_BPOD, ITI)
    EPHYS_trimmed2 = Add_trial(EPHYS_trimmed)
    return EPHYS_trimmed2













def Zetatest_unit(Data, State,unit,Plot=False, window = 0.4,lag=0.2):
    timess = np.asarray(Data.loc[Data[Data.Trialtype==str(State)].groupby('Trial')['time_bin'].idxmin()]['time_bin'])[np.newaxis].T-lag
    Spikes = np.asarray(Data.query('cluster_id=='+ str(unit)).loc[:,'time_bin'])[np.newaxis].T
    if Plot==False:
        dblZetaP = zp.zetatest(Spikes, timess,dblUseMaxDur=window, intResampNum = 5000)[0] 
    else:
        dblZetaP = zp.zetatest(Spikes, timess,dblUseMaxDur=window,boolPlot=True,intResampNum = 5000)[0] 
    return dblZetaP,unit


def Zetatest_all(Data,State,window=0.4,lag=0.2):
    if State in list(Data['Trialtype']):
        p_vals=[]
        for i in Data['cluster_id'].unique():
            temp = Zetatest_unit(Data,State,i,Plot=False,window= window,lag=lag)
            p_vals.append(temp)
        p_vals = pd.DataFrame.from_records(p_vals,columns=['p_val', 'unit'])
        return p_vals
    else: print('FU')
    
