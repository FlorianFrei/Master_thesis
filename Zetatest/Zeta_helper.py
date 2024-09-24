# -*- coding: utf-8 -*-
"""
Created on Tue Sep 24 17:43:06 2024

@author: Florian Freitag
"""

import numpy as np
import pandas as pd
import zetapy as zp

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
    