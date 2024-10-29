# -*- coding: utf-8 -*-
"""
Created on Tue Feb 20 13:24:50 2024

@author: Florian Freitag
"""
import numpy as np
import pandas as pd
from Zeta_helper import*
import itertools
import matplotlib.pyplot as plt



#%%

Data = pd.read_csv("C:/Users/deepl/Desktop/com_Data/wrangled_Data/Zetadata.csv")

#%% for single tries

# =============================================================================
# Data3 = Data.loc[(Data['Day'] == 1) & (Data['Mouse'] == "M7")]
#  
# timess = np.asarray(Data3.loc[Data3[Data3.Trialtype=="Airpuff"].groupby('Trial')['time_bin'].idxmin()]['time_bin'])[np.newaxis].T-0.5
# Spikes = np.asarray(Data3.query('cluster_id==101').loc[:,'time_bin'])[np.newaxis].T
# zp.zetatest(Spikes, timess,boolPlot=True)
#  #%% for single tries
# Data3 = Data.loc[(Data['Day'] == 1) & (Data['Mouse'] == "M7")]
# Zetatest_unit(Data3.loc[Data3['Behv']=='PC'], 'Airpuff2', 176,window=2,lag=0.1,Plot=True)
# Zetatest_unit(Data3.loc[Data3['Behv']=='A2L'], 'Airpuff', 176,window=2,lag=0.1,Plot=True)
# Zetatest_unit(Data3.loc[Data3['Behv']=='W2T'], 'HIT', 176,window=2,lag=0.5,Plot=True)
# Zetatest_unit(Data3.loc[Data3['Behv']=='W2T'], 'Audio', 176,window=2,lag=0.5,Plot=True)
# Zetatest_unit(Data3.loc[Data3['Behv']=='W2T'], 'LeftReward', 176,window=2,lag=0.5,Plot=True)
# =============================================================================
#     
# =============================================================================
#%% Zeta_all
W2T = ZetaMass(Data,'W2T','HIT',lag = 0.01, window = 0.3)
A2L = ZetaMass(Data,'A2L','Airpuff',lag=0.01, window = 0.3)
PC = ZetaMass(Data,'PC','Airpuff2',lag=0.01, window = 0.3)
HIT = ZetaMass(Data,'W2T','HIT',lag=0.3, window = 0.3)

#%%
basefolder = "C:/Users/deepl/Desktop/com_Data/wrangled_Data/Zeta_Tests"
W2T.to_csv(basefolder + str('/W2T.csv'))
A2L.to_csv(basefolder + str('/A2L.csv'))
PC.to_csv(basefolder + str('/PC.csv'))
HIT.to_csv(basefolder + str('/HIT.csv'))