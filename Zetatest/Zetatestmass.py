# -*- coding: utf-8 -*-
"""
Created on Tue Feb 20 13:24:50 2024

@author: Florian Freitag
"""
import numpy as np
import pandas as pd
from load_mat import*
import itertools
import matplotlib.pyplot as plt



#%%

Data = pd.read_csv('C:/Users/deepl/Desktop/ffOE/master_code/data/Zetadata2.csv')

#%%
Data3 = Data.loc[(Data['Day'] == 1) & (Data['Mouse'] == "M4")]


timess = np.asarray(Data3.loc[Data3[Data3.Trialtype=="Airpuff"].groupby('Trial')['time_bin'].idxmin()]['time_bin'])[np.newaxis].T-0.5
Spikes = np.asarray(Data3.query('cluster_id==10').loc[:,'time_bin'])[np.newaxis].T
   
zp.zetatest(Spikes, timess,boolPlot=True)
#%%
Data3 = Data.loc[(Data['Day'] == 3) & (Data['Mouse'] == "M6")]
Zetatest_unit(Data3.loc[Data3['Behv']=='PC'], 'Airpuff2', 19,window=None,lag=0.1,Plot=True)
Zetatest_unit(Data3.loc[Data3['Behv']=='A2L'], 'Airpuff', 19,window=None,lag=0.1,Plot=True)
Zetatest_unit(Data3.loc[Data3['Behv']=='W2T'], 'HIT', 19,window=1,lag=0.5,Plot=True)
#Zetatest_unit(Data3.loc[Data3['Behv']=='W2T'], 'Audio', 10,window=None,lag=0.5,Plot=True)
#Zetatest_unit(Data3.loc[Data3['Behv']=='W2T'], 'LeftReward', 10,window=None,lag=0.5,Plot=True)

#%%


  
    
#%%
output_data1 = []

for iteration in range(1, 2):  
    print(iteration)
    # Generate all unique combinations of 'Mouse' and 'Day' from the available data
    unique_combinations = list(itertools.product(Data['Mouse'].unique(), Data['Day'].unique()))
    for Mouse, Day in unique_combinations:
        Data2 = Data.loc[(Data['Day'] == Day) & (Data['Mouse'] == Mouse)]
        ZTd1 = Zetatest_all(Data2.loc[Data2['Behv'] == 'W2T'], 'HIT', window=None, lag=0.5)
        if ZTd1 is not None:
            for index, row in ZTd1.iterrows():
                output_data1.append([iteration, Day, Mouse, row['p_val'], row['unit']])
        else:
            print(f"No data found for Day: {Day}, Mouse: {Mouse}")

W2T = pd.DataFrame(output_data1, columns=['Iteration', 'Day', 'Mouse', 'p_val', 'unit'])
W2T.to_csv('C:/Users/deepl/Desktop/ffOE/master_code/data/koen3/W2T.csv')   




output_data1 = []

for iteration in range(1, 2):  
    print(iteration)
    # Generate all unique combinations of 'Mouse' and 'Day' from the available data
    unique_combinations = list(itertools.product(Data['Mouse'].unique(), Data['Day'].unique()))
    for Mouse, Day in unique_combinations:
        Data2 = Data.loc[(Data['Day'] == Day) & (Data['Mouse'] == Mouse)]
        ZTd1 = Zetatest_all(Data2.loc[Data2['Behv'] == 'A2L'], 'Airpuff', window=None, lag=0.1)
        if ZTd1 is not None:
            for index, row in ZTd1.iterrows():
                output_data1.append([iteration, Day, Mouse, row['p_val'], row['unit']])
        else:
            print(f"No data found for Day: {Day}, Mouse: {Mouse}")

A2L = pd.DataFrame(output_data1, columns=['Iteration', 'Day', 'Mouse', 'p_val', 'unit'])
A2L.to_csv('C:/Users/deepl/Desktop/ffOE/master_code/data/koen3/A2L.csv')  



output_data1 = []

for iteration in range(1, 2):  
    print(iteration)
    # Generate all unique combinations of 'Mouse' and 'Day' from the available data
    unique_combinations = list(itertools.product(Data['Mouse'].unique(), Data['Day'].unique()))
    for Mouse, Day in unique_combinations:
        Data2 = Data.loc[(Data['Day'] == Day) & (Data['Mouse'] == Mouse)]
        ZTd1 = Zetatest_all(Data2.loc[Data2['Behv'] == 'PC'], 'Airpuff2', window=None, lag=0.1)
        if ZTd1 is not None:
            for index, row in ZTd1.iterrows():
                output_data1.append([iteration, Day, Mouse, row['p_val'], row['unit']])
        else:
            print(f"No data found for Day: {Day}, Mouse: {Mouse}")

PC = pd.DataFrame(output_data1, columns=['Iteration', 'Day', 'Mouse', 'p_val', 'unit'])
PC.to_csv('C:/Users/deepl/Desktop/ffOE/master_code/data/koen3/PC.csv')   

