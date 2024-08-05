

# -*- coding: utf-8 -*-
"""
Created on Fri Aug  2 17:12:42 2024

@author: AG_Larkum
"""


import os

import deeplabcut

def getsubfolders(folder):
    ''' returns list of subfolders '''
    return [os.path.join(folder,p) for p in os.listdir(folder) if os.path.isdir(os.path.join(folder,p))]
    


shuffle=1


projectpath='C:/Users/AG_Larkum/Documents/Liv/DLC_whisker/whisker_rec-liv-2024-04-12'
config=os.path.join(projectpath,'config.yaml')

basepath='C:/Users/AG_Larkum/Desktop/try' #data'

'''

Imagine that the data (here: videos of 3 different types) are in subfolders:
    /January/January29 ..
    /February/February1
    /February/February2
    
    etc.

'''

subfolders=getsubfolders(basepath)
for subfolder in subfolders: #this would be January, February etc. in the upper example
    print("Starting analyze data in:", subfolder)
    subsubfolders=getsubfolders(subfolder)
    for subsubfolder in subsubfolders: #this would be Febuary1, etc. in the upper example...
        print("Starting analyze data in:", subsubfolder)
        for vtype in ['.mp4','.mp4v','.m4v']:
            deeplabcut.analyze_videos(config,[subsubfolder],shuffle=shuffle,videotype=vtype,save_as_csv=True)
            
        