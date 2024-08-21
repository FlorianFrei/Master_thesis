# -*- coding: utf-8 -*-
"""
Created on Mon Aug 19 16:25:11 2024

@author: Florian Freitag
"""

import os
import pandas as pd 

def resturcture_csv(root_dir):
    """
    Renames CSV files in a directory structure, including full path in the new filename.

    Args:
        root_dir: The root directory containing the CSV files.
    """

    for root, dirs, files in os.walk(root_dir):
        for file in files:
            if file.endswith(".csv"):
                file_path = os.path.join(root, file)
                
                df = pd.read_csv(file_path, header=None)  # Load without header to treat all rows equally

                # 1. Delete the first row
                df = df.drop(index=0).reset_index(drop=True)
                
                # 2. Combine the second and third rows to make the new header
                new_header = df.iloc[0].fillna('') + '_' + df.iloc[1].fillna('')  # Combine rows 0 and 1
                df = df[2:]  # Drop the now combined rows
                df.columns = new_header  # Set the new header
                
                df.to_csv(file_path, index=False)
                print(f"Processed and saved: {file_path}")
               
resturcture_csv("D:/try")
