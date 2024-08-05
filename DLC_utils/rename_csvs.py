# -*- coding: utf-8 -*-
"""
Created on Mon Aug  5 12:06:49 2024

@author: AG_Larkum
"""

import os
import re

def rename_csv_files(root_dir):
    """
    Renames CSV files in a directory structure, including full path in the new filename.

    Args:
        root_dir: The root directory containing the CSV files.
    """

    for root, dirs, files in os.walk(root_dir):
        for file in files:
            if file.endswith(".csv"):
                # Extract video number
                match = re.search(r'\d{2,3}DLC_', file)
                if match:
                    video_number = match.group()[:-4]  # Remove 'DLC_' from the match

                    # Construct the new filename with full path
                    relative_path = os.path.relpath(root, root_dir).replace(os.sep, '_')
                    new_name = f"{relative_path}_{video_number}.csv"
                    new_path = os.path.join(root, new_name)

                    # Check for existing file and increment counter if necessary
                    counter = 1
                    while os.path.exists(new_path):
                        new_name = f"{relative_path}_{video_number}_{counter}.csv"
                        new_path = os.path.join(root, new_name)
                        counter += 1

                    old_path = os.path.join(root, file)
                    print(f"Old path: {old_path}")
                    print(f"New path: {new_path}")

                    if os.path.exists(old_path):
                        try:
                            os.rename(old_path, new_path)
                            print(f"Renamed {old_path} to {new_path}")
                        except OSError as e:
                            print(f"Error renaming {old_path}: {e}")
                    else:
                        print(f"File not found: {old_path}")

if __name__ == "__main__":
    root_directory = "C:/Users/AG_Larkum/Desktop/try"  # Replace with your root directory
    rename_csv_files(root_directory)
