# -*- coding: utf-8 -*-
"""
Created on Fri Aug  2 16:45:29 2024

@author: AG_Larkum
"""

import os
import re

def rename_videos(root_dir):
    """
    Renames videos in a directory structure with multiple nested folders, including full path in the new filename.

    Args:
        root_dir: The root directory containing the video files.
    """

    for root, dirs, files in os.walk(root_dir):
        for file in files:
            if file.endswith(".mp4v"):
                video_number = re.search(r'\d+(?=\.mp4)', file).group()  # Capture digits before .mp4 extension

                video_number = video_number.zfill(2)  # Pad with zeros

                # Construct the new filename with full path
                relative_path = os.path.relpath(root, root_dir)
                new_name_parts = relative_path.split(os.sep) + [f"{video_number}.mp4"]
                new_path = os.path.join(root_dir, *new_name_parts)

                # Check for existing file and increment counter if necessary
                counter = 1
                while os.path.exists(new_path):
                    new_name_parts[-2] += f"_{counter}"
                    new_path = os.path.join(root_dir, *new_name_parts)
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
    rename_videos(root_directory)
