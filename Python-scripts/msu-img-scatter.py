#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  7 22:22:24 2020

@author: Warren.Sink
"""

import os
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.pyplot import imread
from PIL import Image

os.chdir("/Users/Warren.Sink/github/facial-polyphenism/CSVs")

df = pd.read_csv('msu_extract_nonnosenorm.csv')
df_y = pd.read_csv('msu_extract_nonnosenorm_y.csv')
df_x = pd.read_csv('msu_extract_nonnosenorm_x.csv')

os.chdir("/Users/Warren.Sink/Desktop/OneDrive-Van-Andel-Institute/msu-images-extract")
source_dir = "/Users/Warren.Sink/Desktop/OneDrive-Van-Andel-Institute/msu-images-extract"
target_dir = "/Users/Warren.Sink/Desktop/OneDrive-Van-Andel-Institute/msu-images-scatter"
directory = os.fsencode("/Users/Warren.Sink/Desktop/OneDrive-Van-Andel-Institute/msu-images-extract")

source_dir = "/Users/Warren.Sink/Desktop/OneDrive-Van-Andel-Institute/msu-images-extract"
filenames = []
for idx, file in enumerate(os.listdir(source_dir)):
    filename = os.fsdecode(file)
    if filename.endswith(".jpg"): 
        filenames.append(filename)
        continue
    else:
        continue    
    
filenames.sort()

def scatterp_img(filenames, landmark_x, landmark_y):
    for idx, IMAGE in enumerate(filenames):
     counter = idx + 1
     img = imread(IMAGE)
     im = Image.open(IMAGE)
     width, height = im.size
     plt.scatter(x = landmark_x.iloc[:,counter], y = landmark_y.iloc[:,counter], c='b', s=15)
     plt.imshow(img, zorder=0, extent=[0, width, 0, height])
     plt.savefig(fname = target_dir + "/" + IMAGE[2:6] + ".pdf")
     plt.show()
     














