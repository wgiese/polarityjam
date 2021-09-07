import pandas as pd
import yaml
import numpy as np
import scipy.ndimage
import time, os, sys
import skimage.io
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['figure.dpi'] = 300
from cellpose import utils
from matplotlib import cm
from skimage.filters import threshold_otsu
from skimage.measure import label, regionprops, regionprops_table
import math
from scipy.ndimage import gaussian_filter
import sys
from cellpose import models, io, plot
import os
import glob


sys.path.append("../feature_extraction/")
import functions

if len(sys.argv) > 1:
    parameter_file = sys.argv[1]
else:
    parameter_file = "../feature_extraction/local/parameters.yml"
   
print("Read parameters from: %s" % parameter_file)

parameters = functions.read_parameters(parameter_file)
input_path = parameters["input_path"]



file_list = glob.glob(input_path + "*.tif")
merged_properties_df = pd.DataFrame()

for ind, filepath in enumerate(file_list):
 
    output_path = parameters['output_path']
    
    print("processing", filepath.split("/")[-1].split(".")[0])
    filename = filepath.split("/")[-1]
  
    img = functions.read_image(parameters,filepath)
    img_seg = functions.get_image_for_segmentation(parameters,img)

    fig, ax = plt.subplots(1,2)
    ax[0].imshow(img_seg[0,:,:])
    ax[1].imshow(img_seg[1,:,:])
    plt.savefig(output_path + filename + "-seg.png")

    cellpose_mask = functions.get_cellpose_segmentation(parameters, img_seg)
    print("Number of cell labels: ")
    print(np.max(cellpose_mask))
    #nuclei_mask = functions.get_nuclei_mask(parameters, img, cellpose_mask)
    #golgi_mask = functions.get_golgi_mask(parameters, img, cellpose_mask)
    properties_df = functions.get_features_from_cellpose_seg(parameters, img, cellpose_mask, filename)
    
    
    print(properties_df.head())
    #properties_df.to_csv("image_%s.csv" % ind)
    properties_df.to_csv(output_path + filename.split(".")[0] + ".csv")
    
    if len(merged_properties_df.index) < 10:
        merged_properties_df = properties_df.copy()
    else:
        merged_properties_df.append(properties_df, ignore_index=True)
        
    merged_properties_df.to_csv(output_path + "merged.csv")


