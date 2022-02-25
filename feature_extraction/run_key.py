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
    parameter_file = "base/parameters.yml"
   
print("Read parameters from: %s" % parameter_file)

parameters = functions.read_parameters(parameter_file)
key_file = pd.read_csv(parameters['key_file'])
output_path_base = parameters['output_path']    
summary_df = pd.DataFrame()
summary_ind = 0


for index, row in key_file.iterrows():
        
    input_path = parameters['input_path'] + str(row['folder_name'])
    #subfolder = "/results_%sdyn_%sh_%s/" % (row["shear_stress"], row['flow_exposure_time'], row['treatment'])
    subfolder = str(row["short_name"]) + "/"
    output_path = output_path_base + subfolder    
    if not os.path.exists(output_path):
        os.mkdir(output_path)

    print("input path: %s" % input_path)
    print("output path: %s" % output_path)
    


    file_list = glob.glob(input_path + "*.tif")
    merged_properties_df = pd.DataFrame()
    
    for filepath in file_list:    
        print(filepath)
    

    for ind, filepath in enumerate(file_list):
     
        #output_path = parameters['output_path']
        
        print("processing", filepath.split("/")[-1].split(".")[0])
        filename = filepath.split("/")[-1]
      
        img = functions.read_image(parameters,filepath)
        img_seg = functions.get_image_for_segmentation(parameters,img)
        if parameters["channel_nucleus"] >= 0:
            fig, ax = plt.subplots(1,2)
            ax[0].imshow(img_seg[0,:,:])
            ax[1].imshow(img_seg[1,:,:])
            plt.savefig(output_path + filename + "-seg.png")
        else:
            fig, ax = plt.subplots()
            ax.imshow(img_seg[:,:])
            plt.savefig(output_path + filename + "-seg.png")

        cellpose_mask = functions.get_cellpose_segmentation(parameters, img_seg)
        print("Number of cell labels: ")
        print(len(np.unique(cellpose_mask)))
        #nuclei_mask = functions.get_nuclei_mask(parameters, img, cellpose_mask)
        #golgi_mask = functions.get_golgi_mask(parameters, img, cellpose_mask)
        properties_df = functions.get_features_from_cellpose_seg(parameters, img, cellpose_mask, filename, output_path)
        
        properties_df["condition"] = row["short_name"]

        print(properties_df.head())
        #properties_df.to_csv("image_%s.csv" % ind)
        properties_df.to_csv(output_path + filename.split(".")[0] + ".csv")
        
        if len(merged_properties_df.index) < 10:
            merged_properties_df = properties_df.copy()
        else:
            merged_properties_df = merged_properties_df.append(properties_df, ignore_index=True)

        summary_df.at[summary_ind, "folder_name"] = row["folder_name"] 
        summary_df.at[summary_ind, "short_name"] = row["short_name"] 
        summary_df.at[summary_ind, "filepath"] = filepath 
        summary_df.at[summary_ind, "cell_number"] = len(np.unique(cellpose_mask))

        summary_df.to_csv(output_path_base + "summary_table" + ".csv")
        
        summary_ind += 1
           
        #merged_properties_df.to_csv(output_path + "merged.csv")


