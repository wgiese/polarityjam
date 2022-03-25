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

sys.path.append("../feature_extraction/")
import functions

### read parameters

if len(sys.argv) > 0:
    parameter_file = sys.argv[1] 
else:
    print("please provide a parameters file")

print("Read parameters from: %s" % parameter_file)
parameters = functions.read_parameters(parameter_file)

print("Parameters:")
print(parameters)

filename = parameters['image_filename']
input_path = parameters['input_path']
output_path = parameters['output_path']


filepath = input_path + filename
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
print(np.max(cellpose_mask))
#nuclei_mask = functions.get_nuclei_mask(parameters, img, cellpose_mask)
#golgi_mask = functions.get_golgi_mask(parameters, img, cellpose_mask)
properties_df = functions.get_features_from_cellpose_seg(parameters, img, cellpose_mask, filename, output_path)


print(properties_df.head())
properties_df.to_csv(output_path + filename.split(".")[0] + ".csv")


