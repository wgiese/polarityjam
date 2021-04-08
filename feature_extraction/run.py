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

### read parameters

if len(sys.argv) > 0:
    parameter_file = sys.argv[1] 
else:
    print("please provide a parameters file")

def read_parameters(parameter_file):

    parameters = dict()

    with open("base/parameters.yml") as file:
        parameters = yaml.load(file, Loader=yaml.FullLoader)
    with open(parameter_file) as file:
        parameters_local = yaml.load(file, Loader=yaml.FullLoader)

    # overwrite global parameters with local setting
    for key in parameters_local:
        parameters[key] = parameters_local[key]
    
    return parameters

parameters = read_parameters(parameter_file)

print("Parameters:")
print(parameters)

filename = parameters['image_path']
img = skimage.io.imread(filename)
output_path = parameters['output_path']
output_filename = parameters["output_filename"]
output_filepath = output_path + output_filename

# extract channels

im_junction = img[:,:,parameters["channel_junction"]]
im_nucleus = img[:,:,parameters["channel_nucleus"]]
im_golgi = img[:,:,parameters["channel_golgi"]]

# use junction and nucleus channel for cellpose segmentation 
im_seg = np.array([im_junction, im_nucleus])


model = models.Cellpose(gpu=True, model_type='cyto')

channels = [1,2]

masks, flows, styles, diams = model.eval(im_seg, diameter=100, channels=channels)

io.masks_flows_to_seg(im_seg , masks, flows, diams, output_filename, channels)

#read cellpose segmentation or feature extraction

cellpose_seg = np.load(output_filename + "_seg.npy", allow_pickle=True)
mask = cellpose_seg.item()['masks']
for key in cellpose_seg.item():
    print(key)

img_golgi_blur = gaussian_filter(img[:,:,0], sigma= 3)
img_nuclei_blur = gaussian_filter(img[:,:,1], sigma=3)

nuclei_mask = np.where(img_nuclei_blur > threshold_otsu(img_nuclei_blur), True, False)
golgi_mask = np.where(img_golgi_blur > threshold_otsu(img_golgi_blur), True, False)

nuclei_label = nuclei_mask*mask
golgi_label = golgi_mask*mask
print(np.max(mask))

# feature extraction

single_cell_props = pd.DataFrame()
counter = 0

for label in range(1,np.max(mask)-1):
    
    single_cell_mask = np.where(mask ==label, 1, 0)
    single_nucleus_mask = np.where(nuclei_label==label, 1, 0)
    single_golgi_mask = np.where(golgi_label ==label, 1, 0)    
    
    area_cell_px2 = len(single_cell_mask[single_cell_mask==1])
    area_golgi_px2 = len(single_golgi_mask[single_golgi_mask==1])
    area_nucleus_px2 = len(single_nucleus_mask[single_nucleus_mask==1])    
    
    if (area_nucleus_px2) < 10 or (area_golgi_px2 < 10):
        continue
                
    #print(len(single_cell_mask[single_cell_mask==1]))
    #print(len(single_nucleus_mask[single_nucleus_mask==1]))
    #print(len(single_golgi_mask[single_golgi_mask==1]))
      
    regions = regionprops(single_cell_mask)
    for props in regions:
        x_cell, y_cell = props.centroid
        orientation = props.orientation
        minor_axis_length = props.minor_axis_length
        major_axis_length = props.major_axis_length
                
    single_cell_props.at[counter, "label"] = label
    single_cell_props.at[counter, "X_cell"] = x_cell
    single_cell_props.at[counter, "Y_cell"] = y_cell
    single_cell_props.at[counter, "shape_orientation"] = x_cell
    single_cell_props.at[counter, "major_axis_length"] = major_axis_length
    single_cell_props.at[counter, "minor_axis_length"] = minor_axis_length
    
    regions = regionprops(single_nucleus_mask)
    for props in regions:
        x_nucleus, y_nucleus = props.centroid
    
    single_cell_props.at[counter, "X_nuc"] = x_nucleus
    single_cell_props.at[counter, "Y_nuc"] = y_nucleus
    
    regions = regionprops(single_golgi_mask)
    for props in regions:
        x_golgi, y_golgi = props.centroid
    
    single_cell_props.at[counter, "X_golgi"] = x_golgi
    single_cell_props.at[counter, "Y_golgi"] = y_golgi
    
    distance2 = (x_golgi - x_nucleus)**2
    distance2 += (y_golgi - y_nucleus)**2
    
    single_cell_props.at[counter, "distance"] = np.sqrt(distance2)
    
    vec_x = x_golgi - x_nucleus
    vec_y = y_golgi - y_nucleus
    angle_rad = np.arctan2(vec_y, vec_x)
    
    
    
    counter += 1
    
single_cell_props.to_csv(output_filepath + "_.csv") 

# visualiztion of nuclei-golgi vectors

fig, ax = plt.subplots()

nuclei_mask_ = np.where(nuclei_mask == True, 70, 0)
golgi_mask_ = np.where(golgi_mask == True, 20, 0)

ax.imshow(im_junction, cmap=plt.cm.gray, alpha = 1.0)
ax.imshow(np.ma.masked_where(nuclei_mask_ == 0, nuclei_mask_),  plt.cm.gist_rainbow, vmin=0, vmax=100, alpha=1.0)

ax.imshow(np.ma.masked_where(golgi_mask_ == 0, golgi_mask_),  plt.cm.gist_rainbow, vmin=0, vmax=100, alpha=1.0)

for index, row in single_cell_props.iterrows():
    ax.plot( row["Y_nuc"], row["X_nuc"], '.g', markersize=1)
    ax.plot( row["Y_golgi"], row["X_golgi"], '.m', markersize=1)
    ax.arrow(row["Y_nuc"], row["X_nuc"], row["Y_golgi"]- row["Y_nuc"],row["X_golgi"]- row["X_nuc"], color = 'red', width = 2)

ax.set_xlim(0,im_junction.shape[0])
ax.set_ylim(0,im_junction.shape[1])
plt.savefig(output_filepath + "_nuclei_golgi_vector.pdf")
plt.savefig(output_filepath + "_nuclei_golgi_vector.png")

# visualiztion of nuclei-golgi vectors

fig, ax = plt.subplots()
ax.imshow(im_junction, cmap=plt.cm.gray)

regions = regionprops(mask)
for props in regions:
    y0, x0 = props.centroid
    orientation = props.orientation
    x1 = x0 + math.cos(orientation) * 0.5 * props.minor_axis_length
    y1 = y0 - math.sin(orientation) * 0.5 * props.minor_axis_length
    x2 = x0 - math.sin(orientation) * 0.5 * props.major_axis_length
    y2 = y0 - math.cos(orientation) * 0.5 * props.major_axis_length

    ax.plot((x0, x1), (y0, y1), '--r', linewidth=0.5)
    ax.plot((x0, x2), (y0, y2), '--r', linewidth=0.5)
    ax.plot(x0, y0, '.b', markersize=5)

ax.set_xlim(0,im_junction.shape[0])
ax.set_ylim(0,im_junction.shape[1])
plt.savefig(output_filepath + "_cellshape_orientation.pdf")
plt.savefig(output_filepath + "_cellshape_orientation.png")


