import os
import yaml
import skimage.io
import numpy as np
import cellpose.utils
import cellpose.io
import cellpose.models
import scipy.ndimage as ndi
import skimage.filters 
import skimage.measure
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl

def read_parameters(parameter_file):

    parameters = dict()

    if os.path.exists("../feature_extraction/base/parameters.yml"):
        with open("../feature_extraction/base/parameters.yml") as file:
            parameters = yaml.load(file, Loader=yaml.FullLoader)

    #with open("base/parameters.yml") as file:
    #    parameters = yaml.load(file, Loader=yaml.FullLoader)
    with open(parameter_file) as file:
        parameters_local = yaml.load(file, Loader=yaml.FullLoader)

    # overwrite global parameters with local setting
    for key in parameters_local:
        parameters[key] = parameters_local[key]
    
    return parameters
    
def read_image(parameters, filename):

    img = skimage.io.imread(filename)
    return img

def get_image_for_segmentation(parameters, img):

    print(parameters["channel_junction"])
    print(img.shape)
    im_junction = img[:,:,int(parameters["channel_junction"])]
    im_nucleus = img[:,:,int(parameters["channel_nucleus"])]
    im_seg = np.array([im_junction, im_nucleus])

    return im_seg

def get_cellpose_segmentation(parameters, im_seg):

    model = cellpose.models.Cellpose(gpu=True, model_type='cyto')

    channels = [1,2]

    masks, flows, styles, diams = model.eval(im_seg, diameter=parameters["estimated_cell_diameter"], channels=channels)  

    return masks

def get_nuclei_mask(parameters, img, cellpose_mask):

    img_nuclei_blur = ndi.gaussian_filter(img[:,:,parameters["channel_nucleus"]], sigma=3)

    nuclei_mask = np.where(img_nuclei_blur > skimage.filters.threshold_otsu(img_nuclei_blur), True, False)

    nuclei_label = nuclei_mask*cellpose_mask

    return nuclei_label

def get_nuclei_cellpose(parameters, img, cellpose_mask):

    model = cellpose.models.Cellpose(gpu=True, model_type='nuclei')

    channels = [0,0]

    masks, flows, styles, diams = model.eval(img, diameter=parameters["estimated_cell_diameter"], channels=channels)
    #masks, flows, styles, diams = model.eval(img, channels=channels)

    nuclei_mask = np.where(masks > 0, True, False)

    nuclei_label = nuclei_mask*cellpose_mask

    return nuclei_label

def get_golgi_mask(parameters, img, cellpose_mask):

    img_golgi_blur = ndi.gaussian_filter(img[:,:,parameters["channel_golgi"]], sigma=3)
    
    golgi_mask = np.where(img_golgi_blur > skimage.filters.threshold_otsu(img_golgi_blur), True, False)

    golgi_label = golgi_mask*cellpose_mask
    
    return golgi_label


def get_features_from_cellpose_seg(parameters, img, cell_mask, filename):

    #nuclei_mask = get_nuclei_mask(parameters, img, cell_mask)
    nuclei_mask = get_nuclei_cellpose(parameters, img, cell_mask)
    golgi_mask = get_golgi_mask(parameters, img, cell_mask)
    if parameters["channel_expression_marker"] >= 0:
        im_marker = img[:,:,parameters["channel_expression_marker"]]

    single_cell_props = pd.DataFrame()
    counter = 0

    for label in range(1,np.max(cell_mask)-1):
        
        single_cell_mask = np.where(cell_mask ==label, 1, 0)
        single_nucleus_mask = np.where(nuclei_mask==label, 1, 0)
        single_golgi_mask = np.where(golgi_mask ==label, 1, 0)    
        
        area_cell_px2 = len(single_cell_mask[single_cell_mask==1])
        area_golgi_px2 = len(single_golgi_mask[single_golgi_mask==1])
        area_nucleus_px2 = len(single_nucleus_mask[single_nucleus_mask==1])    
        
        if (area_nucleus_px2) < 10 or (area_golgi_px2 < 10):
            continue
                    
        #print(len(single_cell_mask[single_cell_mask==1]))
        #print(len(single_nucleus_mask[single_nucleus_mask==1]))
        #print(len(single_golgi_mask[single_golgi_mask==1]))
          
        regions = skimage.measure.regionprops(single_cell_mask)
        for props in regions:
            x_cell, y_cell = props.centroid
            orientation = props.orientation
            minor_axis_length = props.minor_axis_length
            major_axis_length = props.major_axis_length
            area = props.area
            perimeter = props.perimeter     


        regions = skimage.measure.regionprops(single_nucleus_mask)
        for props in regions:
            x_nucleus, y_nucleus = props.centroid
            orientation_nuc = props.orientation
            minor_axis_length_nuc = props.minor_axis_length
            major_axis_length_nuc = props.major_axis_length
            area_nuc = props.area
            perimeter_nuc = props.perimeter     


        single_cell_props.at[counter, "label"] = label
        single_cell_props.at[counter, "X_cell"] = x_cell
        single_cell_props.at[counter, "Y_cell"] = y_cell
        single_cell_props.at[counter, "shape_orientation"] = orientation
        single_cell_props.at[counter, "flow_shape_alignment"] = np.sin(orientation) # assumes flow from left to right anlong x-axis
        single_cell_props.at[counter, "major_axis_length"] = major_axis_length
        single_cell_props.at[counter, "minor_axis_length"] = minor_axis_length
        single_cell_props.at[counter, "area"] = area
        single_cell_props.at[counter, "perimeter"] = perimeter
        
        single_cell_props.at[counter, "X_nuc"] = x_nucleus
        single_cell_props.at[counter, "Y_nuc"] = y_nucleus
        single_cell_props.at[counter, "shape_orientation_nuc"] = orientation_nuc
        single_cell_props.at[counter, "flow_shape_alignment_nuc"] = np.sin(orientation_nuc) # assumes flow from left to right anlong x-axis
        single_cell_props.at[counter, "major_axis_length"] = major_axis_length_nuc
        single_cell_props.at[counter, "minor_axis_length"] = minor_axis_length_nuc
        single_cell_props.at[counter, "area"] = area_nuc
        single_cell_props.at[counter, "perimeter"] = perimeter_nuc
 


        if parameters["channel_expression_marker"] >= 0:
            regions = skimage.measure.regionprops(single_cell_mask, intensity_image=im_marker)
            for props in regions:
                mean_expression = props.mean_intensity
            single_cell_props.at[counter, "mean_expression"] = mean_expression
            
            regions = skimage.measure.regionprops(single_nucleus_mask, intensity_image=im_marker)
            for props in regions:
                mean_expression_nuc = props.mean_intensity
            single_cell_props.at[counter, "mean_expression_nuc"] = mean_expression_nuc

            cytosol_mask = np.logical_xor(single_cell_mask.astype(bool), single_nucleus_mask.astype(bool))

            regions = skimage.measure.regionprops(cytosol_mask.astype(int), intensity_image=im_marker)
            for props in regions:
                mean_expression_cyt = props.mean_intensity
            single_cell_props.at[counter, "mean_expression_cyt"] = mean_expression_cyt

        
        regions = skimage.measure.regionprops(single_golgi_mask)
        for props in regions:
            x_golgi, y_golgi = props.centroid
        
        single_cell_props.at[counter, "X_golgi"] = x_golgi
        single_cell_props.at[counter, "Y_golgi"] = y_golgi
        
        distance2 = (x_golgi - x_nucleus)**2
        distance2 += (y_golgi - y_nucleus)**2
        
        single_cell_props.at[counter, "distance"] = np.sqrt(distance2)
        
        vec_x = x_golgi - x_nucleus
        vec_y = y_golgi - y_nucleus
        angle_rad_ = np.arctan2(vec_x, vec_y)
        
        angle_rad = angle_rad_
        
        if (angle_rad_ < 0.0):
            angle_rad = 2.0*np.pi + angle_rad_
        
        single_cell_props.at[counter, "vec_X"] = vec_x
        single_cell_props.at[counter, "vec_Y"] = vec_y
        single_cell_props.at[counter, "angle_rad"] = angle_rad
        single_cell_props.at[counter, "flow_alignment"] = np.sin(angle_rad)
        single_cell_props.at[counter, "angle_deg"] = 180.0*angle_rad/np.pi   
        
        counter += 1

    im_junction = img[:,:,int(parameters["channel_junction"])]
    im_marker = img[:,:,int(parameters["channel_expression_marker"])]
    
    plot_polarity(parameters, im_junction, [cell_mask, nuclei_mask, golgi_mask], single_cell_props, filename)
    plot_marker(parameters, im_junction, [cell_mask, nuclei_mask, golgi_mask], single_cell_props, filename)

    return single_cell_props

def plot_polarity(parameters, im_junction, masks, single_cell_props, filename):
    
    output_path = parameters['output_path']
    output_filename = parameters["output_filename"]
    output_filepath = output_path + output_filename

    fig, ax = plt.subplots()
    
    cell_mask = masks[0] 
    nuclei_mask = masks[1].astype(bool)
    golgi_mask = masks[2].astype(bool)

    nuclei_mask_ = np.where(nuclei_mask == True, 70, 0)
    golgi_mask_ = np.where(golgi_mask == True, 1, 0)

    ax.imshow(im_junction, cmap=plt.cm.gray, alpha = 1.0)
    ax.imshow(cell_mask, cmap=plt.cm.Set3, alpha = 0.25)
    ax.imshow(np.ma.masked_where(nuclei_mask_ == 0, nuclei_mask_),  plt.cm.gist_rainbow, vmin=0, vmax=100, alpha = 0.5)
    ax.imshow(np.ma.masked_where(golgi_mask_ == 0, golgi_mask_),  plt.cm.gist_rainbow, vmin=0, vmax=100, alpha = 0.5)

    for index, row in single_cell_props.iterrows():
        ax.plot( row["Y_nuc"], row["X_nuc"], '.g', markersize=1)
        ax.plot( row["Y_golgi"], row["X_golgi"], '.m', markersize=1)
        ax.arrow(row["Y_nuc"], row["X_nuc"], row["Y_golgi"]- row["Y_nuc"],row["X_golgi"]- row["X_nuc"], color = 'white', width = 2)

    ax.set_xlim(0,im_junction.shape[0])
    ax.set_ylim(0,im_junction.shape[1])
    plt.savefig(output_path + filename + "_nuclei_golgi_vector.pdf")
    plt.savefig(output_path + filename + "_nuclei_golgi_vector.png")

    return 0


def plot_marker(parameters, im_marker, masks, single_cell_props, filename):
    
    output_path = parameters['output_path']
    output_filename = parameters["output_filename"]
    output_filepath = output_path + output_filename

    fig, ax = plt.subplots(1,3, figsize=(30,10))
    
    cell_mask = masks[0] 
    nuclei_mask = masks[1].astype(bool)

    nuclei_mask_ = np.where(nuclei_mask == True, 70, 0)

    ax[0].imshow(im_marker, cmap=plt.cm.gray, alpha = 1.0)
    ax[1].imshow(im_marker, cmap=plt.cm.gray, alpha = 1.0)
    ax[1].imshow(cell_mask, cmap=plt.cm.Set3, alpha = 0.25)

    ax[2].imshow(im_marker, cmap=plt.cm.gray, alpha = 1.0)
    ax[2].imshow(nuclei_mask,  cmap = plt.cm.Set3, alpha = 0.5)

    for index, row in single_cell_props.iterrows():
        ax[2].text( row["Y_nuc"], row["X_nuc"], str(np.round(row["mean_expression_nuc"],1)), color = "w", fontsize=6)
        ax[1].text( row["Y_cell"], row["X_cell"], str(np.round(row["mean_expression"],1)), color = "w", fontsize=6)
    #    ax.plot( row["Y_golgi"], row["X_golgi"], '.m', markersize=1)
    #    ax.arrow(row["Y_nuc"], row["X_nuc"], row["Y_golgi"]- row["Y_nuc"],row["X_golgi"]- row["X_nuc"], color = 'white', width = 2)

    #ax.set_xlim(0,im_marker.shape[0])
    #ax.set_ylim(0,im_marker.shape[1])
    plt.savefig(output_path + filename + "_marker_expression.pdf")
    plt.savefig(output_path + filename + "_marker_expression.png")

    return 0

