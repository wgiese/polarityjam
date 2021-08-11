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


def get_golgi_mask(parameters, img, cellpose_mask):

    img_golgi_blur = ndi.gaussian_filter(img[:,:,parameters["channel_golgi"]], sigma=3)
    
    golgi_mask = np.where(img_golgi_blur > skimage.filters.threshold_otsu(img_golgi_blur), True, False)

    golgi_label = golgi_mask*cellpose_mask
    
    return golgi_label


def get_features_from_cellpose_seg(parameters, img, cell_mask):

    nuclei_mask = get_nuclei_mask(parameters, img, cell_mask)
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

        single_cell_props.at[counter, "label"] = label
        single_cell_props.at[counter, "X_cell"] = x_cell
        single_cell_props.at[counter, "Y_cell"] = y_cell
        single_cell_props.at[counter, "shape_orientation"] = orientation
        single_cell_props.at[counter, "major_axis_length"] = major_axis_length
        single_cell_props.at[counter, "minor_axis_length"] = minor_axis_length
        single_cell_props.at[counter, "area"] = area
        single_cell_props.at[counter, "perimeter"] = perimeter
        
        if parameters["channel_expression_marker"] >= 0:
            regions = skimage.measure.regionprops(single_cell_mask, intensity_image=im_marker)
            for props in regions:
                mean_expression = props.mean_intensity
            single_cell_props.at[counter, "mean_expression"] = mean_expression

        regions = skimage.measure.regionprops(single_nucleus_mask)
        for props in regions:
            x_nucleus, y_nucleus = props.centroid
        
        single_cell_props.at[counter, "X_nuc"] = x_nucleus
        single_cell_props.at[counter, "Y_nuc"] = y_nucleus
        
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
        single_cell_props.at[counter, "angle_deg"] = 180.0*angle_rad/np.pi   
        
        counter += 1


    return single_cell_props

