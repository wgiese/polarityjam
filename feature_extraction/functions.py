import os
import sys
import yaml
import skimage.io
import numpy as np
import cellpose.utils
import cellpose.io
import cellpose.models
import scipy.ndimage as ndi
import skimage.filters 
import skimage.measure
import skimage.segmentation
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import math
import numpy as np
import tifffile as tiff
from skimage.measure import regionprops
from skimage.future.graph import RAG
import networkx as nw
import plot_fcts


def read_parameters(parameter_file):

    parameters = dict()

    if sys.platform.startswith("win"):
        if os.path.exists("..\\feature_extraction\\base\\parameters.yml"):
            with open("../feature_extraction/base/parameters.yml") as file:
                parameters = yaml.load(file, Loader=yaml.FullLoader)
    else:
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

    img_ = skimage.io.imread(filename)

    if len(img_.shape) <= 2:
        img_ = np.array([img_,img_]) 
        

    if img_.shape[0] < min(img_.shape[1],img_.shape[2]):
        print("Warning: channel is on the first dimension of the image.")
        img = img_.reshape(img_.shape[1], img_.shape[2], img_.shape[0])
    else:
        img = img_

    return img

def get_image_for_segmentation(parameters, img):

    print(parameters["channel_junction"])
    print(img.shape)
    im_junction = img[:,:,int(parameters["channel_junction"])]
    if parameters["channel_nucleus"] >= 0:
    	im_nucleus = img[:,:,int(parameters["channel_nucleus"])]
    	im_seg = np.array([im_junction, im_nucleus])
    else:
    	return im_junction

    return im_seg

def get_cellpose_segmentation(parameters, im_seg):

    model = cellpose.models.Cellpose(gpu=parameters['use_gpu'], model_type='cyto')
    if parameters["channel_nucleus"] >= 0:
    	channels = [1,2]
    else:
    	channels = [0,0]

    masks, flows, styles, diams = model.eval(im_seg, diameter=parameters["estimated_cell_diameter"], channels=channels)  

    if parameters["clear_border"]:
        masks = skimage.segmentation.clear_border(masks)

    return masks

def get_nuclei_mask(parameters, img, cellpose_mask):

    img_nuclei_blur = ndi.gaussian_filter(img[:,:,parameters["channel_nucleus"]], sigma=3)

    nuclei_mask = np.where(img_nuclei_blur > skimage.filters.threshold_otsu(img_nuclei_blur), True, False)

    nuclei_label = nuclei_mask*cellpose_mask

    return nuclei_label

def get_nuclei_cellpose(parameters, img, cellpose_mask):

    model = cellpose.models.Cellpose(gpu=parameters["use_gpu"], model_type='nuclei')

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


def get_features_from_cellpose_seg(parameters, img, cell_mask, filename, output_path):

    if parameters["channel_nucleus"] >= 0: 
        nuclei_mask = get_nuclei_mask(parameters, img, cell_mask)
    #nuclei_mask = get_nuclei_cellpose(parameters, img, cell_mask)
    if parameters["channel_golgi"] >= 0:
        golgi_mask = get_golgi_mask(parameters, img, cell_mask)
    
    if parameters["channel_expression_marker"] >= 0:
        im_marker = img[:,:,parameters["channel_expression_marker"]]

    single_cell_props = pd.DataFrame()
    counter = 0

    for label in np.unique(cell_mask):
#initialize graph - no features assosciated with nodes
        graph_nf = orientation_graph_nf(cell_mask)
        if label == 0:
            continue
        #in range(1,np.max(cell_mask)-1):
        
        single_cell_mask = np.where(cell_mask ==label, 1, 0)
        if parameters["channel_nucleus"] >= 0:
            single_nucleus_mask = np.where(nuclei_mask==label, 1, 0)
            if len(single_nucleus_mask[single_nucleus_mask==1]) < parameters["min_nucleus_size"]:
                continue
            if parameters["channel_golgi"] >= 0:
                single_golgi_mask = np.where(golgi_mask ==label, 1, 0)    
                if len(single_golgi_mask[single_golgi_mask==1]) < parameters["min_golgi_size"]:
                    continue

        #area_golgi_px2 = len(single_golgi_mask[single_golgi_mask==1])
        #area_nucleus_px2 = len(single_nucleus_mask[single_nucleus_mask==1])    
        
        #if (area_nucleus_px2) < 10 or (area_golgi_px2 < 10):
        #    continue
                    
        #print(len(single_cell_mask[single_cell_mask==1]))
        #print(len(single_nucleus_mask[single_nucleus_mask==1]))
        #print(len(single_golgi_mask[single_golgi_mask==1]))
        
        if len(single_cell_mask[single_cell_mask==1]) < parameters["min_cell_size"]:
            continue
  
        regions = skimage.measure.regionprops(single_cell_mask)
        
        for props in regions:
            x_cell, y_cell = props.centroid
            # note, the values of orientation from props are in [-pi/2,pi/2] with zero along the y-axis
            orientation = np.pi/2.0 - props.orientation
            minor_axis_length = props.minor_axis_length
            major_axis_length = props.major_axis_length
            eccentricity = props.eccentricity
            area = props.area
            perimeter = props.perimeter     

            neighbours = len(list(graph_nf.neighbors(label)))
        
        single_cell_props.at[counter, "filename"] = filename
        single_cell_props.at[counter, "label"] = label
        single_cell_props.at[counter, "X_cell"] = x_cell
        single_cell_props.at[counter, "Y_cell"] = y_cell
        single_cell_props.at[counter, "shape_orientation"] = orientation
        single_cell_props.at[counter, "flow_shape_alignment"] = np.sin(orientation) # assumes flow from left to right anlong x-axis
        single_cell_props.at[counter, "major_axis_length"] = major_axis_length
        single_cell_props.at[counter, "minor_axis_length"] = minor_axis_length
        single_cell_props.at[counter, "eccentricity"] = eccentricity
        single_cell_props.at[counter, "major_to_minor_ratio"] = major_axis_length/minor_axis_length
        single_cell_props.at[counter, "area"] = area
        single_cell_props.at[counter, "perimeter"] = perimeter
        single_cell_props.at[counter, "n_neighbors"] = neighbours
        

         
        if parameters["channel_nucleus"] >= 0:
            regions = skimage.measure.regionprops(single_nucleus_mask)
            for props in regions:
                x_nucleus, y_nucleus = props.centroid
                # note, the values of orientation from props are in [-pi/2,pi/2] with zero along the y-axis
                orientation_nuc = np.pi/2.0 - props.orientation
                minor_axis_length_nuc = props.minor_axis_length
                major_axis_length_nuc = props.major_axis_length
                area_nuc = props.area
                perimeter_nuc = props.perimeter     
                eccentricity_nuc = props.eccentricity

            
            single_cell_props.at[counter, "X_nuc"] = x_nucleus
            single_cell_props.at[counter, "Y_nuc"] = y_nucleus
            single_cell_props.at[counter, "shape_orientation_nuc"] = orientation_nuc
            single_cell_props.at[counter, "flow_shape_alignment_nuc"] = np.sin(orientation_nuc) # assumes flow from left to right anlong x-axis
            single_cell_props.at[counter, "major_axis_length_nuc"] = major_axis_length_nuc
            single_cell_props.at[counter, "minor_axis_length_nuc"] = minor_axis_length_nuc
            single_cell_props.at[counter, "area"] = area_nuc
            single_cell_props.at[counter, "perimeter"] = perimeter_nuc
            single_cell_props.at[counter, "eccentricity_nuc"] = eccentricity_nuc
            single_cell_props.at[counter, "major_to_minor_ratio_nuc"] = major_axis_length/minor_axis_length

        ### compute marker polarity

        if parameters["channel_expression_marker"] >= 0:
            regions = skimage.measure.regionprops(single_cell_mask, intensity_image=im_marker)
            for props in regions:
                mean_expression = props.mean_intensity
                area_cell = props.area
                x_weighted, y_weighted = props.weighted_centroid
                x_cell, y_cell = props.centroid
            
            single_cell_props.at[counter, "mean_expression"] = mean_expression
            single_cell_props.at[counter, "sum_expression"] = mean_expression*area_cell
            single_cell_props.at[counter, "X_weighted"] = x_weighted
            single_cell_props.at[counter, "Y_weighted"] = y_weighted
             
            distance2 = (x_cell - x_weighted)**2
            distance2 += (y_cell - y_weighted)**2
            
            single_cell_props.at[counter, "maker_vec_norm"] = np.sqrt(distance2)
            
            vec_x = x_weighted - x_cell
            vec_y = y_weighted - y_cell
            angle_rad_ = np.arctan2(vec_x, vec_y)
            
            angle_rad = angle_rad_
            
            if (angle_rad_ < 0.0):
                angle_rad = 2.0*np.pi + angle_rad_
            
            single_cell_props.at[counter, "marker_polarity_rad"] = angle_rad
            single_cell_props.at[counter, "marker_polarity_deg"] = 180.0*angle_rad/np.pi   
        


            
            if parameters["channel_nucleus"] >= 0:
                regions = skimage.measure.regionprops(single_nucleus_mask, intensity_image=im_marker)
                for props in regions:
                    mean_expression_nuc = props.mean_intensity
                    area_nuc = props.area
                single_cell_props.at[counter, "mean_expression_nuc"] = mean_expression_nuc
                single_cell_props.at[counter, "sum_expression_nuc"] = mean_expression_nuc*area_nuc

                cytosol_mask = np.logical_xor(single_cell_mask.astype(bool), single_nucleus_mask.astype(bool))

                regions = skimage.measure.regionprops(cytosol_mask.astype(int), intensity_image=im_marker)
                for props in regions:
                    mean_expression_cyt = props.mean_intensity
                    area_cyt = props.area
                single_cell_props.at[counter, "mean_expression_cyt"] = mean_expression_cyt
                single_cell_props.at[counter, "sum_expression_cyt"] = mean_expression_cyt*area_cyt
            
            membrane_mask = get_outline_from_mask(single_cell_mask.astype(bool), parameters["membrane_thickness"])

            regions = skimage.measure.regionprops(membrane_mask.astype(int), intensity_image=im_marker)
            for props in regions:
                mean_expression_mem = props.mean_intensity
                area_mem = props.area
            single_cell_props.at[counter, "mean_expression_mem"] = mean_expression_mem
            single_cell_props.at[counter, "sum_expression_mem"] = mean_expression_mem*area_mem
 

            

   
        if (parameters["channel_nucleus"] >= 0) and (parameters["channel_golgi"] >= 0):
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

        #rag.nodes["label"][feature_of_interest] = single_cell_props.at[counter, feature_of_interest]
        counter += 1

    im_junction = img[:,:,int(parameters["channel_junction"])]
    im_marker = img[:,:,int(parameters["channel_expression_marker"])]

    if parameters["plot_polarity"] and (parameters["channel_golgi"] >= 0):
        plot_fcts.plot_polarity(parameters, im_junction, [cell_mask, nuclei_mask, golgi_mask], single_cell_props, filename, output_path)
    if parameters["plot_marker"] and (parameters["channel_nucleus"] >= 0):
        plot_fcts.plot_marker(parameters, im_marker, [cell_mask, nuclei_mask], single_cell_props, filename, output_path)
        plot_fcts.plot_marker_polarity(parameters, im_marker, [cell_mask], single_cell_props, filename, output_path)
    if parameters["plot_marker"] and (parameters["channel_nucleus"] < 0):
        plot_fcts.plot_marker(parameters, im_marker, [cell_mask], single_cell_props, filename, output_path)
        plot_fcts.plot_marker_polarity(parameters, im_marker, [cell_mask], single_cell_props, filename, output_path)
    if parameters["plot_alignment"] and (parameters["channel_nucleus"] >= 0):
        plot_fcts.plot_alignment(parameters, im_junction, [cell_mask, nuclei_mask], single_cell_props, filename, output_path)
    if parameters["plot_marker"] and (parameters["channel_nucleus"] < 0):
        plot_fcts.plot_alignment(parameters, im_junction, [cell_mask], single_cell_props, filename, output_path)
    if parameters["plot_ratio_method"]:
        plot_fcts.plot_ratio_method(parameters, im_junction, [cell_mask], single_cell_props, filename, output_path)
######################    
    '''
    weihgts = psy.lib.weights.W.from_networkx(rag)

    moron_eye_feature_list = []

    rag_labels = list(rag.nodes)
    moran_keys = weihgts.neighbors.keys()

    for nombritas in zip(rag_labels,moran_keys):
        #print(nombritas)
        #print(len(list(rag.neighbors(nombritas[0]))))
        #print(len(weihgts.neighbors[nombritas[1]]))

        feature2append = rag.nodes[nombritas[0]]
        single_feature = (feature2append[feature_of_interest])

        moron_eye_feature_list.append(single_feature)

    mi = psy.explore.esda.Moran(moron_eye_feature_list, weihgts, two_tailed=False)
    print("%.3f"%mi.I)
    print(mi.EI)
    print("%f"%mi.p_norm)
    '''    
    return single_cell_props
######################
def get_outline_from_mask(mask, width = 1):
 
    dilated_mask = ndi.morphology.binary_dilation(mask.astype(bool), iterations = width)
    eroded_mask = ndi.morphology.binary_erosion(mask.astype(bool), iterations = width)
    outline_mask = np.logical_xor(dilated_mask, eroded_mask)
    
    return outline_mask
def orientation_graph_nf(img):

    rag = RAG(img.astype("int"))
    return(rag)
def orientation_graph(img):

    rag = RAG(img.astype("int"))
    rag.remove_node(0)

    regions = regionprops(img.astype("int"))
    for region in regions:
        rag.nodes[region['label']]['orientation'] = region['orientation']
        rag.nodes[region['label']]['area'] = region['area']
        rag.nodes[region['label']]['polarity'] = region['major_axis_length'] / region['minor_axis_length']
        rag.nodes[region['label']]['aspect_ratio'] = (region["bbox"][2] - region["bbox"][0]) / (region["bbox"][3] - region["bbox"][1])
    return(rag)
def remove_edges(mask):
    segments = mask.astype("int")
    end_1,end_2 = mask.shape[0],mask.shape[1]

    start_1, start_2 = 0,0

    the_start = np.empty(mask.shape[0])
    the_start.fill(int(start_1))
    the_end = np.empty(mask.shape[0])
    the_end.fill(int(end_1 - 1 ))
    lower_right_arr = np.asarray(range(start_1,end_1))

    #list of points with 0, 0 - max
    roof = np.asarray([the_start,lower_right_arr])
    #list of of points with max 0-max
    floor = np.asarray([the_end,lower_right_arr])
    #list of point 0-max, 0
    left_wall = np.asarray([lower_right_arr,the_start])
    #list of point 0-max, max
    right_wall = np.asarray([lower_right_arr,the_end])

    concat_arr = np.hstack([left_wall,right_wall,roof,floor]).astype("int")
    x_indexed = segments[(concat_arr[1,:]),(concat_arr[0,:])]

    for elemet in np.unique(x_indexed):
        segments[segments == elemet] = 0
    return(segments)
def remove_islands(frame_graph, mask):
    list_of_islands = []
    #Get list of islands - nodes with no neighbkluors,
    #remove nodes with neighbours
    #frame_graph = orientation_graph(dat_no_edges[i,:,:])
    for nodez in frame_graph.nodes:
        if len(list(frame_graph.neighbors(nodez))) == 0:
            list_of_islands.append(nodez)
    
    print(list_of_islands)
    #remove islands from image and graph
    for elemet in np.unique(list_of_islands):
        frame_graph.remove_node(elemet)
        mask[:,:][mask[:,:] == elemet] = 0
    return(frame_graph,mask)
