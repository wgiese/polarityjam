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
import math
import numpy as np
import tifffile as tiff
from skimage.measure import regionprops
from skimage.future.graph import RAG
import networkx as nw
import plot_fcts


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

    if img.shape[0] < min(img.shape[1],img.shape[2]):
        print("Warning: channel is on the first dimension of the image.")
        img = img.reshape(np.shape[1], np.shape[2], np.shape[0])

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

    model = cellpose.models.Cellpose(gpu=True, model_type='cyto')
    if parameters["channel_nucleus"] >= 0:
    	channels = [1,2]
    else:
    	channels = [0,0]

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

    if parameters["channel_nucleus"] >= 0: 
        nuclei_mask = get_nuclei_mask(parameters, img, cell_mask)
    #nuclei_mask = get_nuclei_cellpose(parameters, img, cell_mask)
    if parameters["channel_golgi"] >= 0:
        golgi_mask = get_golgi_mask(parameters, img, cell_mask)
    
    if parameters["channel_expression_marker"] >= 0:
        im_marker = img[:,:,parameters["channel_expression_marker"]]

    single_cell_props = pd.DataFrame()
    counter = 0

    for label in range(1,np.max(cell_mask)-1):
        
        single_cell_mask = np.where(cell_mask ==label, 1, 0)
        if parameters["channel_nucleus"] >= 0:
            single_nucleus_mask = np.where(nuclei_mask==label, 1, 0)
        if parameters["channel_golgi"] >= 0:
            single_golgi_mask = np.where(golgi_mask ==label, 1, 0)    
        
        #area_golgi_px2 = len(single_golgi_mask[single_golgi_mask==1])
        #area_nucleus_px2 = len(single_nucleus_mask[single_nucleus_mask==1])    
        
        #if (area_nucleus_px2) < 10 or (area_golgi_px2 < 10):
        #    continue
                    
        #print(len(single_cell_mask[single_cell_mask==1]))
        #print(len(single_nucleus_mask[single_nucleus_mask==1]))
        #print(len(single_golgi_mask[single_golgi_mask==1]))
          
        regions = skimage.measure.regionprops(single_cell_mask)
        for props in regions:
            x_cell, y_cell = props.centroid
            # note, the values of orientation from props are in [-pi/2,pi/2] with zero along the y-axis
            orientation = np.pi/2.0 - props.orientation
            minor_axis_length = props.minor_axis_length
            major_axis_length = props.major_axis_length
            area = props.area
            perimeter = props.perimeter     


        
        single_cell_props.at[counter, "filename"] = filename
        single_cell_props.at[counter, "label"] = label
        single_cell_props.at[counter, "X_cell"] = x_cell
        single_cell_props.at[counter, "Y_cell"] = y_cell
        single_cell_props.at[counter, "shape_orientation"] = orientation
        single_cell_props.at[counter, "flow_shape_alignment"] = np.sin(orientation) # assumes flow from left to right anlong x-axis
        single_cell_props.at[counter, "major_axis_length"] = major_axis_length
        single_cell_props.at[counter, "minor_axis_length"] = minor_axis_length
        single_cell_props.at[counter, "area"] = area
        single_cell_props.at[counter, "perimeter"] = perimeter
        

         
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


            
            single_cell_props.at[counter, "X_nuc"] = x_nucleus
            single_cell_props.at[counter, "Y_nuc"] = y_nucleus
            single_cell_props.at[counter, "shape_orientation_nuc"] = orientation_nuc
            single_cell_props.at[counter, "flow_shape_alignment_nuc"] = np.sin(orientation_nuc) # assumes flow from left to right anlong x-axis
            single_cell_props.at[counter, "major_axis_length_nuc"] = major_axis_length_nuc
            single_cell_props.at[counter, "minor_axis_length_nuc"] = minor_axis_length_nuc
            single_cell_props.at[counter, "area"] = area_nuc
            single_cell_props.at[counter, "perimeter"] = perimeter_nuc
     


        if parameters["channel_expression_marker"] >= 0:
            regions = skimage.measure.regionprops(single_cell_mask, intensity_image=im_marker)
            for props in regions:
                mean_expression = props.mean_intensity
                area_cell = props.area
                #x_weighted, y_weighted = props.centroid_weighted

            single_cell_props.at[counter, "mean_expression"] = mean_expression
            single_cell_props.at[counter, "sum_expression"] = mean_expression*area_cell
            #single_cell_porps.at[counter, "X_weighted"] = x_weighted
            #single_cell_porps.at[counter, "Y_weighted"] = y_weighted
 
            
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
 
            ### compute marker polarity

            

   
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
        
        counter += 1

    im_junction = img[:,:,int(parameters["channel_junction"])]
    im_marker = img[:,:,int(parameters["channel_expression_marker"])]
    
    if parameters["plot_polarity"] and (parameters["channel_golgi"] >= 0):
        plot_fcts.plot_polarity(parameters, im_junction, [cell_mask, nuclei_mask, golgi_mask], single_cell_props, filename)
    if parameters["plot_marker"] and (parameters["channel_nucleus"] >= 0):
        plot_fcts.plot_marker(parameters, im_marker, [cell_mask, nuclei_mask], single_cell_props, filename)
    if parameters["plot_marker"] and (parameters["channel_nucleus"] < 0):
        plot_fcts.plot_marker(parameters, im_marker, [cell_mask], single_cell_props, filename)
    if parameters["plot_alignment"] and (parameters["channel_nucleus"] >= 0):
        plot_fcts.plot_alignment(parameters, im_junction, [cell_mask, nuclei_mask], single_cell_props, filename)
    if parameters["plot_marker"] and (parameters["channel_nucleus"] < 0):
        plot_fcts.plot_alignment(parameters, im_junction, [cell_mask], single_cell_props, filename)

    return single_cell_props

def get_outline_from_mask(mask, width = 1):
 
    dilated_mask = ndi.morphology.binary_dilation(mask.astype(bool), iterations = width)
    eroded_mask = ndi.morphology.binary_erosion(mask.astype(bool), iterations = width)
    outline_mask = np.logical_xor(dilated_mask, eroded_mask)
    
    return outline_mask

def plot_polarity(parameters, im_junction, masks, single_cell_props, filename):
    
    output_path = parameters['output_path']
    output_filename = parameters["output_filename"]
    output_filepath = output_path + output_filename

    fig, ax = plt.subplots(figsize=(10,10))
    
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
        if parameters["show_polarity_angles"]:
            ax.text(row["Y_cell"], row["X_cell"], str(int(np.round(row["angle_deg"],0))), color = "yellow", fontsize = 6)
    

    ax.set_xlim(0,im_junction.shape[0])
    ax.set_ylim(0,im_junction.shape[1])
    plt.savefig(output_path + filename + "_nuclei_golgi_vector.pdf")
    plt.savefig(output_path + filename + "_nuclei_golgi_vector.png")

    return 0


def plot_marker(parameters, im_marker, masks, single_cell_props, filename):
    
    output_path = parameters['output_path']
    output_filename = parameters["output_filename"]
    output_filepath = output_path + output_filename

    sub_figs = len(masks) + 1
    fig, ax = plt.subplots(1, sub_figs, figsize=(10*sub_figs,10))
    
    cell_mask = masks[0]
    if sub_figs > 2: 
        nuclei_mask = masks[1].astype(bool)

        nuclei_mask_ = np.where(nuclei_mask == True, 70, 0)

    for i in range(sub_figs):
        ax[i].imshow(im_marker, cmap=plt.cm.gray, alpha = 1.0)

    #ax[1].imshow(im_marker, cmap=plt.cm.gray, alpha = 1.0)
    #ax[2].imshow(im_marker, cmap=plt.cm.gray, alpha = 1.0)


    #outline_nuc = get_outline_from_mask(nuclei_mask, parameters["outline_width"])
    #outline_nuc_ = np.where(outline_nuc == True, 30, 0)
    #ax[0].imshow(np.ma.masked_where(outline_nuc_ == 0, outline_nuc_),  plt.cm.Wistia, vmin=0, vmax=100, alpha = 0.75)
 
    if sub_figs > 2:
        outline_nuc = get_outline_from_mask(nuclei_mask, parameters["outline_width"])
        outline_nuc_ = np.where(outline_nuc == True, 30, 0)
        ax[sub_figs - 1].imshow(np.ma.masked_where(outline_nuc_ == 0, outline_nuc_),  plt.cm.Wistia, vmin=0, vmax=100, alpha = 0.75)

    outlines_cell = np.zeros((im_marker.shape[0], im_marker.shape[1]))
    outlines_mem = np.zeros((im_marker.shape[0], im_marker.shape[1]))

    for label in range(1,np.max(cell_mask)+1):
        single_cell_mask = np.where(cell_mask ==label, 1, 0)
        outline_cell = get_outline_from_mask(single_cell_mask, parameters["outline_width"])
        outline_cell_ = np.where(outline_cell == True, 30, 0)
        outlines_cell += outline_cell_
        
        outline_mem = get_outline_from_mask(single_cell_mask, parameters["membrane_thickness"])
        outline_mem_ = np.where(outline_mem == True, 30, 0)
        outlines_mem += outline_mem_
    
        #ax[1].imshow(np.ma.masked_where(outline_nuc_ == 0, outline_nuc_),  plt.cm.Wistia, vmin=0, vmax=100, alpha = 0.75)
        #ax[1].imshow(np.ma.masked_where(outline_cell_ == 0, outline_cell_),  plt.cm.Wistia, vmin=0, vmax=100, alpha = 0.75)
    
    outlines_cell_ = np.where(outlines_cell > 0, 30, 0)
    ax[0].imshow(np.ma.masked_where(outlines_cell_ == 0, outlines_cell_),  plt.cm.Wistia, vmin=0, vmax=100, alpha = 0.5)
        
    outlines_mem_ = np.where(outlines_mem > 0, 30, 0)
    ax[1].imshow(np.ma.masked_where(outlines_mem_ == 0, outlines_mem_),  plt.cm.Wistia, vmin=0, vmax=100, alpha = 0.5)
    
    '''
    for label in range(1,np.max(cell_mask)+1):
        single_cell_mask = np.where(cell_mask ==label, 1, 0)
        outline_cell = get_outline_from_mask(single_cell_mask, parameters["outline_width"])
        outline_cell_ = np.where(outline_cell == True, 30, 0)
        #ax[1].imshow(np.ma.masked_where(outline_nuc_ == 0, outline_nuc_),  plt.cm.Wistia, vmin=0, vmax=100, alpha = 0.75)
        #ax[1].imshow(np.ma.masked_where(outline_cell_ == 0, outline_cell_),  plt.cm.Wistia, vmin=0, vmax=100, alpha = 0.75)
        
        outline_mem = get_outline_from_mask(single_cell_mask, parameters["membrane_thickness"])
        outline_mem_ = np.where(outline_mem == True, 30, 0)
        #ax[2].imshow(np.ma.masked_where(outline_mem_ == 0, outline_mem_),  plt.cm.Wistia, vmin=0, vmax=100, alpha = 0.75)
    '''


    for index, row in single_cell_props.iterrows():
        
        ax[0].text( row["Y_cell"], row["X_cell"], str(np.round(row["mean_expression"],1)), color = "w", fontsize=6)
        ax[1].text( row["Y_cell"], row["X_cell"], str(np.round(row["mean_expression_mem"],1)), color = "w", fontsize=6)
        if sub_figs > 2:
            ax[sub_figs-1].text( row["Y_nuc"], row["X_nuc"], str(np.round(row["mean_expression_nuc"],1)), color = "w", fontsize=6)

    #    ax.plot( row["Y_golgi"], row["X_golgi"], '.m', markersize=1)
    #    ax.arrow(row["Y_nuc"], row["X_nuc"], row["Y_golgi"]- row["Y_nuc"],row["X_golgi"]- row["X_nuc"], color = 'white', width = 2)
    

    ax[0].set_title("mean intensity cell")
    ax[1].set_title("mean intensity membrane")
    if sub_figs > 2:
       ax[sub_figs - 1].set_title("mean intensity nucleus")

    #ax.set_xlim(0,im_marker.shape[0])
    #ax.set_ylim(0,im_marker.shape[1])
    plt.savefig(output_path + filename + "_marker_expression.pdf")
    plt.savefig(output_path + filename + "_marker_expression.png")
    
    #print("plotted expression")
        
    return 0

def plot_alignment(parameters, im_junction, masks, single_cell_props, filename):

    output_path = parameters['output_path']
    output_filename = parameters["output_filename"]
    output_filepath = output_path + output_filename

    sub_figs = len(masks)
    fig, ax = plt.subplots(1, sub_figs)

    cell_mask = masks[0]
    if sub_figs > 1:
        nuclei_mask = masks[1].astype(bool)

    #regions = regionprops(mask)
    #for props in regions:
    #    y0, x0 = props.centroid
    #    orientation = props.orientation
    #    x1 = x0 + math.cos(orientation) * 0.5 * props.minor_axis_length
    #    y1 = y0 - math.sin(orientation) * 0.5 * props.minor_axis_length
    #    x2 = x0 - math.sin(orientation) * 0.5 * props.major_axis_length
    #    y2 = y0 - math.cos(orientation) * 0.5 * props.major_axis_length

    #ax.plot((x0, x1), (y0, y1), '--r', linewidth=0.5)
    #ax.plot((x0, x2), (y0, y2), '--r', linewidth=0.5)
    #ax.plot(x0, y0, '.b', markersize=5)

    #ax[0].imshow(im_marker, cmap=plt.cm.gray, alpha = 1.0)
    
    if sub_figs > 1:
        ax[0].imshow(im_junction, cmap=plt.cm.gray, alpha = 1.0)
        ax[0].imshow(cell_mask, cmap=plt.cm.Set3, alpha = 0.25)

        ax[1].imshow(im_junction, cmap=plt.cm.gray, alpha = 1.0)
        #ax[1].imshow(nuclei_mask,  cmap = plt.cm.Set3, alpha = 0.5)

        nuclei_mask_ = np.where(nuclei_mask == True, 70, 0)
        ax[1].imshow(np.ma.masked_where(nuclei_mask_ == 0, nuclei_mask_),  plt.cm.gist_rainbow, vmin=0, vmax=100, alpha = 0.5)
    else:
        ax.imshow(im_junction, cmap=plt.cm.gray, alpha = 1.0)
        ax.imshow(cell_mask, cmap=plt.cm.Set3, alpha = 0.25)




    for index, row in single_cell_props.iterrows():
        orientation = row['shape_orientation']
        x0 = row['X_cell']
        y0 = row['Y_cell']
        
        #x1 = x0 + math.cos(orientation) * 0.5 * row['major_axis_length']
        #y1 = y0 + math.sin(orientation) * 0.5 * row['major_axis_length']
        #x2 = x0 + math.sin(orientation) * 0.5 * row['minor_axis_length']
        #y2 = y0 - math.cos(orientation) * 0.5 * row['minor_axis_length']
        
        x1 = x0 + math.sin(orientation) * 0.5 * row['major_axis_length']
        y1 = y0 + math.cos(orientation) * 0.5 * row['major_axis_length']
        x2 = x0 + math.cos(orientation) * 0.5 * row['minor_axis_length']
        y2 = y0 - math.sin(orientation) * 0.5 * row['minor_axis_length']
        
        if sub_figs > 1:
            ax[0].plot((y0, y1), (x0, x1), '--r', linewidth=0.5)
            ax[0].plot((y0, y2), (x0, x2), '--r', linewidth=0.5)
            ax[0].plot(y0, x0, '.b', markersize=5)
            orientation_degree = 180.0*orientation/np.pi
            ax[0].text( y0, x0, str(int(np.round(orientation_degree,0))), color = "yellow", fontsize=4)
        else:
            ax.plot((y0, y1), (x0, x1), '--r', linewidth=0.5)
            ax.plot((y0, y2), (x0, x2), '--r', linewidth=0.5)
            ax.plot(y0, x0, '.b', markersize=5)
            orientation_degree = 180.0*orientation/np.pi
            ax.text( y0, x0, str(int(np.round(orientation_degree,0))), color = "yellow", fontsize=4)
               

        
        if sub_figs > 1:
            orientation = row['shape_orientation_nuc']
            x0 = row['X_nuc']
            y0 = row['Y_nuc']
            #x1 = x0 + math.cos(orientation) * 0.5 * row['major_axis_length_nuc']
            #y1 = y0 + math.sin(orientation) * 0.5 * row['major_axis_length_nuc']
            #x2 = x0 + math.sin(orientation) * 0.5 * row['minor_axis_length_nuc']
            #y2 = y0 - math.cos(orientation) * 0.5 * row['minor_axis_length_nuc']

            x1 = x0 + math.sin(orientation) * 0.5 * row['major_axis_length_nuc']
            y1 = y0 + math.cos(orientation) * 0.5 * row['major_axis_length_nuc']
            x2 = x0 + math.cos(orientation) * 0.5 * row['minor_axis_length_nuc']
            y2 = y0 - math.sin(orientation) * 0.5 * row['minor_axis_length_nuc']
            
            ax[1].plot((y0, y1), (x0, x1), '--r', linewidth=0.5)
            ax[1].plot((y0, y2), (x0, x2), '--r', linewidth=0.5)
            ax[1].plot(y0, x0, '.b', markersize=5)

            orientation_degree = 180.0*orientation/np.pi
            ax[1].text( y0, x0, str(int(np.round(orientation_degree,0))), color = "yellow", fontsize=4)

        #ax[1].text( row["Y_cell"], row["X_cell"], str(np.round(row["mean_expression"],1)), color = "w", fontsize=6)

    if sub_figs > 1:
        ax[0].set_title("alignment cell shape")
        ax[0].set_xlim(0,im_junction.shape[0])
        ax[0].set_ylim(0,im_junction.shape[1])
    
        ax[1].set_title("alignment nuclei shape")
        ax[1].set_xlim(0,im_junction.shape[0])
        ax[1].set_ylim(0,im_junction.shape[1])
    else:
        ax.set_title("alignment cell shape")
        ax.set_xlim(0,im_junction.shape[0])
        ax.set_ylim(0,im_junction.shape[1])
    
    

    plt.tight_layout()
    plt.savefig(output_path + filename + "_alignment.pdf")
    plt.savefig(output_path + filename + "_alignment.png")

    return 0


