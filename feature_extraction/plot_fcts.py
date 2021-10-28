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

    fig, ax = plt.subplots(1,3, figsize=(30,10))
    
    cell_mask = masks[0] 
    nuclei_mask = masks[1].astype(bool)

    nuclei_mask_ = np.where(nuclei_mask == True, 70, 0)

    #ax[0].imshow(im_marker, cmap=plt.cm.gray, alpha = 1.0)
    ax[0].imshow(im_marker, cmap=plt.cm.gray, alpha = 1.0)
    #ax[1].imshow(cell_mask, cmap=plt.cm.Set3, alpha = 0.25)

    ax[1].imshow(im_marker, cmap=plt.cm.gray, alpha = 1.0)
    ax[2].imshow(im_marker, cmap=plt.cm.gray, alpha = 1.0)
    #ax[2].imshow(nuclei_mask,  cmap = plt.cm.Set3, alpha = 0.5)


    #outline_nuc = get_outline_from_mask(nuclei_mask, parameters["outline_width"])
    #outline_nuc_ = np.where(outline_nuc == True, 30, 0)
    #ax[0].imshow(np.ma.masked_where(outline_nuc_ == 0, outline_nuc_),  plt.cm.Wistia, vmin=0, vmax=100, alpha = 0.75)
 
    outline_nuc = get_outline_from_mask(nuclei_mask, parameters["outline_width"])
    outline_nuc_ = np.where(outline_nuc == True, 30, 0)
    ax[0].imshow(np.ma.masked_where(outline_nuc_ == 0, outline_nuc_),  plt.cm.Wistia, vmin=0, vmax=100, alpha = 0.75)

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
    ax[1].imshow(np.ma.masked_where(outlines_cell_ == 0, outlines_cell_),  plt.cm.Wistia, vmin=0, vmax=100, alpha = 0.5)
        
    outlines_mem_ = np.where(outlines_mem > 0, 30, 0)
    ax[2].imshow(np.ma.masked_where(outlines_mem_ == 0, outlines_mem_),  plt.cm.Wistia, vmin=0, vmax=100, alpha = 0.5)
    
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
        ax[0].text( row["Y_nuc"], row["X_nuc"], str(np.round(row["mean_expression_nuc"],1)), color = "w", fontsize=6)
        ax[1].text( row["Y_cell"], row["X_cell"], str(np.round(row["mean_expression_cyt"],1)), color = "w", fontsize=6)
        ax[2].text( row["Y_cell"], row["X_cell"], str(np.round(row["mean_expression_mem"],1)), color = "w", fontsize=6)
    #    ax.plot( row["Y_golgi"], row["X_golgi"], '.m', markersize=1)
    #    ax.arrow(row["Y_nuc"], row["X_nuc"], row["Y_golgi"]- row["Y_nuc"],row["X_golgi"]- row["X_nuc"], color = 'white', width = 2)
    

    ax[0].set_title("mean intensity nucleus")
    ax[1].set_title("mean intensity cytosol")
    ax[2].set_title("mean intensity membrane")

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

    fig, ax = plt.subplots(1,2)

    cell_mask = masks[0]
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
    ax[0].imshow(im_junction, cmap=plt.cm.gray, alpha = 1.0)
    ax[0].imshow(cell_mask, cmap=plt.cm.Set3, alpha = 0.25)

    ax[1].imshow(im_junction, cmap=plt.cm.gray, alpha = 1.0)
    #ax[1].imshow(nuclei_mask,  cmap = plt.cm.Set3, alpha = 0.5)

    nuclei_mask_ = np.where(nuclei_mask == True, 70, 0)
    ax[1].imshow(np.ma.masked_where(nuclei_mask_ == 0, nuclei_mask_),  plt.cm.gist_rainbow, vmin=0, vmax=100, alpha = 0.5)



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
        

        ax[0].plot((y0, y1), (x0, x1), '--r', linewidth=0.5)
        ax[0].plot((y0, y2), (x0, x2), '--r', linewidth=0.5)
        ax[0].plot(y0, x0, '.b', markersize=5)

        orientation_degree = 180.0*orientation/np.pi
        ax[0].text( y0, x0, str(int(np.round(orientation_degree,0))), color = "yellow", fontsize=4)
        
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

    ax[0].set_title("alignment cell shape")
    ax[1].set_title("alignment nuclei shape")
    plt.tight_layout()
    ax[0].set_xlim(0,im_junction.shape[0])
    ax[0].set_ylim(0,im_junction.shape[1])
    ax[1].set_xlim(0,im_junction.shape[0])
    ax[1].set_ylim(0,im_junction.shape[1])
    plt.savefig(output_path + filename + "_alignment.pdf")
    plt.savefig(output_path + filename + "_alignment.png")

    return 0


