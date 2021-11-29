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

def plot_polarity(parameters, im_junction, masks, single_cell_props, filename, output_path):
    
    #output_path = parameters['output_path']
    #output_filename = parameters["output_filename"]
    #output_filepath = output_path + output_filename

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
    plt.savefig(output_path + filename + "_nuclei_golgi_vector.svg")
    plt.savefig(output_path + filename + "_nuclei_golgi_vector.png")

    return 0


def plot_marker(parameters, im_marker, masks, single_cell_props, filename, output_path):
    
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

    #for label in range(1,np.max(cell_mask)+1):
    for label in np.unique(cell_mask):
        
        if label == 0:
            continue

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

def plot_marker_polarity(parameters, im_marker, masks, single_cell_props, filename, output_path):

    #sub_figs = len(masks) + 1
    #fig, ax = plt.subplots(1, sub_figs, figsize=(10*sub_figs,10))
    fig, ax = plt.subplots(1, figsize=(10,10))

    
    cell_mask = masks[0]
    ax.imshow(im_marker, cmap=plt.cm.gray, alpha = 1.0)
    
    outlines_cell = np.zeros((im_marker.shape[0], im_marker.shape[1]))

    #for label in range(1,np.max(cell_mask)+1):
    for label in np.unique(cell_mask):
        
        if label == 0:
            continue

        single_cell_mask = np.where(cell_mask ==label, 1, 0)
        outline_cell = get_outline_from_mask(single_cell_mask, parameters["outline_width"])
        outline_cell_ = np.where(outline_cell == True, 30, 0)
        outlines_cell += outline_cell_
        
    outlines_cell_ = np.where(outlines_cell > 0, 30, 0)
    ax.imshow(np.ma.masked_where(outlines_cell_ == 0, outlines_cell_),  plt.cm.Wistia, vmin=0, vmax=100, alpha = 0.5)
     

    for index, row in single_cell_props.iterrows():
        ax.plot( row["Y_cell"], row["X_cell"], '.m', markersize=1)
        ax.plot( row["Y_weighted"], row["X_weighted"], '.m', markersize=1)
        ax.arrow(row["Y_cell"], row["X_cell"], row["Y_weighted"]- row["Y_cell"],row["X_weighted"]- row["X_cell"], color = 'white', width = 2)
        #if parameters["show_polarity_angles"]:
        #    ax.text(row["Y_cell"], row["X_cell"], str(int(np.round(row["angle_deg"],0))), color = "yellow", fontsize = 6)

    ax.set_title("marker polarity")

    #ax.set_xlim(0,im_marker.shape[0])
    #ax.set_ylim(0,im_marker.shape[1])
    plt.savefig(output_path + filename + "_marker_polarity.pdf")
    plt.savefig(output_path + filename + "_marker_polarity.png")
 
    return 0


def plot_alignment(parameters, im_junction, masks, single_cell_props, filename, output_path):

    sub_figs = len(masks)
    fig, ax = plt.subplots(1, sub_figs, figsize=(sub_figs*5,5))

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
    
    cell_eccentricity = np.zeros((cell_mask.shape[0], cell_mask.shape[1]))

    for index, row in single_cell_props.iterrows():
        label = int(row['label'])
        if label == 0:
            continue
        single_cell_mask = np.where(cell_mask ==label, 1, 0)*row['eccentricity']
        cell_eccentricity += single_cell_mask

    if sub_figs > 1:
        ax[0].imshow(im_junction, cmap=plt.cm.gray, alpha = 1.0)
        #ax[0].imshow(cell_mask, cmap=plt.cm.Set3, alpha = 0.25)
        #ax[0].imshow(cell_eccentricity, cmap=plt.cm.bwr, alpha = 0.25)
        
        cax_0 = ax[0].imshow(np.ma.masked_where(cell_mask == 0, cell_eccentricity), cmap=plt.cm.bwr,  vmin=0, vmax=1.0, alpha = 0.5)
        cbar = fig.colorbar(cax_0, ax=ax[0], shrink = 0.3)#, extend='both')
        cbar.set_label('eccentricity')
        ax[1].imshow(im_junction, cmap=plt.cm.gray, alpha = 1.0)
        #ax[1].imshow(nuclei_mask,  cmap = plt.cm.Set3, alpha = 0.5)

        
        nuclei_eccentricity = np.zeros((cell_mask.shape[0], cell_mask.shape[1]))
        for index, row in single_cell_props.iterrows():
            label = int(row['label'])
            if label == 0:
                continue
            single_cell_mask = np.where(cell_mask ==label, True, 0)*row['eccentricity']
            single_nuclei_mask_ = np.logical_and(single_cell_mask, nuclei_mask)
            #np.where(nuclei_mask == label, 1, 0)*row['eccentricity_nuc']
            single_nuclei_mask = np.where(single_nuclei_mask_ == True, 1, 0)*row['eccentricity_nuc']
            nuclei_eccentricity += single_nuclei_mask

        print(nuclei_eccentricity)
        print(np.max(nuclei_eccentricity))
        print(np.min(nuclei_eccentricity))

        nuclei_mask_ = np.where(nuclei_mask == True, 70, 0)
        ##nuclei_mask_ = np.where(nuclei_mask == True, 70, 0)
        #ax[1].imshow(np.ma.masked_where(nuclei_mask_ == 0, nuclei_mask_),  plt.cm.gist_rainbow, vmin=0, vmax=1.0, alpha = 0.5)
        cax_1 = ax[1].imshow(np.ma.masked_where(nuclei_mask_ == 0, nuclei_eccentricity),  plt.cm.bwr, vmin=0, vmax=1.0, alpha = 0.5)
        #fig.colorbar(ax=ax[0])
        cbar = fig.colorbar(cax_1, ax=ax[1], shrink = 0.3)#, extend='both')
        cbar.set_label('eccentricity')

    else:
        ax.imshow(im_junction, cmap=plt.cm.gray, alpha = 1.0)
        #ax.imshow(cell_mask, cmap=plt.cm.Set3, alpha = 0.25)

        #ax.imshow(cell_eccentricity, cmap=plt.cm.Wistia, alpha = 0.25)

        ax.imshow(np.ma.masked_where(cell_mask == 0, cell_eccentricity), cmap=plt.cm.bwr, alpha = 0.25)
        #fig.colorbar(ax=ax)

    for index, row in single_cell_props.iterrows():
        orientation = row['shape_orientation']
        x0 = row['X_cell']
        y0 = row['Y_cell']
        
        #x1 = x0 + math.cos(orientation) * 0.5 * row['major_axis_length']
        #y1 = y0 + math.sin(orientation) * 0.5 * row['major_axis_length']
        #x2 = x0 + math.sin(orientation) * 0.5 * row['minor_axis_length']
        #y2 = y0 - math.cos(orientation) * 0.5 * row['minor_axis_length']
        
        x1_major = x0 + math.sin(orientation) * 0.5 * row['major_axis_length']
        y1_major = y0 + math.cos(orientation) * 0.5 * row['major_axis_length']
        x2_major = x0 - math.sin(orientation) * 0.5 * row['major_axis_length']
        y2_major = y0 - math.cos(orientation) * 0.5 * row['major_axis_length']

        x1_minor = x0 + math.cos(orientation) * 0.5 * row['minor_axis_length']
        y1_minor = y0 - math.sin(orientation) * 0.5 * row['minor_axis_length']
        x2_minor = x0 - math.cos(orientation) * 0.5 * row['minor_axis_length']
        y2_minor = y0 + math.sin(orientation) * 0.5 * row['minor_axis_length']

        if sub_figs > 1:
            ax[0].plot((y1_major, y2_major), (x1_major, x2_major), '--w', linewidth=0.5)
            ax[0].plot((y1_minor, y2_minor), (x1_minor, x2_minor), '--w', linewidth=0.5)
            ax[0].plot(y0, x0, '.b', markersize=5)
            orientation_degree = 180.0*orientation/np.pi
            ax[0].text( y0, x0, str(int(np.round(orientation_degree,0))), color = "yellow", fontsize=4)
        else:
            ax.plot((y1_major, y2_major), (x1_major, x2_major), '--w', linewidth=0.5)
            ax.plot((y1_minor, y2_minor), (x1_minor, x2_minor), '--w', linewidth=0.5)
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

            #x1 = x0 + math.sin(orientation) * 0.5 * row['major_axis_length_nuc']
            #y1 = y0 + math.cos(orientation) * 0.5 * row['major_axis_length_nuc']
            #x2 = x0 + math.cos(orientation) * 0.5 * row['minor_axis_length_nuc']
            #y2 = y0 - math.sin(orientation) * 0.5 * row['minor_axis_length_nuc']
         
            x1_major = x0 + math.sin(orientation) * 0.5 * row['major_axis_length']
            y1_major = y0 + math.cos(orientation) * 0.5 * row['major_axis_length']
            x2_major = x0 - math.sin(orientation) * 0.5 * row['major_axis_length']
            y2_major = y0 - math.cos(orientation) * 0.5 * row['major_axis_length']

            x1_minor = x0 + math.cos(orientation) * 0.5 * row['minor_axis_length']
            y1_minor = y0 - math.sin(orientation) * 0.5 * row['minor_axis_length']
            x2_minor = x0 - math.cos(orientation) * 0.5 * row['minor_axis_length']
            y2_minor = y0 + math.sin(orientation) * 0.5 * row['minor_axis_length']

            ax[1].plot((y1_major, y2_major), (x1_major, x2_major), '--w', linewidth=0.5)
            ax[1].plot((y1_minor, y2_minor), (x1_minor, x2_minor), '--w', linewidth=0.5)

            #ax[1].plot((y0, y1), (x0, x1), '--w', linewidth=0.5)
            #ax[1].plot((y0, y2), (x0, x2), '--w', linewidth=0.5)
            #ax[1].plot(y0, x0, '.b', markersize=5)

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
    plt.savefig(output_path + filename + "_alignment.svg")
    plt.savefig(output_path + filename + "_alignment.png")

    return 0

def plot_ratio_method(parameters, im_junction, masks, single_cell_props, filename, output_path):

    sub_figs = len(masks)
    fig, ax = plt.subplots(1, sub_figs)

    cell_mask = masks[0]

    ax.imshow(im_junction, cmap=plt.cm.gray, alpha = 1.0)
    ax.imshow(cell_mask, cmap=plt.cm.Set3, alpha = 0.25)

    outlines_cell = np.zeros((im_junction.shape[0], im_junction.shape[1]))

    #for label in range(1,np.max(cell_mask)+1):

    for label in np.unique(cell_mask):
        
        if label == 0:
            continue

        single_cell_mask = np.where(cell_mask ==label, 1, 0)
        outline_cell = get_outline_from_mask(single_cell_mask, parameters["membrane_thickness"])
        outline_cell_ = np.where(outline_cell == True, 30, 0)
        outlines_cell += outline_cell_

    outlines_cell_ = np.where(outlines_cell > 0, 30, 0)
    ax.imshow(np.ma.masked_where(outlines_cell_ == 0, outlines_cell_),  plt.cm.bwr, vmin=0, vmax=100, alpha = 0.5)


    for index, row in single_cell_props.iterrows():
        orientation = row['shape_orientation']
        x0 = row['X_cell']
        y0 = row['Y_cell']
        
        
        x1 = x0 + math.sin(np.pi/4.0) * 0.5 * row['major_axis_length']
        y1 = y0 + math.cos(np.pi/4.0) * 0.5 * row['major_axis_length']
        x2 = x0 + math.cos(np.pi/4.0) * 0.5 * row['major_axis_length']
        y2 = y0 - math.sin(np.pi/4.0) * 0.5 * row['major_axis_length']
        
        ax.plot((y0, y1), (x0, x1), '--r', linewidth=0.5)
        ax.plot((y0, y2), (x0, x2), '--r', linewidth=0.5)
        ax.plot(y0, x0, '.b', markersize=5)
        
        x1 = x0 - math.sin(np.pi/4.0) * 0.5 * row['major_axis_length']
        y1 = y0 - math.cos(np.pi/4.0) * 0.5 * row['major_axis_length']
        x2 = x0 - math.cos(np.pi/4.0) * 0.5 * row['major_axis_length']
        y2 = y0 + math.sin(np.pi/4.0) * 0.5 * row['major_axis_length']
        
        ax.plot((y0, y1), (x0, x1), '--r', linewidth=0.5)
        ax.plot((y0, y2), (x0, x2), '--r', linewidth=0.5)
        ax.plot(y0, x0, '.b', markersize=5)
               
        ax.set_title("ratio method")
        ax.set_xlim(0,im_junction.shape[0])
        ax.set_ylim(0,im_junction.shape[1])
    
    plt.tight_layout()
    plt.savefig(output_path + filename + "_ratio_method.pdf")
    plt.savefig(output_path + filename + "_ratio_method.png")

    return 0


