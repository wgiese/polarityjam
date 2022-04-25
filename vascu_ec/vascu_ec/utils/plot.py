import math
from pathlib import Path

import cmocean as cm
import matplotlib as mpl
import tifffile
import numpy as np
from matplotlib import pyplot as plt
from skimage.future import graph
from skimage.measure import label, regionprops

from vascu_ec.utils.rag import orientation_graph_nf
from vascu_ec.utils.seg import get_outline_from_mask
from vascu_ec.vascu_ec_logging import get_logger

# for figure plot resolution
FIGURE_DPI = 300
FONTSIZE_TEXT_ANNOTATIONS = 3
MARKERSIZE = 2
ALPHA_MASKS = 0.5

def set_figure_dpi():
    mpl.rcParams['figure.dpi'] = FIGURE_DPI


def plot_seg_channels(seg_img, output_path, filename):
    """"""
    output_path = Path(output_path)

    if len(seg_img.shape) > 2:
        fig, ax = plt.subplots(1, 2)
        ax[0].imshow(seg_img[0, :, :])
        ax[0].set_title("junction channel")
        ax[1].imshow(seg_img[1, :, :])
        ax[1].set_title("nuclei channel")
        plt.savefig(str(output_path.joinpath(filename + "-seg.png")))
        plt.close(fig)
    else:
        fig, ax = plt.subplots()
        ax.imshow(seg_img[:, :])
        plt.savefig(str(output_path.joinpath(filename + "-seg.png")))
        plt.close(fig)


def plot_cellpose_masks(seg_img, cellpose_mask, output_path, filename):
    """"""
    output_path = Path(output_path)

    if len(seg_img.shape) > 2:
        fig, ax = plt.subplots(1, 3)
        ax[0].imshow(seg_img[0, :, :])
        ax[0].set_title("junction channel")
        ax[1].imshow(seg_img[1, :, :])
        ax[1].set_title("nuclei channel")
        ax[2].imshow(seg_img[0, :, :])
        ax[2].imshow(cellpose_mask, cmap=plt.cm.Set3, alpha=0.5)
        ax[2].set_title("cellpose segmentation")
        plt.savefig(str(output_path.joinpath(filename + "-cellpose-seg.png")))
        plt.close(fig)
    else:
        fig, ax = plt.subplots(1, 2)
        ax[0].imshow(seg_img[:, :])
        ax[1].imshow(cellpose_mask, cmap=plt.cm.Set3, alpha=0.5)
        plt.savefig(str(output_path.joinpath(filename + "-cellpose-seg.png")))
        plt.close(fig)


def plot_organelle_polarity(parameters, im_junction, cell_mask, nuclei_mask, golgi_mask, single_cell_props, base_filename,
                  output_path):
    """ function to plot nuclei-golgi polarity vectors 
    
    parameters  :   dict
                    user defined parameters
    im_junction :   numpy.array (2-dim), float
                    channel containing the junction staining (used for segmentation)
    cell_mask   :   numpy.array (2-dim), int
                    same dimension as im_junction, contains cell masks
    nuclei_mask :   numpy.array (2-dim), int
                    same dimension as im_junction, contains nuclei masks
    golgi_mask  :   numpy.array (2-dim), int
                    same dimension as im_junction, contains golgi masks
    single_cell_props : pandas data frame
    base_filename : string 
                    base_filename for plots
    output_path :   string
                    output path for plots 
    """

    width = parameters["graphics_width"]
    fig, ax = plt.subplots(figsize=(width, width))

    nuclei_mask = nuclei_mask.astype(bool)
    golgi_mask = golgi_mask.astype(bool)

    cell_angle = np.zeros((cell_mask.shape[0], cell_mask.shape[1]))

    for index, row in single_cell_props.iterrows():
        row_label = int(row['label'])
        if row_label == 0:
            continue
        single_cell_mask = np.where(cell_mask == row_label, 1, 0) * row['angle_deg']
        cell_angle += single_cell_mask

    ax.imshow(im_junction, cmap=plt.cm.gray, alpha=1.0)
    cax = ax.imshow(
        np.ma.masked_where(cell_mask == 0, cell_angle), cmap=cm.cm.phase, vmin=0, vmax=360, alpha=0.5
    )
    color_bar = fig.colorbar(cax, ax=ax, shrink=0.3)  # , extend='both')
    color_bar.set_label("polarity angle")
    color_bar.ax.set_yticks([0,90,180,270,360])    

    zero = np.zeros((im_junction.shape[0], im_junction.shape[1]))
    rgb_golgi = np.dstack((golgi_mask.astype(int)*256,zero,zero,golgi_mask.astype(float)*0.5))
    rgb_nuclei = np.dstack((zero,zero,nuclei_mask.astype(int)*256,nuclei_mask.astype(float)*0.5))

    ax.imshow(rgb_nuclei)
    ax.imshow(rgb_golgi)

    for index, row in single_cell_props.iterrows():
        ax.plot(row["Y_nuc"], row["X_nuc"], '.g', markersize=MARKERSIZE)
        ax.plot(row["Y_golgi"], row["X_golgi"], '.m', markersize=MARKERSIZE)
        ax.arrow(row["Y_nuc"], row["X_nuc"], row["Y_golgi"] - row["Y_nuc"], row["X_golgi"] - row["X_nuc"],
                 color='white', width=2)
        if parameters["show_polarity_angles"]:
            ax.text(row["Y_cell"], row["X_cell"], str(int(np.round(row["angle_deg"], 0))), color="yellow", fontsize=6)

    ax.set_xlim(0, im_junction.shape[1])
    ax.set_ylim(0, im_junction.shape[0])
    ax.invert_yaxis()
    ax.axis('off')

    if "pdf" in parameters["graphics_output_format"]:
        plt.savefig(str(Path(output_path).joinpath(base_filename + "_nuclei_golgi_vector.pdf")))
    if "svg" in parameters["graphics_output_format"]:
        plt.savefig(str(Path(output_path).joinpath(base_filename + "_nuclei_golgi_vector.svg")))
    if "png" in parameters["graphics_output_format"]:
        plt.savefig(str(Path(output_path).joinpath(base_filename + "_nuclei_golgi_vector.png")))
    plt.close(fig)

    if "tif" in parameters["graphics_output_format"]:
        tifffile.imsave(str(Path(output_path).joinpath(base_filename + "_nuclei_golgi_polarity.tif")), np.ma.masked_where(cell_mask == 0, cell_angle))


def _get_outline_and_membrane_thickness(im_marker, cell_mask, parameters):
    outlines_cell = np.zeros((im_marker.shape[0], im_marker.shape[1]))
    outlines_mem = np.zeros((im_marker.shape[0], im_marker.shape[1]))

    for cell_label in np.unique(cell_mask):
        # exclude background
        if cell_label == 0:
            continue

        single_cell_mask = np.where(cell_mask == cell_label, 1, 0)
        outline_cell = get_outline_from_mask(single_cell_mask, parameters["outline_width"])
        outline_cell_ = np.where(outline_cell == True, 30, 0)  # todo: magic number 30?
        outlines_cell += outline_cell_

        outline_mem = get_outline_from_mask(single_cell_mask, parameters["membrane_thickness"])
        outline_mem_ = np.where(outline_mem == True, 30, 0)  # todo: magic number 30?
        outlines_mem += outline_mem_

    return [outlines_cell, outlines_mem]


def plot_marker_expression(parameters, im_marker, cell_mask, single_cell_dataset, filename, output_path,
                           nuclei_mask=None):
    number_sub_figs = 2  # mean intensity cell, mean intensity membrane
    if nuclei_mask is not None:
        nuclei_mask = nuclei_mask.astype(bool)
        number_sub_figs = 3  # (optional) mean intensity nucleus

    fig, ax = plt.subplots(1, number_sub_figs, figsize=(10 * number_sub_figs, 10))

    # plot marker intensity for all subplots
    for i in range(number_sub_figs):
        ax[i].imshow(im_marker, cmap=plt.cm.gray, alpha=1.0)

    outlines_cell, outlines_mem = _get_outline_and_membrane_thickness(im_marker, cell_mask, parameters)

    # cell and membrane outline
    outlines_cell_ = np.where(outlines_cell > 0, 30, 0)  # todo: what is that threshold?
    ax[0].imshow(np.ma.masked_where(outlines_cell_ == 0, outlines_cell_), plt.cm.Wistia, vmin=0, vmax=100, alpha=0.5)

    outlines_mem_ = np.where(outlines_mem > 0, 30, 0)  # todo: what is that threshold?
    ax[1].imshow(np.ma.masked_where(outlines_mem_ == 0, outlines_mem_), plt.cm.Wistia, vmin=0, vmax=100, alpha=0.5)

    # nuclei marker intensity
    if nuclei_mask is not None:
        outline_nuc = get_outline_from_mask(nuclei_mask, parameters["outline_width"])
        outline_nuc_ = np.where(outline_nuc == True, 30, 0)  # todo: 30 ?
        ax[2].imshow(
            np.ma.masked_where(outline_nuc_ == 0, outline_nuc_), plt.cm.Wistia, vmin=0, vmax=100, alpha=0.75
        )  # always last axis

    # plot mean expression value of cell and membrane as text
    for index, row in single_cell_dataset.iterrows():
        ax[0].text(row["Y_cell"], row["X_cell"], str(np.round(row["mean_expression"], 1)), color="w", fontsize=7)
        ax[1].text(row["Y_cell"], row["X_cell"], str(np.round(row["mean_expression_mem"], 1)), color="w", fontsize=7)
        if nuclei_mask is not None:
            ax[2].text(
                row["Y_nuc"], row["X_nuc"], str(np.round(row["mean_expression_nuc"], 1)), color="w", fontsize=7
            )

    # set title
    ax[0].set_title("mean intensity cell")
    ax[1].set_title("mean intensity membrane")
    if nuclei_mask is not None:
        ax[2].set_title("mean intensity nucleus")

    # save output
    if "pdf" in parameters["graphics_output_format"]:
        plt.savefig(str(Path(output_path).joinpath(filename + "_marker_expression.pdf")))
    if "svg" in parameters["graphics_output_format"]:
        plt.savefig(str(Path(output_path).joinpath(filename + "_marker_expression.svg")))
    if "png" in parameters["graphics_output_format"]:
        plt.savefig(str(Path(output_path).joinpath(filename + "_marker_expression.png")))

    plt.close(fig)


def plot_marker_polarity(parameters, im_marker, cell_mask, single_cell_props, filename, output_path):
    fig, ax = plt.subplots(1, figsize=(10, 10))
    ax.imshow(im_marker, cmap=plt.cm.gray, alpha=1.0)

    outlines_cell = np.zeros((im_marker.shape[0], im_marker.shape[1]))

    for cell_label in np.unique(cell_mask):
        # exclude background
        if cell_label == 0:
            continue

        single_cell_mask = np.where(cell_mask == cell_label, 1, 0)
        outline_cell = get_outline_from_mask(single_cell_mask, parameters["outline_width"])
        outline_cell_ = np.where(outline_cell == True, 30, 0)
        outlines_cell += outline_cell_

    outlines_cell_ = np.where(outlines_cell > 0, 30, 0)
    ax.imshow(np.ma.masked_where(outlines_cell_ == 0, outlines_cell_), plt.cm.Wistia, vmin=0, vmax=100, alpha=0.5)

    for index, row in single_cell_props.iterrows():
        ax.plot(row["Y_cell"], row["X_cell"], '.m', markersize=1)
        ax.plot(row["Y_weighted"], row["X_weighted"], '.m', markersize=1)
        ax.arrow(row["Y_cell"], row["X_cell"], row["Y_weighted"] - row["Y_cell"], row["X_weighted"] - row["X_cell"],
                 color='white', width=2)

    ax.set_title("marker polarity")
    if "pdf" in parameters["graphics_output_format"]:
        plt.savefig(str(Path(output_path).joinpath(filename + "_marker_polarity.pdf")))
    if "svg" in parameters["graphics_output_format"]:
        plt.savefig(str(Path(output_path).joinpath(filename + "_marker_polarity.svg")))
    if "png" in parameters["graphics_output_format"]:
        plt.savefig(str(Path(output_path).joinpath(filename + "_marker_polarity.png")))

    plt.close(fig)


def plot_shape_props(parameters, im_junction, single_cell_props, filename, output_path, cell_mask, feature_to_plot, nuclei_mask=None):
    """ function to plot cell (and optionally nuclei) orientation/alignment 
    
    parameters  :   dict
                    user defined parameters
    im_junction :   numpy.array (2-dim), float
                    channel containing the junction staining (used for segmentation)
    masks       :   numpy.array (2-dim), int
                    same dimension as im_junction, contains cell masks
    """

    print("Plot: ", feature_to_plot)

    number_sub_figs = 1

    if nuclei_mask is not None:
        nuclei_mask = nuclei_mask.astype(bool)
        number_sub_figs = 2

    fig, ax = plt.subplots(1, number_sub_figs, figsize=(number_sub_figs * 5, 5))


    if feature_to_plot == 'shape_orientation':
        v_min = 0.0         
        v_max = 180.0
        yticks = [0.0,45.0,90.0,135.0,180.0]
        plot_title = "cell shape orientation"
        color_bar_label = "cell shape orientation (degree)"    
    else:
        v_min = 0.0
        v_max = 1.0
        yticks = [0.0,0.5,1.0]
        plot_title = "cell elongation"
        color_bar_label = "eccentricity"    

    cell_eccentricity = np.zeros((cell_mask.shape[0], cell_mask.shape[1]))

    for index, row in single_cell_props.iterrows():
        row_label = int(row['label'])
        if row_label == 0:
            continue
        if feature_to_plot == 'shape_orientation':     
            single_cell_mask = np.where(cell_mask == row_label, 1, 0) * row['shape_orientation']*180.0/np.pi
        else:
            single_cell_mask = np.where(cell_mask == row_label, 1, 0) * row['eccentricity']
        cell_eccentricity += single_cell_mask

    if nuclei_mask is not None:
        ax[0].imshow(im_junction, cmap=plt.cm.gray, alpha=1.0)

        if feature_to_plot == 'shape_orientation':   
            cax_0 = ax[0].imshow(
                np.ma.masked_where(cell_mask == 0, cell_eccentricity), cmap=cm.cm.phase, vmin=v_min, vmax=v_max, alpha=0.5
            )
        else:
            cax_0 = ax[0].imshow(
                np.ma.masked_where(cell_mask == 0, cell_eccentricity), cmap=plt.cm.bwr, vmin=v_min, vmax=v_max, alpha=0.5
            )
        
        color_bar = fig.colorbar(cax_0, ax=ax[0], shrink=0.3)  # , extend='both')
        color_bar.set_label(feature_to_plot)
        color_bar.ax.set_yticks(yticks)

        #color_bar = fig.colorbar(cax_0, ax=ax[0], shrink=0.3)
        #color_bar.set_label('eccentricity')
 
        ax[1].imshow(im_junction, cmap=plt.cm.gray, alpha=1.0)

        get_logger().info("Calculating nuclei eccentricity...")
        nuclei_eccentricity = np.zeros((cell_mask.shape[0], cell_mask.shape[1]))
        for index, row in single_cell_props.iterrows():
            row_label = int(row['label'])

            # exclude background
            if row_label == 0:
                continue

            single_cell_mask = np.where(cell_mask == row_label, True, 0) 
            single_nuclei_mask_ = np.logical_and(single_cell_mask, nuclei_mask)
            if feature_to_plot == 'shape_orientation':     
                single_nuclei_mask = np.where(single_nuclei_mask_ == True, 1, 0) * row['shape_orientation_nuc']*180.0/np.pi
            else:
                single_nuclei_mask = np.where(single_nuclei_mask_ == True, 1, 0) * row['eccentricity_nuc']
            #single_nuclei_mask = np.where(single_nuclei_mask_ == True, 1, 0) * row['eccentricity_nuc']
            nuclei_eccentricity += single_nuclei_mask

        get_logger().info("Maximal nuclei eccentricity: %s" % str(np.max(nuclei_eccentricity)))
        get_logger().info("Minimal nuclei eccentricity: %s" % str(np.min(nuclei_eccentricity)))

        nuclei_mask_ = np.where(nuclei_mask == True, 70, 0)  # todo: threshold?

        if feature_to_plot == 'shape_orientation': 
            cax_1 = ax[1].imshow(np.ma.masked_where(
                nuclei_mask_ == 0, nuclei_eccentricity), cmap =cm.cm.phase, vmin=v_min, vmax=v_max, alpha=0.75
            )
        else:  
            cax_1 = ax[1].imshow(np.ma.masked_where(
                nuclei_mask_ == 0, nuclei_eccentricity), cmap =plt.cm.bwr, vmin=v_min, vmax=v_max, alpha=0.5
            )

        #color_bar = fig.colorbar(cax_1, ax=ax[1], shrink=0.3)  # , extend='both')
        #color_bar.set_label('eccentricity')
        color_bar = fig.colorbar(cax_1, ax=ax[1], shrink=0.3)  # , extend='both')
        color_bar.set_label(color_bar_label)
        color_bar.ax.set_yticks(yticks)


    else:
        ax.imshow(im_junction, cmap=plt.cm.gray, alpha=1.0)

        if feature_to_plot == 'shape_orientation':
            v_min = 0.0         
            v_max = 180.0
            yticks = [0.0,45.0,90.0,135.0,180.0]
        else:
            v_min = 0.0
            v_max = 1.0
            yticks = [0.0,0.5,1.0]
        
        if feature_to_plot == 'shape_orientation':   
            cax = ax.imshow(
                np.ma.masked_where(cell_mask == 0, cell_eccentricity), cmap=cm.cm.phase, vmin=v_min, vmax=v_max, alpha=0.75
            )
        else:
            cax = ax.imshow(
                np.ma.masked_where(cell_mask == 0, cell_eccentricity), cmap=plt.cm.bwr, vmin=v_min, vmax=v_max, alpha=0.75
            )

        #ax.imshow(np.ma.masked_where(cell_mask == 0, cell_eccentricity), cmap=plt.cm.bwr, alpha=0.25)
        
        color_bar = fig.colorbar(cax, ax=ax, shrink=0.3)  # , extend='both')
        color_bar.set_label(color_bar_label)
        color_bar.ax.set_yticks(yticks)

    for index, row in single_cell_props.iterrows():
        orientation = row['shape_orientation']
        x0 = row['X_cell']
        y0 = row['Y_cell']

        x1_major = x0 + math.sin(orientation) * 0.5 * row['major_axis_length']
        y1_major = y0 - math.cos(orientation) * 0.5 * row['major_axis_length']
        x2_major = x0 - math.sin(orientation) * 0.5 * row['major_axis_length']
        y2_major = y0 + math.cos(orientation) * 0.5 * row['major_axis_length']

        x1_minor = x0 - math.cos(orientation) * 0.5 * row['minor_axis_length']
        y1_minor = y0 - math.sin(orientation) * 0.5 * row['minor_axis_length']
        x2_minor = x0 + math.cos(orientation) * 0.5 * row['minor_axis_length']
        y2_minor = y0 + math.sin(orientation) * 0.5 * row['minor_axis_length']

        if nuclei_mask is not None:
            # plot orientation degree
            ax[0].plot((y1_major, y2_major), (x1_major, x2_major), '--w', linewidth=0.5)
            ax[0].plot((y1_minor, y2_minor), (x1_minor, x2_minor), '--w', linewidth=0.5)
            ax[0].plot(y0, x0, '.b', markersize=MARKERSIZE)
            orientation_degree = 180.0 * orientation / np.pi
            ax[0].text(y0, x0, str(int(np.round(orientation_degree, 0))), color="yellow", fontsize=FONTSIZE_TEXT_ANNOTATIONS)

            # plot orientation degree nucleus
            orientation = row['shape_orientation_nuc']
            x0 = row['X_nuc']
            y0 = row['Y_nuc']

            x1_major = x0 + math.sin(orientation) * 0.5 * row['major_axis_length_nuc']
            y1_major = y0 - math.cos(orientation) * 0.5 * row['major_axis_length_nuc']
            x2_major = x0 - math.sin(orientation) * 0.5 * row['major_axis_length_nuc']
            y2_major = y0 + math.cos(orientation) * 0.5 * row['major_axis_length_nuc']

            x1_minor = x0 - math.cos(orientation) * 0.5 * row['minor_axis_length_nuc']
            y1_minor = y0 - math.sin(orientation) * 0.5 * row['minor_axis_length_nuc']
            x2_minor = x0 + math.cos(orientation) * 0.5 * row['minor_axis_length_nuc']
            y2_minor = y0 + math.sin(orientation) * 0.5 * row['minor_axis_length_nuc']

            ax[1].plot((y1_major, y2_major), (x1_major, x2_major), '--w', linewidth=0.5)
            ax[1].plot((y1_minor, y2_minor), (x1_minor, x2_minor), '--w', linewidth=0.5)
            ax[1].plot(y0, x0, '.b', markersize=MARKERSIZE)
            orientation_degree = 180.0 * orientation / np.pi
            ax[1].text(y0, x0, str(int(np.round(orientation_degree, 0))), color="yellow", fontsize=FONTSIZE_TEXT_ANNOTATIONS)
        else:
            ax.plot((y1_major, y2_major), (x1_major, x2_major), '--w', linewidth=0.5)
            ax.plot((y1_minor, y2_minor), (x1_minor, x2_minor), '--w', linewidth=0.5)
            ax.plot(y0, x0, '.b', markersize=MARKERSIZE)
            if feature_to_plot == "shape_orientation":
                orientation_degree = 180.0 * orientation / np.pi
                ax.text(y0, x0, str(int(np.round(orientation_degree, 0))), color="yellow", fontsize=FONTSIZE_TEXT_ANNOTATIONS)
            else:
                ax.text(y0, x0, str(np.round(row["eccentricity"], 2)), color="yellow", fontsize=FONTSIZE_TEXT_ANNOTATIONS)

    # set title and ax limits
    if nuclei_mask is not None:
        ax[0].set_title(plot_title)
        ax[0].set_xlim(0, im_junction.shape[1])
        ax[0].set_ylim(0, im_junction.shape[0])
        ax[0].invert_yaxis()
        ax[0].axis('off')
        
        ax[1].set_title(plot_title)
        ax[1].set_xlim(0, im_junction.shape[1])
        ax[1].set_ylim(0, im_junction.shape[0])
        ax[1].invert_yaxis()
        ax[1].axis('off')
    else:
        ax.set_title(plot_title)
        ax.set_xlim(0, im_junction.shape[1])
        ax.set_ylim(0, im_junction.shape[0])
        ax.invert_yaxis()
        ax.axis('off')

    # save to disk
    plt.tight_layout()
    
    if "pdf" in parameters["graphics_output_format"]:
        plt.savefig(str(Path(output_path).joinpath(filename + "_" + feature_to_plot + ".pdf")))
    if "svg" in parameters["graphics_output_format"]:
        plt.savefig(str(Path(output_path).joinpath(filename + "_" + feature_to_plot + ".svg")))
    if "png" in parameters["graphics_output_format"]:
        plt.savefig(str(Path(output_path).joinpath(filename + "_" + feature_to_plot + ".png")))
    plt.close(fig)


def plot_ratio_method(parameters, im_junction, cell_mask, single_cell_props, filename, output_path):
    fig, ax = plt.subplots(1, 1)

    ax.imshow(im_junction, cmap=plt.cm.gray, alpha=1.0)
    ax.imshow(cell_mask, cmap=plt.cm.Set3, alpha=0.25)

    outlines_cell = np.zeros((im_junction.shape[0], im_junction.shape[1]))

    for cell_label in np.unique(cell_mask):

        # exclude background
        if cell_label == 0:
            continue

        single_cell_mask = np.where(cell_mask == cell_label, 1, 0)
        outline_cell = get_outline_from_mask(single_cell_mask, parameters["membrane_thickness"])
        outline_cell_ = np.where(outline_cell == True, 30, 0)
        outlines_cell += outline_cell_

    outlines_cell_ = np.where(outlines_cell > 0, 30, 0)
    ax.imshow(np.ma.masked_where(outlines_cell_ == 0, outlines_cell_), plt.cm.bwr, vmin=0, vmax=100, alpha=0.5)

    for index, row in single_cell_props.iterrows():
        x0 = row['X_cell']
        y0 = row['Y_cell']

        x1 = x0 + math.sin(np.pi / 4.0) * 0.5 * row['major_axis_length']
        y1 = y0 + math.cos(np.pi / 4.0) * 0.5 * row['major_axis_length']
        x2 = x0 + math.cos(np.pi / 4.0) * 0.5 * row['major_axis_length']
        y2 = y0 - math.sin(np.pi / 4.0) * 0.5 * row['major_axis_length']

        ax.plot((y0, y1), (x0, x1), '--r', linewidth=0.5)
        ax.plot((y0, y2), (x0, x2), '--r', linewidth=0.5)
        ax.plot(y0, x0, '.b', markersize=5)

        x1 = x0 - math.sin(np.pi / 4.0) * 0.5 * row['major_axis_length']
        y1 = y0 - math.cos(np.pi / 4.0) * 0.5 * row['major_axis_length']
        x2 = x0 - math.cos(np.pi / 4.0) * 0.5 * row['major_axis_length']
        y2 = y0 + math.sin(np.pi / 4.0) * 0.5 * row['major_axis_length']

        ax.plot((y0, y1), (x0, x1), '--r', linewidth=0.5)
        ax.plot((y0, y2), (x0, x2), '--r', linewidth=0.5)
        ax.plot(y0, x0, '.b', markersize=5)

        ax.set_title("ratio method")
        ax.set_xlim(0, im_junction.shape[0])
        ax.set_ylim(0, im_junction.shape[1])

    plt.tight_layout()
    
    if "pdf" in parameters["graphics_output_format"]:
        plt.savefig(str(Path(output_path).joinpath(filename + "_ratio_method.pdf")))
    if "svg" in parameters["graphics_output_format"]:
        plt.savefig(str(Path(output_path).joinpath(filename + "_ratio_method.svg")))
    if "png" in parameters["graphics_output_format"]:
        plt.savefig(str(Path(output_path).joinpath(filename + "_ratio_method.png")))

    plt.close(fig)

## TODO: remove this function, has been integrated in plot_shape_props
#def plot_cyclical(parameters, input_image, cell_masks_approx, filename, output_path):
#
#
#    get_logger().info("Plotting shape orientation on a cyclic scale")
#
#
#    #lab_image = label(input_image)
#    lab_image = label(cell_masks_approx)
#    regions = regionprops(lab_image)
#    orient_list = []
#    for region in regions:
#        cell_masks_approx[lab_image == region.label] = region.orientation * (180 / np.pi) + 90
#        orient_list.append(np.round(math.degrees(region.orientation)))
#
#    #cell_masks_approx[cell_masks_approx == 0] = np.nan
#    cell_masks_approx= np.ma.masked_where(cell_masks_approx == 0, cell_masks_approx)
#
#    masked_data = (lab_image != 0)
#    t_f = (masked_data * 1)
#    masked_data[np.where(t_f == 0), np.where(t_f == 0)] = np.nan
#
#    fig = plt.figure(figsize=(20, 20))
#    cm_phase = plt.imshow(cell_masks_approx, cmap=cm.cm.phase)  # todo: color scale OK? Put on top of file.
#    plt.imshow(np.where(masked_data == 0, 1, np.nan), cmap='binary', vmin=0, vmax=1)
#    #plt.imshow(np.where(masked_data == 0, 1, np.nan), cmap='binary', vmin=0, vmax=1)
#    plt.colorbar(cm_phase)
#
#    plt.savefig(str(Path(output_path).joinpath(filename + "_cyclic.pdf")))
#    plt.savefig(str(Path(output_path).joinpath(filename + "_cyclic.png")))
#    plt.close(fig)
#
#    #return cm_phase


def plot_adjacency_matrix(label_image, intensity_image):
    # todo: needed?
    rag = orientation_graph_nf(label_image)
    out = graph.draw_rag(label_image, rag, intensity_image, node_color="#ffde00")
    return out


def plot_dataset(parameters, img, properties_dataset, output_path, filename, cell_mask_rem_island, nuclei_mask,
                 golgi_mask, im_marker):
    """Plots the properties dataset"""
    get_logger().info("Plotting data...")
    im_junction = img[:, :, int(parameters["channel_junction"])]

    # TODO: adapt name plot polarity in parameter files
    if parameters["plot_polarity"] and nuclei_mask is not None and golgi_mask is not None:
        plot_organelle_polarity(
            parameters,
            im_junction,
            cell_mask_rem_island,
            nuclei_mask,
            golgi_mask,
            properties_dataset,
            filename, output_path,
        )
    if parameters["plot_marker"] and im_marker is not None:
        plot_marker_expression(
            parameters, im_marker, cell_mask_rem_island, properties_dataset, filename, output_path,
            nuclei_mask=nuclei_mask
        )
        plot_marker_polarity(
            parameters, im_marker, cell_mask_rem_island, properties_dataset, filename, output_path
        )
    # TODO: rename alignment to orientation
    if parameters["plot_alignment"]:
        plot_shape_props(
            parameters,
            im_junction,
            properties_dataset,
            filename,
            output_path,
            cell_mask_rem_island,
            "eccentricity",
            nuclei_mask=nuclei_mask,
        )
    if parameters["plot_ratio_method"]:
        plot_ratio_method(
            parameters,
            im_junction,
            cell_mask_rem_island,
            properties_dataset,
            filename,
            output_path
        )
    if parameters["plot_cyclic_orientation"]:
        plot_shape_props(
            parameters,
            im_junction,
            properties_dataset,
            filename,
            output_path,
            cell_mask_rem_island,
            "shape_orientation",
            nuclei_mask=nuclei_mask,
        )

        #plot_cyclical(
        #    parameters,
        #    im_junction,
        #    cell_mask_rem_island,
        #    filename,
        #    output_path
        #)
