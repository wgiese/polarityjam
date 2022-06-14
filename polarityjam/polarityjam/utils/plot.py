import math
from pathlib import Path

import cmocean as cm
import matplotlib as mpl
import numpy as np
import tifffile
from matplotlib import pyplot as plt
from skimage.future import graph

from polarityjam.polarityjam_logging import get_logger
from polarityjam.utils.masks import get_single_cell_mask, get_outline_from_mask
from polarityjam.utils.rag import orientation_graph_nf

# for figure plot resolution
FIGURE_DPI = 300
FONTSIZE_TEXT_ANNOTATIONS = 3
MARKERSIZE = 2
ALPHA_MASKS = 0.5


def save_current_fig(graphics_output_format, output_path, filename, filename_suffix, image=None):
    if "pdf" in graphics_output_format:
        plt.savefig(str(Path(output_path).joinpath(filename + filename_suffix + ".pdf")))
    if "svg" in graphics_output_format:
        plt.savefig(str(Path(output_path).joinpath(filename + filename_suffix + ".svg")))
    if "png" in graphics_output_format:
        plt.savefig(str(Path(output_path).joinpath(filename + filename_suffix + ".png")))
    if "tif" in graphics_output_format and image is not None:
        tifffile.imwrite(str(Path(output_path).joinpath(filename + filename_suffix + ".tif")), image)


def set_figure_dpi():
    mpl.rcParams['figure.dpi'] = FIGURE_DPI


def _add_single_cell_polarity_vector(ax, x_pos_p1, y_pos_p1, x_pos_p2, y_pos_p2):
    ax.plot(y_pos_p1, x_pos_p1, '.g', markersize=MARKERSIZE)
    ax.plot(y_pos_p2, x_pos_p2, '.m', markersize=MARKERSIZE)
    ax.arrow(
        y_pos_p1,
        x_pos_p1,
        y_pos_p2 - y_pos_p1,
        x_pos_p2 - x_pos_p1,
        color='white', width=2
    )


def _get_outline_and_membrane_thickness(im_marker, cell_mask, parameters):
    outlines_cell = np.zeros((im_marker.shape[0], im_marker.shape[1]))
    outlines_mem = np.zeros((im_marker.shape[0], im_marker.shape[1]))

    for cell_label in np.unique(cell_mask):
        # exclude background
        if cell_label == 0:
            continue

        single_cell_mask = get_single_cell_mask(cell_label, cell_mask)
        outline_cell = get_outline_from_mask(single_cell_mask, parameters["outline_width"])
        outline_cell_ = np.where(outline_cell == True, 30, 0)  # todo: magic number 30?
        outlines_cell += outline_cell_

        outline_mem = get_outline_from_mask(single_cell_mask, parameters["membrane_thickness"])
        outline_mem_ = np.where(outline_mem == True, 30, 0)  # todo: magic number 30?
        outlines_mem += outline_mem_

    return [outlines_cell, outlines_mem]


def plot_seg_channels(seg_img, output_path, filename):
    """Plots the separate channels from the input file given."""
    output_path = Path(output_path)
    filename_out = str(output_path.joinpath(filename + "_seg.png"))
    if len(seg_img.shape) > 2:
        fig, ax = plt.subplots(1, 2)
        ax[0].imshow(seg_img[0, :, :])
        ax[0].set_title("junction channel")
        ax[1].imshow(seg_img[1, :, :])
        ax[1].set_title("nuclei channel")
    else:
        fig, ax = plt.subplots()
        ax.imshow(seg_img[:, :])
    plt.savefig(filename_out)
    plt.close(fig)


def plot_cellpose_masks(seg_img, cellpose_mask, output_path, filename):
    """Plots the cellpose segmentation output, together with the separate channels from the input image."""
    output_path = Path(output_path)
    filename_out = str(output_path.joinpath(filename + "_cellpose_seg.png"))

    if len(seg_img.shape) > 2:
        fig, ax = plt.subplots(1, 3)
        ax[0].imshow(seg_img[0, :, :])
        ax[0].set_title("junction channel")
        ax[1].imshow(seg_img[1, :, :])
        ax[1].set_title("nuclei channel")
        ax[2].imshow(seg_img[0, :, :])
        ax[2].imshow(cellpose_mask, cmap=plt.cm.Set3, alpha=0.5)
        ax[2].set_title("cellpose segmentation")
    else:
        fig, ax = plt.subplots(1, 2)
        ax[0].imshow(seg_img[:, :])
        ax[0].set_title("junction channel")
        ax[1].imshow(cellpose_mask, cmap=plt.cm.Set3, alpha=0.5)
        ax[1].set_title("cellpose segmentation")
    plt.savefig(filename_out)
    plt.close(fig)


def plot_organelle_polarity(parameters, im_junction, cell_mask, nuclei_mask, golgi_mask, single_cell_props,
                            base_filename, output_path):
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
    # figure and axes
    width = parameters["graphics_width"]
    fig, ax = plt.subplots(figsize=(width, width))

    # base image
    ax.imshow(im_junction, cmap=plt.cm.gray, alpha=1.0)

    # determine polarity_angle
    cell_angle = np.zeros((cell_mask.shape[0], cell_mask.shape[1]))
    for index, row in single_cell_props.iterrows():
        row_label = int(row['label'])
        if row_label == 0:
            continue
        cell_angle += get_single_cell_mask(row_label, cell_mask) * row['angle_deg']
    polarity_angle = np.ma.masked_where(cell_mask == 0, cell_angle)

    # plot polarity angle
    cax = ax.imshow(polarity_angle, cmap=cm.cm.phase, vmin=0, vmax=360, alpha=0.5)
    color_bar = fig.colorbar(cax, ax=ax, shrink=0.3)  # , extend='both')
    color_bar.set_label("polarity angle")
    color_bar.ax.set_yticks([0, 90, 180, 270, 360])

    # plot differently colored golgi (red) and nuclei (blue)
    zero = np.zeros((im_junction.shape[0], im_junction.shape[1]))
    rgb_golgi = np.dstack((golgi_mask.astype(int) * 256, zero, zero, golgi_mask.astype(float) * 0.5))
    rgb_nuclei = np.dstack((zero, zero, nuclei_mask.astype(int) * 256, nuclei_mask.astype(float) * 0.5))
    ax.imshow(rgb_nuclei)
    ax.imshow(rgb_golgi)

    # plot polarity vector
    for index, row in single_cell_props.iterrows():
        _add_single_cell_polarity_vector(ax, row["X_nuc"], row["Y_nuc"], row["X_golgi"], row["Y_golgi"])
        if parameters["show_polarity_angles"]:
            ax.text(row["Y_cell"], row["X_cell"], str(int(np.round(row["angle_deg"], 0))), color="yellow", fontsize=6)

    # set ax limits
    ax.set_xlim(0, im_junction.shape[1])
    ax.set_ylim(0, im_junction.shape[0])
    ax.invert_yaxis()
    ax.axis('off')

    # save output & close
    save_current_fig(
        parameters["graphics_output_format"],
        output_path, base_filename,
        "_nuclei_golgi_vector",
        image=polarity_angle
    )
    plt.close(fig)


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

    # save output & close
    save_current_fig(parameters["graphics_output_format"], output_path, filename, "_marker_expression")
    plt.close(fig)


def plot_marker_polarity(parameters, im_marker, cell_mask, single_cell_props, filename, output_path):
    fig, ax = plt.subplots(1, figsize=(10, 10))
    ax.imshow(im_marker, cmap=plt.cm.gray, alpha=1.0)

    outlines_cell = np.zeros((im_marker.shape[0], im_marker.shape[1]))

    for cell_label in np.unique(cell_mask):
        # exclude background
        if cell_label == 0:
            continue

        single_cell_mask = get_single_cell_mask(cell_label, cell_mask)
        outline_cell = get_outline_from_mask(single_cell_mask, parameters["outline_width"])
        outline_cell_ = np.where(outline_cell == True, 30, 0)
        outlines_cell += outline_cell_

    outlines_cell_ = np.where(outlines_cell > 0, 30, 0)
    ax.imshow(np.ma.masked_where(outlines_cell_ == 0, outlines_cell_), plt.cm.Wistia, vmin=0, vmax=100, alpha=0.5)

    for index, row in single_cell_props.iterrows():
        _add_single_cell_polarity_vector(ax, row["X_cell"], row["Y_cell"], row["X_weighted"], row["Y_weighted"])

    ax.set_title("marker polarity")

    # save output & close
    save_current_fig(parameters["graphics_output_format"], output_path, filename, "_marker_polarity")
    plt.close(fig)


def _add_single_cell_eccentricity_axis(
        ax, y0, x0, orientation, major_axis_length, minor_axis_length, eccentricity
):
    x1_ma, x1_mi, x2_ma, x2_mi, y1_ma, y1_mi, y2_ma, y2_mi = _calc_single_cell_axis_orientation_vector(
        x0, y0, orientation, major_axis_length, minor_axis_length
    )

    ax.plot((y1_ma, y2_ma), (x1_ma, x2_ma), '--w', linewidth=0.5)
    ax.plot((y1_mi, y2_mi), (x1_mi, x2_mi), '--w', linewidth=0.5)
    ax.plot(y0, x0, '.b', markersize=MARKERSIZE)
    ax.text(y0, x0, str(np.round(eccentricity, 2)), color="yellow", fontsize=FONTSIZE_TEXT_ANNOTATIONS)


def _add_cell_eccentricity(fig, ax, im_junction, cell_mask, cell_eccentricity):
    v_min = 0.0
    v_max = 1.0
    yticks = [0.0, 0.5, 1.0]

    ax.imshow(im_junction, cmap=plt.cm.gray, alpha=1.0)

    # show cell_eccentricity everywhere but background label
    cax_0 = ax.imshow(
        np.ma.masked_where(cell_mask == 0, cell_eccentricity), cmap=plt.cm.bwr, vmin=v_min, vmax=v_max,
        alpha=0.5
    )

    # colorbar
    color_bar = fig.colorbar(cax_0, ax=ax, shrink=0.3)
    color_bar.set_label("eccentricity")
    color_bar.ax.set_yticks(yticks)


def _add_nuclei_eccentricity(fig, ax, im_junction, nuclei_mask, nuclei_eccentricity):
    v_min = 0.0
    v_max = 1.0
    yticks = [0.0, 0.5, 1.0]

    ax.imshow(im_junction, cmap=plt.cm.gray, alpha=1.0)

    # show nuclei eccentricity everywhere but background label
    nuclei_mask_ = np.where(nuclei_mask == True, 70, 0)  # todo: threshold?
    cax_1 = ax.imshow(np.ma.masked_where(
        nuclei_mask_ == 0, nuclei_eccentricity), cmap=plt.cm.bwr, vmin=v_min, vmax=v_max, alpha=0.5
    )

    # colorbar
    color_bar = fig.colorbar(cax_1, ax=ax, shrink=0.3)  # , extend='both')
    color_bar.set_label("eccentricity")
    color_bar.ax.set_yticks(yticks)


def _add_title(ax, plot_title, im_junction):
    ax.set_title(plot_title)
    ax.set_xlim(0, im_junction.shape[1])
    ax.set_ylim(0, im_junction.shape[0])
    ax.invert_yaxis()
    ax.axis('off')


def _calc_nuc_eccentricity(single_cell_props, cell_mask, nuclei_mask):
    get_logger().info("Calculating nuclei eccentricity...")
    nuclei_eccentricity = np.zeros((cell_mask.shape[0], cell_mask.shape[1]))
    for index, row in single_cell_props.iterrows():
        row_label = int(row['label'])

        # exclude background
        if row_label == 0:
            continue
        single_cell_mask = np.where(cell_mask == row_label, True, 0)
        single_nuclei_mask_ = np.logical_and(single_cell_mask, nuclei_mask)
        nuclei_eccentricity += np.where(single_nuclei_mask_ == True, 1, 0) * row['eccentricity_nuc']

    get_logger().info("Maximal nuclei eccentricity: %s" % str(np.max(nuclei_eccentricity)))
    get_logger().info("Minimal nuclei eccentricity: %s" % str(np.min(nuclei_eccentricity)))

    return nuclei_eccentricity


def _calc_cell_eccentricity(single_cell_props, cell_mask):
    get_logger().info("Calculating cell eccentricity...")
    cell_eccentricity = np.zeros((cell_mask.shape[0], cell_mask.shape[1]))
    for index, row in single_cell_props.iterrows():
        row_label = int(row['label'])
        if row_label == 0:
            continue
        cell_eccentricity += get_single_cell_mask(row_label, cell_mask) * row['eccentricity']

    get_logger().info("Maximal cell eccentricity: %s" % str(np.max(cell_eccentricity)))
    get_logger().info("Minimal cell eccentricity: %s" % str(np.min(cell_eccentricity)))

    return cell_eccentricity


def _calc_single_cell_axis_orientation_vector(x, y, orientation, major_axis_length, minor_axis_length):
    x1_major = x + math.sin(orientation) * 0.5 * major_axis_length
    y1_major = y - math.cos(orientation) * 0.5 * major_axis_length
    x2_major = x - math.sin(orientation) * 0.5 * major_axis_length
    y2_major = y + math.cos(orientation) * 0.5 * major_axis_length

    x1_minor = x - math.cos(orientation) * 0.5 * minor_axis_length
    y1_minor = y - math.sin(orientation) * 0.5 * minor_axis_length
    x2_minor = x + math.cos(orientation) * 0.5 * minor_axis_length
    y2_minor = y + math.sin(orientation) * 0.5 * minor_axis_length

    return [x1_major, x1_minor, x2_major, x2_minor, y1_major, y1_minor, y2_major, y2_minor]


def plot_eccentricity(parameters, im_junction, single_cell_props, filename, output_path, cell_mask, nuclei_mask=None):
    """ function to plot cell (and optionally nuclei) eccentricity

    parameters  :   dict
                    user defined parameters
    im_junction :   numpy.array (2-dim), float
                    channel containing the junction staining (used for segmentation)
    single_cell_props       :   pandas properties dataset
    filename    :   base filename for the output file
    output_path :   desired output path
    cell_mask   :   cellpose cell mask
    nuclei_mask :   cellpose nuclei mask
    """
    get_logger().info("Prepare plotting: eccentricity")

    # figure
    number_sub_figs = 1
    if nuclei_mask is not None:
        nuclei_mask = nuclei_mask.astype(bool)  # todo ?
        number_sub_figs = 2
    fig, ax = plt.subplots(1, number_sub_figs, figsize=(number_sub_figs * 5, 5))

    # get cell_eccentricity
    cell_eccentricity = _calc_cell_eccentricity(single_cell_props, cell_mask)

    # add cell (and nuclei) eccentricity to the figure
    if nuclei_mask is not None:
        _add_cell_eccentricity(fig, ax[0], im_junction, cell_mask, cell_eccentricity)
        # get nuclei eccentricity
        nuclei_eccentricity = _calc_nuc_eccentricity(single_cell_props, cell_mask, nuclei_mask)
        _add_nuclei_eccentricity(fig, ax[1], im_junction, nuclei_mask, nuclei_eccentricity)
    else:
        _add_cell_eccentricity(fig, ax, im_junction, cell_mask, cell_eccentricity)

    # plot major and minor axis
    for index, row in single_cell_props.iterrows():
        if nuclei_mask is not None:
            # plot orientation degree
            _add_single_cell_eccentricity_axis(
                ax[0],
                row['Y_cell'],
                row['X_cell'],
                row['shape_orientation'],
                row['major_axis_length'],
                row['minor_axis_length'],
                row["eccentricity"]
            )

            # plot orientation degree nucleus
            _add_single_cell_eccentricity_axis(
                ax[1],
                row['Y_nuc'],
                row['X_nuc'],
                row['shape_orientation_nuc'],
                row['major_axis_length_nuc'],
                row['minor_axis_length_nuc'],
                row["eccentricity_nuc"]
            )
        else:
            _add_single_cell_eccentricity_axis(
                ax,
                row['Y_cell'],
                row['X_cell'],
                row['shape_orientation'],
                row['major_axis_length'],
                row['minor_axis_length'],
                row["eccentricity"]
            )

    # set title and ax limits
    if nuclei_mask is not None:
        _add_title(ax[0], "cell elongation", im_junction)
        _add_title(ax[1], "nuclei elongation", im_junction)
    else:
        _add_title(ax, "cell elongation", im_junction)

    # save to disk
    plt.tight_layout()

    # save output & close
    save_current_fig(parameters["graphics_output_format"], output_path, filename, "_eccentricity")
    plt.close(fig)


def _calc_nuc_orientation(single_cell_props, cell_mask, nuclei_mask):
    get_logger().info("Calculating nuclei eccentricity...")
    nuclei_orientation = np.zeros((cell_mask.shape[0], cell_mask.shape[1]))
    for index, row in single_cell_props.iterrows():
        row_label = int(row['label'])

        # exclude background
        if row_label == 0:
            continue
        single_cell_mask = np.where(cell_mask == row_label, True, 0)
        single_nuclei_mask_ = np.logical_and(single_cell_mask, nuclei_mask)
        nuclei_orientation += np.where(single_nuclei_mask_ == True, 1, 0) * row[
            'shape_orientation_nuc'] * 180.0 / np.pi

    get_logger().info("Maximal nuclei orientation: %s" % str(np.max(nuclei_orientation)))
    get_logger().info("Minimal nuclei orientation: %s" % str(np.min(nuclei_orientation)))

    return nuclei_orientation


def _calc_cell_orientation(single_cell_props, cell_mask):
    get_logger().info("Calculating cell orientation...")
    cell_orientation = np.zeros((cell_mask.shape[0], cell_mask.shape[1]))
    for index, row in single_cell_props.iterrows():
        row_label = int(row['label'])
        if row_label == 0:
            continue
        cell_orientation += get_single_cell_mask(row_label, cell_mask) * row[
            'shape_orientation_nuc'] * 180.0 / np.pi

    get_logger().info("Maximal cell orientation: %s" % str(np.max(cell_orientation)))
    get_logger().info("Minimal cell orientation: %s" % str(np.min(cell_orientation)))

    return cell_orientation


def _add_cell_orientation(fig, ax, im_junction, cell_mask, cell_orientation):
    v_min = 0.0
    v_max = 180.0
    yticks = [0.0, 45.0, 90.0, 135.0, 180.0]

    ax.imshow(im_junction, cmap=plt.cm.gray, alpha=1.0)

    # show cell_eccentricity everywhere but background label
    cax = ax.imshow(
        np.ma.masked_where(cell_mask == 0, cell_orientation), cmap=cm.cm.phase, vmin=v_min, vmax=v_max,
        alpha=0.75
    )

    # colorbar
    color_bar = fig.colorbar(cax, ax=ax, shrink=0.3)
    color_bar.set_label("shape orientation (degree)")
    color_bar.ax.set_yticks(yticks)


def _add_nuclei_orientation(fig, ax, im_junction, nuclei_mask, nuclei_orientation):
    v_min = 0.0
    v_max = 180.0
    yticks = [0.0, 45.0, 90.0, 135.0, 180.0]

    ax.imshow(im_junction, cmap=plt.cm.gray, alpha=1.0)

    # show nuclei eccentricity everywhere but background label
    nuclei_mask_ = np.where(nuclei_mask == True, 70, 0)  # todo: threshold?
    cax_1 = ax.imshow(np.ma.masked_where(
        nuclei_mask_ == 0, nuclei_orientation), cmap=cm.cm.phase, vmin=v_min, vmax=v_max, alpha=0.75
    )

    # colorbar
    color_bar = fig.colorbar(cax_1, ax=ax, shrink=0.3)  # , extend='both')
    color_bar.set_label("shape orientation (degree)")
    color_bar.ax.set_yticks(yticks)


def _add_single_cell_orientation_degree_axis(
        ax, y0, x0, orientation, major_axis_length, minor_axis_length,
):
    x1_ma, x1_mi, x2_ma, x2_mi, y1_ma, y1_mi, y2_ma, y2_mi = _calc_single_cell_axis_orientation_vector(
        x0, y0, orientation, major_axis_length, minor_axis_length
    )
    orientation_degree = 180.0 * orientation / np.pi

    ax.plot((y1_ma, y2_ma), (x1_ma, x2_ma), '--w', linewidth=0.5)
    ax.plot((y1_mi, y2_mi), (x1_mi, x2_mi), '--w', linewidth=0.5)
    ax.plot(y0, x0, '.b', markersize=MARKERSIZE)
    ax.text(y0, x0, str(int(np.round(orientation_degree, 0))), color="yellow", fontsize=FONTSIZE_TEXT_ANNOTATIONS)


def plot_orientation(parameters, im_junction, single_cell_props, filename, output_path, cell_mask, nuclei_mask=None):
    """ function to plot cell (and optionally nuclei) orientation

    parameters  :   dict
                    user defined parameters
    im_junction :   numpy.array (2-dim), float
                    channel containing the junction staining (used for segmentation)
    single_cell_props       :   pandas properties dataset
    filename    :   base filename for the output file
    output_path :   desired output path
    cell_mask   :   cellpose cell mask
    nuclei_mask :   cellpose nuclei mask
    """
    get_logger().info("Prepare plotting: orientation")

    # figure
    number_sub_figs = 1
    if nuclei_mask is not None:
        nuclei_mask = nuclei_mask.astype(bool)  # todo ?
        number_sub_figs = 2
    fig, ax = plt.subplots(1, number_sub_figs, figsize=(number_sub_figs * 5, 5))

    # get cell_orientation
    cell_orientation = _calc_cell_orientation(single_cell_props, cell_mask)

    # add cell (and nuclei) orientation to the figure
    if nuclei_mask is not None:
        _add_cell_orientation(fig, ax[0], im_junction, cell_mask, cell_orientation)
        # get nuclei orientation
        nuclei_orientation = _calc_nuc_orientation(single_cell_props, cell_mask, nuclei_mask)
        _add_nuclei_orientation(fig, ax[1], im_junction, nuclei_mask, nuclei_orientation)
    else:
        _add_cell_orientation(fig, ax, im_junction, cell_mask, cell_orientation)

    # plot major and minor axis
    for index, row in single_cell_props.iterrows():
        if nuclei_mask is not None:
            # plot orientation degree
            _add_single_cell_orientation_degree_axis(
                ax[0],
                row['Y_cell'],
                row['X_cell'],
                row['shape_orientation'],
                row['major_axis_length'],
                row['minor_axis_length']
            )

            # plot orientation degree nucleus
            _add_single_cell_orientation_degree_axis(
                ax[1],
                row['Y_nuc'],
                row['X_nuc'],
                row['shape_orientation_nuc'],
                row['major_axis_length_nuc'],
                row['minor_axis_length_nuc']
            )
        else:
            _add_single_cell_orientation_degree_axis(
                ax,
                row['Y_cell'],
                row['X_cell'],
                row['shape_orientation'],
                row['major_axis_length'],
                row['minor_axis_length']
            )

    # set title and ax limits
    if nuclei_mask is not None:
        _add_title(ax[0], "cell shape orientation", im_junction)
        _add_title(ax[1], "nuclei shape orientation", im_junction)
    else:
        _add_title(ax, "cell shape orientation", im_junction)

    # save to disk
    plt.tight_layout()

    # save output & close
    save_current_fig(parameters["graphics_output_format"], output_path, filename, "_shape_orientation")
    plt.close(fig)


def plot_shape_props(parameters, im_junction, single_cell_props, filename, output_path, cell_mask, feature_to_plot,
                     nuclei_mask=None):
    """ function to plot cell (and optionally nuclei) orientation/alignment

    parameters  :   dict
                    user defined parameters
    im_junction :   numpy.array (2-dim), float
                    channel containing the junction staining (used for segmentation)
    masks       :   numpy.array (2-dim), int
                    same dimension as im_junction, contains cell masks
    """
    get_logger().info("Prepare plotting: %s" % feature_to_plot)


def plot_ratio_method(parameters, im_junction, cell_mask, single_cell_props, filename, output_path):
    fig, ax = plt.subplots(1, 1)

    ax.imshow(im_junction, cmap=plt.cm.gray, alpha=1.0)
    ax.imshow(cell_mask, cmap=plt.cm.Set3, alpha=0.25)

    cell_outlines = np.zeros((im_junction.shape[0], im_junction.shape[1]))

    for cell_label in np.unique(cell_mask):
        # exclude background
        if cell_label == 0:
            continue

        single_cell_mask = get_single_cell_mask(cell_label, cell_mask)
        cell_outline = get_outline_from_mask(single_cell_mask, parameters["membrane_thickness"])
        cell_outlines += np.where(cell_outline == True, 30, 0)

    outlines_cell_ = np.where(cell_outlines > 0, 30, 0)
    ax.imshow(np.ma.masked_where(outlines_cell_ == 0, outlines_cell_), plt.cm.bwr, vmin=0, vmax=100, alpha=0.5)

    for index, row in single_cell_props.iterrows():
        x0 = row['X_cell']
        y0 = row['Y_cell']

        # upper
        x1 = x0 + math.sin(np.pi / 4.0) * 0.5 * row['major_axis_length']
        y1 = y0 + math.cos(np.pi / 4.0) * 0.5 * row['major_axis_length']
        x2 = x0 + math.cos(np.pi / 4.0) * 0.5 * row['major_axis_length']
        y2 = y0 - math.sin(np.pi / 4.0) * 0.5 * row['major_axis_length']

        ax.plot((y0, y1), (x0, x1), '--r', linewidth=0.5)
        ax.plot((y0, y2), (x0, x2), '--r', linewidth=0.5)
        ax.plot(y0, x0, '.b', markersize=5)

        # lower
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

    save_current_fig(parameters["graphics_output_format"], output_path, filename, "_ratio_method")
    plt.close(fig)


def plot_adjacency_matrix(label_image, intensity_image):
    # todo: needed?
    rag = orientation_graph_nf(label_image)
    out = graph.draw_rag(label_image, rag, intensity_image, node_color="#ffde00")
    return out


def plot_dataset(parameters, img, properties_ds, output_path, filename, cell_mask, nuclei_mask, golgi_mask, im_marker):
    """Plots the properties dataset"""
    get_logger().info("Plotting...")
    im_junction = img[:, :, int(parameters["channel_junction"])]

    # TODO: adapt name plot polarity in parameter files
    if parameters["plot_polarity"] and nuclei_mask is not None and golgi_mask is not None:
        plot_organelle_polarity(
            parameters,
            im_junction,
            cell_mask,
            nuclei_mask,
            golgi_mask,
            properties_ds,
            filename, output_path,
        )
    if parameters["plot_marker"] and im_marker is not None:
        plot_marker_expression(
            parameters, im_marker, cell_mask, properties_ds, filename, output_path,
            nuclei_mask=nuclei_mask
        )
        plot_marker_polarity(
            parameters, im_marker, cell_mask, properties_ds, filename, output_path
        )
    if parameters["plot_orientation"]:
        plot_eccentricity(
            parameters,
            im_junction,
            properties_ds,
            filename,
            output_path,
            cell_mask,
            nuclei_mask=nuclei_mask,
        )
    if parameters["plot_ratio_method"]:
        plot_ratio_method(
            parameters,
            im_junction,
            cell_mask,
            properties_ds,
            filename,
            output_path
        )
    if parameters["plot_cyclic_orientation"]:
        plot_orientation(
            parameters,
            im_junction,
            properties_ds,
            filename,
            output_path,
            cell_mask,
            nuclei_mask=nuclei_mask,
        )
