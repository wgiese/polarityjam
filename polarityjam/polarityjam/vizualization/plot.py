import json
import math
from pathlib import Path

import cmocean as cm
import matplotlib as mpl
import numpy as np
import tifffile
from matplotlib import pyplot as plt

from polarityjam.model.masks import get_single_cell_mask, get_outline_from_mask
from polarityjam.polarityjam_logging import get_logger
from polarityjam.utils import parameters

# for figure plot resolution  # todo: parameters?
FIGURE_DPI = 300
FONTSIZE_TEXT_ANNOTATIONS = 3
MARKERSIZE = 2
ALPHA_MASKS = 0.5
CELL_OUTLINE_INTENSITY = 30


def save_current_fig(graphics_output_format, output_path, filename, filename_suffix, image=None):
    # prevent text outside figure area
    plt.tight_layout()

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


def _get_outline_and_membrane_thickness(im_marker, cell_mask):
    outlines_cell = np.zeros((im_marker.shape[0], im_marker.shape[1]))
    outlines_mem_accumulated = np.zeros((im_marker.shape[0], im_marker.shape[1]))

    for cell_label in np.unique(cell_mask):
        # exclude background
        if cell_label == 0:
            continue

        single_cell_mask = get_single_cell_mask(cell_label, cell_mask)
        outline_cell = get_outline_from_mask(single_cell_mask, parameters.outline_width)
        outline_cell_ = np.where(outline_cell == True, 1, 0)
        outlines_cell += outline_cell_

        outline_mem = get_outline_from_mask(single_cell_mask, parameters.membrane_thickness)
        outline_mem_ = np.where(outline_mem == True, 1, 0)
        outlines_mem_accumulated += outline_mem_

    return [outlines_cell, outlines_mem_accumulated]


def plot_seg_channels(seg_img, output_path, filename):
    """Plots the separate channels from the input file given."""
    get_logger().info("Plotting: input channels")

    output_path = Path(output_path)
    filename_out = str(output_path.joinpath(filename + "_seg.png"))
    if len(seg_img.shape) > 2:
        fig, ax = plt.subplots(1, 2)
        if not parameters.show_graphics_axis:
            ax[0].axis('off')
            ax[1].axis('off')
        ax[0].imshow(seg_img[0, :, :])
        ax[0].set_title("junction channel")
        ax[1].imshow(seg_img[1, :, :])
        ax[1].set_title("nuclei channel")
    else:
        fig, ax = plt.subplots()
        if not parameters.show_graphics_axis:
            ax.axis('off')
        ax.imshow(seg_img[:, :])
    plt.savefig(filename_out)
    plt.close(fig)


def plot_cellpose_masks(seg_img, cellpose_mask, output_path, filename):
    """Plots the cellpose segmentation output, together with the separate channels from the input image."""
    get_logger().info("Plotting: cellpose masks")

    # figure and axes
    w, h = parameters.graphics_width, parameters.graphics_height

    if len(seg_img.shape) > 2:
        fig, ax = plt.subplots(1, 3, figsize=(3 * w, h))
        ax[0].imshow(seg_img[0, :, :])
        ax[0].set_title("junction channel")
        ax[1].imshow(seg_img[1, :, :])
        ax[1].set_title("nuclei channel")
        ax[2].imshow(seg_img[0, :, :])
        ax[2].imshow(cellpose_mask, cmap=plt.cm.Set3, alpha=0.5)
        ax[2].set_title("cellpose segmentation")
    else:
        fig, ax = plt.subplots(1, 2, figsize=(2 * w, h))
        ax[0].imshow(seg_img[:, :])
        ax[0].set_title("junction channel")
        ax[1].imshow(cellpose_mask, cmap=plt.cm.Set3, alpha=0.5)
        ax[1].set_title("cellpose segmentation")

    if not parameters.show_graphics_axis:
        for ax_ in ax:
            ax_.axis('off')

    # save output & close
    save_current_fig(
        parameters.graphics_output_format,
        output_path, filename,
        "_cellpose_seg",
    )
    plt.close(fig)


def plot_organelle_polarity(im_junction, cell_mask, nuclei_mask, organelle_mask, single_cell_props,
                            base_filename, output_path):
    """ function to plot nuclei-organelle polarity vectors

    im_junction :   numpy.array (2-dim), float
                    channel containing the junction staining (used for segmentation)
    cell_mask   :   numpy.array (2-dim), int
                    same dimension as im_junction, contains cell masks
    nuclei_mask :   numpy.array (2-dim), int
                    same dimension as im_junction, contains nuclei masks
    organelle_mask  :   numpy.array (2-dim), int
                    same dimension as im_junction, contains organelle masks
    single_cell_props : pandas data frame
    base_filename : string 
                    base_filename for plots
    output_path :   string
                    output path for plots 
    """
    get_logger().info("Plotting: organelle polarity")

    # figure and axes
    w, h = parameters.graphics_width, parameters.graphics_height
    fig, ax = plt.subplots(figsize=(w, h))

    # resources image
    ax.imshow(im_junction, cmap=plt.cm.gray, alpha=1.0)

    # determine polarity_angle
    cell_angle = np.zeros((cell_mask.shape[0], cell_mask.shape[1]))
    for index, row in single_cell_props.iterrows():
        row_label = int(row['label'])
        if row_label == 0:
            continue
        cell_angle += get_single_cell_mask(row_label, cell_mask) * row['organelle_orientation_deg']
    polarity_angle = np.ma.masked_where(cell_mask == 0, cell_angle)

    # plot polarity angle
    cax = ax.imshow(polarity_angle, cmap=cm.cm.phase, vmin=0, vmax=360, alpha=0.5)
    color_bar = fig.colorbar(cax, ax=ax, shrink=0.3)  # , extend='both')
    color_bar.set_label("polarity angle")
    color_bar.ax.set_yticks([0, 90, 180, 270, 360])

    # plot differently colored organelle (red) and nuclei (blue)
    zero = np.zeros((im_junction.shape[0], im_junction.shape[1]))
    rgb_organelle = np.dstack((organelle_mask.astype(int) * 256, zero, zero, organelle_mask.astype(float) * 0.5))
    rgb_nuclei = np.dstack((zero, zero, nuclei_mask.astype(int) * 256, nuclei_mask.astype(float) * 0.5))
    ax.imshow(rgb_nuclei)
    ax.imshow(rgb_organelle)

    # plot polarity vector
    for index, row in single_cell_props.iterrows():
        _add_single_cell_polarity_vector(ax, row["nuc_X"], row["nuc_Y"], row["organelle_X"], row["organelle_Y"])
        if parameters.show_polarity_angles:
            ax.text(row["cell_Y"], row["cell_X"], str(int(np.round(row["organelle_orientation_deg"], 0))),
                    color="yellow", fontsize=6)

    # set title and ax limits
    _add_title(ax, "organelle orientation", im_junction, parameters.show_graphics_axis)

    # save output & close
    save_current_fig(
        parameters.graphics_output_format,
        output_path, base_filename,
        "_nuclei_organelle_vector",
        image=polarity_angle
    )
    plt.close(fig)


def plot_nuc_displacement_orientation(im_junction, cell_mask, nuclei_mask, single_cell_props,
                                      base_filename, output_path):
    """ function to plot nuclei-organelle polarity vectors

    im_junction :   numpy.array (2-dim), float
                    channel containing the junction staining (used for segmentation)
    cell_mask   :   numpy.array (2-dim), int
                    same dimension as im_junction, contains cell masks
    nuclei_mask :   numpy.array (2-dim), int
                    same dimension as im_junction, contains nuclei masks
    single_cell_props : pandas data frame
    base_filename : string
                    base_filename for plots
    output_path :   string
                    output path for plots
    """
    get_logger().info("Plotting: marker nucleus polarity")

    # figure and axes
    w, h = parameters.graphics_width, parameters.graphics_height
    fig, ax = plt.subplots(figsize=(w, h))

    # resources image
    ax.imshow(im_junction, cmap=plt.cm.gray, alpha=1.0)

    # determine nucleus polarity_angle
    nuc_angle = np.zeros((cell_mask.shape[0], cell_mask.shape[1]))
    for index, row in single_cell_props.iterrows():
        row_label = int(row['label'])
        if row_label == 0:
            continue
        nuc_angle += get_single_cell_mask(row_label, cell_mask) * row["nuc_displacement_orientation_deg"]
    nuc_polarity_angle = np.ma.masked_where(cell_mask == 0, nuc_angle)

    # plot polarity angle
    cax = ax.imshow(nuc_polarity_angle, cmap=cm.cm.phase, vmin=0, vmax=360, alpha=0.5)
    color_bar = fig.colorbar(cax, ax=ax, shrink=0.3)  # , extend='both')
    color_bar.set_label("polarity angle")
    color_bar.ax.set_yticks([0, 90, 180, 270, 360])

    # plot nuclei (blue)
    zero = np.zeros((im_junction.shape[0], im_junction.shape[1]))
    rgb_nuclei = np.dstack((zero, zero, nuclei_mask.astype(int) * 256, nuclei_mask.astype(float) * 0.5))
    ax.imshow(rgb_nuclei)

    # plot polarity vector
    for index, row in single_cell_props.iterrows():
        _add_single_cell_polarity_vector(ax, row["cell_X"], row["cell_Y"], row["nuc_X"], row["nuc_Y"])
        if parameters.show_polarity_angles:
            ax.text(
                row["nuc_Y"], row["nuc_X"], str(int(np.round(row["nuc_displacement_orientation_deg"], 0))),
                color="yellow", fontsize=6
            )

    # set title and ax limits
    _add_title(ax, "nucleus displacement orientation", im_junction, parameters.show_graphics_axis)

    # save output & close
    save_current_fig(
        parameters.graphics_output_format,
        output_path, base_filename,
        "_nucleus_displacement_orientation",
        image=nuc_polarity_angle
    )
    plt.close(fig)


def plot_marker_nucleus_orientation(im_junction, cell_mask, nuclei_mask, single_cell_props, base_filename,
                                    output_path):
    """ function to plot nuclei-organelle polarity vectors

    im_junction :   numpy.array (2-dim), float
                    channel containing the junction staining (used for segmentation)
    cell_mask   :   numpy.array (2-dim), int
                    same dimension as im_junction, contains cell masks
    nuclei_mask :   numpy.array (2-dim), int
                    same dimension as im_junction, contains nuclei masks
    single_cell_props : pandas data frame
    base_filename : string
                    base_filename for plots
    output_path :   string
                    output path for plots
    """
    get_logger().info("Plotting: marker nucleus polarity")

    # figure and axes
    w, h = parameters.graphics_width, parameters.graphics_height
    fig, ax = plt.subplots(figsize=(w, h))

    # resources image
    ax.imshow(im_junction, cmap=plt.cm.gray, alpha=1.0)

    # determine nucleus polarity_angle
    nuc_angle = np.zeros((cell_mask.shape[0], cell_mask.shape[1]))
    for index, row in single_cell_props.iterrows():
        row_label = int(row['label'])
        if row_label == 0:
            continue
        nuc_angle += get_single_cell_mask(row_label, cell_mask) * row['marker_nucleus_orientation_deg']
    nuc_polarity_angle = np.ma.masked_where(cell_mask == 0, nuc_angle)

    # plot polarity angle
    cax = ax.imshow(nuc_polarity_angle, cmap=cm.cm.phase, vmin=0, vmax=360, alpha=0.5)
    color_bar = fig.colorbar(cax, ax=ax, shrink=0.3)  # , extend='both')
    color_bar.set_label("polarity angle")
    color_bar.ax.set_yticks([0, 90, 180, 270, 360])

    # plot nuclei (blue)
    zero = np.zeros((im_junction.shape[0], im_junction.shape[1]))
    rgb_nuclei = np.dstack((zero, zero, nuclei_mask.astype(int) * 256, nuclei_mask.astype(float) * 0.5))
    ax.imshow(rgb_nuclei)

    # plot polarity vector
    for index, row in single_cell_props.iterrows():
        _add_single_cell_polarity_vector(ax, row["nuc_X"], row["nuc_Y"], row["marker_centroid_X"],
                                         row["marker_centroid_Y"])
        if parameters.show_polarity_angles:
            ax.text(
                row["nuc_Y"], row["nuc_X"], str(int(np.round(row["marker_nucleus_orientation_deg"], 0))),
                color="yellow", fontsize=6
            )

    # set title and ax limits
    _add_title(ax, "marker nucleus orientation", im_junction, parameters.show_graphics_axis)

    # save output & close
    save_current_fig(
        parameters.graphics_output_format,
        output_path, base_filename,
        "_marker_nucleus_orientation",
        image=nuc_polarity_angle
    )
    plt.close(fig)


def plot_marker_expression(im_marker, cell_mask, single_cell_dataset, filename, output_path,
                           nuclei_mask=None):
    get_logger().info("Plotting: marker expression")
    # figure and axes
    number_sub_figs = 2  # mean intensity cell, mean intensity membrane
    if nuclei_mask is not None:
        nuclei_mask = nuclei_mask.astype(bool)
        number_sub_figs = 3  # (optional) mean intensity nucleus

    w, h = parameters.graphics_width, parameters.graphics_height
    fig, ax = plt.subplots(1, number_sub_figs, figsize=(w * number_sub_figs, h))

    # plot marker intensity for all subplots
    for i in range(number_sub_figs):
        ax[i].imshow(im_marker, cmap=plt.cm.gray, alpha=1.0)

    outlines_cell, outlines_mem = _get_outline_and_membrane_thickness(im_marker, cell_mask)

    # cell and membrane outline
    outlines_cell_ = np.where(outlines_cell > 0, CELL_OUTLINE_INTENSITY, 0)
    ax[0].imshow(np.ma.masked_where(outlines_cell_ == 0, outlines_cell_), plt.cm.Wistia, vmin=0, vmax=100, alpha=0.5)

    outlines_mem_ = np.where(outlines_mem > 0, CELL_OUTLINE_INTENSITY, 0)
    ax[1].imshow(np.ma.masked_where(outlines_mem_ == 0, outlines_mem_), plt.cm.Wistia, vmin=0, vmax=100, alpha=0.5)

    # nuclei marker intensity
    if nuclei_mask is not None:
        outline_nuc = get_outline_from_mask(nuclei_mask, parameters.outline_width)
        outline_nuc_ = np.where(outline_nuc == True, CELL_OUTLINE_INTENSITY, 0)
        ax[2].imshow(
            np.ma.masked_where(outline_nuc_ == 0, outline_nuc_), plt.cm.Wistia, vmin=0, vmax=100, alpha=0.75
        )  # always last axis

    # plot mean expression value of cell and membrane as text
    for index, row in single_cell_dataset.iterrows():
        ax[0].text(row["cell_Y"], row["cell_X"], str(np.round(row["marker_mean_expr"], 1)), color="w", fontsize=7)
        ax[1].text(row["cell_Y"], row["cell_X"], str(np.round(row["marker_mean_expression_mem"], 1)), color="w",
                   fontsize=7)
        if nuclei_mask is not None:
            ax[2].text(
                row["nuc_Y"], row["nuc_X"], str(np.round(row["marker_mean_expression_nuc"], 1)), color="w", fontsize=7
            )

    # set title
    ax[0].set_title("mean intensity cell")
    ax[1].set_title("mean intensity membrane")
    if nuclei_mask is not None:
        ax[2].set_title("mean intensity nucleus")

    if not parameters.show_graphics_axis:
        for ax_ in ax:
            ax_.axis('off')

    # save output & close
    save_current_fig(parameters.graphics_output_format, output_path, filename, "_marker_expression")
    plt.close(fig)


def plot_marker_polarity(im_marker, cell_mask, single_cell_props, filename, output_path):
    get_logger().info("Plotting: marker polarity")
    # figure and axes
    w, h = parameters.graphics_width, parameters.graphics_height
    fig, ax = plt.subplots(1, figsize=(w, h))

    # plot marker intensity
    ax.imshow(im_marker, cmap=plt.cm.gray, alpha=1.0)

    # cumulative cell outlines
    outlines_cell_accumulated = np.zeros((im_marker.shape[0], im_marker.shape[1]))
    for cell_label in np.unique(cell_mask):
        # exclude background
        if cell_label == 0:
            continue

        single_cell_mask = get_single_cell_mask(cell_label, cell_mask)
        outline_cell = get_outline_from_mask(single_cell_mask, parameters.outline_width)
        outline_cell_ = np.where(outline_cell == True, 30, 0)
        outlines_cell_accumulated += outline_cell_

    # plot non-cumulative cell outlines
    outlines_cell_ = np.where(outlines_cell_accumulated > 0, CELL_OUTLINE_INTENSITY, 0)
    ax.imshow(np.ma.masked_where(outlines_cell_ == 0, outlines_cell_), plt.cm.Wistia, vmin=0, vmax=100, alpha=0.5)

    # add all polarity vectors
    for index, row in single_cell_props.iterrows():
        _add_single_cell_polarity_vector(ax, row["cell_X"], row["cell_Y"], row["marker_centroid_X"],
                                         row["marker_centroid_Y"])

    ax.set_title("marker polarity")
    if not parameters.show_graphics_axis:
        ax.axis('off')

    # save output & close
    save_current_fig(parameters.graphics_output_format, output_path, filename, "_marker_polarity")
    plt.close(fig)


def plot_corners(im_junction, single_cell_props, filename, output_path):
    w, h = parameters.graphics_width, parameters.graphics_height
    fig, ax = plt.subplots(1, figsize=(w, h))

    # plot marker intensity
    ax.imshow(im_junction, cmap=plt.cm.gray, alpha=1.0)

    for index, row in single_cell_props.iterrows():
        plt.scatter(np.array(json.loads(row["cell_corner_points"]))[:, 0],
                    np.array(json.loads(row["cell_corner_points"]))[:, 1],
                    [4]*len(np.array(json.loads(row["cell_corner_points"]))[:, 1]))

    save_current_fig(parameters.graphics_output_format, output_path, filename, "_cell_corners")


def plot_junction_polarity(im_junction, cell_mask, single_cell_props, filename, output_path):
    get_logger().info("Plotting: junction polarity")
    # figure and axes
    w, h = parameters.graphics_width, parameters.graphics_height
    fig, ax = plt.subplots(1, figsize=(w, h))

    # plot marker intensity
    ax.imshow(im_junction, cmap=plt.cm.gray, alpha=1.0)

    # cumulative cell outlines
    outlines_cell_accumulated = np.zeros((im_junction.shape[0], im_junction.shape[1]))
    for cell_label in np.unique(cell_mask):
        # exclude background
        if cell_label == 0:
            continue

        single_cell_mask = get_single_cell_mask(cell_label, cell_mask)
        outline_cell = get_outline_from_mask(single_cell_mask, parameters.outline_width)
        outline_cell_ = np.where(outline_cell == True, 30, 0)
        outlines_cell_accumulated += outline_cell_

    # plot non-cumulative cell outlines
    outlines_cell_ = np.where(outlines_cell_accumulated > 0, CELL_OUTLINE_INTENSITY, 0)
    ax.imshow(np.ma.masked_where(outlines_cell_ == 0, outlines_cell_), plt.cm.Wistia, vmin=0, vmax=100, alpha=0.5)

    # add all polarity vectors
    for index, row in single_cell_props.iterrows():
        _add_single_cell_polarity_vector(ax, row["cell_X"], row["cell_Y"], row["junction_centroid_X"],
                                         row["junction_centroid_Y"])

    ax.set_title("junction polarity")

    if not parameters.show_graphics_axis:
        ax.axis('off')

    # save output & close
    save_current_fig(parameters.graphics_output_format, output_path, filename, "_junction_polarity")
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


def _add_colorbar(fig, cax, ax, yticks, label):
    color_bar = fig.colorbar(cax, ax=ax, shrink=0.3)
    color_bar.set_label(label)
    color_bar.ax.set_yticks(yticks)


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
    _add_colorbar(fig, cax_0, ax, yticks, "eccentricity")


def _add_nuclei_eccentricity(fig, ax, im_junction, nuclei_mask, nuclei_eccentricity):
    v_min = 0.0
    v_max = 1.0
    yticks = [0.0, 0.5, 1.0]

    ax.imshow(im_junction, cmap=plt.cm.gray, alpha=1.0)

    # show nuclei eccentricity everywhere but background label
    nuclei_mask_ = np.where(nuclei_mask == True, 1, 0)
    cax_1 = ax.imshow(np.ma.masked_where(
        nuclei_mask_ == 0, nuclei_eccentricity), cmap=plt.cm.bwr, vmin=v_min, vmax=v_max, alpha=0.5
    )

    # colorbar
    _add_colorbar(fig, cax_1, ax, yticks, "eccentricity")


def _add_title(ax, plot_title, im_junction, axis_on):
    ax.set_title(plot_title)
    ax.set_xlim(0, im_junction.shape[1])
    ax.set_ylim(0, im_junction.shape[0])
    ax.invert_yaxis()
    if not axis_on:
        ax.axis('off')


def _calc_nuc_eccentricity(single_cell_props, cell_mask, nuclei_mask):
    nuclei_eccentricity = np.zeros((cell_mask.shape[0], cell_mask.shape[1]))
    for index, row in single_cell_props.iterrows():
        row_label = int(row['label'])

        # exclude background
        if row_label == 0:
            continue
        single_cell_mask = np.where(cell_mask == row_label, True, 0)
        single_nuclei_mask_ = np.logical_and(single_cell_mask, nuclei_mask)
        nuclei_eccentricity += np.where(single_nuclei_mask_ == True, 1, 0) * row['nuc_eccentricity']

    get_logger().info("Maximal nuclei eccentricity: %s" % str(np.max(nuclei_eccentricity)))
    get_logger().info("Minimal nuclei eccentricity: %s" % str(np.min(nuclei_eccentricity)))

    return nuclei_eccentricity


def _calc_cell_eccentricity(single_cell_props, cell_mask):
    cell_eccentricity = np.zeros((cell_mask.shape[0], cell_mask.shape[1]))
    for index, row in single_cell_props.iterrows():
        row_label = int(row['label'])
        if row_label == 0:
            continue
        cell_eccentricity += get_single_cell_mask(row_label, cell_mask) * row['cell_eccentricity']

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


def plot_eccentricity(im_junction, single_cell_props, cell_mask, filename, output_path, nuclei_mask=None):
    """ function to plot cell (and optionally nuclei) eccentricity

    im_junction :   numpy.array (2-dim), float
                    channel containing the junction staining (used for segmentation)
    single_cell_props       :   pandas properties dataset
    filename    :   resources filename for the output file
    output_path :   desired output path
    cell_mask   :   cellpose cell mask
    nuclei_mask :   cellpose nuclei mask
    """
    get_logger().info("Plotting: eccentricity")

    # figure and axes
    number_sub_figs = 1
    if nuclei_mask is not None:
        nuclei_mask = nuclei_mask.astype(bool)
        number_sub_figs = 2

    w, h = parameters.graphics_width, parameters.graphics_height
    fig, ax = plt.subplots(1, number_sub_figs, figsize=(number_sub_figs * w, h))

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
                row['cell_Y'],
                row['cell_X'],
                row['cell_shape_orientation'],
                row['cell_major_axis_length'],
                row['cell_minor_axis_length'],
                row["cell_eccentricity"]
            )

            # plot orientation degree nucleus
            _add_single_cell_eccentricity_axis(
                ax[1],
                row['nuc_Y'],
                row['nuc_X'],
                row['nuc_shape_orientation'],
                row['nuc_major_axis_length'],
                row['nuc_minor_axis_length'],
                row["nuc_eccentricity"]
            )
        else:
            _add_single_cell_eccentricity_axis(
                ax,
                row['cell_Y'],
                row['cell_X'],
                row['cell_shape_orientation'],
                row['cell_major_axis_length'],
                row['cell_minor_axis_length'],
                row["cell_eccentricity"]
            )

    # set title and ax limits
    if nuclei_mask is not None:
        _add_title(ax[0], "cell elongation", im_junction, parameters.show_graphics_axis)
        _add_title(ax[1], "nuclei elongation", im_junction, parameters.show_graphics_axis)
    else:
        _add_title(ax, "cell elongation", im_junction, parameters.show_graphics_axis)

    # save output & close
    save_current_fig(parameters.graphics_output_format, output_path, filename, "_eccentricity")
    plt.close(fig)


def _calc_nuc_orientation(single_cell_props, cell_mask, nuclei_mask):
    get_logger().info("Calculating nuclei orientation...")
    nuclei_orientation = np.zeros((cell_mask.shape[0], cell_mask.shape[1]))
    for index, row in single_cell_props.iterrows():
        row_label = int(row['label'])

        # exclude background
        if row_label == 0:
            continue
        single_cell_mask = np.where(cell_mask == row_label, True, 0)
        single_nuclei_mask_ = np.logical_and(single_cell_mask, nuclei_mask)
        nuclei_orientation += np.where(single_nuclei_mask_ == True, 1, 0) * row['nuc_shape_orientation'] * 180.0 / np.pi

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
            'nuc_shape_orientation'] * 180.0 / np.pi

    get_logger().info("Maximal cell orientation: %s" % str(np.max(cell_orientation)))
    get_logger().info("Minimal cell orientation: %s" % str(np.min(cell_orientation)))

    return cell_orientation


def _add_cell_orientation(fig, ax, im_junction, cell_mask, cell_orientation):
    v_min = 0.0
    v_max = 180.0
    yticks = [0.0, 45.0, 90.0, 135.0, 180.0]

    ax.imshow(im_junction, cmap=plt.cm.gray, alpha=1.0)

    # show cell_orientation everywhere but background label
    cax = ax.imshow(
        np.ma.masked_where(cell_mask == 0, cell_orientation), cmap=cm.cm.phase, vmin=v_min, vmax=v_max,
        alpha=0.75
    )

    # colorbar
    _add_colorbar(fig, cax, ax, yticks, "shape orientation (degree)")


def _add_nuclei_orientation(fig, ax, im_junction, nuclei_mask, nuclei_orientation):
    v_min = 0.0
    v_max = 180.0
    yticks = [0.0, 45.0, 90.0, 135.0, 180.0]

    ax.imshow(im_junction, cmap=plt.cm.gray, alpha=1.0)

    # show nuclei orientation everywhere but background label
    nuclei_mask_ = np.where(nuclei_mask == True, 1, 0)
    cax_1 = ax.imshow(np.ma.masked_where(
        nuclei_mask_ == 0, nuclei_orientation), cmap=cm.cm.phase, vmin=v_min, vmax=v_max, alpha=0.75
    )

    # colorbar
    _add_colorbar(fig, cax_1, ax, yticks, "shape orientation (degree)")


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


def plot_orientation(im_junction, single_cell_props, filename, output_path, cell_mask, nuclei_mask=None):
    """ function to plot cell (and optionally nuclei) orientation

    im_junction :   numpy.array (2-dim), float
                    channel containing the junction staining (used for segmentation)
    single_cell_props       :   pandas properties dataset
    filename    :   resources filename for the output file
    output_path :   desired output path
    cell_mask   :   cellpose cell mask
    nuclei_mask :   cellpose nuclei mask
    """
    get_logger().info("Plotting: orientation")

    # figure and axes
    number_sub_figs = 1
    if nuclei_mask is not None:
        nuclei_mask = nuclei_mask.astype(bool)
        number_sub_figs = 2

    w, h = parameters.graphics_width, parameters.graphics_height
    fig, ax = plt.subplots(1, number_sub_figs, figsize=(number_sub_figs * w, h))

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
                row['cell_Y'],
                row['cell_X'],
                row['cell_shape_orientation'],
                row['cell_major_axis_length'],
                row['cell_minor_axis_length']
            )

            # plot orientation degree nucleus
            _add_single_cell_orientation_degree_axis(
                ax[1],
                row['nuc_Y'],
                row['nuc_X'],
                row['nuc_shape_orientation'],
                row['nuc_major_axis_length'],
                row['nuc_minor_axis_length']
            )
        else:
            # plot orientation degree
            _add_single_cell_orientation_degree_axis(
                ax,
                row['cell_Y'],
                row['cell_X'],
                row['cell_shape_orientation'],
                row['cell_major_axis_length'],
                row['cell_minor_axis_length']
            )

    # set title and ax limits
    if nuclei_mask is not None:
        _add_title(ax[0], "cell shape orientation", im_junction, parameters.show_graphics_axis)
        _add_title(ax[1], "nuclei shape orientation", im_junction, parameters.show_graphics_axis)
    else:
        _add_title(ax, "cell shape orientation", im_junction, parameters.show_graphics_axis)

    # save output & close
    save_current_fig(parameters.graphics_output_format, output_path, filename, "_shape_orientation")
    plt.close(fig)


def plot_ratio_method(im_junction, cell_mask, single_cell_props, filename, output_path):
    get_logger().info("Plotting: ratio method")

    # figure and axes
    w, h = parameters.graphics_width, parameters.graphics_height
    fig, ax = plt.subplots(1, 1, figsize=(w, h))

    # show junction and cell mask overlay
    ax.imshow(im_junction, cmap=plt.cm.gray, alpha=1.0)
    ax.imshow(cell_mask, cmap=plt.cm.Set3, alpha=0.25)

    cell_outlines_accumulated = np.zeros((im_junction.shape[0], im_junction.shape[1]))
    for cell_label in np.unique(cell_mask):
        # exclude background
        if cell_label == 0:
            continue

        single_cell_mask = get_single_cell_mask(cell_label, cell_mask)
        cell_outline = get_outline_from_mask(single_cell_mask, parameters.membrane_thickness)
        # accumulates cell outlines. overlapping outlines have a higher value
        cell_outlines_accumulated += np.where(cell_outline == True, 1, 0)

    # overlapping accumulated outlines are ignored and set to 1.
    cell_outlines = np.where(cell_outlines_accumulated > 0, CELL_OUTLINE_INTENSITY, 0)
    ax.imshow(np.ma.masked_where(cell_outlines == 0, cell_outlines), plt.cm.bwr, vmin=0, vmax=100, alpha=0.5)

    # plot major axis around coordinates of each cell
    for index, row in single_cell_props.iterrows():
        x0 = row['cell_X']
        y0 = row['cell_Y']

        # upper
        x1 = x0 + math.sin(np.pi / 4.0) * 0.5 * row['cell_major_axis_length']
        y1 = y0 + math.cos(np.pi / 4.0) * 0.5 * row['cell_major_axis_length']
        x2 = x0 + math.cos(np.pi / 4.0) * 0.5 * row['cell_major_axis_length']
        y2 = y0 - math.sin(np.pi / 4.0) * 0.5 * row['cell_major_axis_length']

        ax.plot((y0, y1), (x0, x1), '--r', linewidth=0.5)
        ax.plot((y0, y2), (x0, x2), '--r', linewidth=0.5)
        ax.plot(y0, x0, '.b', markersize=5)

        # lower
        x1 = x0 - math.sin(np.pi / 4.0) * 0.5 * row['cell_major_axis_length']
        y1 = y0 - math.cos(np.pi / 4.0) * 0.5 * row['cell_major_axis_length']
        x2 = x0 - math.cos(np.pi / 4.0) * 0.5 * row['cell_major_axis_length']
        y2 = y0 + math.sin(np.pi / 4.0) * 0.5 * row['cell_major_axis_length']

        ax.plot((y0, y1), (x0, x1), '--r', linewidth=0.5)
        ax.plot((y0, y2), (x0, x2), '--r', linewidth=0.5)
        ax.plot(y0, x0, '.b', markersize=5)

    ax.set_title("ratio method")
    ax.set_xlim(0, im_junction.shape[0])
    ax.set_ylim(0, im_junction.shape[1])

    if not parameters.show_graphics_axis:
        ax.axis('off')

    # save output & close
    save_current_fig(parameters.graphics_output_format, output_path, filename, "_ratio_method")
    plt.close(fig)


# def plot_adjacency_matrix(label_image, intensity_image):
#    # todo: needed?
#    rag = orientation_graph_nf(label_image)
#    out = graph.draw_rag(label_image, rag, intensity_image, node_color="#ffde00")
#    return out


def plot_dataset(properties_ds, cell_mask, nuclei_mask, organelle_mask, img_marker, img_junction, filename,
                 output_path):
    """Plots the properties dataset"""
    get_logger().info("Plotting...")

    # TODO: adapt name plot polarity in parameter files
    if parameters.plot_polarity and nuclei_mask is not None and organelle_mask is not None:
        plot_organelle_polarity(
            img_junction,
            cell_mask,
            nuclei_mask,
            organelle_mask,
            properties_ds,
            filename,
            output_path
        )
        if nuclei_mask is not None:
            plot_nuc_displacement_orientation(
                img_junction,
                cell_mask,
                nuclei_mask,
                properties_ds,
                filename,
                output_path
            )
    if parameters.plot_marker and img_marker is not None:
        plot_marker_expression(
            img_marker,
            cell_mask,
            properties_ds,
            filename,
            output_path,
            nuclei_mask=nuclei_mask
        )
        plot_marker_polarity(
            img_marker,
            cell_mask,
            properties_ds,
            filename,
            output_path
        )
        if nuclei_mask is not None:
            plot_marker_nucleus_orientation(
                img_junction,
                cell_mask,
                nuclei_mask,
                properties_ds,
                filename,
                output_path
            )
    if parameters.plot_junctions and img_junction is not None:
        plot_junction_polarity(img_junction, cell_mask, properties_ds, filename, output_path)
        plot_corners(img_junction, properties_ds, filename, output_path)

    if parameters.plot_orientation:
        plot_eccentricity(img_junction, properties_ds, cell_mask, filename, output_path,
                          nuclei_mask=nuclei_mask)
    if parameters.plot_ratio_method:
        plot_ratio_method(
            img_junction,
            cell_mask,
            properties_ds,
            filename,
            output_path
        )
    if parameters.plot_cyclic_orientation:
        plot_orientation(
            img_junction,
            properties_ds,
            filename,
            output_path,
            cell_mask,
            nuclei_mask=nuclei_mask
        )
