import math
from pathlib import Path

import cmocean as cm
import matplotlib as mpl
import numpy as np
from matplotlib import pyplot as plt
from skimage.future import graph
from skimage.measure import label, regionprops

from vascu_ec.utils.rag import orientation_graph_nf
from vascu_ec.utils.seg import get_outline_from_mask
from vascu_ec.vascu_ec_logging import get_logger

# for figure plot resolution
FIGURE_DPI = 300


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


def plot_polarity(parameters, im_junction, cell_mask, nuclei_mask, golgi_mask, single_cell_props, filename,
                  output_path):
    """
    parameters  :   dict
                    user defined parameters
    im_junction :   numpy.array (2-dim), float
                    channel containing the junction staining (used for segmentation)
    masks       :   numpy.array (2-dim), int
                    same dimension as im_junction, contains cell masks
    """

    width = parameters["graphics_width"]
    fig, ax = plt.subplots(figsize=(width, width))

    nuclei_mask = nuclei_mask.astype(bool)
    golgi_mask = golgi_mask.astype(bool)

    nuclei_mask_ = np.where(nuclei_mask == True, 70, 0)  # todo: why 70? why doing that at all if it is a mask?
    golgi_mask_ = np.where(golgi_mask == True, 1, 0)

    ax.imshow(im_junction, cmap=plt.cm.gray, alpha=1.0)
    ax.imshow(cell_mask, cmap=plt.cm.Set3, alpha=0.25)
    ax.imshow(np.ma.masked_where(nuclei_mask_ == 0, nuclei_mask_), plt.cm.gist_rainbow, vmin=0, vmax=100, alpha=0.5)
    ax.imshow(np.ma.masked_where(golgi_mask_ == 0, golgi_mask_), plt.cm.gist_rainbow, vmin=0, vmax=100, alpha=0.5)

    for index, row in single_cell_props.iterrows():
        ax.plot(row["Y_nuc"], row["X_nuc"], '.g', markersize=1)
        ax.plot(row["Y_golgi"], row["X_golgi"], '.m', markersize=1)
        ax.arrow(row["Y_nuc"], row["X_nuc"], row["Y_golgi"] - row["Y_nuc"], row["X_golgi"] - row["X_nuc"],
                 color='white', width=2)
        if parameters["show_polarity_angles"]:
            ax.text(row["Y_cell"], row["X_cell"], str(int(np.round(row["angle_deg"], 0))), color="yellow", fontsize=6)

    ax.set_xlim(0, im_junction.shape[0])
    ax.set_ylim(0, im_junction.shape[1])

    if "pdf" in parameters["graphics_output_format"]:
        plt.savefig(str(Path(output_path).joinpath(filename + "_nuclei_golgi_vector.pdf")))
    if "svg" in parameters["graphics_output_format"]:
        plt.savefig(str(Path(output_path).joinpath(filename + "_nuclei_golgi_vector.svg")))
    if "png" in parameters["graphics_output_format"]:
        plt.savefig(str(Path(output_path).joinpath(filename + "_nuclei_golgi_vector.png")))
    plt.close(fig)


def plot_marker_expression(parameters, im_marker, cell_mask, single_cell_dataset, filename, output_path,
                           nuclei_mask=None):
    number_sub_figs = 1
    if nuclei_mask is not None:
        nuclei_mask = nuclei_mask.astype(bool)
        number_sub_figs = 2

    fig, ax = plt.subplots(1, number_sub_figs, figsize=(10 * number_sub_figs, 10))

    if nuclei_mask is not None:
        outline_nuc = get_outline_from_mask(nuclei_mask, parameters["outline_width"])
        outline_nuc_ = np.where(outline_nuc == True, 30, 0)  # todo: 30 ?
        ax[number_sub_figs - 1].imshow(
            np.ma.masked_where(outline_nuc_ == 0, outline_nuc_), plt.cm.Wistia, vmin=0, vmax=100, alpha=0.75
        )

    for i in range(number_sub_figs):
        ax[i].imshow(im_marker, cmap=plt.cm.gray, alpha=1.0)

    outlines_cell = np.zeros((im_marker.shape[0], im_marker.shape[1]))
    outlines_mem = np.zeros((im_marker.shape[0], im_marker.shape[1]))

    for cell_label in np.unique(cell_mask):
        # exclude background
        if cell_label == 0:
            continue

        single_cell_mask = np.where(cell_mask == cell_label, 1, 0)
        outline_cell = get_outline_from_mask(single_cell_mask, parameters["outline_width"])
        outline_cell_ = np.where(outline_cell == True, 30, 0)
        outlines_cell += outline_cell_

        outline_mem = get_outline_from_mask(single_cell_mask, parameters["membrane_thickness"])
        outline_mem_ = np.where(outline_mem == True, 30, 0)
        outlines_mem += outline_mem_

    outlines_cell_ = np.where(outlines_cell > 0, 30, 0)  # todo: what is that threshold?
    ax[0].imshow(np.ma.masked_where(outlines_cell_ == 0, outlines_cell_), plt.cm.Wistia, vmin=0, vmax=100, alpha=0.5)

    outlines_mem_ = np.where(outlines_mem > 0, 30, 0)  # todo: what is that threshold?
    ax[1].imshow(np.ma.masked_where(outlines_mem_ == 0, outlines_mem_), plt.cm.Wistia, vmin=0, vmax=100, alpha=0.5)

    for index, row in single_cell_dataset.iterrows():
        ax[0].text(row["Y_cell"], row["X_cell"], str(np.round(row["mean_expression"], 1)), color="w", fontsize=6)
        ax[1].text(row["Y_cell"], row["X_cell"], str(np.round(row["mean_expression_mem"], 1)), color="w", fontsize=6)
        if nuclei_mask is not None:
            ax[number_sub_figs - 1].text(row["Y_nuc"], row["X_nuc"], str(np.round(row["mean_expression_nuc"], 1)),
                                         color="w",
                                         fontsize=6)

    ax[0].set_title("mean intensity cell")
    ax[1].set_title("mean intensity membrane")
    if nuclei_mask is not None:
        ax[number_sub_figs - 1].set_title("mean intensity nucleus")

    # save output
    plt.savefig(str(Path(output_path).joinpath(filename + "_marker_expression.pdf")))
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

    plt.savefig(str(Path(output_path).joinpath(filename + "_marker_polarity.pdf")))
    plt.savefig(str(Path(output_path).joinpath(filename + "_marker_polarity.png")))
    plt.close(fig)


def plot_alignment(im_junction, single_cell_props, filename, output_path, cell_mask, nuclei_mask=None):
    number_sub_figs = 1

    if nuclei_mask is not None:
        nuclei_mask = nuclei_mask.astype(bool)
        number_sub_figs = 2

    fig, ax = plt.subplots(1, number_sub_figs, figsize=(number_sub_figs * 5, 5))

    cell_eccentricity = np.zeros((cell_mask.shape[0], cell_mask.shape[1]))

    for index, row in single_cell_props.iterrows():
        row_label = int(row['label'])
        if row_label == 0:
            continue
        single_cell_mask = np.where(cell_mask == row_label, 1, 0) * row['eccentricity']
        cell_eccentricity += single_cell_mask

    if nuclei_mask is not None:
        ax[0].imshow(im_junction, cmap=plt.cm.gray, alpha=1.0)

        cax_0 = ax[0].imshow(
            np.ma.masked_where(cell_mask == 0, cell_eccentricity), cmap=plt.cm.bwr, vmin=0, vmax=1.0, alpha=0.5
        )
        color_bar = fig.colorbar(cax_0, ax=ax[0], shrink=0.3)
        color_bar.set_label('eccentricity')
        ax[1].imshow(im_junction, cmap=plt.cm.gray, alpha=1.0)

        get_logger().info("Calculating nuclei eccentricity...")
        nuclei_eccentricity = np.zeros((cell_mask.shape[0], cell_mask.shape[1]))
        for index, row in single_cell_props.iterrows():
            row_label = int(row['label'])

            # exclude background
            if row_label == 0:
                continue

            single_cell_mask = np.where(cell_mask == row_label, True, 0) * row['eccentricity']
            single_nuclei_mask_ = np.logical_and(single_cell_mask, nuclei_mask)
            single_nuclei_mask = np.where(single_nuclei_mask_ == True, 1, 0) * row['eccentricity_nuc']
            nuclei_eccentricity += single_nuclei_mask

        get_logger().info("Maximal nuclei eccentricity: %s" % str(np.max(nuclei_eccentricity)))
        get_logger().info("Minimal nuclei eccentricity: %s" % str(np.min(nuclei_eccentricity)))

        nuclei_mask_ = np.where(nuclei_mask is True, 70, 0)  # todo: threshold?

        cax_1 = ax[1].imshow(np.ma.masked_where(
            nuclei_mask_ == 0, nuclei_eccentricity), plt.cm.bwr, vmin=0, vmax=1.0, alpha=0.5
        )

        color_bar = fig.colorbar(cax_1, ax=ax[1], shrink=0.3)  # , extend='both')
        color_bar.set_label('eccentricity')

    else:
        ax.imshow(im_junction, cmap=plt.cm.gray, alpha=1.0)
        ax.imshow(np.ma.masked_where(cell_mask == 0, cell_eccentricity), cmap=plt.cm.bwr, alpha=0.25)

    for index, row in single_cell_props.iterrows():
        orientation = row['shape_orientation']
        x0 = row['X_cell']
        y0 = row['Y_cell']

        x1_major = x0 + math.sin(orientation) * 0.5 * row['major_axis_length']
        y1_major = y0 + math.cos(orientation) * 0.5 * row['major_axis_length']
        x2_major = x0 - math.sin(orientation) * 0.5 * row['major_axis_length']
        y2_major = y0 - math.cos(orientation) * 0.5 * row['major_axis_length']

        x1_minor = x0 + math.cos(orientation) * 0.5 * row['minor_axis_length']
        y1_minor = y0 - math.sin(orientation) * 0.5 * row['minor_axis_length']
        x2_minor = x0 - math.cos(orientation) * 0.5 * row['minor_axis_length']
        y2_minor = y0 + math.sin(orientation) * 0.5 * row['minor_axis_length']

        if nuclei_mask is not None:
            # plot orientation degree
            ax[0].plot((y1_major, y2_major), (x1_major, x2_major), '--w', linewidth=0.5)
            ax[0].plot((y1_minor, y2_minor), (x1_minor, x2_minor), '--w', linewidth=0.5)
            ax[0].plot(y0, x0, '.b', markersize=5)
            orientation_degree = 180.0 * orientation / np.pi
            ax[0].text(y0, x0, str(int(np.round(orientation_degree, 0))), color="yellow", fontsize=4)

            # plot orientation degree nucleus
            orientation = row['shape_orientation_nuc']
            x0 = row['X_nuc']
            y0 = row['Y_nuc']

            x1_major = x0 + math.sin(orientation) * 0.5 * row['major_axis_length_nuc']
            y1_major = y0 + math.cos(orientation) * 0.5 * row['major_axis_length_nuc']
            x2_major = x0 - math.sin(orientation) * 0.5 * row['major_axis_length_nuc']
            y2_major = y0 - math.cos(orientation) * 0.5 * row['major_axis_length_nuc']

            x1_minor = x0 + math.cos(orientation) * 0.5 * row['minor_axis_length_nuc']
            y1_minor = y0 - math.sin(orientation) * 0.5 * row['minor_axis_length_nuc']
            x2_minor = x0 - math.cos(orientation) * 0.5 * row['minor_axis_length_nuc']
            y2_minor = y0 + math.sin(orientation) * 0.5 * row['minor_axis_length_nuc']

            ax[1].plot((y1_major, y2_major), (x1_major, x2_major), '--w', linewidth=0.5)
            ax[1].plot((y1_minor, y2_minor), (x1_minor, x2_minor), '--w', linewidth=0.5)

            orientation_degree = 180.0 * orientation / np.pi
            ax[1].text(y0, x0, str(int(np.round(orientation_degree, 0))), color="yellow", fontsize=4)
        else:
            ax.plot((y1_major, y2_major), (x1_major, x2_major), '--w', linewidth=0.5)
            ax.plot((y1_minor, y2_minor), (x1_minor, x2_minor), '--w', linewidth=0.5)
            ax.plot(y0, x0, '.b', markersize=5)
            orientation_degree = 180.0 * orientation / np.pi
            ax.text(y0, x0, str(int(np.round(orientation_degree, 0))), color="yellow", fontsize=4)

    # set title and ax limits
    if nuclei_mask is not None:
        ax[0].set_title("alignment cell shape")
        ax[0].set_xlim(0, im_junction.shape[0])
        ax[0].set_ylim(0, im_junction.shape[1])

        ax[1].set_title("alignment nuclei shape")
        ax[1].set_xlim(0, im_junction.shape[0])
        ax[1].set_ylim(0, im_junction.shape[1])
    else:
        ax.set_title("alignment cell shape")
        ax.set_xlim(0, im_junction.shape[0])
        ax.set_ylim(0, im_junction.shape[1])

    # save to disk
    plt.tight_layout()
    plt.savefig(str(Path(output_path).joinpath(filename + "_alignment.pdf")))
    plt.savefig(str(Path(output_path).joinpath(filename + "_alignment.svg")))
    plt.savefig(str(Path(output_path).joinpath(filename + "_alignment.png")))
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
    plt.savefig(str(Path(output_path).joinpath(filename + "_ratio_method.pdf")))
    plt.savefig(str(Path(output_path).joinpath(filename + "_ratio_method.png")))
    plt.close(fig)


def plot_cyclical(input_image, cell_masks_approx):
    lab_image = label(input_image)
    regions = regionprops(lab_image)
    orient_list = []
    for region in regions:
        cell_masks_approx[lab_image == region.label] = region.orientation * (180 / np.pi) + 90
        orient_list.append(np.round(math.degrees(region.orientation)))

    cell_masks_approx[cell_masks_approx == 0] = np.nan

    masked_data = (lab_image != 0)
    t_f = (masked_data * 1)
    masked_data[np.where(t_f == 0), np.where(t_f == 0)] = np.nan

    fig = plt.figure(figsize=(20, 20))
    cm_phase = plt.imshow(cell_masks_approx, cmap=cm.cm.phase)  # todo: color scale OK? Put on top of file.
    plt.imshow(np.where(masked_data == 0, 1, np.nan), cmap='binary', vmin=0, vmax=1)
    plt.colorbar(cm_phase)
    plt.close(fig)

    return cm_phase


def plot_adjacency_matrix(label_image, intensity_image):
    # todo: needed?
    rag = orientation_graph_nf(label_image)
    out = graph.draw_rag(label_image, rag, intensity_image, node_color="#ffde00")
    return out


def plot_dataset(parameters, img, properties_dataset, output_path, filename, cell_mask_rem_island, nuclei_mask,
                 golgi_mask, im_marker):
    """Plots the properties dataset"""
    # todo: close figures
    get_logger().info("Plotting data...")
    im_junction = img[:, :, int(parameters["channel_junction"])]

    if parameters["plot_polarity"] and nuclei_mask is not None and golgi_mask is not None:
        plot_polarity(
            parameters,
            im_junction,
            cell_mask_rem_island,
            nuclei_mask,
            golgi_mask,
            properties_dataset,
            filename, output_path,
        )
    if parameters["plot_marker"] and nuclei_mask is not None and im_marker is not None:
        plot_marker_expression(
            parameters, im_marker, cell_mask_rem_island, properties_dataset, filename, output_path,
            nuclei_mask=nuclei_mask
        )
        plot_marker_polarity(
            parameters, im_marker, cell_mask_rem_island, properties_dataset, filename, output_path
        )
    if parameters["plot_marker"] and nuclei_mask is None and im_marker is not None:
        plot_marker_expression(parameters, im_marker, cell_mask_rem_island, properties_dataset, filename, output_path)
        plot_marker_polarity(
            parameters, im_marker, cell_mask_rem_island, properties_dataset, filename, output_path
        )

    if parameters["plot_alignment"] and nuclei_mask is not None:
        plot_alignment(
            im_junction,
            properties_dataset,
            filename,
            output_path,
            cell_mask_rem_island,
            nuclei_mask=nuclei_mask,
        )
    if parameters["plot_marker"] and nuclei_mask is None:
        plot_alignment(
            im_junction,
            properties_dataset,
            filename,
            output_path,
            cell_mask_rem_island,
        )
    if parameters["plot_ratio_method"]:
        plot_ratio_method(
            parameters,
            im_junction,
            [cell_mask_rem_island],
            properties_dataset,
            filename,
            output_path
        )
