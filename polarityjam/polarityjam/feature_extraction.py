import numpy as np
import pandas as pd
import skimage.filters
import skimage.io
import skimage.measure
import skimage.segmentation

from polarityjam.polarityjam_logging import get_logger
from polarityjam.utils.masks import get_single_cell_membrane_mask, get_single_cell_mask, get_single_cell_organelle_mask, \
    get_single_cell_cytosol_mask, get_single_cell_nucleus_mask, get_organelle_mask, get_nuclei_mask
from polarityjam.utils.moran import Moran
from polarityjam.utils.plot import plot_dataset
from polarityjam.utils.rag import orientation_graph_nf
from polarityjam.utils.weights import W


def get_image_for_segmentation(parameters, img):
    """Returns segmentation (junction, nucleus) if multichannel, else junction only."""

    ch_junction = int(parameters["channel_junction"])
    ch_nucleus = int(parameters["channel_nucleus"])

    get_logger().info("Image shape: %s" % str(img.shape))
    get_logger().info("Junction channel at position: %s" % str(ch_junction))

    # add debug out output here
    im_junction = img[:, :, ch_junction]
    if ch_nucleus >= 0:
        im_nucleus = img[:, :, ch_nucleus]
        im_seg = np.array([im_junction, im_nucleus])
    else:
        return im_junction

    return im_seg


def threshold(parameters, single_cell_mask, single_nucleus_mask=None, single_organelle_mask=None):
    """Thresholds given single_cell_mask. Returns True if falls under threshold."""
    # TODO: check if this can be removed, we already remove small objects from the cellpose mask
    if len(single_cell_mask[single_cell_mask == 1]) < parameters["min_cell_size"]:
        return True

    if single_nucleus_mask is not None:
        if len(single_nucleus_mask[single_nucleus_mask == 1]) < parameters["min_nucleus_size"]:
            return True

    if single_organelle_mask is not None:
        if len(single_organelle_mask[single_organelle_mask == 1]) < parameters["min_organelle_size"]:
            return True

    return False


def get_image_marker(parameters, img):
    """Gets the image marker."""
    if parameters["channel_expression_marker"] >= 0:
        return img[:, :, parameters["channel_expression_marker"]]


def set_single_cell_props(properties_dataset, index, filename, connected_component_label, single_cell_mask, graph_nf):
    single_cell_props = get_single_cell_prop(single_cell_mask)
    neighbours = len(list(graph_nf.neighbors(connected_component_label)))
    fill_single_cell_general_data_frame(
        properties_dataset, index, filename, connected_component_label, single_cell_props, neighbours
    )

    return single_cell_props


def set_single_cell_nucleus_props(properties_dataset, index, single_nucleus_mask):
    regions = skimage.measure.regionprops(single_nucleus_mask)
    nucleus_props = None
    if regions:
        nucleus_props = regions[-1]
        fill_single_nucleus_data_frame(properties_dataset, index, nucleus_props)

    return nucleus_props


def set_single_cell_organelle_props(properties_dataset, index, single_organelle_mask, nucleus_props):
    regions = skimage.measure.regionprops(single_organelle_mask)
    organelle_props = None
    if regions:
        organelle_props = regions[-1]
        fill_single_cell_organelle_data_frame(properties_dataset, index, organelle_props, nucleus_props)

    return organelle_props


def set_single_cell_marker_props(properties_dataset, index, single_cell_mask, im_marker):
    regions = skimage.measure.regionprops(single_cell_mask, intensity_image=im_marker)
    marker_props = None
    if regions:
        marker_props = regions[-1]
        fill_single_cell_marker_polarity(properties_dataset, index, marker_props)

    return marker_props


def set_single_cell_marker_nuclei_props(properties_dataset, index, single_nucleus_mask, im_marker):
    regions = skimage.measure.regionprops(single_nucleus_mask, intensity_image=im_marker)
    marker_nuc_props = None
    if regions:
        marker_nuc_props = regions[-1]
        fill_single_cell_marker_nuclei_data_frame(properties_dataset, index, marker_nuc_props)

    return marker_nuc_props


def set_single_cell_marker_cytosol_props(properties_dataset, index, single_cytosol_mask, im_marker):
    regions = skimage.measure.regionprops(single_cytosol_mask.astype(int), intensity_image=im_marker)
    marker_nuc_cyt_props = None
    if regions:
        marker_nuc_cyt_props = regions[-1]
        fill_single_cell_marker_nuclei_cytosol_data_frame(properties_dataset, index, marker_nuc_cyt_props)

    return marker_nuc_cyt_props


def set_single_cell_marker_membrane_props(properties_dataset, index, single_membrane_mask, im_marker):
    regions = skimage.measure.regionprops(single_membrane_mask.astype(int), intensity_image=im_marker)
    marker_membrane_props = None
    if regions:
        marker_membrane_props = regions[-1]
        fill_single_cell_marker_membrane_data_frame(properties_dataset, index, marker_membrane_props)

    return marker_membrane_props


def get_features_from_cellpose_seg_multi_channel(parameters, img, cell_mask, filename, output_path):
    """Extracts features from a cellpose segmentation based on parameters given."""

    # initialize graph - no features associated with nodes
    rag = orientation_graph_nf(cell_mask)
    rag, cell_mask_rem_island = remove_islands(rag, cell_mask)

    get_logger().info("Number of RAG nodes: %s " % len(list(rag.nodes)))

    # get masks
    nuclei_mask = get_nuclei_mask(parameters, img, cell_mask_rem_island)
    organelle_mask = get_organelle_mask(parameters, img, cell_mask_rem_island)
    im_marker = get_image_marker(parameters, img)

    properties_dataset = pd.DataFrame()

    excluded = 0
    # iterate through each unique segmented cell
    for index, connected_component_label in enumerate(np.unique(cell_mask_rem_island)):
        # ignore background
        if connected_component_label == 0:
            continue

        # get single cell masks
        single_cell_mask = get_single_cell_mask(connected_component_label, cell_mask_rem_island)
        single_nucleus_mask = get_single_cell_nucleus_mask(connected_component_label, nuclei_mask)  # can be None
        single_organelle_mask = get_single_cell_organelle_mask(connected_component_label, organelle_mask)  # can be None
        single_membrane_mask = get_single_cell_membrane_mask(parameters, im_marker, single_cell_mask)  # can be None
        single_cytosol_mask = get_single_cell_cytosol_mask(single_cell_mask, im_marker,
                                                           single_nucleus_mask)  # can be None

        # threshold
        if threshold(
                parameters,
                single_cell_mask,
                single_nucleus_mask=single_nucleus_mask,
                single_organelle_mask=single_organelle_mask
        ):
            get_logger().info("Cell \"%s\" falls under threshold! Removed from RAG..." % connected_component_label)
            excluded += 1
            # remove a cell from the RAG
            rag.remove_node(connected_component_label)
            continue

        # properties for single cell
        set_single_cell_props(properties_dataset, index, filename, connected_component_label, single_cell_mask, rag)

        # properties for nucleus:
        if single_nucleus_mask is not None:
            nucleus_props = set_single_cell_nucleus_props(properties_dataset, index, single_nucleus_mask)

            # properties for organelle
            if nucleus_props and single_organelle_mask is not None:
                set_single_cell_organelle_props(properties_dataset, index, single_organelle_mask, nucleus_props)

        # properties for marker
        if im_marker is not None:
            set_single_cell_marker_props(properties_dataset, index, single_cell_mask, im_marker)
            set_single_cell_marker_membrane_props(properties_dataset, index, single_membrane_mask, im_marker)

            # marker nuclei properties
            if nuclei_mask is not None:
                set_single_cell_marker_nuclei_props(properties_dataset, index, single_nucleus_mask, im_marker)
                set_single_cell_marker_cytosol_props(properties_dataset, index, single_cytosol_mask, im_marker)

        # append feature of interest to the RAG node for being able to do further analysis
        foi = properties_dataset.at[index, parameters["feature_of_interest"]]
        foi_name = str(parameters["feature_of_interest"])
        rag.nodes[connected_component_label.astype('int')][foi_name] = foi

        get_logger().info(
            " ".join(
                str(x) for x in ["Cell %s - feature \"%s\": %s" % (connected_component_label, foi_name, foi)]
            )
        )

    get_logger().info("Excluded cells: %s" % str(excluded))
    get_logger().info("Leftover cells: %s" % str(len(np.unique(cell_mask)) - excluded))

    # morans I analysis based on FOI
    fill_morans_i(properties_dataset, rag, parameters["feature_of_interest"])

    # neighborhood feature analysis
    k_neighbor_dif(rag, parameters["feature_of_interest"], properties_dataset)

    plot_dataset(
        parameters, img, properties_dataset, output_path, filename, cell_mask_rem_island, nuclei_mask, organelle_mask,
        im_marker
    )

    return properties_dataset


def fill_morans_i(properties_dataset, rag, foi):
    get_logger().info("Calculating morans I group statistic...")

    # extract FOI and weights
    morans_feature, morans_weights = morans_data_prep(rag, foi)
    morans_i = run_morans(morans_feature, morans_weights)

    get_logger().info("Morans I value: %s " % morans_i.I)
    get_logger().info("Morans I p norm: %s " % morans_i.p_norm)

    # extend dataset
    properties_dataset["morans_i"] = [morans_i.I] * len(properties_dataset)
    properties_dataset["morans_p_norm"] = [morans_i.p_norm] * len(properties_dataset)


def get_single_cell_prop(single_cell_mask):
    """Gets the single cell properties."""
    # we are looking at a single cell. There is only one region!
    regions = skimage.measure.regionprops(single_cell_mask)
    if len(regions) > 1:
        raise ValueError("Too many regions for a single cell!")
    props = regions[-1]

    return props


def fill_single_cell_marker_membrane_data_frame(properties_dataset, index, marker_membrane_props):
    """Fills the dataset with the single cell marker membrane properties."""
    properties_dataset.at[index, "marker_mean_expression_mem"] = marker_membrane_props.mean_intensity
    properties_dataset.at[
        index, "marker_sum_expression_mem"] = marker_membrane_props.mean_intensity * marker_membrane_props.area


def fill_single_cell_marker_nuclei_cytosol_data_frame(properties_dataset, index, marker_nuc_cyt_props):
    """Fills the dataset with the single cell marker nuclei cytosol properties."""
    properties_dataset.at[index, "marker_mean_expression_cyt"] = marker_nuc_cyt_props.mean_intensity
    properties_dataset.at[
        index, "marker_sum_expression_cyt"] = marker_nuc_cyt_props.mean_intensity * marker_nuc_cyt_props.area


def fill_single_cell_marker_nuclei_data_frame(properties_dataset, index, marker_nuc_props):
    """Fills the dataset with the single cell marker nuclei properties."""
    properties_dataset.at[index, "marker_mean_expression_nuc"] = marker_nuc_props.mean_intensity
    properties_dataset.at[
        index, "marker_sum_expression_nuc"] = marker_nuc_props.mean_intensity * marker_nuc_props.area


def fill_single_cell_organelle_data_frame(dataset, index, organelle_props, nucleus_props):
    """Fills the dataset with the single cell organelle properties."""
    x_organelle, y_organelle = organelle_props.centroid
    angle_rad = compute_marker_centroid_orientation_rad(
        nucleus_props.centroid[0], nucleus_props.centroid[1], x_organelle, y_organelle
    )
    dataset.at[index, "organelle_X"] = x_organelle
    dataset.at[index, "organelle_Y"] = y_organelle

    dataset.at[index, "organelle_distance"] = compute_marker_vector_norm(
        x_organelle, y_organelle, nucleus_props.centroid[0], nucleus_props.centroid[1]
    )
    dataset.at[index, "organelle_orientation_rad"] = angle_rad
    dataset.at[index, "organelle_orientation_deg"] = 180.0 * angle_rad / np.pi


def fill_single_cell_general_data_frame(dataset, index, filename, connected_component_label, props, neighbours):
    """Fills the dataset with the general single cell properties."""
    dataset.at[index, "filename"] = filename
    dataset.at[index, "label"] = connected_component_label
    dataset.at[index, "cell_X"] = props.centroid[0]
    dataset.at[index, "cell_Y"] = props.centroid[1]
    # note, the values of orientation from props are in [-pi/2,pi/2] with zero along the y-axis
    dataset.at[index, "cell_shape_orientation"] = np.pi / 2.0 + props.orientation
    # assumes flow from left to right anlong x-axis
    dataset.at[index, "cell_major_axis_length"] = props.major_axis_length
    dataset.at[index, "cell_minor_axis_length"] = props.minor_axis_length
    dataset.at[index, "cell_eccentricity"] = props.eccentricity
    dataset.at[index, "cell_major_to_minor_ratio"] = props.major_axis_length / props.minor_axis_length
    dataset.at[index, "cell_area"] = props.area
    dataset.at[index, "cell_perimeter"] = props.perimeter
    dataset.at[index, "cell_neighbors"] = neighbours


def fill_single_nucleus_data_frame(dataset, index, props):
    """Fills the dataset with the single cell nucleus properties."""
    dataset.at[index, "nuc_X"] = props.centroid[0]
    dataset.at[index, "nuc_Y"] = props.centroid[1]
    # note, the values of orientation from props are in [-pi/2,pi/2] with zero along the y-axis
    dataset.at[index, "nuc_shape_orientation"] = np.pi / 2.0 + props.orientation
    # assumes flow from left to right anlong x-axis
    dataset.at[index, "nuc_major_axis_length"] = props.major_axis_length
    dataset.at[index, "nuc_minor_axis_length"] = props.minor_axis_length
    dataset.at[index, "nuc_area"] = props.area
    dataset.at[index, "nuc_perimeter"] = props.perimeter
    dataset.at[index, "nuc_eccentricity"] = props.eccentricity
    dataset.at[index, "nuc_major_to_minor_ratio"] = props.major_axis_length / props.minor_axis_length


def fill_single_cell_marker_polarity(dataset, index, props):
    """Fills the dataset with the single cell marker properties."""
    marker_centroid_x, marker_centroid_y = props.weighted_centroid
    cell_x, cell_y = props.centroid
    angle_rad = compute_marker_centroid_orientation_rad(cell_x, cell_y, marker_centroid_x, marker_centroid_y)
    dataset.at[index, "marker_mean_expr"] = props.mean_intensity
    dataset.at[index, "marker_sum_expression"] = props.mean_intensity * props.area
    dataset.at[index, "marker_centroid_X"] = props.weighted_centroid[0]
    dataset.at[index, "marker_centroid_Y"] = props.weighted_centroid[1]
    dataset.at[index, "marker_centroid_orientation_rad"] = angle_rad  # FIXME: no sign correct?
    dataset.at[index, "marker_centroid_orientation_deg"] = 180.0 * angle_rad / np.pi


def compute_marker_vector_norm(cell_x, cell_y, marker_centroid_x, marker_centroid_y):
    """Computes the marker vector norm"""
    distance2 = (cell_x - marker_centroid_x) ** 2
    distance2 += (cell_y - marker_centroid_y) ** 2

    return np.sqrt(distance2)


def compute_marker_centroid_orientation_rad(cell_x, cell_y, marker_centroid_x, marker_centroid_y):
    """Computes the marker polarity radius"""
    vec_x = cell_x - marker_centroid_x
    vec_y = cell_y - marker_centroid_y
    organelle_orientation_rad = np.pi - np.arctan2(vec_x, vec_y)

    return organelle_orientation_rad


def remove_edges(mask):
    """DescribeMe"""
    segments = mask.astype("int")
    end_1, end_2 = mask.shape[0], mask.shape[1]

    start_1, start_2 = 0, 0

    start = np.empty(mask.shape[0])
    end = np.empty(mask.shape[0])

    start.fill(int(start_1))
    end.fill(int(end_1 - 1))

    lower_right_arr = np.asarray(range(start_1, end_1))

    # list of points with 0, 0 - max
    roof = np.asarray([start, lower_right_arr])
    # list of of points with max 0-max
    floor = np.asarray([end, lower_right_arr])
    # list of point 0-max, 0
    left_wall = np.asarray([lower_right_arr, start])
    # list of point 0-max, max
    right_wall = np.asarray([lower_right_arr, end])

    concat_arr = np.hstack([left_wall, right_wall, roof, floor]).astype("int")
    x_indexed = segments[(concat_arr[1, :]), (concat_arr[0, :])]

    for element in np.unique(x_indexed):
        segments[segments == element] = 0
        get_logger().info("Removed edge s: %s" % str(element))

    get_logger().info("Number of removed edges: %s" % len(np.unique(x_indexed)))

    return segments


def remove_islands(frame_graph, mask):
    """Remove unconnected cells (Cells without neighbours)."""
    # Get list of islands - nodes with no neighbours and remove them
    list_of_islands = []
    for nodes in frame_graph.nodes:
        if len(list(frame_graph.neighbors(nodes))) == 0:
            list_of_islands.append(nodes)

    get_logger().info("Removed number of islands: %s" % len(list_of_islands))

    # remove islands from image and graph
    for elemet in np.unique(list_of_islands):
        frame_graph.remove_node(elemet)
        mask[:, :][mask[:, :] == elemet] = 0

    return frame_graph, mask


def morans_data_prep(rag, feature):
    """Retrieves the Weights and the feature of interest from the rag."""
    weights = W.from_networkx(rag)
    # extract the feature of interest from the rag
    morans_features = [rag.nodes[nodes_idx][feature] for nodes_idx in list(rag.nodes)]

    return morans_features, weights


def run_morans(morans_features, weihgts):
    """Run morans I, measure of spatial correlation and significance."""
    mi = Moran(morans_features, weihgts, two_tailed=False)

    return mi


def k_neighbor_dif(graph, foi, properties_dataset):
    """Extracts neighborhood statics from graph."""
    get_logger().info("Calculating first and second nearest neighbor statistic...")

    mean_dif_first_neighbors = []
    median_dif_first_neighbors = []
    var_dif_first_neighbors = []
    range_dif_first_neighbors = []

    mean_dif_second_neighbors = []
    median_dif_second_neighbors = []
    var_dif_second_neighbors = []
    range_dif_second_neighbors = []
    for node in graph.nodes():

        # feature of interest, first and second_nearest_neighbors
        node_foi = graph.nodes[node][foi]
        first_nearest_neighbors, second_nearest_neighbors = second_neighbors(graph, node)

        # feature of interest in comparison with first_nearest neighbor nodes
        foi_first_nearest_diff = [graph.nodes[neighbor][foi] - node_foi for neighbor in first_nearest_neighbors]

        # feature of interest in comparison to second nearest neighbor nodes
        foi_second_nearest_diff = [graph.nodes[neighbor][foi] - node_foi for neighbor in second_nearest_neighbors]

        # append first nearest neighbors statistics
        if len(foi_first_nearest_diff) > 0:
            mean_dif_first_neighbors.append(np.mean(foi_first_nearest_diff))
            median_dif_first_neighbors.append(np.median(foi_first_nearest_diff))
            var_dif_first_neighbors.append(np.std(foi_first_nearest_diff))
            range_dif_first_neighbors.append(abs(np.min(foi_first_nearest_diff) - np.max(foi_first_nearest_diff)))
        else:
            mean_dif_first_neighbors.append(0)
            median_dif_first_neighbors.append(0)
            var_dif_first_neighbors.append(0)
            range_dif_first_neighbors.append(0)

        # append second nearest neighbors statistics
        if len(foi_second_nearest_diff) > 0:
            mean_dif_second_neighbors.append(np.mean(foi_second_nearest_diff))
            median_dif_second_neighbors.append(np.median(foi_second_nearest_diff))
            var_dif_second_neighbors.append(np.std(foi_second_nearest_diff))
            range_dif_second_neighbors.append(abs(np.min(foi_second_nearest_diff) - np.max(foi_second_nearest_diff)))
        else:
            mean_dif_second_neighbors.append(0)
            median_dif_second_neighbors.append(0)
            var_dif_second_neighbors.append(0)
            range_dif_second_neighbors.append(0)

    # fill properties for first nearest neighbors
    properties_dataset["mean_dif_1st_neighbors"] = mean_dif_first_neighbors
    properties_dataset["median_dif_1st_neighbors"] = median_dif_first_neighbors
    properties_dataset["stddev_dif_1st_neighbors"] = var_dif_first_neighbors
    properties_dataset["range_dif_1st_neighbors"] = range_dif_first_neighbors

    # fill properties for second nearest neighbors
    properties_dataset["mean_dif_2nd_neighbors"] = mean_dif_second_neighbors
    properties_dataset["median_dif_2nd_neighbors"] = median_dif_second_neighbors
    properties_dataset["stddev_dif_2nd_neighbors"] = var_dif_second_neighbors
    properties_dataset["range_dif_2nd_neighbors"] = range_dif_second_neighbors


def second_neighbors(rag, node):
    """Return all second and first nearest neighbor nodes of a given node in a rag"""
    first_nearest_list = shared_edges(rag, node)

    # get list of second nearest neighbors
    k_nearest_list_unfiltered = []
    for first_nearest in first_nearest_list:
        k_nearest = list(rag.neighbors(first_nearest))
        k_nearest_list_unfiltered += k_nearest

    first_nearest_set = set(first_nearest_list)
    k_nearest_set = set(k_nearest_list_unfiltered)

    # get the difference of both sets and remove the node of interest from the final set
    second_nearest_list = list(k_nearest_set.difference(first_nearest_set))
    second_nearest_list = list(filter(lambda n: n != node, second_nearest_list))

    return first_nearest_list, second_nearest_list


def shared_edges(rag, node):
    """Return all first nearest neighbor nodes of a given node in a rag"""
    first_nearest = list(set(rag.neighbors(node)))
    return first_nearest
