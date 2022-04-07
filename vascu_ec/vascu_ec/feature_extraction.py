import numpy as np
import pandas as pd
import scipy.ndimage as ndi
import skimage.filters
import skimage.io
import skimage.measure
import skimage.segmentation

from vascu_ec.utils.plot import plot_dataset
from vascu_ec.utils.rag import orientation_graph_nf
from vascu_ec.utils.seg import get_outline_from_mask
from vascu_ec.vascu_ec_logging import get_logger


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


def get_nuclei_mask(parameters, img, cellpose_mask):
    """Gets the nuclei mask."""
    if parameters["channel_nucleus"] >= 0:
        img_nuclei_blur = ndi.gaussian_filter(img[:, :, parameters["channel_nucleus"]], sigma=3)

        nuclei_mask = np.where(img_nuclei_blur > skimage.filters.threshold_otsu(img_nuclei_blur), True, False)

        nuclei_label = nuclei_mask * cellpose_mask

        return nuclei_label
    else:
        return None


def get_golgi_mask(parameters, img, cellpose_mask):
    """Gets the golgi mask."""
    if parameters["channel_golgi"] >= 0:
        img_golgi_blur = ndi.gaussian_filter(img[:, :, parameters["channel_golgi"]], sigma=3)

        golgi_mask = np.where(img_golgi_blur > skimage.filters.threshold_otsu(img_golgi_blur), True, False)

        golgi_label = golgi_mask * cellpose_mask

        return golgi_label
    else:
        return None


def threshold(parameters, single_cell_mask, single_nucleus_mask=None, single_golgi_mask=None):
    """Thresholds given single_cell_mask. Returns True if falls under threshold."""
    #TODO: check if this can be removed, we already remove small objects from the cellpose mask
    if len(single_cell_mask[single_cell_mask == 1]) < parameters["min_cell_size"]:
        return True

    if single_nucleus_mask is not None:
        if len(single_nucleus_mask[single_nucleus_mask == 1]) < parameters["min_nucleus_size"]:
            return True

    if single_golgi_mask is not None:
        if len(single_golgi_mask[single_golgi_mask == 1]) < parameters["min_golgi_size"]:
            return True

    return False


def get_single_cell_nucleus_mask(connected_component_label, nuclei_mask):
    """Gets the single cell nucleus mask."""
    if nuclei_mask is not None:
        return np.where(nuclei_mask == connected_component_label, 1, 0)
    else:
        return None


def get_single_cell_golgi_mask(connected_component_label, golgi_mask):
    """Gets the single cell golgi mask."""
    if golgi_mask is not None:
        return np.where(golgi_mask == connected_component_label, 1, 0)
    else:
        return None


def get_single_cell_mask(connected_component_label, cell_mask_rem_island):
    """Gets the single cell mask."""
    if cell_mask_rem_island is not None:
        return np.where(cell_mask_rem_island == connected_component_label, 1, 0)
    else:
        return None


def get_single_cell_membrane_mask(parameters, im_marker, single_cell_mask):
    """Gets the single cellmembrane mask."""
    if im_marker is not None:
        return get_outline_from_mask(single_cell_mask.astype(bool), parameters["membrane_thickness"])
    else:
        return None


def get_single_cell_cytosol_mask(single_cell_mask, im_marker, single_nucleus_mask):
    """Gets the cytosol mask."""
    if single_nucleus_mask is not None and im_marker is not None:
        return np.logical_xor(single_cell_mask.astype(bool), single_nucleus_mask.astype(bool))
    else:
        return None


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


def set_single_cell_golgi_props(properties_dataset, index, single_golgi_mask, nucleus_props):
    regions = skimage.measure.regionprops(single_golgi_mask)
    golgi_props = None
    if regions:
        golgi_props = regions[-1]
        fill_single_cell_golgi_data_frame(properties_dataset, index, golgi_props, nucleus_props)

    return golgi_props


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
    golgi_mask = get_golgi_mask(parameters, img, cell_mask_rem_island)
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
        single_golgi_mask = get_single_cell_golgi_mask(connected_component_label, golgi_mask)  # can be None
        single_membrane_mask = get_single_cell_membrane_mask(parameters, im_marker, single_cell_mask)  # can be None
        single_cytosol_mask = get_single_cell_cytosol_mask(single_cell_mask, im_marker,
                                                           single_nucleus_mask)  # can be None

        # threshold
        if threshold(
                parameters,
                single_cell_mask,
                single_nucleus_mask=single_nucleus_mask,
                single_golgi_mask=single_golgi_mask
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

            # properties for golgi
            if nucleus_props and single_golgi_mask is not None:
                set_single_cell_golgi_props(properties_dataset, index, single_golgi_mask, nucleus_props)

        # properties for marker
        if im_marker is not None:
            set_single_cell_marker_props(properties_dataset, index, single_cell_mask, im_marker)
            set_single_cell_marker_membrane_props(properties_dataset, index, single_membrane_mask, im_marker)

            # marker nuclei properties
            if nuclei_mask is not None:
                set_single_cell_marker_nuclei_props(properties_dataset, index, single_nucleus_mask, im_marker)
                set_single_cell_marker_cytosol_props(properties_dataset, index, single_cytosol_mask, im_marker)

        # append feature of interest to the RAG node for MAYBE being able to do further analysis
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

    plot_dataset(
        parameters, img, properties_dataset, output_path, filename, cell_mask_rem_island, nuclei_mask, golgi_mask,
        im_marker
    )

    return properties_dataset


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
    properties_dataset.at[index, "mean_expression_mem"] = marker_membrane_props.mean_intensity
    properties_dataset.at[
        index, "sum_expression_mem"] = marker_membrane_props.mean_intensity * marker_membrane_props.area


def fill_single_cell_marker_nuclei_cytosol_data_frame(properties_dataset, index, marker_nuc_cyt_props):
    """Fills the dataset with the single cell marker nuclei cytosol properties."""
    properties_dataset.at[index, "mean_expression_cyt"] = marker_nuc_cyt_props.mean_intensity
    properties_dataset.at[
        index, "sum_expression_cyt"] = marker_nuc_cyt_props.mean_intensity * marker_nuc_cyt_props.area


def fill_single_cell_marker_nuclei_data_frame(properties_dataset, index, marker_nuc_props):
    """Fills the dataset with the single cell marker nuclei properties."""
    properties_dataset.at[index, "mean_expression_nuc"] = marker_nuc_props.mean_intensity
    properties_dataset.at[
        index, "sum_expression_nuc"] = marker_nuc_props.mean_intensity * marker_nuc_props.area


def fill_single_cell_golgi_data_frame(dataset, index, golgi_props, nucleus_props):
    """Fills the dataset with the single cell golgi properties."""
    x_golgi, y_golgi = golgi_props.centroid
    angle_rad, vec_x, vec_y = compute_marker_polarity_rad(
        nucleus_props.centroid[0], nucleus_props.centroid[1], x_golgi, y_golgi
    )
    dataset.at[index, "X_golgi"] = x_golgi
    dataset.at[index, "Y_golgi"] = y_golgi

    dataset.at[index, "distance"] = compute_marker_vector_norm(
        x_golgi, y_golgi, nucleus_props.centroid[0], nucleus_props.centroid[1]
    )
    dataset.at[index, "vec_X"] = vec_x
    dataset.at[index, "vec_Y"] = vec_y
    dataset.at[index, "angle_rad"] = angle_rad
    dataset.at[index, "flow_alignment"] = np.sin(angle_rad)
    dataset.at[index, "angle_deg"] = 180.0 * angle_rad / np.pi


def fill_single_cell_general_data_frame(dataset, index, filename, connected_component_label, props, neighbours):
    """Fills the dataset with the general single cell properties."""
    dataset.at[index, "filename"] = filename
    dataset.at[index, "label"] = connected_component_label
    dataset.at[index, "X_cell"] = props.centroid[0]
    dataset.at[index, "Y_cell"] = props.centroid[1]
    # note, the values of orientation from props are in [-pi/2,pi/2] with zero along the y-axis
    dataset.at[index, "shape_orientation"] = np.pi / 2.0 - props.orientation
    # assumes flow from left to right anlong x-axis
    dataset.at[index, "flow_shape_alignment"] = np.sin(props.orientation)
    dataset.at[index, "major_axis_length"] = props.major_axis_length
    dataset.at[index, "minor_axis_length"] = props.minor_axis_length
    dataset.at[index, "eccentricity"] = props.eccentricity
    dataset.at[index, "major_to_minor_ratio"] = props.major_axis_length / props.minor_axis_length
    dataset.at[index, "area"] = props.area
    dataset.at[index, "perimeter"] = props.perimeter
    dataset.at[index, "n_neighbors"] = neighbours


def fill_single_nucleus_data_frame(dataset, index, props):
    """Fills the dataset with the single cell nucleus properties."""
    dataset.at[index, "X_nuc"] = props.centroid[0]
    dataset.at[index, "Y_nuc"] = props.centroid[1]
    # note, the values of orientation from props are in [-pi/2,pi/2] with zero along the y-axis
    dataset.at[index, "shape_orientation_nuc"] = np.pi / 2.0 - props.orientation
    # assumes flow from left to right anlong x-axis
    dataset.at[index, "flow_shape_alignment_nuc"] = np.sin(np.pi / 2.0 - props.orientation)
    dataset.at[index, "major_axis_length_nuc"] = props.major_axis_length
    dataset.at[index, "minor_axis_length_nuc"] = props.minor_axis_length
    dataset.at[index, "area"] = props.area
    dataset.at[index, "perimeter"] = props.perimeter
    dataset.at[index, "eccentricity_nuc"] = props.eccentricity
    dataset.at[index, "major_to_minor_ratio_nuc"] = props.major_axis_length / props.minor_axis_length


def fill_single_cell_marker_polarity(dataset, index, props):
    """Fills the dataset with the single cell marker properties."""
    x_weighted, y_weighted = props.weighted_centroid
    x_cell, y_cell = props.centroid
    angle_rad, _, _ = compute_marker_polarity_rad(x_cell, y_cell, x_weighted, y_weighted)
    dataset.at[index, "mean_expression"] = props.mean_intensity
    dataset.at[index, "sum_expression"] = props.mean_intensity * props.area
    dataset.at[index, "X_weighted"] = props.weighted_centroid[0]
    dataset.at[index, "Y_weighted"] = props.weighted_centroid[1]
    dataset.at[index, "maker_vec_norm"] = compute_marker_vector_norm(x_cell, y_cell, x_weighted, y_weighted)
    dataset.at[index, "marker_polarity_rad"] = angle_rad  # FIXME: no sign correct?
    dataset.at[index, "marker_polarity_deg"] = 180.0 * angle_rad / np.pi


def compute_marker_vector_norm(x_cell, y_cell, x_weighted, y_weighted):
    """Computes the marker vector norm"""
    distance2 = (x_cell - x_weighted) ** 2
    distance2 += (y_cell - y_weighted) ** 2

    return np.sqrt(distance2)


def compute_marker_polarity_rad(x_cell, y_cell, x_weighted, y_weighted):
    """Computes the marker polarity radius"""
    vec_x = x_weighted - x_cell
    vec_y = y_weighted - y_cell
    angle_rad_ = np.arctan2(vec_x, vec_y)
    angle_rad = angle_rad_

    if angle_rad_ < 0.0:
        angle_rad = 2.0 * np.pi + angle_rad_

    return angle_rad, vec_x, vec_y


def remove_edges(mask):
    """DescribeMe"""
    segments = mask.astype("int")
    end_1, end_2 = mask.shape[0], mask.shape[1]

    start_1, start_2 = 0, 0

    the_start = np.empty(mask.shape[0])
    the_start.fill(int(start_1))
    the_end = np.empty(mask.shape[0])
    the_end.fill(int(end_1 - 1))
    lower_right_arr = np.asarray(range(start_1, end_1))

    # list of points with 0, 0 - max
    roof = np.asarray([the_start, lower_right_arr])
    # list of of points with max 0-max
    floor = np.asarray([the_end, lower_right_arr])
    # list of point 0-max, 0
    left_wall = np.asarray([lower_right_arr, the_start])
    # list of point 0-max, max
    right_wall = np.asarray([lower_right_arr, the_end])

    concat_arr = np.hstack([left_wall, right_wall, roof, floor]).astype("int")
    x_indexed = segments[(concat_arr[1, :]), (concat_arr[0, :])]

    for element in np.unique(x_indexed):
        segments[segments == element] = 0
        get_logger().info("Removed edge s: %s" % str(element))

    get_logger().info("Number of removed edges: %s" % len(np.unique(x_indexed)))

    return segments


def remove_islands(frame_graph, mask):
    """Remove unconnected cells (Cells without neighbours."""
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
