import cellpose.io
import cellpose.models
import cellpose.utils
import numpy as np
import pandas as pd
import scipy.ndimage as ndi
import skimage.filters
import skimage.io
import skimage.measure
import skimage.segmentation
from skimage.future.graph import RAG
from skimage.measure import regionprops

from vascu_ec.logging import get_logger
from vascu_ec.utils import plot


def get_image_for_segmentation(parameters, img):
    """Returns segmentation (junction, nucleus) if multichannel, else junction only."""

    ch_junction = int(parameters["channel_junction"])
    ch_nucleus = int(parameters["channel_nucleus"])

    print(str(ch_junction))
    print(img.shape)

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


def get_nuclei_cellpose(parameters, img, cellpose_mask):
    """Gets the nuclei cellpose mask."""
    model = cellpose.models.Cellpose(gpu=parameters["use_gpu"], model_type='nuclei')

    channels = [0, 0]

    masks, flows, styles, diams = model.eval(img, diameter=parameters["estimated_cell_diameter"], channels=channels)
    # masks, flows, styles, diams = model.eval(img, channels=channels)

    nuclei_mask = np.where(masks > 0, True, False)

    nuclei_label = nuclei_mask * cellpose_mask

    return nuclei_label


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
    if len(single_cell_mask[single_cell_mask == 1]) < parameters["min_cell_size"]:
        return True

    if single_nucleus_mask:
        if len(single_nucleus_mask[single_nucleus_mask == 1]) < parameters["min_nucleus_size"]:
            return True

    if single_golgi_mask:
        if len(single_golgi_mask[single_golgi_mask == 1]) < parameters["min_golgi_size"]:
            return True

    return False


def get_single_cell_nucleus_mask(connected_component_label, nuclei_mask):
    """Gets the single cell nucleus mask."""
    if nuclei_mask:
        return np.where(nuclei_mask == connected_component_label, 1, 0)
    else:
        return None


def get_single_cell_golgi_mask(connected_component_label, golgi_mask):
    """Gets the single cell golgi mask."""
    if golgi_mask:
        return np.where(golgi_mask == connected_component_label, 1, 0)
    else:
        return None


def get_single_cell_mask(connected_component_label, cell_mask_rem_island):
    """Gets the single cell mask."""
    if cell_mask_rem_island:
        return np.where(cell_mask_rem_island == connected_component_label, 1, 0)
    else:
        return None


def get_single_cell_membrane_mask(parameters, im_marker, single_cell_mask):
    """Gets the single cellmembrane mask."""
    if im_marker:
        return get_outline_from_mask(single_cell_mask.astype(bool), parameters["membrane_thickness"])
    else:
        return None


def get_single_cell_cytosol_mask(single_cell_mask, im_marker, single_nucleus_mask):
    """Gets the cytosol mask."""
    if single_nucleus_mask and im_marker:
        return np.logical_xor(single_cell_mask.astype(bool), single_nucleus_mask.astype(bool))
    else:
        return None


def get_image_marker(parameters, img):
    """Gets the image marker."""
    if parameters["channel_expression_marker"] >= 0:
        return img[:, :, parameters["channel_expression_marker"]]


def get_features_from_cellpose_seg_multi_channel(parameters, img, cell_mask, filename, output_path):
    """Extracts features from a cellpose segmentation based on parameters given."""
    rag = orientation_graph_nf(cell_mask)
    rag, cell_mask_rem_island = remove_islands(rag, cell_mask)

    # initialize graph - no features associated with nodes
    graph_nf = orientation_graph_nf(cell_mask_rem_island)
    get_logger().info(list(rag.nodes))

    # get masks
    nuclei_mask = get_nuclei_mask(parameters, img, cell_mask_rem_island)
    golgi_mask = get_golgi_mask(parameters, img, cell_mask_rem_island)
    im_marker = get_image_marker(parameters, img)

    properties_dataset = pd.DataFrame()

    # iterate through each unique segmented cell pixel,
    for index, connected_component_label in enumerate(np.unique(cell_mask_rem_island)):
        # ignore background
        if connected_component_label == 0:
            continue

        # get masks
        single_cell_mask = get_single_cell_mask(connected_component_label, cell_mask_rem_island)
        single_nucleus_mask = get_single_cell_nucleus_mask(connected_component_label, nuclei_mask)  # can be None
        single_golgi_mask = get_single_cell_golgi_mask(connected_component_label, golgi_mask)  # can be None
        single_membrane_mask = get_single_cell_membrane_mask(parameters, im_marker, single_cell_mask)  # can be None
        single_cytosol_mask = get_single_cell_cytosol_mask(single_cell_mask, im_marker, single_nucleus_mask)  # can be None

        # threshold
        if threshold(parameters, single_nucleus_mask, single_golgi_mask, single_cell_mask):
            rag.remove_node(connected_component_label)
            continue

        # properties for single cell
        single_cell_props = get_single_cell_prop(single_cell_mask)
        neighbours = len(list(graph_nf.neighbors(connected_component_label)))
        fill_single_cell_general_data_frame(
            properties_dataset, index, filename, connected_component_label, single_cell_props, neighbours
        )

        # properties for nucleus:
        if single_nucleus_mask:
            regions = skimage.measure.regionprops(single_nucleus_mask)
            nucleus_props = regions[-1]
            fill_single_nucleus_data_frame(single_nucleus_mask, index, nucleus_props)

        # properties for marker
        if im_marker:
            regions = skimage.measure.regionprops(single_cell_mask, intensity_image=im_marker)
            marker_props = regions[-1]
            fill_single_cell_marker_polarity(properties_dataset, index, marker_props)

            # marker nuclei properties
            if nuclei_mask:
                regions = skimage.measure.regionprops(single_nucleus_mask, intensity_image=im_marker)
                marker_nuc_props = regions[-1]
                fill_single_cell_marker_nuclei_data_frame(properties_dataset, index, marker_nuc_props)

                regions = skimage.measure.regionprops(single_cytosol_mask.astype(int), intensity_image=im_marker)
                marker_nuc_cyt_props = regions[-1]
                fill_single_cell_marker_nuclei_cytosol_data_frame(properties_dataset, index, marker_nuc_cyt_props)

            regions = skimage.measure.regionprops(single_membrane_mask.astype(int), intensity_image=im_marker)
            marker_membrane_props = regions[-1]
            fill_single_cell_marker_membrane_data_frame(properties_dataset, index, marker_membrane_props)

        # properties for golgi
        if single_nucleus_mask and single_golgi_mask:
            regions = skimage.measure.regionprops(single_golgi_mask)
            golgi_props = regions[-1]
            fill_single_cell_golgi_data_frame(properties_dataset, index, golgi_props, nucleus_props)

        # rag.nodes["label"][feature_of_interest] = single_cell_props.at[counter, feature_of_interest]
        f2a = properties_dataset.at[index, parameters["feature_of_interest"]]
        foe = str(parameters["feature_of_interest"])
        rag.nodes[connected_component_label.astype('int')][foe] = f2a
        # nw.set_node_attributes(graph_nf, {label.astype('int'):f2a}, parameters["feature_of_interest"])
        get_logger().info(f2a, parameters["feature_of_interest"], connected_component_label,
                          rag.nodes[connected_component_label.astype('int')][foe])
        # print(nx.get_node_attributes(G,parameters["feature_of_interest"]))

    plot_dataset(
        parameters, img, properties_dataset, output_path, filename, cell_mask_rem_island, nuclei_mask, golgi_mask
    )

    return properties_dataset


def plot_dataset(parameters, img, properties_dataset, output_path, filename, cell_mask_rem_island, nuclei_mask,
                 golgi_mask):
    """Plots the properties dataset"""
    im_junction = img[:, :, int(parameters["channel_junction"])]
    im_marker = img[:, :, int(parameters["channel_expression_marker"])]

    if parameters["plot_polarity"] and (parameters["channel_golgi"] >= 0):
        plot.plot_polarity(
            parameters,
            im_junction,
            [cell_mask_rem_island, nuclei_mask, golgi_mask],
            properties_dataset,
            filename, output_path
        )
    if parameters["plot_marker"] and (parameters["channel_nucleus"] >= 0):
        plot.plot_marker(
            parameters, im_marker, [cell_mask_rem_island, nuclei_mask], properties_dataset, filename, output_path
        )
        plot.plot_marker_polarity(
            parameters, im_marker, [cell_mask_rem_island], properties_dataset, filename,
            output_path
        )
    if parameters["plot_marker"] and (parameters["channel_nucleus"] < 0):
        plot.plot_marker(parameters, im_marker, [cell_mask_rem_island], properties_dataset, filename, output_path)
        plot.plot_marker_polarity(
            parameters, im_marker, [cell_mask_rem_island], properties_dataset, filename, output_path
        )
    if parameters["plot_alignment"] and (parameters["channel_nucleus"] >= 0):
        plot.plot_alignment(
            parameters,
            im_junction,
            [cell_mask_rem_island, nuclei_mask],
            properties_dataset,
            filename,
            output_path
        )
    if parameters["plot_marker"] and (parameters["channel_nucleus"] < 0):
        plot.plot_alignment(
            parameters,
            im_junction,
            [cell_mask_rem_island],
            properties_dataset,
            filename,
            output_path
        )
    if parameters["plot_ratio_method"]:
        plot.plot_ratio_method(
            parameters,
            im_junction,
            [cell_mask_rem_island],
            properties_dataset,
            filename,
            output_path
        )


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
    dataset.at[index, "shape_orientation"] = np.sin(np.pi / 2.0 - props.orientation)
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


######################
def get_outline_from_mask(mask, width=1):
    """DescribeMe"""
    dilated_mask = ndi.morphology.binary_dilation(mask.astype(bool), iterations=width)
    eroded_mask = ndi.morphology.binary_erosion(mask.astype(bool), iterations=width)
    outline_mask = np.logical_xor(dilated_mask, eroded_mask)

    return outline_mask


def orientation_graph_nf(img):
    """Gets the RegionAdjecencyGraph for an Image """
    rag = RAG(img.astype("int"))
    rag.remove_node(0)
    return (rag)


def orientation_graph(img):
    """DescribeMe"""
    rag = RAG(img.astype("int"))
    rag.remove_node(0)

    regions = regionprops(img.astype("int"))
    for region in regions:
        rag.nodes[region['label']]['orientation'] = region['orientation']
        rag.nodes[region['label']]['area'] = region['area']
        rag.nodes[region['label']]['polarity'] = region['major_axis_length'] / region['minor_axis_length']
        rag.nodes[region['label']]['aspect_ratio'] = (region["bbox"][2] - region["bbox"][0]) / (
                region["bbox"][3] - region["bbox"][1])
    return (rag)


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

    for elemet in np.unique(x_indexed):
        segments[segments == elemet] = 0
    return (segments)


def remove_islands(frame_graph, mask):
    """Remove unconnected masks"""
    list_of_islands = []
    # Get list of islands - nodes with no neighbours,
    # remove nodes with neighbours
    # frame_graph = orientation_graph(dat_no_edges[i,:,:])
    for nodes in frame_graph.nodes:
        if len(list(frame_graph.neighbors(nodes))) == 0:
            list_of_islands.append(nodes)

    print(list_of_islands)
    # remove islands from image and graph
    for elemet in np.unique(list_of_islands):
        frame_graph.remove_node(elemet)
        mask[:, :][mask[:, :] == elemet] = 0
    return frame_graph, mask
