import numpy as np
import pandas as pd

from polarityjam.polarityjam_logging import get_logger
from polarityjam.utils.compute import compute_reference_target_orientation_rad, compute_marker_vector_norm
from polarityjam.utils.moran import run_morans
from polarityjam.utils.weights import W


class Collector:
    """Collects output"""
    def __init__(self):
        self.dataset = pd.DataFrame()

    def fill_sc_props(self, sc_props, index, filename, connected_component_label):
        self.fill_sc_general_data_frame(index, filename, connected_component_label, sc_props.single_cell_props)

        if sc_props.marker_props:
            self.fill_sc_marker_polarity(index, sc_props.marker_props)

        if sc_props.nucleus_props:
            self.fill_sc_nucleus_data_frame(index, sc_props.nucleus_props)

        if sc_props.organelle_props and sc_props.nucleus_props:
            self.fill_sc_organelle_data_frame(index, sc_props.organelle_props, sc_props.nucleus_props)

        if sc_props.marker_nuc_props:
            self.fill_sc_marker_nuclei_data_frame(index, sc_props.marker_nuc_props)

        if sc_props.marker_nuc_cyt_props:
            self.fill_sc_marker_nuclei_cytosol_data_frame(index, sc_props.marker_nuc_cyt_props)

        if sc_props.marker_membrane_props:
            self.fill_sc_marker_membrane_data_frame(index, sc_props.marker_membrane_props)

        if sc_props.junction_props:
            self.fill_sc_junction_data_frame(index, sc_props.junction_props)

        if sc_props.junction_props:
            self.fill_sc_junction_sec_stat_data_frame(index, sc_props.junction_props)

    def fill_sc_marker_polarity(self, index, props):
        """Fills the dataset with the single cell marker properties."""
        marker_centroid_x, marker_centroid_y = props.weighted_centroid
        cell_x, cell_y = props.centroid
        angle_rad = compute_reference_target_orientation_rad(cell_x, cell_y, marker_centroid_x, marker_centroid_y)
        self.dataset.at[index, "marker_mean_expr"] = props.mean_intensity
        self.dataset.at[index, "marker_sum_expression"] = props.mean_intensity * props.area
        self.dataset.at[index, "marker_centroid_X"] = props.weighted_centroid[0]
        self.dataset.at[index, "marker_centroid_Y"] = props.weighted_centroid[1]
        self.dataset.at[index, "marker_centroid_orientation_rad"] = angle_rad
        self.dataset.at[index, "marker_centroid_orientation_deg"] = 180.0 * angle_rad / np.pi

    def fill_sc_nucleus_data_frame(self, index, props):
        """Fills the dataset with the single cell nucleus properties."""
        self.dataset.at[index, "nuc_X"] = props.centroid[0]
        self.dataset.at[index, "nuc_Y"] = props.centroid[1]
        # compute nucleus displacement
        nucleus_displacement_orientation_rad = compute_reference_target_orientation_rad(
            self.dataset["cell_X"][index], self.dataset["cell_Y"][index], props.centroid[0], props.centroid[1],
        )
        self.dataset.at[index, "nuc_displacement_orientation_rad"] = nucleus_displacement_orientation_rad
        self.dataset.at[
            index, "nuc_displacement_orientation_deg"] = 180.0 * nucleus_displacement_orientation_rad / np.pi
        # note, the values of orientation from props are in [-pi/2,pi/2] with zero along the y-axis
        # TODO: should be nuc_shape_orientation_deg and rad ?
        self.dataset.at[index, "nuc_shape_orientation"] = np.pi / 2.0 + props.orientation
        # assumes flow from left to right anlong x-axis
        self.dataset.at[index, "nuc_major_axis_length"] = props.major_axis_length
        self.dataset.at[index, "nuc_minor_axis_length"] = props.minor_axis_length
        self.dataset.at[index, "nuc_area"] = props.area
        self.dataset.at[index, "nuc_perimeter"] = props.perimeter
        self.dataset.at[index, "nuc_eccentricity"] = props.eccentricity
        self.dataset.at[index, "nuc_major_to_minor_ratio"] = props.major_axis_length / props.minor_axis_length

    def fill_sc_general_data_frame(self, index, filename, connected_component_label, props):
        """Fills the dataset with the general single cell properties."""
        self.dataset.at[index, "filename"] = filename
        self.dataset.at[index, "label"] = connected_component_label
        self.dataset.at[index, "cell_X"] = props.centroid[0]
        self.dataset.at[index, "cell_Y"] = props.centroid[1]
        # note, the values of orientation from props are in [-pi/2,pi/2] with zero along the y-axis
        self.dataset.at[index, "cell_shape_orientation"] = np.pi / 2.0 + props.orientation
        # assumes flow from left to right anlong x-axis
        self.dataset.at[index, "cell_major_axis_length"] = props.major_axis_length
        self.dataset.at[index, "cell_minor_axis_length"] = props.minor_axis_length
        self.dataset.at[index, "cell_eccentricity"] = props.eccentricity
        self.dataset.at[index, "cell_major_to_minor_ratio"] = props.major_axis_length / props.minor_axis_length
        self.dataset.at[index, "cell_area"] = props.area
        self.dataset.at[index, "cell_perimeter"] = props.perimeter

    def fill_sc_organelle_data_frame(self, index, organelle_props, nucleus_props):
        """Fills the dataset with the single cell organelle properties."""
        x_organelle, y_organelle = organelle_props.centroid
        angle_rad = compute_reference_target_orientation_rad(
            nucleus_props.centroid[0], nucleus_props.centroid[1], x_organelle, y_organelle
        )
        self.dataset.at[index, "organelle_X"] = x_organelle
        self.dataset.at[index, "organelle_Y"] = y_organelle

        self.dataset.at[index, "organelle_distance"] = compute_marker_vector_norm(
            x_organelle, y_organelle, nucleus_props.centroid[0], nucleus_props.centroid[1]
        )
        self.dataset.at[index, "organelle_orientation_rad"] = angle_rad
        self.dataset.at[index, "organelle_orientation_deg"] = 180.0 * angle_rad / np.pi

    def fill_sc_marker_nuclei_data_frame(self, index, marker_nuc_props):
        """Fills the dataset with the single cell marker nuclei properties."""
        # compute marker nucleus orientation
        marker_nucleus_orientation_rad = compute_reference_target_orientation_rad(
            self.dataset["nuc_X"][index],
            self.dataset["nuc_Y"][index],
            self.dataset["marker_centroid_X"][index],
            self.dataset["marker_centroid_Y"][index]
        )

        self.dataset.at[index, "marker_mean_expression_nuc"] = marker_nuc_props.mean_intensity
        self.dataset.at[
            index, "marker_sum_expression_nuc"] = marker_nuc_props.mean_intensity * marker_nuc_props.area
        self.dataset.at[index, "marker_nucleus_orientation_rad"] = marker_nucleus_orientation_rad
        self.dataset.at[index, "marker_nucleus_orientation_deg"] = 180.0 * marker_nucleus_orientation_rad / np.pi

    def fill_sc_marker_nuclei_cytosol_data_frame(self, index, marker_nuc_cyt_props):
        """Fills the dataset with the single cell marker nuclei cytosol properties."""
        self.dataset.at[index, "marker_mean_expression_cyt"] = marker_nuc_cyt_props.mean_intensity
        self.dataset.at[
            index, "marker_sum_expression_cyt"] = marker_nuc_cyt_props.mean_intensity * marker_nuc_cyt_props.area

    def fill_sc_marker_membrane_data_frame(self, index, marker_membrane_props):
        """Fills the dataset with the single cell marker membrane properties."""
        self.dataset.at[index, "marker_mean_expression_mem"] = marker_membrane_props.mean_intensity
        self.dataset.at[
            index, "marker_sum_expression_mem"] = marker_membrane_props.mean_intensity * marker_membrane_props.area

    def fill_sc_junction_data_frame(self, index, sc_junction_props):
        """Fills the dataset with the single cell junction properties."""
        junctions_centroid_x, junctions_centroid_y = sc_junction_props.circular_junction_props.weighted_centroid

        self.dataset.at[index, "junction_centroid_X"] = junctions_centroid_x
        self.dataset.at[index, "junction_centroid_Y"] = junctions_centroid_y
        self.dataset.at[index, "junction_perimeter"] = sc_junction_props.interface_perimeter
        self.dataset.at[index, "junction_protein_area"] = sc_junction_props.junction_protein_area_props.area
        # dataset.at[index, "junction_fragmented_perimeter"] = sc_junction_props.junction_fragmented_perimeter
        self.dataset.at[index, "junction_mean_expression"] = sc_junction_props.interface_props.mean_intensity
        self.dataset.at[
            index, "junction_protein_intensity"] = sc_junction_props.interface_props.mean_intensity * sc_junction_props.interface_props.area

    def fill_sc_junction_sec_stat_data_frame(self, index, sc_junction_props):
        # dataset.at[index, "junction_coverage_index"] = sc_junction_props.junction_coverage_index
        self.dataset.at[index, "junction_interface_occupancy"] = sc_junction_props.junction_interface_occupancy
        self.dataset.at[
            index, "junction_intensity_per_interface_area"] = sc_junction_props.junction_intensity_per_interface_area
        self.dataset.at[index, "junction_cluster_density"] = sc_junction_props.junction_cluster_density

    def fill_morans_i(self, rag, foi):
        get_logger().info("Calculating morans I group statistic...")

        # extract FOI and weights
        weights = W.from_networkx(rag)
        # extract the feature of interest from the rag
        morans_features = [rag.nodes[nodes_idx][foi] for nodes_idx in list(rag.nodes)]

        morans_i = run_morans(morans_features, weights)

        get_logger().info("Morans I value: %s " % morans_i.I)
        get_logger().info("Morans I p norm: %s " % morans_i.p_norm)

        # extend dataset
        self.dataset["morans_i"] = [morans_i.I] * len(self.dataset)
        self.dataset["morans_p_norm"] = [morans_i.p_norm] * len(self.dataset)

    def fill_neighborhood_features(self, props):
        self.dataset["neighbors_cell"] = props["num_neighbours"]

        # fill properties for first nearest neighbors
        self.dataset["neighbors_mean_dif_1st"] = props["mean_dif_first_neighbors"]
        self.dataset["neighbors_median_dif_1st"] = props["median_dif_first_neighbors"]
        self.dataset["neighbors_stddev_dif_1st"] = props["var_dif_first_neighbors"]
        self.dataset["neighbors_range_dif_1st"] = props["range_dif_first_neighbors"]

        # fill properties for second nearest neighbors
        self.dataset["neighbors_mean_dif_2nd"] = props["mean_dif_second_neighbors"]
        self.dataset["neighbors_median_dif_2nd"] = props["median_dif_second_neighbors"]
        self.dataset["neighbors_stddev_dif_2nd"] = props["var_dif_second_neighbors"]
        self.dataset["neighbors_range_dif_2nd"] = props["range_dif_second_neighbors"]

def fill_neighborhood_features(dataset, props):
    dataset["neighbors_cell"] = props["num_neighbours"]

    # fill properties for first nearest neighbors
    dataset["neighbors_mean_dif_1st"] = props["mean_dif_first_neighbors"]
    dataset["neighbors_median_dif_1st"] = props["median_dif_first_neighbors"]
    dataset["neighbors_stddev_dif_1st"] = props["var_dif_first_neighbors"]
    dataset["neighbors_range_dif_1st"] = props["range_dif_first_neighbors"]

    # fill properties for second nearest neighbors
    dataset["neighbors_mean_dif_2nd"] = props["mean_dif_second_neighbors"]
    dataset["neighbors_median_dif_2nd"] = props["median_dif_second_neighbors"]
    dataset["neighbors_stddev_dif_2nd"] = props["var_dif_second_neighbors"]
    dataset["neighbors_range_dif_2nd"] = props["range_dif_second_neighbors"]


def fill_single_cell_marker_polarity(dataset, index, props):
    """Fills the dataset with the single cell marker properties."""
    marker_centroid_x, marker_centroid_y = props.weighted_centroid
    cell_x, cell_y = props.centroid
    angle_rad = compute_reference_target_orientation_rad(cell_x, cell_y, marker_centroid_x, marker_centroid_y)
    dataset.at[index, "marker_mean_expr"] = props.mean_intensity
    dataset.at[index, "marker_sum_expression"] = props.mean_intensity * props.area
    dataset.at[index, "marker_centroid_X"] = props.weighted_centroid[0]
    dataset.at[index, "marker_centroid_Y"] = props.weighted_centroid[1]
    dataset.at[index, "marker_centroid_orientation_rad"] = angle_rad
    dataset.at[index, "marker_centroid_orientation_deg"] = 180.0 * angle_rad / np.pi


def fill_single_nucleus_data_frame(dataset, index, props):
    """Fills the dataset with the single cell nucleus properties."""
    dataset.at[index, "nuc_X"] = props.centroid[0]
    dataset.at[index, "nuc_Y"] = props.centroid[1]
    # compute nucleus displacement
    nucleus_displacement_orientation_rad = compute_reference_target_orientation_rad(
        dataset["cell_X"][index], dataset["cell_Y"][index], props.centroid[0], props.centroid[1],
    )
    dataset.at[index, "nuc_displacement_orientation_rad"] = nucleus_displacement_orientation_rad
    dataset.at[index, "nuc_displacement_orientation_deg"] = 180.0 * nucleus_displacement_orientation_rad / np.pi
    # note, the values of orientation from props are in [-pi/2,pi/2] with zero along the y-axis
    # TODO: should be nuc_shape_orientation_deg and rad ?
    dataset.at[index, "nuc_shape_orientation"] = np.pi / 2.0 + props.orientation
    # assumes flow from left to right anlong x-axis
    dataset.at[index, "nuc_major_axis_length"] = props.major_axis_length
    dataset.at[index, "nuc_minor_axis_length"] = props.minor_axis_length
    dataset.at[index, "nuc_area"] = props.area
    dataset.at[index, "nuc_perimeter"] = props.perimeter
    dataset.at[index, "nuc_eccentricity"] = props.eccentricity
    dataset.at[index, "nuc_major_to_minor_ratio"] = props.major_axis_length / props.minor_axis_length


def fill_single_cell_general_data_frame(dataset, index, filename, connected_component_label, props):
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


def fill_single_cell_organelle_data_frame(dataset, index, organelle_props, nucleus_props):
    """Fills the dataset with the single cell organelle properties."""
    x_organelle, y_organelle = organelle_props.centroid
    angle_rad = compute_reference_target_orientation_rad(
        nucleus_props.centroid[0], nucleus_props.centroid[1], x_organelle, y_organelle
    )
    dataset.at[index, "organelle_X"] = x_organelle
    dataset.at[index, "organelle_Y"] = y_organelle

    dataset.at[index, "organelle_distance"] = compute_marker_vector_norm(
        x_organelle, y_organelle, nucleus_props.centroid[0], nucleus_props.centroid[1]
    )
    dataset.at[index, "organelle_orientation_rad"] = angle_rad
    dataset.at[index, "organelle_orientation_deg"] = 180.0 * angle_rad / np.pi


def fill_single_cell_marker_nuclei_data_frame(properties_dataset, index, marker_nuc_props):
    """Fills the dataset with the single cell marker nuclei properties."""
    # compute marker nucleus orientation
    marker_nucleus_orientation_rad = compute_reference_target_orientation_rad(
        properties_dataset["nuc_X"][index],
        properties_dataset["nuc_Y"][index],
        properties_dataset["marker_centroid_X"][index],
        properties_dataset["marker_centroid_Y"][index]
    )

    properties_dataset.at[index, "marker_mean_expression_nuc"] = marker_nuc_props.mean_intensity
    properties_dataset.at[
        index, "marker_sum_expression_nuc"] = marker_nuc_props.mean_intensity * marker_nuc_props.area
    properties_dataset.at[index, "marker_nucleus_orientation_rad"] = marker_nucleus_orientation_rad
    properties_dataset.at[index, "marker_nucleus_orientation_deg"] = 180.0 * marker_nucleus_orientation_rad / np.pi


def fill_single_cell_marker_nuclei_cytosol_data_frame(properties_dataset, index, marker_nuc_cyt_props):
    """Fills the dataset with the single cell marker nuclei cytosol properties."""
    properties_dataset.at[index, "marker_mean_expression_cyt"] = marker_nuc_cyt_props.mean_intensity
    properties_dataset.at[
        index, "marker_sum_expression_cyt"] = marker_nuc_cyt_props.mean_intensity * marker_nuc_cyt_props.area


def fill_single_cell_marker_membrane_data_frame(properties_dataset, index, marker_membrane_props):
    """Fills the dataset with the single cell marker membrane properties."""
    properties_dataset.at[index, "marker_mean_expression_mem"] = marker_membrane_props.mean_intensity
    properties_dataset.at[
        index, "marker_sum_expression_mem"] = marker_membrane_props.mean_intensity * marker_membrane_props.area


def fill_morans_i(properties_dataset, rag, foi):
    get_logger().info("Calculating morans I group statistic...")

    # extract FOI and weights
    weights = W.from_networkx(rag)
    # extract the feature of interest from the rag
    morans_features = [rag.nodes[nodes_idx][foi] for nodes_idx in list(rag.nodes)]

    morans_i = run_morans(morans_features, weights)

    get_logger().info("Morans I value: %s " % morans_i.I)
    get_logger().info("Morans I p norm: %s " % morans_i.p_norm)

    # extend dataset
    properties_dataset["morans_i"] = [morans_i.I] * len(properties_dataset)
    properties_dataset["morans_p_norm"] = [morans_i.p_norm] * len(properties_dataset)


def fill_single_cell_junction_data_frame(
        dataset, index, junction_perimeter, junction_props, junction_protein_area, weighted_centroid
):
    """Fills the dataset with the single cell junction properties."""
    junctions_centroid_x, junctions_centroid_y = weighted_centroid

    dataset.at[index, "junction_centroid_X"] = junctions_centroid_x
    dataset.at[index, "junction_centroid_Y"] = junctions_centroid_y
    dataset.at[index, "junction_perimeter"] = junction_perimeter
    dataset.at[index, "junction_protein_area"] = junction_protein_area
    # dataset.at[index, "junction_fragmented_perimeter"] = junction_fragmented_perimeter
    dataset.at[index, "junction_mean_expression"] = junction_props.mean_intensity
    dataset.at[index, "junction_protein_intensity"] = junction_props.mean_intensity * junction_props.area


def fill_single_cell_junction_sec_stat_data_frame(
        dataset, index, junction_interface_occupancy, junction_intensity_per_interface_area,
        cluster_density
):
    # dataset.at[index, "junction_coverage_index"] = junction_coverage_index
    dataset.at[index, "junction_interface_occupancy"] = junction_interface_occupancy
    dataset.at[index, "junction_intensity_per_interface_area"] = junction_intensity_per_interface_area
    dataset.at[index, "junction_cluster_density"] = cluster_density
