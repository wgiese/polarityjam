import json
from scipy import ndimage as ndi
import numpy as np
import pandas as pd
from polarityjam.compute.compute import otsu_thresh_mask

from polarityjam.model.masks import SingleCellMasksCollection
from polarityjam.model.properties import SingleCellCellProps, SingleCellNucleusProps, SingleCellOrganelleProps, \
    SingleCellMarkerProps, SingleCellMarkerMembraneProps, SingleCellMarkerNucleiProps, SingleCellMarkerCytosolProps, \
    SingleCellJunctionInterfaceProps, SingleCellJunctionProteinProps, SingleCellJunctionProteinCircularProps, \
    SingleCellJunctionProps, SingleCellPropertiesCollection
from polarityjam.utils import parameters


class PropertyCollector:
    """Collects features "as they come" in a large dataset. Not responsible for feature calculation!"""

    def __init__(self):
        self.dataset = pd.DataFrame()

    def collect_sc_props(self, sc_props, index, filename, connected_component_label):
        self.collect_sc_general_props(index, filename, connected_component_label, sc_props.single_cell_props)

        if sc_props.marker_props:
            self.collect_sc_marker_polarity_props(index, sc_props.marker_props)

        if sc_props.nucleus_props:
            self.collect_sc_nucleus_props(index, sc_props.nucleus_props)

        if sc_props.organelle_props:
            self.collect_sc_organelle_props(index, sc_props.organelle_props)

        if sc_props.marker_nuc_props:
            self.collect_sc_marker_nuclei_props(index, sc_props.marker_nuc_props)

        if sc_props.marker_nuc_cyt_props:
            self.collect_sc_marker_nuclei_cytosol_props(index, sc_props.marker_nuc_cyt_props)

        if sc_props.marker_membrane_props:
            self.collect_sc_marker_membrane_props(index, sc_props.marker_membrane_props)

        if sc_props.junction_props:
            self.collect_sc_junction_props(index, sc_props.junction_props)

        if sc_props.junction_props:
            self.collect_sc_junction_sec_stat_props(index, sc_props.junction_props)

    def collect_sc_marker_polarity_props(self, index, props):
        """Fills the dataset with the single cell marker properties."""
        self.dataset.at[index, "marker_mean_expr"] = props.mean_intensity
        self.dataset.at[index, "marker_sum_expression"] = props.mean_intensity * props.area
        self.dataset.at[index, "marker_centroid_X"] = props.weighted_centroid[0]
        self.dataset.at[index, "marker_centroid_Y"] = props.weighted_centroid[1]
        self.dataset.at[index, "marker_centroid_orientation_rad"] = props.marker_centroid_orientation_rad
        self.dataset.at[index, "marker_centroid_orientation_deg"] = props.marker_centroid_orientation_deg

    def collect_sc_nucleus_props(self, index, props):
        """Fills the dataset with the single cell nucleus properties."""
        self.dataset.at[index, "nuc_X"] = props.centroid[0]
        self.dataset.at[index, "nuc_Y"] = props.centroid[1]
        self.dataset.at[index, "nuc_displacement_orientation_rad"] = props.nuc_displacement_orientation_rad
        self.dataset.at[index, "nuc_displacement_orientation_deg"] = props.nuc_displacement_orientation_deg
        self.dataset.at[index, "nuc_shape_orientation"] = props.nuc_shape_orientation
        self.dataset.at[index, "nuc_major_axis_length"] = props.major_axis_length
        self.dataset.at[index, "nuc_minor_axis_length"] = props.minor_axis_length
        self.dataset.at[index, "nuc_area"] = props.area
        self.dataset.at[index, "nuc_perimeter"] = props.perimeter
        self.dataset.at[index, "nuc_eccentricity"] = props.eccentricity
        self.dataset.at[index, "nuc_major_to_minor_ratio"] = props.nuc_major_to_minor_ratio

    def collect_sc_general_props(self, index, filename, connected_component_label, props):
        """Fills the dataset with the general single cell properties."""
        self.dataset.at[index, "filename"] = filename
        self.dataset.at[index, "label"] = connected_component_label
        self.dataset.at[index, "cell_X"] = props.centroid[0]
        self.dataset.at[index, "cell_Y"] = props.centroid[1]
        self.dataset.at[index, "cell_shape_orientation"] = props.cell_shape_orientation
        self.dataset.at[index, "cell_major_axis_length"] = props.major_axis_length
        self.dataset.at[index, "cell_minor_axis_length"] = props.minor_axis_length
        self.dataset.at[index, "cell_eccentricity"] = props.eccentricity
        self.dataset.at[index, "cell_major_to_minor_ratio"] = props.cell_major_to_minor_ratio
        self.dataset.at[index, "cell_area"] = props.area
        self.dataset.at[index, "cell_perimeter"] = props.perimeter
        self.dataset.at[index, "cell_corner_points"] = json.dumps(props.cell_corner_points.tolist())

    def collect_sc_organelle_props(self, index, organelle_props):
        """Fills the dataset with the single cell organelle properties."""
        self.dataset.at[index, "organelle_X"] = organelle_props.centroid[0]
        self.dataset.at[index, "organelle_Y"] = organelle_props.centroid[1]
        self.dataset.at[index, "organelle_distance"] = organelle_props.organelle_distance
        self.dataset.at[index, "organelle_orientation_rad"] = organelle_props.organelle_orientation_rad
        self.dataset.at[index, "organelle_orientation_deg"] = organelle_props.organelle_orientation_deg

    def collect_sc_marker_nuclei_props(self, index, marker_nuc_props):
        """Fills the dataset with the single cell marker nuclei properties."""
        self.dataset.at[index, "marker_mean_expression_nuc"] = marker_nuc_props.mean_intensity
        self.dataset.at[index, "marker_sum_expression_nuc"] = marker_nuc_props.marker_sum_expression_nuc
        self.dataset.at[index, "marker_nucleus_orientation_rad"] = marker_nuc_props.marker_nucleus_orientation_rad
        self.dataset.at[index, "marker_nucleus_orientation_deg"] = marker_nuc_props.marker_nucleus_orientation_deg

    def collect_sc_marker_nuclei_cytosol_props(self, index, marker_nuc_cyt_props):
        """Fills the dataset with the single cell marker nuclei cytosol properties."""
        self.dataset.at[index, "marker_mean_expression_cyt"] = marker_nuc_cyt_props.mean_intensity
        self.dataset.at[index, "marker_sum_expression_cyt"] = marker_nuc_cyt_props.marker_sum_expression_cyt

    def collect_sc_marker_membrane_props(self, index, marker_membrane_props):
        """Fills the dataset with the single cell marker membrane properties."""
        self.dataset.at[index, "marker_mean_expression_mem"] = marker_membrane_props.mean_intensity
        self.dataset.at[index, "marker_sum_expression_mem"] = marker_membrane_props.marker_sum_expression_mem

    def collect_sc_junction_props(self, index, sc_junction_props):
        """Fills the dataset with the single cell junction properties."""
        junctions_centroid_x, junctions_centroid_y = sc_junction_props.sc_junction_protein_circular_props.weighted_centroid

        self.dataset.at[index, "junction_centroid_X"] = junctions_centroid_x
        self.dataset.at[index, "junction_centroid_Y"] = junctions_centroid_y
        self.dataset.at[index, "junction_perimeter"] = sc_junction_props.interface_perimeter
        self.dataset.at[index, "junction_protein_area"] = sc_junction_props.sc_junction_protein_props.area
        # dataset.at[index, "junction_fragmented_perimeter"] = sc_junction_props.junction_fragmented_perimeter
        self.dataset.at[
            index, "junction_mean_expression"] = sc_junction_props.sc_junction_interface_props.mean_intensity
        self.dataset.at[
            index, "junction_protein_intensity"] = sc_junction_props.junction_protein_intensity

    def collect_sc_junction_sec_stat_props(self, index, sc_junction_props):
        # dataset.at[index, "junction_coverage_index"] = sc_junction_props.junction_coverage_index
        self.dataset.at[index, "junction_interface_occupancy"] = sc_junction_props.junction_interface_occupancy
        self.dataset.at[
            index, "junction_intensity_per_interface_area"] = sc_junction_props.junction_intensity_per_interface_area
        self.dataset.at[index, "junction_cluster_density"] = sc_junction_props.junction_cluster_density

    def collect_morans_i_props(self, morans_i):
        # extend dataset
        self.dataset["morans_i"] = [morans_i.I] * len(self.dataset)
        self.dataset["morans_p_norm"] = [morans_i.p_norm] * len(self.dataset)

    def collect_neighborhood_props(self, neighborhood_props_list):
        for index, neighborhood_props in enumerate(neighborhood_props_list):
            index += 1  # offset from pandas df header
            self.dataset.at[index, "neighbors_cell"] = neighborhood_props.num_neighbours

            # fill properties for first nearest neighbors
            self.dataset.at[index, "neighbors_mean_dif_1st"] = neighborhood_props.mean_dif_first_neighbors
            self.dataset.at[index, "neighbors_median_dif_1st"] = neighborhood_props.median_dif_first_neighbors
            self.dataset.at[index, "neighbors_stddev_dif_1st"] = neighborhood_props.var_dif_first_neighbors
            self.dataset.at[index, "neighbors_range_dif_1st"] = neighborhood_props.range_dif_first_neighbors

            # fill properties for second-nearest neighbors
            self.dataset.at[index, "neighbors_mean_dif_2nd"] = neighborhood_props.mean_dif_second_neighbors
            self.dataset.at[index, "neighbors_median_dif_2nd"] = neighborhood_props.median_dif_second_neighbors
            self.dataset.at[index, "neighbors_stddev_dif_2nd"] = neighborhood_props.var_dif_second_neighbors
            self.dataset.at[index, "neighbors_range_dif_2nd"] = neighborhood_props.range_dif_second_neighbors


class SingleCellPropertyCollector:

    def __init__(self):
        pass

    def calc_sc_props(self, sc_masks, im_marker, im_junction):
        """calculates all properties for the single cell"""

        # properties for single cell
        sc_cell_props = self.calc_sc_cell_props(sc_masks.sc_mask.astype(int))

        # init optional properties
        sc_nuc_props = None
        sc_organelle_props = None
        sc_marker_props = None
        sc_marker_membrane_props = None
        sc_marker_nuclei_props = None
        sc_marker_cytosol_props = None
        sc_junction_props = None

        # properties for nucleus:
        if sc_masks.sc_nucleus_mask is not None:
            sc_nuc_props = self.calc_sc_nucleus_props(sc_masks.sc_nucleus_mask.astype(int), sc_cell_props)

            # properties for organelle
            if sc_nuc_props and sc_masks.sc_organelle_mask is not None:
                sc_organelle_props = self.calc_sc_organelle_props(sc_masks.sc_organelle_mask.astype(int), sc_nuc_props)

        # properties for marker
        if im_marker is not None:
            sc_marker_props = self.calc_sc_marker_props(sc_masks.sc_mask.astype(int), im_marker)
            sc_marker_membrane_props = self.calc_sc_marker_membrane_props(sc_masks.sc_membrane_mask.astype(int),
                                                                          im_marker)

            # properties for marker nuclei
            if sc_masks.sc_nucleus_mask is not None:
                sc_marker_nuclei_props = self.calc_sc_marker_nuclei_props(sc_masks.sc_nucleus_mask.astype(int),
                                                                          im_marker, sc_nuc_props, sc_marker_props)
                sc_marker_cytosol_props = self.calc_sc_marker_cytosol_props(sc_masks.sc_cytosol_mask.astype(int),
                                                                            im_marker)

        # properties for junctions
        if im_junction is not None:
            sc_junction_props = self.calc_sc_junction_props(
                sc_masks.sc_mask.astype(int),
                sc_masks.sc_membrane_mask.astype(int),
                sc_masks.sc_junction_protein_area_mask.astype(int),
                im_junction,
                sc_cell_props.minor_axis_length
            )

        return SingleCellPropertiesCollection(
            sc_cell_props,
            sc_nuc_props,
            sc_organelle_props,
            sc_marker_props,
            sc_marker_membrane_props,
            sc_marker_nuclei_props,
            sc_marker_cytosol_props,
            sc_junction_props
        )

    @staticmethod
    def calc_sc_cell_props(sc_mask):
        return SingleCellCellProps(sc_mask)

    @staticmethod
    def calc_sc_nucleus_props(sc_nucleus_maks, sc_props):
        return SingleCellNucleusProps(sc_nucleus_maks, sc_props)

    @staticmethod
    def calc_sc_organelle_props(sc_organelle_mask, nucleus_props):
        return SingleCellOrganelleProps(sc_organelle_mask, nucleus_props)

    @staticmethod
    def calc_sc_marker_props(sc_mask, im_marker):
        return SingleCellMarkerProps(sc_mask, im_marker)

    @staticmethod
    def calc_sc_marker_membrane_props(sc_membrane_mask, im_marker):
        return SingleCellMarkerMembraneProps(sc_membrane_mask, im_marker)

    @staticmethod
    def calc_sc_marker_nuclei_props(sc_nucleus_mask, im_marker, nucleus_props, marker_props):
        return SingleCellMarkerNucleiProps(sc_nucleus_mask, im_marker, nucleus_props, marker_props)

    @staticmethod
    def calc_sc_marker_cytosol_props(sc_cytosol_mask, im_marker):
        return SingleCellMarkerCytosolProps(sc_cytosol_mask, im_marker)

    @staticmethod
    def calc_sc_junction_props(sc_mask, single_membrane_mask, single_junction_protein_area_mask,
                               im_junction, cell_minor_axis_length):

        im_junction_protein_single_cell = otsu_thresh_mask(single_membrane_mask, im_junction)

        sc_junction_interface_props = SingleCellJunctionInterfaceProps(single_membrane_mask, im_junction)

        sc_junction_protein_props = SingleCellJunctionProteinProps(single_junction_protein_area_mask,
                                                                   im_junction_protein_single_cell)

        sc_junction_protein_circular_props = SingleCellJunctionProteinCircularProps(
            im_junction_protein_single_cell,
            cell_minor_axis_length,
            sc_junction_interface_props.centroid
        )

        return SingleCellJunctionProps(sc_junction_interface_props, sc_junction_protein_props,
                                       sc_junction_protein_circular_props, sc_mask)


class SingleCellMaskCollector:

    def __init__(self):
        pass

    def calc_sc_masks(self, masks, connected_component_label, im_junction):
        sc_mask = self.get_sc_mask(masks.cell_mask_rem_island, connected_component_label)

        sc_membrane_mask = self.get_sc_membrane_mask(sc_mask, parameters.membrane_thickness)

        # init optional sc masks
        sc_nucleus_mask = None
        sc_organelle_mask = None
        sc_cytosol_mask = None
        sc_junction_protein_mask = None

        if masks.nuclei_mask is not None:
            sc_nucleus_mask = self.get_sc_nucleus_mask(masks.nuclei_mask, connected_component_label)
            sc_cytosol_mask = self.get_sc_cytosol_mask(sc_mask, sc_nucleus_mask)

        if masks.organelle_mask is not None:
            sc_organelle_mask = self.get_sc_organelle_mask(masks.organelle_mask, connected_component_label)

        if im_junction is not None:
            sc_junction_protein_mask = self.get_sc_junction_protein_mask(
                sc_membrane_mask, im_junction
            )

        return SingleCellMasksCollection(
            connected_component_label,
            sc_mask,
            sc_nucleus_mask,
            sc_organelle_mask,
            sc_membrane_mask,
            sc_cytosol_mask,
            sc_junction_protein_mask
        )

    @staticmethod
    def get_sc_mask(cell_mask_rem_island, connected_component_label):
        """Gets the single cell mask from a mask where each cell has an increasing connected component value."""
        sc_mask = np.where(cell_mask_rem_island == connected_component_label, 1, 0)
        return sc_mask

    @staticmethod
    def get_sc_nucleus_mask(nuclei_mask, connected_component_label):
        """Gets the single cell nucleus mask."""
        sc_nucleus_mask = np.where(nuclei_mask == connected_component_label, 1, 0)
        return sc_nucleus_mask

    @staticmethod
    def get_sc_organelle_mask(organelle_mask, connected_component_label):
        """Gets the single cell organelle mask."""
        sc_organelle_mask = np.where(organelle_mask == connected_component_label, 1, 0)

        return sc_organelle_mask

    @staticmethod
    def get_sc_membrane_mask(sc_mask, membrane_thickness):
        """Gets the single cell membrane mask."""
        sc_membrane_mask = SingleCellMaskCollector._get_outline_from_mask(sc_mask, membrane_thickness)
        return sc_membrane_mask

    @staticmethod
    def get_sc_cytosol_mask(sc_mask, sc_nucleus_mask):
        """Gets the cytosol mask."""
        sc_cytosol_mask = np.logical_xor(sc_mask.astype(bool), sc_nucleus_mask.astype(bool))
        return sc_cytosol_mask

    @staticmethod
    def get_sc_junction_protein_mask(sc_membrane_mask, im_junction):
        single_cell_junction_protein = otsu_thresh_mask(sc_membrane_mask, im_junction)
        return single_cell_junction_protein.astype(bool)

    @staticmethod
    def _get_outline_from_mask(mask, width=1):
        """Computes outline for a mask with a single label"""
        mask = mask.astype(bool)
        dilated_mask = ndi.binary_dilation(mask, iterations=width)
        eroded_mask = ndi.binary_erosion(mask, iterations=width)
        outline_mask = np.logical_xor(dilated_mask, eroded_mask)

        return outline_mask
