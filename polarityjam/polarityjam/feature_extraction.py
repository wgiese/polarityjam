import numpy as np
import pandas as pd

from polarityjam.polarityjam_logging import get_logger
from polarityjam.utils import parameters
from polarityjam.utils.collector import fill_morans_i, Collector
from polarityjam.utils.masks import get_single_cell_membrane_mask, get_single_cell_mask, get_single_cell_organelle_mask, \
    get_single_cell_cytosol_mask, get_single_cell_nucleus_mask, get_organelle_mask, get_nuclei_mask, \
    get_single_junction_protein_mask, Masks, SingleCellMasks
from polarityjam.utils.moran import run_morans
from polarityjam.utils.neighborhood import k_neighbor_dif
from polarityjam.utils.plot import plot_dataset
from polarityjam.utils.properties import set_single_cell_nucleus_props, set_single_cell_organelle_props, \
    set_single_cell_props, set_single_cell_marker_nuclei_props, set_single_cell_marker_cytosol_props, \
    set_single_cell_marker_membrane_props, set_single_cell_marker_props, set_single_cell_junction_props, \
    SingleCellProperties
from polarityjam.utils.rag import orientation_graph_nf, remove_islands
from polarityjam.utils.seg import get_image_marker, get_image_junction


class Extractor:
    def __init__(self, img, cells_mask, filename, output_path):
        self.img = img
        self.img_marker = self.get_image_marker(self.img)
        self.img_junction = self.get_image_junction(self.img)
        self.img_nucleus = self.get_image_nucleus(self.img)
        self.img_organelle = self.get_image_organelle(self.img)
        self.filename = filename
        self.output_path = output_path
        self.collector = Collector()
        self.masks = Masks(cells_mask)

    @staticmethod
    def threshold(single_cell_mask, single_nucleus_mask=None, single_organelle_mask=None):
        """Thresholds given single_cell_mask. Returns True if falls under threshold."""
        # TODO: check if this can be removed, we already remove small objects from the cellpose mask
        if len(single_cell_mask[single_cell_mask == 1]) < parameters.min_cell_size:
            return True

        if single_nucleus_mask is not None:
            if len(single_nucleus_mask[single_nucleus_mask == 1]) < parameters.min_nucleus_size:
                return True

        if single_organelle_mask is not None:
            if len(single_organelle_mask[single_organelle_mask == 1]) < parameters.min_organelle_size:
                return True

        return False

    @staticmethod
    def get_image_marker(img):
        """Gets the image of the marker channel specified in the parameters."""
        if parameters.channel_expression_marker >= 0:
            return img[:, :, parameters.channel_expression_marker]
        return None

    @staticmethod
    def get_image_junction(img):
        """Gets the image of the junction channel specified in the parameters."""
        if parameters.channel_junction >= 0:
            return img[:, :, parameters.channel_junction]
        return None

    @staticmethod
    def get_image_nucleus(img):
        """Gets the image of the nucleus channel specified in the parameters."""
        if parameters.channel_nucleus >= 0:
            return img[:, :, parameters.channel_nucleus]
        return None

    @staticmethod
    def get_image_organelle(img):
        """Gets the image of the organelle channel specified in the parameters."""
        if parameters.channel_organelle >= 0:
            return img[:, :, parameters.channel_organelle]
        return None

    def extract(self, p):
        """ Extracts the features from an input image."""
        # initialize graph - no features associated with nodes
        rag = orientation_graph_nf(self.masks.cell_mask)
        list_of_islands = self.masks.set_cell_mask_rem_island(rag)
        rag = remove_islands(rag, list_of_islands)

        get_logger().info("Number of RAG nodes: %s " % len(list(rag.nodes)))

        self.masks.set_nuclei_mask(self.img_nucleus)
        self.masks.set_organelle_mask(self.img_nucleus)

        excluded = 0
        # iterate through each unique segmented cell
        for index, connected_component_label in enumerate(np.unique(self.masks.cell_mask_rem_island)):
            # ignore background
            if connected_component_label == 0:
                continue

            # get single cell masks
            sc_masks = SingleCellMasks(self.masks, connected_component_label)
            sc_masks.calc_sc_masks(self.img_junction)

            # threshold
            if self.threshold(
                    sc_masks.sc_mask,
                    single_nucleus_mask=sc_masks.sc_nucleus_mask,
                    single_organelle_mask=sc_masks.sc_organelle_mask
            ):
                get_logger().info("Cell \"%s\" falls under threshold! Removed from RAG..." % connected_component_label)
                excluded += 1
                # remove a cell from the RAG
                rag.remove_node(connected_component_label)
                continue

            sc_props = SingleCellProperties()
            sc_props.calc_sc_props(sc_masks, self.masks, self.img_marker, self.img_junction)

            self.collector.collect_sc_props(sc_props, index, self.filename, connected_component_label)

            # append feature of interest to the RAG node for being able to do further analysis
            foi_val = self.collector.dataset.at[index, parameters.feature_of_interest]
            rag.nodes[connected_component_label.astype('int')][parameters.feature_of_interest] = foi_val

            get_logger().info(
                " ".join(
                    str(x) for x in ["Cell %s - feature \"%s\": %s" % (
                        connected_component_label, parameters.feature_of_interest, foi_val
                    )]
                )
            )

        get_logger().info("Excluded cells: %s" % str(excluded))
        get_logger().info("Leftover cells: %s" % str(len(np.unique(self.masks.cell_mask)) - excluded))

        # morans I analysis based on FOI
        morans_i = run_morans(rag, parameters.feature_of_interest)
        self.collector.collect_morans_i_props(morans_i)

        # neighborhood feature analysis based on FOI
        neighborhood_props_list = k_neighbor_dif(rag, parameters.feature_of_interest)
        self.collector.collect_neighborhood_props(neighborhood_props_list)

        plot_dataset(
            p, self.img, self.collector.dataset, self.output_path, self.filename,
            self.masks.cell_mask_rem_island, self.masks.nuclei_mask,
            self.masks.organelle_mask,
            self.img_marker, self.img_junction
        )

        return self.collector.dataset


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
    im_junction = get_image_junction(parameters, img)

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
        single_membrane_mask = get_single_cell_membrane_mask(parameters, im_marker, im_junction,
                                                             single_cell_mask)  # can be None
        single_cytosol_mask = get_single_cell_cytosol_mask(single_cell_mask, im_marker,
                                                           single_nucleus_mask)  # can be None
        single_junction_protein_area_mask, single_junction_protein = get_single_junction_protein_mask(
            single_membrane_mask,
            im_junction
        )  # can be None

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
        set_single_cell_props(properties_dataset, index, filename, connected_component_label, single_cell_mask)

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

            # properties for marker nuclei
            if nuclei_mask is not None:
                set_single_cell_marker_nuclei_props(properties_dataset, index, single_nucleus_mask, im_marker)
                set_single_cell_marker_cytosol_props(properties_dataset, index, single_cytosol_mask, im_marker)

        # properties for junctions
        if im_junction is not None:
            set_single_cell_junction_props(
                properties_dataset,
                index,
                single_membrane_mask,
                im_junction,
                single_junction_protein,
                single_junction_protein_area_mask
            )

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

    # neighborhood feature analysis based on FOI
    k_neighbor_dif(properties_dataset, rag, parameters["feature_of_interest"])

    plot_dataset(
        parameters, img, properties_dataset, output_path, filename, cell_mask_rem_island, nuclei_mask, organelle_mask,
        im_marker, im_junction
    )

    return properties_dataset
