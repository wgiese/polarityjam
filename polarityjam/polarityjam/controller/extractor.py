import numpy as np

from polarityjam.controller.collector import PropertyCollector, SingleCellPropertyCollector
from polarityjam.model.masks import MasksCollection, SingleCellMasksCollection
from polarityjam.model.properties import SingleCellPropertiesCollection
from polarityjam.polarityjam_logging import get_logger
from polarityjam.utils import parameters
from polarityjam.compute.moran import run_morans
from polarityjam.compute.neighborhood import k_neighbor_dif
from polarityjam.vizualization.plot import plot_dataset
from polarityjam.utils.rag import orientation_graph_nf, remove_islands


class Extractor:
    def __init__(self, img, cells_mask, filename, output_path):
        self.img = img
        self.img_marker = self.get_image_marker(self.img)
        self.img_junction = self.get_image_junction(self.img)
        self.img_nucleus = self.get_image_nucleus(self.img)
        self.img_organelle = self.get_image_organelle(self.img)
        self.filename = filename
        self.output_path = output_path
        self.collector = PropertyCollector()
        self.masks = MasksCollection(cells_mask)

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

    def extract(self):
        """ Extracts the features from an input image."""
        # initialize graph - no features associated with nodes
        rag = orientation_graph_nf(self.masks.cell_mask)
        list_of_islands = self.masks.set_cell_mask_rem_island(rag)
        rag = remove_islands(rag, list_of_islands)

        get_logger().info("Number of RAG nodes: %s " % len(list(rag.nodes)))

        self.masks.set_nuclei_mask(self.img_nucleus)
        self.masks.set_organelle_mask(self.img_organelle)

        excluded = 0
        # iterate through each unique segmented cell
        for index, connected_component_label in enumerate(np.unique(self.masks.cell_mask_rem_island)):
            index -= excluded  # index should be continuous

            # ignore background
            if connected_component_label == 0:
                continue

            # get single cell masks
            sc_masks = SingleCellMasksCollection(self.masks, connected_component_label)
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

            sc_props_collection = SingleCellPropertyCollector().calc_sc_props(
                sc_masks, self.img_marker, self.img_junction
            )

            self.collector.collect_sc_props(sc_props_collection, index, self.filename, connected_component_label)

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

        plot_dataset(self.collector.dataset, self.masks.cell_mask_rem_island, self.masks.nuclei_mask,
                     self.masks.organelle_mask, self.img_marker, self.img_junction, self.filename, self.output_path)

        return self.collector.dataset
