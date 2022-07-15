import numpy as np
import skimage.filters
from scipy import ndimage as ndi

from polarityjam.polarityjam_logging import get_logger
from polarityjam.utils import parameters
from polarityjam.compute.compute import otsu_thresh_mask


class SingleCellMasksCollection:

    def __init__(self, masks, connected_component_label):
        self.connected_component_label = connected_component_label
        self.masks = masks

        self.sc_mask = None
        self.sc_nucleus_mask = None
        self.sc_organelle_mask = None
        self.sc_membrane_mask = None
        self.sc_cytosol_mask = None
        self.sc_junction_protein_area_mask = None
        self.sc_junction_protein = None  # todo: put me elsewhere

    def calc_sc_masks(self, im_junction):  # todo: make me a controller
        self.get_sc_mask()
        self.get_sc_nucleus_mask()  # can be None
        self.get_sc_organelle_mask()  # can be None
        self.get_sc_membrane_mask()  # can be None
        self.get_sc_cytosol_mask()  # can be None
        self.sc_junction_protein_area_mask, self.sc_junction_protein = self.get_sc_junction_protein_mask(
            im_junction)  # can be None

    def get_sc_organelle_mask(self):
        """Gets the single cell organelle mask."""
        sc_organelle_mask = None
        if self.masks.organelle_mask is not None:
            sc_organelle_mask = np.where(self.masks.organelle_mask == self.connected_component_label, 1, 0)

        self.sc_organelle_mask = sc_organelle_mask

    def get_sc_nucleus_mask(self):
        """Gets the single cell nucleus mask."""
        sc_nucleus_mask = None
        if self.masks.nuclei_mask is not None:
            sc_nucleus_mask = np.where(self.masks.nuclei_mask == self.connected_component_label, 1, 0)
        self.sc_nucleus_mask = sc_nucleus_mask

    def get_sc_mask(self):
        """Gets the single cell mask from a mask where each cell has an increasing connected component value."""
        sc_mask = None
        if self.masks.cell_mask_rem_island is not None:
            sc_mask = np.where(self.masks.cell_mask_rem_island == self.connected_component_label, 1, 0)
        self.sc_mask = sc_mask

    def get_sc_membrane_mask(self):
        """Gets the single cell membrane mask."""
        if not self.sc_mask is not None:
            self.get_sc_mask()
        sc_membrane_mask = self._get_outline_from_mask(self.sc_mask, parameters.membrane_thickness)
        self.sc_membrane_mask = sc_membrane_mask

    def get_sc_cytosol_mask(self):
        """Gets the cytosol mask."""
        sc_cytosol_mask = None
        if self.sc_nucleus_mask is not None:
            sc_cytosol_mask = np.logical_xor(self.sc_mask.astype(bool), self.sc_nucleus_mask.astype(bool))
        self.sc_cytosol_mask = sc_cytosol_mask

    def get_sc_junction_protein_mask(self, im_junction):
        if im_junction is not None:
            single_cell_junction_protein = otsu_thresh_mask(self.sc_membrane_mask, im_junction)
            return single_cell_junction_protein.astype(bool), single_cell_junction_protein
        return None, None

    @staticmethod
    def _get_outline_from_mask(mask, width=1):
        """Computes outline for a mask with a single label"""
        mask = mask.astype(bool)
        dilated_mask = ndi.binary_dilation(mask, iterations=width)
        eroded_mask = ndi.binary_erosion(mask, iterations=width)
        outline_mask = np.logical_xor(dilated_mask, eroded_mask)

        return outline_mask


class MasksCollection:
    def __init__(self, cell_mask):
        self.cell_mask = cell_mask
        self.cell_mask_rem_island = None
        self.nuclei_mask = None
        self.organelle_mask = None

    def set_cell_mask_rem_island(self, neighborhood_graph):
        """Remove unconnected cells (Cells without neighbours)."""
        # Get list of islands - nodes with no neighbours and remove them
        list_of_islands = []
        for nodes in neighborhood_graph.nodes:
            if len(list(neighborhood_graph.neighbors(nodes))) == 0:
                list_of_islands.append(nodes)

        list_of_islands = np.unique(list_of_islands)

        mask = self.cell_mask

        # remove islands from mask
        for elemet in list_of_islands:
            mask[:, :][mask[:, :] == elemet] = 0

        self.cell_mask_rem_island = mask

        get_logger().info("Removed number of islands: %s" % len(list_of_islands))

        return list_of_islands

    def set_nuclei_mask(self, img_nuclei):
        """Sets the nuclei mask."""
        nuclei_mask = None
        if parameters.channel_nucleus >= 0:
            img_nuclei_blur = ndi.gaussian_filter(img_nuclei, sigma=3)
            nuclei_mask = np.where(img_nuclei_blur > skimage.filters.threshold_otsu(img_nuclei_blur), True, False)
            if self.cell_mask_rem_island is not None:
                nuclei_mask = nuclei_mask * self.cell_mask_rem_island
            self.nuclei_mask = nuclei_mask

        return nuclei_mask

    def set_organelle_mask(self, img_organelle):
        """Set the organelle mask."""
        if self.cell_mask_rem_island is None:
            return None

        organelle_mask = None
        if parameters.channel_organelle >= 0:
            img_organelle_blur = ndi.gaussian_filter(img_organelle, sigma=3)
            organelle_mask_o = np.where(
                img_organelle_blur > skimage.filters.threshold_otsu(img_organelle_blur), True, False
            )
            organelle_mask = organelle_mask_o * self.cell_mask_rem_island
            self.organelle_mask = organelle_mask

        return organelle_mask


def get_single_cell_mask(connected_component_label, cell_mask):
    """Gets the single cell mask from a mask where each cell has an increasing connected component value."""
    return np.where(cell_mask == connected_component_label, 1, 0)  # convert connected_component_label to 1


def get_outline_from_mask(mask, width=1):
    """Computes outline for a mask with a single label"""

    mask = mask.astype(bool)
    dilated_mask = ndi.binary_dilation(mask, iterations=width)
    eroded_mask = ndi.binary_erosion(mask, iterations=width)
    outline_mask = np.logical_xor(dilated_mask, eroded_mask)

    return outline_mask
