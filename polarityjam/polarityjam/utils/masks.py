import numpy as np
import skimage.filters
from cellpose import utils
from scipy import ndimage as ndi

from polarityjam.polarityjam_logging import get_logger
from polarityjam.utils import parameters


class SingleCellMasks:

    def __init__(self):
        self.sc_mask = None
        self.sc_nucleus_mask = None
        self.sc_organelle_mask = None
        self.sc_membrane_mask = None
        self.sc_cytosol_mask = None
        self.sc_junction_protein_area_mask = None
        self.sc_junction_protein = None  # todo: put me elsewhere

    def calc_sc_masks(self, masks, connected_component_label, im_marker, im_junction):
        self.sc_mask = masks.get_sc_mask(connected_component_label)
        self.sc_nucleus_mask = masks.get_sc_nucleus_mask(connected_component_label)  # can be None
        self.sc_organelle_mask = masks.get_sc_organelle_mask(connected_component_label)  # can be None
        self.sc_membrane_mask = masks.get_sc_membrane_mask(im_marker, im_junction,
                                                           self.sc_mask)  # can be None
        self.sc_cytosol_mask = masks.get_sc_cytosol_mask(self.sc_mask, im_marker,
                                                         self.sc_nucleus_mask)  # can be None
        self.sc_junction_protein_area_mask, self.sc_junction_protein = masks.get_sc_junction_protein_mask(
            self.sc_membrane_mask,
            im_junction
        )  # can be None


class Masks:
    def __init__(self, cell_mask):
        self.cell_mask = cell_mask
        self.cell_mask_rem_island = None
        self.nuclei_mask = None
        self.organelle_mask = None

    def remove_islands(self, frame_graph):
        """Remove unconnected cells (Cells without neighbours)."""
        # Get list of islands - nodes with no neighbours and remove them
        list_of_islands = []
        for nodes in frame_graph.nodes:
            if len(list(frame_graph.neighbors(nodes))) == 0:
                list_of_islands.append(nodes)

        get_logger().info("Removed number of islands: %s" % len(list_of_islands))

        mask = self.cell_mask

        # remove islands from image and graph
        for elemet in np.unique(list_of_islands):
            frame_graph.remove_node(elemet)
            mask[:, :][mask[:, :] == elemet] = 0

        self.cell_mask_rem_island = mask

        return frame_graph

    def set_nuclei_mask(self, img):
        """Gets the nuclei mask."""
        if parameters.channel_nucleus >= 0:
            img_nuclei_blur = ndi.gaussian_filter(img[:, :, parameters.channel_nucleus], sigma=3)

            nuclei_mask = np.where(img_nuclei_blur > skimage.filters.threshold_otsu(img_nuclei_blur), True, False)

            nuclei_label = nuclei_mask * self.cell_mask_rem_island

            self.nuclei_mask = nuclei_label
            return True
        return False

    def set_organelle_mask(self, img):
        """Gets the organelle mask."""
        if parameters.channel_organelle >= 0:
            img_organelle_blur = ndi.gaussian_filter(img[:, :, parameters.channel_organelle], sigma=3)

            organelle_mask = np.where(img_organelle_blur > skimage.filters.threshold_otsu(img_organelle_blur), True,
                                      False)
            organelle_label = organelle_mask * self.cell_mask_rem_island

            self.organelle_mask = organelle_label
            return True
        return False

    def get_sc_mask(self, connected_component_label):
        """Gets the single cell mask from a mask where each cell has an increasing connected component value."""
        return np.where(self.cell_mask_rem_island == connected_component_label, 1,
                        0)  # convert connected_component_label to 1

    def get_sc_organelle_mask(self, connected_component_label):
        """Gets the single cell organelle mask."""
        if self.organelle_mask is not None:
            return np.where(self.organelle_mask == connected_component_label, 1, 0)
        return None

    def get_sc_nucleus_mask(self, connected_component_label):
        """Gets the single cell nucleus mask."""
        if self.nuclei_mask is not None:
            return np.where(self.nuclei_mask == connected_component_label, 1, 0)
        return None

    def get_sc_membrane_mask(self, im_marker, im_junctions, single_cell_mask):
        """Gets the single cell membrane mask."""
        if im_marker is not None or im_junctions is not None:
            return get_outline_from_mask(single_cell_mask.astype(bool), parameters.membrane_thickness)
        return None

    def get_sc_cytosol_mask(self, single_cell_mask, im_marker, single_nucleus_mask):
        """Gets the cytosol mask."""
        if single_nucleus_mask is not None and im_marker is not None:
            return np.logical_xor(single_cell_mask.astype(bool), single_nucleus_mask.astype(bool))
        return None

    def get_sc_junction_protein_mask(self, single_membrane_mask, im_junction):
        if im_junction is not None:
            single_cell_junction_protein = get_single_junction_protein(single_membrane_mask, im_junction)
            return single_cell_junction_protein.astype(bool), single_cell_junction_protein
        return None, None

    # todo: move me elsewhere. But avoid circular import error (don't move to feature_extraction.py)
    def get_sc_junction_protein(self, single_membrane_mask, im_junction):
        # otsu threshold membrane (junctions) intensity to get junction protein area
        single_cell_junction_protein = single_membrane_mask * im_junction
        otsu_val = skimage.filters.threshold_otsu(single_cell_junction_protein)
        single_cell_junction_protein[single_cell_junction_protein <= otsu_val] = 0

        return single_cell_junction_protein


def get_single_cell_membrane_mask(parameters, im_marker, im_junctions, single_cell_mask):
    """Gets the single cell membrane mask."""
    if im_marker is not None or im_junctions is not None:
        return get_outline_from_mask(single_cell_mask.astype(bool), parameters["membrane_thickness"])
    return None


def get_single_cell_mask(connected_component_label, cell_mask):
    """Gets the single cell mask from a mask where each cell has an increasing connected component value."""
    return np.where(cell_mask == connected_component_label, 1, 0)  # convert connected_component_label to 1


def get_single_cell_organelle_mask(connected_component_label, organelle_mask):
    """Gets the single cell organelle mask."""
    if organelle_mask is not None:
        return np.where(organelle_mask == connected_component_label, 1, 0)
    return None


def get_single_cell_cytosol_mask(single_cell_mask, im_marker, single_nucleus_mask):
    """Gets the cytosol mask."""
    if single_nucleus_mask is not None and im_marker is not None:
        return np.logical_xor(single_cell_mask.astype(bool), single_nucleus_mask.astype(bool))
    return None


def get_single_cell_nucleus_mask(connected_component_label, nuclei_mask):
    """Gets the single cell nucleus mask."""
    if nuclei_mask is not None:
        return np.where(nuclei_mask == connected_component_label, 1, 0)
    return None


def get_organelle_mask(parameters, img, cellpose_mask):
    """Gets the organelle mask."""
    if parameters["channel_organelle"] >= 0:
        img_organelle_blur = ndi.gaussian_filter(img[:, :, parameters["channel_organelle"]], sigma=3)

        organelle_mask = np.where(img_organelle_blur > skimage.filters.threshold_otsu(img_organelle_blur), True, False)
        organelle_label = organelle_mask * cellpose_mask

        return organelle_label
    return None


def get_single_junction_protein_mask(single_membrane_mask, im_junction):
    if im_junction is not None:
        single_cell_junction_protein = get_single_junction_protein(single_membrane_mask, im_junction)
        return single_cell_junction_protein.astype(bool), single_cell_junction_protein
    return None, None


# todo: move me elsewhere. But avoid circular import error (don't move to feature_extraction.py)
def get_single_junction_protein(single_membrane_mask, im_junction):
    # otsu threshold membrane (junctions) intensity to get junction protein area
    single_cell_junction_protein = single_membrane_mask * im_junction
    otsu_val = skimage.filters.threshold_otsu(single_cell_junction_protein)
    single_cell_junction_protein[single_cell_junction_protein <= otsu_val] = 0

    return single_cell_junction_protein


def get_nuclei_mask(parameters, img, cellpose_mask):
    """Gets the nuclei mask."""
    if parameters["channel_nucleus"] >= 0:
        img_nuclei_blur = ndi.gaussian_filter(img[:, :, parameters["channel_nucleus"]], sigma=3)

        nuclei_mask = np.where(img_nuclei_blur > skimage.filters.threshold_otsu(img_nuclei_blur), True, False)

        nuclei_label = nuclei_mask * cellpose_mask

        return nuclei_label
    return None


def get_outline_from_mask(mask, width=1):
    """Computes outline for a mask with a single label"""

    mask = mask.astype(bool)
    dilated_mask = ndi.binary_dilation(mask, iterations=width)
    eroded_mask = ndi.binary_erosion(mask, iterations=width)
    outline_mask = np.logical_xor(dilated_mask, eroded_mask)

    return outline_mask


def get_outline_with_multiple_labels(mask, width=1):
    """ TODO: not used at the moment, faster for multiple masks/mask labels"""

    outline_list = np.array(utils.outlines_list(mask))
    outlines = np.zeros((mask.shape[0], mask.shape[1]))
    for mask_id, outline_coords in enumerate(outline_list):
        if outline_coords.T.shape[1] < 10:
            outlines[tuple(outline_coords.T)] = mask_id + 1

    outlines_mask = ndi.morphology.binary_dilation(outlines.astype(bool), iterations=width)

    return outlines_mask
