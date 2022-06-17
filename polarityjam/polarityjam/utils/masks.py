import numpy as np
import skimage.filters
from cellpose import utils
from scipy import ndimage as ndi


def get_single_cell_membrane_mask(parameters, im_marker, single_cell_mask):
    """Gets the single cell membrane mask."""
    if im_marker is not None:
        return get_outline_from_mask(single_cell_mask.astype(bool), parameters["membrane_thickness"])
    else:
        return None


def get_single_cell_mask(connected_component_label, cell_mask):
    """Gets the single cell mask from a mask where each cell has an increasing connected component value."""
    return np.where(cell_mask == connected_component_label, 1, 0)  # convert connected_component_label to 1


def get_single_cell_organelle_mask(connected_component_label, organelle_mask):
    """Gets the single cell organelle mask."""
    if organelle_mask is not None:
        return np.where(organelle_mask == connected_component_label, 1, 0)
    else:
        return None


def get_single_cell_cytosol_mask(single_cell_mask, im_marker, single_nucleus_mask):
    """Gets the cytosol mask."""
    if single_nucleus_mask is not None and im_marker is not None:
        return np.logical_xor(single_cell_mask.astype(bool), single_nucleus_mask.astype(bool))
    else:
        return None


def get_single_cell_nucleus_mask(connected_component_label, nuclei_mask):
    """Gets the single cell nucleus mask."""
    if nuclei_mask is not None:
        return np.where(nuclei_mask == connected_component_label, 1, 0)
    else:
        return None


def get_organelle_mask(parameters, img, cellpose_mask):
    """Gets the organelle mask."""
    if parameters["channel_organelle"] >= 0:
        img_organelle_blur = ndi.gaussian_filter(img[:, :, parameters["channel_organelle"]], sigma=3)

        organelle_mask = np.where(img_organelle_blur > skimage.filters.threshold_otsu(img_organelle_blur), True, False)

        organelle_label = organelle_mask * cellpose_mask

        return organelle_label
    else:
        return None


def get_nuclei_mask(parameters, img, cellpose_mask):
    """Gets the nuclei mask."""
    if parameters["channel_nucleus"] >= 0:
        img_nuclei_blur = ndi.gaussian_filter(img[:, :, parameters["channel_nucleus"]], sigma=3)

        nuclei_mask = np.where(img_nuclei_blur > skimage.filters.threshold_otsu(img_nuclei_blur), True, False)

        nuclei_label = nuclei_mask * cellpose_mask

        return nuclei_label
    else:
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
