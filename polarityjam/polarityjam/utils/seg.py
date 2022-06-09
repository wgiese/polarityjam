import pickle
from pathlib import Path

import cellpose.models
import numpy as np
import skimage.segmentation
from cellpose import utils
from scipy import ndimage as ndi
from skimage import morphology

from polarityjam.polarityjam_logging import get_logger


def get_cellpose_model(parameters):
    """Gets the specified cellpose model"""

    if parameters["cp_model_type"] == "custom":
        model = cellpose.models.CellposeModel(gpu=parameters["use_gpu"], pretrained_model=parameters["cp_model_path"])
    else:
        model = cellpose.models.Cellpose(gpu=parameters["use_gpu"], model_type=parameters["cp_model_type"])

    return model


def get_cellpose_segmentation(parameters, im_seg, filepath):
    """Gets the cellpose segmentation"""
    get_logger().info("Calculate cellpose segmentation. This might take some time...")

    model = get_cellpose_model(parameters)
    if parameters["channel_nucleus"] >= 0:
        channels = [1, 2]
    else:
        channels = [0, 0]

    # masks, flows, styles, diams = model.eval(im_seg, channels=channels)

    if parameters["cp_model_type"] == "custom":
        masks, flows, styles = model.eval(im_seg, diameter=parameters["estimated_cell_diameter"], channels=channels)
    else:
        masks, flows, styles, diams = model.eval(im_seg, diameter=parameters["estimated_cell_diameter"],
                                                 channels=channels)

    if parameters["store_segmentation"]:
        segmentation_list = {"masks": masks}
        segmentation, _ = get_segmentation_file_name(parameters, filepath)

        get_logger().info("Storing cellpose segmentation: %s" % segmentation)
        np.save(segmentation, segmentation_list, allow_pickle=True)

    return masks


def get_segmentation_file_name(parameters, filepath):
    stem = Path(filepath).stem

    suffix = "_seg.npy"
    if parameters["manually_annotated_mask"]:
        suffix = parameters["manually_annotated_mask"]
    segmentation = Path(filepath).parent.joinpath(stem + suffix)

    return segmentation, stem


def load_or_get_cellpose_segmentation(parameters, img_seg, filepath):
    get_logger().info("Look up cellpose segmentation...")
    segmentation, _ = get_segmentation_file_name(parameters, filepath)

    if segmentation.exists() and parameters["use_given_mask"]:
        get_logger().info("Load cellpose segmentation...")

        # in case an annotated mask is available
        cellpose_seg = np.load(str(segmentation), allow_pickle=True)
        cellpose_mask = cellpose_seg.item()['masks']

    else:
        cellpose_mask = get_cellpose_segmentation(parameters, img_seg, filepath)

    if parameters["clear_border"]:
        cellpose_mask_clear_border = skimage.segmentation.clear_border(cellpose_mask)
        number_of_cellpose_borders = len(np.unique(cellpose_mask)) - len(np.unique(cellpose_mask_clear_border))
        cellpose_mask = cellpose_mask_clear_border

        get_logger().info("Removed number of cellpose borders: %s" % number_of_cellpose_borders)

        # TODO: remove small objects here
        cellpose_mask_remove_small_objects = morphology.remove_small_objects(
            cellpose_mask, parameters["min_cell_size"], connectivity=2
        )
        number_of_cellpose_small_objects = len(np.unique(cellpose_mask)) - len(
            np.unique(cellpose_mask_remove_small_objects))
        cellpose_mask = cellpose_mask_remove_small_objects

        get_logger().info("Removed number of small objects: %s" % number_of_cellpose_small_objects)

    get_logger().info("Detected number of cellpose labels: %s" % len(np.unique(cellpose_mask)))

    return cellpose_mask


def get_outline_from_mask(mask, width=1):
    """computes outline for a mask with a single label"""

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
