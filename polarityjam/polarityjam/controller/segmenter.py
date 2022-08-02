from pathlib import Path

import cellpose.models
import numpy as np
import skimage
from skimage import morphology

from polarityjam.model.parameter import SegmentationParameter
from polarityjam.polarityjam_logging import get_logger


class Segmenter:

    def __init__(self, params: SegmentationParameter):
        self.params = params

    def get_cellpose_model(self):
        """Gets the specified cellpose model"""

        if self.params.cp_model_type == "custom":
            model = cellpose.models.CellposeModel(gpu=self.params.use_gpu, pretrained_model=self.params.cp_model_path)
        else:
            model = cellpose.models.Cellpose(gpu=self.params.use_gpu, model_type=self.params.cp_model_type)

        return model

    def get_cellpose_segmentation(self, im_seg, filepath):
        """Gets the cellpose segmentation. Expects im_seg to have junction channel first, then nucleus channel."""
        get_logger().info("Calculate cellpose segmentation. This might take some time...")

        model = self.get_cellpose_model()
        if im_seg.ndim > 1:
            channels = [1, 2]
        else:
            channels = [0, 0]

        # masks, flows, styles, diams = model.eval(im_seg, channels=channels)

        if self.params.cp_model_type == "custom":
            masks, flows, styles = model.eval(im_seg, diameter=self.params.estimated_cell_diameter, channels=channels)
        else:
            masks, flows, styles, diams = model.eval(im_seg, diameter=self.params.estimated_cell_diameter,
                                                     channels=channels)

        if self.params.store_segmentation:
            segmentation_list = {"masks": masks}
            segmentation, _ = self.get_segmentation_file_name(filepath)

            get_logger().info("Storing cellpose segmentation: %s" % segmentation)
            np.save(str(segmentation), segmentation_list, allow_pickle=True)

        return masks

    def get_segmentation_file_name(self, filepath):
        stem = Path(filepath).stem

        suffix = "_seg.npy"
        if self.params.manually_annotated_mask:
            suffix = self.params.manually_annotated_mask
        segmentation = Path(filepath).parent.joinpath(stem + suffix)

        return segmentation, stem

    def load_or_get_cellpose_segmentation(self, img_seg, filepath):
        get_logger().info("Look up cellpose segmentation...")
        segmentation, _ = self.get_segmentation_file_name(filepath)

        if segmentation.exists() and self.params.use_given_mask:
            get_logger().info("Load cellpose segmentation...")

            # in case an annotated mask is available
            cellpose_seg = np.load(str(segmentation), allow_pickle=True)
            cellpose_mask = cellpose_seg.item()['masks']

        else:
            cellpose_mask = self.get_cellpose_segmentation(img_seg, filepath)

        if self.params.clear_border:
            cellpose_mask_clear_border = skimage.segmentation.clear_border(cellpose_mask)
            number_of_cellpose_borders = len(np.unique(cellpose_mask)) - len(np.unique(cellpose_mask_clear_border))
            cellpose_mask = cellpose_mask_clear_border

            get_logger().info("Removed number of cellpose borders: %s" % number_of_cellpose_borders)

            # TODO: remove small objects here
            cellpose_mask_remove_small_objects = morphology.remove_small_objects(
                cellpose_mask, self.params.min_cell_size, connectivity=2
            )
            number_of_cellpose_small_objects = len(np.unique(cellpose_mask)) - len(
                np.unique(cellpose_mask_remove_small_objects))
            cellpose_mask = cellpose_mask_remove_small_objects

            get_logger().info("Removed number of small objects: %s" % number_of_cellpose_small_objects)

        get_logger().info("Detected number of cellpose labels: %s" % len(np.unique(cellpose_mask)))

        return cellpose_mask
