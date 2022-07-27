from polarityjam.model.collection import PropertiesCollection
from polarityjam.polarityjam_logging import get_logger
from polarityjam.utils import parameters
from polarityjam.vizualization.plot import plot_organelle_polarity, plot_nuc_displacement_orientation, \
    plot_marker_expression, plot_marker_polarity, plot_marker_nucleus_orientation, plot_junction_polarity, plot_corners, \
    plot_eccentricity, plot_ratio_method, plot_orientation


class Plotter:

    def __init__(self):
        pass

    def plot_organelle_polarity(self, collection, img_name):
        plot_organelle_polarity(
            collection.img_channel_dict[img_name]["junction"],
            collection.masks_dict[img_name].cell_mask_rem_island,
            collection.masks_dict[img_name].nuclei_mask,
            collection.masks_dict[img_name].organelle_mask,
            collection.dataset.loc[collection.dataset["filename"] == img_name],
            img_name,
            collection.out_path_dict[img_name]
        )

    def plot_nuc_displacement_orientation(self, collection, img_name):
        plot_nuc_displacement_orientation(
            collection.img_channel_dict[img_name]["junction"],
            collection.masks_dict[img_name].cell_mask_rem_island,
            collection.masks_dict[img_name].nuclei_mask,
            collection.dataset.loc[collection.dataset["filename"] == img_name],
            img_name,
            collection.out_path_dict[img_name]
        )

    def plot_marker_expression(self, collection, img_name):
        plot_marker_expression(
            collection.img_channel_dict[img_name]["marker"],
            collection.masks_dict[img_name].cell_mask_rem_island,
            collection.dataset.loc[collection.dataset["filename"] == img_name],
            img_name,
            collection.out_path_dict[img_name],
            nuclei_mask=collection.masks_dict[img_name].nuclei_mask
        )

    def plot_marker_polarity(self, collection, img_name):
        plot_marker_polarity(
            collection.img_channel_dict[img_name]["marker"],
            collection.masks_dict[img_name].cell_mask_rem_island,
            collection.dataset.loc[collection.dataset["filename"] == img_name],
            img_name,
            collection.out_path_dict[img_name]
        )

    def plot_marker_nucleus_orientation(self, collection, img_name):
        plot_marker_nucleus_orientation(
            collection.img_channel_dict[img_name]["junction"],
            collection.masks_dict[img_name].cell_mask_rem_island,
            collection.masks_dict[img_name].nuclei_mask,
            collection.dataset.loc[collection.dataset["filename"] == img_name],
            img_name,
            collection.out_path_dict[img_name]
        )

    def plot_junction_polarity(self, collection, img_name):
        plot_junction_polarity(
            collection.img_channel_dict[img_name]["junction"],
            collection.masks_dict[img_name].cell_mask_rem_island,
            collection.dataset.loc[collection.dataset["filename"] == img_name],
            img_name,
            collection.out_path_dict[img_name]
        )

    def plot_corners(self, collection, img_name):
        plot_corners(
            collection.img_channel_dict[img_name]["junction"],
            collection.dataset.loc[collection.dataset["filename"] == img_name],
            img_name,
            collection.out_path_dict[img_name]
        )

    def plot_eccentricity(self, collection, img_name):
        plot_eccentricity(
            collection.img_channel_dict[img_name]["junction"],
            collection.dataset.loc[collection.dataset["filename"] == img_name],
            collection.masks_dict[img_name].cell_mask_rem_island,
            img_name,
            collection.out_path_dict[img_name],
            nuclei_mask=collection.masks_dict[img_name].nuclei_mask
        )

    def plot_ratio_method(self, collection, img_name):
        plot_ratio_method(
            collection.img_channel_dict[img_name]["junction"],
            collection.masks_dict[img_name].cell_mask_rem_island,
            collection.dataset.loc[collection.dataset["filename"] == img_name],
            img_name,
            collection.out_path_dict[img_name]
        )

    def plot_orientation(self, collection, img_name):
        plot_orientation(
            collection.img_channel_dict[img_name]["junction"],
            collection.dataset.loc[collection.dataset["filename"] == img_name],
            img_name,
            collection.out_path_dict[img_name],
            collection.masks_dict[img_name].cell_mask_rem_island,
            nuclei_mask=collection.masks_dict[img_name].nuclei_mask
        )

    def plot_collection(self, collection: PropertiesCollection):
        """Plots the properties dataset"""
        get_logger().info("Plotting...")

        for key in collection.img_channel_dict.keys():

            nuclei_mask = collection.masks_dict[key].nuclei_mask
            organelle_mask = collection.masks_dict[key].organelle_mask
            img_marker = collection.img_channel_dict[key]["marker"]
            img_junction = collection.img_channel_dict[key]["junction"]

            if parameters.plot_polarity and nuclei_mask is not None and organelle_mask is not None:
                self.plot_organelle_polarity(collection, key)
                if nuclei_mask is not None:
                    self.plot_nuc_displacement_orientation(collection, key)

            if parameters.plot_marker and img_marker is not None:
                self.plot_marker_expression(collection, key)
                self.plot_marker_polarity(collection, key)
                if nuclei_mask is not None:
                    self.plot_marker_nucleus_orientation(collection, key)

            if parameters.plot_junctions and img_junction is not None:
                self.plot_junction_polarity(collection, key)
                self.plot_corners(collection, key)

            if parameters.plot_orientation:
                self.plot_eccentricity(collection, key)

            if parameters.plot_ratio_method:
                self.plot_ratio_method(collection, key)

            if parameters.plot_cyclic_orientation:
                self.plot_orientation(collection, key)
