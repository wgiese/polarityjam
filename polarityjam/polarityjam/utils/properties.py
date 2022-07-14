import skimage.measure

from polarityjam.utils.compute import map_single_cell_to_circle
from polarityjam.utils.collector import fill_single_cell_marker_polarity, fill_single_nucleus_data_frame, \
    fill_single_cell_general_data_frame, fill_single_cell_organelle_data_frame, \
    fill_single_cell_marker_nuclei_data_frame, fill_single_cell_marker_nuclei_cytosol_data_frame, \
    fill_single_cell_marker_membrane_data_frame, fill_single_cell_junction_data_frame, \
    fill_single_cell_junction_sec_stat_data_frame
from polarityjam.utils.masks import SingleCellMasks, Masks


class SCJunctionProps:
    def __init__(
            self, interface_props,
            circular_junction_props,
            interface_perimeter,
            junction_protein_area_props,
            junction_fragmented_perimeter,
            junction_interface_occupancy,
            junction_protein_intensity,
            junction_intensity_per_interface_area,
            junction_cluster_density
    ):
        self.interface_props = interface_props
        self.circular_junction_props = circular_junction_props
        self.interface_perimeter = interface_perimeter
        self.junction_protein_area_props = junction_protein_area_props
        self.junction_fragmented_perimeter = junction_fragmented_perimeter
        self.junction_interface_occupancy = junction_interface_occupancy
        self.junction_protein_intensity = junction_protein_intensity
        self.junction_intensity_per_interface_area = junction_intensity_per_interface_area
        self.junction_cluster_density = junction_cluster_density


class NeighborhoodProps:
    def __init__(self):
        self.num_neighbours = None
        self.mean_dif_first_neighbors = None
        self.median_dif_first_neighbors = None
        self.var_dif_first_neighbors = None
        self.range_dif_first_neighbors = None
        self.mean_dif_second_neighbors = None
        self.median_dif_second_neighbors = None
        self.var_dif_second_neighbors = None
        self.range_dif_second_neighbors = None


class SingleCellProperties:

    def __init__(self):
        self.nucleus_props = None
        self.organelle_props = None
        self.single_cell_props = None
        self.marker_nuc_props = None
        self.marker_nuc_cyt_props = None
        self.marker_membrane_props = None
        self.marker_props = None
        self.junction_props = None

    def calc_sc_props(self, sc_masks: SingleCellMasks, masks: Masks, im_marker, im_junction):
        """calculates all properties for the single cell"""

        # properties for single cell
        self.set_single_cell_props(sc_masks.sc_mask)

        # properties for nucleus:
        if sc_masks.sc_nucleus_mask is not None:
            self.set_single_cell_nucleus_props(sc_masks.sc_nucleus_mask)

            # properties for organelle
            if self.nucleus_props and sc_masks.sc_organelle_mask is not None:
                self.set_single_cell_organelle_props(sc_masks.sc_organelle_mask)

        # properties for marker
        if im_marker is not None:
            self.set_single_cell_marker_props(sc_masks.sc_mask, im_marker)
            self.set_single_cell_marker_membrane_props(sc_masks.sc_membrane_mask, im_marker)

            # properties for marker nuclei
            if masks.nuclei_mask is not None:
                self.set_single_cell_marker_nuclei_props(sc_masks.sc_nucleus_mask, im_marker)
                self.set_single_cell_marker_cytosol_props(sc_masks.sc_cytosol_mask, im_marker)

        # properties for junctions
        if im_junction is not None:
            self.set_single_cell_junction_props(
                sc_masks.sc_membrane_mask,
                im_junction,
                sc_masks.sc_junction_protein,
                sc_masks.sc_junction_protein_area_mask,
                2,  # todo: replace me
                2  # todo: replace me
            )

    def set_single_cell_props(self, single_cell_mask):
        self.single_cell_props = get_single_cell_prop(single_cell_mask)

    def set_single_cell_nucleus_props(self, single_nucleus_mask):
        regions = skimage.measure.regionprops(single_nucleus_mask)
        nucleus_props = None
        if regions:
            nucleus_props = regions[-1]

        self.nucleus_props = nucleus_props

    def set_single_cell_organelle_props(self, single_organelle_mask):
        regions = skimage.measure.regionprops(single_organelle_mask)
        organelle_props = None
        if regions:
            organelle_props = regions[-1]

        self.organelle_props = organelle_props

    def set_single_cell_marker_props(self, single_cell_mask, im_marker):
        regions = skimage.measure.regionprops(single_cell_mask, intensity_image=im_marker)
        marker_props = None
        if regions:
            marker_props = regions[-1]

        self.marker_props = marker_props

    def set_single_cell_marker_membrane_props(self, single_membrane_mask, im_marker):
        regions = skimage.measure.regionprops(single_membrane_mask.astype(int), intensity_image=im_marker)
        marker_membrane_props = None
        if regions:
            marker_membrane_props = regions[-1]

        self.marker_membrane_props = marker_membrane_props

    def set_single_cell_marker_nuclei_props(self, single_nucleus_mask, im_marker):
        regions = skimage.measure.regionprops(single_nucleus_mask, intensity_image=im_marker)
        marker_nuc_props = None
        if regions:
            marker_nuc_props = regions[-1]

        self.marker_nuc_props = marker_nuc_props

    def set_single_cell_marker_cytosol_props(self, single_cytosol_mask, im_marker):
        regions = skimage.measure.regionprops(single_cytosol_mask.astype(int), intensity_image=im_marker)
        marker_nuc_cyt_props = None
        if regions:
            marker_nuc_cyt_props = regions[-1]

        self.marker_nuc_cyt_props = marker_nuc_cyt_props

    def set_single_cell_junction_props(
            self, single_membrane_mask, im_junction, single_cell_junction_protein,
            single_junction_protein_area_mask, cell_major_axis_length, cell_minor_axis_length
    ):
        # get junction properties. According to junction mapper: https://doi.org/10.7554/eLife.45413
        interface_props = get_single_cell_prop(single_membrane_mask.astype(int), intensity=im_junction)

        # get circular junction props  # todo: radius = minor_axis/2 ???
        r = (cell_major_axis_length - cell_minor_axis_length) / 2
        n_img = map_single_cell_to_circle(single_cell_junction_protein, interface_props.centroid[0],
                                          interface_props.centroid[1], r)
        circular_junction_props = get_single_cell_prop(n_img.astype(bool).astype(int), intensity=n_img)

        # membrane perimeter
        interface_perimeter = interface_props.perimeter  # todo: use me on dilated mask!!!!!!!

        # single_cell_junction_protein mask, area & perimeter
        junction_protein_area_props = get_single_cell_prop(single_junction_protein_area_mask.astype(int),
                                                           intensity=single_cell_junction_protein)
        # junction_fragmented_perimeter = junction_protein_area_props.perimeter  # todo: this seems not correct! Perimeter calculation on discontinous parts seems not possible.

        # secondary statistics  # interface area = junction_props.area
        # junction_coverage_index = junction_fragmented_perimeter / junction_perimeter  # todo: calc junction_fragmented_perimeter
        junction_interface_occupancy = junction_protein_area_props.area / interface_props.area
        junction_protein_intensity = junction_protein_area_props.mean_intensity * junction_protein_area_props.area
        junction_intensity_per_interface_area = junction_protein_intensity / interface_props.area
        junction_cluster_density = junction_protein_intensity / junction_protein_area_props.area

        self.junction_props = SCJunctionProps(
            interface_props, circular_junction_props, interface_perimeter, junction_protein_area_props, None,
            junction_interface_occupancy, junction_protein_intensity,
            junction_intensity_per_interface_area, junction_cluster_density
        )


def set_single_cell_nucleus_props(properties_dataset, index, single_nucleus_mask):
    regions = skimage.measure.regionprops(single_nucleus_mask)
    nucleus_props = None
    if regions:
        nucleus_props = regions[-1]
        fill_single_nucleus_data_frame(properties_dataset, index, nucleus_props)

    return nucleus_props


def set_single_cell_organelle_props(properties_dataset, index, single_organelle_mask, nucleus_props):
    regions = skimage.measure.regionprops(single_organelle_mask)
    organelle_props = None
    if regions:
        organelle_props = regions[-1]
        fill_single_cell_organelle_data_frame(properties_dataset, index, organelle_props, nucleus_props)

    return organelle_props


def set_single_cell_props(properties_dataset, index, filename, connected_component_label, single_cell_mask):
    single_cell_props = get_single_cell_prop(single_cell_mask)
    fill_single_cell_general_data_frame(
        properties_dataset, index, filename, connected_component_label, single_cell_props
    )

    return single_cell_props


def set_single_cell_marker_nuclei_props(properties_dataset, index, single_nucleus_mask, im_marker):
    regions = skimage.measure.regionprops(single_nucleus_mask, intensity_image=im_marker)
    marker_nuc_props = None
    if regions:
        marker_nuc_props = regions[-1]
        fill_single_cell_marker_nuclei_data_frame(properties_dataset, index, marker_nuc_props)

    return marker_nuc_props


def set_single_cell_marker_cytosol_props(properties_dataset, index, single_cytosol_mask, im_marker):
    regions = skimage.measure.regionprops(single_cytosol_mask.astype(int), intensity_image=im_marker)
    marker_nuc_cyt_props = None
    if regions:
        marker_nuc_cyt_props = regions[-1]
        fill_single_cell_marker_nuclei_cytosol_data_frame(properties_dataset, index, marker_nuc_cyt_props)

    return marker_nuc_cyt_props


def set_single_cell_marker_membrane_props(properties_dataset, index, single_membrane_mask, im_marker):
    regions = skimage.measure.regionprops(single_membrane_mask.astype(int), intensity_image=im_marker)
    marker_membrane_props = None
    if regions:
        marker_membrane_props = regions[-1]
        fill_single_cell_marker_membrane_data_frame(properties_dataset, index, marker_membrane_props)

    return marker_membrane_props


def set_single_cell_marker_props(properties_dataset, index, single_cell_mask, im_marker):
    regions = skimage.measure.regionprops(single_cell_mask, intensity_image=im_marker)
    marker_props = None
    if regions:
        marker_props = regions[-1]
        fill_single_cell_marker_polarity(properties_dataset, index, marker_props)

    return marker_props


def set_single_cell_junction_props(
        properties_dataset, index, single_membrane_mask, im_junction, single_cell_junction_protein,
        single_junction_protein_area_mask
):
    # get junction properties. According to junction mapper: https://doi.org/10.7554/eLife.45413
    interface_props = get_single_cell_prop(single_membrane_mask.astype(int), intensity=im_junction)

    # get circular junction props  # todo: radius = minor_axis/2 ???
    r = (properties_dataset["cell_major_axis_length"].values[0] - properties_dataset["cell_minor_axis_length"].values[
        0]) / 2
    n_img = map_single_cell_to_circle(single_cell_junction_protein, interface_props.centroid[0],
                                      interface_props.centroid[1], r)
    circular_junction_props = get_single_cell_prop(n_img.astype(bool).astype(int), intensity=n_img)

    # membrane perimeter
    interface_perimeter = interface_props.perimeter  # todo: use me on dilated mask!!!!!!!

    # single_cell_junction_protein mask, area & perimeter
    junction_protein_area_props = get_single_cell_prop(single_junction_protein_area_mask.astype(int),
                                                       intensity=single_cell_junction_protein)
    # junction_fragmented_perimeter = junction_protein_area_props.perimeter  # todo: this seems not correct! Perimeter calculation on discontinous parts seems not possible.

    # secondary statistics  # interface area = junction_props.area
    # junction_coverage_index = junction_fragmented_perimeter / junction_perimeter  # todo: calc junction_fragmented_perimeter
    junction_interface_occupancy = junction_protein_area_props.area / interface_props.area
    junction_protein_intensity = junction_protein_area_props.mean_intensity * junction_protein_area_props.area
    junction_intensity_per_interface_area = junction_protein_intensity / interface_props.area
    junction_cluster_density = junction_protein_intensity / junction_protein_area_props.area

    fill_single_cell_junction_data_frame(
        properties_dataset,
        index,
        interface_perimeter,
        interface_props,
        junction_protein_area_props.area,
        circular_junction_props.centroid_weighted
    )
    fill_single_cell_junction_sec_stat_data_frame(
        properties_dataset, index, junction_interface_occupancy,
        junction_intensity_per_interface_area,
        junction_cluster_density
    )


def get_single_cell_prop(single_cell_mask, intensity=None):
    """Gets the single cell properties."""
    # we are looking at a single cell. There is only one region!
    regions = skimage.measure.regionprops(single_cell_mask, intensity_image=intensity)
    if len(regions) > 1:
        raise ValueError("Too many regions for a single cell!")
    props = regions[-1]

    return props
