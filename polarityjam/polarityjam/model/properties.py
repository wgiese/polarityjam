import numpy as np
from skimage.measure._regionprops import RegionProperties

from polarityjam.compute.compute import map_single_cell_to_circle, compute_single_cell_prop, \
    compute_reference_target_orientation_rad, compute_angle_deg, compute_marker_vector_norm, compute_shape_orientation, \
    straight_line_length
from polarityjam.compute.corner import get_corner
from polarityjam.model.masks import SingleCellMasksCollection, MasksCollection
from polarityjam.utils import parameters
from scipy import ndimage as ndi


class SingleCellProps(RegionProperties):
    """Base class for all single cell properties.

    """

    def __init__(self, single_cell_mask, intensity=None):
        if not np.issubdtype(single_cell_mask.dtype, np.integer):
            raise RuntimeError("Only integer images allowed!")

        objects = ndi.find_objects(single_cell_mask)

        if len(objects) > 1:
            raise RuntimeError("Several objects detected in single cell mask! Aborting...")

        sl = objects[0]

        super().__init__(sl, 1, single_cell_mask, intensity, True)


class SingleCellCellProps(SingleCellProps):

    def __init__(self, single_cell_mask):
        super().__init__(single_cell_mask)

    @property
    def cell_shape_orientation(self):
        # note, the values of orientation from props are in [-pi/2,pi/2] with zero along the y-axis
        return compute_shape_orientation(self.orientation)

    @property
    def cell_major_to_minor_ratio(self):
        return self.major_axis_length / self.minor_axis_length


class SingleCellNucleusProps(SingleCellProps):
    def __init__(self, single_nucleus_mask):
        super().__init__(single_nucleus_mask)

    @property
    def nucleus_displacement_orientation_rad(self):
        return compute_reference_target_orientation_rad(
            self.single_cell_props.centroid[0], self.single_cell_props.centroid[1], self.centroid[0], self.centroid[1]
        )

    @property
    def nucleus_displacement_orientation_deg(self):
        return compute_angle_deg(self.nucleus_displacement_orientation_rad())

    @property
    def nuc_shape_orientation(self):
        # note, the values of orientation from props are in [-pi/2,pi/2] with zero along the y-axis
        return compute_shape_orientation(self.orientation)  # TODO: should be nuc_shape_orientation_deg and rad ?

    @property
    def nuc_major_to_minor_ratio(self):
        return self.major_axis_length / self.minor_axis_length


class SingleCellOrganelleProps(SingleCellProps):
    def __init__(self, single_organelle_mask, nucleus_props):
        super().__init__(single_organelle_mask)

        self._nucleus_props = nucleus_props

    @property
    def organelle_distance(self):
        return compute_marker_vector_norm(
            self.centroid[0], self.centroid[1], self._nucleus_props.centroid[0], self._nucleus_props.centroid[1]
        )

    @property
    def organelle_orientation_rad(self):
        return compute_reference_target_orientation_rad(
            self.nucleus_props.centroid[0], self.nucleus_props.centroid[1], self.centroid[0], self.centroid[1]
        )

    @property
    def organelle_orientation_deg(self):
        return compute_angle_deg(self.organelle_orientation_rad)


class SingleCellMarkerProps(SingleCellProps):
    def __init__(self, single_cell_mask, im_marker):
        super().__init__(single_cell_mask, im_marker)

    @property
    def marker_centroid_orientation_rad(self):
        return compute_reference_target_orientation_rad(
            self.centroid[0], self.centroid[1], self.weighted_centroid[0], self.weighted_centroid[1]
        )

    @property
    def marker_centroid_orientation_deg(self):
        return compute_angle_deg(self.marker_centroid_orientation_rad)


class SingleCellMarkerMembraneProps(SingleCellProps):
    def __init__(self, single_membrane_mask, im_marker):
        super().__init__(single_membrane_mask, im_marker)

    @property
    def marker_sum_expression_mem(self):
        return self.mean_intensity * self.area


class SingleCellMarkerNucleiProps(SingleCellProps):
    def __init__(self, single_nucleus_mask, im_marker, nucleus_props, marker_props):
        super().__init__(single_nucleus_mask, im_marker)
        self._nucleus_props = nucleus_props
        self._marker_props = marker_props

    @property
    def marker_nucleus_orientation_rad(self):
        return compute_reference_target_orientation_rad(
            self._nucleus_props.centroid[0],
            self._nucleus_props.centroid[1],
            self._marker_props.centroid[0],
            self._marker_props.centroid[0]
        )

    @property
    def marker_nucleus_orientation_deg(self):
        return compute_angle_deg(self.marker_nucleus_orientation_rad)

    @property
    def marker_sum_expression_nuc(self):
        return self.mean_intensity * self.area


class SingleCellMarkerCytosolProps(SingleCellProps):
    def __init__(self, single_cytosol_mask, im_marker):
        super().__init__(single_cytosol_mask, im_marker)

    @property
    def marker_sum_expression_cyt(self):
        return self.mean_intensity * self.area


class SingleCellJunctionInterfaceProps(SingleCellProps):
    # Based on junction mapper: https://doi.org/10.7554/eLife.45413
    def __init__(self, single_membrane_mask, im_junction):
        super().__init__(single_membrane_mask, im_junction)


class SingleCellJunctionProteinProps(SingleCellProps):
    def __init__(self, single_junction_protein_area_mask,
                 im_junction_protein_single_cell):
        super().__init__(single_junction_protein_area_mask, im_junction_protein_single_cell)


class SingleCellJunctionProteinCircularProps(SingleCellProps):
    def __init__(self, im_junction_protein_single_cell, cell_minor_axis_length, interface_centroid):
        r = cell_minor_axis_length / 2
        circular_img = map_single_cell_to_circle(im_junction_protein_single_cell, interface_centroid[0],
                                                 interface_centroid[1], r)
        circular_junction_protein_single_cell_mask = circular_img.astype(bool).astype(int)

        super().__init__(circular_junction_protein_single_cell_mask, circular_img)


class SingleCellJunctionProps:
    def __init__(self, sc_junction_interface_props, sc_junction_protein_props, sc_junction_protein_circular_props,
                 sc_mask):
        self.sc_mask = sc_mask
        self.sc_junction_interface_props = sc_junction_interface_props
        self.sc_junction_protein_props = sc_junction_protein_props
        self.sc_junction_protein_circular_props = sc_junction_protein_circular_props

    @property
    def straight_line_junction_length(self):
        return straight_line_length(get_corner(self.sc_mask, parameters.dp_epsilon))

    @property
    def interface_perimeter(self):
        return 1  # todo: perimeter on dilated mask

    @property
    def junction_interface_occupancy(self):
        return self.sc_junction_protein_props.area / self.sc_junction_interface_props.area

    @property
    def junction_protein_intensity(self):
        return self.sc_junction_protein_props.mean_intensity * self.sc_junction_protein_props.area

    @property
    def junction_intensity_per_interface_area(self):
        return self.junction_protein_intensity / self.sc_junction_interface_props.area

    @property
    def junction_cluster_density(self):
        return self.junction_protein_intensity / self.sc_junction_protein_props.area


class SingleCellJunctionPropss:
    def __init__(
            self, interface_props,
            circular_junction_props,
            interface_perimeter,
            junction_protein_area_props,
            junction_fragmented_perimeter,
            junction_interface_occupancy,
            junction_protein_intensity,
            junction_intensity_per_interface_area,
            junction_cluster_density,
            corners
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
        self.corners = corners


class NeighborhoodProps:
    def __init__(self):
        self.num_neighbours = None
        self.mean_dif_first_neighbors = 0
        self.median_dif_first_neighbors = 0
        self.var_dif_first_neighbors = 0
        self.range_dif_first_neighbors = 0
        self.mean_dif_second_neighbors = 0
        self.median_dif_second_neighbors = 0
        self.var_dif_second_neighbors = 0
        self.range_dif_second_neighbors = 0


class SingleCellPropertiesCollection:

    def __init__(self):
        self.nucleus_props = None
        self.organelle_props = None
        self.single_cell_props = None
        self.marker_nuc_props = None
        self.marker_nuc_cyt_props = None
        self.marker_membrane_props = None
        self.marker_props = None
        self.junction_props = None

    def set_sc_props(self, sc_masks: SingleCellMasksCollection, masks: MasksCollection, im_marker,
                     im_junction):
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
                sc_masks.sc_mask,
                im_junction,
                sc_masks.sc_junction_protein,
                sc_masks.sc_junction_protein_area_mask,
                self.single_cell_props.major_axis_length,
                self.single_cell_props.minor_axis_length
            )

    def set_single_cell_props(self, single_cell_mask):
        p = compute_single_cell_prop(single_cell_mask)
        # note, the values of orientation from props are in [-pi/2,pi/2] with zero along the y-axis
        p.cell_shape_orientation = compute_shape_orientation(p.orientation)
        p.cell_major_to_minor_ratio = p.major_axis_length / p.minor_axis_length

        self.single_cell_props = p

    def set_single_cell_nucleus_props(self, single_nucleus_mask):
        p = compute_single_cell_prop(single_nucleus_mask)

        # compute nucleus displacement
        angle_rad = compute_reference_target_orientation_rad(
            self.single_cell_props.centroid[0], self.single_cell_props.centroid[1], p.centroid[0], p.centroid[1]
        )
        p.nucleus_displacement_orientation_rad = angle_rad
        p.nucleus_displacement_orientation_deg = compute_angle_deg(angle_rad)

        # note, the values of orientation from props are in [-pi/2,pi/2] with zero along the y-axis
        p.nuc_shape_orientation = compute_shape_orientation(
            p.orientation)  # TODO: should be nuc_shape_orientation_deg and rad ?

        p.nuc_major_to_minor_ratio = p.major_axis_length / p.minor_axis_length

        self.nucleus_props = p

    def set_single_cell_organelle_props(self, single_organelle_mask):
        p = compute_single_cell_prop(single_organelle_mask)

        x_organelle, y_organelle = p.centroid
        angle_rad = compute_reference_target_orientation_rad(
            self.nucleus_props.centroid[0], self.nucleus_props.centroid[1], x_organelle, y_organelle
        )

        p.organelle_distance = compute_marker_vector_norm(
            x_organelle, y_organelle, self.nucleus_props.centroid[0], self.nucleus_props.centroid[1]

        )
        p.organelle_orientation_rad = angle_rad
        p.organelle_orientation_deg = compute_angle_deg(angle_rad)

        self.organelle_props = p

    def set_single_cell_marker_props(self, single_cell_mask, im_marker):
        p = compute_single_cell_prop(single_cell_mask, intensity=im_marker)
        marker_centroid_x, marker_centroid_y = p.weighted_centroid
        cell_x, cell_y = p.centroid
        angle_rad = compute_reference_target_orientation_rad(cell_x, cell_y, marker_centroid_x, marker_centroid_y)

        p.marker_centroid_orientation_rad = angle_rad
        p.marker_centroid_orientation_deg = compute_angle_deg(angle_rad)

        self.marker_props = p

    def set_single_cell_marker_membrane_props(self, single_membrane_mask, im_marker):
        p = compute_single_cell_prop(single_membrane_mask.astype(int), intensity=im_marker)

        p.marker_sum_expression_mem = p.mean_intensity * p.area

        self.marker_membrane_props = p

    def set_single_cell_marker_nuclei_props(self, single_nucleus_mask, im_marker):
        p = compute_single_cell_prop(single_nucleus_mask, intensity=im_marker)

        # compute marker nucleus orientation
        angle_rad = compute_reference_target_orientation_rad(
            self.nucleus_props.centroid[0],
            self.nucleus_props.centroid[1],
            self.marker_props.centroid[0],
            self.marker_props.centroid[0]
        )

        p.marker_sum_expression_nuc = p.mean_intensity * p.area
        p.marker_nucleus_orientation_rad = angle_rad
        p.marker_nucleus_orientation_deg = compute_angle_deg(angle_rad)

        self.marker_nuc_props = p

    def set_single_cell_marker_cytosol_props(self, single_cytosol_mask, im_marker):
        p = compute_single_cell_prop(single_cytosol_mask.astype(int), intensity=im_marker)

        p.marker_sum_expression_cyt = p.mean_intensity * p.area

        self.marker_nuc_cyt_props = p

    def set_single_cell_junction_props(
            self, single_membrane_mask, sc_mask, im_junction, single_cell_junction_protein,
            single_junction_protein_area_mask, cell_major_axis_length, cell_minor_axis_length
    ):
        # get corners
        corners = get_corner(sc_mask, parameters.dp_epsilon)

        # straight-line junction length
        straight_line_junction_length = straight_line_length(corners)

        # get junction properties. According to junction mapper: https://doi.org/10.7554/eLife.45413
        interface_props = compute_single_cell_prop(single_membrane_mask.astype(int), intensity=im_junction)

        # get circular junction props  # todo: radius = minor_axis/2 ???
        r = cell_minor_axis_length / 2
        n_img = map_single_cell_to_circle(single_cell_junction_protein, interface_props.centroid[0],
                                          interface_props.centroid[1], r)
        circular_junction_props = compute_single_cell_prop(n_img.astype(bool).astype(int), intensity=n_img)

        # membrane perimeter
        interface_perimeter = interface_props.perimeter  # todo: use me on dilated mask!!!!!!!

        # single_cell_junction_protein mask, area & perimeter
        junction_protein_area_props = compute_single_cell_prop(single_junction_protein_area_mask.astype(int),
                                                               intensity=single_cell_junction_protein)
        # junction_fragmented_perimeter = junction_protein_area_props.perimeter  # todo: this seems not correct! Perimeter calculation on discontinous parts seems not possible.

        # secondary statistics  # interface area = junction_props.area
        # junction_coverage_index = junction_fragmented_perimeter / junction_perimeter  # todo: calc junction_fragmented_perimeter
        junction_interface_occupancy = junction_protein_area_props.area / interface_props.area
        junction_protein_intensity = junction_protein_area_props.mean_intensity * junction_protein_area_props.area
        junction_intensity_per_interface_area = junction_protein_intensity / interface_props.area
        junction_cluster_density = junction_protein_intensity / junction_protein_area_props.area

        junctions_centroid_x, junctions_centroid_y = circular_junction_props.weighted_centroid

        junction_protein_area = junction_protein_area_props.area

        junction_centroid_X = junctions_centroid_x
        junction_centroid_Y = junctions_centroid_y

        # junction_protein_intensity = interface_props.mean_intensity * interface_props.area

        self.junction_props = SingleCellJunctionPropss(
            interface_props,
            circular_junction_props,
            interface_perimeter,
            junction_protein_area_props,
            None,
            junction_interface_occupancy,
            junction_protein_intensity,
            junction_intensity_per_interface_area,
            junction_cluster_density,
            corners
        )
