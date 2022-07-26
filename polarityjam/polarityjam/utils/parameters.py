## cellpose/segmentation parameters
manually_annotated_mask = ""  # in case a manually annotated cellpose mask exists specify suffix here
store_segmentation = False
use_given_mask = True

cp_model_type = "cyto"  # "cyto", "cyto2", "custom"
cp_model_path = ""
estimated_cell_diameter = 100

use_gpu = False
clear_border = True
min_cell_size = 50
min_nucleus_size = 10
min_organelle_size = 10

dp_epsilon = 5

## spatial statistics
feature_of_interest = "cell_area"

## plot options
plot_junctions = True
plot_polarity = True
plot_orientation = True
plot_marker = True
plot_ratio_method = False
plot_cyclic_orientation = False

outline_width = 2
show_polarity_angles = True
show_graphics_axis = False

graphics_output_format = ['png', 'pdf']  # 'png', 'pdf', 'svg'
dpi = 300
graphics_width = 5  # figure width in inches
graphics_height = 5  # figure height in inches

# channel parameters
channel_junction = 0
channel_nucleus = 1
channel_organelle = 2
channel_expression_marker = -1
membrane_thickness = 5
