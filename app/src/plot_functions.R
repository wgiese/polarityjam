
transform_2_axial <- function(feature_2_axial) {
	feature_2_axial <- unlist(feature_2_axial)
	feature_transformed <- feature_2_axial
	for (i in 1:length(feature_2_axial)) {
		if( feature_2_axial[i] < pi/2.0) {
      			feature_transformed[i]  = feature_2_axial[i] + pi
    		}
    		else {
    			feature_transformed[i]  = feature_2_axial[i] 
    		}
	}
	return(feature_transformed)
}


#rose_plot_2_axial <- function(feature_2_axial)

#rose_plot_circular <- function(feature_circular)

