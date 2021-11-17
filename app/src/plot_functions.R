
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

rose_plot_circular <- function(parameters, input, statistics, feature_circular) {
  
  bin_size = 360/input$bins
  
  #feature <- parameters[input$feature_select][[1]][1]
  plot_title <- parameters[input$feature_select][[1]][3]
  
  p <- ggplot() +
    geom_histogram(aes(x = feature_circular, y = ..ncount..),
                   breaks = seq(0, 360, bin_size),
                   colour = "black",
                   fill = "black",
                   alpha = 0.5) +
    ggtitle(plot_title) +
    theme(axis.text.x = element_text(size = 18)) +
    coord_polar(start = -pi/2.0, direction = -1) +
    scale_x_continuous(limits = c(0, 360),
                       breaks = (c(0, 90, 180, 270))) +
    theme_minimal(base_size = 14) +
    xlab(sprintf("number of cells = : %s \n condition: %s", length(feature_circular), input$exp_condition)) +
    ylab("polarity index") 
  #theme(axis.text.y=element_blank()) +
  
  if (input$area_scaled) {
    p <- p + scale_y_sqrt()
  }
  
  p <- p + geom_segment(data = statistics, aes(x=angle_mean_deg, y=0, xend=angle_mean_deg, yend=polarity_index, size = 1.5, color="red", lineend = "butt"), arrow = arrow()) + theme(legend.position = "none") 
  return(p)
  
}


rose_plot_2_axial <- function(parameters, input, feature_circular) {
  
  #bin_size = 360/input$bins
  plot_title <- parameters[input$feature_select][[1]][3]
  
  p <- ggplot() +
    geom_histogram(aes(x = feature_circular, y = ..ncount..),
                   breaks = seq(0, 360, bin_size),
                   colour = "black",
                   fill = "black",
                   alpha = 0.5) +
    ggtitle(plot_title) +
    theme(axis.text.x = element_text(size = 18)) +
    coord_polar(start = -pi/2.0, direction = -1) +
    scale_x_continuous(limits = c(0, 360),
                       breaks = (c(0, 90, 180, 270))) +
    theme_minimal(base_size = 14) +
    xlab(sprintf("number of cells = : %s \n condition: %s", length(feature_circular), input$exp_condition)) #+
    #ylab("polarity index") 
  #theme(axis.text.y=element_blank()) +
  
  if (input$area_scaled) {
    p <- p + scale_y_sqrt()
  }
  
  return(p)
  
}

linear_histogram <- function(parameters, input, feature_linear) {
  
  bin_size = max(feature_linear)/input$bins
  plot_title <- parameters[input$feature_select][[1]][3]
  
  p <- ggplot() +
    geom_histogram(aes(x = feature_linear, y = ..ncount..),
                   breaks = seq(0, max(feature_linear), bin_size),
                   colour = "black",
                   fill = "black",
                   alpha = 0.5) +
    ggtitle(plot_title) +
    theme(axis.text.x = element_text(size = 18)) +
    theme_minimal(base_size = 14) +
    xlab(sprintf("number of cells = : %s \n condition: %s", length(feature_linear), input$exp_condition)) 

  return(p)
  
}



