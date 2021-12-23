
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

rose_plot_circular <- function(parameters, input, statistics, feature_circular, plot_title, text_size = 24) {
  
    bin_size = 360/input$bins
  
    #feature <- parameters[input$feature_select][[1]][1]
    #plot_title <- parameters[input$feature_select][[1]][3]

    polarity_index <- signif(statistics[1,"polarity_index"], digits=3)
    p_value_ <- signif(statistics[1,"rayleigh_test"], digits=3)
    if (statistics[1,"rayleigh_test"] < 0.001) 
        p_value <- "p < 0.001"
    else
        p_value <- paste0("p = ", toString(p_value_))
 
 
    p <- ggplot() +
        geom_histogram(aes(x = feature_circular, y = ..ncount..),
                   breaks = seq(0, 360, bin_size),
                   colour = "black",
                   fill = "black",
                   alpha = 0.5) +
        #geom_density(aes(x = feature_circular)) +
        ggtitle(plot_title) +
        theme(axis.text.x = element_text(size = 18)) +
        coord_polar(start = -pi/2.0, direction = -1) +
        scale_x_continuous(limits = c(0, 360),
                       breaks = (c(0, 90, 180, 270))) +
        theme_minimal(base_size = text_size) +
        xlab(sprintf("number of cells = : %s \n polarity index: %s, %s, \n condition: %s", length(feature_circular), polarity_index, p_value, input$exp_condition)) +
        ylab("polarity index") 
  #theme(axis.text.y=element_blank()) +
  
    if (input$area_scaled) {
            p <- p + scale_y_sqrt()
    }

 
    p <- p + geom_segment(data = statistics, aes(x=mean, y=0, xend=mean, yend=polarity_index, size = 1.5, color="red", lineend = "butt"), arrow = arrow()) + theme(legend.position = "none") 
    
    return(p)
  
}

compare_plot_circular <- function(parameters, input, statistics, feature_circular_1, feature_circular_2, plot_title, text_size = 24) {
  
    bin_size = 360/input$bins_comparison
  
    #feature <- parameters[input$feature_select][[1]][1]
    #plot_title <- parameters[input$feature_select][[1]][3]

    polarity_index <- signif(statistics[1,"polarity_index"], digits=3)
    p_value_ <- signif(statistics[1,"rayleigh_test"], digits=3)
    if (statistics[1,"rayleigh_test"] < 0.001) 
        p_value <- "p < 0.001"
    else
        p_value <- paste0("p = ", toString(p_value_))
 
    
    p <- ggplot() 

    if (input$kde_comparison) {
    p <- p + geom_histogram(aes(x = feature_circular_1, y = ..density..),#, y = ..ncount..),
                   breaks = seq(0, 360, bin_size),
                   colour = "black",
                   fill = "blue",
                   alpha = 0.5) +
        geom_histogram(aes(x = feature_circular_2, y = ..density..),#, y = ..ncount..),
                   breaks = seq(0, 360, bin_size),
                   colour = "black",
                   fill = "red",
                   alpha = 0.5)
    }
    
    if (input$histogram_comparison) {
        p <- p + geom_density(aes(x = feature_circular_1, y = ..density.. ),
                    colour = "black",
                    fill = "red",
                    alpha = 0.5) +
        geom_density(aes(x = feature_circular_2, y = ..density.. ),
                    colour = "black",
                    fill = "blue",
                    alpha = 0.5) 
    }
    p <- p + ggtitle(plot_title) +
        theme(axis.text.x = element_text(size = 18)) +
        coord_polar(start = -pi/2.0, direction = -1) +
        scale_x_continuous(limits = c(0, 360),
                       breaks = (c(0, 90, 180, 270))) +
        theme_minimal(base_size = text_size) +
        xlab(sprintf("number of cells = : %s \n polarity index: %s, %s, \n condition: %s", length(feature_circular_1), polarity_index, p_value, input$exp_condition)) +
        ylab("polarity index") 
  #theme(axis.text.y=element_blank()) +
  
    if (input$area_scaled) {
            p <- p + scale_y_sqrt()
    }

 
    #p <- p + geom_segment(data = statistics, aes(x=mean, y=0, xend=mean, yend=polarity_index, size = 1.5, color="red", lineend = "butt"), arrow = arrow()) + theme(legend.position = "none") 
    
    return(p)
  
}




rose_plot_2_axial <- function(parameters, input, statistics, feature_circular, plot_title, text_size = 24) {
  
    bin_size = 360/input$bins
    #plot_title <- parameters[input$feature_select][[1]][3]
  
    polarity_index <- signif(statistics[1,"polarity_index"], digits=3)
    p_value_ <- signif(statistics[1,"rayleigh_test"], digits=3)
    if (statistics[1,"rayleigh_test"] < 0.001)
        p_value <- "p < 0.001"
    else
        p_value <- paste0("p = ", toString(p_value_))
 
    p <- ggplot() +
        geom_histogram(aes(x = feature_circular, y = ..ncount..),
                   breaks = seq(0, 360, bin_size),
                   colour = "black",
                   fill = "black",
                   alpha = 0.5) +
#        geom_density(aes(x = feature_circular)) +
        ggtitle(plot_title) +
        theme(axis.text.x = element_text(size = 18)) +
        coord_polar(start = -pi/2.0, direction = -1) +
        scale_x_continuous(limits = c(0, 360),
                       breaks = (c(0, 90, 180, 270))) +
        theme_minimal(base_size = text_size) +
        xlab(sprintf("number of cells = : %s \n polarity index: %s, %s, \n condition: %s" , length(feature_circular), polarity_index, p_value, input$exp_condition)) +
        ylab("polarity index") 
  #theme(axis.text.y=element_blank()) +
  
    if (input$area_scaled) {
        p <- p + scale_y_sqrt()
    }

    if (input$left_axial) {
        print(statistics[1,"mean"]) 
        if( statistics[1,"mean"] < 90.0) {
      	    statistics[1,"mean"] = statistics[1,"mean"] + 180.0
    	}
    }
 
 
    p <- p + geom_segment(data = statistics, aes(x=mean, y=0, xend=mean, yend=polarity_index, size = 1.5, color="red", lineend = "butt"), arrow = arrow()) + theme(legend.position = "none") 
  return(p)
  
}

linear_histogram <- function(parameters, input, statistics, feature_linear, plot_title) {
  
    bin_size = max(feature_linear)/input$bins
    #plot_title <- parameters[input$feature_select][[1]][3]
  
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

    p <- p + geom_vline(data = statistics, aes(xintercept=mean),  col="red", size = 2) 
    return(p)
  
}



