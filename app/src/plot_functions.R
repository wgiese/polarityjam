
transform_2_axial <- function(input, feature_2_axial) {

    feature_2_axial <- unlist(feature_2_axial)
	feature_transformed <- feature_2_axial

    if (input$hemi_rose_options == "left") {

	    for (i in 1:length(feature_2_axial)) {
	        if( feature_2_axial[i] < pi/2.0) {
      	        feature_transformed[i]  = feature_2_axial[i] + pi
    	    }
    	    else {
    	        feature_transformed[i]  = feature_2_axial[i] 
    	    }
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
    if(input$stats_method == "V-Test")
        p_value_ <- signif(statistics[1,"v_test"], digits=3)
    

    if (statistics[1,"rayleigh_test"] < 0.001) 
        p_value <- "P < 0.001"
    else
        p_value <- paste0("P = ", toString(p_value_))
    
    
    if (input$stats_method == "Watson's Test")
        p_value <-  statistics[1,"watson_test"] 
    if (input$stats_method == "Rao's Test")
        p_value <-  statistics[1,"rao_test"] 


    p <- ggplot()

    colour_fill <- "black"
    colour <- "black"
    
    if (input$select_colourmap == "Okabe_Ito") {
        colour_n <- input$select_colour %% length(Okabe_Ito)
        colour_fill <- Okabe_Ito[colour_n]
        colour <- Okabe_Ito[colour_n]
    }
    
    alpha <- 0.5
    if (input$adjust_alpha == TRUE) {
        alpha <- input$alpha_fill
        if (input$outline != "colour")
            colour <- input$outline
    }


    if (input$kde_plot) {
      p <-  p + 
        geom_density(aes(x = feature_circular, y = ..count../max(..count..) ),
                     colour = colour,
                     fill = colour_fill,
                     alpha = alpha)
    }
    
    if (input$histogram_plot) {
      p <- p + 
          geom_histogram(aes(x = feature_circular, y = ..ncount..),
                        breaks = seq(0, 360, bin_size),
                        colour = colour,
                        fill = colour_fill,
                        alpha = alpha)
    }    
    
    if (input$scatter_plot) {
        print("plot points")
        p + geom_point(aes(x = feature_circular, y = 1))
    }

    p <- p + ggtitle(plot_title) +
        theme(plot.title = element_text(size = 18, face = "bold")) +
        theme(axis.text.x = element_text(size = 18)) +
        coord_polar(start = -pi/2.0, direction = -1) +
        scale_x_continuous(limits = c(0, 360),
                       breaks = (c(0, 90, 180, 270))) +
        scale_y_continuous(limits = c(0, 1.1)) +
        theme_minimal(base_size = text_size) +
        #xlab(sprintf("number of cells = : %s \n polarity index: %s, %s, \n condition: %s", length(feature_circular), polarity_index, p_value, input$exp_condition)) +
        ylab("polarity index") 
        #theme(axis.text.y=element_blank()) +
    

    if (input$stats_method != "None") {
        p<- p + xlab(sprintf("number of cells = : %s \n polarity index: %s, %s \n condition: %s", length(feature_circular), polarity_index, p_value, input$exp_condition)) 
    }
    else {
        p<- p + xlab(sprintf("number of cells = : %s \n polarity index: %s \n condition: %s", length(feature_circular), polarity_index, input$exp_condition)) 
    } 


    if (input$scatter_plot) {
        print("plot points")
        p <- p + geom_point(aes(x = feature_circular, y = 1))

    }


 
    if (input$area_scaled) {
            p <- p + scale_y_sqrt(limits = c(0,sqrt(1.1))) ##+ scale_y_continuous(limits = c(0, sqrt(1.1)))
    }

 
    p <- p + geom_segment(data = statistics, aes(x=mean, y=0, xend=mean, yend=polarity_index, size = 1.5, color="red"), arrow = arrow()) + theme(legend.position = "none") 
    
    if (input$ci_plot) {
        if (input$ci_method == "95% CI of the mean") {
            p <- p + geom_segment(data = statistics, aes(x=ci_95_lower_limit, y=0, xend=ci_95_lower_limit, yend=1), size = 1.5, color="red",linetype = "dashed", arrow = NULL) + theme(legend.position = "none") 
            p <- p + geom_segment(data = statistics, aes(x=ci_95_upper_limit, y=0, xend=ci_95_upper_limit, yend=1), size = 1.5, color="red",linetype = "dashed", arrow = NULL) + theme(legend.position = "none") 
        }
        if (input$ci_method == "90% CI of the mean") {
            p <- p + geom_segment(data = statistics, aes(x=ci_90_lower_limit, y=0, xend=ci_90_lower_limit, yend=1), size = 1.5, color="red",linetype = "dashed", arrow = NULL) + theme(legend.position = "none") 
            p <- p + geom_segment(data = statistics, aes(x=ci_90_upper_limit, y=0, xend=ci_90_upper_limit, yend=1), size = 1.5, color="red",linetype = "dashed", arrow = NULL) + theme(legend.position = "none") 
        }
        if (input$ci_method == "50% CI of the mean") {
            p <- p + geom_segment(data = statistics, aes(x=ci_50_lower_limit, y=0, xend=ci_50_lower_limit, yend=1), size = 1.5, color="red",linetype = "dashed", arrow = NULL) + theme(legend.position = "none") 
            p <- p + geom_segment(data = statistics, aes(x=ci_50_upper_limit, y=0, xend=ci_50_upper_limit, yend=1), size = 1.5, color="red",linetype = "dashed", arrow = NULL) + theme(legend.position = "none") 
        }
        if (input$ci_method == "circular standard deviation") {
            p <- p + geom_segment(data = statistics, aes(x=std_circ_low_lim, y=0, xend=std_circ_low_lim, yend=1), size = 1.5, color="red",linetype = "dashed", arrow = NULL) + theme(legend.position = "none") 
            p <- p + geom_segment(data = statistics, aes(x=std_circ_up_lim, y=0, xend=std_circ_up_lim, yend=1), size = 1.5, color="red",linetype = "dashed", arrow = NULL) + theme(legend.position = "none") 
        }
        if (input$ci_method == "angular standard deviation") {
            p <- p + geom_segment(data = statistics, aes(x=std_ang_low_lim, y=0, xend=std_ang_low_lim, yend=1), size = 1.5, color="red",linetype = "dashed", arrow = NULL) + theme(legend.position = "none") 
            p <- p + geom_segment(data = statistics, aes(x=std_ang_up_lim, y=0, xend=std_ang_up_lim, yend=1), size = 1.5, color="red",linetype = "dashed", arrow = NULL) + theme(legend.position = "none") 
        }
    }

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

    if (input$histogram_comparison) {
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
    
    if (input$kde_comparison) {
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
 
    if (input$hemi_rose_options == "all") {
        n <- length(feature_circular)
        feature_circular_ <- numeric(2*n)
        for (i in 1:n) {
	        feature_circular_[i]  = feature_circular[i] 
    	    feature_circular_[i+n]  = feature_circular[i] + 180.0
	    }
    }
    else {
        feature_circular_ <- feature_circular 
    }




    p <- ggplot() +
        geom_histogram(aes(x = feature_circular_, y = ..ncount..),
                   breaks = seq(0, 360, bin_size),
                   colour = "black",
                   fill = "black",
                   alpha = 0.5) +
#        geom_density(aes(x = feature_circular)) +
        ggtitle(plot_title) +
        theme(plot.title = element_text(size = 10, face = "bold")) +
        theme(axis.text.x = element_text(size = 18)) +
        coord_polar(start = -pi/2.0, direction = -1) +
        scale_x_continuous(limits = c(0, 360),
                       breaks = (c(0, 90, 180, 270))) +
        theme_minimal(base_size = text_size) +
#        xlab(sprintf("number of cells = : %s \n polarity index: %s, %s, \n condition: %s" , length(feature_circular), polarity_index, p_value, input$exp_condition)) +
        ylab("polarity index") 
  #theme(axis.text.y=element_blank()) +
 
    if (input$stats_method != "None") {
        p<- p + xlab(sprintf("number of cells = : %s \n polarity index: %s, %s \n condition: %s", length(feature_circular), polarity_index, p_value, input$exp_condition)) 
    }
    else {
        p<- p + xlab(sprintf("number of cells = : %s \n polarity index: %s \n condition: %s", length(feature_circular), polarity_index, input$exp_condition)) 
    } 


 
    if (input$area_scaled) {
        p <- p + scale_y_sqrt()
    }

    if (input$hemi_rose_options == "left") {
    #if (input$left_axial) {
        print(statistics[1,"mean"]) 
        if( statistics[1,"mean"] < 90.0) {
      	    statistics[1,"mean"] = statistics[1,"mean"] + 180.0
    	}
        if( statistics[1,"ci_95_lower_limit"] < 90.0) {
      	    statistics[1,"ci_95_lower_limit"] = statistics[1,"ci_95_lower_limit"] + 180.0
    	}
        if( statistics[1,"ci_95_upper_limit"] < 90.0) {
      	    statistics[1,"ci_95_upper_limit"] = statistics[1,"ci_95_upper_limit"] + 180.0
    	}
        if( statistics[1,"ci_90_lower_limit"] < 90.0) {
      	    statistics[1,"ci_90_lower_limit"] = statistics[1,"ci_90_lower_limit"] + 180.0
    	}
        if( statistics[1,"ci_50_upper_limit"] < 90.0) {
      	    statistics[1,"ci_50_upper_limit"] = statistics[1,"ci_50_upper_limit"] + 180.0
    	}
   
    }
    
    

 
    p <- p + geom_segment(data = statistics, aes(x=mean, y=0, xend=mean, yend=polarity_index, size = 1.5, color="red", lineend = "butt"), arrow = NULL) + theme(legend.position = "none") 
 
    if (input$ci_plot) {
        #p <- p + geom_segment(data = statistics, aes(x=ci_lower_limit, y=0, xend=ci_lower_limit, yend=1), size = 1.5, color="red",linetype = "dashed", arrow = NULL) + theme(legend.position = "none") 
        #p <- p + geom_segment(data = statistics, aes(x=ci_upper_limit, y=0, xend=ci_upper_limit, yend=1), size = 1.5, color="red",linetype = "dashed", arrow = NULL) + theme(legend.position = "none") 
        if (input$ci_method == "95% CI of the mean") {
            p <- p + geom_segment(data = statistics, aes(x=ci_95_lower_limit, y=0, xend=ci_95_lower_limit, yend=1), size = 1.5, color="red",linetype = "dashed", arrow = NULL) + theme(legend.position = "none") 
            p <- p + geom_segment(data = statistics, aes(x=ci_95_upper_limit, y=0, xend=ci_95_upper_limit, yend=1), size = 1.5, color="red",linetype = "dashed", arrow = NULL) + theme(legend.position = "none") 
        }
        if (input$ci_method == "90% CI of the mean") {
            p <- p + geom_segment(data = statistics, aes(x=ci_90_lower_limit, y=0, xend=ci_90_lower_limit, yend=1), size = 1.5, color="red",linetype = "dashed", arrow = NULL) + theme(legend.position = "none") 
            p <- p + geom_segment(data = statistics, aes(x=ci_90_upper_limit, y=0, xend=ci_90_upper_limit, yend=1), size = 1.5, color="red",linetype = "dashed", arrow = NULL) + theme(legend.position = "none") 
        }
        if (input$ci_method == "50% CI of the mean") {
            p <- p + geom_segment(data = statistics, aes(x=ci_50_lower_limit, y=0, xend=ci_50_lower_limit, yend=1), size = 1.5, color="red",linetype = "dashed", arrow = NULL) + theme(legend.position = "none") 
            p <- p + geom_segment(data = statistics, aes(x=ci_50_upper_limit, y=0, xend=ci_50_upper_limit, yend=1), size = 1.5, color="red",linetype = "dashed", arrow = NULL) + theme(legend.position = "none") 
        }
        if (input$ci_method == "circular standard deviation") {
            p <- p + geom_segment(data = statistics, aes(x=std_circ_low_lim, y=0, xend=std_circ_low_lim, yend=1), size = 1.5, color="red",linetype = "dashed", arrow = NULL) + theme(legend.position = "none") 
            p <- p + geom_segment(data = statistics, aes(x=std_circ_up_lim, y=0, xend=std_circ_up_lim, yend=1), size = 1.5, color="red",linetype = "dashed", arrow = NULL) + theme(legend.position = "none") 
        }
        if (input$ci_method == "angular standard deviation") {
            p <- p + geom_segment(data = statistics, aes(x=std_ang_low_lim, y=0, xend=std_ang_low_lim, yend=1), size = 1.5, color="red",linetype = "dashed", arrow = NULL) + theme(legend.position = "none") 
            p <- p + geom_segment(data = statistics, aes(x=std_ang_up_lim, y=0, xend=std_ang_up_lim, yend=1), size = 1.5, color="red",linetype = "dashed", arrow = NULL) + theme(legend.position = "none") 
        }
    }

    
    if (input$hemi_rose_options == "all") {
    
      	statistics[1,"mean"] = statistics[1,"mean"] + 180.0
      	statistics[1,"ci_95_lower_limit"] = statistics[1,"ci_95_lower_limit"] + 180.0
      	statistics[1,"ci_95_upper_limit"] = statistics[1,"ci_95_upper_limit"] + 180.0
        statistics[1,"ci_90_lower_limit"] = statistics[1,"ci_90_lower_limit"] + 180.0
      	statistics[1,"ci_90_upper_limit"] = statistics[1,"ci_90_upper_limit"] + 180.0
        statistics[1,"ci_50_lower_limit"] = statistics[1,"ci_50_lower_limit"] + 180.0
      	statistics[1,"ci_50_upper_limit"] = statistics[1,"ci_50_upper_limit"] + 180.0

        p <- p + geom_segment(data = statistics, aes(x=mean, y=0, xend=mean, yend=polarity_index, size = 1.5, color="red", lineend = "butt"), arrow = NULL) + theme(legend.position = "none") 
 
        if (input$ci_plot) {
            if (input$ci_method == "95% CI of the mean") {
                p <- p + geom_segment(data = statistics, aes(x=ci_95_lower_limit, y=0, xend=ci_95_lower_limit, yend=1), size = 1.5, color="red",linetype = "dashed", arrow = NULL) + theme(legend.position = "none") 
                p <- p + geom_segment(data = statistics, aes(x=ci_95_upper_limit, y=0, xend=ci_95_upper_limit, yend=1), size = 1.5, color="red",linetype = "dashed", arrow = NULL) + theme(legend.position = "none") 
            }
            if (input$ci_method == "90% CI of the mean") {
                p <- p + geom_segment(data = statistics, aes(x=ci_90_lower_limit, y=0, xend=ci_90_lower_limit, yend=1), size = 1.5, color="red",linetype = "dashed", arrow = NULL) + theme(legend.position = "none") 
                p <- p + geom_segment(data = statistics, aes(x=ci_90_upper_limit, y=0, xend=ci_90_upper_limit, yend=1), size = 1.5, color="red",linetype = "dashed", arrow = NULL) + theme(legend.position = "none") 
            }
            if (input$ci_method == "50% CI of the mean") {
                p <- p + geom_segment(data = statistics, aes(x=ci_50_lower_limit, y=0, xend=ci_50_lower_limit, yend=1), size = 1.5, color="red",linetype = "dashed", arrow = NULL) + theme(legend.position = "none") 
                p <- p + geom_segment(data = statistics, aes(x=ci_50_upper_limit, y=0, xend=ci_50_upper_limit, yend=1), size = 1.5, color="red",linetype = "dashed", arrow = NULL) + theme(legend.position = "none") 
            }
            if (input$ci_method == "circular standard deviation") {
                p <- p + geom_segment(data = statistics, aes(x=std_circ_low_lim, y=0, xend=std_circ_low_lim, yend=1), size = 1.5, color="red",linetype = "dashed", arrow = NULL) + theme(legend.position = "none") 
                p <- p + geom_segment(data = statistics, aes(x=std_circ_up_lim, y=0, xend=std_circ_up_lim, yend=1), size = 1.5, color="red",linetype = "dashed", arrow = NULL) + theme(legend.position = "none") 
            }
            if (input$ci_method == "angular standard deviation") {
                p <- p + geom_segment(data = statistics, aes(x=std_ang_low_lim, y=0, xend=std_ang_low_lim, yend=1), size = 1.5, color="red",linetype = "dashed", arrow = NULL) + theme(legend.position = "none") 
                p <- p + geom_segment(data = statistics, aes(x=std_ang_up_lim, y=0, xend=std_ang_up_lim, yend=1), size = 1.5, color="red",linetype = "dashed", arrow = NULL) + theme(legend.position = "none") 
            }


            #p <- p + geom_segment(data = statistics, aes(x=ci_lower_limit, y=0, xend=ci_lower_limit, yend=1), size = 1.5, color="red",linetype = "dashed", arrow = NULL) + theme(legend.position = "none") 
            #p <- p + geom_segment(data = statistics, aes(x=upper_limit, y=0, xend=upper_limit, yend=1), size = 1.5, color="red",linetype = "dashed", arrow = NULL) + theme(legend.position = "none") 
        }


    }

    return(p)
  
}

compare_plot_2_axial <- function(parameters, input, statistics, feature_circular, plot_title, text_size = 24) {
  
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
 
 
    p <- p + geom_segment(data = statistics, aes(x=mean, y=0, xend=mean, yend=polarity_index, size = 1.5, color="red", lineend = "butt"), arrow = NULL) + theme(legend.position = "none") 
  return(p)
  
}

linear_histogram <- function(parameters, input, statistics, feature_linear, plot_title, text_size = 24) {
  
    
    bin_size = max(feature_linear)/input$bins
    #plot_title <- parameters[input$feature_select][[1]][3]
  
    p <- ggplot()

    if (input$histogram_plot) {
        p <- p + geom_histogram(aes(x = feature_linear, y = ..density..),
                   breaks = seq(0, max(feature_linear), bin_size),
                   colour = "black",
                   fill = "black",
                   alpha = 0.5)
    }

    if (input$kde_plot) {
        p <- p + geom_density(aes(x = feature_linear, y = ..density..),
                   colour = "black",
                   fill = "black",
                   alpha = 0.5)

    }

    p <- p + ggtitle(plot_title) +
        theme(axis.text.x = element_text(size = 18)) +
        theme_minimal(base_size = text_size) +
        xlab(sprintf("number of cells = : %s \n condition: %s", length(feature_linear), input$exp_condition)) 
 
    p <- p + geom_vline(data = statistics, aes(xintercept=mean),  col="red", size = 2) 
    return(p)
  
}




linear_histogram <- function(parameters, input, statistics, feature_linear, plot_title, text_size = 24) {
  
    
    bin_size = max(feature_linear)/input$bins
    #plot_title <- parameters[input$feature_select][[1]][3]
  
    p <- ggplot()

    if (input$histogram_plot) {
        p <- p + geom_histogram(aes(x = feature_linear, y = ..density..),
                   breaks = seq(0, max(feature_linear), bin_size),
                   colour = "black",
                   fill = "black",
                   alpha = 0.5)
    }

    if (input$kde_plot) {
        p <- p + geom_density(aes(x = feature_linear, y = ..density..),
                   colour = "black",
                   fill = "black",
                   alpha = 0.5)

    }

    p <- p + ggtitle(plot_title) +
        theme(axis.text.x = element_text(size = 18)) +
        theme_minimal(base_size = text_size) +
        xlab(sprintf("number of cells = : %s \n condition: %s", length(feature_linear), input$exp_condition)) 
 
    p <- p + geom_vline(data = statistics, aes(xintercept=mean),  col="red", size = 2) 
    return(p)
  
}


compare_plot_linear <- function(parameters, input, statistics, feature_linear_1, feature_linear_2, plot_title, text_size = 24) {
    
    #plot_title <- parameters[input$feature_select][[1]][3]
  
    p <- ggplot()
    max_1 <- max(feature_linear_1)
    max_2 <- max(feature_linear_2)
 
    #df_1 <- as.data
   
    bin_size = max(max_1,max_2)/input$bins

    if (input$histogram_plot) {
        p <- p + geom_histogram(aes(x = feature_linear_1, y = ..density..),
                   breaks = seq(0, max(max_1,max_2), bin_size),
                   colour = "black",
                   fill = "black",
                   alpha = 0.5) +
             geom_histogram(aes(x = feature_linear_2, y = ..density..),
                   breaks = seq(0, max(max_1,max_2), bin_size),
                   colour = "black",
                   fill = "black",
                   alpha = 0.5)
    }

    if (input$kde_plot) {
        p <- p + geom_density(aes(x = feature_linear_1, y = ..density..),
                   colour = "black",
                   fill = "black",
                   alpha = 0.5) +
            geom_density(aes(x = feature_linear_1, y = ..density..),
                   colour = "black",
                   fill = "black",
                   alpha = 0.5)
    }

    p <- p + ggtitle(plot_title) +
        theme(axis.text.x = element_text(size = 18)) +
        theme_minimal(base_size = text_size) +
        xlab(sprintf("number of cells = : %s \n condition: %s", length(feature_linear_1), input$exp_condition)) 
 
    p <- p + geom_vline(data = statistics, aes(xintercept=mean),  col="red", size = 2) 
    return(p)
  
}


