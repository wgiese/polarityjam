plot_circular_circular <- function(correlation_data, input, parameters, plot_nr = 0, text_size = 24) {
  
  feature_1 <- parameters[input$feature_select_1][[1]][1]
  feature_2 <- parameters[input$feature_select_2][[1]][1]
  feature_1_values <- unlist(correlation_data[feature_1])
  feature_1_values_ <- correlation_data[feature_1] * 180.0 / pi
  feature_2_values <- unlist(correlation_data[feature_2])
  feature_2_values_ <- correlation_data[feature_2] * 180.0 / pi
  
  feature_1_name <- parameters[input$feature_select_1][[1]][3]
  feature_2_name <- parameters[input$feature_select_2][[1]][3]
  
  mode_1 <- parameters[input$feature_select_1][[1]][2]
  mode_2 <- parameters[input$feature_select_2][[1]][2]
  
  conditions <- correlation_data[input$condition_col]
  
  mean_dir_1 <- compute_mean_circular(feature_1_values, mode_1)
  mean_dir_2 <- compute_mean_circular(feature_2_values, mode_2)
  
  res <- compute_correlation(feature_1_values, mode_1, feature_2_values, mode_2) 
  
  print("Mean directions")
  print(mean_dir_1)
  print(mean_dir_2)
  
  
  if (input$center_corr_plot == TRUE) {
  
    if (mode_1 == "directional") {
      if ((mean_dir_1 < pi / 2.0) | (mean_dir_1 > 3.0 * pi / 2.0)) {
        for (i in 1:length(feature_2_values)) {
          if (feature_1_values[i] > pi) {
            correlation_data[i, feature_1] <- correlation_data[i, feature_1] - 2.0 * pi
          }
        }
        feature_1_values_ <- correlation_data[feature_1] * 180.0 / pi
      }
    } else if (mode_1 == "undirectional") {
      if ((mean_dir_1 < pi / 4.0) | (mean_dir_1 > 3.0 * pi / 4.0)) {
        for (i in 1:length(feature_2_values)) {
          if (feature_1_values[i] > pi/2.0) {
            correlation_data[i, feature_1] <- correlation_data[i, feature_1] -  pi
          }
        }
        feature_1_values_ <- correlation_data[feature_1] * 180.0 / pi
      }
    }
    
    if (mode_1 == "directional") {
      if ((mean_dir_2 < pi / 2.0) | (mean_dir_2 > 3.0 * pi / 2.0)) {
        for (i in 1:length(feature_2_values)) {
          if (feature_2_values[i] > pi) {
            correlation_data[i, feature_2] <- correlation_data[i, feature_2] - 2.0 * pi
          }
        }
        feature_2_values_ <- correlation_data[feature_2] * 180.0 / pi
      }
    } else if  (mode_2 == "undirectional") {
      if ((mean_dir_2 < pi / 4.0) | (mean_dir_2 > 3.0 * pi / 4.0)) {
        for (i in 1:length(feature_2_values)) {
          if (feature_2_values[i] > pi/2.0) {
            correlation_data[i, feature_2] <- correlation_data[i, feature_2] - pi
          }
        }
        feature_2_values_ <- correlation_data[feature_2] * 180.0 / pi
      }
    }
  }
  
  reg_coeff <- signif(res$r, digits = 3)
  
  plot_df <- as.data.frame(c(feature_1_values_, feature_2_values_, conditions))
  
  p_value_ <- signif(res$p.value, digits = 3)
  
  if (p_value_ < 0.001) {
    p_value <- "P < 0.001"
  } else if (p_value_ < 0.01) {
    p_value <- "P < 0.01"
  } else {
    p_value <- p_value_
  }
  
  colnames(plot_df) <- c("x", "y", "condition")
  p <- ggplot(plot_df, aes(x = x, y = y, color = condition)) +
    geom_point(size = input$marker_size_corr) +
    theme_minimal(base_size = text_size) # theme_bw()
  p <- p + theme(aspect.ratio = 3 / 3)
  p <- p + ggtitle(sprintf("number of cells = : %s \n r = %s, p-value: %s", length(feature_1_values), reg_coeff, p_value))
  p <- p + xlab(feature_1_name) + ylab(feature_2_name)
  
}

compute_mean_circular <- function(feature_values, mode = "directional") {
  mean_dir <- 0.0
  if ( mode == "directional") {
    mean_dir <- circ.mean(feature_values)
    
    if (mean_dir < 0.0) {
      mean_dir <- mean_dir + 2.0 * pi
    }
    
  }
  else {
    feature_values_ <- feature_values
    
    for (i in 1:length(feature_values)) {
        feature_values_[i] <- 2.0 *feature_values[i]
    }
    
    mean_dir <- circ.mean(feature_values_)/2.0
    
    if (mean_dir < 0.0) {
      mean_dir <- mean_dir + pi
    }
    
  }
  
  return(mean_dir)
}


compute_correlation <- function(feature_1_values, mode_1 = "directional", feature_2_values, mode_2 = "directional") {
  
  feature_1_values_ <- feature_1_values
  feature_2_values_ <- feature_2_values
  
  if ( mode_1 == "undirectional") {
    for (i in 1:length(feature_1_values)) {
      feature_1_values_[i] <- 2.0 *feature_1_values[i]
    }
  }
  
  if ( mode_2 == "undirectional") {
    for (i in 1:length(feature_2_values)) {
      feature_2_values_[i] <- 2.0 *feature_2_values[i]
    }
  }
  
  res <- NULL
  
  if ((mode_1 != "linear") && (mode_2 != "linear")) {
    
    res <- circ.cor(feature_1_values_, feature_2_values_, test = TRUE)
  
  }
  
  #if ((mode_1 == "linear") && (mode_2 != "linear")) {
  #  
  #  r_sx <- cor(method= "pearson")
  #
  #}
  
  if ((mode_1 == "linear") && (mode_2 == "linear")) {
    
    res <- cor(feature_1_values_, feature_2_values_, method= "pearson")
    
  }
  #res <- circ.cor(feature_1_values, feature_2_values, test = TRUE)
  
  
  
  return(res)
  
}
  
  
  