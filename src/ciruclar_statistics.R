

compute_polarity_index <- function(polarity_data) {
  print(polarity_data$angle_deg)
  
  sin_sum <- 0.0
  cos_sum <- 0.0
  #polarity_index <- 0.0
  for (angle in polarity_data$angle_rad) {
    sin_sum <- sin_sum + sin(angle)
    cos_sum <- cos_sum + cos(angle)
  }
  sin_mean <- sin_sum/length(polarity_data$angle_rad)
  cos_mean <- cos_sum/length(polarity_data$angle_rad)
  polarity_index <- sqrt(sin_mean*sin_mean + cos_mean*cos_mean)
  angle_mean_rad <- atan2(sin_mean, cos_mean)
  angle_mean_deg <- angle_mean_rad*180.0/3.1415926
  if (angle_mean_rad < 0.0) {
    angle_mean_deg <- 360.0 + angle_mean_rad*180.0/3.1415926
  }
    

  #values <- c(polarity_index, sin_mean, cos_mean, angle_mean_deg)
  #names(values) <- c("polarity_index", "sin_mean", "cos_mean", "angle_mean_deg")
  
  values <- data.frame("polarity_index" = polarity_index, "sin_mean" = sin_mean, "cos_mean" = cos_mean, "angle_mean_deg" = angle_mean_deg )
  
  return(values)
}
