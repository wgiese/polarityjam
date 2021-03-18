

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
  polarity_index2 <- sin_mean*sin_mean + cos_mean*cos_mean
  
  return(sqrt(polarity_index2))
}