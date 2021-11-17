

compute_polarity_index <- function(polarity_data) {
  print(polarity_data)
  
  sin_sum <- 0.0
  cos_sum <- 0.0
  #polarity_index <- 0.0
  for (i in 1:length(polarity_data)) {
    angle <- polarity_data[i]
    sin_sum <- sin_sum + sin(angle)
    cos_sum <- cos_sum + cos(angle)
  }
  sin_mean <- sin_sum/length(polarity_data)
  cos_mean <- cos_sum/length(polarity_data)
  polarity_index <- sqrt(sin_mean*sin_mean + cos_mean*cos_mean)
  angle_mean_rad <- atan2(sin_mean, cos_mean)
  angle_mean_deg <- angle_mean_rad*180.0/3.1415926
  if (angle_mean_rad < 0.0) {
    angle_mean_deg <- 360.0 + angle_mean_rad*180.0/3.1415926
  }
    
  #against_flow <- polarity_data[(polarity_data 150*pi/180.0),]
  #against_flow <- against_flow [(against_flow < 210*pi/180.0),]
  #with_flow <- polarity_data[((polarity_data > 330*pi/180.0) | (polarity_data < 30*pi/180.0)),]
  
  #print("Delme!")
  #print(against_flow)
  #print(with_flow)
  #print("number of rows:")
  #print(nrow(polarity_data))
  #print(nrow(against_flow))
  #print(nrow(with_flow))
  
  #signed_polarity_index <- (nrow(against_flow) - nrow(with_flow))/(nrow(against_flow) + nrow(with_flow))

  #values <- c(polarity_index, sin_mean, cos_mean, angle_mean_deg)
  #names(values) <- c("polarity_index", "sin_mean", "cos_mean", "angle_mean_deg")
  
    
  #values <- data.frame("polarity_index" = polarity_index, "signed_polarity_index" = signed_polarity_index, "sin_mean" = sin_mean, "cos_mean" = cos_mean, "angle_mean_deg" = angle_mean_deg )
  
  values <- data.frame("polarity_index" = polarity_index,  "sin_mean" = sin_mean, "cos_mean" = cos_mean, "angle_mean_deg" = angle_mean_deg )
  
  return(values)
}

