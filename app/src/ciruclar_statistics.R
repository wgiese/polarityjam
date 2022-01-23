

compute_circular_statistics <- function(data, feature, parameters) {

    circular_data <- unlist(data[feature])
  
    sin_sum <- 0.0
    cos_sum <- 0.0
    polarity_index <- 0.0
    
    for (i in 1:length(circular_data)) {
        angle <- circular_data[i]
        sin_sum <- sin_sum + sin(angle)
        cos_sum <- cos_sum + cos(angle)
    }
    sin_mean <- sin_sum/length(circular_data)
    cos_mean <- cos_sum/length(circular_data)
    polarity_index <- sqrt(sin_mean*sin_mean + cos_mean*cos_mean)
    angle_mean_rad <- atan2(sin_mean, cos_mean)
    angle_mean_deg <- angle_mean_rad*180.0/3.1415926
    if (angle_mean_rad < 0.0) {
        angle_mean_deg <- 360.0 + angle_mean_rad*180.0/3.1415926
    }

    
    rayleigh_test_res <- r.test(circular_data)
    #rayleigh_test_res <- r.test(results_df$angle_deg, degree = TRUE)
    rayleigh_test_mu_res <- v0.test(circular_data, mu0 = pi)
    #rayleigh_test_mu_res <- v0.test(results_df$angle_deg, mu0 = 180.0, degree = TRUE)
    
    rayleigh_test <- rayleigh_test_res$p.value
    rayleigh_test_mu <- rayleigh_test_mu_res$p.value

    ci_res <- vm.bootstrap.ci(circular_data)
    print("Confidence interval:")
    print(ci_res$mu.ci)
    print(str(ci_res$mu.ci))
    print(ci_res$mu.ci[[1]])
    lower_limit <- 180.0*ci_res$mu.ci[[1]]/3.1415926
    upper_limit <- 180.0*ci_res$mu.ci[[2]]/3.1415926

    print(ci_res$mu.ci[[2]])
    #against_flow <- polarity_data[(polarity_data 150*pi/180.0),]
    #against_flow <- against_flow [(against_flow < 210*pi/180.0),]
    #with_flow <- polarity_data[((polarity_data > 330*pi/180.0) | (polarity_data < 30*pi/180.0)),]
  
    #signed_polarity_index <- (nrow(against_flow) - nrow(with_flow))/(nrow(against_flow) + nrow(with_flow))

    #values <- c(polarity_index, sin_mean, cos_mean, angle_mean_deg)
    #names(values) <- c("polarity_index", "sin_mean", "cos_mean", "angle_mean_deg")
    
    #values <- data.frame("polarity_index" = polarity_index, "signed_polarity_index" = signed_polarity_index, "sin_mean" = sin_mean, "cos_mean" = cos_mean, "angle_mean_deg" = angle_mean_deg )
  
    values <- data.frame( "polarity_index" = polarity_index, 
                          "mean" = angle_mean_deg,
                          "rayleigh_test" = rayleigh_test,
                          "rayleigh_test_mu" = rayleigh_test_mu,
                          "lower_limit" = lower_limit,
                          "upper_limit" = upper_limit)
  
    return(values)
}

comparison_circular_statistics <- function(data_1, data_2, feature, parameters) {
  
  circular_data_1 <- unlist(data_1[feature])
  circular_data_2 <- unlist(data_2[feature])
  
  values_1 <- compute_circular_statistics(data_1, feature, parameters)
  values_2 <- compute_circular_statistics(data_2, feature, parameters)
  
  #values <- data.frame( "polarity_index" = polarity_index, 
  #                      "mean" = angle_mean_deg,
  #                      "rayleigh_test" = rayleigh_test,
  #                      "rayleigh_test_mu" = rayleigh_test_mu)
  
  return(values_1)
}


compute_2_axial_statistics <- function(data, feature, parameters) {
  
    p_axial_data <- unlist(data[feature])
    circular_data <- unlist(data[feature])
  
    sin_sum <- 0.0
    cos_sum <- 0.0
    polarity_index <- 0.0
    
    for (i in 1:length(circular_data)) {
        circular_data[i] <- 2.0*p_axial_data[i]
        angle <- circular_data[i]
        sin_sum <- sin_sum + sin(angle)
        cos_sum <- cos_sum + cos(angle)
    }
    sin_mean <- sin_sum/length(circular_data)
    cos_mean <- cos_sum/length(circular_data)
    
    polarity_index <- sqrt(sin_mean*sin_mean + cos_mean*cos_mean)
    angle_mean_rad <- atan2(sin_mean, cos_mean)/2.0
    angle_mean_deg <- angle_mean_rad*180.0/3.1415926
    if (angle_mean_rad < 0.0) {
        angle_mean_deg <- 180.0 + angle_mean_rad*180.0/3.1415926
    }

    rayleigh_test_res <- r.test(circular_data)
    #rayleigh_test_res <- r.test(results_df$angle_deg, degree = TRUE)
    #rayleigh_test_mu_res <- v0.test(circular_data, mu0 = pi)
    #rayleigh_test_mu_res <- v0.test(results_df$angle_deg, mu0 = 180.0, degree = TRUE)
    
    rayleigh_test <- rayleigh_test_res$p.value
    #rayleigh_test_mu <- rayleigh_test_mu_res$p.value
    ci_res <- vm.bootstrap.ci(circular_data)
    

    
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
  
  values <- data.frame( "polarity_index" = polarity_index, 
                        "rayleigh_test" = rayleigh_test,
                        "mean" = angle_mean_deg )

  return(values)
}


compute_linear_statistics <- function(data, feature, parameters) {

    data <- unlist(data[feature])

    mean_ <- mean(data)
    std_ <- sd(data)
    median_ <- median(data)

    values <- data.frame("mean" = mean_,  "std" = std_, "median" = median_ )

}

