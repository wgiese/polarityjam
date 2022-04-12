

compute_circular_statistics <- function(data, feature, parameters) {

    circular_data <- unlist(data[feature])
  
    sin_sum <- 0.0
    cos_sum <- 0.0
    polarity_index <- 0.0
    std_circular <- 0.0    
    std_angular <- 0.0    

    for (i in 1:length(circular_data)) {
        angle <- circular_data[i]
        sin_sum <- sin_sum + sin(angle)
        cos_sum <- cos_sum + cos(angle)
    }
    sin_mean <- sin_sum/length(circular_data)
    cos_mean <- cos_sum/length(circular_data)
    polarity_index <- sqrt(sin_mean*sin_mean + cos_mean*cos_mean)
    std_angular <- sqrt(2.0*(1.0-polarity_index))*180.0/pi
    std_circular <- sqrt(-2.0*log(polarity_index))*180.0/pi


    print("Polarity index")
    print(polarity_index)    
    print("STD from circular")
    print(std_angular)
    print(std_circular)

    angle_mean_rad <- atan2(sin_mean, cos_mean)
    angle_mean_deg <- angle_mean_rad*180.0/pi
    if (angle_mean_rad < 0.0) {
        angle_mean_deg <- 360.0 + angle_mean_rad*180.0/pi
    }
    
    print("STD:")
    std_circ_up_lim <- angle_mean_deg + std_circular
    if ( std_circ_up_lim > 360.0)
        std_circ_up_lim <- std_circ_up_lim - 360.0    
    print(std_circ_up_lim)
    std_circ_low_lim <- angle_mean_deg - std_circular
    if ( std_circ_low_lim < 0.0)
        std_circ_low_lim <- std_circ_low_lim + 360.0    
    print(std_circ_low_lim)
    
    std_ang_up_lim <- angle_mean_deg + std_angular
    if ( std_ang_up_lim > 360.0)
        std_ang_up_lim <- std_ang_up_lim - 360.0    
    std_ang_low_lim <- angle_mean_deg - std_angular
    if ( std_ang_low_lim < 0.0)
        std_ang_low_lim <- std_ang_low_lim + 360.0    
 
    
    rayleigh_test_res <- r.test(circular_data)
    #rayleigh_test_res <- r.test(results_df$angle_deg, degree = TRUE)
    watson_res <- capture.output(watson.test(circular_data, alpha = 0, dist = "vonmises"))
    #v_test_res <- v0.test(circular_data, mu0 = pi)
    rao_res <- capture.output(rao.spacing.test(circular_data, alpha = 0))     

    v_test_res <- v0.test(circular_data, mu0 = pi)
    #rayleigh_test_mu_res <- v0.test(results_df$angle_deg, mu0 = 180.0, degree = TRUE)
    if (input$stats_method %in% c('V-Test')) {
        v_test_res <- v0.test(circular_data, mu0 = pi*input$cond_mean_direction/180.0)
        #av_test_res <- v0.test(circular_data, mu0 = pi*input$cond_mean_direction/180.0)
    }


    rayleigh_test <- rayleigh_test_res$p.value
    v_test <- v_test_res$p.value
   # rayleigh_test_mu <- rayleigh_test_mu_res$p.value

    ci_95_res <- vm.bootstrap.ci(circular_data, alpha = 0.05)
    #print("Confidence interval:")
    #print(ci_res$mu.ci)
    #print(str(ci_res$mu.ci))
    #print(ci_res$mu.ci[[1]])
    ci_95_lower_limit <- transform_rad_degrees(ci_95_res$mu.ci[[1]],-pi,pi,0.0,360.0)
    ci_95_upper_limit <- transform_rad_degrees(ci_95_res$mu.ci[[2]],-pi,pi,0.0,360.0)


    ci_90_res <- vm.bootstrap.ci(circular_data, alpha = 0.1)
    ci_90_lower_limit <- transform_rad_degrees(ci_90_res$mu.ci[[1]],-pi,pi,0.0,360.0)
    ci_90_upper_limit <- transform_rad_degrees(ci_90_res$mu.ci[[2]],-pi,pi,0.0,360.0)

    ci_50_res <- vm.bootstrap.ci(circular_data, alpha = 0.5)
    ci_50_lower_limit <- transform_rad_degrees(ci_50_res$mu.ci[[1]],-pi,pi,0.0,360.0)
    ci_50_upper_limit <- transform_rad_degrees(ci_50_res$mu.ci[[2]],-pi,pi,0.0,360.0)

#against_flow <- polarity_data[(polarity_data 150*pi/180.0),]
    #against_flow <- against_flow [(against_flow < 210*pi/180.0),]
    #with_flow <- polarity_data[((polarity_data > 330*pi/180.0) | (polarity_data < 30*pi/180.0)),]
  
    #signed_polarity_index <- (nrow(against_flow) - nrow(with_flow))/(nrow(against_flow) + nrow(with_flow))

    #values <- c(polarity_index, sin_mean, cos_mean, angle_mean_deg)
    #names(values) <- c("polarity_index", "sin_mean", "cos_mean", "angle_mean_deg")
    
    #values <- data.frame("polarity_index" = polarity_index, "signed_polarity_index" = signed_polarity_index, "sin_mean" = sin_mean, "cos_mean" = cos_mean, "angle_mean_deg" = angle_mean_deg )
  
    values <- data.frame( "polarity_index" = polarity_index, 
                          "mean" = angle_mean_deg,
                          "std_angular" = std_angular,
                          "std_circ_up_lim" = std_circ_up_lim,
                          "std_circ_low_lim" = std_circ_low_lim,
                          "std_circular" = std_circular,
                          "std_ang_up_lim" = std_ang_up_lim,
                          "std_ang_low_lim" = std_ang_low_lim,
                          "rayleigh_test" = rayleigh_test,
                          "v_test" = v_test,
                          "watson_test" = watson_res[5], 
                          "rao_test" = rao_res[5], 
                          "ci_95_lower_limit" = ci_95_lower_limit,
                          "ci_95_upper_limit" = ci_95_upper_limit,
                          "ci_90_lower_limit" = ci_90_lower_limit,
                          "ci_90_upper_limit" = ci_90_upper_limit,
                          "ci_50_lower_limit" = ci_50_lower_limit,
                          "ci_50_upper_limit" = ci_50_upper_limit
    )
  
    return(values)
}

transform_rad_degrees <- function(value, source_low, source_high, target_low, target_high) {
    
    res <- 180.0*value/pi
    if (value < target_low) {
        res <- 180.0*value/pi + target_high
    }    
    if (value > target_high) {
        res <- 180.0*value/pi - target_high
    }    

    return(res)

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
    std_angular <- sqrt(2.0*(1.0-polarity_index))*180.0/pi
    std_circular <- sqrt(-2.0*log(polarity_index))*180.0/pi    
    angle_mean_rad <- atan2(sin_mean, cos_mean)/2.0
    angle_mean_deg <- angle_mean_rad*180.0/3.1415926
    if (angle_mean_rad < 0.0) {
        angle_mean_deg <- 180.0 + angle_mean_rad*180.0/pi
    }
    
    

    print("mean")
    print(angle_mean_deg)
    rayleigh_test_res <- r.test(circular_data)
    #rayleigh_test_res <- r.test(results_df$angle_deg, degree = TRUE)
    #rayleigh_test_mu_res <- v0.test(circular_data, mu0 = pi)
    #rayleigh_test_mu_res <- v0.test(results_df$angle_deg, mu0 = 180.0, degree = TRUE)
    
    rayleigh_test <- rayleigh_test_res$p.value
    rayleigh_test_mu_res <- v0.test(circular_data, mu0 = pi)
    #rayleigh_test_mu <- rayleigh_test_mu_res$p.value
    ci_res <- vm.bootstrap.ci(circular_data)
    print("Confidence interval:")
    print(ci_res$mu.ci)
    print(str(ci_res$mu.ci))
    
    ci_95_res <- vm.bootstrap.ci(circular_data, alpha = 0.05)
    ci_95_lower_limit <- transform_rad_degrees(ci_95_res$mu.ci[[1]]/2.0,-pi/2.0,pi/2.0,0.0,180.0)
    ci_95_upper_limit <- transform_rad_degrees(ci_95_res$mu.ci[[2]]/2.0,-pi/2.0,pi/2.0,0.0,180.0)

    ci_90_res <- vm.bootstrap.ci(circular_data, alpha = 0.1)
    ci_90_lower_limit <- transform_rad_degrees(ci_90_res$mu.ci[[1]]/2.0,-pi/2.0,pi/2.0,0.0,180.0)
    ci_90_upper_limit <- transform_rad_degrees(ci_90_res$mu.ci[[2]]/2.0,-pi/2.0,pi/2.0,0.0,180.0)

    ci_50_res <- vm.bootstrap.ci(circular_data, alpha = 0.5)
    ci_50_lower_limit <- transform_rad_degrees(ci_50_res$mu.ci[[1]]/2.0,-pi/2.0,pi/2.0,0.0,180.0)
    ci_50_upper_limit <- transform_rad_degrees(ci_50_res$mu.ci[[2]]/2.0,-pi/2.0,pi/2.0,0.0,180.0)




    ci_lower_limit <- 90.0*ci_res$mu.ci[[1]]/pi
    ci_upper_limit <- 90.0*ci_res$mu.ci[[2]]/pi

    if (ci_res$mu.ci[[1]] < 0.0) {
        ci_lower_limit <- 180.0 + 90.0*ci_res$mu.ci[[1]]/pi
    }

    if (ci_res$mu.ci[[2]] < 0.0) {
        ci_upper_limit <- 180.0 + 90.0*ci_res$mu.ci[[2]]/pi
    }

    print("STD:")
    std_circ_up_lim <- angle_mean_deg + std_circular
    if ( std_circ_up_lim > 180.0)
        std_circ_up_lim <- std_circ_up_lim - 180.0    
    print(std_circ_up_lim)
    std_circ_low_lim <- angle_mean_deg - std_circular
    if ( std_circ_low_lim < 0.0)
        std_circ_low_lim <- std_circ_low_lim + 180.0    
    print(std_circ_low_lim)
    
    std_ang_up_lim <- angle_mean_deg + std_angular
    if ( std_ang_up_lim > 180.0)
        std_ang_up_lim <- std_ang_up_lim - 180.0    
    std_ang_low_lim <- angle_mean_deg - std_angular
    if ( std_ang_low_lim < 0.0)
        std_ang_low_lim <- std_ang_low_lim + 180.0

 
  #signed_polarity_index <- (nrow(against_flow) - nrow(with_flow))/(nrow(against_flow) + nrow(with_flow))

  #values <- c(polarity_index, sin_mean, cos_mean, angle_mean_deg)
  #names(values) <- c("polarity_index", "sin_mean", "cos_mean", "angle_mean_deg")
  
    
  #values <- data.frame("polarity_index" = polarity_index, "signed_polarity_index" = signed_polarity_index, "sin_mean" = sin_mean, "cos_mean" = cos_mean, "angle_mean_deg" = angle_mean_deg )
  
  values <- data.frame( "polarity_index" = polarity_index, 
                        "rayleigh_test" = rayleigh_test,
                        "mean" = angle_mean_deg,
                        "std_angular" = std_angular,
                        "std_circ_up_lim" = std_circ_up_lim,
                        "std_circ_low_lim" = std_circ_low_lim,
                        "std_circular" = std_circular,
                        "std_ang_up_lim" = std_ang_up_lim,
                        "std_ang_low_lim" = std_ang_low_lim,
                        "ci_95_lower_limit" = ci_95_lower_limit,
                        "ci_95_upper_limit" = ci_95_upper_limit,
                        "ci_90_lower_limit" = ci_90_lower_limit,
                        "ci_90_upper_limit" = ci_90_upper_limit,
                        "ci_50_lower_limit" = ci_50_lower_limit,
                        "ci_50_upper_limit" = ci_50_upper_limit
                    )
  return(values)
}


compute_linear_statistics <- function(data, feature, parameters) {

    data <- unlist(data[feature])

    mean_ <- mean(data)
    std_ <- sd(data)
    median_ <- median(data)

    values <- data.frame("mean" = mean_,  "std" = std_, "median" = median_ )

}

