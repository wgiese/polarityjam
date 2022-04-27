# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Polarity JaM: Shiny app for plotting and comparing polarity data (beta 0.2)
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Takes spreadsheet type data as input with circular and non-circular features
# Visualization of circular and non-circular distributions
# Comparative circular statistics 
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# MIT License
#
# Copyright (c) 2021 Wolfgang Giese
# electronic mail address: wolfgang #dot# giese #at# mdc #minus# berlin #dot# de
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


options(shiny.maxRequestSize = 30*1024^2)

library(shiny)
library(shinyFiles)
library(circular)
#library(CircMLE)
library(ggplot2)
library(shape)
library(shinyWidgets)
library(tools)
library(grid)
library(lattice)
library(gridExtra)
library(FNN)
library(tidyverse)
library(CircStats)
library(readxl)
library(fs)
library(rjson)


#From Paul Tol: https://personal.sron.nl/~pault/
Tol_bright <- c('#EE6677', '#228833', '#4477AA', '#CCBB44', '#66CCEE', '#AA3377', '#BBBBBB')

Tol_muted <- c('#88CCEE', '#44AA99', '#117733', '#332288', '#DDCC77', '#999933','#CC6677', '#882255', '#AA4499', '#DDDDDD')

Tol_light <- c('#BBCC33', '#AAAA00', '#77AADD', '#EE8866', '#EEDD88', '#FFAABB', '#99DDFF', '#44BB99', '#DDDDDD')

#From Color Universal Design (CUD): https://jfly.uni-koeln.de/color/
Okabe_Ito <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")


# Create a reactive object here that we can share between all the sessions.
vals <- reactiveValues(count=0)

###### UI: User interface #########

ui <- navbarPage("Polarity JaM - a web app for visualizing cell polarity, junction and morphology data (beta 0.2)",

### Panel 0: Data preparation

    tabPanel("Data preparation",
        sidebarLayout(
            sidebarPanel(
                shinyDirButton("dir", "Input directory", "Upload"),
                verbatimTextOutput("dir", placeholder = TRUE),
                actionButton("refreshStack", "Refresh"),
                selectInput("dataset_merged", "Choose a dataset:",
                            choices = c("merged_file")),
                downloadButton("downloadProcessedData", "Download")
            ),
 
            #conditionalPanel(condition = "input.tidyInput==true",
            #    selectInput("x_var", "Select variable for x-axis", choices = ""),
            #    selectInput("y_var", "Select variable for y-axis", choices = ""),
            #    selectInput("g_var", "Identifier of samples", choices = ""),
            #    selectInput("c_var", "Identifier of conditions", choices = ""),
            #    selectInput("filter_column", "Filter based on this parameter:", choices = ""),
            #    selectInput("remove_these_conditions", "Deselect these conditions:", "", multiple = TRUE)
            #),
               
            # Show a plot of the generated distribution
            mainPanel(
                tabsetPanel(
                    tabPanel("Table", tableOutput("merged_stack"))
                )
            )
        )
    ),
 
                 
### Panel A: Image stack histogram
                 
    tabPanel("Image stack analysis",
        sidebarLayout(
            sidebarPanel(
                fileInput("stackData", "Upload data file",
                            accept = c( "text/csv",
                            "text/comma-separated-values,text/plain",
                            ".csv",".xlsx")),
                            tags$hr(),
                checkboxInput("header_correlation", "File upload", TRUE),
#                shinyDirButton("dir", "Input directory", "Upload"),
#                verbatimTextOutput("dir", placeholder = TRUE),
#                actionButton("refreshStack", "Refresh"),
                sliderInput("bins",
                            "Number of bins:",
                            min = 4,
                            max = 36,
                            value = 12),
                sliderInput("min_eccentricity",
                            "Mininum eccentricity",
                            min = 0,
                            max = 1,
                            step = 0.1,
                            value = 0.0),
                sliderInput("min_nuclei_golgi_dist",
                            "Minimum nuclei golgi distance",
                            min = 0,
                            max = 10,
                            step = 1,
                            value = 0),
                selectInput("feature_select", "Choose a feature:",
                            choices = c("nuclei_golgi_polarity","major_axis_shape_orientation",
                            "major_axis_nucleus_orientation","eccentricity","major_over_minor_ratio",
                            "mean_expression","marker_polarity","area","perimeter")),
                selectInput("stats_method", "Choose a stats test", 
                            choices = c("Rayleigh uniform", "V-Test", "Rao's Test", "Watson's Test", "None")),
                conditionalPanel(
                   # condition = "input.stats_method %in% c('V-Test')",
                    condition = "input.stats_method == 'V-Test'",
                    numericInput("cond_mean_direction",
                            "Conditional mean direction", value = 180),
                    NULL
                ),    
                textInput("exp_condition", "Exp. condition", "condition A"),
 
                checkboxInput("ci_plot", "Confidence interval (CI)", TRUE),
                conditionalPanel(
                    condition = "input.ci_plot == true",
                    selectInput("ci_method", "CI method", 
                                choices = c("95% CI of the mean","90% CI of the mean","50% CI of the mean", 
                                "circular standard deviation", "angular standard deviation"))
                ),
                checkboxInput("kde_plot", "KDE plot", FALSE),
                checkboxInput("scatter_plot", "Scatter plot", FALSE),
                checkboxInput("histogram_plot", "Histogram plot", TRUE),
                checkboxInput("area_scaled", "area scaled histogram", TRUE),
                #checkboxInput("left_axial", "hemirose on left", FALSE),
                
                selectInput("plot_mode", "Choose data modality:",
                            choices = c("circular", "semicircular", "linear")),

                conditionalPanel(
                    condition = "input.plot_mode == 'semicircular'",
                    selectInput("hemi_rose_options", "Hemirose plot options:",
                    choices = c("up","down","left","right","all"))
                ),

                selectInput("dataset", "Choose a dataset:",
                            choices = c("statistics_file","merged_plot_file","multi_plot_file")),
                selectInput("image_file_format", "Choose image file format:",
                            choices = c(".pdf",".eps",".png")),
                downloadButton("downloadData", "Download")
            ),
                
            # Show a plot of the generated distribution
            mainPanel(
                tabsetPanel(
#                    tabPanel("Table", tableOutput("merged_stack")),
                    tabPanel("Plot", plotOutput("merged_plot", height = "860px")),
                    tabPanel("MultiPlot", plotOutput("multi_dist_plot", height = "860px")),
                    tabPanel("Statistics", tableOutput("merged_statistics"))
                )
            )
        )
    ),
                 

            
### Panel B: Correlation analysis

    tabPanel("Correlation analysis",
        sidebarLayout(
            sidebarPanel(
                fileInput("correlationData", "Upload data file",
                            accept = c( "text/csv",
                            "text/comma-separated-values,text/plain",
                            ".csv",".xlsx")),
                            tags$hr(),
                checkboxInput("header_correlation", "File upload", TRUE),
                selectInput("feature_select_1", "Choose a feature 1:",
                            choices = c("nuclei_golgi_polarity","major_axis_shape_orientation","major_axis_nucleus_orientation","eccentricity","mean_expression","area","perimeter")),
                selectInput("feature_select_2", "Choose a feature 2:",
                            choices = c("nuclei_golgi_polarity","major_axis_shape_orientation","major_axis_nucleus_orientation","eccentricity","mean_expression","area","perimeter")),
                selectInput("datasetSingleImage", "Download:",
                            choices = c("results_file","statistics_file","orientation_plot", "rose_histogram")),
                            tags$hr(),
                checkboxInput("header_image", "File upload", TRUE),
                downloadButton("downloadCorrelationData", "Download")
            ),
            mainPanel(
                tabsetPanel(
                tabPanel("Plot", plotOutput("correlation_plot", height = "1000px"))#,
                #tabPanel("Statistics", tableOutput("singleImageStatistics"))
                )
            )
        )
    ),
      
### Panel C: Comparison statistics

    tabPanel("Compare",                    
        sidebarLayout(
            sidebarPanel(
                fileInput("condition_1", "Condition 1",
                            accept = c( "text/csv",
                            "text/comma-separated-values,text/plain",
                            ".csv")),    
                tags$hr(),
                checkboxInput("header_cond1", "File upload", TRUE),
                fileInput("condition_2", "Condition 2",
                            accept = c( "text/csv",
                            "text/comma-separated-values,text/plain",
                            ".csv")), 
                tags$hr(),
                checkboxInput("header_cond2", "File upload", TRUE),
                sliderInput("bins_comparison",
                            "Number of bins:",
                            min = 1,
                            max = 30,
                            value = 12),
                selectInput("feature_comparison", "Choose a feature:",
                            choices = c("nuclei_golgi_polarity","major_axis_shape_orientation",
                            "major_axis_nucleus_orientation","eccentricity","major_over_minor_ratio",
                            "mean_expression","marker_polarity","area","perimeter")),
                checkboxInput("kde_comparison", "KDE plot", FALSE),
                checkboxInput("histogram_comparison", "Histogram plot", TRUE),
                checkboxInput("split_view_comparison", "Split view", TRUE),
            ),
            mainPanel(
                #tabPanel("Plot", plotOutput("comparison_plot", height = "1000px")),
                tabsetPanel(
                    tabPanel("Plot", plotOutput("comparison_plot", height = "1000px")),
                    tabPanel("CDF Plot", plotOutput("CDFPlot")), 
                    tabPanel("Statistics", tableOutput("comparison_statistics"))
                )
            )
        )
    )
)



loadVolumes <- function (exclude = "") 
{
  "This function returns all volumes in a named list.
  TODO: check if built-in function would also work"
  
  osSystem <- Sys.info()["sysname"]
  if (osSystem == "Darwin") {
    volumes <- dir_ls("/Volumes")
    names(volumes) <- basename(volumes)
  }
  else if (osSystem == "Linux") {
    volumes <- c(Computer = "/")
    media <- c(media = '/media/')
    home <- c(home = '~')
    #media <- dir_ls("/media")
    #names(media) <- basename(media)
    volumes <- c(volumes,media,home)
  }
  else if (osSystem == "Windows") {
    wmic <- paste0(Sys.getenv("SystemRoot"), "\\System32\\Wbem\\WMIC.exe")
    if (!file.exists(wmic)) {
      message("\nThe wmic program does not seem to be in the default location")
      message("Please report this problem and include output from the command")
      message("'where wmic' to https://github.com/thomasp85/shinyFiles/issues")
      volumes <- Sys.getenv("HOMEDRIVE")
      volNames <- ""
    }
    else {
      volumes <- system(paste(wmic, "logicaldisk get Caption"), 
                        intern = T)
      volumes <- sub(" *\\r$", "", volumes)
      keep <- !tolower(volumes) %in% c("caption", 
                                       "")
      volumes <- volumes[keep]
      volNames <- system(paste(wmic, "logicaldisk get VolumeName"), 
                         intern = T)
      volNames <- sub(" *\\r$", "", volNames)
      volNames <- volNames[keep]
      volNames <- paste0(volNames, ifelse(volNames == "", 
                                          "", " "))
    }
    volNames <- paste0(volNames, "(", volumes, ")")
    names(volumes) <- volNames
    volumes <- gsub(":$", ":/", volumes)
  }
  else {
    stop("unsupported OS")
  }
  if (!is.null(exclude)) {
    volumes <- volumes[!names(volumes) %in% exclude]
  }
  volumes
}

# Define server logic
server <- function(input, output, session) {
  
  
  ### Panel A
  
  # function to choose directory for data from image stack
  volumes <- loadVolumes()
  print("Volumes")
  print(volumes)
  shinyDirChoose(
    input,
    'dir',
    roots = volumes#,
  )
  
  dir <- reactive(input$dir)
  stack_data_info <- reactiveValues(datapath = getwd())
  
  output$dir <- renderText({
    print(stack_data_info$datapath)
    stack_data_info$datapath
  })
  
  observeEvent(ignoreNULL = TRUE,
               eventExpr = {
                 input$dir
               },
               handlerExpr = {
               if (!"path" %in% names(dir())) return()

                rootdir <- normalizePath(paste(volumes[input$dir[[2]]][[1]]))
                stack_data_info$datapath <- file.path(rootdir,paste(unlist(dir()$path[-1]), collapse = .Platform$file.sep))
                 
                 #home <- normalizePath("~")
                 #stack_data_info$datapath <-
                 #  file.path(paste(unlist(dir()$path[-1]), collapse = .Platform$file.sep))
 
               })
  
  

  mergedStack <- reactive({
    "
    reactive function that reads all csv files 
    from the directory given in stack_data_info$datapath 
    and combines them into one data frame
    "
   
    inFileStackData <- input$stackData

    if (is.null(inFileStackData)) {
        datapath <- stack_data_info$datapath 
        
        # get list of files in the datapath
        file_list <- list.files(datapath)
        print("File_list")
        print(file_list)

        counter <- 1
        plist <- list()
        results_all_df <-data.frame()
        tag <- FALSE
        
        for (file_name in file_list[1:length(file_list)]){

          if (file_ext(file_name) == "csv") {
            results_df <- read.csv(paste0(datapath,"/",file_name))
            results_df <- cbind(results_df,  data.frame("filename"=rep(file_name, nrow(results_df))) )
            results_all_df  <- rbind(results_all_df, results_df )
          }
        }
        
        if (length(results_all_df) > 1) {
          results_all_df$datapath <- datapath
          results_all_df$experimental_condition <- input$exp_condition
        }
    } else {    
        results_all_df <- read.csv(inFileStackData$datapath, header = input$header_correlation)
    }
    
    
    results_all_df
  })
  # end of merged stack function
  
  
  output$merged_stack <- renderTable({
    "
    function that the merged stack of polarity data and angles in table format
    "
    mergedStack()
  })
  
  
  mergedStatistics <- reactive({
    "
    reactive function that reads a stack of spreadsheet and returns a data frame 
    with descriptive statistics including circular mean, circular standard deviation 
    and nearest neighbours for the merged stack of data
    "
    
    results_df <- mergedStack()
    
    #inFileStackData <- input$stackData

    #if (!is.null(inFileStackData))
    #   results_df <- read.csv(inFileStackData$datapath, header = input$header_correlation)

    print("Data Frame:")
    print(head(results_df))
    
    #threshold <- input$max_golgi_nuclei_distance
    #if ("distance" %in% colnames(results_df)){
    #  results_df <- subset(results_df, results_df$distance < threshold)
    #}
    
    source(file = paste0(getwd(),"/src/circular_statistics.R"), local=T)
    parameters <- fromJSON(file = "parameters/parameters.json")

    feature <- parameters[input$feature_select][[1]][1]

    threshold <- input$min_nuclei_golgi_dist
    if ("distance" %in% colnames(results_df)){
      results_df <- subset(results_df, results_df$distance > threshold)
    }

    statistics_df <- as.data.frame(matrix(nrow=5,ncol=3))
    colnames(statistics_df) <- c("entity", "value") #, "comment")

    if (parameters[input$feature_select][[1]][2] == "axial") {
      
        x_data <- unlist(results_df[feature])*180.0/pi
        statistics <- compute_circular_statistics(results_df, feature, parameters)
        print("Statistics")
        print(statistics)

        p_value <- signif(statistics[1,"rayleigh_test"], digits = 3)
        #if (statistics[1,"rayleigh_test"] < 0.001)
        #    p_value <- "p < 0.001"
        p_value_mu <- signif(statistics[1,"rayleigh_test_mu"], digits = 3)
        #if (statistics[1,"rayleigh_test_mu"] < 0.001)
        #    p_value_mu <- "p < 0.001"
     

        ind <-1
        statistics_df[ind,1] <- "number of cells"
        statistics_df[ind,2] <- nrow(results_df)
        
        ind <- ind + 1
        statistics_df[ind,1] <- "mean (degree)"
        statistics_df[ind,2] <- signif(statistics[1,"mean"], digits = 3)
        
        ind <- ind + 1
        statistics_df[ind,1] <- "polarity index"
        statistics_df[ind,2] <- signif(statistics[1,"polarity_index"], digits = 3)
        
        ind <- ind + 1
        statistics_df[ind,1] <- "angular standard deviation"
        statistics_df[ind,2] <- signif(statistics[1,"std_angular"])
        statistics_df[ind,3] <- "angular standard deviation, takes values in [0,sqrt(2)], see https://doi.org/10.18637/jss.v031.i10 for more info."

        ind <- ind + 1
        statistics_df[ind,1] <- "circular standard deviation"
        statistics_df[ind,2] <- signif(statistics[1,"std_circular"], digits = 3)
        statistics_df[ind,3] <- "circular standard deviation, takes values in [0,inf], see https://doi.org/10.18637/jss.v031.i10 for more info."

        ind <- ind + 1
        statistics_df[ind,1] <- "95% confidence interval of the mean, lower limit: "
        statistics_df[ind,2] <- signif(statistics[1,"ci_lower_limit"], digits = 3)

        ind <- ind + 1
        statistics_df[ind,1] <- "95% confidence interval of the mean, upper limit: "
        statistics_df[ind,2] <- signif(statistics[1,"ci_upper_limit"], digits = 3)

        ind <- ind + 1
        statistics_df[ind,1] <- "Rayleigh test, p-value:"
        statistics_df[ind,2] <- p_value
        
        ind <- ind + 1
        statistics_df[ind,1] <- "Rayleigh test, p-value (cond. mean = 180): "
        statistics_df[ind,2] <- p_value_mu





        #entity <- c("cells", "circular sample mean (degree)", "polarity index", "Rayleigh test (p-value)", "Rayleigh test with mu=180 (p-value)")
        #values <- c( nrow(results_df), statistics[1,"mean"],  statistics[1,"polarity_index"], statistics[1,"rayleigh_test"], 0.0) # statistics[1,"rayleigh_test_mu"])
        
        #print("entity")
        #print(entity)
        #print("value")
        #print(values)
        
        #statistics_df <- data.frame(entity,values)

    }
    else if (parameters[input$feature_select][[1]][2] == "2-axial") {
        statistics <- compute_2_axial_statistics(results_df, feature, parameters)

        p_value <- signif(statistics[1,"rayleigh_test"], digits = 3)
        
        statistics_df[1,1] <- "cells"
        statistics_df[1,2] <- nrow(results_df)
        statistics_df[2,1] <- "mean (degree)"
        statistics_df[2,2] <- signif(statistics[1,"mean"], digits = 3)
        statistics_df[3,1] <- "polarity index"
        statistics_df[3,2] <- signif(statistics[1,"polarity_index"], digits = 3)
        statistics_df[4,1] <- "Rayleigh test, p-value:"
        statistics_df[4,2] <- p_value

    } else {
      
        statistics <- compute_linear_statistics(results_df, feature, parameters)
    
        statistics_df[1,1] <- "cells"
        statistics_df[1,2] <- nrow(results_df)
        statistics_df[2,1] <- "mean"
        statistics_df[2,2] <- signif(statistics[1,"mean"], digits = 3)
        statistics_df[3,1] <- "standard deviation"
        statistics_df[3,2] <- signif(statistics[1,"std"], digits = 3)
        statistics_df[4,1] <- "median"
        statistics_df[4,2] <- signif(statistics[1,"median"], digits = 3)

    }
    
    #values <- compute_polarity_index(results_df)
    #print(values)
    #polarity_index <- values[["polarity_index"]]
    ##signed_polarity_index <- values[["signed_polarity_index"]]
    #angle_mean_deg <- values[["angle_mean_deg"]]
    
    #angle_degree <- conversion.circular(results_df$angle_deg, units = "degrees", zero = 0, modulo = "2pi")
    
    #variance_degree  <- var(angle_degree)
    #mean_degree <- mean.circular(angle_degree)
    #rayleigh_test_res <- r.test(results_df$angle_deg, degree = TRUE)
    #rayleigh_test_mu_res <- v0.test(results_df$angle_deg, mu0 = 180.0, degree = TRUE)
    #rayleigh_test <- rayleigh_test_res$p.value
    #rayleigh_test_mu <- rayleigh_test_mu_res$p.value
    ###print(struct(rayleight_test))
    ##print(struct(rayleight_test_mu))
    #sd_degree  <- sd(angle_degree)
    #median_degree  <- median(angle_degree)
    
    ##entity <- c("nucleus-golgi pairs", "circular sample mean (degree)",  "circular standard deviation (degree)", "circular median (degree)", "polarity index")
    ##value <- c(nrow(results_df), angle_mean_deg , sd_degree , median_degree, polarity_index)
    
    #print("Values in data frame")
    #print(angle_mean_deg)
    #print(polarity_index)
    #print(rayleigh_test)
    #print(rayleigh_test_mu)
    
    #angle_mean_deg <- statistics
    #polarity_index <- 0
    #rayleigh_test <- 0
    #rayleigh_test_mu <- 0
    
    #entity <- c("cells", "circular sample mean (degree)", "polarity index", "Rayleigh test (p-value)", "Rayleigh test with mu=180 (p-value)")
    
    #value <- c(nrow(results_df), angle_mean_deg ,  polarity_index, rayleigh_test, rayleigh_test_mu)
    
    
    #entity <- c("cells", "circular sample mean (degree)", "polarity index", "Rayleigh test (p-value)", "Rayleigh test with mu=180 (p-value)")
    #value <- c(nrow(results_df), angle_mean_deg ,  polarity_index, rayleigh_test, rayleigh_test_mu)

    #statistics_df <- data.frame(entity,values)
    
    statistics_df
    
  })
  
  output$merged_statistics <- renderTable({
    "
    function that shows the descriptive statistics of the merged data stack in table format
    "
    
    statistics_df <- mergedStatistics()
    statistics_df 
  }, digits = 3)
  
  merged_plot <- reactive({
    
    source(file = paste0(getwd(),"/src/plot_functions.R"), local=T)
    source(file = paste0(getwd(),"/src/circular_statistics.R"), local=T)
    
    parameters <- fromJSON(file = "parameters/parameters.json")
    text_size <- as.integer(parameters["text_size_merged_plot"])
    
    results_all_df <- mergedStack()
 
    #inFileStackData <- input$stackData

    #if (!is.null(inFileStackData))
    #    results_all_df <- read.csv(inFileStackData$datapath, header = input$header_correlation)

   
    for(i in 1:nrow(results_all_df)) {
      row <- results_all_df[i,]
      a <- row$major_axis_length
      b <- row$minor_axis_length
      
      eccentricity <- sqrt(1.0 - b*b/(a*a))
      results_all_df[i,"eccentricity"] = eccentricity 
    }
    
    threshold <- input$min_nuclei_golgi_dist
    if ("distance" %in% colnames(results_all_df)){
      results_all_df <- subset(results_all_df, results_all_df$distance > threshold)
    }
    
    bin_size = 360/input$bins
    exp_condition <- input$exp_condition
    datapath <- stack_data_info$datapath 

    feature <- parameters[input$feature_select][[1]][1]
    
    print("Feature:")
    print(feature)

    if (parameters[input$feature_select][[1]][2] == "axial") {
      
        print("Axial feature!")
      
      
        x_data <- unlist(results_all_df[feature])*180.0/pi
        statistics <- compute_circular_statistics(results_all_df, feature, parameters)
        plot_title <- parameters[input$feature_select][[1]][3]
        p <- rose_plot_circular(parameters, input, statistics, x_data, plot_title, text_size)
      
    }
    else if (parameters[input$feature_select][[1]][2] == "2-axial") {
      
      x_data <- results_all_df[feature]        
      statistics <- compute_2_axial_statistics(results_all_df, feature, parameters)
      #if (input$left_axial) {
      #  x_data <- unlist(transform_2_axial(input,x_data))*180.0/pi
      #} else {
      #  x_data <- unlist(results_all_df[feature])*180.0/pi
      #}
        x_data <- unlist(transform_2_axial(input,x_data))*180.0/pi

        plot_title <- parameters[input$feature_select][[1]][3]
        p <- rose_plot_2_axial(parameters, input, statistics, x_data, plot_title, text_size)
      
    } else {
      
        x_data <- unlist(results_all_df[feature])
        statistics <- compute_linear_statistics(results_all_df, feature, parameters)
        plot_title <- parameters[input$feature_select][[1]][3]
        p <- linear_histogram(parameters, input, statistics, x_data, plot_title, text_size)
    }
    
    p

  })  
  
  output$merged_plot <- renderPlot({
    
    p <-merged_plot()
    p
  })  
  
  multi_plot <- reactive({
    
    source(file = paste0(getwd(),"/src/plot_functions.R"), local=T)
    source(file = paste0(getwd(),"/src/circular_statistics.R"), local=T)
    
    parameters <- fromJSON(file = "parameters/parameters.json")
    text_size <- 12
    
    datapath <- stack_data_info$datapath
    print(datapath)
    
    file_list <- list.files(datapath)
    print(file_list)
    
    i <- 1
    angle_dists <- list()
    file_names <- list()
    polarity_indices <- list()
    angle_mean_degs <- list()
    
    results_all_df <- mergedStack()
    
    for(row_nr in 1:nrow(results_all_df)) {
      row <- results_all_df[row_nr,]
      a <- row$major_axis_length
      b <- row$minor_axis_length
      
      eccentricity <- sqrt(1.0 - b*b/(a*a))
      results_all_df[row_nr,"eccentricity"] = eccentricity 
    }
    
    threshold <- input$min_eccentricity
    if ("eccentricity" %in% colnames(results_all_df)){
      results_all_df <- subset(results_all_df, results_all_df$eccentricity> threshold)
    }
    threshold <- input$min_nuclei_golgi_dist
    if ("distance" %in% colnames(results_all_df)){
      results_all_df <- subset(results_all_df, results_all_df$distance > threshold)
    }
    
    feature <- parameters[input$feature_select][[1]][1]
    
    plist <- vector('list', length(unique(results_all_df$filename)))
    
    for(file_name in unique(results_all_df$filename)) {
      results_df <- subset(results_all_df, results_all_df$filename == file_name )
      
      #values <- compute_polarity_index(results_df)
      
      
      #print(values)
      #polarity_index <- values[["polarity_index"]]
      #angle_mean_deg <- values[["angle_mean_deg"]]
      
      x <- unlist(results_df[feature])
      angle_dists[[i]] <- x

      file_names[[i]] <- file_name
         
      #polarity_indices[[i]] <- polarity_index
      #angle_mean_degs[[i]] <- angle_mean_deg
      i <- i+1
      
    }
    
    
    
    n <- length(angle_dists)
    nCol <- floor(sqrt(n))
    
    bin_size = 360/input$bins
    
    plotseries <- function(i){
        
        angle_dist <- angle_dists[[i]]
        file_name <- file_names[[i]]
        #polarity_index <- polarity_indices[[i]]
        #angle_mean_deg <- angle_mean_degs[[i]]
      
        results_df <- subset(results_all_df, results_all_df$filename == file_name)
      
        plot_title <- file_name
        if (nchar(file_name) > 15) {
            plot_title <- paste0("image #",toString(i))
            print(paste0("filename: ",file_name," too long, will be replaced by",plot_title))
        }        


      if (parameters[input$feature_select][[1]][2] == "axial") {
        
        statistics <- compute_circular_statistics(results_df, feature, parameters)
        #statistics <- compute_polarity_index(unlist(results_df[feature]))
        x_data <- unlist(results_df[feature])*180.0/pi
        print(paste0("Length of filename", toString(i)))
        
                p <- rose_plot_circular(parameters, input, statistics, x_data, plot_title, text_size)
        
      }
      else if (parameters[input$feature_select][[1]][2] == "2-axial") {
        
        x_data <- results_df[feature]        
        #print(x_data)
        statistics <- compute_2_axial_statistics(results_df, feature, parameters)
        #if (input$left_axial) {
        x_data <- unlist(transform_2_axial(input,x_data))*180.0/pi
        #} else {
        #  x_data <- unlist(results_df[feature])*180.0/pi
        #}
        #plot_title <- file_name
        p <- rose_plot_2_axial(parameters, input, statistics, x_data, plot_title, text_size)
        
      } else {
        
        x_data <- unlist(results_df[feature])
        statistics <- compute_linear_statistics(results_df, feature, parameters)
        #plot_title <- file_name
        p <- linear_histogram(parameters, input, statistics, x_data, plot_title)
      }
      
      

    }
    
    
    myplots <- lapply(1:length(angle_dists), plotseries)
    
    #print(myplots)
    grid.arrange(grobs = myplots, nrow = nCol) #, widths = list(10,10))
    
  })  
  
  output$multi_dist_plot <- renderPlot({
    
      multi_plot()
    
  })  
  
    output$downloadProcessedData <- downloadHandler(
        filename = function() {
            filename <- "merged_file.csv"
            return(filename)
        },
        content = function(file) {
            return(write.csv(mergedStack(), file, row.names = FALSE))
        }
    )
  
    output$downloadData <- downloadHandler(
   
        filename = function() {
            filename <- "merged_file.csv"
            if (input$dataset == "statistics_file"){
                filename <- "statistics_file.csv"
                print("Download merged_file.csv")
            }
            if (input$dataset == "merged_plot_file"){
                filename <- paste0("merge_plot", input$image_file_format)
            }
            if (input$dataset == "multi_plot_file"){
                filename <- paste0("multi_plot", input$image_file_format)
            }
            return(filename)
        },
        content = function(file) {
           
            parameters <- fromJSON(file = "parameters/parameters.json")
            width_ <- as.double(parameters["pdf_figure_size_inches"])

            if (input$dataset == "statistics_file"){
                return(write.csv(mergedStatistics(), file, row.names = FALSE))
            }
            else if ((input$dataset == "multi_plot_file") && (input$image_file_format == ".pdf")) {
                #pdf(file, width=14, height=14)
                pdf(file, family = "ArialMT", width = width_, height =width_, pointsize = 18)
                p <- multi_plot()
                plot(p)
                dev.off()
            }
            else if ((input$dataset == "multi_plot_file") && (input$image_file_format == ".png")) {
                png(file, width=960, height=960)
                p <- multi_plot()
                plot(p)
                dev.off()
            }
            else if ((input$dataset == "multi_plot_file") && (input$image_file_format == ".eps")) {
                eps(file, width=14, height=14)
                p <- multi_plot()
                plot(p)
                dev.off()
            }
            else if ((input$dataset == "merged_plot_file") && (input$image_file_format == ".pdf")){
                print("Saving merge pdf")
                pdf(file, family = "ArialMT", width = width_, height =width_, pointsize = 18)
                p <- merged_plot()
                plot(p)
                dev.off()
            }
            else if ((input$dataset == "merged_plot_file") && (input$image_file_format == ".png")){
                png(file, width=960,height=960)
                p <- merged_plot()
                plot(p)
                dev.off()
            }
            else if ((input$dataset == "merged_plot_file") && (input$image_file_format == ".eps")){
                eps(file, width=width_,height=width_)
                p <- merged_plot()
                plot(p)
                dev.off()
            }
            else (
                return(write.csv(mergedStack(), file, row.names = FALSE))
            )
        }
    )
    
    output$downloadDataSingleImage <- downloadHandler(
      filename = function() {
        if (input$datasetSingleImage == "results_file"){
          filename <- "results_file.csv"
        }
        if (input$datasetSingleImage == "statistics_file"){
          filename <- "statistics_file.csv"
        }
        if (input$datasetSingleImage == "rose_histogram"){
          filename <- "rose_histogram.pdf"
        }
        if (input$datasetSingleImage == "orientation_plot"){
          filename <- "vector_plot.pdf"
        }
        filename
      },
      content = function(file) {
        if (input$datasetSingleImage == "results_file"){
          write.csv(resultSingleImage(), file, row.names = FALSE)
        }
        if (input$datasetSingleImage == "statistics_file"){
          write.csv(singleImageStatistics(), file, row.names = FALSE)
        }
        if (input$datasetSingleImage == "rose_histogram"){
          pdf(file, width=7,height=7)
          p <- rose_histogram_single_image()
          plot(p)
          dev.off()
        }
        if (input$datasetSingleImage == "orientation_plot"){
          pdf(file, width=7,height=7)
          p <- vectorPlot()
          plot(p)
          dev.off()
        }
      }
    )
    
    
    ### Panel B
    
    plot_correlation <- reactive({
      
      parameters <- fromJSON(file = "parameters/parameters.json")
      
      inFileCorrelationData <- input$correlationData
      
      if (is.null(inFileCorrelationData))
        return(NULL)
      
      print(inFileCorrelationData$datapath)
      correlation_data <- read.csv(inFileCorrelationData$datapath, header = input$header_correlation)
      
      feature_1 <- parameters[input$feature_select_1][[1]][1]
      feature_2 <- parameters[input$feature_select_2][[1]][1]
      
      plot_df <- as.data.frame(c(correlation_data[feature_1], correlation_data[feature_2]))
      colnames(plot_df) <- c("x","y")
      p <-ggplot(plot_df, aes(x=x, y=y)) + geom_point()
      p
      
    })
    
    output$correlation_plot <- renderPlot({
      
      p <- plot_correlation()
      p
    })  
    
    
    ### Panel C
    
    output$comparison_plot <- renderPlot({
      
        source(file = paste0(getwd(),"/src/plot_functions.R"), local=T)
        source(file = paste0(getwd(),"/src/circular_statistics.R"), local=T)
    
        parameters <- fromJSON(file = "parameters/parameters.json")
        text_size <- as.integer(parameters["text_size_merged_plot"])
        
        inFileCondition_1 <- input$condition_1
        inFileCondition_2 <- input$condition_2
      
        if (is.null(inFileCondition_1))
            return(NULL)
        if (is.null(inFileCondition_2))
            return(NULL)
      
        print(inFileCondition_1$datapath)
        print(inFileCondition_2$datapath)
        cond1_data <- read.csv(inFileCondition_1$datapath, header = input$header_cond1)
        cond2_data <- read.csv(inFileCondition_2$datapath, header = input$header_cond2)
      
        feature <- parameters[input$feature_select][[1]][1]

        bin_size <- 360.0/input$bins_comparison
      #bin_size <- 20.0
      
      p <- ggplot() +
        geom_histogram(aes(cond1_data$angle_deg),
                       breaks = seq(0, 360, bin_size),
                       colour = "black",
                       fill = "green", alpha = 0.2) +
        ggtitle("cellular orientation") +
        theme(axis.text.x = element_text(size = 18)) +
        coord_polar(start = -pi/2.0, direction = -1) +
        scale_x_continuous(limits = c(0, 360),
                           breaks = (c(0, 90, 180, 270))) +
        theme_minimal(base_size = 14) +
        xlab(paste0("n = ", length(cond1_data$angle_deg))) +
        ylab("") +
        theme(axis.text.y=element_blank()) +
        scale_y_sqrt()
      
      p <- p + geom_histogram(aes(cond2_data$angle_deg),
                              breaks = seq(0, 360, bin_size),
                              colour = "black",
                              fill = "red", alpha = 0.2)
      
      
      p1 <- ggplot() +
        geom_histogram(aes(cond1_data$angle_deg),
                       breaks = seq(0, 360, bin_size),
                       colour = "white",
                       fill = "blue",
                       alpha = 0.5) +
        geom_histogram(aes(cond2_data$angle_deg),
                       breaks = seq(0, 360, bin_size),
                       colour = "white",
                       fill = "red",
                       alpha = 0.5) +
        ggtitle(paste0(" ")) +
        theme(axis.text.x = element_text(size = 18)) +
        coord_polar(start = -pi/2.0, direction = -1) +
        scale_x_continuous(limits = c(0, 360),
                           breaks = (c(0, 90, 180, 270))) +
        theme_minimal(base_size = 14) +
        xlab("") +
        ylab("") +
        theme(axis.text.y=element_blank())+
        scale_y_sqrt()
      
      #print(wilcox.test(cond1_data$angle_rad, cond2_data$angle_rad, paired=FALSE)$p.value)

        p2 <- ggplot()
        feature <- parameters[input$feature_comparison][[1]][1]
        
        if (parameters[input$feature_comparison][[1]][2] == "axial") {
      
            cond1_feature <- unlist(cond1_data[feature])*180.0/pi
            cond2_feature <- unlist(cond2_data[feature])*180.0/pi
            x_data <- list(cond1_feature, cond2_feature)            
            condition_data <- list()
            condition_data[[1]] <- cond1_data
            condition_data[[2]] <- cond2_data
            #condition_data <- list(cond1_data, cond2_data)            

            if (input$split_view_comparison) 
            { 
                
                plotseries <- function(i){
                   
                    #print(x_data)
                    #print(x_data[[i]]) 
                    
                    statistics <- compute_circular_statistics(cond1_data, feature, parameters)
                    plot_title <- parameters[input$feature_select][[1]][3]
                    p <- rose_plot_circular(parameters, input, statistics, x_data[[i]], plot_title, text_size)
                    
                    #cond_data <- condition_data[[i]]
                    #x_data <- cond_data[feature]
                    #statistics <- compute_circular_statistics(cond_data, feature, parameters)
                    #plot_title <- parameters[input$feature_select][[1]][3]
                    #p <- rose_plot_circular(parameters, input, statistics, x_data, plot_title, text_size)

                }
                myplots <- lapply(1:2, plotseries)
                p2 <- grid.arrange(grobs = myplots, ncol = 2) #, widths = list(10,10))
            }
            else {
                statistics <- compute_circular_statistics(cond1_data, feature, parameters)
                plot_title <- parameters[input$feature_select][[1]][3]
                p2 <- compare_plot_circular(parameters, input, statistics, cond1_feature, cond2_feature, plot_title)
            }
        }
        else if (parameters[input$feature_comparison][[1]][2] == "2-axial") {
            
            cond1_feature <- unlist(cond1_data[feature])
            cond2_feature <- unlist(cond2_data[feature])
            
            #if (input$left_axial) {
            cond1_feature <- unlist(transform_2_axial(input, cond1_feature))*180.0/pi
            cond2_feature <- unlist(transform_2_axial(input, cond2_feature))*180.0/pi
            #} else {
            #    cond1_feature <- unlist(cond1_data[feature])*180.0/pi
            #    cond2_feature <- unlist(cond2_data[feature])*180.0/pi
            #}

            # x_data <- list(cond1_feature, cond2_feature)            
            
            plotseries <- function(i) {
                
                ## print(x_data)
                # print(x_data[[i]]) 
                # x_data <- unlist(cond1_data[feature])*180.0/pi
                # statistics <- compute_circular_statistics(cond1_data, feature, parameters)
                # plot_title <- parameters[input$feature_select][[1]][3]
                # p <- rose_plot_2_axial(parameters, input, statistics, x_data[[i]], plot_title, text_size)
              
                plot_title <- parameters[input$feature_select][[1]][3]
                
                if (i==1) {
                    statistics <- compute_2_axial_statistics(cond1_data, feature, parameters)
                    p <- rose_plot_2_axial(parameters, input, statistics, cond1_feature, plot_title, text_size)
                } else {
                    statistics <- compute_2_axial_statistics(cond2_data, feature, parameters)
                    p <- rose_plot_2_axial(parameters, input, statistics, cond2_feature, plot_title, text_size)

                }  
                # x_data <- unlist(cond1_data[feature])*180.0/pi
                # statistics <- compute_circular_statistics(cond1_data, feature, parameters)
                # plot_title <- parameters[input$feature_select][[1]][3]
                # p <- rose_plot_2_axial(parameters, input, statistics, x_data[[i]], plot_title, text_size)
            }
            myplots <- lapply(1:2, plotseries)
            p2 <- grid.arrange(grobs = myplots, ncol = 2) 

        }
        else if (parameters[input$feature_comparison][[1]][2] == "linear") {
      
            cond1_feature <- unlist(cond1_data[feature])
            cond2_feature <- unlist(cond2_data[feature])
            
            if (input$split_view_comparison) 
            {
                
                plot_title <- parameters[input$feature_select][[1]][3]
                plotseries <- function(i) {
                
                    if (i==1) {
                        statistics <- compute_linear_statistics(cond1_data, feature, parameters)
                        p <- linear_histogram(parameters, input, statistics, cond1_feature, plot_title, text_size)
                    } else {
                        statistics <- compute_linear_statistics(cond2_data, feature, parameters)
                        p <- linear_histogram(parameters, input, statistics, cond2_feature, plot_title, text_size)
                    }  
            }
            myplots <- lapply(1:2, plotseries)
            p2 <- grid.arrange(grobs = myplots, ncol = 2) 
 
            
            } else {   

                statistics <- compute_linear_statistics(cond1_data, feature, parameters)
                plot_title <- parameters[input$feature_select][[1]][3]
                p2 <- compare_plot_linear(parameters, input, statistics, cond1_feature, cond2_feature, plot_title)
            }

        }



      p2
      #print(as.circular(cond1_data$angle_rad,type=radians))
      #circ.plot(cond1_data$angle_rad, cex =2.0)
    }) 
    
    
    comparisonStatistics <- reactive({
      "
      reactive function that reads a stack of spreadsheet and returns a data frame 
      with descriptive statistics including circular mean, circular standard deviation 
      and nearest neighbours for the merged stack of data
      "
      
      inFileCondition_1 <- input$condition_1
      inFileCondition_2 <- input$condition_2
      
      if (is.null(inFileCondition_1))
        return(NULL)
      if (is.null(inFileCondition_2))
        return(NULL)
      
      print(inFileCondition_1$datapath)
      print(inFileCondition_2$datapath)
      cond1_data <- read.csv(inFileCondition_1$datapath, header = input$header_cond1)
      cond2_data <- read.csv(inFileCondition_2$datapath, header = input$header_cond2)
      
      print(watson.two(cond1_data$angle_rad, cond2_data$angle_rad, alpha=0.05, plot=TRUE))
      #watson.two(data1, data2)
      watson.two(cond1_data$angle_rad, cond2_data$angle_rad)
      #print("Structure")
      #print(str(watson.two(cond1_data$angle_rad, cond2_data$angle_rad)))
      
      #statistics_df <- data.frame()
      #statistics_df
      
    })
    
    output$CDFPlot <- renderPlot({
      "
      function that shows the descriptive statistics in table format
      "
      inFileCondition_1 <- input$condition_1
      inFileCondition_2 <- input$condition_2
      
      if (is.null(inFileCondition_1))
        return(NULL)
      if (is.null(inFileCondition_2))
        return(NULL)
      
      print(inFileCondition_1$datapath)
      print(inFileCondition_2$datapath)
      cond1_data <- read.csv(inFileCondition_1$datapath, header = input$header_cond1)
      cond2_data <- read.csv(inFileCondition_2$datapath, header = input$header_cond2)
      
      watson.two(cond1_data$angle_rad, cond2_data$angle_rad, alpha=0.05, plot=TRUE)
    })
    
    
    output$comparison_statistics <- renderText({
      "
    function that shows the descriptive statistics in table format
    "
      
      statistics_df <- comparisonStatistics()
      statistics_df 
    })
        
}

# Run the application 
shinyApp(ui = ui, server = server)
