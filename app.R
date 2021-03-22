#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# PolApp: Shiny app for plotting and comparing polarity data (beta 0.1)
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# author:  Wolfgang Giese
# e-mail address: wolfgang #dot# giese #at# mdc-berlin #dot# de


options(shiny.maxRequestSize = 30*1024^2)

library(shiny)
library(shinyFiles)
library(circular)
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
#library(ggimage) 
library(fs)


# Create a reactive object here that we can share between all the sessions.
vals <- reactiveValues(count=0)

###### UI: User interface #########

ui <- navbarPage("PolApp - a web app for visualizing cell polarity data (beta 0.1)",
                 
                ### Panel A: Image stack histogram
                 
                tabPanel("Image stack histogram",
                          sidebarLayout(
                            sidebarPanel(
                              shinyDirButton("dir", "Input directory", "Upload"),
                              verbatimTextOutput("dir", placeholder = TRUE),
                              sliderInput("bins",
                                          "Number of bins:",
                                          min = 1,
                                          max = 30,
                                          value = 10),
                              sliderInput("max_golgi_nuclei_distance",
                                          "Max. distance nuclei to golgi",
                                          min = 0,
                                          max = 50,
                                          step = 0.1,
                                          value = 30),
                              textInput("exp_condition", "Exp. condition", "condition A"),
                              selectInput("dataset", "Choose a dataset:",
                                          choices = c("merged_file","statistics_file")),
                              downloadButton("downloadData", "Download")
                            ),
                            # Show a plot of the generated distribution
                            mainPanel(
                              tabsetPanel(
                                tabPanel("Table", tableOutput("merged_stack")),
                                tabPanel("Plot", plotOutput("merged_plot", height = "860px")),
                                tabPanel("MultiPlot", plotOutput("multi_dist_plot", height = "860px")),
                                
                                tabPanel("Statistics", tableOutput("merged_statistics"))
                              )
                            ))),
                 

            
           ### Panel B1: Single image analysis

           tabPanel("Single image analysis",
                    sidebarLayout(
                      sidebarPanel(
                        fileInput("SingleImageNuclei", "Nuclei positions",
                                  accept = c( "text/csv",
                                              "text/comma-separated-values,text/plain",
                                              ".csv",".xlsx")
                        ),    
                        tags$hr(),
                        checkboxInput("header_nuclei", "File upload", TRUE),
                        fileInput("SingleImageGolgi", "Golgi positions",
                                  accept = c( "text/csv",
                                              "text/comma-separated-values,text/plain",
                                              ".csv",".xlsx")
                        ), 
                        tags$hr(),
                        checkboxInput("header_golgi", "File upload", TRUE),
                        knobInput("flow_direction",
                                  "Flow direction (degree):",
                                  angleOffset = 90,
                                  angleArc = 360,
                                  min = 0,
                                  max = 360,
                                  value = 0),
                        sliderInput("max_golgi_nuclei_distance_single_image",
                                    "Max. distance nuclei to golgi",
                                    min = 0,
                                    max = 50,
                                    step = 0.1,
                                    value = 30),
                        checkboxInput("showNucleiLabels", "Show NucleusIDs", FALSE),
                        checkboxInput("showGolgiLabels", "Show GolgiIDs", FALSE),
                        selectInput("datasetSingleImage", "Choose a dataset:",
                                    choices = c("results_file","statistics_file","orientation_plot", "rose_histogram")),
                        fileInput("SingleImage", "Microscopic image",
                                  accept = c(".png", ".tiff", ".jpeg" )
                        ), 
                        tags$hr(),
                        checkboxInput("header_image", "File upload", TRUE),
                        sliderInput("scale_image",
                                    "Scale image",
                                    min = 0,
                                    max = 5.0,
                                    step = 0.01,
                                    value = 1.0),
                        # Button
                        downloadButton("downloadDataSingleImage", "Download")
                      ),
                      mainPanel(
                      tabsetPanel(
                        tabPanel("Plot", plotOutput("merge_nuclei_golgi", height = "1000px")),
                        tabPanel("Statistics", tableOutput("singleImageStatistics"))
                      )
                      ))),
           
           
           ### Panel B2 (remove): "Histogram"
           ### statistics for single image histograms
           
           tabPanel("Histogram",
                    sidebarLayout(
                      sidebarPanel(
                        #shinyFilesButton("file", "File select", "Please select a file", multiple = FALSE),
                        fileInput("file1", "Choose CSV File",
                                  accept = c( "text/csv",
                                              "text/comma-separated-values,text/plain",
                                              ".csv", ".xlsx")
                        ),      
                        tags$hr(),
                        checkboxInput("header", "File upload", TRUE),
                        sliderInput("bins",
                                    "Number of bins:",
                                    min = 1,
                                    max = 30,
                                    value = 10),
                        sliderInput("max_golgi_nuclei_distance",
                                    "Max. distance nuclei to golgi",
                                    min = 0,
                                    max = 5,
                                    step = 0.1,
                                    value = 5)
                      ),
                      # Show a plot of the generated distribution
                      mainPanel(
                        tabsetPanel(
                          tabPanel("Plot", plotOutput("distPlot")), 
                          tabPanel("OrientationPlot", plotOutput("orientationPlot")),
                          tabPanel("Table", tableOutput("contents")),
                          tabPanel("Summary", verbatimTextOutput("summary"))
                        )
                      ))),
           
           ### Panel C: Compare 
           ### statistical comparison of two histograms
           
           tabPanel("Compare",                    
                    sidebarLayout(
                      sidebarPanel(
                        fileInput("condition_1", "Condition 1",
                         accept = c( "text/csv",
                                     "text/comma-separated-values,text/plain",
                                     ".csv")
                         ),    
                         tags$hr(),
                         checkboxInput("header_cond1", "File upload", TRUE),
                         fileInput("condition_2", "Condition 2",
                                   accept = c( "text/csv",
                                               "text/comma-separated-values,text/plain",
                                               ".csv")
                         ), 
                         tags$hr(),
                         checkboxInput("header_cond2", "File upload", TRUE)
                      ),
                      mainPanel(
                        #tabPanel("Plot", plotOutput("comparison_plot", height = "1000px")),
                        
                        tabsetPanel(
                          tabPanel("Plot", plotOutput("comparison_plot", height = "1000px")),
                          tabPanel("CDF Plot", plotOutput("CDFPlot")), 
                          tabPanel("Statistics", tableOutput("comparison_statistics"))
                        )
                      )
                    ))
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
    
    print("Data Frame:")
    print(head(results_df))
    
    threshold <- input$max_golgi_nuclei_distance
    if ("distance" %in% colnames(results_df)){
      results_df <- subset(results_df, results_df$distance < threshold)
    }
    
    source(file = paste0(getwd(),"/src/ciruclar_statistics.R"), local=T)

    values <- compute_polarity_index(results_df)
    print(values)
    polarity_index <- values[["polarity_index"]]
    angle_mean_deg <- values[["angle_mean_deg"]]
    
    angle_degree <- conversion.circular(results_df$angle_deg, units = "degrees", zero = 0, modulo = "2pi")
    
    variance_degree  <- var(angle_degree)
    mean_degree <- mean.circular(angle_degree)
    sd_degree  <- sd(angle_degree)
    median_degree  <- median(angle_degree)
    
    entity <- c("nucleus-golgi pairs", "circular sample mean (degree)",  "circular standard deviation (degree)", "circular median (degree)", "polarity index")
    value <- c(nrow(results_df), angle_mean_deg , sd_degree , median_degree, polarity_index)
    statistics_df <- data.frame(entity,value)
    
    statistics_df
    
  })
  
  output$merged_statistics <- renderTable({
    "
    function that shows the descriptive statistics of the merged data stack in table format
    "
    
    statistics_df <- mergedStatistics()
    statistics_df 
  })
  
  output$merged_plot <- renderPlot({
    
    results_all_df <- mergedStack()
    
    threshold <- input$max_golgi_nuclei_distance
    if ("distance" %in% colnames(results_all_df)){
      results_all_df <- subset(results_all_df, results_all_df$distance < threshold)
    }
    
    bin_size = 360/input$bins
    exp_condition <- input$exp_condition
    datapath <- stack_data_info$datapath 

    source(file = paste0(getwd(),"/src/ciruclar_statistics.R"), local=T)
    values <- compute_polarity_index(results_all_df)
    print(values)
    
    
    #polarity_index <- values[["polarity_index"]]
    #angle_mean_deg <- values[["angle_mean_deg"]]
    #values <- data.frame(values)
    #print(values)
    
    p <- ggplot() +
      geom_histogram(aes(x = results_all_df$angle_deg, y = ..ncount..),
                     breaks = seq(0, 360, bin_size),
                     colour = "black",
                     fill = "grey80") +
      ggtitle("cellular orientation") +
      theme(axis.text.x = element_text(size = 18)) +
      coord_polar(start = -pi/2.0, direction = -1) +
      scale_x_continuous(limits = c(0, 360),
                         breaks = (c(0, 90, 180, 270))) +
      theme_minimal(base_size = 14) +
      xlab(sprintf("number of cells = : %s \n condition: %s", length(results_all_df$angle_rad), exp_condition)) +
      ylab("") +
      #theme(axis.text.y=element_blank()) +
      scale_y_sqrt()

  p <- p + geom_segment(data = values, aes(x=angle_mean_deg, y=0, xend=angle_mean_deg, yend=polarity_index, size = 2, color="red", lineend = "butt", ))+
        #scale_linetype_manual("segment legend",values=c("segment legend"=2)) +
        #theme(legend.title=element_blank())
        theme(legend.position = "none")
  p

  })  
  
  
  
  output$multi_dist_plot <- renderPlot({
    
    datapath <- stack_data_info$datapath
    print(datapath)
    
    file_list <- list.files(datapath)
    print(file_list)
    
    
    i <- 1
    angle_dists <- list()
    file_names <- list()
    
    results_all_df <- mergedStack()
    
    plist <- vector('list', length(unique(results_all_df$filename)))
    
    for(file_name in unique(results_all_df$filename)) {
      results_df <- subset(results_all_df, results_all_df$filename == file_name )
      
      x <- results_df$angle_deg
      angle_dists[[i]] <- x
      file_names[[i]] <- file_name
      i <- i+1
      
    }
    
    
    
    n <- length(angle_dists)
    nCol <- floor(sqrt(n))
    
    bin_size = 360/input$bins
    
    plotseries <- function(i){
      angle_dist <- angle_dists[[i]]
      file_name <- file_names[[i]]
      ggplot() +
        geom_histogram(aes(angle_dist),
                       breaks = seq(0, 360, bin_size),
                       colour = "black",
                       fill = "grey80") +
        ggtitle(file_name) +
        #theme(axis.text.x = element_text(size = 18)) +
        coord_polar(start = -pi/2.0, direction = -1) +
        scale_x_continuous(limits = c(0, 360),
                           breaks = (c(0, 90, 180, 270))) +
        #theme_minimal(base_size = 14) +
        xlab(paste0("n = ", length(angle_dist))) +
        ylab("") +
        theme(axis.text.y=element_blank()) +
        scale_y_sqrt()
    }
    
    #plotseries(2)
    
    myplots <- lapply(1:length(angle_dists), plotseries)
    
    #plotOutput(outputId, height = "400px")
    
    print(myplots)
    #myplots <- lapply(colnames(data2), plot_data_column, data = data2)
    grid.arrange(grobs = myplots, nrow = nCol) #, widths = list(10,10))
    
  })  
  

  ### Panel B1
  
  resultSingleImage <- reactive({
    "
    reactive function that reads a csv file with nuclei postitions and a csv file with golgi positions
    and returns a data frame with nucleus-golgi polarity vectors and polarity angles in degrees and radians
    with respect to a given flow direction
    "

    inFileNuclei <- input$SingleImageNuclei
    inFileGolgi <- input$SingleImageGolgi
    
    
    angle <- input$flow_direction*pi/180.0
    
    flow_vec = c(cos(angle),
                 sin(angle))
    
    if (is.null(inFileNuclei))
      return(NULL)
    if (is.null(inFileGolgi))
      return(NULL)
    
    print(inFileNuclei$datapath)
    print(inFileGolgi$datapath)
    if (file_ext(inFileNuclei$datapath) == "csv"){
      print(input$header_nuclei)
      nuclei_data <- read.csv(inFileNuclei$datapath, header = input$header_nuclei)
      golgi_data <- read.csv(inFileGolgi$datapath, header = input$header_golgi)
    }
    
    if (file_ext(inFileNuclei$datapath) == "xlsx"){
      #nuclei_data <- read.xlsx(inFileNuclei$datapath, 1,header = input$header_nuclei)
      #golgi_data <- read.xlsx(inFileGolgi$datapath, 1,header = input$header_nuclei)
      nuclei_data <- data.frame(read_excel(inFileNuclei$datapath, 1))#,header = input$header_nuclei)
      golgi_data <- data.frame(read_excel(inFileGolgi$datapath, 1))#,header = input$header_nuclei)
    }
    
    print(head(nuclei_data))
    print(head(golgi_data))
    
    find_nearest_neighbour <- function(x,y,position_data) {
      min_dist2 = (position_data[1,'X'] -x)**2
      min_dist2 = min_dist2 + (position_data[1,'Y'] -y)**2	
      min_X = position_data[1,'X']
      min_Y = position_data[1,'Y']
      
      index <- 1
      
      for(k in 2:nrow(position_data)) {
        dist2 = (position_data[k,'X'] -x)**2
        dist2 = dist2 + (position_data[k,'Y'] -y)**2
        if (dist2 < min_dist2) {
          min_dist2 = dist2
          min_X = position_data[k,'X']
          min_Y = position_data[k,'Y']
          index <- k
        }
      }
      return(index)
    }
    
    angle <- function(v1,v2){
      dot.prod <- v1%*%v2
      norm.v1 <- norm(v1,type="2")
      norm.v2 <- norm(v2,type="2")
      #theta <- acos(dot.prod / (norm.v1 * norm.v2))
      theta = atan2(v2[2], v2[1]) - atan2(v1[2], v1[1]);
      if (theta < 0) { theta = theta + 2 * pi;}
      as.numeric(theta)
    }
    
    results_df <- data.frame(X_nuclei=double(),
                             Y_nuclei=double(),
                             X_golgi=double(),
                             Y_golgi=double(),
                             distance=double(),
                             vector_X=double(),
                             vector_Y=double(),
                             angle_rad=double(),
                             angle_deg=double())
    
    
    for(i in 1:nrow(nuclei_data)) {
      row <- nuclei_data[i,]
      golgi_pos <- find_nearest_neighbour(row$X,row$Y,golgi_data)
      nuclei_pos <- find_nearest_neighbour(golgi_data[golgi_pos,'X'],golgi_data[golgi_pos,'Y'],nuclei_data)
      
      if (nuclei_pos == i) {
        golgi_X <- golgi_data[golgi_pos,'X']
        golgi_Y <- golgi_data[golgi_pos,'Y']
        nucleus_X <- nuclei_data[nuclei_pos,'X']
        nucleus_Y <- nuclei_data[nuclei_pos,'Y']
        
        distance = sqrt((golgi_X - nucleus_X)**2 + (golgi_Y - nucleus_Y)**2)
        vector_X = golgi_X - nucleus_X
        vector_Y = golgi_Y - nucleus_Y
        cell_vec = c(vector_X,vector_Y)
        angle_rad = angle(cell_vec,flow_vec)
        angle_deg = 180.0*angle_rad/pi
        results_df = rbind(results_df, c(nuclei_pos,golgi_pos,nucleus_X,nucleus_Y,golgi_X,golgi_Y,distance,vector_X,vector_Y,angle_rad,angle_deg))  
      }	
    }
    
    
    colnames(results_df) <- c("nucleus_ID","golgi_ID","X_nuclei","Y_nuclei","X_golgi","Y_golgi","distance","vector_X","vector_Y","angle_rad","angle_deg")
    
    results_df
  })
  
  singleImageStatistics <- reactive({ 
    "
    reactive function that reads a spreadsheet with nuclei postitions and a spreadsheet with golgi positions
    and returns a data frame with descriptive statistics including circular mean, circular standard deviation 
    and nearest neighbours
    "
    
    #TODO: remove code duplication here
    
    inFileNuclei <- input$SingleImageNuclei
    inFileGolgi <- input$SingleImageGolgi
    
    nuclei_data <- read.csv(inFileNuclei$datapath, header = input$header_nuclei)
    golgi_data <- read.csv(inFileGolgi$datapath, header = input$header_golgi)
    
    print(head(nuclei_data))
    print(head(golgi_data))
    
    #res <- get.knn(nuclei_data %>% select(X,Y), k=1)
    res <- get.knn(cbind(nuclei_data$X,nuclei_data$Y), k=1)
    mean_nn <- mean(res$nn.dist)
    
    res <- get.knn(cbind(nuclei_data$X,nuclei_data$Y), k=4)
    mean_4nn <- mean(res$nn.dist)

    
    results_df <- resultSingleImage()
    
    threshold <- input$max_golgi_nuclei_distance_single_image
    results_df <- subset(results_df, results_df$distance < threshold)
    
    print(results_df)
    
    #TODO: check if input is in intervall [0,360]
    
    angle_degree <- conversion.circular(results_df$angle_deg, units = "degrees", zero = pi)
    
    variance_degree  <- var(angle_degree)
    mean_degree <- mean(angle_degree)
    sd_degree  <- sd(angle_degree)
    median_degree  <- median(angle_degree)

    entity <- c("number nuclei positions", "number golgi positions", "nucleus-golgi pairs","mean nearest neighbour distances nuclei", "mean 4 nearest neighbour distances nuclei",
                "circular sample mean (degree)",  "circular standard deviation (degree)", "circular median (degree)")
    value <- c(nrow(nuclei_data), nrow(golgi_data), nrow(results_df), mean_nn, mean_4nn, mean_degree , sd_degree , median_degree)
    statistics_df <- data.frame(entity,value)
    
    statistics_df
    
    })
  
  
  output$singleImageStatistics <- renderTable({
    "
    function that shows the descriptive statistics in table format
    "
    
    statistics_df <- singleImageStatistics()
    
  })
  
  vectorPlot <- reactive({
    
    threshold <- input$max_golgi_nuclei_distance_single_image
    
    
    
    results_df <- resultSingleImage()
    
    results_df <- subset(results_df, distance < threshold)

    x_max <-  max(c(max(results_df$X_golgi),results_df$X_nuclei))*1.1
    y_max <-  max(c(max(results_df$Y_golgi),results_df$Y_nuclei))*1.1
    
    scale_image <- input$scale_image
    
    p4 <- ggplot(results_df) +
      geom_segment(aes(x = X_nuclei, y = y_max - Y_nuclei, xend =  X_golgi, yend = y_max - Y_golgi), arrow=arrow(length = unit(0.1,"cm")), size=0.5, color="red")+
      ggtitle("nuclei-golgi vectors") +
      xlim(0,x_max*scale_image) +
      ylim(0,x_max*scale_image) +
      xlab("x") +
      ylab("y") +
      theme_bw()

    inFileImage <- input$SingleImage
    if (!is.null(inFileImage)) {
      #p4 <- p4 +
      #  xlim(0,500) +
      #  ylim(0,500)
      img = inFileImage$datapath 
      p4 <- ggbackground(p4, img)
      

    }
    
    if (input$showGolgiLabels) 
      p4 <- p4 +  geom_text(aes( x = X_golgi, y = y_max - Y_golgi, label = golgi_ID) ) 
    if (input$showNucleiLabels) 
      p4 <- p4 + geom_text(aes( x = X_nuclei, y = y_max - Y_nuclei, label = nucleus_ID) )
    
    p4
    
  })  
  
  
  rose_histogram_single_image <- reactive({
    
    results_df <- resultSingleImage()
    
    threshold <- input$max_golgi_nuclei_distance_single_image
    
    results_df <- subset(results_df, results_df$distance < threshold)
    
    bin_size = 360/16
    
    #threshold <- input$max_golgi_nuclei_distance_single_image
    
    #results_df <- subset(results_df, results_df$distance < threshold)
    
    p1 <- ggplot() +
      geom_histogram(aes(results_df$angle_deg),
                     breaks = seq(0, 360, bin_size),
                     colour = "black",
                     fill = "grey80") +
      ggtitle("cellular orientation w.r.t flow") +
      theme(axis.text.x = element_text(size = 18)) +
      coord_polar(start = -pi/2.0, direction = -1) +
      scale_x_continuous(limits = c(0, 360),
                         breaks = (c(0, 90, 180, 270))) +
      theme_minimal(base_size = 14) +
      xlab(paste0("n = ", length(results_df$angle_rad))) +
      ylab("") +
      theme(axis.text.y=element_blank()) +
      scale_y_sqrt()
    
    p1
  })
  

    output$distPlot <- renderPlot({
        
        inFile <- input$file1
        
        if (is.null(inFile))
            return(NULL)
        
        print(list.files(inFile$datapath))
        
        results_df <- read.csv(inFile$datapath, header = input$header)
        
        
        results_df <- subset( results_df, results_df$distance < input$max_golgi_nuclei_distance)
        
        
        x <- results_df$angle_deg

        bin_size = 360/input$bins
        
        ggplot() +
            geom_histogram(aes(results_df$angle_deg),
                           breaks = seq(0, 360, bin_size),
                           colour = "black",
                           fill = "grey80") +
            ggtitle("cellular orientation") +
            theme(axis.text.x = element_text(size = 18)) +
            coord_polar(start = -pi/2.0, direction = -1) +
            scale_x_continuous(limits = c(0, 360),
                               breaks = (c(0, 90, 180, 270))) +
            theme_minimal(base_size = 14) +
            xlab(paste0("n = ", length(results_df$angle_rad))) +
            ylab("") +
            theme(axis.text.y=element_blank()) +
            scale_y_sqrt()
    })
    
    
    output$orientationPlot <- renderPlot({
      
        inFile <- input$file1
        
        if (is.null(inFile))
            return(NULL)
        
        results_df <- read.csv(inFile$datapath, header = input$header)

        scale_factor = 1.4453
        
        plot(NA,NA, xlim=c(0,10), ylim=c(0,10))
        
        Arrows(25,  975, 125, 975, lwd=10, arr.type="triangle", arr.width = 1.0, arr.col = "red", col="red")
        
        for(i in 1:nrow(results_df)){
          row <- results_df[i,]
          Arrows(row$X_nuclei*scale_factor,  row$Y_nuclei*scale_factor, row$X_golgi*scale_factor, row$Y_golgi*scale_factor, lwd=2, arr.type="triangle",arr.col = "yellow")
          
        }
    })
    
    output$merge_nuclei_golgi <- renderPlot({
      
      angle <- input$flow_direction*pi/180.0
      
      flow_vec = c(cos(angle),
                   sin(angle))
    
      results_df <- resultSingleImage()
      
      if (is.null(results_df))
        return(NULL)
      
      #bin_size = 360/input$bins
      bin_size = 360/16
      
      threshold <- input$max_golgi_nuclei_distance_single_image
      
      results_df <- subset(results_df, results_df$distance < threshold)

      
      p1 <- rose_histogram_single_image()
      
      p2 <- ggplot() +
          geom_segment(aes(x = 0.0, y = 0.0, xend =  flow_vec[1], yend = flow_vec[2], colour = "segment"), arrow=arrow(), size=2, color="blue")+
          ggtitle("flow direction") +
          xlim(-1.0,1.0) +
          ylim(-1.0,1.0) +
          theme_bw()
      
      #img <- png::readPNG("/media/wgiese/DATA/Projects/Olya/polarity_analysis/PolApp/arch+1_golgi.png")
      
      
      

      
      x    <- results_df$distance
      #bins <- seq(min(x), max(x), length.out = input$bins + 1)
      bins <- seq(min(x), max(x), length.out = 10)
      
      # draw the histogram with the specified number of bins
      #p3 <- hist(x, breaks = bins, col = 'darkgray', border = 'white')
      
      
      
      p3 <- ggplot(data=results_df, aes(results_df$distance)) + 
            geom_histogram(binwidth=2) + 
            geom_vline(xintercept=threshold, col="red") +
            xlab("nuclei-golgi distance") +
            ylab("count") +
            theme_bw()
            
      
      p4 <- vectorPlot() 
      
      #p4 <- ggplot() +
      #  geom_segment(aes(x = 0.0, y = 0.0, xend =  flow_vec[1], yend = flow_vec[2], colour = "segment"), arrow=arrow(), size=2, color="blue")+
      #  ggtitle("flow direction") +
      #  xlim(-1.0,1.0) +
      #  ylim(-1.0,1.0)
      
      
      #p3 <- ggplot(data=nuclei_data, aes(nuclei_data$Area)) + geom_histogram()
      
      #p4 <- ggplot(data=golgi_data, aes(golgi_data$Area)) + geom_histogram()
      
      grid.arrange(p1, p2, p3, p4, nrow = 2)
      
    })
    
    #print(merged_stack_df)
    
    output$downloadData <- downloadHandler(
      filename = function() {
        filename <- "merged_file.csv"
        if (input$dataset == "statistics_file"){
          filename <- "statistics_file.csv"
        }
        return(filename)
      },
      content = function(file) {
        if (input$dataset == "statistics_file"){
          return(write.csv(mergedStatistics(), file, row.names = FALSE))
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
    
    
    
    
    
    output$summary <- renderText({ 
        
        inFile <- input$file1
        
        if (is.null(inFile))
            return(NULL)
        
        results_df <- read.csv(inFile$datapath, header = input$header)
        
        cos_sum <- 0.0
        sin_sum <- 0.0
        mean_angle <- 0.0
        
        for (row in 1:nrow(results_df)) {
            #print(row)
            cos_sum <- cos_sum + cos(results_df[row,"angle_rad"]) 
            sin_sum <- sin_sum + sin(results_df[row,"angle_rad"]) 
        }
        
        if ( nrow(results_df) > 0) {
            #cos_sum <- cos_sum/nrow(results_df)
            #sin_sum <- sin_sum/nrow(results_df)
            mean_angle <- atan2(sin_sum/nrow(results_df),cos_sum/nrow(results_df))
        }
        
        print(angular.variance(results_df$angle_rad))
        #print(mean.circular(as.circular(results_df$angle_rad))*180.0/3.14)
        
        
        
        output1 <- paste0("Mean angle of the distribution ",mean_angle*180.0/3.14 ,"\n")
        output2 <- paste0("Variance angle of the distribution ",angular.variance(results_df$angle_rad)*180.0/3.14 ,"\n")
        
        
        
        
        
        output3 <- "For details of the computation see (\"Mean_of_circular_quantities\")\n"
        output4 <- "https://en.wikipedia.org/wiki/Mean_of_circular_quantities \n https://en.wikipedia.org/wiki/Directional_statistics"
        print(paste0(output1,output2,output3,output4))
        
    })
    
    
    output$contents <- renderTable({
        # input$file1 will be NULL initially. After the user selects
        # and uploads a file, it will be a data frame with 'name',
        # 'size', 'type', and 'datapath' columns. The 'datapath'
        # column will contain the local filenames where the data can
        # be found.
        inFile <- input$file1
        
        if (is.null(inFile))
            return(NULL)
        
        read.csv(inFile$datapath, header = input$header)
    })
    
    ### Panel C
    
    output$comparison_plot <- renderPlot({
      
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
      
      bin_size = 16
      
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
      
      p1
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
