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


options(shiny.maxRequestSize = 30 * 1024^2)

library(shiny)
library(shinyFiles)
library(shinycssloaders)
library(circular)
# library(CircMLE)
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
#library(fs)
library(rjson)
library(optparse)

option_list <- list(
  make_option(c("-p", "--port"), type = "integer", default = 8888)
)
opt <- parse_args(OptionParser(option_list = option_list))

# Discussion of color palettes https://thenode.biologists.com/data-visualization-with-flying-colors/research/ and more examples of use see https://huygens.science.uva.nl/PlotTwist/
# Color palettes Paul Tol: https://personal.sron.nl/~pault/

# From Paul Tol: https://personal.sron.nl/~pault/
Tol_bright <- c("#EE6677", "#228833", "#4477AA", "#CCBB44", "#66CCEE", "#AA3377", "#BBBBBB")

Tol_muted <- c("#88CCEE", "#44AA99", "#117733", "#332288", "#DDCC77", "#999933", "#CC6677", "#882255", "#AA4499", "#DDDDDD")

Tol_light <- c("#BBCC33", "#AAAA00", "#77AADD", "#EE8866", "#EEDD88", "#FFAABB", "#99DDFF", "#44BB99", "#DDDDDD")

# From Color Universal Design (CUD): https://jfly.uni-koeln.de/color/
Okabe_Ito <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")


# Create a reactive object here that we can share between all the sessions.
vals <- reactiveValues(count = 0)

###### UI: User interface #########

ui <- navbarPage(
  "Polarity JaM - a web app for visualizing cell polarity, junction and morphology data",

  ### Panel 0: Data preparation

  tabPanel(
    "Data preparation",
    sidebarLayout(
      sidebarPanel(

        # radioButtons("data_upload_form", "Data from:", choices = list("example 1", "single file", "folder", "key file"), selected = "example 1"),
        radioButtons("data_upload_form", "Data from:", choices = list("example 1", "upload data"), selected = "example 1"),
        conditionalPanel(
          condition = "input.data_upload_form == 'upload data'",
          checkboxInput("terms_of_use", "I agree to terms of use", FALSE),
        ),
        conditionalPanel(
          condition = "input.terms_of_use == true",
          radioButtons("data_upload_source", "Data from:", choices = list("single file", "folder", "key file"), selected = "single file"),
        ),
        conditionalPanel(
          condition = "(input.data_upload_source == 'single file') & (input.terms_of_use == true)",
          fileInput("stackData", "Upload data file",
            accept = c(
              "text/csv",
              "text/comma-separated-values,text/plain",
              ".csv", ".xlsx"
            )
          ),
          tags$hr(),
          checkboxInput("header_correlation", "File upload", TRUE),
        ),
        conditionalPanel(
          condition = "(input.data_upload_source == 'folder') & (input.terms_of_use == true)",
          shinyDirButton("dir", "Input directory", "Upload"),
          verbatimTextOutput("dir", placeholder = TRUE),
          actionButton("refreshStack", "Refresh"),
        ),
        conditionalPanel(
          condition = "(input.data_upload_source == 'key file') & (input.terms_of_use == true)",
          fileInput("keyData", "Upload catalogue",
            accept = c(
              "text/csv",
              "text/comma-separated-values,text/plain",
              ".csv", ".xlsx"
            )
          ),
          tags$hr(),
          checkboxInput("header_correlation_key", "File upload", TRUE),
          shinyDirButton("dir_key", "Input directory", "Upload"),
          verbatimTextOutput("dir_key", placeholder = TRUE),
          actionButton("refreshStack", "Refresh"),
        ),

        # in case of a very large spreadsheet file, the rows can be subsampled. In this case only every n-th row is selected.
        checkboxInput("subsample_data", "Subsample data", FALSE),
        conditionalPanel(
          condition = "input.subsample_data == true",
          numericInput("subsample_n", "Select every n-th row:", value = 1, min = 1, max = 50, step = 1)
        ),

        # 
        selectInput("sample_col", "Identifier of samples", choices = ""),
        selectInput("condition_col", "Identifier of conditions", choices = ""),
        
        #TODO: needs to be implemented
        selectInput("filter_column", "Filter based on this parameter:", choices = ""),
        #TODO: needs to be implemented
        selectInput("remove_these_conditions", "Deselect these conditions:", "", multiple = TRUE),
        selectInput("dataset_merged", "Choose a dataset:",
          choices = c("merged_file")
        ),
        downloadButton("downloadProcessedData", "Download")
      ),

      # TODO: Add Terms of Use text
      mainPanel(
        tabsetPanel(
          tabPanel("Data", htmlOutput("terms_of_use_text"), tableOutput("merged_stack"))
        )
      )
    )
  ),


  ### Panel A: Image stack analysis

  tabPanel(
    "Image stack analysis",
    sidebarLayout(
      sidebarPanel(
        selectInput("feature_select", "Choose a feature:", choices = ""),
        selectInput("stats_method", "Choose a stats test",
          choices = c("None", "Rayleigh uniform", "V-Test", "Rao's Test", "Watson's Test")
        ),
        conditionalPanel(
          condition = "input.stats_method == 'V-Test'",
          numericInput("cond_mean_direction",
            "Conditional mean direction",
            value = 180
          ),
          NULL
        ),
        checkboxInput("ci_plot", "Confidence interval (CI)", TRUE),
        conditionalPanel(
          condition = "input.ci_plot == true",
          selectInput("ci_method", "CI method",
            choices = c(
              "95% CI of the mean", "90% CI of the mean", "50% CI of the mean",
              "circular standard deviation", "angular standard deviation"
            )
          )
        ),
        checkboxInput("histogram_plot", "Histogram plot", TRUE),
        conditionalPanel(
          condition = "input.histogram_plot == true",
          sliderInput("bins",
            "Number of bins:",
            min = 4,
            max = 36,
            value = 12
          ),
        ),
        checkboxInput("scatter_plot", "Scatter plot", FALSE),
        checkboxInput("kde_plot", "KDE plot", FALSE),
        checkboxInput("area_scaled", "area scaled histogram", TRUE),
        # checkboxInput("left_directional", "hemirose on left", FALSE),

        checkboxInput("filter_data", "filter data", FALSE),
        conditionalPanel(
          condition = "input.filter_data == true",
          sliderInput("min_eccentricity",
            "Mininum eccentricity",
            min = 0,
            max = 1,
            step = 0.1,
            value = 0.0
          ),
          sliderInput("min_nuclei_golgi_dist",
            "Minimum nuclei golgi distance",
            min = 0,
            max = 10,
            step = 1,
            value = 0
          ),
        ),
        selectInput("plot_mode", "Choose data modality:",
          choices = c("circular", "semicircular", "linear")
        ),
        conditionalPanel(
          condition = "input.plot_mode == 'semicircular'",
          selectInput("hemi_rose_options", "Hemirose plot options:",
            choices = c("up", "down", "left", "right", "all")
          )
        ),
        selectInput("select_colormap", "Choose a color scheme",
          choices = c("gray", "Okabe_Ito", "Tol_bright", "Tol_muted", "Tol_light")
        ),
        conditionalPanel(
          condition = "input.select_colormap != 'gray'",
          numericInput("select_color", "Select a color from color scheme:", value = 1, min = 1, max = 10, step = 1),
        ),
        checkboxInput("adjust_alpha", "adjust transparency", FALSE),
        conditionalPanel(
          condition = "input.adjust_alpha == true",
          numericInput("alpha_fill", "set alpha fill:", value = 0.5, min = 0.0, max = 1.0, step = 0.1),
          selectInput("outline", "choose outline style:", choice = c("color", "white", "black"))
        ),
        numericInput("plot_height_A", "Height (# pixels): ", value = 720),
        numericInput("plot_width_A", "Width (# pixels):", value = 1280),
        selectInput("dataset", "Choose a dataset:",
          choices = c("statistics_file", "merged_plot_file", "multi_plot_file")
        ),
        # selectInput("image_file_format", "Choose image file format:",
        #            choices = c(".pdf",".eps",".png")),
        downloadButton("downloadData", "Download")
      ),

      # Show a plot of the generated distribution
      mainPanel(
        tabsetPanel(
          #                    tabPanel("Table", tableOutput("merged_stack")),
          tabPanel(
            "Plot", downloadButton("downloadMultiPlotPDF", "Download pdf-file"),
            downloadButton("downloadMultiPlotEPS", "Download eps-file"),
            downloadButton("downloadMultiPlotSVG", "Download svg-file"),
            downloadButton("downloadMultiPlotPNG", "Download png-file"),
            div(`data-spy` = "affix", `data-offset-top` = "10", withSpinner(plotOutput("multi_dist_plot", height = "120%"))),
            # textOutput("parameter_error"),
            NULL,
          ),
          tabPanel(
            "MergedPlot", downloadButton("downloadMergedPlotPDF", "Download pdf-file"),
            downloadButton("downloadMergedPlotEPS", "Download eps-file"),
            downloadButton("downloadMergedPlotSVG", "Download svg-file"),
            downloadButton("downloadMergedPlotPNG", "Download png-file"),
            div(`data-spy` = "affix", `data-offset-top` = "10", withSpinner(plotOutput("merged_plot", height = "120%"))),
            textOutput("parameter_error"),
            NULL,
          ),
          tabPanel("Statistics", tableOutput("merged_statistics"))
        )
      )
    )
  ),



  ### Panel B: Correlation analysis

  tabPanel(
    "Correlation analysis",
    sidebarLayout(
      sidebarPanel(
        #                fileInput("correlationData", "Upload data file",
        #                            accept = c( "text/csv",
        #                            "text/comma-separated-values,text/plain",
        #                            ".csv",".xlsx")),
        #                            tags$hr(),
        #                checkboxInput("header_correlation", "File upload", TRUE),
        selectInput("feature_select_1", "Choose a feature:", choices = ""),
        selectInput("feature_select_2", "Choose a feature:", choices = ""),
        # selectInput("feature_select_1", "Choose a feature 1:",
        #            choices = c("organelle_orientation","major_axis_shape_orientation","major_axis_nucleus_orientation","eccentricity","mean_expression","area","perimeter")),
        # selectInput("feature_select_2", "Choose a feature 2:",
        #            choices = c("organelle_orientation","major_axis_shape_orientation","major_axis_nucleus_orientation","eccentricity","mean_expression","area","perimeter")),
        selectInput("datasetSingleImage", "Download:",
          choices = c("results_file", "statistics_file", "orientation_plot", "rose_histogram")
        ),
        # tags$hr(),
        selectInput("corr_plot_option", "Choose a plot option:",
          choices = c("correlation plot", "spoke plot")
        ),
        conditionalPanel(
          condition = "input.corr_plot_option == 'correlation plot'",
          checkboxInput("center_corr_plot", "center correlation plot", TRUE),
        ),
        conditionalPanel(
          condition = "input.corr_plot_option == 'spoke plot'",
          numericInput("spoke_subsample_n", "Subsample every n-th row:", value = 1, min = 1, max = 50, step = 1)
        ),
        numericInput("text_size_corr", "text size", value = 24, min = 4, max = 50, step = 1),
        numericInput("marker_size_corr", "marker size", value = 3, min = 1, max = 20, step = 1),
        numericInput("plot_height_corr", "Height (# pixels): ", value = 600),
        numericInput("plot_width_corr", "Width (# pixels):", value = 800),
        checkboxInput("header_image", "File upload", TRUE),
        downloadButton("downloadCorrelationData", "Download")
      ),
      mainPanel(
        tabsetPanel(
          tabPanel(
            "Plot", downloadButton("downloadPlotPDF", "Download pdf-file"),
            downloadButton("downloadPlotEPS", "Download eps-file"),
            downloadButton("downloadPlotSVG", "Download svg-file"),
            downloadButton("downloadPlotPNG", "Download png-file"),
            div(`data-spy` = "affix", `data-offset-top` = "10", withSpinner(plotOutput("correlation_plot", height = "120%"))),
            NULL,
          ),
          tabPanel("Statistics", tableOutput("correlation_statistics"))
          # plotOutput("correlation_plot", height = "1000px")),#,
          # tabPanel("Spoke Plot", plotOutput("spoke_plot", height = "1000px"))#,
          # tabPanel("Statistics", tableOutput("singleImageStatistics"))
        )
      )
    )
  ),

  ### Panel C: Comparison statistics

  tabPanel(
    "Compare",
    sidebarLayout(
      sidebarPanel(
        #                fileInput("control_condition", "Control condition",
        #                            accept = c( "text/csv",
        #                            "text/comma-separated-values,text/plain",
        #                            ".csv")),
        #                tags$hr(),
        #                checkboxInput("header_cond1", "File upload", TRUE),

        #                fileInput("condition_2", "Condition 2",
        #                            accept = c( "text/csv",
        #                            "text/comma-separated-values,text/plain",
        #                            ".csv")),
        #                tags$hr(),
        #                checkboxInput("header_cond2", "File upload", TRUE),
        #                sliderInput("bins_comparison",
        #                            "Number of bins:",
        #                            min = 1,
        #                            max = 30,
        #                            value = 12),


        selectInput("control_condition", "control condition", choices = ""),
        selectInput("feature_comparison", "Choose a feature:",
          choices = c(
            "organelle_orientation", "major_axis_shape_orientation",
            "major_axis_nucleus_orientation", "eccentricity", "major_over_minor_ratio",
            "mean_expression", "marker_polarity", "area", "perimeter"
          )
        ),
        checkboxInput("kde_comparison", "KDE plot", FALSE),
        checkboxInput("histogram_comparison", "Histogram plot", TRUE),
        #                checkboxInput("split_view_comparison", "Split view", TRUE),
      ),
      mainPanel(
        # tabPanel("Plot", plotOutput("comparison_plot", height = "1000px")),
        tabsetPanel(
          tabPanel("Plot", plotOutput("comparison_plot", height = "1000px")),
          tabPanel("CDF Plot", plotOutput("CDFPlot")),
          tabPanel("Statistics", tableOutput("comparison_statistics"))
        )
      )
    )
  ),

  ### Panel D: Terms of Use

  tabPanel(
    "Terms of Use",
    sidebarLayout(
      sidebarPanel(
        checkboxInput("terms_of_use_all", "I agree to the terms of use", FALSE),
      ),
      mainPanel(
        tabsetPanel(
          tabPanel("Text", htmlOutput("terms_of_use_text_all"))
        )
      )
    )
  )
)


loadVolumes <- function(exclude = "") {
  "This function returns all volumes in a named list.
  TODO: check if built-in function would also work."

  osSystem <- Sys.info()["sysname"]
  if (osSystem == "Darwin") {
    volumes <- dir_ls("/Volumes")
    names(volumes) <- basename(volumes)
  } else if (osSystem == "Linux") {
    volumes <- c(Computer = "/")
    media <- c(media = "/media/")
    home <- c(home = "~")
    # media <- dir_ls("/media")
    # names(media) <- basename(media)
    volumes <- c(volumes, media, home)
  } else if (osSystem == "Windows") {
    wmic <- paste0(Sys.getenv("SystemRoot"), "\\System32\\Wbem\\WMIC.exe")
    if (!file.exists(wmic)) {
      message("\nThe wmic program does not seem to be in the default location")
      message("Please report this problem and include output from the command")
      message("'where wmic' to https://github.com/thomasp85/shinyFiles/issues")
      volumes <- Sys.getenv("HOMEDRIVE")
      volNames <- ""
    } else {
      volumes <- system(paste(wmic, "logicaldisk get Caption"),
        intern = T
      )
      volumes <- sub(" *\\r$", "", volumes)
      keep <- !tolower(volumes) %in% c(
        "caption",
        ""
      )
      volumes <- volumes[keep]
      volNames <- system(paste(wmic, "logicaldisk get VolumeName"),
        intern = T
      )
      volNames <- sub(" *\\r$", "", volNames)
      volNames <- volNames[keep]
      volNames <- paste0(volNames, ifelse(volNames == "",
        "", " "
      ))
    }
    volNames <- paste0(volNames, "(", volumes, ")")
    names(volumes) <- volNames
    volumes <- gsub(":$", ":/", volumes)
  } else {
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
    "dir",
    roots = volumes # ,
  )
  shinyDirChoose(
    input,
    "dir_key",
    roots = volumes # ,
  )

  dir <- reactive(input$dir)
  stack_data_info <- reactiveValues(datapath = getwd())

  output$dir <- renderText({
    print(stack_data_info$datapath)
    stack_data_info$datapath
  })

  observeEvent(
    ignoreNULL = TRUE,
    eventExpr = {
      input$dir
    },
    handlerExpr = {
      if (!"path" %in% names(dir())) {
        return()
      }

      rootdir <- normalizePath(paste(volumes[input$dir[[2]]][[1]]))
      stack_data_info$datapath <- file.path(rootdir, paste(unlist(dir()$path[-1]), collapse = .Platform$file.sep))

      # home <- normalizePath("~")
      # stack_data_info$datapath <-
      #  file.path(paste(unlist(dir()$path[-1]), collapse = .Platform$file.sep))
    }
  )
  dir_key <- reactive(input$dir_key)
  key_data_info <- reactiveValues(datapath = getwd())

  output$dir_key <- renderText({
    print(key_data_info$datapath)
    key_data_info$datapath
  })

  observeEvent(
    ignoreNULL = TRUE,
    eventExpr = {
      input$dir_key
    },
    handlerExpr = {
      if (!"path" %in% names(dir_key())) {
        return()
      }

      rootdir_key <- normalizePath(paste(volumes[input$dir_key[[2]]][[1]]))
      key_data_info$datapath <- file.path(rootdir_key, paste(unlist(dir_key()$path[-1]), collapse = .Platform$file.sep))

      # home <- normalizePath("~")
      # stack_data_info$datapath <-
      #  file.path(paste(unlist(dir()$path[-1]), collapse = .Platform$file.sep))
    }
  )

  observe({
    var_names <- colnames(mergedStack())
    print("var_names")
    print(var_names)
    if (length(var_names) > 0) {
      var_list <- c("none", var_names)
    } else {
      print("var_name is not a list")
      var_list <- c("none")
    }
    print("var_list")
    print(var_list)
    # }
 
    updateSelectInput(session, "sample_col", choices = var_list, selected = "label")
    updateSelectInput(session, "condition_col", choices = var_list, selected = "filename")
    updateSelectInput(session, "feature_select", choices = var_list, selected = "cell_shape_orientation")
    updateSelectInput(session, "feature_select_1", choices = var_list, selected = "cell_shape_orientation")
    updateSelectInput(session, "feature_select_2", choices = var_list, selected = "nuc_shape_orientation")
    updateSelectInput(session, "feature_comparison", choices = var_list, selected = "nuclei_golgi_polarity")
    updateSelectInput(session, "filter_column", choices = var_list, selected="none")
    updateSelectInput(session, "filter_column", choices = var_list, selected="none")
  })


  mergedStack <- reactive({
    "
    reactive function that reads all csv files 
    from the directory given in stack_data_info$datapath 
    and combines them into one data frame
    "

    inFileStackData <- input$stackData

    if (input$data_upload_form == "example 1") {
      
      results_all_df <- read.csv("example_1/example_1.csv", header = TRUE)
    
    } else if ((input$data_upload_source == "single file") & !is.null(inFileStackData)) {
    
      results_all_df <- read.csv(inFileStackData$datapath, header = input$header_correlation)
    
    } else if (input$data_upload_source == "folder") {
      datapath <- stack_data_info$datapath

      file_list <- list.files(datapath)
      print("File_list")
      print(file_list)

      counter <- 1
      plist <- list()
      results_all_df <- data.frame()
      tag <- FALSE

      for (file_name in file_list[1:length(file_list)]) {
        if (file_ext(file_name) == "csv") {
          results_df <- read.csv(paste0(datapath, "/", file_name))
          results_df <- cbind(results_df, data.frame("filename" = rep(file_name, nrow(results_df))))
          results_all_df <- rbind(results_all_df, results_df)
        }
      }

      if (length(results_all_df) > 1) {
        results_all_df$datapath <- datapath
        results_all_df$experimental_condition <- input$exp_condition
      }
    } else if ((input$data_upload_source == "key file") & !is.null(input$keyData)) {
      results_all_df <- data.frame()
      print("key data path")
      print(input$keyData$datapath)
      key_file <- read.csv(input$keyData$datapath, header = input$header_correlation)
      print(key_file)

      print("data_root")
      data_root <- key_data_info$datapath
      print(data_root)


      for (i in 1:nrow(key_file)) {


        # data_path <- key_file[i,"feature_table"]
        folder_path <- key_file[i, "folder_name"]
        short_name <- key_file[i, "short_name"]
        data_path <- paste0(data_root, "/", short_name, "/")

        print("Data path")
        print(data_path)

        if (dir.exists(data_path)) {
          file_list <- list.files(data_path)
          print("File_list")
          print(file_list)

          counter <- 1
          plist <- list()
          tag <- FALSE

          for (file_name in file_list[1:length(file_list)]) {
            print("File name")
            print(file_list)

            if (file_ext(file_name) == "csv") {
              results_df <- read.csv(paste0(data_path, "/", file_name))
              results_df <- cbind(results_df, data.frame("short_name" = rep(short_name, nrow(results_df))))
              results_all_df <- rbind(results_all_df, results_df)
            }
          }
        } else {
          print("Directory does not exist")
        }



        # if (length(results_all_df) > 1) {
        #    results_all_df$datapath <- datapath
        #    results_all_df$experimental_condition <- input$exp_condition
        # }


        # print("data_path")
        # print(data_path)

        # data_path <- paste0(data_path, "merged_table_", short_name, ".csv")
        # print("file_path")
        # print(data_path)

        # split_path <- SplitPath(input$keyData)$dirname
        # print(split_path)
        # results_df <- read.csv(data_path)
        # results_df <- cbind(results_df,  data.frame("filename"=rep(data_path, nrow(results_df))) )
        # results_all_df  <- rbind(results_all_df, results_df )
      }
    } else {
      results_all_df <- data.frame()
      # datapath = "../test_data/stack_EC_microscopy/120821 BSA #01.csv"
      # results_all_df <- read.csv(inFileStackData$datapath, header = input$header_correlation)
    }

    if (input$subsample_data) {
      N <- nrow(results_all_df) %/% input$subsample_n
      if (nrow(results_all_df) > N) {
        results_all_df <- results_all_df[sample(nrow(results_all_df), N), ]
      }
    }
  
    # TODO: remove this computation and move to python PolarityJaM
    for(i in 1:nrow(results_all_df)) {
      row <- results_all_df[i,]
      A <- row$cell_area
      P <- row$cell_perimeter
    
      cell_circularity <- 4*A/(pi*P*P)
      results_all_df[i,"cell_circularity"] = cell_circularity
      
      A <- row$nuc_area
      P <- row$nuc_perimeter
      nuc_circularity <- 4*A/(pi*P*P)
      results_all_df[i,"nuc_circularity"] = nuc_circularity
    } 


    # if (is.null(inFileStackData)) {
    #    datapath <- stack_data_info$datapath
    #
    #    # get list of files in the datapath
    #    file_list <- list.files(datapath)
    #    print("File_list")
    #    print(file_list)
    #
    #    counter <- 1
    #    plist <- list()
    #    results_all_df <-data.frame()
    #    tag <- FALSE
    #
    #    for (file_name in file_list[1:length(file_list)]){
    #
    #      if (file_ext(file_name) == "csv") {
    #        results_df <- read.csv(paste0(datapath,"/",file_name))
    #        results_df <- cbind(results_df,  data.frame("filename"=rep(file_name, nrow(results_df))) )
    #        results_all_df  <- rbind(results_all_df, results_df )
    #      }
    #    }
    #
    #    if (length(results_all_df) > 1) {
    #      results_all_df$datapath <- datapath
    #      results_all_df$experimental_condition <- input$exp_condition
    #    }
    # } else {
    #    results_all_df <- read.csv(inFileStackData$datapath, header = input$header_correlation)
    # }

#    sample_col <- input$sample_col
#    print(sample_col)
#    results_all_df <- results_all_df[!is.na(results_all_df[sample_col]),]
    results_all_df <- na.omit(results_all_df)
    results_all_df
  })
  # end of merged stack function


  output$terms_of_use_text <- renderText({
    "
    function that the merged stack of polarity data and angles in table format
    "
    if ((input$data_upload_form == "upload data") & (input$terms_of_use == FALSE)) {
      # if ((input$data_upload_form == "upload data")) {
      # HTML("Dear user, data upload is currently not possible in the online version. Please download the Rshiny app from <a href='https://github.com/wgiese/polarityjam'>polaritjam</a>! on your computer and run this app locally. </p>")
      # HTML("<p>If you enjoyed this tool, please consider <a href='https://www.gofundme.com/f/fantasy-football-mental-health-initiative?utm_medium=copy_link&utm_source=customer&utm_campaign=p_lico+share-sheet'>donating to the Fantasy Football Mental Health Initiative</a>!</p>")
      HTML("<p>  <font size='+2'> Terms of Use </font><br>
           Text </p>")
    } else {

    }
  })

  output$terms_of_use_text_all <- renderText({
    "
    function that the merged stack of polarity data and angles in table format
    "

    # if ((input$data_upload_form == "upload data")) {
    # HTML("Dear user, data upload is currently not possible in the online version. Please download the Rshiny app from <a href='https://github.com/wgiese/polarityjam'>polaritjam</a>! on your computer and run this app locally. </p>")
    # HTML("<p>If you enjoyed this tool, please consider <a href='https://www.gofundme.com/f/fantasy-football-mental-health-initiative?utm_medium=copy_link&utm_source=customer&utm_campaign=p_lico+share-sheet'>donating to the Fantasy Football Mental Health Initiative</a>!</p>")
    HTML("<p>  <font size='+2'> Terms of Use </font><br>
           Text </p>")
  })


  output$merged_stack <- renderTable({
    "
    function that the merged stack of polarity data and angles in table format
    "
    if ((input$data_upload_form == "upload data") & (input$terms_of_use == FALSE)) {
      # data.frame( "Info" = c("Dear user, data upload is currently not possible in the online version.",
      #                       "Please download the Rshiny app from \n and run locally."))
    } else {
      mergedStack()
    }
  })


  mergedStatistics <- reactive({
    "
    reactive function that reads a stack of spreadsheet and returns a data frame 
    with descriptive statistics including circular mean, circular standard deviation 
    and nearest neighbours for the merged stack of data

    TODO: rework threhsolding

    "

    results_df <- mergedStack()


    print("Data Frame in merged statistics:")
    print(head(results_df))

    source(file = paste0(getwd(), "/src/circular_statistics.R"), local = T)
    parameters <- fromJSON(file = "parameters/parameters.json")

    condition_col <- input$condition_col
    condition_list <- unlist(unique(results_df[condition_col]))

    feature <- parameters[input$feature_select][[1]][1]

    #    threshold <- input$min_nuclei_golgi_dist
    #    if ("organelle_distance" %in% colnames(results_df)){
    #      results_df <- subset(results_df, results_df$distance > threshold)
    #    }

    statistics_df <- as.data.frame(matrix(ncol = length(condition_list) + 2, nrow = 0))
    cols <- c("entity")
    for (condition in condition_list) {
      cols <- c(cols, condition)
    }
    cols <- c(cols, "description")


    colnames(statistics_df) <- cols # c("entity", "value") #, "comment")

    print("Colnames")
    print(colnames(statistics_df))
    # print("Feature property")
    # print(parameters[input$feature_select][[1]][2])
    if (parameters[input$feature_select][[1]][2] == "directional") {

      for (condition in condition_list) {
        condition_data <- subset(results_df, results_df[condition_col] == condition)
        print("Condition subset: ")
        print(head(condition_data))

        x_data <- unlist(condition_data[feature]) * 180.0 / pi
        statistics <- compute_circular_statistics(condition_data, feature, parameters)
        print("Statistics")
        print(statistics)

        p_value <- signif(statistics[1, "rayleigh_test"], digits = 3)
        # if (statistics[1,"rayleigh_test"] < 0.001)
        #    p_value <- "p < 0.001"
        p_value_mu <- signif(statistics[1, "v_test"], digits = 3)
        # if (statistics[1,"rayleigh_test_mu"] < 0.001)
        #    p_value_mu <- "p < 0.001"
        ind <- 1
        statistics_df[ind, 1] <- "number of cells"
        statistics_df[ind, condition] <- nrow(condition_data)

        ind <- ind + 1
        statistics_df[ind, 1] <- "mean (degree)"
        statistics_df[ind, condition] <- signif(statistics[1, "mean"], digits = 3)

        ind <- ind + 1
        statistics_df[ind, 1] <- "polarity index"
        statistics_df[ind, condition] <- signif(statistics[1, "polarity_index"], digits = 3)

        ind <- ind + 1
        statistics_df[ind, 1] <- "signed polarity index (mean = 180)"
        statistics_df[ind, condition] <- signif(statistics[1, "signed_polarity_index"], digits = 3)

        ind <- ind + 1
        statistics_df[ind, 1] <- "angular standard deviation"
        statistics_df[ind, condition] <- signif(statistics[1, "std_angular"], digits = 3)
        # statistics_df[ind,3] <- "angular standard deviation, takes values in [0,sqrt(2)], see https://doi.org/10.18637/jss.v031.i10 for more info."

        ind <- ind + 1
        statistics_df[ind, 1] <- "circular standard deviation"
        statistics_df[ind, condition] <- signif(statistics[1, "std_circular"], digits = 3)
        # statistics_df[ind,3] <- "circular standard deviation, takes values in [0,inf], see https://doi.org/10.18637/jss.v031.i10 for more info."

        ind <- ind + 1
        statistics_df[ind, 1] <- "95% confidence interval of the mean, lower limit: "
        statistics_df[ind, condition] <- signif(statistics[1, "ci_95_lower_limit"], digits = 3)

        ind <- ind + 1
        statistics_df[ind, 1] <- "95% confidence interval of the mean, upper limit: "
        statistics_df[ind, condition] <- signif(statistics[1, "ci_95_upper_limit"], digits = 3)

        ind <- ind + 1
        statistics_df[ind, 1] <- "Rayleigh test, p-value:"
        statistics_df[ind, condition] <- p_value

        ind <- ind + 1
        statistics_df[ind, 1] <- "V-test p-value (cond. mean = 180): "
        statistics_df[ind, condition] <- p_value_mu
      }
    } else if (parameters[input$feature_select][[1]][2] == "undirectional") {
      #statistics <- compute_undirectional_statistics(results_df, feature, parameters)

      #p_value <- signif(statistics[1, "rayleigh_test"], digits = 3)

      #statistics_df[1, 1] <- "cells"
      #statistics_df[1, 2] <- nrow(results_df)
      #statistics_df[2, 1] <- "mean (degree)"
      #statistics_df[2, 2] <- signif(statistics[1, "mean"], digits = 3)
      #statistics_df[3, 1] <- "polarity index"
      #statistics_df[3, 2] <- signif(statistics[1, "polarity_index"], digits = 3)
      #statistics_df[4, 1] <- "Rayleigh test, p-value:"
      #statistics_df[4, 2] <- p_value
      for (condition in condition_list) {
        condition_data <- subset(results_df, results_df[condition_col] == condition)
        print("Condition subset: ")
        print(head(condition_data))
        
        x_data <- unlist(condition_data[feature]) * 180.0 / pi
        statistics <- compute_undirectional_statistics(condition_data, feature, parameters)
        print("Statistics")
        print(statistics)
        
        p_value <- signif(statistics[1, "rayleigh_test"], digits = 3)
        # if (statistics[1,"rayleigh_test"] < 0.001)
        #    p_value <- "p < 0.001"

        ind <- 1
        statistics_df[ind, 1] <- "number of cells"
        statistics_df[ind, condition] <- nrow(condition_data)
        
        ind <- ind + 1
        statistics_df[ind, 1] <- "mean (degree)"
        statistics_df[ind, condition] <- signif(statistics[1, "mean"], digits = 3)
        
        ind <- ind + 1
        statistics_df[ind, 1] <- "polarity index"
        statistics_df[ind, condition] <- signif(statistics[1, "polarity_index"], digits = 3)
        
        #ind <- ind + 1
        #statistics_df[ind, 1] <- "signed polarity index (mean = 180)"
        #statistics_df[ind, condition] <- signif(statistics[1, "signed_polarity_index"], digits = 3)
        
        ind <- ind + 1
        statistics_df[ind, 1] <- "angular standard deviation"
        statistics_df[ind, condition] <- signif(statistics[1, "std_angular"], digits = 3)
        # statistics_df[ind,3] <- "angular standard deviation, takes values in [0,sqrt(2)], see https://doi.org/10.18637/jss.v031.i10 for more info."
        
        ind <- ind + 1
        statistics_df[ind, 1] <- "circular standard deviation"
        statistics_df[ind, condition] <- signif(statistics[1, "std_circular"], digits = 3)
        # statistics_df[ind,3] <- "circular standard deviation, takes values in [0,inf], see https://doi.org/10.18637/jss.v031.i10 for more info."
        
        ind <- ind + 1
        statistics_df[ind, 1] <- "95% confidence interval of the mean, lower limit: "
        statistics_df[ind, condition] <- signif(statistics[1, "ci_95_lower_limit"], digits = 3)
        
        ind <- ind + 1
        statistics_df[ind, 1] <- "95% confidence interval of the mean, upper limit: "
        statistics_df[ind, condition] <- signif(statistics[1, "ci_95_upper_limit"], digits = 3)
        
        ind <- ind + 1
        statistics_df[ind, 1] <- "Rayleigh test, p-value:"
        statistics_df[ind, condition] <- p_value
        
      }
      
    } else {
      statistics <- compute_linear_statistics(results_df, feature, parameters)

      statistics_df[1, 1] <- "cells"
      statistics_df[1, 2] <- nrow(results_df)
      statistics_df[2, 1] <- "mean"
      statistics_df[2, 2] <- signif(statistics[1, "mean"], digits = 3)
      statistics_df[3, 1] <- "standard deviation"
      statistics_df[3, 2] <- signif(statistics[1, "std"], digits = 3)
      statistics_df[4, 1] <- "median"
      statistics_df[4, 2] <- signif(statistics[1, "median"], digits = 3)
    }

    statistics_df
  })

  output$merged_statistics <- renderTable(
    {
      "
    function that shows the descriptive statistics of the merged data stack in table format
    "

      statistics_df <- mergedStatistics()
      statistics_df
    },
    digits = 3
  )

  merged_plot <- reactive({
    source(file = paste0(getwd(), "/src/plot_functions.R"), local = T)
    source(file = paste0(getwd(), "/src/circular_statistics.R"), local = T)

    parameters <- fromJSON(file = "parameters/parameters.json")
    text_size <- as.integer(parameters["text_size_merged_plot"])

    results_all_df <- mergedStack()

    # inFileStackData <- input$stackData

    # if (!is.null(inFileStackData))
    #    results_all_df <- read.csv(inFileStackData$datapath, header = input$header_correlation)

    # TODO: make filtering conditional
    # for(i in 1:nrow(results_all_df)) {
    #  row <- results_all_df[i,]
    #  a <- row$major_axis_length
    #  b <- row$minor_axis_length
    #
    #  eccentricity <- sqrt(1.0 - b*b/(a*a))
    #  results_all_df[i,"cell_eccentricity"] = eccentricity
    # }

    # threshold <- input$min_nuclei_golgi_dist
    # if ("organelle_distance" %in% colnames(results_all_df)){
    #  results_all_df <- subset(results_all_df, results_all_df$distance > threshold)
    # }

    # print("In merged_plot")
    # print(head(results_all_df))

    bin_size <- 360 / input$bins
    exp_condition <- input$exp_condition
    datapath <- stack_data_info$datapath

    feature <- parameters[input$feature_select][[1]][1]

    print("Feature:")
    print(feature)

    if (parameters[input$feature_select][[1]][2] == "directional") {
      print("directional feature!")


      x_data <- unlist(results_all_df[feature]) * 180.0 / pi
      statistics <- compute_circular_statistics(results_all_df, feature, parameters)
      plot_title <- parameters[input$feature_select][[1]][3]
      p <- rose_plot_circular(parameters, input, statistics, x_data, plot_title, 0, text_size)
    } else if (parameters[input$feature_select][[1]][2] == "undirectional") {
      x_data <- results_all_df[feature]
      statistics <- compute_undirectional_statistics(results_all_df, feature, parameters)
      # if (input$left_directional) {
      #  x_data <- unlist(transform_undirectional(input,x_data))*180.0/pi
      # } else {
      #  x_data <- unlist(results_all_df[feature])*180.0/pi
      # }
      x_data <- unlist(transform_undirectional(input, x_data)) * 180.0 / pi

      plot_title <- parameters[input$feature_select][[1]][3]
      p <- rose_plot_undirectional(parameters, input, statistics, x_data, plot_title, 0, text_size)
    } else {
      x_data <- unlist(results_all_df[feature])
      statistics <- compute_linear_statistics(results_all_df, feature, parameters)
      plot_title <- parameters[input$feature_select][[1]][3]
      p <- linear_histogram(parameters, input, statistics, x_data, plot_title, 0, text_size, min(x_data), max(x_data))
    }

    p
  })

  width_A <- reactive({
    input$plot_width_A
  })
  height_A <- reactive({
    input$plot_height_A
  })

  output$merged_plot <- renderPlot(width = width_A, height = height_A, {
    parameters <- fromJSON(file = "parameters/parameters.json")
    # parameters[input$feature_select][[1]][1]

    if (input$feature_select %in% names(parameters)) {
      p <- merged_plot()
      p
    } else {

    }

    # if (input$feature_select == "filename") {
    #
    #        }
    #        else {
    #            p <-merged_plot()
    #            p
    #        }
  })

  output$parameter_error <- renderText({
    parameters <- fromJSON(file = "parameters/parameters.json")
    # parameters[input$feature_select][[1]][1]
    # if (input$feature_select == "filename") {
    #    print("Plotting of this parameter is not supported.")
    # }
    # else {
    #
    # }

    if (input$feature_select %in% names(parameters)) {

    } else {
      print("Plotting of this parameter is not supported.")
    }
  })


  multi_plot <- reactive({
    source(file = paste0(getwd(), "/src/plot_functions.R"), local = T)
    source(file = paste0(getwd(), "/src/circular_statistics.R"), local = T)

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

    #   for(row_nr in 1:nrow(results_all_df)) {
    #       row <- results_all_df[row_nr,]
    # a <- row$major_axis_length
    # b <- row$minor_axis_length
    #
    # eccentricity <- sqrt(1.0 - b*b/(a*a))
    # results_all_df[row_nr,"cell_eccentricity"] = eccentricity
    # }

    # threshold <- input$min_eccentricity
    # if ("cell_eccentricity" %in% colnames(results_all_df)){
    #  results_all_df <- subset(results_all_df, results_all_df$eccentricity> threshold)
    # }
    # threshold <- input$min_nuclei_golgi_dist
    # if ("orgenelle_distance" %in% colnames(results_all_df)){
    #  results_all_df <- subset(results_all_df, results_all_df$distance > threshold)
    # }

    feature <- parameters[input$feature_select][[1]][1]
    condition_col <- input$condition_col

    condition_list <- unlist(unique(results_all_df[condition_col]))
    # plist <- vector('list', length(unique(results_all_df$filename)))
    plist <- vector("list", length(condition_list))
    print("length of plot list")
    print(plist)
    print("list of unique entries")
    print(unlist(unique(results_all_df[condition_col])))

    x_lim <- c(min(results_all_df[feature]), max(results_all_df[feature]))
    # for(file_name in unique(results_all_df$filename)) {
    #  results_df <- subset(results_all_df, results_all_df$filename == file_name )
    for (file_name in condition_list) {
      results_df <- subset(results_all_df, results_all_df[condition_col] == file_name)

      # values <- compute_polarity_index(results_df)


      # print(values)
      # polarity_index <- values[["polarity_index"]]
      # angle_mean_deg <- values[["angle_mean_deg"]]

      x <- unlist(results_df[feature])
      angle_dists[[i]] <- x


      # if (parameters[input$feature_select][[1]][2] == "linear") {
      #
      #    if (x_lim[0] > min(x))
      #        x_lim[0] <- min(x)
      #    if (x_lim[1] < max(x))
      #        x_lim[1] <- max(x)
      #
      #       }

      file_names[[i]] <- file_name

      # polarity_indices[[i]] <- polarity_index
      # angle_mean_degs[[i]] <- angle_mean_deg
      i <- i + 1
    }



    n <- length(angle_dists)
    nCol <- floor(sqrt(n))

    bin_size <- 360 / input$bins

    plotseries <- function(i) {
      angle_dist <- angle_dists[[i]]
      file_name <- file_names[[i]]
      # polarity_index <- polarity_indices[[i]]
      # angle_mean_deg <- angle_mean_degs[[i]]

      # results_df <- subset(results_all_df, results_all_df$filename == file_name)
      results_df <- subset(results_all_df, results_all_df[condition_col] == file_name)

      plot_title <- file_name

      if (nchar(file_name) > 37) {
        max_fl <- 17
        file_name_end <- substr(file_name, nchar(file_name) - max_fl + 1, nchar(file_name))
        file_name_start <- substr(file_name, 1, max_fl)
        plot_title <- paste0(file_name_start, "...", file_name_end)
      }
      # if (nchar(file_name) > 15) {
      #    plot_title <- paste0("image #",toString(i))
      #    print(paste0("filename: ",file_name," too long, will be replaced by",plot_title))
      # }


      if (parameters[input$feature_select][[1]][2] == "directional") {
        statistics <- compute_circular_statistics(results_df, feature, parameters)
        # statistics <- compute_polarity_index(unlist(results_df[feature]))
        x_data <- unlist(results_df[feature]) * 180.0 / pi
        print(paste0("Length of filename", toString(i)))

        p <- rose_plot_circular(parameters, input, statistics, x_data, plot_title, i, text_size)
      } else if (parameters[input$feature_select][[1]][2] == "undirectional") {
        x_data <- results_df[feature]
        # print(x_data)
        statistics <- compute_undirectional_statistics(results_df, feature, parameters)
        # if (input$left_directional) {
        x_data <- unlist(transform_undirectional(input, x_data)) * 180.0 / pi
        # } else {
        #  x_data <- unlist(results_df[feature])*180.0/pi
        # }
        # plot_title <- file_name
        p <- rose_plot_undirectional(parameters, input, statistics, x_data, plot_title, i, text_size)
      } else {
        x_data <- unlist(results_df[feature])
        statistics <- compute_linear_statistics(results_df, feature, parameters)
        # plot_title <- file_name
        # p <- linear_histogram(parameters, input, statistics, x_data,  plot_title, i, text_size, x_lim[0], x_lim[1])
        p <- linear_histogram(parameters, input, statistics, x_data, plot_title, i, text_size, min(results_all_df[feature]), max(results_all_df[feature]))
      }
    }


    myplots <- lapply(1:length(angle_dists), plotseries)

    # print(myplots)
    grid.arrange(grobs = myplots, nrow = nCol) # , widths = list(10,10))
  })

  output$multi_dist_plot <- renderPlot(width = width_A, height = height_A, {
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
      if (input$dataset == "statistics_file") {
        filename <- "statistics_file.csv"
        print("Download merged_file.csv")
      }
      if (input$dataset == "merged_plot_file") {
        filename <- paste0("merge_plot", input$image_file_format)
      }
      if (input$dataset == "multi_plot_file") {
        filename <- paste0("multi_plot", input$image_file_format)
      }
      return(filename)
    },
    content = function(file) {
      parameters <- fromJSON(file = "parameters/parameters.json")
      width_ <- as.double(parameters["pdf_figure_size_inches"])

      if (input$dataset == "statistics_file") {
        return(write.csv(mergedStatistics(), file, row.names = FALSE))
      } else if ((input$dataset == "multi_plot_file") && (input$image_file_format == ".pdf")) {
        # pdf(file, width=14, height=14)
        pdf(file, family = "ArialMT", width = width_, height = width_, pointsize = 18)
        p <- multi_plot()
        plot(p)
        dev.off()
      } else if ((input$dataset == "multi_plot_file") && (input$image_file_format == ".png")) {
        png(file, width = 960, height = 960)
        p <- multi_plot()
        plot(p)
        dev.off()
      } else if ((input$dataset == "multi_plot_file") && (input$image_file_format == ".eps")) {
        eps(file, width = 14, height = 14)
        p <- multi_plot()
        plot(p)
        dev.off()
      } else if ((input$dataset == "merged_plot_file") && (input$image_file_format == ".pdf")) {
        print("Saving merge pdf")
        pdf(file, family = "ArialMT", width = width_, height = width_, pointsize = 18)
        p <- merged_plot()
        plot(p)
        dev.off()
      } else if ((input$dataset == "merged_plot_file") && (input$image_file_format == ".png")) {
        png(file, width = 960, height = 960)
        p <- merged_plot()
        plot(p)
        dev.off()
      } else if ((input$dataset == "merged_plot_file") && (input$image_file_format == ".eps")) {
        eps(file, width = width_, height = width_)
        p <- merged_plot()
        plot(p)
        dev.off()
      } else {
        (
          return(write.csv(mergedStack(), file, row.names = FALSE))
        )
      }
    }
  )

  output$downloadDataSingleImage <- downloadHandler(
    filename = function() {
      if (input$datasetSingleImage == "results_file") {
        filename <- "results_file.csv"
      }
      if (input$datasetSingleImage == "statistics_file") {
        filename <- "statistics_file.csv"
      }
      if (input$datasetSingleImage == "rose_histogram") {
        filename <- "rose_histogram.pdf"
      }
      if (input$datasetSingleImage == "orientation_plot") {
        filename <- "vector_plot.pdf"
      }
      filename
    },
    content = function(file) {
      if (input$datasetSingleImage == "results_file") {
        write.csv(resultSingleImage(), file, row.names = FALSE)
      }
      if (input$datasetSingleImage == "statistics_file") {
        write.csv(singleImageStatistics(), file, row.names = FALSE)
      }
      if (input$datasetSingleImage == "rose_histogram") {
        pdf(file, width = 7, height = 7)
        p <- rose_histogram_single_image()
        plot(p)
        dev.off()
      }
      if (input$datasetSingleImage == "orientation_plot") {
        pdf(file, width = 7, height = 7)
        p <- vectorPlot()
        plot(p)
        dev.off()
      }
    }
  )


  # download for merged plot

  output$downloadMergedPlotPDF <- downloadHandler(
    filename <- function() {
      paste(paste0("PolarityJaM_", input$featue_select,"_Merged_"), Sys.time(), ".pdf", sep = "")
    },
    content <- function(file) {
      pdf(file, width = input$plot_width_A / 72, height = input$plot_height_A / 72)
      plot(merged_plot())
      dev.off()
    },
    contentType = "application/pdf" # MIME type of the image
  )

  output$downloadMergedPlotSVG <- downloadHandler(
    filename <- function() {
      paste("PolarityJaM_Merged_", Sys.time(), ".svg", sep = "")
    },
    content <- function(file) {
      svg(file, width = input$plot_width_A / 72, height = input$plot_height_A / 72)
      plot(merged_plot())
      dev.off()
    },
    contentType = "application/svg" # MIME type of the image
  )

  output$downloadMergedPlotEPS <- downloadHandler(
    filename <- function() {
      paste("PolarityJaM_Merged_", Sys.time(), ".eps", sep = "")
    },
    content <- function(file) {
      cairo_ps(file, width = input$plot_width_A / 72, height = input$plot_height_A / 72)
      plot(merged_plot())
      dev.off()
    },
    contentType = "application/eps" # MIME type of the image
  )

  output$downloadMergedPlotPNG <- downloadHandler(
    filename <- function() {
      paste("PolarityJaM_", input$feature_select, "_Merged_", Sys.time(), ".png", sep = "")
    },
    content <- function(file) {
      png(file, width = input$plot_width_A * 4, height = input$plot_height_A * 4, res = 300)
      # if (input$data_form != "dataaspixel") plot(plot_data())
      # else plot(plot_map())
      plot(merged_plot())
      dev.off()
    },
    contentType = "application/png" # MIME type of the image
  )

  # download for multi plot
  # TODO: check why multi plot files have a grid when downloaded, while no grid is displayed in the app


  output$downloadMultiPlotPDF <- downloadHandler(
    filename <- function() {
      paste("PolarityJaM_Multi_", Sys.time(), ".pdf", sep = "")
    },
    content <- function(file) {
      pdf(file, width = input$plot_width_A / 72, height = input$plot_height_A / 72)
      plot(multi_plot())
      dev.off()
    },
    contentType = "application/pdf" # MIME type of the image
  )

  output$downloadMultiPlotSVG <- downloadHandler(
    filename <- function() {
      paste("PolarityJaM_Multi_", Sys.time(), ".svg", sep = "")
    },
    content <- function(file) {
      svg(file, width = input$plot_width_A / 72, height = input$plot_height_A / 72)
      plot(multi_plot())
      dev.off()
    },
    contentType = "application/svg" # MIME type of the image
  )

  output$downloadMultiPlotEPS <- downloadHandler(
    filename <- function() {
      paste("PolarityJaM_Multi_", Sys.time(), ".eps", sep = "")
    },
    content <- function(file) {
      cairo_ps(file, width = input$plot_width_A / 72, height = input$plot_height_A / 72)
      plot(multi_plot())
      dev.off()
    },
    contentType = "application/eps" # MIME type of the image
  )

  output$downloadMultiPlotPNG <- downloadHandler(
    filename <- function() {
      paste("PolarityJaM_", input$feature_select, "_Multi_", Sys.time(), ".png", sep = "")
    },
    content <- function(file) {
      png(file, width = input$plot_width_A * 4, height = input$plot_height_A * 4, res = 300)
      # if (input$data_form != "dataaspixel") plot(plot_data())
      # else plot(plot_map())
      plot(multi_plot())
      dev.off()
    },
    contentType = "application/png" # MIME type of the image
  )


  ### Panel B

  plot_correlation <- reactive({
    
    
    parameters <- fromJSON(file = "parameters/parameters.json")
    source(file = paste0(getwd(), "/src/plot_functions.R"), local = T)
    source(file = paste0(getwd(), "/src/circular_statistics.R"), local = T)
    source(file = paste0(getwd(), "/src/circular_correlations.R"), local = T)
    
    text_size <- input$text_size_corr

    # inFileCorrelationData <- input$correlationData

    # if (is.null(inFileCorrelationData))
    #    return(NULL)

    # print(inFileCorrelationData$datapath)
    # correlation_data <- read.csv(inFileCorrelationData$datapath, header = input$header_correlation)

    correlation_data <- mergedStack() # read.csv(inFileCorrelationData$datapath, header = input$header_correlation)

    
    
    feature_1 <- parameters[input$feature_select_1][[1]][1]
    feature_2 <- parameters[input$feature_select_2][[1]][1]
    
    
    feature_1_values <- unlist(correlation_data[feature_1])
    feature_1_values_ <- correlation_data[feature_1] * 180.0 / pi
    feature_2_values <- unlist(correlation_data[feature_2])
    feature_2_values_ <- correlation_data[feature_2] * 180.0 / pi

    # feature_1_values_sin <- sin(unlist(correlation_data[feature_1]))
    # feature_2_values_sin <- sin(unlist(correlation_data[feature_2]))

    feature_1_name <- parameters[input$feature_select_1][[1]][3]
    feature_2_name <- parameters[input$feature_select_2][[1]][3]

    conditions <- correlation_data[input$condition_col]
    
    print("feature_values")

    #if (parameters[input$feature_select_1][[1]][2] != "linear" && parameters[input$feature_select_2][[1]][2] != "linear") {
      
    p <- plot_circular_circular(correlation_data, input, parameters, plot_nr = 0, text_size = 24) 
      
 
  })
  
  output$correlation_statistics <- renderTable(
    {
      "
    function that shows the descriptive statistics of the merged data stack in table format
    "
      
      correlation_data <- mergedStack()
      
      source(file = paste0(getwd(), "/src/circular_correlations.R"), local = T)
      parameters <- fromJSON(file = "parameters/parameters.json")
      
      condition_col <- input$condition_col
      condition_list <- unlist(unique(correlation_data[condition_col]))
      
      feature_1 <- parameters[input$feature_select_1][[1]][1]
      feature_2 <- parameters[input$feature_select_2][[1]][1]
      
      feature_1_values <- unlist(correlation_data[feature_1])
      
      feature_2_values <- unlist(correlation_data[feature_2])
      
      
      feature_1_name <- parameters[input$feature_select_1][[1]][3]
      feature_2_name <- parameters[input$feature_select_2][[1]][3]
      
      mode_1 <- parameters[input$feature_select_1][[1]][2]
      mode_2 <- parameters[input$feature_select_2][[1]][2]
      
      
      #    threshold <- input$min_nuclei_golgi_dist
      #    if ("organelle_distance" %in% colnames(results_df)){
      #      results_df <- subset(results_df, results_df$distance > threshold)
      #    }
      
      statistics_df <- as.data.frame(matrix(ncol = length(condition_list) + 3, nrow = 0))
      cols <- c("entity")
      for (condition in condition_list) {
        cols <- c(cols, condition)
      }
      cols <- c(cols, "all", "description")
      
      
      colnames(statistics_df) <- cols # c("entity", "value") #, "comment")
      
      res <- compute_correlation(feature_1_values, mode_1, feature_2_values, mode_2) 
      
      ind <- 1
      statistics_df[ind, 1] <- "pearson r-value"
      if ( (mode_1 == "linear") | (mode_2 == "linear") ) {
        statistics_df[ind, "all"] <- res
      } else {
        statistics_df[ind, "all"] <- res$r
      }
      
      
      print("Colnames")
      print(colnames(statistics_df))

        for (condition in condition_list) {
          condition_data <- subset(correlation_data, correlation_data[condition_col] == condition)
          
          feature_1_values <- unlist(condition_data[feature_1])
          feature_2_values <- unlist(condition_data[feature_2])
          
          res <- compute_correlation(feature_1_values, mode_1, feature_2_values, mode_2) 
          
          ind <- 1
          statistics_df[ind, 1] <- "pearson r-value"
          if ( (mode_1 == "linear") | (mode_2 == "linear") ) {
            statistics_df[ind, condition] <- res
          } else {
            statistics_df[ind, condition] <- res$r
          }
          
        }
      statistics_df
      #statistics_df <- mergedStatistics()
      #statistics_df
    },
    digits = 3
  )

  width <- reactive({
    input$plot_width_corr
  })
  height <- reactive({
    input$plot_height_corr
  })

  output$correlation_plot <- renderPlot(width = width, height = height, {
    if (input$corr_plot_option == "spoke plot") {
      p <- spoke_plot_correlation()
    } else {
      p <- plot_correlation()
    }
    p
  })

  spoke_plot_correlation <- reactive({
    parameters <- fromJSON(file = "parameters/parameters.json")

    text_size <- input$text_size_corr

    correlation_data <- mergedStack()

    if (input$spoke_subsample_n > 1) {
      N <- nrow(correlation_data) %/% input$spoke_subsample_n
      if (nrow(correlation_data) > N) {
        correlation_data <- correlation_data[sample(nrow(correlation_data), N), ]
      }
    }


    feature_1 <- parameters[input$feature_select_1][[1]][1]
    feature_2 <- parameters[input$feature_select_2][[1]][1]
    feature_1_values <- unlist(correlation_data[feature_1])
    feature_2_values <- unlist(correlation_data[feature_2])

    feature_1_name <- parameters[input$feature_select_1][[1]][3]
    feature_2_name <- parameters[input$feature_select_2][[1]][3]

    # res = circ.cor(feature_1_values, feature_2_values, test=TRUE)

    # reg_coeff <- res$r
    # p_value <- res$p.value

    feature_1_values_deg <- unlist(correlation_data[feature_1]) * 180.0 / pi
    feature_2_values_deg <- unlist(correlation_data[feature_2]) * 180.0 / pi

    feature_1_x_a <- list()
    feature_1_y_a <- list()
    feature_1_x_b <- list()
    feature_1_y_b <- list()


    feature_2_x_a <- list()
    feature_2_y_a <- list()
    feature_2_x_b <- list()
    feature_2_y_b <- list()


    print("until here")
    for (i in 1:length(feature_1_values)) {
      #    print(i)
      #    print(feature_1_values[i])
      feature_1_x_a[i] <- 0.5 * cos(feature_1_values[i])
      feature_1_y_a[i] <- 0.5 * sin(feature_1_values[i])
      feature_1_x_b[i] <- 0.5 * cos(feature_1_values[i] + pi)
      feature_1_y_b[i] <- 0.5 * sin(feature_1_values[i] + pi)


      # dist_a <- abs(feature_1_values[i] - feature_2_values[i])
      # dist_b <- abs(feature_1_values[i] - feature_2_values[i] - pi)

      dist_a <- (cos(feature_1_values[i]) - cos(feature_2_values[i])) * (cos(feature_1_values[i]) - cos(feature_2_values[i]))
      dist_a <- dist_a + (sin(feature_1_values[i]) - sin(feature_2_values[i])) * (sin(feature_1_values[i]) - sin(feature_2_values[i]))

      dist_b <- (cos(feature_1_values[i]) - cos(feature_2_values[i] + pi)) * (cos(feature_1_values[i]) - cos(feature_2_values[i] + pi))
      dist_b <- dist_b + (sin(feature_1_values[i]) - sin(feature_2_values[i] + pi)) * (sin(feature_1_values[i]) - sin(feature_2_values[i] + pi))


      if (dist_a < dist_b) {
        feature_2_x_a[i] <- cos(feature_2_values[i])
        feature_2_y_a[i] <- sin(feature_2_values[i])
        feature_2_x_b[i] <- cos(feature_2_values[i] + pi)
        feature_2_y_b[i] <- sin(feature_2_values[i] + pi)
      } else {
        feature_2_x_a[i] <- cos(feature_2_values[i] + pi)
        feature_2_y_a[i] <- sin(feature_2_values[i] + pi)
        feature_2_x_b[i] <- cos(feature_2_values[i])
        feature_2_y_b[i] <- sin(feature_2_values[i])
      }
    }


    f_1 <- data.frame(x1 = unlist(feature_1_x_a), y1 = unlist(feature_1_y_a))
    f_2 <- data.frame(x2 = unlist(feature_2_x_a), y2 = unlist(feature_2_y_a))
    f_1 <- data.frame(x1 = unlist(feature_1_x_b), y1 = unlist(feature_1_y_b))
    f_2 <- data.frame(x2 = unlist(feature_2_x_b), y2 = unlist(feature_2_y_b))


    print(head(f_1))

    p <- ggplot()
    # p <- p + geom_point(aes(x = list(feature_2_x), y = list(feature_2_y), size = 3))
    # p <- p + geom_point(aes(x = list(feature_1_x), y = list(feature_1_y), size = 3))
    # p <- p + geom_point(aes(x = x1, y = y1, size = 3))
    # p <- p + geom_point(aes(x = x2, y = y2, size = 3))
    p <- p + geom_point(aes(x = unlist(feature_2_x_a), y = unlist(feature_2_y_a), size = 3))
    p <- p + geom_point(aes(x = unlist(feature_1_x_a), y = unlist(feature_1_y_a), size = 3))
    p <- p + geom_point(aes(x = unlist(feature_2_x_b), y = unlist(feature_2_y_b), size = 3))
    p <- p + geom_point(aes(x = unlist(feature_1_x_b), y = unlist(feature_1_y_b), size = 3))



    p <- p + geom_segment(aes(x = unlist(feature_1_x_a), y = unlist(feature_1_y_a), xend = unlist(feature_2_x_a), yend = unlist(feature_2_y_a), size = 0.1, color = "red"))
    p <- p + geom_segment(aes(x = unlist(feature_1_x_b), y = unlist(feature_1_y_b), xend = unlist(feature_2_x_b), yend = unlist(feature_2_y_b), size = 0.1, color = "red"))

    p <- p + xlim(-1.0, 1.0)
    p <- p + ylim(-1.0, 1.0)

    p <- p + xlab(feature_1_name) + ylab(feature_2_name)
    p <- p + theme(aspect.ratio = 3 / 3)
    p <- p + geom_point(color = "black", size = input$marker_size_corr)
    # p <- p + theme_minimal(panel.background = element_blank(), base_size = text_size)
    p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))
    # p <- ggplot()
    # p <- p + geom_point(aes(x = feature_2_values_deg, y = 1, size = 3))
    # p <- p + geom_point(aes(x = feature_1_values_deg, y = 0.5, size = 3))

    # p <- p + geom_segment(aes(x=feature_1_values_deg, y=0.5, xend=feature_2_values_deg, yend=1.0, size = 0.1, color="red"))

    #        p <- p + ggtitle("spoke_plot") +
    #            theme(plot.title = element_text(size = 18, face = "bold")) +
    #            theme(axis.text.x = element_text(size = 18)) +
    #            coord_polar(start = -pi/2.0, direction = -1) +
    #            scale_x_continuous(limits = c(0, 180),
    #                       breaks = (c(0, 45, 90, 135))) +
    #            scale_x_continuous(limits = c(0, 360),
    #                       breaks = (c(0, 90, 180, 270))) +
    #            scale_y_continuous(limits = c(0, 1.1)) +
    #            theme_minimal(base_size = text_size)
    p
  })

  output$spoke_plot <- renderPlot({
    p <- spoke_plot_correlation()
    p
  })

  output$downloadPlotPDF <- downloadHandler(
    filename <- function() {
      paste("PolarityJaM_Correlation_", Sys.time(), ".pdf", sep = "")
    },
    content <- function(file) {
      pdf(file, width = input$plot_width_corr / 72, height = input$plot_height_corr / 72)
      plot(plot_correlation())
      dev.off()
    },
    contentType = "application/pdf" # MIME type of the image
  )

  output$downloadPlotSVG <- downloadHandler(
    filename <- function() {
      paste("PolarityJaM_Correlation_", Sys.time(), ".svg", sep = "")
    },
    content <- function(file) {
      svg(file, width = input$plot_width_corr / 72, height = input$plot_height_corr / 72)
      plot(plot_correlation())
      dev.off()
    },
    contentType = "application/svg" # MIME type of the image
  )

  output$downloadPlotEPS <- downloadHandler(
    filename <- function() {
      paste("PolarityJaM_Correlation_", Sys.time(), ".eps", sep = "")
    },
    content <- function(file) {
      cairo_ps(file, width = input$plot_width_corr / 72, height = input$plot_height_corr / 72)
      plot(plot_correlation())
      dev.off()
    },
    contentType = "application/eps" # MIME type of the image
  )

  output$downloadPlotPNG <- downloadHandler(
    filename <- function() {
      paste("PolarityJaM_Correlation_", Sys.time(), ".png", sep = "")
    },
    content <- function(file) {
      png(file, width = input$plot_width_corr * 4, height = input$plot_height_corr * 4, res = 300)
      # if (input$data_form != "dataaspixel") plot(plot_data())
      # else plot(plot_map())
      plot(plot_correlation())
      dev.off()
    },
    contentType = "application/png" # MIME type of the image
  )



  ### Panel C

  observe({
    condition_col <- input$condition_col
    condition_col <- "filename"

    data <- mergedStack()

    # if (is.data.frame(data)) {
    if (length(colnames(data)) > 3) {
      print("Unique condition names: ")
      # print(data[condition_col])
      condition_list <- unique(data[condition_col])
      print(condition_list)
      updateSelectInput(session, "control_condition", choices = condition_list, selected = "filename")
    }
  })



  comparison_plot <- reactive({
    source(file = paste0(getwd(), "/src/plot_functions.R"), local = T)
    source(file = paste0(getwd(), "/src/circular_statistics.R"), local = T)

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

    #        for(row_nr in 1:nrow(results_all_df)) {
    #            row <- results_all_df[row_nr,]
    #      a <- row$major_axis_length
    #      b <- row$minor_axis_length
    #
    #      eccentricity <- sqrt(1.0 - b*b/(a*a))
    #      results_all_df[row_nr,"eccentricity"] = eccentricity
    #    }

    #    threshold <- input$min_eccentricity
    #    if ("eccentricity" %in% colnames(results_all_df)){
    #      results_all_df <- subset(results_all_df, results_all_df$eccentricity> threshold)
    #    }
    #    threshold <- input$min_nuclei_golgi_dist
    #    if ("distance" %in% colnames(results_all_df)){
    #      results_all_df <- subset(results_all_df, results_all_df$distance > threshold)
    #    }

    # feature <- parameters[input$feature_select][[1]][1]
    feature <- parameters[input$feature_comparison][[1]][1]
    condition_col <- input$condition_col

    condition_list <- unlist(unique(results_all_df[condition_col]))
    # plist <- vector('list', length(unique(results_all_df$filename)))
    plist <- vector("list", length(condition_list))
    print("length of plot list")
    print(plist)
    print("list of unique entries")
    print(unlist(unique(results_all_df[condition_col])))

    x_lim <- c(min(results_all_df[feature]), max(results_all_df[feature]))
    # for(file_name in unique(results_all_df$filename)) {
    #  results_df <- subset(results_all_df, results_all_df$filename == file_name )
    for (file_name in condition_list) {
      results_df <- subset(results_all_df, results_all_df[condition_col] == file_name)

      # values <- compute_polarity_index(results_df)


      # print(values)
      # polarity_index <- values[["polarity_index"]]
      # angle_mean_deg <- values[["angle_mean_deg"]]

      x <- unlist(results_df[feature])
      angle_dists[[i]] <- x


      # if (parameters[input$feature_select][[1]][2] == "linear") {
      #
      #    if (x_lim[0] > min(x))
      #        x_lim[0] <- min(x)
      #    if (x_lim[1] < max(x))
      #        x_lim[1] <- max(x)
      #
      #       }

      file_names[[i]] <- file_name

      # polarity_indices[[i]] <- polarity_index
      # angle_mean_degs[[i]] <- angle_mean_deg
      i <- i + 1
    }



    n <- length(angle_dists)
    nCol <- floor(sqrt(n))

    bin_size <- 360 / input$bins

    plotseries <- function(i) {
      angle_dist <- angle_dists[[i]]
      file_name <- file_names[[i]]
      # polarity_index <- polarity_indices[[i]]
      # angle_mean_deg <- angle_mean_degs[[i]]

      # results_df <- subset(results_all_df, results_all_df$filename == file_name)
      results_df <- subset(results_all_df, results_all_df[condition_col] == file_name)

      plot_title <- file_name

      if (nchar(file_name) > 37) {
        max_fl <- 17
        file_name_end <- substr(file_name, nchar(file_name) - max_fl + 1, nchar(file_name))
        file_name_start <- substr(file_name, 1, max_fl)
        plot_title <- paste0(file_name_start, "...", file_name_end)
      }
      # if (nchar(file_name) > 15) {
      #    plot_title <- paste0("image #",toString(i))
      #    print(paste0("filename: ",file_name," too long, will be replaced by",plot_title))
      # }


      if (parameters[input$feature_select][[1]][2] == "directional") {
        statistics <- compute_circular_statistics(results_df, feature, parameters)
        # statistics <- compute_polarity_index(unlist(results_df[feature]))
        x_data <- unlist(results_df[feature]) * 180.0 / pi
        print(paste0("Length of filename", toString(i)))

        p <- rose_plot_circular(parameters, input, statistics, x_data, plot_title, i, text_size)
      } else if (parameters[input$feature_select][[1]][2] == "undirectional") {
        x_data <- results_df[feature]
        # print(x_data)
        statistics <- compute_undirectional_statistics(results_df, feature, parameters)
        # if (input$left_directional) {
        x_data <- unlist(transform_undirectional(input, x_data)) * 180.0 / pi
        # } else {
        #  x_data <- unlist(results_df[feature])*180.0/pi
        # }
        # plot_title <- file_name
        p <- rose_plot_undirectional(parameters, input, statistics, x_data, plot_title, i, text_size)
      } else {
        x_data <- unlist(results_df[feature])
        statistics <- compute_linear_statistics(results_df, feature, parameters)
        # plot_title <- file_name
        # p <- linear_histogram(parameters, input, statistics, x_data,  plot_title, i, text_size, x_lim[0], x_lim[1])
        p <- linear_histogram(parameters, input, statistics, x_data, plot_title, i, text_size, min(results_all_df[feature]), max(results_all_df[feature]))
      }
    }


    myplots <- lapply(1:length(angle_dists), plotseries)

    # print(myplots)
    grid.arrange(grobs = myplots, nrow = nCol) # , widths = list(10,10))
  })


  output$comparison_plot <- renderPlot({

    # source(file = paste0(getwd(),"/src/plot_functions.R"), local=T)
    # source(file = paste0(getwd(),"/src/circular_statistics.R"), local=T)

    # parameters <- fromJSON(file = "parameters/parameters.json")
    # text_size <- as.integer(parameters["text_size_merged_plot"])
    comparison_plot()
  })


  comparisonStatistics <- reactive({
    "
        reactive function that reads a stack of spreadsheet and returns a data frame 
        with descriptive statistics including circular mean, circular standard deviation 
        and nearest neighbours for the merged stack of data
        "

    data <- mergedStack()

    condition_col <- input$condition_col
    control_condition <- input$control_condition
    condition_list <- unlist(unique(data[condition_col]))

    source(file = paste0(getwd(), "/src/plot_functions.R"), local = T)
    source(file = paste0(getwd(), "/src/circular_statistics.R"), local = T)

    parameters <- fromJSON(file = "parameters/parameters.json")


    feature <- parameters[input$feature_comparison][[1]][1]

    res <- data.frame(matrix(ncol = length(condition_list) + 1, nrow = 0))
    cols <- c("Test", "Control")

    for (condition in condition_list) {
      if (condition == control_condition) {
        next
      }
      cols <- c(cols, condition)
    }
    # cols <- append(c("Test"),condition_list)

    colnames(res) <- cols
    res[1, "Test"] <- "WatsonU2"
    res[1, "Control"] <- control_condition

    print("data frame")
    print(res)
    # colnames(res) <- append(c("Test"),condition_list)

    # res[1,"Test"] <- "WatsonU2"


    for (condition in condition_list) {
      if (condition == control_condition) {
        next
      }

      print("Control condition:")
      print(control_condition)
      print("Condition:")
      print(condition)


      condition_data <- subset(data, data[condition_col] == condition)
      control_data <- subset(data, data[condition_col] == control_condition)
      print("Output Watson test object")
      # print(watson.two(condition_data$organelle_orientation_rad, control_data$organelle_orientation_rad, alpha=0.05, plot=TRUE))

      # print("Struct of Watson test object")
      # print(str(watson.two(condition_data$organelle_orientation_rad, control_data$organelle_orientation_rad, alpha=0.05, plot=TRUE)))
      condition_values <- unlist(condition_data[feature])
      control_values <- unlist(control_data[feature])
      # watson1 <- watson.two.test(condition_data$organelle_orientation_rad, control_data$organelle_orientation_rad)
      watson1 <- watson.two.test(condition_values, control_values)
      out <- capture.output(watson.two.test(condition_values, control_values))
      print(out)
      p_value <- out[5]
      res[1, condition] <- p_value
      print("The resulting data frame:")
      print(res)
      # watson1 <- array(as.matrix(unlist(watson1)), dim=c(5, 1))
      # res <- if (watson1[1,]>0.187) (0.04) else (0.06)  # Critical Value: 0.187 #this is just to have number higher and one lower than 0.05 (see coding below)
      # print("Extract")
      # print(watson1)
      # print(res)
    }

    # inFileCondition_1 <- input$condition_1
    # inFileCondition_2 <- input$condition_2

    # if (is.null(inFileCondition_1))
    #    return(NULL)
    # if (is.null(inFileCondition_2))
    #    return(NULL)

    # print(inFileCondition_1$datapath)
    # print(inFileCondition_2$datapath)
    # cond1_data <- read.csv(inFileCondition_1$datapath, header = input$header_cond1)
    # cond2_data <- read.csv(inFileCondition_2$datapath, header = input$header_cond2)

    # print(watson.two(cond1_data$organelle_orientation_rad, cond2_data$organelle_orientation_rad, alpha=0.05, plot=TRUE))

    # watson.two(data1, data2)
    # watson.two(cond1_data$organelle_orientation_rad, cond2_data$organelle_orientation_rad)

    # print("Structure")
    # print(str(watson.two(cond1_data$organelle_orientation_rad, cond2_data$organelle_orientation_rad)))

    # statistics_df <- data.frame()
    # statistics_df
    res
  })

  output$CDFPlot <- renderPlot({
    "
      function that shows the descriptive statistics in table format
      "
    inFileCondition_1 <- input$condition_1
    inFileCondition_2 <- input$condition_2

    if (is.null(inFileCondition_1)) {
      return(NULL)
    }
    if (is.null(inFileCondition_2)) {
      return(NULL)
    }

    print(inFileCondition_1$datapath)
    print(inFileCondition_2$datapath)
    cond1_data <- read.csv(inFileCondition_1$datapath, header = input$header_cond1)
    cond2_data <- read.csv(inFileCondition_2$datapath, header = input$header_cond2)

    watson.two(cond1_data$organelle_orientation_rad, cond2_data$organelle_orientation_rad, alpha = 0.05, plot = TRUE)
  })


  output$comparison_statistics <- renderTable({
    "
    function that shows the descriptive statistics in table format
    "

    statistics_df <- comparisonStatistics()
    statistics_df
  })
}

# Run the application
shinyApp(ui = ui, server = server, options = list(port = opt$p, host = "127.0.0.1"))
