# Metabolite Abundance Across Time Shiny App
# Main application file

# Load required libraries
library(shiny)
library(shinydashboard)
library(shinyWidgets)
library(DT)
library(plotly)
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(readxl)
library(heatmaply)
library(corrplot)
library(FactoMineR)
library(factoextra)

# Source modules
source("modules/data_import.R")
source("modules/data_validation.R")
source("modules/outlier_detection.R")
source("modules/normalization.R")
source("modules/visualization.R")
source("modules/export.R")
source("utils/helper_functions.R")

# Define UI
ui <- dashboardPage(
  dashboardHeader(title = "Metabolomics Time Series Analysis"),
  
  dashboardSidebar(
    sidebarMenu(
      menuItem("Data Import", tabName = "import", icon = icon("upload")),
      menuItem("Data Quality", tabName = "quality", icon = icon("check-circle")),
      menuItem("Outlier Detection", tabName = "outliers", icon = icon("exclamation-triangle")),
      menuItem("Normalization", tabName = "normalization", icon = icon("balance-scale")),
      menuItem("Visualizations", tabName = "visualizations", icon = icon("chart-line")),
      menuItem("Comparison", tabName = "comparison", icon = icon("columns")),
      menuItem("Export", tabName = "export", icon = icon("download"))
    )
  ),
  
  dashboardBody(
    tags$head(
      tags$link(rel = "stylesheet", type = "text/css", href = "custom.css")
    ),
    
    tabItems(
      # Data Import Tab
      tabItem(tabName = "import",
        fluidRow(
          box(
            title = "Data Import", status = "primary", solidHeader = TRUE, width = 12,
            dataImportUI("data_import")
          )
        )
      ),
      
      # Data Quality Tab
      tabItem(tabName = "quality",
        fluidRow(
          box(
            title = "Data Quality Assessment", status = "primary", solidHeader = TRUE, width = 12,
            dataValidationUI("data_validation")
          )
        )
      ),
      
      # Outlier Detection Tab
      tabItem(tabName = "outliers",
        fluidRow(
          box(
            title = "Outlier Detection", status = "primary", solidHeader = TRUE, width = 12,
            outlierDetectionUI("outlier_detection")
          )
        )
      ),
      
      # Normalization Tab
      tabItem(tabName = "normalization",
        fluidRow(
          box(
            title = "Data Normalization", status = "primary", solidHeader = TRUE, width = 12,
            normalizationUI("normalization")
          )
        )
      ),
      
      # Visualizations Tab
      tabItem(tabName = "visualizations",
        fluidRow(
          box(
            title = "Data Visualizations", status = "primary", solidHeader = TRUE, width = 12,
            visualizationUI("visualization")
          )
        )
      ),
      
      # Comparison Tab
      tabItem(tabName = "comparison",
        fluidRow(
          box(
            title = "Raw vs Normalized Comparison", status = "primary", solidHeader = TRUE, width = 12,
            h3("Comparative Analysis"),
            p("Compare raw and normalized data side by side"),
            # Comparison interface will be implemented
            plotlyOutput("comparison_plot", height = "600px")
          )
        )
      ),
      
      # Export Tab
      tabItem(tabName = "export",
        fluidRow(
          box(
            title = "Export Results", status = "primary", solidHeader = TRUE, width = 12,
            exportUI("export")
          )
        )
      )
    )
  )
)

# Define server logic
server <- function(input, output, session) {
  # Reactive values to store data throughout the app
  values <- reactiveValues(
    raw_data = NULL,
    cleaned_data = NULL,
    normalized_data = NULL,
    outliers = NULL,
    statistics = NULL
  )
  
  # Call server modules
  imported_data <- dataImportServer("data_import")
  
  # Update values when data is imported
  observe({
    if (!is.null(imported_data())) {
      values$raw_data <- imported_data()
    }
  })
  
  # Data validation server
  dataValidationServer("data_validation", reactive(values$raw_data))
  
  # Outlier detection server
  outlier_results <- outlierDetectionServer("outlier_detection", reactive(values$raw_data))
  
  # Update cleaned data after outlier removal
  observe({
    if (!is.null(outlier_results())) {
      values$cleaned_data <- outlier_results()$cleaned_data
      values$outliers <- outlier_results()$outliers
    }
  })
  
  # Normalization server
  normalized_results <- normalizationServer("normalization", reactive(values$cleaned_data))
  
  # Update normalized data
  observe({
    if (!is.null(normalized_results())) {
      values$normalized_data <- normalized_results()
    }
  })
  
  # Visualization server
  visualizationServer("visualization", 
                     raw_data = reactive(values$cleaned_data),
                     normalized_data = reactive(values$normalized_data))
  
  # Export server
  exportServer("export", 
               raw_data = reactive(values$cleaned_data),
               normalized_data = reactive(values$normalized_data),
               statistics = reactive(values$statistics))
  
  # Comparison plot
  output$comparison_plot <- renderPlotly({
    req(values$cleaned_data, values$normalized_data)
    
    # Create comparison visualization
    create_comparison_plot(values$cleaned_data, values$normalized_data)
  })
}

# Run the application
shinyApp(ui = ui, server = server)
