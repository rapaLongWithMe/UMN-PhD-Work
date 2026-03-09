# Data Import Module
# Handles file upload, validation, and preview functionality

dataImportUI <- function(id) {
  ns <- NS(id)
  
  tagList(
    fluidRow(
      column(6,
        h4("Upload mzMine Data"),
        fileInput(ns("file"), "Choose CSV or Excel File",
                 multiple = FALSE,
                 accept = c(".csv", ".xlsx", ".xls")),
        
        conditionalPanel(
          condition = paste0("output['", ns("file_uploaded"), "']"),
          h5("File Information"),
          verbatimTextOutput(ns("file_info"))
        )
      ),
      
      column(6,
        h4("Data Format Requirements"),
        div(
          p(strong("Required Columns:")),
          tags$ul(
            tags$li("Metabolite ID/Name"),
            tags$li("Sample ID"),
            tags$li("Time point"),
            tags$li("Abundance values"),
            tags$li("Replicate information")
          ),
          p(strong("Optional Columns:")),
          tags$ul(
            tags$li("Sample groups"),
            tags$li("Experimental conditions"),
            tags$li("Additional metadata")
          )
        )
      )
    ),
    
    hr(),
    
    conditionalPanel(
      condition = paste0("output['", ns("data_loaded"), "']"),
      fluidRow(
        column(12,
          h4("Data Preview"),
          DTOutput(ns("data_preview")),
          br(),
          actionButton(ns("confirm_data"), "Confirm Data Import", 
                      class = "btn-success", icon = icon("check"))
        )
      )
    ),
    
    conditionalPanel(
      condition = paste0("output['", ns("data_confirmed"), "']"),
      fluidRow(
        column(12,
          div(class = "alert alert-success",
              icon("check-circle"), " Data successfully imported!")
        )
      )
    )
  )
}

dataImportServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    # Reactive values
    values <- reactiveValues(
      raw_data = NULL,
      confirmed = FALSE
    )
    
    # File upload handling
    observeEvent(input$file, {
      req(input$file)
      
      ext <- tools::file_ext(input$file$datapath)
      
      # Read file based on extension
      if (ext == "csv") {
        values$raw_data <- read_csv(input$file$datapath, show_col_types = FALSE)
      } else if (ext %in% c("xlsx", "xls")) {
        values$raw_data <- read_excel(input$file$datapath)
      } else {
        showNotification("Unsupported file format", type = "warning")
        return()
      }
      
      # Basic validation
      if (nrow(values$raw_data) == 0) {
        showNotification("File appears to be empty", type = "warning")
        values$raw_data <- NULL
        return()
      }
      
      showNotification("File uploaded successfully", type = "default")
    })
    
    # File information output
    output$file_info <- renderText({
      req(input$file, values$raw_data)
      
      paste(
        "Filename:", input$file$name,
        "\nSize:", format(input$file$size, units = "Mb", digits = 2),
        "\nRows:", nrow(values$raw_data),
        "\nColumns:", ncol(values$raw_data)
      )
    })
    
    # Data preview
    output$data_preview <- renderDT({
      req(values$raw_data)
      
      datatable(
        head(values$raw_data, 100),
        options = list(
          scrollX = TRUE,
          pageLength = 10,
          dom = 'tip'
        ),
        class = 'cell-border stripe'
      )
    })
    
    # Confirm data import
    observeEvent(input$confirm_data, {
      req(values$raw_data)
      values$confirmed <- TRUE
      showNotification("Data import confirmed", type = "message")
    })
    
    # Output flags for conditional panels
    output$file_uploaded <- reactive({
      !is.null(input$file)
    })
    outputOptions(output, "file_uploaded", suspendWhenHidden = FALSE)
    
    output$data_loaded <- reactive({
      !is.null(values$raw_data)
    })
    outputOptions(output, "data_loaded", suspendWhenHidden = FALSE)
    
    output$data_confirmed <- reactive({
      values$confirmed
    })
    outputOptions(output, "data_confirmed", suspendWhenHidden = FALSE)
    
    # Return the imported data
    reactive({
      if (values$confirmed) {
        values$raw_data
      } else {
        NULL
      }
    })
  })
}
