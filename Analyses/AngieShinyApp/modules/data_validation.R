# Data Validation Module
# Handles data quality assessment and validation

dataValidationUI <- function(id) {
  ns <- NS(id)
  
  tagList(
    conditionalPanel(
      condition = paste0("output['", ns("has_data"), "']"),
      fluidRow(
        column(6,
          h4("Data Structure"),
          DTOutput(ns("data_structure")),
          br(),
          h4("Column Information"),
          DTOutput(ns("column_info"))
        ),
        
        column(6,
          h4("Data Quality Metrics"),
          valueBoxOutput(ns("total_rows"), width = 12),
          valueBoxOutput(ns("total_columns"), width = 12),
          valueBoxOutput(ns("missing_values"), width = 12),
          valueBoxOutput(ns("completeness"), width = 12)
        )
      ),
      
      hr(),
      
      fluidRow(
        column(12,
          h4("Missing Values Analysis"),
          plotlyOutput(ns("missing_plot"), height = "400px")
        )
      ),
      
      hr(),
      
      fluidRow(
        column(12,
          h4("Data Quality Report"),
          downloadButton(ns("download_report"), "Download Quality Report", 
                        class = "btn-primary"),
          br(), br(),
          verbatimTextOutput(ns("quality_summary"))
        )
      )
    ),
    
    conditionalPanel(
      condition = paste0("!output['", ns("has_data"), "']"),
      div(class = "alert alert-info",
          icon("info-circle"), " Please import data first to view quality assessment.")
    )
  )
}

dataValidationServer <- function(id, data) {
  moduleServer(id, function(input, output, session) {
    
    # Check if data exists
    output$has_data <- reactive({
      !is.null(data()) && nrow(data()) > 0
    })
    outputOptions(output, "has_data", suspendWhenHidden = FALSE)
    
    # Data structure table
    output$data_structure <- renderDT({
      req(data())
      
      structure_df <- data.frame(
        Column = names(data()),
        Type = sapply(data(), class),
        Sample_Value = sapply(names(data()), function(x) {
          val <- data()[[x]][1]
          if (is.na(val)) "NA" else as.character(val)
        }),
        stringsAsFactors = FALSE
      )
      
      datatable(
        structure_df,
        options = list(
          pageLength = 15,
          dom = 't'
        ),
        rownames = FALSE
      )
    })
    
    # Column information
    output$column_info <- renderDT({
      req(data())
      
      # Create column information in a simpler way
      col_names <- names(data())
      col_info <- data.frame(
        Column = col_names,
        Missing = sapply(col_names, function(x) sum(is.na(data()[[x]]))),
        Complete = sapply(col_names, function(x) sum(!is.na(data()[[x]]))),
        Unique = sapply(col_names, function(x) n_distinct(data()[[x]], na.rm = TRUE)),
        stringsAsFactors = FALSE
      )
      
      datatable(
        col_info,
        options = list(
          pageLength = 15,
          dom = 'tip'
        ),
        rownames = FALSE
      )
    })
    
    # Value boxes
    output$total_rows <- renderValueBox({
      valueBox(
        value = ifelse(is.null(data()), 0, nrow(data())),
        subtitle = "Total Rows",
        icon = icon("table"),
        color = "blue"
      )
    })
    
    output$total_columns <- renderValueBox({
      valueBox(
        value = ifelse(is.null(data()), 0, ncol(data())),
        subtitle = "Total Columns",
        icon = icon("columns"),
        color = "green"
      )
    })
    
    output$missing_values <- renderValueBox({
      missing_count <- if (is.null(data())) 0 else sum(is.na(data()))
      valueBox(
        value = missing_count,
        subtitle = "Missing Values",
        icon = icon("exclamation-triangle"),
        color = if (missing_count > 0) "red" else "green"
      )
    })
    
    output$completeness <- renderValueBox({
      if (is.null(data())) {
        completeness <- 0
      } else {
        total_cells <- nrow(data()) * ncol(data())
        missing_cells <- sum(is.na(data()))
        completeness <- round((1 - missing_cells / total_cells) * 100, 1)
      }
      
      valueBox(
        value = paste0(completeness, "%"),
        subtitle = "Data Completeness",
        icon = icon("check-circle"),
        color = if (completeness >= 90) "green" else if (completeness >= 70) "yellow" else "red"
      )
    })
    
    # Missing values plot
    output$missing_plot <- renderPlotly({
      req(data())
      
      missing_data <- data() %>%
        summarise_all(~sum(is.na(.))) %>%
        gather(key = "Column", value = "Missing_Count") %>%
        mutate(Percentage = Missing_Count / nrow(data()) * 100) %>%
        arrange(desc(Missing_Count))
      
      p <- ggplot(missing_data, aes(x = reorder(Column, Missing_Count), y = Missing_Count)) +
        geom_col(fill = "steelblue", alpha = 0.7) +
        coord_flip() +
        labs(
          title = "Missing Values by Column",
          x = "Column",
          y = "Number of Missing Values"
        ) +
        theme_minimal() +
        theme(
          plot.title = element_text(size = 14, face = "bold"),
          axis.text = element_text(size = 10)
        )
      
      ggplotly(p, tooltip = c("x", "y"))
    })
    
    # Quality summary
    output$quality_summary <- renderText({
      req(data())
      
      total_rows <- nrow(data())
      total_cols <- ncol(data())
      total_cells <- total_rows * total_cols
      missing_cells <- sum(is.na(data()))
      completeness <- round((1 - missing_cells / total_cells) * 100, 1)
      
      # Check for required columns (basic metabolomics structure)
      required_patterns <- c("metabolite", "sample", "time", "abundance", "replicate")
      col_names_lower <- tolower(names(data()))
      found_columns <- sapply(required_patterns, function(pattern) {
        any(grepl(pattern, col_names_lower, ignore.case = TRUE))
      })
      
      summary_text <- paste0(
        "DATA QUALITY SUMMARY\n",
        "===================\n\n",
        "Dimensions: ", total_rows, " rows × ", total_cols, " columns\n",
        "Missing values: ", missing_cells, " (", round(missing_cells/total_cells*100, 2), "%)\n",
        "Data completeness: ", completeness, "%\n\n",
        "COLUMN STRUCTURE ASSESSMENT:\n",
        "============================\n"
      )
      
      for (i in seq_along(required_patterns)) {
        status <- if (found_columns[i]) "✓ FOUND" else "✗ MISSING"
        summary_text <- paste0(summary_text, 
                              "Required pattern '", required_patterns[i], "': ", status, "\n")
      }
      
      summary_text <- paste0(summary_text, "\n",
        "RECOMMENDATIONS:\n",
        "===============\n")
      
      if (completeness < 90) {
        summary_text <- paste0(summary_text, 
                              "• Review missing values before proceeding with analysis\n")
      }
      
      if (!all(found_columns)) {
        summary_text <- paste0(summary_text,
                              "• Verify column names match expected metabolomics data format\n")
      }
      
      if (completeness >= 90 && all(found_columns)) {
        summary_text <- paste0(summary_text,
                              "• Data appears to be in good condition for analysis\n")
      }
      
      summary_text
    })
    
    # Download handler for quality report
    output$download_report <- downloadHandler(
      filename = function() {
        paste0("data_quality_report_", Sys.Date(), ".txt")
      },
      content = function(file) {
        req(data())
        
        # Generate comprehensive quality report
        report_content <- paste0(
          "METABOLOMICS DATA QUALITY REPORT\n",
          "Generated on: ", Sys.time(), "\n",
          "================================\n\n",
          "DATASET OVERVIEW:\n",
          "Rows: ", nrow(data()), "\n",
          "Columns: ", ncol(data()), "\n",
          "Total cells: ", nrow(data()) * ncol(data()), "\n",
          "Missing values: ", sum(is.na(data())), "\n",
          "Completeness: ", round((1 - sum(is.na(data())) / (nrow(data()) * ncol(data()))) * 100, 2), "%\n\n",
          "COLUMN DETAILS:\n"
        )
        
        # Add column-specific information
        for (col in names(data())) {
          missing_count <- sum(is.na(data()[[col]]))
          unique_count <- n_distinct(data()[[col]], na.rm = TRUE)
          
          report_content <- paste0(report_content,
            "- ", col, ": ", class(data()[[col]])[1], 
            " (Missing: ", missing_count, 
            ", Unique: ", unique_count, ")\n")
        }
        
        writeLines(report_content, file)
      }
    )
  })
}
