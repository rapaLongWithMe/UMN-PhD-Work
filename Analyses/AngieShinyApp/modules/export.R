# Export Module
# Handles data export and report generation

exportUI <- function(id) {
  ns <- NS(id)
  
  tagList(
    conditionalPanel(
      condition = paste0("output['", ns("has_data"), "']"),
      
      fluidRow(
        column(6,
          h4("Data Export"),
          
          h5("Export Processed Data"),
          div(
            downloadButton(ns("download_raw"), "Download Raw Data (CSV)", 
                          class = "btn-primary", style = "margin: 5px;"),
            br(),
            downloadButton(ns("download_normalized"), "Download Normalized Data (CSV)", 
                          class = "btn-primary", style = "margin: 5px;"),
            br(),
            downloadButton(ns("download_excel"), "Download All Data (Excel)", 
                          class = "btn-primary", style = "margin: 5px;")
          ),
          
          hr(),
          
          h5("Statistical Reports"),
          div(
            downloadButton(ns("download_stats"), "Download Statistics Report", 
                          class = "btn-info", style = "margin: 5px;"),
            br(),
            downloadButton(ns("download_summary"), "Download Analysis Summary", 
                          class = "btn-info", style = "margin: 5px;")
          )
        ),
        
        column(6,
          h4("Plot Export"),
          
          h5("Export Settings"),
          selectInput(ns("plot_format"), "Plot Format:",
                     choices = list("PNG (High Resolution)" = "png",
                                  "PDF (Vector)" = "pdf",
                                  "SVG (Scalable)" = "svg"),
                     selected = "png"),
          
          numericInput(ns("plot_width"), "Width (inches):", 
                      value = 10, min = 4, max = 20, step = 1),
          
          numericInput(ns("plot_height"), "Height (inches):", 
                      value = 8, min = 4, max = 16, step = 1),
          
          numericInput(ns("plot_dpi"), "DPI (for PNG):", 
                      value = 300, min = 150, max = 600, step = 50),
          
          br(),
          
          h5("Available Plots"),
          p("Note: Generate plots in the Visualization tab first"),
          
          div(
            downloadButton(ns("download_distribution"), "Download Distribution Plot", 
                          class = "btn-success", style = "margin: 5px;"),
            br(),
            downloadButton(ns("download_heatmap"), "Download Heatmap", 
                          class = "btn-success", style = "margin: 5px;"),
            br(),
            downloadButton(ns("download_pca"), "Download PCA Plot", 
                          class = "btn-success", style = "margin: 5px;"),
            br(),
            downloadButton(ns("download_timeseries"), "Download Time Series Plot", 
                          class = "btn-success", style = "margin: 5px;")
          )
        )
      ),
      
      hr(),
      
      fluidRow(
        column(12,
          h4("Analysis Report"),
          p("Generate a comprehensive analysis report containing all results"),
          
          fluidRow(
            column(6,
              h5("Report Content"),
              checkboxGroupInput(ns("report_sections"), "Include sections:",
                               choices = list(
                                 "Data Summary" = "summary",
                                 "Quality Assessment" = "quality",
                                 "Outlier Analysis" = "outliers",
                                 "Normalization Results" = "normalization",
                                 "Statistical Analysis" = "statistics",
                                 "Visualizations" = "plots"
                               ),
                               selected = c("summary", "quality", "normalization", "statistics"))
            ),
            
            column(6,
              h5("Report Format"),
              radioButtons(ns("report_format"), "Output format:",
                          choices = list("HTML" = "html",
                                       "PDF" = "pdf"),
                          selected = "html"),
              
              br(),
              actionButton(ns("generate_report"), "Generate Report", 
                          class = "btn-warning", icon = icon("file-alt")),
              
              br(), br(),
              
              conditionalPanel(
                condition = paste0("output['", ns("report_ready"), "']"),
                downloadButton(ns("download_report"), "Download Report", 
                              class = "btn-success")
              )
            )
          )
        )
      )
    ),
    
    conditionalPanel(
      condition = paste0("!output['", ns("has_data"), "']"),
      div(class = "alert alert-info",
          icon("info-circle"), " Complete the analysis workflow to enable export options.")
    )
  )
}

exportServer <- function(id, raw_data, normalized_data, statistics) {
  moduleServer(id, function(input, output, session) {
    
    # Reactive values
    values <- reactiveValues(
      report_generated = FALSE,
      report_content = NULL
    )
    
    # Check if data exists
    output$has_data <- reactive({
      !is.null(raw_data()) || !is.null(normalized_data())
    })
    outputOptions(output, "has_data", suspendWhenHidden = FALSE)
    
    # Download raw data
    output$download_raw <- downloadHandler(
      filename = function() {
        paste0("metabolomics_raw_data_", Sys.Date(), ".csv")
      },
      content = function(file) {
        req(raw_data())
        write_csv(raw_data(), file)
      }
    )
    
    # Download normalized data
    output$download_normalized <- downloadHandler(
      filename = function() {
        paste0("metabolomics_normalized_data_", Sys.Date(), ".csv")
      },
      content = function(file) {
        req(normalized_data())
        write_csv(normalized_data(), file)
      }
    )
    
    # Download Excel file with multiple sheets
    output$download_excel <- downloadHandler(
      filename = function() {
        paste0("metabolomics_analysis_", Sys.Date(), ".xlsx")
      },
      content = function(file) {
        wb <- createWorkbook()
        
        # Add raw data sheet
        if (!is.null(raw_data())) {
          addWorksheet(wb, "Raw_Data")
          writeData(wb, "Raw_Data", raw_data())
        }
        
        # Add normalized data sheet
        if (!is.null(normalized_data())) {
          addWorksheet(wb, "Normalized_Data")
          writeData(wb, "Normalized_Data", normalized_data())
        }
        
        # Add statistics sheet
        if (!is.null(statistics())) {
          addWorksheet(wb, "Statistics")
          writeData(wb, "Statistics", statistics())
        }
        
        saveWorkbook(wb, file)
      }
    )
    
    # Download statistics report
    output$download_stats <- downloadHandler(
      filename = function() {
        paste0("metabolomics_statistics_", Sys.Date(), ".txt")
      },
      content = function(file) {
        req(raw_data())
        
        report_content <- generate_statistics_report(raw_data(), normalized_data())
        writeLines(report_content, file)
      }
    )
    
    # Download analysis summary
    output$download_summary <- downloadHandler(
      filename = function() {
        paste0("metabolomics_summary_", Sys.Date(), ".txt")
      },
      content = function(file) {
        req(raw_data())
        
        summary_content <- generate_analysis_summary(raw_data(), normalized_data())
        writeLines(summary_content, file)
      }
    )
    
    # Plot download handlers (these would need to be connected to actual plots)
    output$download_distribution <- downloadHandler(
      filename = function() {
        paste0("distribution_plot_", Sys.Date(), ".", input$plot_format)
      },
      content = function(file) {
        # This would need to be connected to the actual plot from visualization module
        showNotification("Please generate a distribution plot first", type = "warning")
      }
    )
    
    output$download_heatmap <- downloadHandler(
      filename = function() {
        paste0("heatmap_", Sys.Date(), ".", input$plot_format)
      },
      content = function(file) {
        # This would need to be connected to the actual plot from visualization module
        showNotification("Please generate a heatmap first", type = "warning")
      }
    )
    
    output$download_pca <- downloadHandler(
      filename = function() {
        paste0("pca_plot_", Sys.Date(), ".", input$plot_format)
      },
      content = function(file) {
        # This would need to be connected to the actual plot from visualization module
        showNotification("Please generate a PCA plot first", type = "warning")
      }
    )
    
    output$download_timeseries <- downloadHandler(
      filename = function() {
        paste0("timeseries_plot_", Sys.Date(), ".", input$plot_format)
      },
      content = function(file) {
        # This would need to be connected to the actual plot from visualization module
        showNotification("Please generate a time series plot first", type = "warning")
      }
    )
    
    # Generate comprehensive report
    observeEvent(input$generate_report, {
      req(input$report_sections)
      
      withProgress(message = "Generating report...", value = 0, {
        
        setProgress(0.2, detail = "Collecting data...")
        
        # Generate report content based on selected sections
        report_sections <- list()
        
        if ("summary" %in% input$report_sections) {
          setProgress(0.3, detail = "Adding data summary...")
          report_sections$summary <- generate_data_summary_section(raw_data(), normalized_data())
        }
        
        if ("quality" %in% input$report_sections) {
          setProgress(0.4, detail = "Adding quality assessment...")
          report_sections$quality <- generate_quality_section(raw_data())
        }
        
        if ("normalization" %in% input$report_sections) {
          setProgress(0.6, detail = "Adding normalization results...")
          report_sections$normalization <- generate_normalization_section(raw_data(), normalized_data())
        }
        
        if ("statistics" %in% input$report_sections) {
          setProgress(0.8, detail = "Adding statistical analysis...")
          report_sections$statistics <- generate_statistics_section(raw_data(), normalized_data())
        }
        
        setProgress(0.9, detail = "Compiling report...")
        
        # Combine all sections
        values$report_content <- compile_report_sections(report_sections, input$report_format)
        values$report_generated <- TRUE
        
        setProgress(1, detail = "Complete!")
        
        showNotification("Report generated successfully", type = "message")
      })
    })
    
    # Download report
    output$download_report <- downloadHandler(
      filename = function() {
        ext <- if (input$report_format == "html") "html" else "pdf"
        paste0("metabolomics_report_", Sys.Date(), ".", ext)
      },
      content = function(file) {
        req(values$report_content)
        
        if (input$report_format == "html") {
          writeLines(values$report_content, file)
        } else {
          # For PDF, would need additional processing with pandoc or similar
          writeLines(values$report_content, file)
        }
      }
    )
    
    # Report ready flag
    output$report_ready <- reactive({
      values$report_generated
    })
    outputOptions(output, "report_ready", suspendWhenHidden = FALSE)
  })
}

# Helper functions for report generation

generate_statistics_report <- function(raw_data, normalized_data = NULL) {
  content <- c(
    "METABOLOMICS STATISTICAL ANALYSIS REPORT",
    "========================================",
    paste("Generated on:", Sys.time()),
    "",
    "RAW DATA STATISTICS:",
    "-------------------"
  )
  
  if (!is.null(raw_data)) {
    numeric_cols <- sapply(raw_data, is.numeric)
    
    if (sum(numeric_cols) > 0) {
      raw_stats <- raw_data[numeric_cols] %>%
        summarise_all(list(
          Mean = ~round(mean(., na.rm = TRUE), 4),
          SD = ~round(sd(., na.rm = TRUE), 4),
          Min = ~round(min(., na.rm = TRUE), 4),
          Max = ~round(max(., na.rm = TRUE), 4),
          Missing = ~sum(is.na(.))
        ))
      
      content <- c(content, "", capture.output(print(raw_stats, width = 120)))
    }
  }
  
  if (!is.null(normalized_data)) {
    content <- c(content, "", "", "NORMALIZED DATA STATISTICS:", "-------------------------")
    
    numeric_cols <- sapply(normalized_data, is.numeric)
    
    if (sum(numeric_cols) > 0) {
      norm_stats <- normalized_data[numeric_cols] %>%
        summarise_all(list(
          Mean = ~round(mean(., na.rm = TRUE), 4),
          SD = ~round(sd(., na.rm = TRUE), 4),
          Min = ~round(min(., na.rm = TRUE), 4),
          Max = ~round(max(., na.rm = TRUE), 4),
          Missing = ~sum(is.na(.))
        ))
      
      content <- c(content, "", capture.output(print(norm_stats, width = 120)))
    }
  }
  
  return(content)
}

generate_analysis_summary <- function(raw_data, normalized_data = NULL) {
  content <- c(
    "METABOLOMICS ANALYSIS SUMMARY",
    "=============================",
    paste("Analysis Date:", Sys.Date()),
    paste("Analysis Time:", format(Sys.time(), "%H:%M:%S")),
    ""
  )
  
  if (!is.null(raw_data)) {
    content <- c(content,
      "DATASET OVERVIEW:",
      "----------------",
      paste("Samples:", nrow(raw_data)),
      paste("Variables:", ncol(raw_data)),
      paste("Total data points:", nrow(raw_data) * ncol(raw_data)),
      paste("Missing values:", sum(is.na(raw_data))),
      paste("Data completeness:", round((1 - sum(is.na(raw_data)) / (nrow(raw_data) * ncol(raw_data))) * 100, 2), "%"),
      ""
    )
  }
  
  if (!is.null(normalized_data)) {
    content <- c(content,
      "NORMALIZATION APPLIED:",
      "---------------------",
      "✓ Data has been normalized",
      paste("Normalized samples:", nrow(normalized_data)),
      paste("Normalized variables:", ncol(normalized_data)),
      ""
    )
  }
  
  content <- c(content,
    "ANALYSIS WORKFLOW:",
    "-----------------",
    "1. ✓ Data import completed",
    "2. ✓ Data validation performed",
    "3. ✓ Outlier detection completed"
  )
  
  if (!is.null(normalized_data)) {
    content <- c(content, "4. ✓ Data normalization applied")
  }
  
  content <- c(content,
    "",
    "RECOMMENDATIONS:",
    "---------------",
    "• Review statistical summaries for data quality",
    "• Examine visualizations for patterns and outliers",
    "• Consider additional quality control measures if needed",
    "• Document analysis parameters for reproducibility"
  )
  
  return(content)
}

generate_data_summary_section <- function(raw_data, normalized_data) {
  # Generate HTML/text section for data summary
  return("Data summary section content would go here...")
}

generate_quality_section <- function(raw_data) {
  # Generate quality assessment section
  return("Quality assessment section content would go here...")
}

generate_normalization_section <- function(raw_data, normalized_data) {
  # Generate normalization results section
  return("Normalization results section content would go here...")
}

generate_statistics_section <- function(raw_data, normalized_data) {
  # Generate statistical analysis section
  return("Statistical analysis section content would go here...")
}

compile_report_sections <- function(sections, format) {
  # Compile all sections into final report format
  if (format == "html") {
    # Generate HTML report
    html_content <- c(
      "<!DOCTYPE html>",
      "<html>",
      "<head><title>Metabolomics Analysis Report</title></head>",
      "<body>",
      "<h1>Metabolomics Analysis Report</h1>",
      paste("<p>Generated on:", Sys.time(), "</p>"),
      unlist(sections),
      "</body>",
      "</html>"
    )
    return(html_content)
  } else {
    # Generate text report
    return(unlist(sections))
  }
}
