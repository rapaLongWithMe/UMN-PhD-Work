# Outlier Detection Module
# Handles identification and management of outlier replicates

outlierDetectionUI <- function(id) {
  ns <- NS(id)
  
  tagList(
    conditionalPanel(
      condition = paste0("output['", ns("has_data"), "']"),
      fluidRow(
        column(4,
          h4("Outlier Detection Methods"),
          
          selectInput(ns("outlier_method"), "Detection Method:",
                     choices = list(
                       "Z-Score" = "zscore",
                       "IQR Method" = "iqr", 
                       "Modified Z-Score (MAD)" = "mad"
                     ),
                     selected = "zscore"),
          
          conditionalPanel(
            condition = paste0("input['", ns("outlier_method"), "'] == 'zscore'"),
            numericInput(ns("zscore_threshold"), "Z-Score Threshold:", 
                        value = 3, min = 1, max = 5, step = 0.1)
          ),
          
          conditionalPanel(
            condition = paste0("input['", ns("outlier_method"), "'] == 'iqr'"),
            numericInput(ns("iqr_multiplier"), "IQR Multiplier:", 
                        value = 1.5, min = 1, max = 3, step = 0.1)
          ),
          
          conditionalPanel(
            condition = paste0("input['", ns("outlier_method"), "'] == 'mad'"),
            numericInput(ns("mad_threshold"), "MAD Threshold:", 
                        value = 3, min = 1, max = 5, step = 0.1)
          ),
          
          br(),
          actionButton(ns("detect_outliers"), "Detect Outliers", 
                      class = "btn-primary", icon = icon("search")),
          
          br(), br(),
          
          conditionalPanel(
            condition = paste0("output['", ns("outliers_detected"), "']"),
            h5("Outlier Summary"),
            valueBoxOutput(ns("outlier_count"), width = 12),
            valueBoxOutput(ns("outlier_percentage"), width = 12)
          )
        ),
        
        column(8,
          conditionalPanel(
            condition = paste0("output['", ns("outliers_detected"), "']"),
            h4("Outlier Visualization"),
            plotlyOutput(ns("outlier_plot"), height = "400px"),
            
            br(),
            
            h4("Detected Outliers"),
            DTOutput(ns("outlier_table")),
            
            br(),
            
            fluidRow(
              column(6,
                actionButton(ns("remove_outliers"), "Remove Selected Outliers", 
                           class = "btn-warning", icon = icon("trash"))
              ),
              column(6,
                actionButton(ns("keep_all"), "Keep All Data", 
                           class = "btn-success", icon = icon("check"))
              )
            )
          )
        )
      ),
      
      hr(),
      
      conditionalPanel(
        condition = paste0("output['", ns("data_processed"), "']"),
        fluidRow(
          column(12,
            div(class = "alert alert-success",
                icon("check-circle"), " Outlier detection completed!"),
            
            h4("Final Dataset Summary"),
            verbatimTextOutput(ns("final_summary"))
          )
        )
      )
    ),
    
    conditionalPanel(
      condition = paste0("!output['", ns("has_data"), "']"),
      div(class = "alert alert-info",
          icon("info-circle"), " Please import data first to detect outliers.")
    )
  )
}

outlierDetectionServer <- function(id, data) {
  moduleServer(id, function(input, output, session) {
    
    # Reactive values
    values <- reactiveValues(
      outliers = NULL,
      cleaned_data = NULL,
      processing_complete = FALSE
    )
    
    # Check if data exists
    output$has_data <- reactive({
      !is.null(data()) && nrow(data()) > 0
    })
    outputOptions(output, "has_data", suspendWhenHidden = FALSE)
    
    # Detect outliers
    observeEvent(input$detect_outliers, {
      req(data())
      
      withProgress(message = "Detecting outliers...", value = 0, {
        
        # Identify numeric columns (abundance data)
        numeric_cols <- sapply(data(), is.numeric)
        
        if (sum(numeric_cols) == 0) {
          showNotification("No numeric columns found for outlier detection", type = "warning")
          return()
        }
        
        setProgress(0.3, detail = "Analyzing data distribution...")
        
        # Apply selected outlier detection method
        outlier_indices <- c()
        
        if (input$outlier_method == "zscore") {
          outlier_indices <- detect_outliers_zscore(data(), input$zscore_threshold)
        } else if (input$outlier_method == "iqr") {
          outlier_indices <- detect_outliers_iqr(data(), input$iqr_multiplier)
        } else if (input$outlier_method == "mad") {
          outlier_indices <- detect_outliers_mad(data(), input$mad_threshold)
        }
        
        setProgress(0.7, detail = "Compiling results...")
        
        if (length(outlier_indices) > 0) {
          values$outliers <- data()[outlier_indices, ]
          values$outliers$outlier_row <- outlier_indices
        } else {
          values$outliers <- data.frame()
        }
        
        setProgress(1, detail = "Complete!")
        
        showNotification(
          paste("Detected", nrow(values$outliers), "potential outliers"),
          type = "message"
        )
      })
    })
    
    # Outlier count value box
    output$outlier_count <- renderValueBox({
      outlier_count <- if (is.null(values$outliers)) 0 else nrow(values$outliers)
      valueBox(
        value = outlier_count,
        subtitle = "Detected Outliers",
        icon = icon("exclamation-triangle"),
        color = if (outlier_count > 0) "yellow" else "green"
      )
    })
    
    # Outlier percentage value box
    output$outlier_percentage <- renderValueBox({
      if (is.null(data()) || is.null(values$outliers)) {
        percentage <- 0
      } else {
        percentage <- round((nrow(values$outliers) / nrow(data())) * 100, 1)
      }
      
      valueBox(
        value = paste0(percentage, "%"),
        subtitle = "Percentage of Data",
        icon = icon("percent"),
        color = if (percentage > 10) "red" else if (percentage > 5) "yellow" else "green"
      )
    })
    
    # Outlier visualization
    output$outlier_plot <- renderPlotly({
      req(values$outliers)
      
      if (nrow(values$outliers) == 0) {
        # No outliers found
        p <- ggplot() + 
          annotate("text", x = 0.5, y = 0.5, label = "No outliers detected", 
                   size = 6, color = "darkgreen") +
          theme_void() +
          xlim(0, 1) + ylim(0, 1)
        
        return(ggplotly(p))
      }
      
      # Create visualization showing outliers vs normal data
      numeric_cols <- sapply(data(), is.numeric)
      if (sum(numeric_cols) > 0) {
        first_numeric_col <- names(data())[numeric_cols][1]
        
        plot_data <- data()
        plot_data$is_outlier <- 1:nrow(plot_data) %in% values$outliers$outlier_row
        plot_data$row_number <- 1:nrow(plot_data)
        
        p <- ggplot(plot_data, aes(x = row_number, y = .data[[first_numeric_col]], 
                                   color = is_outlier)) +
          geom_point(alpha = 0.7, size = 2) +
          scale_color_manual(values = c("FALSE" = "steelblue", "TRUE" = "red"),
                           labels = c("Normal", "Outlier"),
                           name = "Data Type") +
          labs(
            title = paste("Outlier Detection Results -", first_numeric_col),
            x = "Sample Index",
            y = first_numeric_col
          ) +
          theme_minimal() +
          theme(
            plot.title = element_text(size = 14, face = "bold"),
            legend.position = "top"
          )
        
        ggplotly(p, tooltip = c("x", "y", "colour"))
      }
    })
    
    # Outlier table
    output$outlier_table <- renderDT({
      req(values$outliers)
      
      if (nrow(values$outliers) == 0) {
        return(datatable(data.frame(Message = "No outliers detected")))
      }
      
      display_outliers <- values$outliers %>%
        select(-outlier_row) %>%
        slice_head(n = 100)  # Limit display for performance
      
      datatable(
        display_outliers,
        options = list(
          pageLength = 10,
          scrollX = TRUE,
          dom = 'tip'
        ),
        selection = 'multiple',
        rownames = FALSE
      ) %>%
        formatStyle(columns = 1:ncol(display_outliers), 
                   backgroundColor = "#ffeeee")
    })
    
    # Remove outliers
    observeEvent(input$remove_outliers, {
      req(data(), values$outliers)
      
      if (nrow(values$outliers) > 0) {
        outlier_rows <- values$outliers$outlier_row
        values$cleaned_data <- data()[-outlier_rows, ]
        values$processing_complete <- TRUE
        
        showNotification(
          paste("Removed", length(outlier_rows), "outliers from dataset"),
          type = "message"
        )
      }
    })
    
    # Keep all data
    observeEvent(input$keep_all, {
      req(data())
      
      values$cleaned_data <- data()
      values$processing_complete <- TRUE
      
      showNotification("Proceeding with complete dataset", type = "message")
    })
    
    # Final summary
    output$final_summary <- renderText({
      req(values$cleaned_data)
      
      original_rows <- nrow(data())
      final_rows <- nrow(values$cleaned_data)
      removed_rows <- original_rows - final_rows
      
      paste0(
        "OUTLIER PROCESSING SUMMARY\n",
        "==========================\n",
        "Original dataset: ", original_rows, " rows\n",
        "Outliers detected: ", ifelse(is.null(values$outliers), 0, nrow(values$outliers)), "\n",
        "Outliers removed: ", removed_rows, "\n",
        "Final dataset: ", final_rows, " rows\n",
        "Data retention: ", round((final_rows / original_rows) * 100, 1), "%"
      )
    })
    
    # Output flags for conditional panels
    output$outliers_detected <- reactive({
      !is.null(values$outliers)
    })
    outputOptions(output, "outliers_detected", suspendWhenHidden = FALSE)
    
    output$data_processed <- reactive({
      values$processing_complete
    })
    outputOptions(output, "data_processed", suspendWhenHidden = FALSE)
    
    # Return processed data
    reactive({
      if (values$processing_complete) {
        list(
          cleaned_data = values$cleaned_data,
          outliers = values$outliers
        )
      } else {
        NULL
      }
    })
  })
}

# Helper functions for outlier detection
detect_outliers_zscore <- function(data, threshold = 3) {
  numeric_cols <- sapply(data, is.numeric)
  
  if (sum(numeric_cols) == 0) return(integer(0))
  
  outlier_rows <- c()
  
  for (col in names(data)[numeric_cols]) {
    values <- data[[col]]
    z_scores <- abs((values - mean(values, na.rm = TRUE)) / sd(values, na.rm = TRUE))
    outliers <- which(z_scores > threshold & !is.na(z_scores))
    outlier_rows <- c(outlier_rows, outliers)
  }
  
  unique(outlier_rows)
}

detect_outliers_iqr <- function(data, multiplier = 1.5) {
  numeric_cols <- sapply(data, is.numeric)
  
  if (sum(numeric_cols) == 0) return(integer(0))
  
  outlier_rows <- c()
  
  for (col in names(data)[numeric_cols]) {
    values <- data[[col]]
    Q1 <- quantile(values, 0.25, na.rm = TRUE)
    Q3 <- quantile(values, 0.75, na.rm = TRUE)
    IQR <- Q3 - Q1
    
    lower_bound <- Q1 - multiplier * IQR
    upper_bound <- Q3 + multiplier * IQR
    
    outliers <- which((values < lower_bound | values > upper_bound) & !is.na(values))
    outlier_rows <- c(outlier_rows, outliers)
  }
  
  unique(outlier_rows)
}

detect_outliers_mad <- function(data, threshold = 3) {
  numeric_cols <- sapply(data, is.numeric)
  
  if (sum(numeric_cols) == 0) return(integer(0))
  
  outlier_rows <- c()
  
  for (col in names(data)[numeric_cols]) {
    values <- data[[col]]
    median_val <- median(values, na.rm = TRUE)
    mad_val <- mad(values, na.rm = TRUE)
    
    if (mad_val == 0) next  # Skip if MAD is 0
    
    modified_z <- 0.6745 * (values - median_val) / mad_val
    outliers <- which(abs(modified_z) > threshold & !is.na(modified_z))
    outlier_rows <- c(outlier_rows, outliers)
  }
  
  unique(outlier_rows)
}
