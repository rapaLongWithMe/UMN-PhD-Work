# Normalization Module
# Handles data normalization with multiple strategies

normalizationUI <- function(id) {
  ns <- NS(id)
  
  tagList(
    conditionalPanel(
      condition = paste0("output['", ns("has_data"), "']"),
      fluidRow(
        column(4,
          h4("Normalization Settings"),
          
          selectInput(ns("norm_method"), "Normalization Method:",
                     choices = list(
                       "Z-Score (Standard)" = "zscore",
                       "Min-Max Scaling" = "minmax",
                       "Robust Scaling (Median + MAD)" = "robust",
                       "Log Transformation" = "log",
                       "Log + Z-Score" = "log_zscore"
                     ),
                     selected = "zscore"),
          
          conditionalPanel(
            condition = paste0("input['", ns("norm_method"), "'] == 'log' || input['", ns("norm_method"), "'] == 'log_zscore'"),
            radioButtons(ns("log_base"), "Log Base:",
                        choices = list("Natural (ln)" = "natural",
                                     "Base 2" = "base2",
                                     "Base 10" = "base10"),
                        selected = "natural")
          ),
          
          hr(),
          
          h5("Preview Options"),
          selectInput(ns("preview_metabolite"), "Metabolite to Preview:",
                     choices = NULL),
          
          br(),
          actionButton(ns("preview_norm"), "Preview Normalization", 
                      class = "btn-info", icon = icon("eye")),
          
          br(), br(),
          actionButton(ns("apply_norm"), "Apply Normalization", 
                      class = "btn-success", icon = icon("check"))
        ),
        
        column(8,
          conditionalPanel(
            condition = paste0("output['", ns("preview_available"), "']"),
            h4("Normalization Preview"),
            plotlyOutput(ns("preview_plot"), height = "400px")
          ),
          
          conditionalPanel(
            condition = paste0("output['", ns("normalization_applied"), "']"),
            h4("Normalization Results"),
            
            fluidRow(
              column(6,
                h5("Before Normalization"),
                DTOutput(ns("before_stats"))
              ),
              column(6,
                h5("After Normalization"),
                DTOutput(ns("after_stats"))
              )
            ),
            
            br(),
            
            h5("Distribution Comparison"),
            plotlyOutput(ns("comparison_plot"), height = "400px")
          )
        )
      )
    ),
    
    conditionalPanel(
      condition = paste0("!output['", ns("has_data"), "']"),
      div(class = "alert alert-info",
          icon("info-circle"), " Please complete outlier detection first.")
    )
  )
}

normalizationServer <- function(id, data) {
  moduleServer(id, function(input, output, session) {
    
    # Reactive values
    values <- reactiveValues(
      normalized_data = NULL,
      preview_data = NULL,
      normalization_applied = FALSE
    )
    
    # Check if data exists
    output$has_data <- reactive({
      !is.null(data()) && nrow(data()) > 0
    })
    outputOptions(output, "has_data", suspendWhenHidden = FALSE)
    
    # Update metabolite choices for preview
    observe({
      req(data())
      
      # Identify likely metabolite columns (numeric columns)
      numeric_cols <- names(data())[sapply(data(), is.numeric)]
      
      updateSelectInput(session, "preview_metabolite",
                       choices = numeric_cols,
                       selected = numeric_cols[1])
    })
    
    # Preview normalization
    observeEvent(input$preview_norm, {
      req(data(), input$preview_metabolite, input$norm_method)
      
      withProgress(message = "Generating preview...", value = 0, {
        
        setProgress(0.3, detail = "Calculating normalization...")
        
        # Get the column to preview
        col_data <- data()[[input$preview_metabolite]]
        
        # Apply normalization to preview column
        normalized_col <- apply_normalization_method(col_data, input$norm_method, input$log_base)
        
        setProgress(0.7, detail = "Creating visualization...")
        
        # Store preview data
        values$preview_data <- data.frame(
          Original = col_data,
          Normalized = normalized_col,
          Index = 1:length(col_data)
        )
        
        setProgress(1, detail = "Complete!")
      })
    })
    
    # Apply normalization to full dataset
    observeEvent(input$apply_norm, {
      req(data(), input$norm_method)
      
      withProgress(message = "Applying normalization...", value = 0, {
        
        # Identify numeric columns
        numeric_cols <- sapply(data(), is.numeric)
        
        if (sum(numeric_cols) == 0) {
          showNotification("No numeric columns found for normalization", type = "warning")
          return()
        }
        
        setProgress(0.2, detail = "Preparing data...")
        
        # Start with original data
        normalized_data <- data()
        
        setProgress(0.4, detail = "Normalizing columns...")
        
        # Apply normalization to each numeric column
        for (col in names(data())[numeric_cols]) {
          setProgress(0.4 + 0.4 * which(names(data())[numeric_cols] == col) / sum(numeric_cols),
                     detail = paste("Normalizing", col, "..."))
          
          normalized_data[[col]] <- apply_normalization_method(
            data()[[col]], input$norm_method, input$log_base
          )
        }
        
        setProgress(0.9, detail = "Finalizing...")
        
        values$normalized_data <- normalized_data
        values$normalization_applied <- TRUE
        
        setProgress(1, detail = "Complete!")
        
        showNotification(
          paste("Normalization applied using", input$norm_method, "method"),
          type = "message"
        )
      })
    })
    
    # Preview plot
    output$preview_plot <- renderPlotly({
      req(values$preview_data)
      
      # Create before/after comparison
      plot_data <- values$preview_data %>%
        gather(key = "Type", value = "Value", Original, Normalized) %>%
        filter(!is.na(Value) & is.finite(Value))
      
      p <- ggplot(plot_data, aes(x = Value, fill = Type)) +
        geom_histogram(alpha = 0.7, bins = 30, position = "identity") +
        facet_wrap(~Type, scales = "free") +
        scale_fill_manual(values = c("Original" = "steelblue", "Normalized" = "orange")) +
        labs(
          title = paste("Normalization Preview:", input$preview_metabolite),
          x = "Value",
          y = "Frequency"
        ) +
        theme_minimal() +
        theme(
          plot.title = element_text(size = 14, face = "bold"),
          strip.text = element_text(size = 12, face = "bold")
        )
      
      ggplotly(p)
    })
    
    # Before normalization statistics
    output$before_stats <- renderDT({
      req(data(), values$normalized_data)
      
      numeric_cols <- sapply(data(), is.numeric)
      numeric_data <- data()[numeric_cols]
      
      if (ncol(numeric_data) > 0) {
        col_names <- names(numeric_data)
        stats <- data.frame(
          Metabolite = col_names,
          Mean = round(sapply(numeric_data, function(x) mean(x, na.rm = TRUE)), 4),
          SD = round(sapply(numeric_data, function(x) sd(x, na.rm = TRUE)), 4),
          Min = round(sapply(numeric_data, function(x) min(x, na.rm = TRUE)), 4),
          Max = round(sapply(numeric_data, function(x) max(x, na.rm = TRUE)), 4),
          stringsAsFactors = FALSE
        )
        
        datatable(
          stats,
          options = list(
            pageLength = 10,
            dom = 'tip',
            scrollX = TRUE
          ),
          rownames = FALSE
        ) %>%
          formatRound(columns = 2:5, digits = 4)
      }
    })
    
    # After normalization statistics
    output$after_stats <- renderDT({
      req(values$normalized_data)
      
      numeric_cols <- sapply(values$normalized_data, is.numeric)
      numeric_data <- values$normalized_data[numeric_cols]
      
      if (ncol(numeric_data) > 0) {
        col_names <- names(numeric_data)
        stats <- data.frame(
          Metabolite = col_names,
          Mean = round(sapply(numeric_data, function(x) mean(x, na.rm = TRUE)), 4),
          SD = round(sapply(numeric_data, function(x) sd(x, na.rm = TRUE)), 4),
          Min = round(sapply(numeric_data, function(x) min(x, na.rm = TRUE)), 4),
          Max = round(sapply(numeric_data, function(x) max(x, na.rm = TRUE)), 4),
          stringsAsFactors = FALSE
        )
        
        datatable(
          stats,
          options = list(
            pageLength = 10,
            dom = 'tip',
            scrollX = TRUE
          ),
          rownames = FALSE
        ) %>%
          formatRound(columns = 2:5, digits = 4)
      }
    })
    
    # Distribution comparison plot
    output$comparison_plot <- renderPlotly({
      req(data(), values$normalized_data)
      
      # Use first numeric column for comparison
      numeric_cols <- names(data())[sapply(data(), is.numeric)]
      
      if (length(numeric_cols) > 0) {
        first_col <- numeric_cols[1]
        
        comparison_data <- data.frame(
          Original = data()[[first_col]],
          Normalized = values$normalized_data[[first_col]]
        ) %>%
          gather(key = "Type", value = "Value") %>%
          filter(!is.na(Value) & is.finite(Value))
        
        p <- ggplot(comparison_data, aes(x = Value, fill = Type)) +
          geom_density(alpha = 0.7) +
          scale_fill_manual(values = c("Original" = "steelblue", "Normalized" = "orange")) +
          labs(
            title = paste("Distribution Comparison:", first_col),
            x = "Value",
            y = "Density"
          ) +
          theme_minimal() +
          theme(
            plot.title = element_text(size = 14, face = "bold"),
            legend.position = "top"
          )
        
        ggplotly(p)
      }
    })
    
    # Output flags for conditional panels
    output$preview_available <- reactive({
      !is.null(values$preview_data)
    })
    outputOptions(output, "preview_available", suspendWhenHidden = FALSE)
    
    output$normalization_applied <- reactive({
      values$normalization_applied
    })
    outputOptions(output, "normalization_applied", suspendWhenHidden = FALSE)
    
    # Return normalized data
    reactive({
      if (values$normalization_applied) {
        values$normalized_data
      } else {
        NULL
      }
    })
  })
}

# Helper function to apply normalization methods
apply_normalization_method <- function(data, method, log_base = "natural") {
  
  # Remove non-finite values for calculation
  finite_data <- data[is.finite(data) & !is.na(data)]
  
  if (length(finite_data) == 0) {
    return(data)  # Return original if no finite values
  }
  
  result <- data
  
  if (method == "zscore") {
    # Z-score normalization: (x - mean) / sd
    mean_val <- mean(finite_data)
    sd_val <- sd(finite_data)
    if (sd_val > 0) {
      result <- (data - mean_val) / sd_val
    }
    
  } else if (method == "minmax") {
    # Min-max scaling: (x - min) / (max - min)
    min_val <- min(finite_data)
    max_val <- max(finite_data)
    if (max_val > min_val) {
      result <- (data - min_val) / (max_val - min_val)
    }
    
  } else if (method == "robust") {
    # Robust scaling: (x - median) / MAD
    median_val <- median(finite_data)
    mad_val <- mad(finite_data)
    if (mad_val > 0) {
      result <- (data - median_val) / mad_val
    }
    
  } else if (method == "log") {
    # Log transformation
    if (all(finite_data > 0)) {
      if (log_base == "natural") {
        result <- log(data)
      } else if (log_base == "base2") {
        result <- log2(data)
      } else if (log_base == "base10") {
        result <- log10(data)
      }
    } else {
      # Add small constant to handle zeros/negatives
      min_positive <- min(finite_data[finite_data > 0])
      offset <- min_positive / 2
      
      if (log_base == "natural") {
        result <- log(data + offset)
      } else if (log_base == "base2") {
        result <- log2(data + offset)
      } else if (log_base == "base10") {
        result <- log10(data + offset)
      }
    }
    
  } else if (method == "log_zscore") {
    # Log transformation followed by z-score
    if (all(finite_data > 0)) {
      if (log_base == "natural") {
        log_data <- log(data)
      } else if (log_base == "base2") {
        log_data <- log2(data)
      } else if (log_base == "base10") {
        log_data <- log10(data)
      }
    } else {
      # Add small constant to handle zeros/negatives
      min_positive <- min(finite_data[finite_data > 0])
      offset <- min_positive / 2
      
      if (log_base == "natural") {
        log_data <- log(data + offset)
      } else if (log_base == "base2") {
        log_data <- log2(data + offset)
      } else if (log_base == "base10") {
        log_data <- log10(data + offset)
      }
    }
    
    # Apply z-score to log-transformed data
    finite_log <- log_data[is.finite(log_data) & !is.na(log_data)]
    if (length(finite_log) > 0) {
      mean_val <- mean(finite_log)
      sd_val <- sd(finite_log)
      if (sd_val > 0) {
        result <- (log_data - mean_val) / sd_val
      }
    }
  }
  
  return(result)
}
