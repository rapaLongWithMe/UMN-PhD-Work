# Visualization Module
# Handles all data visualization functionality

visualizationUI <- function(id) {
  ns <- NS(id)
  
  tagList(
    conditionalPanel(
      condition = paste0("output['", ns("has_data"), "']"),
      
      fluidRow(
        column(3,
          h4("Visualization Controls"),
          
          selectInput(ns("plot_type"), "Plot Type:",
                     choices = list(
                       "Distribution Plots" = "distribution",
                       "Heatmap" = "heatmap",
                       "PCA Analysis" = "pca",
                       "Time Series" = "timeseries"
                     ),
                     selected = "distribution"),
          
          hr(),
          
          # Distribution plot controls
          conditionalPanel(
            condition = paste0("input['", ns("plot_type"), "'] == 'distribution'"),
            h5("Distribution Settings"),
            selectInput(ns("dist_metabolite"), "Metabolite:",
                       choices = NULL),
            radioButtons(ns("dist_plot_type"), "Plot Type:",
                        choices = list("Box Plot" = "box",
                                     "Violin Plot" = "violin",
                                     "Histogram" = "histogram"),
                        selected = "box"),
            checkboxInput(ns("show_points"), "Show individual points", value = TRUE)
          ),
          
          # Heatmap controls
          conditionalPanel(
            condition = paste0("input['", ns("plot_type"), "'] == 'heatmap'"),
            h5("Heatmap Settings"),
            checkboxInput(ns("cluster_rows"), "Cluster rows", value = TRUE),
            checkboxInput(ns("cluster_cols"), "Cluster columns", value = TRUE),
            selectInput(ns("color_scale"), "Color Scale:",
                       choices = list("Viridis" = "viridis",
                                    "Blue-Red" = "RdBu",
                                    "Heat" = "heat.colors"),
                       selected = "viridis"),
            numericInput(ns("heatmap_rows"), "Max rows to display:", 
                        value = 50, min = 10, max = 500)
          ),
          
          # PCA controls
          conditionalPanel(
            condition = paste0("input['", ns("plot_type"), "'] == 'pca'"),
            h5("PCA Settings"),
            selectInput(ns("pca_type"), "PCA Plot:",
                       choices = list("Samples" = "samples",
                                    "Variables" = "variables",
                                    "Biplot" = "biplot"),
                       selected = "samples"),
            numericInput(ns("pc_x"), "PC X-axis:", value = 1, min = 1, max = 10),
            numericInput(ns("pc_y"), "PC Y-axis:", value = 2, min = 1, max = 10)
          ),
          
          # Time series controls
          conditionalPanel(
            condition = paste0("input['", ns("plot_type"), "'] == 'timeseries'"),
            h5("Time Series Settings"),
            selectInput(ns("ts_metabolites"), "Metabolites:",
                       choices = NULL,
                       multiple = TRUE),
            checkboxInput(ns("show_trend"), "Show trend lines", value = TRUE),
            checkboxInput(ns("show_confidence"), "Show confidence intervals", value = FALSE)
          ),
          
          hr(),
          
          actionButton(ns("generate_plot"), "Generate Plot", 
                      class = "btn-primary", icon = icon("chart-line")),
          
          br(), br(),
          
          conditionalPanel(
            condition = paste0("output['", ns("plot_generated"), "']"),
            downloadButton(ns("download_plot"), "Download Plot", 
                          class = "btn-success")
          )
        ),
        
        column(9,
          conditionalPanel(
            condition = paste0("output['", ns("plot_generated"), "']"),
            h4("Visualization"),
            
            conditionalPanel(
              condition = paste0("input['", ns("plot_type"), "'] == 'heatmap'"),
              plotlyOutput(ns("heatmap_plot"), height = "600px")
            ),
            
            conditionalPanel(
              condition = paste0("input['", ns("plot_type"), "'] != 'heatmap'"),
              plotlyOutput(ns("main_plot"), height = "600px")
            ),
            
            br(),
            
            # Show PCA summary for PCA plots
            conditionalPanel(
              condition = paste0("input['", ns("plot_type"), "'] == 'pca' && output['", ns("plot_generated"), "']"),
              h4("PCA Summary"),
              DTOutput(ns("pca_summary"))
            )
          )
        )
      )
    ),
    
    conditionalPanel(
      condition = paste0("!output['", ns("has_data"), "']"),
      div(class = "alert alert-info",
          icon("info-circle"), " Please complete data normalization first.")
    )
  )
}

visualizationServer <- function(id, raw_data, normalized_data) {
  moduleServer(id, function(input, output, session) {
    
    # Reactive values
    values <- reactiveValues(
      current_plot = NULL,
      pca_results = NULL,
      plot_generated = FALSE
    )
    
    # Current data (use normalized if available, otherwise raw)
    current_data <- reactive({
      if (!is.null(normalized_data())) {
        normalized_data()
      } else if (!is.null(raw_data())) {
        raw_data()
      } else {
        NULL
      }
    })
    
    # Check if data exists
    output$has_data <- reactive({
      !is.null(current_data()) && nrow(current_data()) > 0
    })
    outputOptions(output, "has_data", suspendWhenHidden = FALSE)
    
    # Update metabolite choices
    observe({
      req(current_data())
      
      numeric_cols <- names(current_data())[sapply(current_data(), is.numeric)]
      
      updateSelectInput(session, "dist_metabolite",
                       choices = numeric_cols,
                       selected = numeric_cols[1])
      
      updateSelectInput(session, "ts_metabolites",
                       choices = numeric_cols,
                       selected = numeric_cols[1:min(3, length(numeric_cols))])
    })
    
    # Generate plot
    observeEvent(input$generate_plot, {
      req(current_data(), input$plot_type)
      
      withProgress(message = "Generating visualization...", value = 0, {
        
        if (input$plot_type == "distribution") {
          setProgress(0.3, detail = "Creating distribution plot...")
          values$current_plot <- create_distribution_plot(
            current_data(), input$dist_metabolite, input$dist_plot_type, input$show_points
          )
          
        } else if (input$plot_type == "heatmap") {
          setProgress(0.3, detail = "Creating heatmap...")
          values$current_plot <- create_heatmap_plot(
            current_data(), input$cluster_rows, input$cluster_cols, 
            input$color_scale, input$heatmap_rows
          )
          
        } else if (input$plot_type == "pca") {
          setProgress(0.3, detail = "Performing PCA analysis...")
          pca_result <- perform_pca_analysis(current_data())
          values$pca_results <- pca_result
          
          setProgress(0.7, detail = "Creating PCA plot...")
          values$current_plot <- create_pca_plot(
            pca_result, input$pca_type, input$pc_x, input$pc_y
          )
          
        } else if (input$plot_type == "timeseries") {
          setProgress(0.3, detail = "Creating time series plot...")
          values$current_plot <- create_timeseries_plot(
            current_data(), input$ts_metabolites, input$show_trend, input$show_confidence
          )
        }
        
        setProgress(1, detail = "Complete!")
        values$plot_generated <- TRUE
      })
    })
    
    # Main plot output
    output$main_plot <- renderPlotly({
      req(values$current_plot, input$plot_type != "heatmap")
      values$current_plot
    })
    
    # Heatmap plot output
    output$heatmap_plot <- renderPlotly({
      req(values$current_plot, input$plot_type == "heatmap")
      values$current_plot
    })
    
    # PCA summary table
    output$pca_summary <- renderDT({
      req(values$pca_results, input$plot_type == "pca")
      
      if (!is.null(values$pca_results$variance_explained)) {
        summary_data <- data.frame(
          Component = paste0("PC", 1:length(values$pca_results$variance_explained)),
          Variance_Explained = round(values$pca_results$variance_explained, 3),
          Cumulative_Variance = round(cumsum(values$pca_results$variance_explained), 3)
        )
        
        datatable(
          summary_data,
          options = list(
            pageLength = 10,
            dom = 't'
          ),
          rownames = FALSE
        ) %>%
          formatPercentage(columns = c("Variance_Explained", "Cumulative_Variance"), digits = 1)
      }
    })
    
    # Download handler
    output$download_plot <- downloadHandler(
      filename = function() {
        paste0("metabolomics_", input$plot_type, "_", Sys.Date(), ".png")
      },
      content = function(file) {
        req(values$current_plot)
        
        if (input$plot_type == "heatmap") {
          # Special handling for heatmap
          ggsave(file, values$current_plot, width = 12, height = 8, dpi = 300)
        } else {
          # For plotly objects, convert to static image
          orca(values$current_plot, file, width = 1200, height = 800)
        }
      }
    )
    
    # Output flag for conditional panels
    output$plot_generated <- reactive({
      values$plot_generated
    })
    outputOptions(output, "plot_generated", suspendWhenHidden = FALSE)
  })
}

# Helper functions for plot creation

create_distribution_plot <- function(data, metabolite, plot_type, show_points) {
  req(data, metabolite)
  
  plot_data <- data.frame(
    Value = data[[metabolite]],
    Sample = 1:nrow(data)
  ) %>%
    filter(!is.na(Value) & is.finite(Value))
  
  if (plot_type == "box") {
    p <- ggplot(plot_data, aes(x = "", y = Value)) +
      geom_boxplot(fill = "steelblue", alpha = 0.7)
    
    if (show_points) {
      p <- p + geom_jitter(width = 0.2, alpha = 0.5)
    }
    
  } else if (plot_type == "violin") {
    p <- ggplot(plot_data, aes(x = "", y = Value)) +
      geom_violin(fill = "steelblue", alpha = 0.7)
    
    if (show_points) {
      p <- p + geom_jitter(width = 0.2, alpha = 0.5)
    }
    
  } else if (plot_type == "histogram") {
    p <- ggplot(plot_data, aes(x = Value)) +
      geom_histogram(fill = "steelblue", alpha = 0.7, bins = 30)
  }
  
  p <- p +
    labs(
      title = paste("Distribution of", metabolite),
      x = if (plot_type == "histogram") metabolite else "",
      y = if (plot_type == "histogram") "Frequency" else metabolite
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      axis.title = element_text(size = 12)
    )
  
  ggplotly(p)
}

create_heatmap_plot <- function(data, cluster_rows, cluster_cols, color_scale, max_rows) {
  req(data)
  
  # Select numeric columns
  numeric_data <- data[sapply(data, is.numeric)]
  
  if (ncol(numeric_data) == 0) {
    return(NULL)
  }
  
  # Limit rows for performance
  if (nrow(numeric_data) > max_rows) {
    numeric_data <- numeric_data[1:max_rows, ]
  }
  
  # Remove rows/columns with all NA
  numeric_data <- numeric_data[rowSums(!is.na(numeric_data)) > 0, ]
  numeric_data <- numeric_data[, colSums(!is.na(numeric_data)) > 0]
  
  if (nrow(numeric_data) == 0 || ncol(numeric_data) == 0) {
    return(NULL)
  }
  
  # Create heatmap
  heatmaply(
    as.matrix(numeric_data),
    colors = if (color_scale == "viridis") viridis::viridis(100) 
             else if (color_scale == "RdBu") RColorBrewer::brewer.pal(11, "RdBu")
             else heat.colors(100),
    dendrogram = if (cluster_rows && cluster_cols) "both" 
                else if (cluster_rows) "row" 
                else if (cluster_cols) "column" 
                else "none",
    main = "Metabolite Abundance Heatmap",
    margins = c(60, 100, 40, 20)
  )
}

perform_pca_analysis <- function(data) {
  req(data)
  
  # Select numeric columns
  numeric_data <- data[sapply(data, is.numeric)]
  
  if (ncol(numeric_data) < 2) {
    return(NULL)
  }
  
  # Remove rows with missing values
  complete_data <- numeric_data[complete.cases(numeric_data), ]
  
  if (nrow(complete_data) < 3) {
    return(NULL)
  }
  
  # Perform PCA
  pca_result <- prcomp(complete_data, scale. = TRUE, center = TRUE)
  
  # Calculate variance explained
  variance_explained <- (pca_result$sdev^2) / sum(pca_result$sdev^2)
  
  list(
    pca = pca_result,
    data = complete_data,
    variance_explained = variance_explained
  )
}

create_pca_plot <- function(pca_results, plot_type, pc_x, pc_y) {
  req(pca_results)
  
  if (plot_type == "samples") {
    # Scores plot
    scores_data <- data.frame(
      PC_X = pca_results$pca$x[, pc_x],
      PC_Y = pca_results$pca$x[, pc_y],
      Sample = 1:nrow(pca_results$pca$x)
    )
    
    p <- ggplot(scores_data, aes(x = PC_X, y = PC_Y, text = paste("Sample:", Sample))) +
      geom_point(size = 3, alpha = 0.7, color = "steelblue") +
      labs(
        title = "PCA Scores Plot",
        x = paste0("PC", pc_x, " (", round(pca_results$variance_explained[pc_x] * 100, 1), "%)"),
        y = paste0("PC", pc_y, " (", round(pca_results$variance_explained[pc_y] * 100, 1), "%)")
      ) +
      theme_minimal() +
      theme(plot.title = element_text(size = 14, face = "bold"))
    
  } else if (plot_type == "variables") {
    # Loadings plot
    loadings_data <- data.frame(
      PC_X = pca_results$pca$rotation[, pc_x],
      PC_Y = pca_results$pca$rotation[, pc_y],
      Variable = rownames(pca_results$pca$rotation)
    )
    
    p <- ggplot(loadings_data, aes(x = PC_X, y = PC_Y, text = Variable)) +
      geom_point(size = 3, alpha = 0.7, color = "orange") +
      labs(
        title = "PCA Loadings Plot",
        x = paste0("PC", pc_x, " (", round(pca_results$variance_explained[pc_x] * 100, 1), "%)"),
        y = paste0("PC", pc_y, " (", round(pca_results$variance_explained[pc_y] * 100, 1), "%)")
      ) +
      theme_minimal() +
      theme(plot.title = element_text(size = 14, face = "bold"))
  }
  
  ggplotly(p, tooltip = "text")
}

create_timeseries_plot <- function(data, metabolites, show_trend, show_confidence) {
  req(data, metabolites)
  
  # Try to identify time column
  time_cols <- names(data)[grepl("time|day|hour|week", names(data), ignore.case = TRUE)]
  
  if (length(time_cols) == 0) {
    # Use row number as time if no time column found
    time_data <- 1:nrow(data)
    time_label <- "Sample Index"
  } else {
    time_data <- data[[time_cols[1]]]
    time_label <- time_cols[1]
  }
  
  # Create long format data
  plot_data <- data.frame(
    Time = rep(time_data, length(metabolites)),
    Value = unlist(data[metabolites]),
    Metabolite = rep(metabolites, each = nrow(data))
  ) %>%
    filter(!is.na(Value) & is.finite(Value))
  
  p <- ggplot(plot_data, aes(x = Time, y = Value, color = Metabolite)) +
    geom_line(alpha = 0.7, size = 1) +
    geom_point(alpha = 0.7, size = 2)
  
  if (show_trend) {
    p <- p + geom_smooth(method = "loess", se = show_confidence, alpha = 0.3)
  }
  
  p <- p +
    labs(
      title = "Metabolite Time Series",
      x = time_label,
      y = "Abundance"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      legend.position = "bottom"
    )
  
  ggplotly(p)
}
