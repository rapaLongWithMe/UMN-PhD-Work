# Helper Functions for Metabolomics Shiny App
# Utility functions used across modules

# Function to create comparison plot for raw vs normalized data
create_comparison_plot <- function(raw_data, normalized_data) {
  if (is.null(raw_data) || is.null(normalized_data)) {
    return(NULL)
  }
  
  # Select first numeric column for comparison
  numeric_cols <- names(raw_data)[sapply(raw_data, is.numeric)]
  
  if (length(numeric_cols) == 0) {
    return(NULL)
  }
  
  first_col <- numeric_cols[1]
  
  # Create comparison data
  comparison_data <- data.frame(
    Raw = raw_data[[first_col]],
    Normalized = normalized_data[[first_col]],
    Sample = 1:nrow(raw_data)
  ) %>%
    filter(!is.na(Raw) & !is.na(Normalized) & 
           is.finite(Raw) & is.finite(Normalized)) %>%
    gather(key = "Type", value = "Value", Raw, Normalized)
  
  # Create side-by-side box plots
  p <- ggplot(comparison_data, aes(x = Type, y = Value, fill = Type)) +
    geom_boxplot(alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.3, size = 1) +
    scale_fill_manual(values = c("Raw" = "steelblue", "Normalized" = "orange")) +
    labs(
      title = paste("Raw vs Normalized Comparison:", first_col),
      x = "Data Type",
      y = "Value",
      fill = "Data Type"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      legend.position = "bottom"
    )
  
  ggplotly(p)
}

# Function to validate metabolomics data structure
validate_metabolomics_data <- function(data) {
  if (is.null(data) || nrow(data) == 0) {
    return(list(valid = FALSE, message = "No data provided"))
  }
  
  # Check for required column patterns
  col_names <- tolower(names(data))
  
  required_patterns <- list(
    metabolite = c("metabolite", "compound", "feature", "peak"),
    sample = c("sample", "subject", "id"),
    abundance = c("abundance", "intensity", "concentration", "level"),
    time = c("time", "timepoint", "day", "hour", "week"),
    replicate = c("replicate", "rep", "biological", "technical")
  )
  
  missing_patterns <- c()
  
  for (pattern_name in names(required_patterns)) {
    patterns <- required_patterns[[pattern_name]]
    found <- any(sapply(patterns, function(p) any(grepl(p, col_names))))
    
    if (!found) {
      missing_patterns <- c(missing_patterns, pattern_name)
    }
  }
  
  if (length(missing_patterns) > 0) {
    message <- paste("Missing expected column patterns:", 
                    paste(missing_patterns, collapse = ", "))
    return(list(valid = FALSE, message = message))
  }
  
  # Check for numeric data
  numeric_cols <- sum(sapply(data, is.numeric))
  if (numeric_cols == 0) {
    return(list(valid = FALSE, message = "No numeric columns found"))
  }
  
  return(list(valid = TRUE, message = "Data structure appears valid"))
}

# Function to calculate coefficient of variation
calculate_cv <- function(x, na.rm = TRUE) {
  if (na.rm) {
    x <- x[!is.na(x)]
  }
  
  if (length(x) == 0 || mean(x) == 0) {
    return(NA)
  }
  
  (sd(x) / mean(x)) * 100
}

# Function to identify likely metabolite columns
identify_metabolite_columns <- function(data) {
  # Look for numeric columns that might represent metabolite abundances
  numeric_cols <- names(data)[sapply(data, is.numeric)]
  
  # Filter out columns that are likely identifiers or metadata
  metadata_patterns <- c("id", "index", "row", "time", "day", "hour", "week", 
                        "replicate", "rep", "group", "treatment", "condition")
  
  likely_metabolites <- numeric_cols[!grepl(paste(metadata_patterns, collapse = "|"), 
                                           tolower(numeric_cols))]
  
  return(likely_metabolites)
}

# Function to generate quality control metrics
generate_qc_metrics <- function(data) {
  if (is.null(data) || nrow(data) == 0) {
    return(NULL)
  }
  
  metrics <- list(
    total_samples = nrow(data),
    total_variables = ncol(data),
    numeric_variables = sum(sapply(data, is.numeric)),
    missing_values = sum(is.na(data)),
    completeness = 1 - (sum(is.na(data)) / (nrow(data) * ncol(data))),
    numeric_completeness = NA
  )
  
  # Calculate completeness for numeric columns only
  numeric_data <- data[sapply(data, is.numeric)]
  if (ncol(numeric_data) > 0) {
    metrics$numeric_completeness <- 1 - (sum(is.na(numeric_data)) / 
                                        (nrow(numeric_data) * ncol(numeric_data)))
  }
  
  return(metrics)
}

# Function to prepare data for PCA
prepare_pca_data <- function(data) {
  # Select only numeric columns
  numeric_data <- data[sapply(data, is.numeric)]
  
  if (ncol(numeric_data) < 2) {
    return(NULL)
  }
  
  # Remove columns with zero variance
  zero_var_cols <- apply(numeric_data, 2, function(x) var(x, na.rm = TRUE) == 0)
  if (any(zero_var_cols)) {
    numeric_data <- numeric_data[, !zero_var_cols]
  }
  
  # Remove rows with too many missing values (>50%)
  missing_threshold <- ncol(numeric_data) * 0.5
  complete_rows <- rowSums(is.na(numeric_data)) <= missing_threshold
  
  if (sum(complete_rows) < 3) {
    return(NULL)
  }
  
  prepared_data <- numeric_data[complete_rows, ]
  
  # For remaining missing values, use column means for imputation
  for (col in names(prepared_data)) {
    missing_indices <- is.na(prepared_data[[col]])
    if (any(missing_indices)) {
      prepared_data[[col]][missing_indices] <- mean(prepared_data[[col]], na.rm = TRUE)
    }
  }
  
  return(prepared_data)
}

# Function to calculate normalization effectiveness
assess_normalization_effectiveness <- function(raw_data, normalized_data) {
  if (is.null(raw_data) || is.null(normalized_data)) {
    return(NULL)
  }
  
  # Compare coefficient of variation before and after normalization
  numeric_cols <- names(raw_data)[sapply(raw_data, is.numeric)]
  
  assessment <- data.frame(
    Metabolite = character(),
    CV_Raw = numeric(),
    CV_Normalized = numeric(),
    CV_Improvement = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (col in numeric_cols) {
    if (col %in% names(normalized_data)) {
      cv_raw <- calculate_cv(raw_data[[col]])
      cv_norm <- calculate_cv(normalized_data[[col]])
      
      cv_improvement <- if (is.na(cv_raw) || is.na(cv_norm)) NA else cv_raw - cv_norm
      
      assessment <- rbind(assessment, data.frame(
        Metabolite = col,
        CV_Raw = cv_raw,
        CV_Normalized = cv_norm,
        CV_Improvement = cv_improvement,
        stringsAsFactors = FALSE
      ))
    }
  }
  
  return(assessment)
}

# Function to format numbers for display
format_number <- function(x, digits = 3) {
  if (is.na(x) || !is.finite(x)) {
    return("N/A")
  }
  
  if (abs(x) < 0.001) {
    return(formatC(x, format = "e", digits = digits))
  } else {
    return(round(x, digits))
  }
}

# Function to create sample data for testing
create_sample_metabolomics_data <- function(n_samples = 50, n_metabolites = 20, n_timepoints = 5) {
  # Set seed for reproducibility
  set.seed(42)
  
  sample_data <- data.frame(
    Sample_ID = paste0("Sample_", 1:n_samples),
    Time_Point = rep(1:n_timepoints, length.out = n_samples),
    Replicate = rep(1:3, length.out = n_samples),
    Group = rep(c("Control", "Treatment"), length.out = n_samples)
  )
  
  # Generate metabolite data with some patterns and noise
  for (i in 1:n_metabolites) {
    metabolite_name <- paste0("Metabolite_", LETTERS[((i-1) %% 26) + 1], 
                             sprintf("%02d", ceiling(i/26)))
    
    # Create base abundance with time trend
    base_abundance <- 1000 + i * 50
    time_effect <- sample_data$Time_Point * runif(1, -50, 50)
    group_effect <- ifelse(sample_data$Group == "Treatment", 
                          runif(1, -200, 200), 0)
    
    # Add noise
    noise <- rnorm(n_samples, 0, base_abundance * 0.1)
    
    # Combine effects
    abundance <- base_abundance + time_effect + group_effect + noise
    
    # Ensure positive values
    abundance <- pmax(abundance, base_abundance * 0.1)
    
    sample_data[[metabolite_name]] <- abundance
  }
  
  # Introduce some missing values randomly
  metabolite_cols <- grep("Metabolite_", names(sample_data))
  missing_indices <- sample(length(unlist(sample_data[metabolite_cols])), 
                           size = length(unlist(sample_data[metabolite_cols])) * 0.02)
  
  sample_data[metabolite_cols][missing_indices] <- NA
  
  return(sample_data)
}

# Function to save analysis parameters
save_analysis_parameters <- function(parameters, filename = NULL) {
  if (is.null(filename)) {
    filename <- paste0("analysis_parameters_", Sys.Date(), ".json")
  }
  
  # Convert parameters to JSON and save
  # This would require jsonlite package
  # jsonlite::write_json(parameters, filename, pretty = TRUE)
  
  return(filename)
}

# Function to load analysis parameters
load_analysis_parameters <- function(filename) {
  if (!file.exists(filename)) {
    return(NULL)
  }
  
  # Load parameters from JSON
  # This would require jsonlite package
  # parameters <- jsonlite::read_json(filename)
  
  return(NULL)  # Placeholder
}
