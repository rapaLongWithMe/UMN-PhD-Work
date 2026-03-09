# Simple Sample Data Generator for Metabolomics Shiny App

library(dplyr)

# Function to create simple test data
create_test_data <- function(n_samples = 60, n_metabolites = 20) {
  set.seed(42)  # For reproducibility
  
  # Create basic metadata
  sample_data <- data.frame(
    Sample_ID = paste0("Sample_", sprintf("%03d", 1:n_samples)),
    Time_Point = rep(c(0, 2, 4, 8, 12, 24), length.out = n_samples),
    Replicate = rep(1:3, length.out = n_samples),
    Group = rep(c("Control", "Treatment"), length.out = n_samples),
    stringsAsFactors = FALSE
  )
  
  # Generate metabolite names
  metabolite_names <- paste0("Metabolite_", LETTERS[1:n_metabolites])
  
  # Generate abundance data for each metabolite
  for (i in seq_along(metabolite_names)) {
    metabolite_name <- metabolite_names[i]
    
    # Base abundance varies by metabolite
    base_abundance <- runif(1, 500, 3000)
    
    # Time effect (some metabolites change over time)
    time_effect <- sample_data$Time_Point * runif(1, -20, 50)
    
    # Group effect (treatment differences)
    group_effect <- ifelse(sample_data$Group == "Treatment", 
                          runif(1, -300, 300), 0)
    
    # Random biological variation
    bio_variation <- rnorm(n_samples, 0, base_abundance * 0.2)
    
    # Combine effects
    abundance <- base_abundance + time_effect + group_effect + bio_variation
    
    # Ensure positive values
    abundance <- pmax(abundance, base_abundance * 0.1)
    
    # Add to dataset
    sample_data[[metabolite_name]] <- round(abundance, 2)
  }
  
  # Add some missing values (2% of data)
  total_cells <- n_samples * n_metabolites
  missing_count <- round(total_cells * 0.02)
  
  for (i in 1:missing_count) {
    row_idx <- sample(1:n_samples, 1)
    col_idx <- sample(metabolite_names, 1)
    sample_data[row_idx, col_idx] <- NA
  }
  
  # Add a few outliers
  outlier_count <- round(total_cells * 0.005)
  
  for (i in 1:outlier_count) {
    row_idx <- sample(1:n_samples, 1)
    col_idx <- sample(metabolite_names, 1)
    
    current_val <- sample_data[row_idx, col_idx]
    if (!is.na(current_val)) {
      # Make it much higher or lower
      multiplier <- sample(c(5, 7, 10, 0.1), 1)
      sample_data[row_idx, col_idx] <- current_val * multiplier
    }
  }
  
  return(sample_data)
}

# Generate sample datasets
cat("Generating sample datasets...\n")

# Small dataset for quick testing
small_data <- create_test_data(n_samples = 30, n_metabolites = 10)
write.csv(small_data, "sample_data_small.csv", row.names = FALSE)
cat("✓ Created sample_data_small.csv (30 samples, 10 metabolites)\n")

# Medium dataset for typical use
medium_data <- create_test_data(n_samples = 60, n_metabolites = 20)
write.csv(medium_data, "sample_data_medium.csv", row.names = FALSE)
cat("✓ Created sample_data_medium.csv (60 samples, 20 metabolites)\n")

# Large dataset for stress testing
large_data <- create_test_data(n_samples = 120, n_metabolites = 40)
write.csv(large_data, "sample_data_large.csv", row.names = FALSE)
cat("✓ Created sample_data_large.csv (120 samples, 40 metabolites)\n")

# Display summary of medium dataset
cat("\nSample of medium dataset:\n")
cat("Dimensions:", nrow(medium_data), "rows ×", ncol(medium_data), "columns\n")
cat("Columns:", paste(names(medium_data)[1:8], collapse = ", "), "...\n")
cat("Missing values:", sum(is.na(medium_data)), "\n")
completeness <- round((1 - sum(is.na(medium_data)) / (nrow(medium_data) * ncol(medium_data))) * 100, 1)
cat("Completeness:", completeness, "%\n")

cat("\nSample data generation completed successfully!\n")
cat("Files ready for testing the Shiny app.\n")
