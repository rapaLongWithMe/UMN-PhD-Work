# Create Sample Metabolomics Data
# This script generates sample data for testing the Shiny app

library(dplyr)

# Function to create sample metabolomics dataset
create_sample_data <- function(n_samples = 60, n_metabolites = 25, n_timepoints = 6) {
  set.seed(123)  # For reproducibility
  
  # Create sample metadata
  samples_per_timepoint <- n_samples %/% n_timepoints
  actual_samples <- samples_per_timepoint * n_timepoints
  
  sample_data <- data.frame(
    Sample_ID = paste0("S", sprintf("%03d", 1:actual_samples)),
    Time_Point = rep(1:n_timepoints, each = samples_per_timepoint),
    Replicate = rep(1:3, length.out = actual_samples),
    Group = rep(c("Control", "Treatment_A", "Treatment_B"), length.out = actual_samples),
    Batch = rep(1:4, length.out = actual_samples),
    stringsAsFactors = FALSE
  )
  
  # Generate metabolite data with realistic patterns
  metabolite_names <- c(
    # Amino acids
    "Alanine", "Glycine", "Leucine", "Isoleucine", "Valine", "Serine", "Threonine",
    # Organic acids
    "Lactate", "Pyruvate", "Citrate", "Succinate", "Fumarate", "Malate",
    # Sugars
    "Glucose", "Fructose", "Sucrose", "Ribose", "Xylose",
    # Fatty acids
    "Palmitate", "Oleate", "Stearate", "Linoleate",
    # Nucleotides
    "ATP", "ADP", "GTP", "CTP"
  )[1:n_metabolites]
  
  # Initialize metabolite data
  for (i in 1:length(metabolite_names)) {
    metabolite <- metabolite_names[i]
    
    # Base abundance varies by metabolite type
    if (grepl("ine$", metabolite)) {
      # Amino acids: moderate abundance
      base_abundance <- runif(1, 800, 2000)
    } else if (grepl("ate$", metabolite)) {
      # Organic acids: variable abundance
      base_abundance <- runif(1, 500, 3000)
    } else if (metabolite %in% c("Glucose", "Fructose", "Sucrose")) {
      # Major sugars: high abundance
      base_abundance <- runif(1, 2000, 8000)
    } else if (grepl("TP$", metabolite)) {
      # Nucleotides: low to moderate abundance
      base_abundance <- runif(1, 200, 1000)
    } else {
      # Others: variable
      base_abundance <- runif(1, 400, 2500)
    }
    
    # Time effects (some metabolites change over time)
    time_effect_strength <- runif(1, -50, 100)
    time_effects <- sample_data$Time_Point * time_effect_strength
    
    # Group effects (treatment differences)
    group_effects <- case_when(
      sample_data$Group == "Control" ~ 0,
      sample_data$Group == "Treatment_A" ~ runif(1, -300, 400),
      sample_data$Group == "Treatment_B" ~ runif(1, -200, 300),
      TRUE ~ 0
    )
    
    # Batch effects (systematic variation)
    batch_effects <- case_when(
      sample_data$Batch == 1 ~ 0,
      sample_data$Batch == 2 ~ runif(1, -100, 100),
      sample_data$Batch == 3 ~ runif(1, -80, 80),
      sample_data$Batch == 4 ~ runif(1, -120, 120),
      TRUE ~ 0
    )
    
    # Biological variation
    bio_variation <- rnorm(actual_samples, 0, base_abundance * 0.15)
    
    # Technical variation
    tech_variation <- rnorm(actual_samples, 0, base_abundance * 0.05)
    
    # Combine all effects
    abundance <- base_abundance + time_effects + group_effects + 
                batch_effects + bio_variation + tech_variation
    
    # Ensure positive values
    abundance <- pmax(abundance, base_abundance * 0.1)
    
    # Add to dataset
    sample_data[[metabolite]] <- round(abundance, 2)
  }
  
  # Introduce some missing values (realistic data often has gaps)
  metabolite_cols <- metabolite_names
  n_missing <- round(length(metabolite_cols) * actual_samples * 0.02)  # 2% missing
  
  for (missing_count in 1:n_missing) {
    row_idx <- sample(1:actual_samples, 1)
    col_idx <- sample(metabolite_cols, 1)
    sample_data[row_idx, col_idx] <- NA
  }
  
  # Add some extreme outliers (simulate contamination or technical issues)
  n_outliers <- round(length(metabolite_cols) * actual_samples * 0.005)  # 0.5% outliers
  
  for (outlier_count in 1:n_outliers) {
    row_idx <- sample(1:actual_samples, 1)
    col_idx <- sample(metabolite_cols, 1)
    
    current_val <- sample_data[row_idx, col_idx]
    if (!is.na(current_val)) {
      # Make it 5-10x higher or 10x lower
      multiplier <- sample(c(5, 6, 7, 8, 9, 10, 0.1), 1)
      sample_data[row_idx, col_idx] <- current_val * multiplier
    }
  }
  
  return(sample_data)
}

# Generate and save sample datasets

# Small dataset for quick testing
small_data <- create_sample_data(n_samples = 30, n_metabolites = 10, n_timepoints = 3)
write.csv(small_data, "sample_data_small.csv", row.names = FALSE)

# Medium dataset for typical use
medium_data <- create_sample_data(n_samples = 60, n_metabolites = 25, n_timepoints = 6)
write.csv(medium_data, "sample_data_medium.csv", row.names = FALSE)

# Large dataset for stress testing
large_data <- create_sample_data(n_samples = 120, n_metabolites = 50, n_timepoints = 8)
write.csv(large_data, "sample_data_large.csv", row.names = FALSE)

cat("Sample datasets created:\n")
cat("- sample_data_small.csv (30 samples, 10 metabolites)\n")
cat("- sample_data_medium.csv (60 samples, 25 metabolites)\n")
cat("- sample_data_large.csv (120 samples, 50 metabolites)\n")
cat("\nThese files can be used to test the Shiny app functionality.\n")

# Display summary of medium dataset
cat("\nSample of medium dataset:\n")
cat("Dimensions:", nrow(medium_data), "rows ×", ncol(medium_data), "columns\n")
cat("Columns:", paste(names(medium_data)[1:8], collapse = ", "), "...\n")
cat("Missing values:", sum(is.na(medium_data)), "\n")
cat("Completeness:", round((1 - sum(is.na(medium_data)) / (nrow(medium_data) * ncol(medium_data))) * 100, 1), "%\n")
