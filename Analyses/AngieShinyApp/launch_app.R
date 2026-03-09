# Launch Script for Metabolomics Shiny App
# Run this script to start the application

cat("Starting Metabolomics Time Series Analysis App...\n")
cat("==============================================\n\n")

# Check if we're in the right directory
if (!file.exists("app.R")) {
  stop("Error: app.R not found. Please make sure you're in the correct directory.")
}

# Check for sample data files
sample_files <- c("sample_data_small.csv", "sample_data_medium.csv", "sample_data_large.csv")
available_samples <- sample_files[file.exists(sample_files)]

if (length(available_samples) > 0) {
  cat("Available sample datasets for testing:\n")
  for (file in available_samples) {
    cat(paste0("  - ", file, "\n"))
  }
  cat("\n")
} else {
  cat("No sample datasets found. You can generate them by running:\n")
  cat("  source('generate_test_data.R')\n\n")
}

# Launch the app
cat("Launching app...\n")
cat("The app will open in your default web browser.\n")
cat("To stop the app, press Ctrl+C (or Cmd+C on Mac) in this console.\n\n")

# Run the Shiny app
shiny::runApp(
  appDir = ".",
  launch.browser = TRUE,
  host = "127.0.0.1",
  port = 3838
)
