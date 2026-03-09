# Simple R script to start the Shiny app
# Usage: R -f start_app.R

cat("🧬 Metabolomics Time Series Analysis App\n")
cat("========================================\n\n")

# Check current directory
current_dir <- getwd()
cat("📁 Current directory:", current_dir, "\n")

# Check for required files
if (!file.exists("app.R")) {
  stop("❌ Error: app.R not found. Please run this from the AngieShinyApp directory.")
}

# Check for sample data
sample_files <- c("sample_data_small.csv", "sample_data_medium.csv", "sample_data_large.csv")
available_samples <- sample_files[file.exists(sample_files)]

if (length(available_samples) > 0) {
  cat("📊 Available sample datasets:\n")
  for (file in available_samples) {
    cat("   •", file, "\n")
  }
} else {
  cat("⚠️  No sample datasets found\n")
}

cat("\n🚀 Starting Shiny app...\n")
cat("📖 Instructions:\n")
cat("   • App will open in your web browser\n")
cat("   • Try uploading sample_data_medium.csv\n") 
cat("   • Press Ctrl+C to stop\n\n")

# Try different ports if 3838 is busy
ports_to_try <- c(3838, 3839, 3840, 8080)

for (port in ports_to_try) {
  cat("🔌 Trying port", port, "...\n")
  
  tryCatch({
    shiny::runApp(
      appDir = ".",
      launch.browser = TRUE,
      host = "127.0.0.1", 
      port = port
    )
    break  # If successful, exit the loop
  }, error = function(e) {
    if (grepl("address already in use", e$message)) {
      cat("   Port", port, "is busy, trying next...\n")
    } else {
      cat("   Error:", e$message, "\n")
      break
    }
  })
}
