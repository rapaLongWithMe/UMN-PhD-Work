# Package Installation Script for Metabolomics Shiny App
# Run this script to install all required packages

# List of required packages
required_packages <- c(
  # Core Shiny packages
  "shiny",
  "shinydashboard", 
  "shinyWidgets",
  "DT",
  
  # Data processing
  "dplyr",
  "tidyr",
  "readr",
  "readxl",
  "openxlsx",
  
  # Statistics
  "stats",
  "psych",
  "corrplot",
  
  # Visualization
  "ggplot2",
  "plotly",
  "heatmaply",
  "pheatmap",
  "RColorBrewer",
  "viridis",
  
  # PCA and multivariate analysis
  "FactoMineR",
  "factoextra",
  "cluster",
  
  # Additional utilities
  "scales",
  "gridExtra",
  "cowplot",
  "htmltools",
  "htmlwidgets"
)

# Function to install missing packages
install_missing_packages <- function(packages) {
  # Set CRAN mirror
  options(repos = c(CRAN = "https://cran.rstudio.com/"))
  
  # Get list of installed packages
  installed <- rownames(installed.packages())
  
  # Find missing packages
  missing <- packages[!packages %in% installed]
  
  if (length(missing) > 0) {
    cat("Installing missing packages:", paste(missing, collapse = ", "), "\n")
    
    # Install from CRAN
    install.packages(missing, dependencies = TRUE)
    
    cat("Package installation completed!\n")
  } else {
    cat("All required packages are already installed.\n")
  }
}

# Install missing packages
install_missing_packages(required_packages)

# Load all packages to verify installation
cat("\nVerifying package installation...\n")
failed_packages <- c()

for (pkg in required_packages) {
  tryCatch({
    library(pkg, character.only = TRUE)
    cat("✓", pkg, "\n")
  }, error = function(e) {
    cat("✗", pkg, "- FAILED\n")
    failed_packages <<- c(failed_packages, pkg)
  })
}

if (length(failed_packages) > 0) {
  cat("\nFailed to load packages:", paste(failed_packages, collapse = ", "), "\n")
  cat("Please install these packages manually.\n")
} else {
  cat("\nAll packages loaded successfully! ✓\n")
  cat("You can now run the Shiny app with: shiny::runApp()\n")
}
