# Quick Start Guide

## Running the Metabolomics Shiny App

### Option 1: Easy Launch (Recommended)
1. Open R or RStudio
2. Set your working directory to this folder:
   ```r
   setwd("/Users/rebecca/projects/AngieShinyApp")
   ```
3. Run the launch script:
   ```r
   source("launch_app.R")
   ```

### Option 2: Manual Launch
1. Open R or RStudio
2. Set your working directory to this folder
3. Run:
   ```r
   shiny::runApp()
   ```

## Testing the App

### Using Sample Data
The app comes with three sample datasets:

1. **sample_data_small.csv** - 30 samples, 10 metabolites (quick testing)
2. **sample_data_medium.csv** - 60 samples, 20 metabolites (typical usage)  
3. **sample_data_large.csv** - 120 samples, 40 metabolites (stress testing)

### App Workflow
1. **Data Import Tab**: Upload one of the sample CSV files
2. **Data Quality Tab**: Review the data quality metrics
3. **Outlier Detection Tab**: Try different outlier detection methods
4. **Normalization Tab**: Apply normalization and see the preview
5. **Visualizations Tab**: Generate different plot types
6. **Export Tab**: Download your results

## Troubleshooting

### If the app won't start:
1. Make sure all packages are installed:
   ```r
   source("install_packages.R")
   ```

2. Check you're in the right directory:
   ```r
   getwd()  # Should show: /Users/rebecca/projects/AngieShinyApp
   list.files()  # Should show app.R and other files
   ```

3. Try loading the app manually:
   ```r
   source("app.R")
   ```

### If you get package errors:
Run the package installation script:
```r
source("install_packages.R")
```

### If data won't upload:
- Make sure the file is CSV or Excel format
- Check that column names don't have special characters
- Try with the provided sample data first

## Need Help?
- Check the console for error messages
- Look at the README.md for detailed documentation
- Try with smaller sample datasets first
