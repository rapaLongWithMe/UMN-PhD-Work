#!/bin/bash
# Shell script to launch the Metabolomics Shiny App

echo "🧬 Metabolomics Time Series Analysis App"
echo "========================================"
echo ""

# Check if we're in the right directory
if [ ! -f "app.R" ]; then
    echo "❌ Error: app.R not found."
    echo "Please run this script from the AngieShinyApp directory."
    exit 1
fi

echo "📁 Found app.R in current directory"
echo "🚀 Launching Shiny app..."
echo ""
echo "📖 Instructions:"
echo "   - The app will open in your web browser"  
echo "   - Try uploading sample_data_medium.csv to test"
echo "   - Press Ctrl+C to stop the app"
echo ""
echo "⏳ Starting app..."

# Launch the R script
/opt/homebrew/bin/Rscript -e "
setwd('$(pwd)')
cat('✓ Working directory set to:', getwd(), '\n')
cat('✓ Loading Shiny app...\n')
shiny::runApp(launch.browser = TRUE, host = '127.0.0.1', port = 3838)
"
