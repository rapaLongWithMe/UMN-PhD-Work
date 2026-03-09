# Metabolite Abundance Across Time Shiny App

A comprehensive R Shiny application for processing, analyzing, and visualizing time series metabolomics data from mzMine output.

## Features

### Data Management
- **Import Support**: CSV and Excel file formats from mzMine
- **Data Validation**: Automatic validation of required columns and data quality assessment
- **Quality Control**: Comprehensive data quality metrics and missing value analysis

### Statistical Analysis
- **Outlier Detection**: Multiple methods (Z-score, IQR, Modified Z-score)
- **Normalization**: Various strategies (Z-score, Min-max, Robust scaling, Log transformation)
- **Descriptive Statistics**: Mean, standard deviation, coefficient of variation

### Visualizations
- **Distribution Plots**: Box plots, violin plots, histograms
- **Heatmaps**: Interactive heatmaps with hierarchical clustering
- **PCA Analysis**: Principal component analysis with interactive plots
- **Time Series**: Temporal visualizations with trend analysis

### Export Capabilities
- **Data Export**: Raw and normalized data in CSV/Excel formats
- **Plot Export**: High-resolution plots (PNG, PDF, SVG)
- **Reports**: Comprehensive analysis reports

## Installation

### Prerequisites
- R (version 4.0 or higher)
- RStudio (recommended)

### Package Installation

1. Clone or download this repository
2. Open R/RStudio and set the working directory to the app folder
3. Run the package installation script:

```r
source("install_packages.R")
```

### Required Packages
The app requires the following R packages:
- `shiny`, `shinydashboard`, `shinyWidgets` - Core Shiny framework
- `dplyr`, `tidyr`, `readr`, `readxl` - Data processing
- `ggplot2`, `plotly`, `heatmaply` - Visualization
- `FactoMineR`, `factoextra` - PCA analysis
- `DT`, `openxlsx` - Data tables and export

## Usage

### Running the App

```r
# Launch the app
shiny::runApp()
```

The app will open in your default web browser.

### Workflow

1. **Data Import**: Upload your mzMine output file (CSV or Excel)
2. **Data Quality**: Review data quality metrics and validation results
3. **Outlier Detection**: Identify and optionally remove outlier samples
4. **Normalization**: Apply appropriate normalization strategy
5. **Visualizations**: Generate various plots for data exploration
6. **Comparison**: Compare raw vs normalized data side-by-side
7. **Export**: Download processed data and publication-ready figures

### Data Format Requirements

Your input file should contain the following columns:
- **Metabolite ID/Name**: Unique identifier for each metabolite
- **Sample ID**: Unique identifier for each sample
- **Time point**: Temporal information
- **Abundance values**: Metabolite abundance measurements
- **Replicate information**: Biological or technical replicate identifiers

Optional metadata columns:
- Sample groups
- Experimental conditions
- Additional annotations

### Example Data Structure

| Sample_ID | Time_Point | Replicate | Group | Metabolite_A01 | Metabolite_B02 | ... |
|-----------|------------|-----------|-------|----------------|----------------|-----|
| Sample_1  | 1          | 1         | Control | 1234.5        | 987.2          | ... |
| Sample_2  | 1          | 2         | Control | 1198.7        | 1023.4         | ... |
| Sample_3  | 2          | 1         | Treatment | 1456.8      | 876.5          | ... |

## Features in Detail

### Outlier Detection Methods

1. **Z-Score Method**: Identifies samples with z-scores exceeding a threshold (default: 3)
2. **IQR Method**: Uses interquartile range to identify outliers (default multiplier: 1.5)
3. **Modified Z-Score**: Uses median absolute deviation for robust outlier detection

### Normalization Strategies

1. **Z-Score Normalization**: (x - mean) / standard deviation
2. **Min-Max Scaling**: (x - min) / (max - min)
3. **Robust Scaling**: (x - median) / MAD
4. **Log Transformation**: Natural log, log2, or log10
5. **Log + Z-Score**: Log transformation followed by z-score normalization

### Visualization Options

- **Distribution Analysis**: Compare distributions across replicates and conditions
- **Heatmaps**: Hierarchical clustering of samples and metabolites
- **PCA Plots**: Sample scores, variable loadings, and biplots
- **Time Series**: Temporal patterns with optional trend lines and confidence intervals

## Performance Considerations

- **Large Datasets**: Supports up to 50,000 metabolites and 1,000 samples
- **Memory Management**: Efficient handling of large datasets
- **Interactive Plots**: Responsive visualizations with zoom and pan capabilities
- **Export Quality**: Publication-ready outputs at 300+ DPI

## Troubleshooting

### Common Issues

1. **Package Installation Errors**
   - Ensure you have the latest version of R
   - Install packages individually if batch installation fails
   - Check for system dependencies (especially for spatial packages)

2. **Data Import Issues**
   - Verify file format (CSV or Excel)
   - Check column names and data types
   - Ensure no special characters in column names

3. **Memory Issues**
   - Reduce dataset size by sampling
   - Close other R sessions
   - Increase available memory if possible

4. **Plot Generation Errors**
   - Ensure data contains numeric columns
   - Check for infinite or extreme values
   - Verify sufficient data points for analysis

### Getting Help

If you encounter issues:
1. Check the console for error messages
2. Verify your data format matches requirements
3. Try with a smaller dataset first
4. Consult R documentation for specific packages

## Development

### Project Structure

```
AngieShinyApp/
├── app.R                 # Main application file
├── install_packages.R    # Package installation script
├── modules/              # Shiny modules
│   ├── data_import.R
│   ├── data_validation.R
│   ├── outlier_detection.R
│   ├── normalization.R
│   ├── visualization.R
│   └── export.R
├── utils/                # Utility functions
│   └── helper_functions.R
├── www/                  # Web assets
│   └── custom.css
├── plan.md              # Technical plan
├── prd.md               # Product requirements
└── README.md            # This file
```

### Contributing

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Test thoroughly
5. Submit a pull request

### Future Enhancements

- Integration with metabolomics databases (HMDB, KEGG)
- Machine learning models for metabolite prediction
- Pathway analysis capabilities
- Multi-omics data integration
- Cloud deployment options

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Citation

If you use this software in your research, please cite:

```
Metabolite Abundance Across Time Shiny App
GitHub Repository: https://github.com/[your-username]/AngieShinyApp
```

## Acknowledgments

- Built with R Shiny framework
- Uses ggplot2 and plotly for visualizations
- Statistical methods implemented with base R and specialized packages
- UI design inspired by modern web applications

## Contact

For questions, suggestions, or bug reports, please:
- Open an issue on GitHub
- Contact the development team
- Consult the documentation

---

**Version**: 1.0.0  
**Last Updated**: December 2024
