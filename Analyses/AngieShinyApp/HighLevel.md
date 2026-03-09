# Metabolite Abundance Across Time Shiny App

Workflow that will take in raw abundances from mzmine and will clean, normalize, and visualize time series based metabolomics data.

## Standardized Data

Take raw metabolomics data in relative abundance of candidate metabolites. large numbers and you'll have some number of metabolites and to standarize them you'll have multiple replicates. You want to first look through those replicates and identify outliers.

## Calculate Stastics

Includes calculating the mean, std. dev, and coeffient of variation. Once you have those stats there will be two ways to use them. The two was are Visualized Spread and Normalize the Data

### Visualized Spread

We're going to make a series of plots that will show how spread the data is. So we'll be making distributions for each metabolite across replicates. We'll also be making heatmaps to look at the spread and principal component analysis.

### Normalize the Data

We're going to use our previously calculated statistics to normalize the raw abundances or each metabolite. 

#### Visualized Spread second run

Instead of using the raw abundances we're going to use the normlized abundances for each metabolite using the same plots. The purpose of this is then to compare the normalized to the raw distributions to choose which normalization strategy works best for the given data set.

## End Goal

Create publish ready final figure data.

## Tech Stack

- R
- Import of data from mzml data files
- TBD: Graphing Visualization