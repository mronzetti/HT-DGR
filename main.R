# ╦ ╦╔╦╗ ╔╦╗╔═╗╦═╗
# ╠═╣ ║───║║║ ╦╠╦╝
# ╩ ╩ ╩  ═╩╝╚═╝╩╚═
#
# Code for processing a csv file into 
#
# Written by Michael Ronzetti, NIH/NCATS 2022
library(tidyverse)
library(ggthemes)
source('functions.R')

# Include paths to raw data file, platemap file, and columns to be filtered (don't recommend touching that).
filePath <- 'data/sample_refolding.csv'
platemapPath <- 'data/platemap.csv'
columnFilter <- c('Well','Tm','Slope','Initial')

# First, we import the data from data files exported from Roche TSA Analysis software.
# Requires: Roche .csv file
df.Raw <- importData(filePath = filePath, columnFilter = columnFilter)

# Assign identifiers from the platemap for buffer and salt conditions
df.Raw <- plate_assignment(df.Raw, platemap_file = platemapPath)

# Change any NA values to be blank
df.Raw <- cleanBlank(df.Raw)

# Calculate HT-DGR scores for the entire plate
df.Raw <- calculateDGR(df.Raw)

# Plot out distributions of scores
plotSpreads(df.Raw)

# Merge replicates and perform mean, stdev, n calculations
df.Agg <- mergeReplicates(df.Raw)

# Plot out all graphs for each buffer condition.
DGRPlots(df.Agg)
