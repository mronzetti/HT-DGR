# ╦ ╦╔╦╗ ╔╦╗╔═╗╦═╗
# ╠═╣ ║───║║║ ╦╠╦╝
# ╩ ╩ ╩  ═╩╝╚═╝╩╚═
#
# Code for processing a csv file into 
#
# Written by Michael Ronzetti, NIH/NCATS 2022

library(tidyverse)
source('functions.R')

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