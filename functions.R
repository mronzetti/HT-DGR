# ╦ ╦╔╦╗ ╔╦╗╔═╗╦═╗
# ╠═╣ ║───║║║ ╦╠╦╝
# ╩ ╩ ╩  ═╩╝╚═╝╩╚═
#
# Code for processing a csv file into
#
# Written by Michael Ronzetti, NIH/NCATS 2022

library(tidyverse)
library(ggplot2)
library(drc)
library(spatialEco)
library(readxl)

#' importData
#'
#' Importer for the Roche LightCycler exported data
#' Need to specify if these are raw curves or the first derivative export
#'
#' @param file
#' df path
#' @return
#' df containing cleaned Roche Thermal Cycler data.
#' @export
#'
#' @examples
#' importData(file = './firstderiv.txt', type = 'FirstDeriv')
importData <-
  function(filePath = '',
           columnFilter = '') {
    # Bit of rough error-checking for the function.
    if (filePath == '') {
      stop('No file path provided for function.')
    }
    if (is_empty(columnFilter) == TRUE) {
      stop('No column filtering chr vector supplied.')
    }
    # Read in csv file and clean columns
    df <- read.csv(file = filePath)
    df <- df %>%
      dplyr::select(any_of(columnFilter))
    message('TSA File imported and cleaned.')
    return(df)
  }

#' plate_assignment
#' Assign compound ids and concentration from platemap
#'
#' @param df
#' @param platemap_file
#'
#' @return
#' @export
#'
#' @examples
plate_assignment <- function(df, platemap_file) {
  platemapParse <- read_xlsx(platemapPath)
  combined.df <- full_join(df,platemapParse,by='Well')
}