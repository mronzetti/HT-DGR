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
    df <- df %>%
      mutate(Tm = as.numeric(Tm)) %>%
      mutate(Slope = as.numeric(Slope)) %>%
      mutate(Initial = as.numeric(Initial))
    message('TSA File imported and cleaned.')
    return(df)
  }

#' Clean NAs and convert to blank
#'
#' @param df 
#'
#' @return
#' @export
#'
#' @examples
cleanBlank <- function(df) {
  df$buffer.Name[is.na(df$buffer.Name)] <- 'Blank'
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
  df.platemap <- read.csv(file = platemap_file)
  combined.df <- full_join(df,df.platemap,by='Well')
  return(combined.df)
}


#' Calculate HT-DGR Score
#'
#' @param df 
#' df containing Tm, Slope, and Initial values
#' 
#' @return
#' @export
#'
#' @examples
calculateDGR <- function(df) {
  df <- df %>%
    mutate(., DGR.Score = (Tm*Slope)/Initial)
  return(df)
}