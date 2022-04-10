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
library(cowplot)
library(ggpubr)

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
  df$buffer.Name[df$buffer.Name == ''] <- 'Blank'
  df$buffer.NaCl.M[is.na(df$buffer.NaCl.M)] <- 0
  df$buffer.M[is.na(df$buffer.M)] <- 0
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
  combined.df <- full_join(df.platemap, df, by = 'Well')
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
    mutate(., DGR.Score = (Tm * Slope) / Initial)
  return(df)
}

plotSpreads <- function(df) {
  # Plot Tm vs. DGR Score
  spread.Tm.DGR <-
    ggplot(df,
           aes(
             x = Tm,
             y = DGR.Score,
             fill = buffer.Name,
             size = buffer.NaCl.M
           )) +
    geom_point(shape = 21) +
    theme_clean() +
    labs(title = 'Tm vs. HT-DGR Score',
         fill = 'Buffer Component',
         size = 'NaCl [M]') +
    theme(
      axis.title = element_text(size = 15),
      axis.text = element_text(size = 12),
      title = element_text(size = 14)
    )
  
  print(spread.Tm.DGR)
  ggsave(filename = './figs/spread_TM_DGR.png',
         dpi = 'retina',
         scale = 1.5)
  
  # Plot Initial vs. DGR Score
  spread.Initial.DGR <-
    ggplot(df,
           aes(
             x = Initial,
             y = DGR.Score,
             fill = buffer.Name,
             size = buffer.NaCl.M
           )) +
    geom_point(shape = 21) +
    theme_clean() +
    labs(title = 'Initial Signal vs. HT-DGR Score',
         fill = 'Buffer Component',
         size = 'NaCl [M]') +
    theme(
      axis.title = element_text(size = 15),
      axis.text = element_text(size = 12),
      title = element_text(size = 14)
    )
  
  print(spread.Initial.DGR)
  ggsave(filename = './figs/spread_Initial_DGR.png',
         dpi = 'retina',
         scale = 1.5)
  
  # Plot Slope vs. DGR Score
  spread.Slope.DGR <-
    ggplot(df,
           aes(
             x = Slope,
             y = DGR.Score,
             fill = buffer.Name,
             size = buffer.NaCl.M
           )) +
    geom_point(shape = 21) +
    theme_clean() +
    labs(title = 'First Derivative Max Slope vs. HT-DGR Score',
         fill = 'Buffer Component',
         size = 'NaCl [M]') +
    theme(
      axis.title = element_text(size = 15),
      axis.text = element_text(size = 12),
      title = element_text(size = 14)
    )
  
  print(spread.Slope.DGR)
  ggsave(filename = './figs/spread_Slope_DGR.png',
         dpi = 'retina',
         scale = 1.5)
  
  # Plot Initial vs. Tm Score
  spread.Initial.Tm <-
    ggplot(df,
           aes(
             x = Tm,
             y = Initial,
             fill = buffer.Name,
             size = buffer.NaCl.M
           )) +
    geom_point(shape = 21) +
    theme_clean() +
    labs(title = 'Initial Signal vs. Tm',
         fill = 'Buffer Component',
         size = 'NaCl [M]') +
    theme(
      axis.title = element_text(size = 15),
      axis.text = element_text(size = 12),
      title = element_text(size = 14)
    )
  
  print(spread.Initial.Tm)
  ggsave(filename = './figs/spread_Initial_Tm.png',
         dpi = 'retina',
         scale = 1.5)
  
  # Plot Tm vs. Slope Score
  spread.Slope.Tm <-
    ggplot(df,
           aes(
             x = Tm,
             y = Slope,
             fill = buffer.Name,
             size = buffer.NaCl.M
           )) +
    geom_point(shape = 21) +
    theme_clean() +
    labs(title = 'First Derivative Max Slope vs. Tm',
         fill = 'Buffer Component',
         size = 'NaCl [M]') +
    theme(
      axis.title = element_text(size = 15),
      axis.text = element_text(size = 12),
      title = element_text(size = 14)
    )
  
  print(spread.Slope.Tm)
  ggsave(filename = './figs/spread_Slope_Tm.png',
         dpi = 'retina',
         scale = 1.5)
}