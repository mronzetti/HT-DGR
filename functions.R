# ╦ ╦╔╦╗ ╔╦╗╔═╗╦═╗
# ╠═╣ ║───║║║ ╦╠╦╝
# ╩ ╩ ╩  ═╩╝╚═╝╩╚═
#
# Code for processing a csv file into
#
# Written by Michael Ronzetti, NIH/NCATS 2022

library(tidyverse)
library(ggplot2)
library(cowplot)
library(ggthemes)
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
  df$buffer.pH[is.na(df$buffer.pH)] <- 0
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

#' mergeReplicates
#'
#' Merges the replicates across a plate and outputs mean, StDev, and N
#'
#' @param df
#'
#' @return
#' @export
#'
#' @examples
mergeReplicates <- function(df) {
  df <- df %>% dplyr::select(-Well)
  df.Agg <- df %>%
    dplyr::group_by(buffer.Name, buffer.M, buffer.pH, buffer.NaCl.M) %>%
    dplyr::summarise_all(funs(mean, sd, n()))
  return(df.Agg)
}

DGRPlots <- function(df) {
  bufferList <- unique(df$buffer.Name)
  for (i in 1:length(bufferList)) {
    # First, make temp df to store the current buffer
    df.temp <- df %>% filter(buffer.Name == bufferList[i])
    
    # Then, plot out graphs for Tm, Slope, Initial, and DGR Score
    plot.Tm <-
      ggplot(data = df.temp, aes(x = buffer.NaCl.M, y = Tm_mean)) +
      geom_errorbar(aes(ymin = Tm_mean - Tm_sd, ymax = Tm_mean + Tm_sd)) +
      geom_point(
        shape = 21,
        size = 4,
        color = 'black',
        fill = '#CC6677'
      ) +
      theme_clean() +
      labs(title = 'Tm vs. NaCl Concentration',
           x = 'NaCl [M]',
           y = 'Tm') +
      theme(
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 16)
      )
    plot.Tm
    
    plot.Slope <-
      ggplot(data = df.temp, aes(x = buffer.NaCl.M, y = Slope_mean)) +
      geom_errorbar(aes(ymin = Slope_mean - Slope_sd, ymax = Slope_mean +
                          Slope_sd)) +
      geom_point(
        shape = 21,
        size = 4,
        color = 'black',
        fill = '#999933'
      ) +
      theme_clean() +
      labs(title = 'Slope vs. NaCl Concentration',
           x = 'NaCl [M]',
           y = 'Slope') +
      theme(
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 16)
      )
    plot.Slope
    
    plot.Initial <-
      ggplot(data = df.temp, aes(x = buffer.NaCl.M, y = Initial_mean)) +
      geom_errorbar(aes(ymin = Initial_mean - Initial_sd, ymax = Initial_mean +
                          Initial_sd)) +
      geom_point(
        shape = 21,
        size = 4,
        color = 'black',
        fill = '#44AA99'
      ) +
      theme_clean() +
      labs(title = 'Initial vs. NaCl Concentration',
           x = 'NaCl [M]',
           y = 'Initial') +
      theme(
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 16)
      )
    plot.Initial
    
    plot.DGR <-
      ggplot(data = df.temp, aes(x = buffer.NaCl.M, y = DGR.Score_mean)) +
      geom_errorbar(aes(
        ymin = DGR.Score_mean - DGR.Score_sd,
        ymax = DGR.Score_mean + DGR.Score_sd
      )) +
      geom_point(
        shape = 21,
        size = 4,
        color = 'black',
        fill = '#AA4499'
      ) +
      theme_clean() +
      labs(title = 'DGR Score vs. NaCl Concentration',
           x = 'NaCl [M]',
           y = 'DGR Score') +
      theme(
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 16)
      )
    plot.DGR
    
    # Now, assemble in a cowplot
    plot.row1 <- plot_grid(plot.Tm, plot.Slope)
    plot.row2 <- plot_grid(plot.Initial, plot.DGR)
    plot.title <- ggdraw() +
      draw_label(
        paste('Buffer: ', bufferList[i], sep = ''),
        fontface = 'bold',
        size = 18,
        x = 0,
        hjust = 0
      ) +
      theme(# add margin on the left of the drawing canvas,
        # so title is aligned with left edge of first plot
        plot.margin = margin(0, 0, 0, 7))
    plot.full <- plot_grid(
      plot.title,
      plot.row1,
      plot.row2,
      ncol = 1,
      rel_heights = c(0.15, 1, 1)
    )
    print(plot.full)
    ggsave(filename = paste('./figs/',bufferList[[i]],'.png',sep=''), dpi = 'retina', scale = 1.25)
  }
}