## Introduction
## Author: Chau Quang (quang.c@wehi.edu.au)
## Purpose: To produce IRR plots by antigen
## References: Takaaki Yuguchi's Stata code (xi_xtnbreg from excel.do)

setwd("~/Downloads/data/inspire/R")

# Load required packages --------------------------------------------------
library(ggplot2)
library(dplyr)
library(grid)
library(patchwork)
library(tidyverse)

# Import and clean data ---------------------------------------------------
metadata <- read.csv("data/InSPIRE_metadata.csv")
stata_irr <- read.csv("data/InSPIRE_stata_irr_results.csv")

# Forest plot -------------------------------------------------------------
forestplot_irr <- function(data, antigen = NULL, antibody = NULL) {
  
  df_plot <- data %>%
    mutate(
      Antibody = stringr::str_replace_all(Antibody, "Fcg", "Fc\u03B3"),
      Antigen  = stringr::str_replace_all(Antigen,  "Fcg", "Fc\u03B3"),
      sig_size = ifelse(Significance == "ns", 3, 4)
    )
  
  if (!is.null(antibody)) antibody <- stringr::str_replace_all(antibody, "Fcg", "Fc\u03B3")
  if (!is.null(antigen))  antigen  <- stringr::str_replace_all(antigen,  "Fcg", "Fc\u03B3")
  
  if (!is.null(antigen)) {
    df_plot <- df_plot %>% filter(Antigen == antigen) %>% mutate(y_variable = reorder(Antibody, IRR))
    y_label  <- "Antibody"
    title_text <- antigen
  }
  if (!is.null(antibody)) {
    df_plot <- df_plot %>% filter(Antibody == antibody) %>% mutate(y_variable = reorder(Antigen, IRR))
    y_label  <- "Antigen"
    title_text <- antibody
  }
  
  if (nrow(df_plot) == 0) {
    message("No data found for the specified filters.")
    return(NULL)
  }
  
  # format p-values
  df_plot <- df_plot %>% 
    mutate(
      Significance = ifelse(is.na(Significance) | Significance == "", "", Significance),
      p_label = ifelse(P < 0.001, "p < 0.001", paste0("p = ", as.character(round(P, 3)))),
      p_fontface = ifelse(P < 0.05, "bold", "plain"),
      sig_vjust = ifelse(Significance == "ns", -1.4, -0.4)
    )
  
  max_uci <- max(df_plot$UCI, na.rm = TRUE)
  
  break_step <- if (max_uci > 2) 0.5 else 0.1
  x_breaks <- seq(0, ceiling(max_uci / break_step) * break_step, by = break_step)
  
  p <- ggplot(df_plot, aes(x = IRR, y = y_variable)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey20", alpha = 0.7) +
    
    geom_errorbar(aes(xmin = LCI, xmax = UCI), width = 0.2, color = "black", linewidth = 0.7) +
    geom_point(color = "black", size = 3) +
    
    geom_text(
      aes(label = Significance, size = sig_size, vjust = sig_vjust),
      hjust = 0.5, color = "grey5", fontface = "bold", show.legend = FALSE
    ) +
    scale_size_identity() +
    
    scale_x_continuous(
      limits = c(0, max_uci),
      breaks = x_breaks,
      labels = function(x) ifelse(x == 0, "0", sprintf("%.1f", x)),
      expand = c(0, 0)
    ) +
    labs(title = title_text, x = "Incidence Rate Ratio (95% CI)", y = y_label) +
    theme_minimal() +
    theme(
      axis.text.x   = element_text(color = "black"),
      axis.text.y   = element_text(size = 10, face = "plain", color = "black"),
      axis.title.y  = element_text(face = "bold"),
      plot.title    = element_text(size = 12, face = "bold", hjust = 0.5),
      panel.grid    = element_blank(),
      axis.line     = element_line(color = "black"),
      axis.ticks.x  = element_line(color = ifelse(x_breaks == 0, "transparent", "black")),
      axis.ticks.y  = element_line(color = "black"),
      axis.ticks.length = unit(4, "pt"),
      plot.margin   = margin(t = 10, r = 110, b = 10, l = 10)
    ) +
    coord_cartesian(clip = "off")
  
  p_annotations <- purrr::pmap(df_plot, function(p_label, p_fontface, y_variable, ...) {
    annotation_custom(
      grob = grid::textGrob(label = p_label, hjust = 0, gp = grid::gpar(fontface = p_fontface, fontsize = 9, col = "black")),
      xmin = max_uci * 1.05, xmax = max_uci * 1.05, ymin = y_variable, ymax = y_variable
    )
  })
  
  return(p + p_annotations)
}

# Plots
forestplot_irr(stata_irr, antibody = "IgG")
forestplot_irr(stata_irr, antibody = "IgG1")
forestplot_irr(stata_irr, antibody = "IgG3")
forestplot_irr(stata_irr, antibody = "C1q")
forestplot_irr(stata_irr, antibody = "FcgRI")
forestplot_irr(stata_irr, antibody = "FcgRII")
forestplot_irr(stata_irr, antibody = "FcgRIII")

## Experimental
# Generate multipanel (2x4)
forestplot_irr_multipanel <- function(data, ncols = 2, nrows = 4) {
  
  # Normalize "Fcg" to unicode formatting
  unique_antibodies <- data %>%
    mutate(Antibody = stringr::str_replace_all(Antibody, "Fcg", "Fc\u03B3")) %>%
    pull(Antibody) %>%
    unique()
  
  # Generate plots for every antibody feature
  plot_list <- lapply(unique_antibodies, function(ab) {
    forestplot_irr(data = data, antibody = ab)
  })
  
  # Clean out any NULL items in case an antibody returned no rows
  plot_list <- compact(plot_list)
  
  # Wrap using patchwork rules
  multipanel_plot <- wrap_plots(plot_list, ncol = ncols, nrow = nrows)
  
  return(multipanel_plot)
}

forestplot_irr_multipanel(stata_irr)

  