## Introduction
## Author: Chau Quang (quang.c@wehi.edu.au)
## Purpose: To produce multipanel (4x2) RAU, MFI boxplots of Taka's InSPIRE data

setwd("~/Downloads/data/inspire/R")

# Load required packages --------------------------------------------------
library(ggplot2)

# Import data (.csv) ------------------------------------------------------
data_IgG  <- transform(read.csv("data/InSPIRE_IgG_grouped.csv"),  assay = "IgG")
data_IgG1 <- transform(read.csv("data/InSPIRE_IgG1_grouped.csv"), assay = "IgG1")
data_IgG3 <- transform(read.csv("data/InSPIRE_IgG3_grouped.csv"), assay = "IgG3")
data_C1q <- transform(read.csv("data/InSPIRE_C1q_grouped.csv"), assay = "C1q")
data_FcgRI <- transform(read.csv("data/InSPIRE_FcgRI_grouped.csv"), assay = "FcgRI")
data_FcgRII <- transform(read.csv("data/InSPIRE_FcgRII_grouped.csv"), assay = "FcgRII")
data_FcgRIII <- transform(read.csv("data/InSPIRE_FcgRIII_grouped.csv"), assay = "FcgRIII")

## Debug datasets: tertiling performed via Stata using XOA IDs as input only (Pk standards, VBDR negative samples removed)
data_IgG  <- transform(read.csv("data/testing/InSPIRE_IgG_grouped_xoa.csv"),  assay = "IgG")
data_IgG1 <- transform(read.csv("data/testing/InSPIRE_IgG1_grouped_xoa.csv"), assay = "IgG1")
data_IgG3 <- transform(read.csv("data/testing/InSPIRE_IgG3_grouped_xoa.csv"), assay = "IgG3")
data_C1q <- transform(read.csv("data/testing/InSPIRE_C1q_grouped_xoa.csv"), assay = "C1q")
data_FcgRI <- transform(read.csv("data/testing/InSPIRE_FcgRI_grouped_xoa.csv"), assay = "FcgRI")
data_FcgRII <- transform(read.csv("data/testing/InSPIRE_FcgRII_grouped_xoa.csv"), assay = "FcgRII")
data_FcgRIII <- transform(read.csv("data/testing/InSPIRE_FcgRIII_grouped_xoa.csv"), assay = "FcgRIII")

# Boxplot function --------------------------------------------------------
boxplot_grouped <- function(data, var_name, show_points = TRUE, show_x = TRUE) {
  
  set.seed(42)
  
  antigen <- gsub("_RAU|_MFI", "", var_name)
  assay <- data$assay[1]
  measure <- sub(".*_", "", var_name)
  group_col <- gsub("_RAU|_MFI", "_Group", var_name)
  
  # Clean + rename groups
  data <- data |>
    dplyr::mutate(
      group = factor(
        .data[[group_col]],
        levels = c(1,2,3),
        labels = c("G1","G2","G3")
      )
    ) |>
    dplyr::filter(
      !is.na(.data[[var_name]]),
      .data[[var_name]] > 0,
      !is.na(group)
    )
  
  # y-axis limits
  yvals <- data[[var_name]]
  ymin <- 10^floor(log10(min(yvals)))
  ymax <- 10^ceiling(log10(max(yvals)))
  
  # Toggle x-axis title label
  x_label <- if (show_x) "Group" else NULL
  
  # Adjust y-axis label for FcR plots
  y_label <- if (grepl("FcgR", assay)) {
    receptor <- sub("FcgR", "", assay)
    
    paste0(antigen, " Fc\u03B3R", receptor, " (", measure, ")")
  } else {
    paste0(antigen, " ", assay, " (", measure, ")")
  }
  
  # Boxplot
  p <- ggplot(
    data,
    aes(x = group,
       y = .data[[var_name]],
       fill = group)) +
    labs(x = x_label, 
         y = y_label) +
    theme_bw() +
    theme(
      legend.position = "none",
      panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
      panel.grid.minor = element_line(color = "grey92"),
      axis.line = element_line(color = "black", linewidth = 0.8),
      axis.title.y = element_text(face = "bold"),
      axis.text = element_text(size = 9.5),
      axis.ticks = element_line(linewidth = 0.8, color = "black"),
      plot.margin = margin(5,5,5,5)) +
    scale_fill_manual(values = c(
      "G1" = "yellowgreen",
      "G2" = "darkcyan",
      "G3" = "firebrick2")) +
    scale_y_log10(
      limits = c(ymin, ymax),
      breaks = 10^(floor(log10(min(yvals))):ceiling(log10(max(yvals)))),
      labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    annotation_logticks(sides = "l")
  
  # Toggle individual data points
  if (show_points) {
    
    p <- p +
      geom_boxplot(
        width = 0.5,
        staplewidth = 0.4,
        linewidth = 0.7,
        outlier.shape = NA
      ) +
      geom_jitter(
        aes(colour = group),
        width = 0.15,
        alpha = 0.4,
        size = 2
      ) +
      scale_colour_manual(values = c(
        "G1" = "darkolivegreen",
        "G2" = "deepskyblue4",
        "G3" = "brown"
      ))
    
  } else {
    
    p <- p +
      geom_boxplot(
        aes(colour = group),
        width = 0.5,
        staplewidth = 0.4,
        linewidth = 0.7
      ) +
      scale_colour_manual(values = c(
        "G1" = "darkgreen",
        "G2" = "royalblue4",
        "G3" = "red4"
      ))
  }
  
  return(p)
}

# Test plots
boxplot_grouped(data_IgG, "Pv41_RAU", show_points = FALSE)
boxplot_grouped(data_IgG1, "Pv41_RAU", show_points = FALSE)
boxplot_grouped(data_FcgRI, "PvEXP1_MFI", show_points = FALSE)

# Multipanel (4x2) --------------------------------------------------------
library(patchwork)

# Define antigen groups for plotting
antigens_RAU <- c("Pv12_RAU", "Pv32_RAU", "Pv41_RAU", "Pv50_RAU", 
              "PvEXP1_RAU", "PvMSP10_RAU", "PvRBP2p1_RAU")
antigens_MFI <- c("Pv12_MFI", "Pv32_MFI", "Pv41_MFI", "Pv50_MFI", 
                  "PvEXP1_MFI", "PvMSP10_MFI", "PvRBP2p1_MFI")

# RAU
IgG_plots <- lapply(antigens_RAU, function(x) boxplot_grouped(data_IgG, x, show_points = FALSE, show_x = FALSE))
IgG1_plots <- lapply(antigens_RAU, function(x) boxplot_grouped(data_IgG1, x, show_points = FALSE, show_x = FALSE))

# MFI
IgG3_plots <- lapply(antigens_MFI, function(x) boxplot_grouped(data_IgG3, x, show_points = FALSE, show_x = FALSE))
C1q_plots <- lapply(antigens_MFI, function(x) boxplot_grouped(data_C1q, x, show_points = FALSE, show_x = FALSE))
FcgRI_plots <- lapply(antigens_MFI, function(x) boxplot_grouped(data_FcgRI, x, show_points = FALSE, show_x = FALSE))
FcgRII_plots <- lapply(antigens_MFI, function(x) boxplot_grouped(data_FcgRII, x, show_points = FALSE, show_x = FALSE))
FcgRIII_plots <- lapply(antigens_MFI, function(x) boxplot_grouped(data_FcgRIII, x, show_points = FALSE, show_x = FALSE))

# Multipanel plots
wrap_plots(IgG_plots, ncol = 4, nrow = 2)
wrap_plots(IgG1_plots, ncol = 4, nrow = 2)
wrap_plots(IgG3_plots, ncol = 4, nrow = 2)
wrap_plots(C1q_plots, ncol = 4, nrow = 2)
wrap_plots(FcgRI_plots, ncol = 4, nrow = 2)
wrap_plots(FcgRII_plots, ncol = 4, nrow = 2)
wrap_plots(FcgRIII_plots, ncol = 4, nrow = 2)
