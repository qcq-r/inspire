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

## Negative control VBDR datasets
neg_IgG <- transform(read.csv("data/VBDR_IgG_RAU.csv"),  assay = "IgG")
neg_IgG1 <- transform(read.csv("data/VBDR_IgG1_RAU.csv"),  assay = "IgG1")
neg_IgG3 <- transform(read.csv("data/VBDR_IgG3_MFI.csv"),  assay = "IgG3")
neg_C1q <- transform(read.csv("data/VBDR_C1q_MFI-Blank.csv"),  assay = "C1q")
neg_FcgRI <- transform(read.csv("data/VBDR_FcgRI_MFI-Blank.csv"),  assay = "FcgRI")
neg_FcgRII <- transform(read.csv("data/VBDR_FcgRII_MFI.csv"),  assay = "FcgRII")
neg_FcgRIII <- transform(read.csv("data/VBDR_FcgRIII_MFI.csv"),  assay = "FcgRIII")

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
                  "PvMSP10_RAU", "PvRBP2p1_RAU", "PvEXP1_RAU")
antigens_MFI <- c("Pv12_MFI", "Pv32_MFI", "Pv41_MFI", "Pv50_MFI", 
                  "PvMSP10_MFI", "PvRBP2p1_MFI", "PvEXP1_MFI")

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

# Include VBDR negative controls ------------------------------------------
boxplot_grouped <- function(data, var_name, neg_data = NULL, show_points = TRUE, show_x = TRUE, show_neg_ci = FALSE, ci_type = "median") {
  
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
  
  # 1. Extract negative control values if provided
  neg_median <- NULL
  neg_lower <- NULL
  neg_upper <- NULL
  
  if (!is.null(neg_data)) {
    neg_col <- paste0(antigen, "_", assay)
    if (neg_col %in% names(neg_data)) {
      neg_vals <- neg_data[[neg_col]]
      neg_vals <- neg_vals[!is.na(neg_vals) & neg_vals > 0] # Filter out missing/zeros
      if (length(neg_vals) > 0) {
        neg_median <- median(neg_vals)
        
        if (show_neg_ci) {
          if (ci_type == "iqr") {
            # Interquartile Range (25th to 75th percentile)
            quants <- quantile(neg_vals, probs = c(0.25, 0.75))
            neg_lower <- quants[1]
            neg_upper <- quants[2]
          } else if (ci_type == "median") {
            # Statistical 95% CI for the median (Non-parametric)
            tryCatch({
              wt <- wilcox.test(neg_vals, conf.int = TRUE, conf.level = 0.95)
              neg_lower <- wt$conf.int[1]
              neg_upper <- wt$conf.int[2]
            }, error = function(e) {
              neg_lower <- neg_median
              neg_upper <- neg_median
            })
          } else if (ci_type == "percentile") {
            # 95% Data Range (2.5th to 97.5th percentile)
            quants <- quantile(neg_vals, probs = c(0.025, 0.975))
            neg_lower <- quants[1]
            neg_upper <- quants[2]
          } else if (ci_type == "mean") {
            # Standard 95% CI for the Mean
            mean_val <- mean(neg_vals)
            se_val <- sd(neg_vals) / sqrt(length(neg_vals))
            margin <- qt(0.975, df = length(neg_vals) - 1) * se_val
            neg_lower <- mean_val - margin
            neg_upper <- mean_val + margin
          }
        }
      }
    }
  }
  
  # 2. Adjust y-axis limits safely to handle bounds without clipping
  yvals <- data[[var_name]]
  if (!is.null(neg_median)) {
    if (show_neg_ci && !is.null(neg_lower) && !is.null(neg_upper)) {
      yvals <- c(yvals, neg_lower, neg_upper)
    } else {
      yvals <- c(yvals, neg_median)
    }
  }
  yvals <- yvals[yvals > 0 & !is.na(yvals)]
  
  ymin <- 10^floor(log10(min(yvals)))
  ymax <- 10^ceiling(log10(max(yvals)))
  
  x_label <- if (show_x) "Group" else NULL
  
  y_label <- if (grepl("FcgR", assay)) {
    receptor <- sub("FcgR", "", assay)
    paste0(antigen, " Fc\u03B3R", receptor, " (", measure, ")")
  } else {
    paste0(antigen, " ", assay, " (", measure, ")")
  }
  
  # Base Boxplot
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
  
  # 3. Add Negative Control Reference Bands (Drawn underneath the boxplots)
  if (!is.null(neg_median)) {
    if (show_neg_ci && !is.null(neg_lower) && !is.null(neg_upper)) {
      p <- p + annotate("rect", xmin = -Inf, xmax = Inf, 
                        ymin = max(ymin, neg_lower), 
                        ymax = neg_upper, 
                        fill = "grey92", alpha = 0.7)
    }
    p <- p + geom_hline(yintercept = neg_median, linetype = "dotted", color = "grey40", linewidth = 0.8)
  }
  
  # 4. Toggle individual data points & Boxplots
  if (show_points) {
    p <- p +
      geom_boxplot(width = 0.5, staplewidth = 0.4, linewidth = 0.7, outlier.shape = NA) +
      geom_jitter(aes(colour = group), width = 0.15, alpha = 0.4, size = 2) +
      scale_colour_manual(values = c("G1" = "darkolivegreen", "G2" = "deepskyblue4", "G3" = "brown"))
  } else {
    p <- p +
      geom_boxplot(aes(colour = group), width = 0.5, staplewidth = 0.4, linewidth = 0.7) +
      scale_colour_manual(values = c("G1" = "darkgreen", "G2" = "royalblue4", "G3" = "red4"))
  }
  
  return(p)
}

IgG_plots <- lapply(antigens_RAU, function(x) {
  boxplot_grouped(
    data = data_IgG, 
    var_name = x, 
    neg_data = neg_IgG, 
    show_points = FALSE, 
    show_x = FALSE, 
    show_neg_ci = TRUE, 
    ci_type = "iqr")})

IgG1_plots <- lapply(antigens_RAU, function(x) {
  boxplot_grouped(
    data = data_IgG1, 
    var_name = x, 
    neg_data = neg_IgG1, 
    show_points = FALSE, 
    show_x = FALSE, 
    show_neg_ci = TRUE, 
    ci_type = "iqr")})

# MFI Layouts
IgG3_plots <- lapply(antigens_MFI, function(x) {
  boxplot_grouped(data_IgG3, 
                  var_name = x, 
                  neg_data = neg_IgG3, 
                  show_points = FALSE, 
                  show_x = FALSE, 
                  show_neg_ci = TRUE, 
                  ci_type = "iqr")})
C1q_plots <- lapply(antigens_MFI, function(x) {
  boxplot_grouped(data_C1q, 
                  var_name = x, 
                  neg_data = neg_C1q, 
                  show_points = FALSE, 
                  show_x = FALSE, 
                  show_neg_ci = TRUE, 
                  ci_type = "iqr")})
FcgRI_plots <- lapply(antigens_MFI, function(x) {
  boxplot_grouped(data_FcgRI, 
                  var_name = x, 
                  neg_data = neg_FcgRI, 
                  show_points = FALSE, 
                  show_x = FALSE, 
                  show_neg_ci = TRUE, 
                  ci_type = "iqr")})
FcgRII_plots <- lapply(antigens_MFI, function(x) {
  boxplot_grouped(data_FcgRII, 
                  var_name = x, 
                  neg_data = neg_FcgRII, 
                  show_points = FALSE, 
                  show_x = FALSE, 
                  show_neg_ci = TRUE, 
                  ci_type = "iqr")})
FcgRIII_plots <- lapply(antigens_MFI, function(x) {
  boxplot_grouped(data_FcgRIII, 
                  var_name = x, 
                  neg_data = neg_FcgRIII, 
                  show_points = FALSE, 
                  show_x = FALSE, 
                  show_neg_ci = TRUE, 
                  ci_type = "iqr")})

# Multipanel plots
wrap_plots(IgG_plots, ncol = 4, nrow = 2)
wrap_plots(IgG1_plots, ncol = 4, nrow = 2)
wrap_plots(IgG3_plots, ncol = 4, nrow = 2)
wrap_plots(C1q_plots, ncol = 4, nrow = 2)
wrap_plots(FcgRI_plots, ncol = 4, nrow = 2)
wrap_plots(FcgRII_plots, ncol = 4, nrow = 2)
wrap_plots(FcgRIII_plots, ncol = 4, nrow = 2)
