## Introduction
## Author: Chau Quang (quang.c@wehi.edu.au)
## Purpose: To produce PCA plots of Taka's InSPIRE data

setwd("~/Downloads/data/inspire/R")

# Load required packages --------------------------------------------------
library(glmnet)
library(progress)
library(tidyverse)
library(car)

# Import data (.csv) ------------------------------------------------------
data_IgG  <- transform(read.csv("data/InSPIRE_IgG_grouped.csv"),  assay = "IgG")
data_IgG1 <- transform(read.csv("data/InSPIRE_IgG1_grouped.csv"), assay = "IgG1")
data_IgG3 <- transform(read.csv("data/InSPIRE_IgG3_grouped.csv"), assay = "IgG3")
data_C1q <- transform(read.csv("data/InSPIRE_C1q_grouped.csv"), assay = "C1q")
data_FcgRI <- transform(read.csv("data/InSPIRE_FcgRI_grouped.csv"), assay = "FcgRI")
data_FcgRII <- transform(read.csv("data/InSPIRE_FcgRII_grouped.csv"), assay = "FcgRII")
data_FcgRIII <- transform(read.csv("data/InSPIRE_FcgRIII_grouped.csv"), assay = "FcgRIII")

# Combine data for each antigen -----------------------------------------------
data_collated <- list(IgG   = data_IgG,
                      IgG1  = data_IgG1,
                      IgG3  = data_IgG3,
                      C1q   = data_C1q,
                      FcgRI = data_FcgRI,
                      FcgRII = data_FcgRII,
                      FcgRIII = data_FcgRIII)

data_collated <- imap(data_collated, function(df, assay_name) {
  df %>%
    select(-assay) %>%
    rename_with(function(col) {
      ifelse(col == "ID", col, paste0(col, "_", assay_name))
    })
}) %>%
  reduce(full_join, by = "ID")

# Function: Pre-process data for feature selection ----------------------------------
pca_preprocess <- function(data) {
  # Remove non-numeric data and convert to numeric matrix
  group_cols <- colnames(data)[grep("_Group", colnames(data))]
  feat_cols  <- colnames(data)[!colnames(data) %in% c("ID", group_cols)]
  data_matrix <- as.matrix(data[, feat_cols])
  
  # Right-shift to remove negative values
  min_vals <- apply(data_matrix, 2, min, na.rm = TRUE) 
  min_vals[min_vals > 0] <- 0
  shifted_data <- t(t(data_matrix) - min_vals)
  
  # log transform and normalise (min 0 - max 1)
  log_data <- log10(shifted_data + 1)
  norm_data <- scale(log_data)

  return(list(norm = norm_data, log = log_data))
}

# Function: LASSO feature selection -------------------------------------------------------
feature_selection <- function(norm_data, matrix_Y, tablelabel, iterations = 100) {
  set.seed(42)
  
  # Rename FcgR
  tablelabel <- gsub("FcgRI", "Fc\u03b3RI", tablelabel)
  tablelabel <- gsub("FcgRII", "Fc\u03b3RII", tablelabel)
  tablelabel <- gsub("FcgRIII", "Fc\u03b3RIII", tablelabel)
  
  all_selected_features <- c()
  num_features_per_run  <- c()
  
  # Progress bar
  pb <- progress_bar$new(
    format = "  Processing [:bar] :percent | ETA: :eta",
    total = iterations, 
    clear = FALSE, 
    width = 60
  )
  
  for (i in 1:iterations) {
    # 10-fold Cross Validation Multi-Response LASSO
    cv_fit <- cv.glmnet(norm_data, matrix_Y, nfolds = 10, alpha = 1, family = "mgaussian")
    
    # Extract coefficients
    coef_list <- coef(cv_fit, s = "lambda.min")
    
    selected_idx <- c()
    for (group_matrix in coef_list) {
      active_vars <- which(group_matrix[-1, 1] != 0)
      selected_idx <- union(selected_idx, active_vars)
    }
    
    all_selected_features <- c(all_selected_features, tablelabel[selected_idx])
    num_features_per_run  <- c(num_features_per_run, length(selected_idx))
    
    pb$tick()
  }
  
  cat("\nLoop complete! Generating plots...\n")
  
  # Tally distributions
  feature_counts <- table(all_selected_features)
  n_counts       <- table(num_features_per_run)
  
  # Plot 1: Features selected
  max_label_length <- max(nchar(names(feature_counts)))
  bottom_margin    <- max(5, max_label_length * 0.45)
  old_mar_p1       <- par(mar = c(bottom_margin, 4, 4, 2) + 0.1)
  
  barplot(sort(feature_counts), las = 2, col = "purple", main = "Features Selected", cex.names = 0.6)
  par(old_mar_p1)
  
  # Plot 2: # of features selected
  barplot(n_counts, col = "navy", main = "Number of Features Selected (Frequency)")
  
  number_features <- as.numeric(names(n_counts)[which.max(n_counts)])
  top_features    <- names(sort(feature_counts, decreasing = TRUE))[1:number_features]
  
  return(top_features)
}

# Function: Run PCA ---------------------------------------------------------
pca_analysis <- function(log_data, top_features, tablelabel, matrix_Y, ellipse = TRUE) {
  set.seed(42)
  # Rename FcgR
  top_features <- gsub("FcgRI", "Fc\u03b3RI", top_features)
  top_features <- gsub("FcgRII", "Fc\u03b3RII", top_features)
  top_features <- gsub("FcgRIII", "Fc\u03b3RIII", top_features)
  tablelabel   <- gsub("FcgRI", "Fc\u03b3RI", tablelabel)
  tablelabel   <- gsub("FcgRII", "Fc\u03b3RII", tablelabel)
  tablelabel   <- gsub("FcgRIII", "Fc\u03b3RIII", tablelabel)
  
  # Extract top features
  idx_final <- match(top_features, tablelabel)
  final_log_data_fs <- log_data[, idx_final, drop = FALSE]
  
  # Run PCA on scaled data
  pca_res <- prcomp(scale(final_log_data_fs), scale. = FALSE)
  
  # Calculate variance
  vars <- pca_res$sdev^2
  explained_var <- (vars / sum(vars)) * 100
  total_variance_2d <- sum(explained_var[1:2])
  cat("Total Variance (PC1+PC2):", round(total_variance_2d, 2), "%\n")
  
  # Average features for each sample and define groups
  avg_group <- round(rowMeans(matrix_Y))
  unique_classes <- sort(unique(avg_group))
  legend_labels <- paste0("G", unique_classes)
  
  master_palette <- c("1" = "#DD5129FF", "2" = "#FAB255FF", "3" = "#43B284FF")
  plot_colors <- master_palette[as.character(avg_group)]
  legend_colors <- master_palette[as.character(unique_classes)]
  
  # Colour antibody features
  ab_types <- sub(".*_", "", top_features)
  unique_ab_types <- sort(unique(ab_types))
  
  ab_palette_pool <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#A65628", "#F781BF", "#999999")
  ab_palette <- ab_palette_pool[1:length(unique_ab_types)]
  names(ab_palette) <- unique_ab_types
  feature_colors <- ab_palette[ab_types]
  
## Plot 1: PCA Plot
  plot(pca_res$x[,1], pca_res$x[,2], col = plot_colors, pch = 19, cex = 1.5,
       main = "PCA Scores",
       xlab = paste0("PC1 (", round(explained_var[1], 2), "%)"),
       ylab = paste0("PC2 (", round(explained_var[2], 2), "%)"))
  
  abline(h = 0, v = 0, lty = 3, col = "gray50")
  
  ### Toggle 95% confidence ellipse
  if (ellipse) {
    library(car)
    for (g in unique_classes) {
      idx <- which(avg_group == g)
      if (length(idx) >= 3) {
        dataEllipse(pca_res$x[idx, 1], pca_res$x[idx, 2], 
                    levels = 0.95, add = TRUE, plot.points = FALSE, 
                    col = master_palette[as.character(g)], 
                    fill = TRUE, fill.alpha = 0.15, 
                    lty = 2, lwd = 1, center.pch = FALSE)
      }
    }
  }
  legend("topleft", legend = legend_labels, col = legend_colors, pch = 19)
  
## Plot 2: PCA Loading Plot
  old_mar_p2 <- par(mar = c(5, 4, 4, 8) + 0.1)
  
  plot(pca_res$rotation[,1], pca_res$rotation[,2], col = feature_colors, pch = 16, cex = 1.2,
       main = "PCA Features (Loadings)",
       xlab = paste0("PC1 (", round(explained_var[1], 2), "%)"),
       ylab = paste0("PC2 (", round(explained_var[2], 2), "%)"))
  
  abline(h = 0, v = 0, lty = 3, col = "gray50")
  
  text(pca_res$rotation[,1], pca_res$rotation[,2], labels = top_features, 
       pos = 3, col = feature_colors, cex = 0.55)
  
  legend(x = par("usr")[2] + (par("usr")[2] - par("usr")[1])*0.04, 
         y = par("usr")[4], 
         legend = unique_ab_types, col = ab_palette, 
         pch = 16, bty = "n", cex = 0.8, xpd = TRUE)
  
  par(old_mar_p2)
  
## Plot 3: PCA Biplot
  x_lim <- range(pca_res$x[,1]) * 1.1
  y_lim <- range(pca_res$x[,2]) * 1.1
  
  plot(pca_res$x[,1], pca_res$x[,2], col = plot_colors, pch = 19, cex = 1.2,
       xlim = x_lim, ylim = y_lim, main = "Overlay Biplot",
       xlab = paste0("PC1 (", round(explained_var[1], 2), "%)"),
       ylab = paste0("PC2 (", round(explained_var[2], 2), "%)"))
  
  abline(h = 0, v = 0, lty = 3, col = "gray50")
  
  # Toggle ellipse
  if (ellipse) {
    for (g in unique_classes) {
      idx <- which(avg_group == g)
      if (length(idx) >= 3) {
        dataEllipse(pca_res$x[idx, 1], pca_res$x[idx, 2], 
                    levels = 0.95, add = TRUE, plot.points = FALSE, 
                    col = master_palette[as.character(g)], 
                    fill = TRUE, fill.alpha = 0.15, 
                    lty = 2, lwd = 1, center.pch = FALSE)
      }
    }
  }
  
  arrow_scale <- min(max(abs(x_lim)) / max(abs(pca_res$rotation[,1])), 
                     max(abs(y_lim)) / max(abs(pca_res$rotation[,2]))) * 0.7
  
  arrows(0, 0, pca_res$rotation[,1] * arrow_scale, pca_res$rotation[,2] * arrow_scale, 
         col = adjustcolor(feature_colors, alpha.f = 0.6), length = 0.07, lwd = 1)
  
  text(pca_res$rotation[,1] * arrow_scale * 1.1, pca_res$rotation[,2] * arrow_scale * 1.1, 
       labels = top_features, col = feature_colors, cex = 0.5)
  
## Plot 4: Horizontal PC1 loadings barplot
  pc1_loadings <- pca_res$rotation[, 1]
  names(pc1_loadings) <- top_features
  
  # Sort from highest positive (top of barplot) to lowest negative (bottom of barplot)
  sorted_indices <- order(pc1_loadings, decreasing = FALSE) 
  sorted_loadings <- pc1_loadings[sorted_indices]
  sorted_colors <- feature_colors[sorted_indices]
  
  # Dynamically adjust left margin based on your longest feature label name lengths
  max_label_length <- max(nchar(names(sorted_loadings)))
  left_margin_size <- max(4, max_label_length * 0.45) 
  
  # Expand margins: c(bottom, left, top, right)
  old_mar_p4 <- par(mar = c(5, left_margin_size, 4, 2) + 0.1)
  
  # Draw the horizontal barplot
  bp <- barplot(sorted_loadings, 
                horiz = TRUE, 
                las = 1,              # Make feature labels read horizontally
                col = sorted_colors,  # Match colors to antibody types
                border = NA,          # Remove messy block borders for ~50 bars
                cex.names = 0.45,     # Tiny text size so 50 labels fit perfectly
                main = "PC1 Feature Loadings (Contributions)",
                xlab = "Loading Weight")
  
  # Vertical dotted baseline at 0
  abline(v = 0, lty = 3, col = "gray30", lwd = 1.2)
  
  # Reset margins
  par(old_mar_p4)
}

# Load data into PCA processing pipeline ----------------------------------
## Identify antibody features from cleaned column names
colnames(data_collated) <- gsub("_RAU_|_MFI_", "_", colnames(data_collated))

## Assign grouping variable
group_cols <- colnames(data_collated)[grep("_Group", colnames(data_collated))]
## Assign antibody feature names
tablelabel <- colnames(data_collated)[!colnames(data_collated) %in% c("ID", group_cols)]

## Extract grouped data matrix
matrix_Y   <- as.matrix(data_collated[, group_cols])
## Retrieve log-transformed, normalised data
datasets   <- pca_preprocess(data_collated)

# Run feature selection
top_features <- feature_selection(datasets$norm, matrix_Y, tablelabel, iterations = 100)

# Generate final PCA plots
pca_analysis(datasets$log, top_features, tablelabel, matrix_Y, ellipse=TRUE)
