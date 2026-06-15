## Introduction
## Author: Chau Quang (quang.c@wehi.edu.au)
## Purpose: To produce correlation heatmaps between Ab features of Taka's InSPIRE data

setwd("~/Downloads/data/inspire/R")

# Load required packages --------------------------------------------------
library(ggplot2)
library(corrplot)
library(dplyr)
library(purrr)
library(tidyr)

# Import data (.csv) ------------------------------------------------------
data_IgG  <- transform(read.csv("data/stata/InSPIRE_IgG_grouped.csv"),  assay = "IgG")
data_IgG1 <- transform(read.csv("data/stata/InSPIRE_IgG1_grouped.csv"), assay = "IgG1")
data_IgG3 <- transform(read.csv("data/stata/InSPIRE_IgG3_grouped.csv"), assay = "IgG3")
data_C1q <- transform(read.csv("data/stata/InSPIRE_C1q_grouped.csv"), assay = "C1q")
data_FcgRI <- transform(read.csv("data/stata/InSPIRE_FcgRI_grouped.csv"), assay = "FcgRI")
data_FcgRII <- transform(read.csv("data/stata/InSPIRE_FcgRII_grouped.csv"), assay = "FcgRII")
data_FcgRIII <- transform(read.csv("data/stata/InSPIRE_FcgRIII_grouped.csv"), assay = "FcgRIII")

# Combine antibody data for each antigen ----------------------------------
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

clean_df <- function(df, prefix) {
  df %>%
    select(ID, starts_with(prefix)) %>%
    {
      ## select columns
      keep_main <- select(.,
                          matches(paste0("^", prefix, "_RAU_IgG$")),
                          matches(paste0("^", prefix, "_RAU_IgG1$")),
                          matches(paste0("^", prefix, "_MFI_IgG3$")),
                          matches(paste0("^", prefix, "_MFI_C1q$")),
                          matches(paste0("^", prefix, "_MFI_FcgRI$")),
                          matches(paste0("^", prefix, "_MFI_FcgRII$")),
                          matches(paste0("^", prefix, "_MFI_FcgRIII$"))
      )
      
      ## reorder columns
      keep_group <- select(., contains("_Group_"))
      
      bind_cols(
        select(., ID),
        keep_main,
        keep_group
      )
    }
}

# Filtered datasets
df_Pv41 <- clean_df(data_collated, "Pv41")
df_Pv50 <- clean_df(data_collated, "Pv50")
df_Pv12 <- clean_df(data_collated, "Pv12")
df_PvMSP10 <- clean_df(data_collated, "PvMSP10")
df_Pv32 <- clean_df(data_collated, "Pv32")
df_PvRBP2p1 <- clean_df(data_collated, "PvRBP2p1")
df_PvEXP1 <- clean_df(data_collated, "PvEXP1")

# Plot Spearman's correlation heatmap ---------------------------------
plot_corr_heatmap <- function(df, method = "spearman", show_sig = TRUE) {
  
  ## prep matrix
  df_clean <- df %>%
    dplyr::select(-ID, -dplyr::contains("_Group_")) %>%
    dplyr::rename_with(~ sub(".*_(IgG1|IgG3|IgG|C1q|FcgRI|FcgRII|FcgRIII)$", "\\1", .))
  
  mat <- as.matrix(df_clean)
  
  ## remove variables with zero variance
  sds <- apply(mat, 2, sd, na.rm = TRUE)
  zero_var <- sds == 0 | is.na(sds)
  
  removed_val <- names(sds[zero_var])
  if (length(removed_val) > 0) {
    message("Removed: ", paste(removed_val, collapse = ", "))
  }
  
  mat <- mat[, !zero_var, drop = FALSE]
  
  ## rename FcRs
  colnames(mat) <- gsub("FcgRI",  "Fc\u03B3RI",  colnames(mat))
  colnames(mat) <- gsub("FcgRII", "Fc\u03B3RII", colnames(mat))
  colnames(mat) <- gsub("FcgRIII","Fc\u03B3RIII",colnames(mat))
  
  ## calculate correlation matrix
  cor_mat <- cor(mat, use = "pairwise.complete.obs", method = method)
  
  ## calculate p-value matrix
  n <- ncol(mat)
  
  p_mat <- matrix(NA, n, n)
  colnames(p_mat) <- colnames(mat)
  rownames(p_mat) <- colnames(mat)
  
  for (i in seq_len(n)) {
    for (j in i:n) {
      
      x <- mat[, i]
      y <- mat[, j]
      
      if (sd(x, na.rm=TRUE)==0 || sd(y, na.rm=TRUE)==0) next
      
      pval <- suppressWarnings(cor.test(x, y, method=method)$p.value)
      
      p_mat[i, j] <- pval
      p_mat[j, i] <- pval
    }
  }
  
  ## remove diagonals
  diag(cor_mat) <- NA
  
  ## colour scale (r=-1 to r=1)
  col_fun <- colorRampPalette(c("#053061", "#2166AC", "#4393C3", "#92C5DE","#F7F7F7","#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F"))(200)
  
  ## correlation heatmap
  corrplot::corrplot(
    cor_mat,
    method = "color",
    col = col_fun,
    type = "lower",
    diag = FALSE,
    
    # Ab feature labels
    font = 2,
    tl.srt = 0,
    tl.col = "black",
    tl.cex = 0.9,
    tl.offset = 0.9,
    
    # show coefficients
    addCoef.col = "black",
    
    # significance
    p.mat = if (show_sig) p_mat else NULL,
    sig.level = c(0.001, 0.01, 0.05),
    insig = "label_sig",
    pch.cex = 1.5,
    pch.col = "black"
  )
  
  ## colour legend title
  mtext("Spearman's correlation coefficient (r)", 
        side = 1, line = 3.5, cex = 0.9)
  
  ## plot title
  title(main = bquote(bold(.(unique(sub("_.*", "", colnames(df)[2]))))), 
        line = 3)
}

# Purpose: To manually adjust y-axis position of significance asterisks.
# Instructions: Run trace(corrplot, edit=TRUE) once at the start of the session as below. To plot correctly call plot_corr_heatmap() again.
trace(corrplot, edit=TRUE) # @Line 446 copy& paste: text(pos.pNew[, 1][sig.locs], pos.pNew[, 2][sig.locs]-0.3

# List of correlation heatmaps -------------------------------------------------------
## Export > Save Image > PNG; Width 609 px, Height 565 px; Use device pixel ratio
plot_corr_heatmap(df_Pv41)
plot_corr_heatmap(df_Pv50)
plot_corr_heatmap(df_Pv12)
plot_corr_heatmap(df_PvMSP10)
plot_corr_heatmap(df_Pv32)
plot_corr_heatmap(df_PvRBP2p1)
plot_corr_heatmap(df_PvEXP1)

# Heatmap by antibody feature ---------------------------------------------
plot_corr_heatmap_feature <- function(df, feature = "IgG",
                                      method = "spearman",
                                      show_sig = TRUE) {
  
  ## select columns for chosen feature
  suffix <- switch(
    feature,
    "IgG"     = "_RAU_IgG$",
    "IgG1"    = "_RAU_IgG1$",
    "IgG3"    = "_MFI_IgG3$",
    "C1q"     = "_MFI_C1q$",
    "FcgRI"   = "_MFI_FcgRI$",
    "FcgRII"  = "_MFI_FcgRII$",
    "FcgRIII" = "_MFI_FcgRIII$"
  )
  
  df_clean <- df %>%
    select(ID, matches(suffix))
  
  ## rename columns to antigen names
  antigen_names <- names(df_clean)[-1]
  
  antigen_names <- antigen_names |>
    sub("_RAU_IgG$", "", x = _) |>
    sub("_RAU_IgG1$", "", x = _) |>
    sub("_MFI_IgG3$", "", x = _) |>
    sub("_MFI_C1q$", "", x = _) |>
    sub("_MFI_FcgRI$", "", x = _) |>
    sub("_MFI_FcgRII$", "", x = _) |>
    sub("_MFI_FcgRIII$", "", x = _)
  
  names(df_clean)[-1] <- antigen_names
  
  mat <- as.matrix(select(df_clean, -ID))
  
  ## remove zero variance antigens
  sds <- apply(mat, 2, sd, na.rm = TRUE)
  
  mat <- mat[, sds > 0 & !is.na(sds), drop = FALSE]
  
  ## calculate correlation matrix
  cor_mat <- cor(mat, use = "pairwise.complete.obs", method = method)
  
  ## calculate p-value matrix
  n <- ncol(mat)
  
  p_mat <- matrix(NA, n, n)
  colnames(p_mat) <- colnames(mat)
  rownames(p_mat) <- colnames(mat)
  
  for (i in seq_len(n)) {
    for (j in i:n) {
      
      x <- mat[, i]
      y <- mat[, j]
      
      if (sd(x, na.rm=TRUE)==0 || sd(y, na.rm=TRUE)==0) next
      
      pval <- suppressWarnings(cor.test(x, y, method=method)$p.value)
      
      p_mat[i, j] <- pval
      p_mat[j, i] <- pval
    }
  }
  
  ## remove diagonals
  diag(cor_mat) <- NA
  
  ## colour scale (r=-1 to r=1)
  col_fun <- colorRampPalette(c("#053061", "#2166AC", "#4393C3", "#92C5DE","#F7F7F7","#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F"))(200)
  
  ## correlation heatmap
  corrplot::corrplot(
    cor_mat,
    method = "color",
    col = col_fun,
    type = "lower",
    diag = FALSE,
    
    # Ab feature labels
    font = 2,
    tl.srt = 45,
    tl.col = "black",
    tl.cex = 0.8,
    tl.offset = 0.5,
    
    # show coefficients
    addCoef.col = "black",
    
    # significance
    p.mat = if (show_sig) p_mat else NULL,
    sig.level = c(0.001, 0.01, 0.05),
    insig = "label_sig",
    pch.cex = 1.5,
    pch.col = "black"
  )
  
  ## colour legend title
  mtext("Spearman's correlation coefficient (r)", 
        side = 1, line = 3.5, cex = 0.9)
  
  ## plot title
  ## feature title
  feature_title <- feature
  
  feature_title <- gsub("FcgRIII", "Fc\u03B3RIII", feature_title)
  feature_title <- gsub("FcgRII",  "Fc\u03B3RII",  feature_title)
  feature_title <- gsub("FcgRI",   "Fc\u03B3RI",   feature_title)
  
  title(
    main = bquote(bold(.(feature_title))),
    line = 1
  )
}

plot_corr_heatmap_feature(data_collated, "IgG")
plot_corr_heatmap_feature(data_collated, "IgG1")
plot_corr_heatmap_feature(data_collated, "IgG3")
plot_corr_heatmap_feature(data_collated, "C1q")
plot_corr_heatmap_feature(data_collated, "FcgRI")
plot_corr_heatmap_feature(data_collated, "FcgRII")
plot_corr_heatmap_feature(data_collated, "FcgRIII")

save.image(file = "inspire.RData")
