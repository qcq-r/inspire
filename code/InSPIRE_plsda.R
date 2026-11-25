## Introduction
## Author: Chau Quang (quang.c@wehi.edu.au)
## Purpose: Experimental code to produce PLS-DA plots using pv500 infection from Taka's InSPIRE data
## NOTE: Code runs but is depracated due to low classification (BER)

setwd("~/Downloads/data/inspire/R")

# Load required packages --------------------------------------------------
library(tidyverse)
library(mixOmics)

library(glmnet)
library(progress)

# Import data (.csv) ------------------------------------------------------
data_IgG  <- transform(read.csv("data/InSPIRE_IgG_grouped.csv"),  assay = "IgG")
data_IgG1 <- transform(read.csv("data/InSPIRE_IgG1_grouped.csv"), assay = "IgG1")
data_IgG3 <- transform(read.csv("data/InSPIRE_IgG3_grouped.csv"), assay = "IgG3")
data_C1q <- transform(read.csv("data/InSPIRE_C1q_grouped.csv"), assay = "C1q")
data_FcgRI <- transform(read.csv("data/InSPIRE_FcgRI_grouped.csv"), assay = "FcgRI")
data_FcgRII <- transform(read.csv("data/InSPIRE_FcgRII_grouped.csv"), assay = "FcgRII")
data_FcgRIII <- transform(read.csv("data/InSPIRE_FcgRIII_grouped.csv"), assay = "FcgRIII")

metadata <- read.csv("data/InSPIRE_metadata.csv") %>%
  rename("ID" = "xoa_id")

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
    dplyr::select(-assay) %>%
    rename_with(function(col) {
      ifelse(col == "ID", col, paste0(col, "_", assay_name))
    })
}) %>%
  reduce(full_join, by = "ID")
colnames(data_collated) <- gsub("_RAU_|_MFI_", "_", colnames(data_collated))
data_collated_metadata <- merge(metadata, data_collated, by = "ID")

# Retrieve Pv infection (>500 parasites/uL) status per sample
infection_status <- data_collated_metadata %>%
  group_by(ID) %>%
  summarise(
    Max_Infection = if(all(is.na(pv_500))) NA else max(pv_500, na.rm = TRUE)
  ) %>%
  filter(!is.na(Max_Infection)) %>% 
  mutate(Ever_Infected = factor(ifelse(Max_Infection >= 1, "Infected", "Uninfected"))) %>%
  arrange(ID)

# Define functions for pre-processing and LASSO feature selection -------------------------------------------------------
pca_preprocess <- function(data) {
  # Remove non-numeric data and convert to numeric matrix
  group_cols <- colnames(data)[grep("_Group", colnames(data))]
  feat_cols  <- colnames(data)[!colnames(data) %in% c("ID", group_cols)]
  data_matrix <- as.matrix(data[, feat_cols])
  
  # Right-shift to remove negative values
  min_vals <- apply(data_matrix, 2, min, na.rm = TRUE) 
  min_vals[min_vals > 0] <- 0
  shifted_data <- t(t(data_matrix) - min_vals)
  
  # log transform and z-score
  log_data <- log10(shifted_data + 1)
  norm_data <- scale(log_data)
  
  return(list(norm = norm_data, log = log_data))
}

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

# Load data into PCA processing pipeline ----------------------------------
## Log-transform and normalise data
datasets <- pca_preprocess(data_collated)

## Extract X normalised data matrix
X_matrix <- datasets$norm
colnames(X_matrix) <- gsub("_RAU_|_MFI_", "_", colnames(X_matrix))

## Extract Y outcome matrix
Y_binary_pls <- infection_status$Ever_Infected
Y_binary_pls <- droplevels(Y_binary_pls)

## Assign tertiles as grouping variable
group_cols <- colnames(data_collated)[grep("_Group", colnames(data_collated))]

## Assign antibody feature names
tablelabel <- colnames(data_collated)[!colnames(data_collated) %in% c("ID", group_cols)]

## Run feature selection
top_features <- feature_selection(datasets$norm, matrix_Y, tablelabel, iterations = 100)

## Subset X matrix to LASSO-selected features only
X_lasso <- X_matrix[, colnames(X_matrix) %in% top_features, drop = FALSE]

## Run PLS-DA on selected features vs infection outcome
plsda_res <- plsda(X_lasso, Y_binary_pls, ncomp = 2)
print(plsda_res)

# Plots --------------------------------------------
# Plot 1: PLSDA scores
plotIndiv(plsda_res, 
          comp = c(1,2),
          group = Y_binary_pls,
          ind.names = FALSE,
          ellipse = TRUE,        
          legend = TRUE,
          col = c("#DD5129FF", "#43B284FF"), 
          title = "PLS-DA")

# Plot 2: Antibody feature contribution to PC1 or PC2
plotLoadings(plsda_res, comp = 1, 
             method = 'median', 
             contrib = 'max', 
             size.name = 0.6,
             title = "Contribution to PC1")

plotLoadings(plsda_res, comp = 2, 
             method = 'median', 
             contrib = 'max', 
             size.name = 0.6,
             title = "Contribution to PC2")

# 8. Run 10-fold cross validation -------------------------------------
perf_plsda <- perf(plsda_res, validation = "Mfold", folds = 10, nrepeat = 10)

# View classification error rates
print(perf_plsda$error.rate)

#> print(perf_plsda$error.rate)
#$overall
#max.dist centroids.dist mahalanobis.dist
#comp1 0.3967568      0.4324324        0.4324324
#comp2 0.3924324      0.4362162        0.4681081

#$BER
#max.dist centroids.dist mahalanobis.dist
#comp1 0.4965938      0.4487368        0.4487368
#comp2 0.4814731      0.4581951        0.4843326