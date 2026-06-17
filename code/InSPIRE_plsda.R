## Introduction
## Author: Chau Quang (quang.c@wehi.edu.au)
## Purpose: Experimental code to produce PLS-DA plots using pv500 infection from Taka's InSPIRE data
## NOTE: Code runs but is depracated due to low classification (BER)

setwd("~/Downloads/data/inspire/R")

# Load required packages --------------------------------------------------
library(tidyverse)
library(mixOmics)

# Retrieve Pv infection (>500 parasites/uL) status per sample
infection_status <- data_collated_metadata %>%
  group_by(ID) %>%
  summarise(
    Max_Infection = if(all(is.na(pv_500))) NA else max(pv_500, na.rm = TRUE)
  ) %>%
  filter(!is.na(Max_Infection)) %>% 
  mutate(Ever_Infected = factor(ifelse(Max_Infection >= 1, "Infected", "Uninfected"))) %>%
  arrange(ID)

# Log-transform and normalise data
datasets_pls <- pca_preprocess(data_collated)

# Extract X normalised data matrix
X_matrix <- datasets_pls$norm
colnames(X_matrix) <- gsub("_RAU_|_MFI_", "_", colnames(X_matrix))

# Extract Y outcome matrix
Y_binary_pls <- infection_status$Ever_Infected
Y_binary_pls <- droplevels(Y_binary_pls)

# Diagnostics
cat("Number of patient rows in X:", nrow(X_matrix), "\n")
cat("Number of target entries in Y:", length(X_matrix), "\n")
print(table(Y_binary_pls, useNA = "ifany")) 

# Subset the matrix down to your robust LASSO features from PCA code
X_lasso <- X_matrix[, colnames(X_matrix) %in% top_features, drop = FALSE]

plsda_res <- plsda(X_lasso, Y_binary, ncomp = 2)
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
