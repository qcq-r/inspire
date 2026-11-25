setwd("~/Downloads/data/inspire/R")
library(tidyverse)

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

# Import metadata
metadata <- read.csv("data/InSPIRE_metadata.csv")

# Filter metadata for IDs in dataset
metadata_filtered <- metadata %>%
  filter(xoa_id %in% data_collated$ID) %>%
  select(xoa_id, ageenyr, sex, village, pvivaxprev, avg_itn) %>%
  distinct(xoa_id, .keep_all = TRUE)

# Generate summary table
summary_table <- data.frame(
  Characteristic = c(
    "N",
    "Age, median (IQR)",
    "Male",
    "Female",
    "Sunuhu",
    "Ilaita",
    "P. vivax positive",
    "Average ITN use"
  ),
  Value = c(
    nrow(metadata_filtered),
    paste0(
      median(metadata_filtered$ageenyr, na.rm = TRUE),
      " (",
      quantile(metadata_filtered$ageenyr, 0.25, na.rm = TRUE),
      "-",
      quantile(metadata_filtered$ageenyr, 0.75, na.rm = TRUE),
      ")"
    ),
    paste0(sum(metadata_filtered$sex == "1"), " (",
           round(100 * mean(metadata_filtered$sex == "1"), 1), "%)"),
    paste0(sum(metadata_filtered$sex == "2"), " (",
           round(100 * mean(metadata_filtered$sex == "2"), 1), "%)"),
    paste0(sum(metadata_filtered$village %in% c("Sunuhu 1", "Sunuhu 2")), " (",
           round(100 * mean(metadata_filtered$village %in% c("Sunuhu 1", "Sunuhu 2")), 1), "%)"),
    paste0(sum(metadata_filtered$village %in% paste("Ilaita", 1:7)), " (",
           round(100 * mean(metadata_filtered$village %in% paste("Ilaita", 1:7)), 1), "%)"),
    paste0(sum(metadata_filtered$pvivaxprev == "1"), " (",
           round(100 * mean(metadata_filtered$pvivaxprev == "1"), 1), "%)"),
    paste0(round(mean(metadata_filtered$avg_itn), 3))
  )
)

