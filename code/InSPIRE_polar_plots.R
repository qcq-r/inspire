## Introduction
## Author: Chau Quang (quang.c@wehi.edu.au)
## Purpose: To produce polar (circular) plots showing the relative magnitude of antibody features of Taka's InSPIRE data

setwd("~/Downloads/data/inspire/R")

# Load required packages --------------------------------------------------
library(ggplot2)
library(tidyverse)
library(datawizard)

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

clean_df <- function(df, prefix) {
  df %>%
    select(ID, starts_with(prefix)) %>%
    {
      keep_main <- select(.,
                          matches(paste0("^", prefix, "_RAU_IgG$")),
                          matches(paste0("^", prefix, "_RAU_IgG1$")),
                          matches(paste0("^", prefix, "_MFI_IgG3$")),
                          matches(paste0("^", prefix, "_MFI_C1q$")),
                          matches(paste0("^", prefix, "_MFI_FcgRI$")),
                          matches(paste0("^", prefix, "_MFI_FcgRII$")),
                          matches(paste0("^", prefix, "_MFI_FcgRIII$"))
      ) %>%
        setNames(c(
          "IgG", "IgG1", "IgG3", "C1q",
          "FcgRI", "FcgRII", "FcgRIII"
        ))
      
      keep_group <- select(.,
                           matches(paste0("^", prefix, "_Group_"))
      ) %>%
        setNames(c(
          "IgG_Group", "IgG1_Group", "IgG3_Group",
          "C1q_Group", "FcgRI_Group", "FcgRII_Group",
          "FcgRIII_Group"
        ))
      
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

# Normalise antibody data -------------------------------------------------
lognorm <- function(data) {
  
  logdata <- log(data[, 2:8])
  
  #Check variables included
  logdata <- logdata[, sapply(logdata, function(x) length(unique(x)) > 1)]
  
  vars <- c("IgG", "IgG1", "IgG3", "C1q", "FcgRI", "FcgRII", "FcgRIII")
  vars <- intersect(vars, names(logdata))
  
  setdiff(c("IgG", "IgG1", "IgG3", "C1q", "FcgRI", "FcgRII", "FcgRIII"),
          names(logdata)) |>
    (\(x) if(length(x)) message(paste("Warning:", x, "was removed")))()
  
  logdata %>%
    as.data.frame() %>%
    normalize(select = vars) %>%
    round(3) %>%
    mutate(ID = data$ID) %>%
    relocate(ID, .before = everything())
}

df_Pv41_norm <- lognorm(df_Pv41)
df_Pv50_norm <- lognorm(df_Pv50)
df_Pv12_norm <- lognorm(df_Pv12)
df_PvMSP10_norm <- lognorm(df_PvMSP10)
df_Pv32_norm <- lognorm(df_Pv32)
df_PvRBP2p1_norm <- lognorm(df_PvRBP2p1)
df_PvEXP1_norm <- lognorm(df_PvEXP1)

norm_med <- function(data, id_col = "ID") {
  
  data %>%
    select(-all_of(id_col)) %>%
    summarise(across(everything(), ~ median(.x, na.rm = TRUE))) %>%
    pivot_longer(
      cols = everything(),
      names_to = "Feature",
      values_to = "Y"
    )
}

Pv41_norm_med <- norm_med(df_Pv41_norm)
Pv50_norm_med <- norm_med(df_Pv50_norm)
Pv12_norm_med <- norm_med(df_Pv12_norm)
PvMSP10_norm_med <- norm_med(df_PvMSP10_norm)
Pv32_norm_med <- norm_med(df_Pv32_norm)
PvRBP2p1_norm_med <- norm_med(df_PvRBP2p1_norm)
PvEXP1_norm_med <- norm_med(df_PvEXP1_norm)

# Polar plots -------------------------------------------------------------
polar_plot <- function(data) {
  antigen <- sub("_norm_med$", "", deparse(substitute(data)))
  data$Feature <- gsub("FcgRI",  "Fc\u03B3RI",  data$Feature)
  data$Feature <- gsub("FcgRII", "Fc\u03B3RII", data$Feature)
  data$Feature <- gsub("FcgRIII","Fc\u03B3RIII",data$Feature)
  
  polarplot <- ggplot(data, aes(Feature, Y)) +
    geom_col(aes(fill = Feature), position = position_dodge(width = 1), width = 1)
  polarplot + coord_polar() +
    ylim(0, 1) +
    ggtitle(antigen) +
    theme_bw() +
    theme(axis.text.x = element_text(size=13, face="bold"),
          axis.text.y = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line = element_line(size = 0),
          legend.position = "none",
          plot.title = element_text(size = 20, face = "bold", hjust = 0.5)) +
    annotate('text', x = 3, y = c(0.25,0.5, 0.75, 1.0), label = c('0.25', '0.5', '0.75', '1.0'), size = 4)
  
}

# Export > width 581 height 590, use device pixel ratio
polar_plot(Pv41_norm_med)
polar_plot(Pv50_norm_med)
polar_plot(Pv12_norm_med)
polar_plot(PvMSP10_norm_med)
polar_plot(Pv32_norm_med)
polar_plot(PvRBP2p1_norm_med)
polar_plot(PvEXP1_norm_med)
