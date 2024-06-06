#!/usr/bin/env Rscript


# Description: plot results from coloc and SuSie.

# Clear workspace ---------------------------------------------------------
rm(list = ls())
setwd("./modules/06_colocalization/")
# Load libraries ----------------------------------------------------------
suppressPackageStartupMessages(library("tidyverse"))
library(ggplot2)
library(data.table)

# Load data ---------------------------------------------------------------

coloc_res <- fread("./results/00_coloc_results.csv")
view(coloc_res)
# Wrangle data ------------------------------------------------------------

coloc_res_longer <- coloc_res %>% pivot_longer(-gene_name,
    values_to = "H4",
    names_to = "phenotype"
)

coloc_res_longer <- coloc_res_longer %>% mutate(
    phenotype = case_when(
        phenotype == "Abeta" ~ "Cerebrospinal fluid amyloid beta levels",
        phenotype == "AD" ~ "Alzheimers",
        phenotype == "ptau" ~ "Cerebrospinal fluid p-tau levels",
        phenotype == "DE" ~ "Dementia",
        phenotype == "Lewy" ~ "Dementia with Lewy bodies",
        phenotype == "PD" ~ "Parkinson's disease"
    )
)


# Create intervals
coloc_res_longer <- coloc_res_longer %>%
    mutate(
        H4 = case_when(
            # H4 >= 1 ~ "+1",
            H4 >= 0.8 ~ "0.8 - 1",
            H4 >= 0.5 & H4 < 0.8 ~ "0.5 - 0.8",
            H4 >= 0.25 & H4 < 0.5 ~ "0.25 - 0.5"
            # H4 >= 1 ~ 1.5,
            # H4 >= 0.8 ~ 1,
            # H4 >= 0.5 & H4 < 0.8 ~ 0.8,
            # H4 < 0.5 ~ 0.3
        )
    ) %>%
    drop_na()


# PLot data ------------------------------------------------------------

coloc_res_longer %>% ggplot(aes(gene_name, phenotype,
    color = H4, size = H4
)) +
    geom_point(alpha = 0.8) +
    scale_colour_manual(values = c(
        "0.25 - 0.5" = "#2C7BB6",
        "0.5 - 0.8" = "#7FD07D",
        "0.8 - 1" = "#EC7200"
    )) +
    #   theme(legend.position="bottom") +
    scale_size_discrete() +
    scale_x_discrete(guide = guide_axis(angle = -45)) +
    guides(color = guide_legend(), size = guide_legend()) +
    # labs(guides = "Dose (mg)")+
    theme_minimal() +
    theme(axis.title = element_blank()) +
    ggsave(
        file = "../../data/99_results/coloc.png",
        units = "cm",
        width = 20,
        height = 10
    )

# Write data ---------------------------------------------------------------
