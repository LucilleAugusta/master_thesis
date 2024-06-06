#!/usr/bin/env Rscript


# Wrangle outcome from metal
# Filter out `+-` and `-+` cases
# Filter Pval>0.05

# Clear workspace ---------------------------------------------------------
rm(list = ls())

# Load libraries ----------------------------------------------------------
suppressPackageStartupMessages(library("tidyverse"))


# Load data ---------------------------------------------------------------
out_data <- read.table(
    file = "data/alzheimers_data/metal_AD1.txt",
    header = TRUE, sep = "\t"
)

out_data_WIDE <- read.table(
    file = "data/alzheimers_data/metal_AD_WIDE1.txt",
    header = TRUE, sep = "\t"
)

# Wrangle data ------------------------------------------------------------

out_data <- out_data %>% filter(
    P.value < 0.05,
    !Direction %in% c("-+", "+-")
)

out_data %>%
    group_by(Direction) %>%
    summarise(n())


out_data_WIDE <- out_data_WIDE %>% filter(
    P.value < 0.05,
    !Direction %in% c("-+", "+-")
)

# Write data ---------------------------------------------------------------
write.table(out_data,
    file = "data/alzheimers_data/out_data_AD.csv",
    sep = "\t", row.names = FALSE, quote = FALSE
)

write.table(out_data_WIDE,
    file = "data/alzheimers_data/out_data_AD_WIDE.csv",
    sep = "\t", row.names = FALSE, quote = FALSE
)
