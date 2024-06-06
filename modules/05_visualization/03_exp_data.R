#!/usr/bin/env Rscript

# Reads exp_data.csv and perform harmonisation and MR.

# Clear workspace ---------------------------------------------------------
rm(list = ls())

# Load libraries ----------------------------------------------------------
suppressPackageStartupMessages(library("tidyverse"))
suppressWarnings(library("TwoSampleMR"))
# library("MRInstruments")
# library("ieugwasr")
# library("xlsx")

# Load data ---------------------------------------------------------------
## Exposure data
# exposure_file <- "data/01_decode_data/exp_data_PAV.csv"
exposure_file <- "data/01_decode_data/exp_data_distinct_PAV.csv"

exp_data <- read_exposure_data(
    exposure_file,
    clump = FALSE,
    sep = "\t",
    # phenotype_col = "Phenotype",
    snp_col = "SNP",
    beta_col = "Beta",
    se_col = "SE",
    eaf_col = "effectAlleleFreq",
    effect_allele_col = "effectAllele",
    other_allele_col = "otherAllele",
    pval_col = "Pval",
    # units_col = "units",
    # ncase_col = "ncase",
    # ncontrol_col = "ncontrol",
    samplesize_col = "N",
    id_col = "SeqId",
    gene_col = "gene_prot",
    min_pval = 1e-200,
    log_pval = FALSE,
    chr_col = "Chrom",
    pos_col = "Pos"
) %>% filter(
    pval.exposure < 5e-8
)



## Outcome data

# outcome_file = "data/alzheimers_data/out_data_AD.csv"
outcome_file <- "data/alzheimers_data/out_data_AD_WIDE.csv"

out_data <- read_outcome_data(
    snps = exp_data$SNP,
    filename = outcome_file,
    sep = "\t",
    snp_col = "MarkerName",
    beta_col = "Effect",
    se_col = "StdErr",
    effect_allele_col = "Allele1",
    other_allele_col = "Allele2",
    # eaf_col = "a1_freq",
    pval_col = "P.value",
    # units_col = "Units",
    # gene_col = "Gene",
    # samplesize_col = "n"
)
dim(out_data) # 425 x 12


## Plot ----------------------------------------------------------
ggplot(
    exp_data,
    aes(
        x = gene.exposure,
        y = SNP
    )
) +
    geom_col()

exp_data %>%
    group_by(gene.exposure) %>%
    summarise(n = n()) %>% # filter(gene.exposure== "CD177") %>%
    ggplot(aes(x = gene.exposure, y = n)) +
    # geom_histogram(bins = 1000 )
    geom_col() +
    theme_minimal() +
    theme(axis.text.x = element_blank())

exp_data %>%
    group_by(gene.exposure) %>%
    summarise(n = n()) %>%
    arrange(desc(n))
