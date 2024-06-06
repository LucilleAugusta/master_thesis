# Clear workspace ---------------------------------------------------------
rm(list = ls())

# Colocalization analyses of AD associated genes
#
# *(a) install required packages
# 1. load data
# 2. make out-of sample LD matrix
# 3. run coloc & susieR
# -------------------------------
#
# Configs
#
# setwd("./modules/06_colocalization/")
# *a
# library(remotes)
# install_github("chr1swallace/coloc@main",build_vignettes=TRUE)
# install.packages("susieR")

library(coloc)
library(susieR)
library(data.table)
library(tidyverse)
library(optparse)

# if (!require(Rfast)) {
#   install.packages("Rfast")
# }
# library("Rfast")

# -------------------------------
#
# 1. Load data
#

# get prot_file
# Collect file list from coding regions
file_list_df <- list.files("../../data/01_decode_data/coding_region") %>%
  as_tibble_col(., column_name = "file_name")

# Extract gene name from file name
file_list_df <- file_list_df %>% mutate(
  seq_id    = str_extract(file_name, "\\d+_\\d+"),
  gene_name = str_extract(file_name, "\\d+_\\d+\\_?[a-zA-Z0-9]+"),
  gene_name = str_extract(gene_name, "[a-zA-Z0-9]+$")
)

# Define Snps to do for the analysis

significant_snps <- fread(file = "../../data/02_analysis/01_mr/MR_results_single_PAV.csv")

# Bonferroni correction, (0.05/ n proteins)
n_gene <- significant_snps %>%
  distinct(gene.exposure) %>%
  dim()
bonferroni <- 0.05 / n_gene[1]


# filter significant genes
significant_snps <- significant_snps %>%
  dplyr::filter(
    p < bonferroni & F_stat >= 10
  ) %>%
  dplyr::select(id.exposure, SNP, gene.exposure, chr.exposure, pos.exposure) # extract sample size

significant_snps <- significant_snps %>%
  dplyr::filter(gene.exposure == "PILRA")


for (i in significant_snps$id.exposure) {
  # -------------------------------
  #
  # 1. load data for i entry
  #

  # ad_fn <- "../../data/alzheimers_data/out_data_AD_WIDE.csv"
  # fread(ad_fn, data.table = F) -> ad_data

  # # Extract file name based on seq.id
  load_file_name <- file_list_df %>%
    dplyr::filter(seq_id == i) %>%
    pull(file_name)
  p_fn <- paste("../../data/01_decode_data/coding_region", load_file_name, sep = "/")
  print(p_fn)
  # fread(p_fn, data.table = F) -> p_data

  # Meta data
  pqtl <- significant_snps %>%
    dplyr::filter(id.exposure == i) %>%
    pull(SNP)
  protein_name <- significant_snps %>%
    dplyr::filter(id.exposure == i) %>%
    pull(gene.exposure)
  print(protein_name)

  # -------------------------------
  #
  # 2. make out-of sample LD matrix for
  #  each protein sumstat

  # source("make_ld_matrices.R")

  # -------------------------------
  #
  # 3. extract outcome data for 6 traits, based on the current cnp_list
  #
  print("Getting outcome data")
  source("coloc_add_traits.R")

  # -------------------------------
  #
  # 4. run coloc & susieR
  #
  print("Doing analysis ")
  # source("coloc_2.R")
  source("coloc_2_copy.R") # only AD and SUSIE
}


# -------------------------------
#
# 4. Save results
#
# source(parseSusieColocResults.R)
