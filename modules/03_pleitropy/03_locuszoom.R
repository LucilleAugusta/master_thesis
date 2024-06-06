#!/usr/bin/env Rscript


# Description: Plotting Locuszoom for N significant genes.

# Clear workspace ---------------------------------------------------------
rm(list = ls())

# Load libraries ----------------------------------------------------------
suppressPackageStartupMessages(library("tidyverse"))
library(locuszoomr)
library("EnsDb.Hsapiens.v75") # not sure if this is correct
library(patchwork)
# Load data ---------------------------------------------------------------
file <- "data/01_decode_data/exp_data_distinct.csv"
file <- "data/01_decode_data/exp_data_distinct_PAV.csv"
exp_data <- read.table(file,
    header = TRUE
)
exp_data <- exp_data %>%
    dplyr::select(
        Chrom, Pos, SNP,
        otherAllele, effectAllele,
        Pval, Beta, SE, gene_prot
    ) %>%
    # dplyr::filter(gene_prot == "PSG5") %>%
    mutate(Chrom = str_remove(Chrom, "chr"))

# Outcome data
# outcome_file <- "data/alzheimers_data/out_data_AD_WIDE.csv"
# out_data <- read.table(outcome_file,
#     header = TRUE
# )
# Load Kunlke to get chrom
# kunkle <- read.table(
#     "data/alzheimers_data/finngen_R10_G6_AD_WIDE",
#     header = FALSE, sep = "\t"
# )

# kunkle_2 <- kunkle[,c(1,2,5)]
# colnames(kunkle_2) <- c("Chrom", "pos", "rsids")

# left_join <- left_join(out_data, kunkle_2,
#     by = c("MarkerName" = "rsids")
# )

# write.table(left_join,
#     file = "data/alzheimers_data/out_data_AD_WIDE_chr.csv",
#     sep = "\t", row.names = FALSE, quote = FALSE
# )
# out_data <- left_join

out_data <- read.table("data/alzheimers_data/out_data_AD_WIDE_chr.csv",
    header = TRUE
)

out_data <- out_data %>% rename(MarkerName = "SNP")

# Get significant genes on chr 19.
# res_single <- read.csv(
#     file = "data/02_analysis/01_mr/MR_results_single_PAV.csv",
#     header = TRUE, sep = "\t"
# )

# # Bonferroni correction, (0.05/ n proteins)
# n_gene <- res_single %>%
#     distinct(gene.exposure) %>%
#     dim()
# bonferroni <- 0.05 / n_gene[1]

# # filter significant genes
# res_single <- res_single %>%
#     dplyr::filter(p < bonferroni &
#         F_stat >= 10 &
#         chr.exposure == "chr19")




# Wrangle data ------------------------------------------------------------
# DO for gene
# my_loc <- locus(
#     data = data, gene = "PSG5", flank = 1e6,
#     ens_db = "EnsDb.Hsapiens.v75"
# )
# summary(my_loc)
# head(my_loc)
# locus_plot(my_loc)


# Plot for SNPs ------------------------------------------------------------
SNP_hit <- "rs1859788" # significant SNP for PILRA

pos_hit <- exp_data %>%
    dplyr::filter(SNP == SNP_hit) %>%
    dplyr::slice(1:1) %>%
    pull(Pos)

gene_hit <- exp_data %>%
    dplyr::filter(SNP == SNP_hit) %>%
    dplyr::slice(1:1) %>%
    pull(gene_prot)
range <- 1000000

# Select SNPs +/- 1Mb from the hit
# data_SNP <- exp_data %>% dplyr::filter(
#     Pos >= pos_hit - range & Pos <= pos_hit + range
# ) # why are all SNPs in the dataset not plotted?
# data_SNP %>% dplyr::filter( Chrom == "19")
# the function is only plotting those present on the same chr...

# Token for LDlink API
# LAG_token <- "ae43a7ce032a"
LAG_token <- "22a7b8338892"

### EXPOSURE ###
gene_loc_exp <- locus(
    data = exp_data, gene = gene_hit,
    # index_snp = SNP_hit,
    flank = range,
    ens_db = "EnsDb.Hsapiens.v75"
)
gene_loc_exp <- link_LD(gene_loc_exp,
    token = LAG_token
)

snp_loc_exp <- locus(
    data = exp_data,
    index_snp = SNP_hit, flank = range,
    ens_db = "EnsDb.Hsapiens.v75"
)
snp_loc_exp <- link_LD(snp_loc_exp,
    token = LAG_token
)
### OUTCOME ###
gene_loc_out <- locus(
    data = out_data, p = "P.value",
    gene = gene_hit,
    flank = range,
    ens_db = "EnsDb.Hsapiens.v75"
)
gene_loc_out <- link_LD(gene_loc_out,
    token = LAG_token
)
p_gene_loc_out <- locus_ggplot(gene_loc_out,
    labels = SNP_hit,
    highlight = gene_hit
)

p_gene_loc_out
# ggtitle( "Out data") +
ggsave(file = "./locus_gene_AD.png")



#############
p_snp_loc_out <- locus_ggplot(snp_loc_out,
    labels = SNP_hit,
    highlight = gene_hit
)


### GENERATE PLOTS ###

# summary(snp_loc_exp)
# head(snp_loc_exp)

p_snp_loc_exp <- locus_ggplot(snp_loc_exp,
    labels = SNP_hit,
    highlight = gene_hit
)
# %>%
#     ggsave(file = "./locus_snp.png")

p_gene_loc_exp <- locus_ggplot(gene_loc_exp,
    labels = SNP_hit,
    highlight = gene_hit
)
# %>%
#     ggsave(file = "./locus_gene.png")


# Write data ---------------------------------------------------------------
