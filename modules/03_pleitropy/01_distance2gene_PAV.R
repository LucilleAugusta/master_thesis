#!/usr/bin/env Rscript

# Clear workspace ---------------------------------------------------------
rm(list = ls())


# Load libraries ----------------------------------------------------------
suppressPackageStartupMessages(library("tidyverse"))

# Load data ---------------------------------------------------------------

### Get input from command line
args <- commandArgs(trailingOnly = TRUE)
snp_file <- args[1]
print(snp_file)
# snp_file <- "data/01_decode_data/associated_SNPs/pleiotropy_clump/rs1064725.txt"

# Get gene list for x SNP
gene_list <- read.table(
    file = snp_file,
    col.names = c(
        "Chrom", "Pos", "Name",
        "rsids", "effectAllele", "otherAllele",
        "Beta", "Pval", "minus_log10_pval", "SE", "N", "ImpMAF"
    ),
    colClasses = c("character")
) %>%
    separate(
        col = Chrom,
        into = c("file_name", "Chr"),
        sep = ":"
    ) %>%
    select(
        !c(
            Beta, Pval, minus_log10_pval,
            SE, N, ImpMAF
        )
    )

# head(gene_list)
## Get coding regions from ref. genome
coding_region <- read.table(
    file = "data/resources/pc_gencodeGenes_filtered.csv"
)

# Get significant SNPs from MR
significant_snps <- read.csv(
    "data/02_analysis/01_mr/MR_results_single_significant_PAV.csv",
    sep = "\t",
    header = TRUE
) %>% select(SNP, gene.exposure)


# Wrangle data ------------------------------------------------------------
gene_list <- gene_list %>% mutate(
    seq_id    = str_extract(file_name, "\\d+_\\d+"),
    gene_name = str_extract(file_name, "\\d+_\\d+\\_?[a-zA-Z0-9]+"),
    gene_name = str_extract(gene_name, "[a-zA-Z0-9]+$")
)

# cis or trans?
gene_list <- coding_region %>%
    select(gene_name, start_pos, end_pos) %>%
    left_join(gene_list, .,
        by = "gene_name"
    ) %>%
    mutate(Pos = as.numeric(Pos))

gene_list <- gene_list %>% mutate(
    cis_gene_start = case_when(
        start_pos >= Pos - 1000000 &
            start_pos <= Pos + 1000000 ~ "TRUE",
        TRUE ~ "FALSE"
    ),
    cis_gene_end = case_when(
        end_pos >= Pos - 1000000 &
            end_pos <= Pos + 1000000 ~ "TRUE",
        TRUE ~ "FALSE"
    )
)

# Add gene for MR-hit
gene_list <- gene_list %>%
    left_join(., significant_snps,
        by = c("rsids" = "SNP")
    ) %>%
    rename("MR_hit" = "gene.exposure") %>%
    relocate(MR_hit, .after = gene_name)


# Write data ---------------------------------------------------------------
write.table(gene_list,
    file = snp_file,
    sep = "\t", row.names = FALSE, quote = FALSE
)

# gene_list  %>% select(file_name, Chr,Pos, gene_name, start_pos, end_pos, cis_gene_start, cis_gene_end)  %>% view()
