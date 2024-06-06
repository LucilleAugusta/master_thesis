#!/usr/bin/env Rscript

# Load pQTLs from deCODE study
# Rename by alternative somascan
# Add EAF col

# Clear workspace ---------------------------------------------------------
rm(list = ls())

# Load libraries ----------------------------------------------------------
suppressPackageStartupMessages(library("tidyverse"))
suppressWarnings(library("TwoSampleMR"))
library("MRInstruments")
# library("ieugwasr")

# Load data ---------------------------------------------------------------
## Exposure data(deCODE)
exposure_file <- "data/01_decode_data/final_pqtls_PAV.csv"
exp_data <- read.table(exposure_file,
    header = TRUE
)

## Get coding regions from ref. genome
coding_region <- read.table(
    file = "data/resources/pc_gencodeGenes_filtered.csv"
)

# Load SomaScan gene names
translate_table <- read.table(
    file = "data/01_decode_data/03_translate_genename/translate_somascan_genes.csv",
    sep = "\t",
    header = TRUE
)

# Load assocvariants.annotated.txt
assocvariants_df <- read.table(
    file = "/home/lag/master_thesis/data/01_decode_data/assocvariants.annotated_filtered.txt",
    header = TRUE
)

# Wrangle data ------------------------------------------------------------


## Handle genes til alt. gene names
# Remove genes with alternative names
exp_data_temp <- exp_data %>%
    filter(!str_detect(
        file_name,
        "region2"
    ))

# Wrangle and replace 'gene_prot' with alternaitve gene name
translate_table <- translate_table %>%
    mutate(SeqId = str_replace(
        SeqId, "-", "_"
    ))

exp_data_alt <- exp_data %>%
    filter(str_detect(
        file_name,
        "region2"
    )) %>%
    left_join(., select(translate_table, gene_name, EntrezGeneName, SeqId),
        by = c(
            "gene_prot" = "gene_name",
            "SeqId" = "SeqId"
        )
    ) %>%
    select(!gene_prot) %>%
    rename(gene_prot = EntrezGeneName)

exp_data <- rbind(exp_data_temp, exp_data_alt)

# Join start and end pos from coding_region
exp_data <- coding_region %>%
    select(gene_name, start_pos, end_pos) %>%
    left_join(exp_data, .,
        by = c("gene_prot" = "gene_name")
    )

# exp_data_1 <- exp_data_1  %>%  mutate(
#     diff_1 = region_end - region_start,
#     diff_2 = end_pos - start_pos
#     ) # end is always bigger than start

# Make cis_snp label according to ref38 (coding_region)
exp_data <- exp_data %>% mutate(
    cis_snp = case_when(
        Pos >= start_pos - 1000000 &
            Pos <= end_pos + 1000000 ~ "TRUE",
        TRUE ~ "FALSE"
    )
)

# Count cis_snp
exp_data %>%
    group_by(cis_snp) %>%
    summarise(n())

# Investigate how "far" FALSE_snps are from cut-off

# exp_data %>%
#     filter(cis_snp == FALSE) %>%
#     select(
#         gene_prot, SeqId,
#         position, Pos,
#         region_start, region_end,
#         start_pos, end_pos
#     ) %>% mutate(diff = case_when(
#         Pos >= start_pos - 1000000 ~ start_pos - 1000000 - Pos,
#             Pos <= end_pos + 1000000 ~ Pos - 1000000 - end_pos))  %>%
#     head(20)


## Add EAF to exposure data
exp_data <- left_join(exp_data,
    assocvariants_df,
    by = "Name"
)

# exp_data <- exp_data %>% select(!rowname)

# exp_data <- read.table("data/01_decode_data/exp_data.csv",
#     header = TRUE
# )


# Make sure SNP with the lowest Pval is kept, when more are present
exp_data_distinct <- exp_data %>%
    group_by(SNP) %>%
    slice_min(Pval, n = 1) %>%
    ungroup()



# Write data ------------------------------------------------------------
write.table(exp_data,
    file = "data/01_decode_data/exp_data_PAV.csv",
    sep = "\t", row.names = FALSE, quote = FALSE
)

write.table(exp_data_distinct,
    file = "data/01_decode_data/exp_data_distinct_PAV.csv",
    sep = "\t", row.names = FALSE, quote = FALSE
)
