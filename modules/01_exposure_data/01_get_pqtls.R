#!/usr/bin/env Rscript


# Find pQTLs from deCODE study and extract them from download files.

# Clear workspace ---------------------------------------------------------
rm(list = ls())

# Load libraries ----------------------------------------------------------
suppressPackageStartupMessages(library("tidyverse"))
library(readxl)

# Load data ---------------------------------------------------------------
decode_data <- read_excel(
    "data/alzheimers_data/41588_2021_978_MOESM4_ESM.xlsx",
    sheet = 2,
    skip = 2
)

# Wrangle data ------------------------------------------------------------

## Wrangle colnames for decode pQTLs
decode_data <- decode_data %>% select(
    "pQTL_ID \r\n(global)",
    "region_ID \r\n(prot.)",
    "gene\r\n (prot.)",
    "chr\r\n(prot.)",
    "SeqId",
    "variant",
    "LD \r\nclass\r\nsize", "chr\r\n(var.)",
    "pos\r\n(var.)",
    "region\r\nstart",
    "region\r\nend",
    "Amin",
    "Amaj",
    "MAF\r\n(%)",
    "cis/\r\ntrans",
    "Any\r\nPAV\r\nsame\r\ngene"
)

decode_data <- decode_data %>% rename(
    pQTLID_global = "pQTL_ID \r\n(global)",
    regionID      = "region_ID \r\n(prot.)",
    gene_prot     = "gene\r\n (prot.)",
    chr_prot      = "chr\r\n(prot.)",
    SNP           = "variant",
    LD_class_size = "LD \r\nclass\r\nsize",
    chr_var       = "chr\r\n(var.)",
    position      = "pos\r\n(var.)",
    region_start  = "region\r\nstart",
    region_end    = "region\r\nend",
    MAF           = "MAF\r\n(%)",
    type          = "cis/\r\ntrans",
    PAV           = "Any\r\nPAV\r\nsame\r\ngene"
)

## Filter data

decode_data_filter <- decode_data %>% filter(
    type == "cis"
)

# Collect file list from coding regions
file_list_df <- list.files("data/01_decode_data/coding_region") %>%
    as_tibble_col(., column_name = "file_name")

# Extract gene name from file name
file_list_df <- file_list_df %>% mutate(
    seq_id    = str_extract(file_name, "\\d+_\\d+"),
    gene_name = str_extract(file_name, "\\d+_\\d+\\_?[a-zA-Z0-9]+"),
    gene_name = str_extract(gene_name, "[a-zA-Z0-9]+$")
)

# Check if gene exists
# (gene_list$gene_prot %in% file_list_df$gene_name)


# For each gene in decode_excel corresponding file is opened and rsIDs is saved
gene_list <- decode_data_filter %>% distinct(gene_prot)
final_table <- tibble()

for (i in gene_list$gene_prot) {
    # If gene are in decode gene list, all files with that gene is opened.
    if ((i %in% file_list_df$gene_name)) {
        file_entry <- file_list_df %>% filter(gene_name == i)

        for (j in 1:nrow(file_entry)) {
            # Extract current file info
            current_file <- file_entry[j, ]

            # Load current file
            gene_file <- read.table(
                paste("data/01_decode_data/coding_region",
                    current_file$file_name,
                    sep = "/"
                ),
                header = TRUE
            )

            # Save seq_id and file_name
            gene_file <- gene_file %>%
                add_column(
                    file_name = current_file$file_name,
                    seq_id = current_file$seq_id
                )

            # Get gene rows from decode list
            gene_decode <- decode_data_filter %>% filter(gene_prot == i)

            # match those variants
            gene_complete <- inner_join(gene_decode, gene_file,
                by = c(
                    "SNP" = "rsids",
                    "SeqId" = "seq_id"
                )
            )

            # save in df
            final_table <- rbind(final_table, gene_complete)
        }
    }
}

# view(final_table)

# Write data ---------------------------------------------------------------
write.table(final_table,
    file = "/home/lag/master_thesis/data/01_decode_data/final_pqtls_PAV.csv",
    sep = "\t",
    quote = FALSE
)
