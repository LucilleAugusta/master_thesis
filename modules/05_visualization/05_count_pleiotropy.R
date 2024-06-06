#
# Create bar plot of count number of times each IV is in deCODE
#

# Clear workspace ---------------------------------------------------------
rm(list = ls())
# Load libraries ----------------------------------------------------------
suppressPackageStartupMessages(library("tidyverse"))
library(ggplot2)
library(data.table)
library(stringr)

# Load data ---------------------------------------------------------------

count_df <- fread("data/01_decode_data/associated_SNPs/pleiotropy_PAV/count_cis_gene.txt")


# Wrangle data ---------------------------------------------------------------
colnames(count_df) <- c("file_path", "gene_name", "total_count", "cis_count", "extra")

# Create regex to extract rs ID
pattern <- "rs[0-9]+"
count_df <- count_df %>% mutate(SNP = str_extract(file_path, pattern))
count_df <- count_df %>% mutate(name = paste0(gene_name, ", ", SNP))

# Plot data ---------------------------------------------------------------

ggplot(
    count_df,
    aes(
        x = reorder(gene_name, -cis_count),
        y = cis_count,
        fill = cis_count
    )
) +
    geom_bar(stat = "identity") +
    geom_hline(
        aes(yintercept = 5),
        col = "gray", linetype = "dashed"
    ) +
    # scale_x_discrete(guide = guide_axis(angle = -45)) +
    scale_fill_gradient(high = "#EC7200", low = "#EAC959") +
    labs(
        # = "cis-count"
        y = "Occurrence of SNP (+/- 1 MB)"
    ) +
    coord_flip() +
    theme_minimal() +
    theme(
        axis.title.y = element_blank(),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
    ) +
    ggsave(
        file = "./data/99_results/cis-count.png",
        units = "cm",
        width = 15,
        height = 10
    )
