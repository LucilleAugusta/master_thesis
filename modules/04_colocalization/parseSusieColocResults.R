#
#
# 1. load results from susieR
# 2. Extract H4
# 3. Save H4 in a matrix
# 4. load coloc resutls
# 5. Extract H4
# --------------------------------------------
#
# Susie
#

setwd("./modules/06_colocalization/")
library(tidyverse)
library(data.table)

susie.fns <- list.files(
    path = "results_2", pattern = "_susie_res.rds",
    full.names = TRUE, recursive = TRUE
)


trait_names <- unique(gsub("(\\w+)_\\w+_susie_res.rds", "\\1", basename(susie.fns)))
protein_names <- gsub("\\w+_(\\w+)_susie_res.rds", "\\1", basename(susie.fns))
protein_names <- unique(protein_names)


susie_h4max <- matrix(
    nrow = length(protein_names), ncol = length(trait_names),
    dimnames = list(protein_names, trait_names)
)
rownames(susie_h4max) <- protein_names
colnames(susie_h4max) <- trait_names

for (susie.fn in susie.fns) {
    # susie.fn <-  susie.fns[2]
    susie_res <- readRDS(susie.fn)

    trait <- gsub("(\\w+)_\\w+_susie_res.rds", "\\1", basename(susie.fn))
    protein <- gsub("\\w+_(\\w+)_susie_res.rds", "\\1", basename(susie.fn))

    # D1 <- susie_res[[2]]
    # D2 <- susie_res[[3]]
    susie_res <- susie_res[[1]]
    h4max <- susie_res$summary[["PP.H4.abf"]] %>% max()
    susie_h4max[protein, trait] <- h4max
}
susie_h4max


# write.table(
#     susie_h4max,
#     file = "results/susie_h4max_bakk.tsv",
#     sep = "\t",
#     quote = FALSE,
#     row.names = TRUE,
#     col.names = TRUE)


# --------------------------------------------
#
# Coloc
#


coloc.fns <- list.files(
    path = "results", pattern = "_coloc.csv",
    full.names = TRUE, recursive = TRUE
)


trait_names <- unique(gsub("(\\w+)_\\w+_coloc.csv", "\\1", basename(coloc.fns)))
protein_names <- gsub("\\w+_(\\w+)_coloc.csv", "\\1", basename(coloc.fns))
protein_names <- unique(protein_names)

coloc_h4max <- matrix(
    nrow = length(protein_names), ncol = length(trait_names),
    dimnames = list(protein_names, trait_names)
)
rownames(coloc_h4max) <- protein_names
colnames(coloc_h4max) <- trait_names

for (coloc.fn in coloc.fns) {
    # coloc.fn <-  coloc.fns[2]
    coloc_res <- fread(coloc.fn)

    trait <- gsub("(\\w+)_\\w+_coloc.csv", "\\1", basename(coloc.fn))
    protein <- gsub("\\w+_(\\w+)_coloc.csv", "\\1", basename(coloc.fn))

    h4max <- coloc_res$PP.H4.abf
    coloc_h4max[protein, trait] <- h4max
}
coloc_h4max

coloc_h4max <- as_tibble(coloc_h4max, rownames = "gene_name")

write.table(
    coloc_h4max,
    file = "./results/00_coloc_results.csv",
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
)
