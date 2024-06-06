#
#
# 2. make out-of sample LD matrix for
#  each protein sumstat
#
# ------------------------------------

library(glue)
library(data.table)
library(tidyverse)

# p_fn <- "../../data/01_decode_data/coding_region/9574_11_BIN1_BIN1_region.csv"
# pqtl <- "rs35103166"
# protein_name <- gsub(".+_(.+)_region\\.csv$", "\\1", basename(p_fn))

# fread(p_fn, data.table = F) -> p_data

# check +- 1Mb
pqtl_position <- p_data %>%
    filter(rsids == pqtl) %>%
    pull(Pos)
lower_bound <- pqtl_position - 1000000
upper_bound <- pqtl_position + 1000000

# Filter the data frame
filtered_df <- p_data %>% filter(Pos >= lower_bound & Pos <= upper_bound)

rids_fn <- paste0("data/", protein_name, "_rsids.txt")
mat_fn <- paste0("data/", protein_name, "_ldmat.txt")


write.table(
    filtered_df$rsid,
    rids_fn,
    row.names = F,
    col.names = F,
    quote = F
)


# reference
g1k <- "../../data/resources/1kg.v3/EUR"

# Construct the command
cmd <- glue("plink --bfile {g1k} --cow --r2 square --extract {rids_fn} --write-snplist --out {mat_fn}")

# Run the command
system(cmd)



#################################################
# EOF # EOF # EOF # EOF # EOF # EOF # EOF # EOF #
#################################################
