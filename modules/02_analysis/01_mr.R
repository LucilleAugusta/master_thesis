#!/usr/bin/env Rscript

# Reads exp_data.csv and perform harmonisation and MR.

# Clear workspace ---------------------------------------------------------
rm(list = ls())

# Load libraries ----------------------------------------------------------
suppressPackageStartupMessages(library("tidyverse"))
suppressWarnings(library("TwoSampleMR"))
library("MRInstruments")
library("ieugwasr")
library("xlsx")

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

# exp_data  %>%  distinct(gene.exposure)

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

# Wrangle data ------------------------------------------------------------

## Clump ------------------------------------------------------------
# Run clumping locally via. ieugwasr
exp_data_clump <- ld_clump(
    tibble(
        rsid = exp_data$SNP,
        pval = exp_data$pval.exposure,
        id = exp_data$exposure
    ),
    clump_p = 5e-8,
    clump_kb = 1000,
    clump_r2 = 0.01,
    plink_bin = genetics.binaRies::get_plink_binary(),
    # plink_bin = "/opt/intomics/software/bin/plink",
    bfile = "data/resources/1kg.v3/EUR"
)

# Join with initial exp_dat to keep all info
exp_data_clump <- inner_join(exp_data, exp_data_clump,
    by = c(
        "SNP" = "rsid",
        "pval.exposure" = "pval",
        "exposure" = "id"
    )
)


## Harmonise data ---------------------------------------------------------
harm_data <- harmonise_data(
    exposure_dat = exp_data_clump,
    outcome_dat = out_data,
    action = 2 # Try to infer positive strand alleles, using allele frequencies for palindromes (default, conservative)
    # action = 3
)
dim(harm_data) # 144 x 34
harm_data %>% distinct(gene.exposure)
# view(harm_data)

## MR ------------------------------------------------------------

mr_res <- mr_singlesnp(
    harm_data,
    parameters = default_parameters(),
    all_method = c("mr_ivw", "mr_ivw_fe", "mr_egger_regression", "mr_ivw_mre"),
    single_method = "mr_wald_ratio"
)


# Wrange data to keep gene name and harm info

mr_res_full <- left_join(mr_res, harm_data,
    by = c(
        "SNP",
        "id.outcome", "outcome",
        "id.exposure", "exposure"
    )
)

# Investigate n SNPs per. exposure
mr_res_full_2 <- mr_res_full %>%
    filter(!is.na(b))

# How many exposures are having +1 SNP?
mr_res_full_2 %>%
    filter(
        SNP == "All - Inverse variance weighted"
    ) %>%
    dim()

harm_data %>%
    filter(mr_keep == "TRUE") %>%
    group_by(id.exposure) %>%
    summarise(n()) %>%
    arrange() %>%
    view()
# 3600_2 (CHIT1) and 4535_50 (BST1) have to SNPs.

# colnames(mr_res_full)

## Calc -Log10(Pval)
# mr_res_full <- mr_res_full %>% mutate(
#     minusLog_pval = -log10(p)
# )

## Create df for single SNP analysis.
res_single <- mr_res_full %>%
    filter(
        SNP != "All - Inverse variance weighted",
        SNP != "All - MR Egger",
        SNP != "All - Inverse variance weighted (fixed effects)",
        SNP != "All - Inverse variance weighted (multiplicative random effects)",
        # samplesize != "Invalid Number"
        mr_keep == "TRUE"
    )

# view(res_single)

# res_single %>%
#     dplyr::select(
#         id.exposure, gene.exposure, p,
#         b, beta.exposure, beta.outcome
#     ) %>%
#     arrange(p) %>%
#     head(30)


## Delta method
# delta se
cov <- 0
res_single <- res_single %>%
    mutate(
        se.delta = sqrt((se.outcome^2 / beta.exposure^2)
        + (beta.outcome^2 / beta.exposure^4) * se.exposure^2
            - 2 * (beta.outcome / beta.exposure^3) * 0)
    ) %>%
    relocate(se.delta, .after = se)


# F stats (and place it after MR)
res_single <- res_single %>%
    group_by(id.exposure) %>%
    mutate(
        k = n(), # number of IV's / gene (k)
        variance_R2 = 2 * beta.exposure^2 * eaf.exposure * (1 - eaf.exposure),
        F_stat = (variance_R2 * (samplesize.exposure - 1 - k)) / ((1 - variance_R2) * k),
    ) %>%
    ungroup() %>%
    dplyr::select(!c(k, variance_R2))

dim(res_single)

res_single %>%
    dplyr::select(
        id.exposure, gene.exposure, SNP,
        p, b, F_stat
    ) %>%
    arrange(p) %>%
    head(20)

## Filter significant genes ------------------------------------------------------------------


# Bonferroni correction, (0.05/ n proteins)
n_gene <- res_single %>%
    distinct(gene.exposure) %>%
    dim()
bonferroni <- 0.05 / n_gene[1]
threshold <- bonferroni

# filter significant genes
res_single %>%
    filter(p < threshold &
        F_stat >= 10) %>%
    summarise(n())

res_single_sig <- res_single %>%
    filter(p < threshold &
        F_stat >= 10)
#     %>%
# select(gene.exposure)

# Get min and max beta + CI
# res_single <-fread("data/02_analysis/01_mr/MR_results_single_PAV.csv")
exp(res_single$b[112] + 1.96 * res_single$se.delta[112])

res_single <- res_single %>%
    mutate(
        CI_up_OR = exp(b + 1.96 * se.delta),
        CI_lo_OR = exp(b - 1.96 * se.delta)
    )

res_single <- res_single %>% mutate(
    OR = round(exp(b), 2),
    conf = paste0("[", round(CI_lo_OR, 2), " : ", round(CI_up_OR, 2), "]")
)




# Outcome and exposure alleles are the same - collapse cols
# res_single_sig %>% mutate(
#     TEST = as.character(other_allele.outcome == other_allele.exposure))  %>%
#     group_by(TEST)  %>%
#     summarise(n())

# Write data ---------------------------------------------------------------
## Write MR results
# write.table(res_single,
#     file = "data/02_analysis/01_mr/mr_results_meta_WIDE_single.csv",
#     sep = "\t", row.names = FALSE, quote = FALSE
# )

write.table(res_single,
    file = "data/02_analysis/01_mr/MR_results_single_PAV.csv",
    sep = "\t", row.names = FALSE, quote = FALSE
)
res_single_excel <- res_single %>%
    select(
        SNP, gene.exposure, chr.exposure, pos.exposure,
        b, p, se.delta, OR, conf,
        effect_allele.exposure, other_allele.exposure,
        beta.exposure, beta.outcome,
        pval.exposure, pval.outcome, F_stat
    )
write.xlsx(as.data.frame(res_single_excel),
    file = "data/99_results/tables/MR_results_single_PAV.xlsx",
    row.names = FALSE,
    col.names = TRUE
)

## Write significant SNPs
# write.table(res_single_sig,
#     file = "data/02_analysis/01_mr/mr_significant_meta_WIDE_clump.csv",
#     sep = "\t", row.names = FALSE, quote = FALSE
# )

# Write filtered table to excel
colnames(res_single_sig)
res_single_sig_excel <- res_single_sig %>%
    select(
        SNP, gene.exposure, chr.exposure, pos.exposure,
        b, p, se.delta,
        effect_allele.exposure, other_allele.exposure,
        beta.exposure, beta.outcome,
        pval.exposure, pval.outcome
    )

# res_single_sig_excel  %>%  view()
write.xlsx(as.data.frame(res_single_sig_excel),
    file = "data/99_results/tables/MR_results_single_sig_PAV.xlsx",
    row.names = FALSE,
    col.names = TRUE
)
