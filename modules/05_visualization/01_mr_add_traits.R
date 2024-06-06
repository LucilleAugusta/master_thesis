#!/usr/bin/env Rscript

# Reads the significant SNPs and perform harmonisation and MR for xxx data.

# Clear workspace ---------------------------------------------------------
rm(list = ls())

# Load libraries ----------------------------------------------------------
suppressPackageStartupMessages(library("tidyverse"))
suppressWarnings(library("TwoSampleMR"))
library("MRInstruments")
library("data.table")
library("ggforce")
library("ieugwasr")

# Load data ---------------------------------------------------------------
## Exposure data ----------------------------------------------------------

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


# significant_snps <- read.csv(
#     "data/02_analysis/01_mr/mr_significant_meta_WIDE_clump.csv",
#     sep = "\t",
#     header = TRUE
# ) %>%
#     filter(SNP != "rs429358")

# exp_data <- exp_data %>% filter(SNP %in% significant_snps$SNP)

significant_snps <- fread(file = "data/02_analysis/01_mr/MR_results_single_PAV.csv")

# Bonferroni correction, (0.05/ n proteins)
n_gene <- significant_snps %>%
    distinct(gene.exposure) %>%
    dim()
bonferroni <- 0.05 / n_gene[1]


# filter significant genes
significant_snps <- significant_snps %>% filter(
    p < bonferroni & F_stat >= 10
)

data <- exp_data_clump %>% filter(SNP %in% significant_snps$SNP)

## Outcome data ---------------------------------------------------------------
# Check if API is avaiblable
# ieugwasr::api_status()
out_data <- extract_outcome_data(
    snps = data$SNP,
    outcomes = c(
        "ebi-a-GCST90029013", # Educational attainment (years of education)
        "ebi-a-GCST006572", # cognitive performance
        "ebi-a-GCST90001390", # Dementia with Lewy bodies
        "ebi-a-GCST90013868",
        "ebi-a-GCST90014023",
        "ebi-a-GCST90018910",
        "ebi-a-GCST009979",
        "ebi-a-GCST90018894",
        "ieu-b-110"
    )
)

## Finngen data
out_finngen_dementia <- read_outcome_data(
    snps = data$SNP,
    filename = "data/studied_phenotypes/finngen_R10_F5_DEMENTIA",
    sep = "\t",
    snp_col = "rsids",
    beta_col = "beta",
    se_col = "sebeta",
    effect_allele_col = "alt",
    other_allele_col = "ref",
    eaf_col = "af_alt",
    pval_col = "pval",
    chr_col = "#chrom",
    pos_col = "pos"
    # units_col = "Units",
    # gene_col = "Gene",
    # samplesize_col = "n"
) %>% mutate(
    outcome = "Dementia",
    originalname.outcome = "Dementia"
)
dim(out_finngen_dementia) # 14 x 15

out_finngen_FTD <- read_outcome_data(
    snps = data$SNP,
    filename = "data/studied_phenotypes/finngen_R10_FTD",
    sep = "\t",
    snp_col = "rsids",
    beta_col = "beta",
    se_col = "sebeta",
    effect_allele_col = "alt",
    other_allele_col = "ref",
    eaf_col = "af_alt",
    pval_col = "pval",
    chr_col = "#chrom",
    pos_col = "pos"
    # units_col = "Units",
    # gene_col = "Gene",
    # samplesize_col = "n"
) %>% mutate(
    outcome = "Frontotemporal dementia",
    originalname.outcome = "Frontotemporal dementia"
)

## CSF
### Cerebrospinal fluid amyloid beta 42 levels
CSF_amyloid <- read.table("data/studied_phenotypes/CSF_adding_rsids_2_sumstat/data_cerebrospinal_fluid_biomarkers/GCST90129599_buildGRCh38.harm.tsv",
    header = TRUE,
    sep = "\t"
)
CSF_amyloid <- CSF_amyloid %>% rename(
    chr = CHR,
    pos = BP,
    # effect_allele = A1,
    # other_allele = A2,
    effect_allele = A2,
    other_allele = A1,
    eaf = FRQ,
    pval = P
)

CSF_amyloid %>%
    group_by(N) %>%
    summarise(n())
# add beta and SE
## Beta = z / sqrt(2p(1− p)(n + z^2)) and
## SE =1 / sqrt(2p(1− p)(n + z^2))
## https://www-nature-com.proxy.findit.cvt.dk/articles/ng.3538

CSF_amyloid <- CSF_amyloid %>% mutate(
    beta = Z / sqrt(2 * eaf * (1 - eaf) * (N + Z^2)),
    se = 1 / sqrt(2 * eaf * (1 - eaf) * (N + Z^2))
)


out_CSF_amyloid <- format_data(
    CSF_amyloid,
    type = "outcome",
    snps = data$SNP,
    header = TRUE,
    snp_col = "SNP",
    beta_col = "beta",
    se_col = "se",
    eaf_col = "eaf",
    effect_allele_col = "effect_allele",
    other_allele_col = "other_allele",
    pval_col = "pval",
    chr_col = "chr",
    pos_col = "pos",
) %>% mutate(
    outcome = "Cerebrospinal fluid amyloid beta 42 levels",
    originalname.outcome = "Cerebrospinal fluid amyloid beta 42 levels"
)

### Cerebrospinal fluid p-tau levels

CSF_ptau <- read.table("data/studied_phenotypes/CSF_adding_rsids_2_sumstat/data_cerebrospinal_fluid_biomarkers/GCST90129600_buildGRCh38.harm.tsv",
    header = TRUE,
    sep = "\t"
)
CSF_ptau <- CSF_ptau %>% rename(
    chr = CHR,
    pos = BP,
    # effect_allele = A1,
    # other_allele = A2,
    effect_allele = A2,
    other_allele = A1,
    eaf = FRQ,
    pval = P
)

# add beta and SE
## Beta = z / sqrt(2p(1− p)(n + z^2)) and
## SE =1 / sqrt(2p(1− p)(n + z^2))
## https://www-nature-com.proxy.findit.cvt.dk/articles/ng.3538

CSF_ptau <- CSF_ptau %>% mutate(
    beta = Z / sqrt(2 * eaf * (1 - eaf) * (N + Z^2)),
    se = 1 / sqrt(2 * eaf * (1 - eaf) * (N + Z^2))
)

# CSF_ptau %>% filter(SNP %in% data$SNP)

out_CSF_ptau <- format_data(
    CSF_ptau,
    type = "outcome",
    snps = data$SNP,
    header = TRUE,
    snp_col = "SNP",
    beta_col = "beta",
    se_col = "se",
    eaf_col = "eaf",
    effect_allele_col = "effect_allele",
    other_allele_col = "other_allele",
    pval_col = "pval",
    chr_col = "chr",
    pos_col = "pos",
) %>% mutate(
    outcome = "Cerebrospinal fluid p-tau levels",
    originalname.outcome = "Cerebrospinal fluid p-tau levels"
)

# Wrangle data ------------------------------------------------------------
## Combine finngen + CSF w. other outcome phenotypes
out_finngen <- full_join(out_finngen_dementia, out_finngen_FTD)

out_CSF <- full_join(out_CSF_amyloid, out_CSF_ptau)

out_finngen_CSF <- full_join(out_finngen, out_CSF) %>%
    mutate(
        pos.outcome = as.character(pos.outcome),
        chr.outcome = as.character(chr.outcome)
    )

out_data_2 <- full_join(out_data, out_finngen_CSF,
    by = join_by(
        SNP, beta.outcome, se.outcome, pval.outcome,
        eaf.outcome, effect_allele.outcome, other_allele.outcome,
        outcome, id.outcome, originalname.outcome, mr_keep.outcome,
        data_source.outcome, pos == pos.outcome, chr == chr.outcome
    )
)


## Harmonise data
harm_data <- harmonise_data(
    exposure_dat = data,
    # outcome_dat = c(out_finngen_FTD, out_finngen_dementia),
    outcome_dat = out_data_2,
    action = 2 # Try to infer positive strand alleles, using allele frequencies for palindromes (default, conservative)
    # action = 3
)
dim(harm_data) # 211 x 45


## MR ------------------------------------------------------------

mr_res <- mr_singlesnp(
    harm_data,
    parameters = default_parameters(),
    # all_method = c("mr_ivw", "mr_ivw_fe", "mr_egger_regression", "mr_ivw_mre"),
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

# ## Create df for single SNP analysis.
res_single <- mr_res_full %>%
    filter(
        SNP != "All - Inverse variance weighted",
        SNP != "All - MR Egger"
        # mr_keep == "TRUE"
    )

## Plot ----------------------------------------------------------
### Create df to plot
significant_snps_reduced <- significant_snps %>%
    dplyr::select(id.exposure, gene.exposure, b, p) %>%
    mutate(originalname.outcome = "Alzheimers")

res_single_reduced <- res_single %>%
    dplyr::select(
        id.exposure, gene.exposure, b, p,
        originalname.outcome
    )
res_single_combined <- bind_rows(
    significant_snps_reduced,
    res_single_reduced
)
# create significant col
# p<0.05*, p<0.01**, p<0.001***
res_single_combined <- res_single_combined %>% mutate(stars = case_when(
    p < 0.001 ~ "***",
    p < 0.01 & p >= 0.001 ~ "**",
    p < 0.05 & p >= 0.01 ~ "*",
    p >= 0.05 ~ ""
))

# significant_snps %>%
#     group_by(chr.exposure, gene.exposure, id.exposure, SNP) %>%
#     summarise(n())

# Renaming
res_single_combined <- res_single_combined %>% mutate(
    originalname.outcome = case_when(
        originalname.outcome == "Coronary artery disease (SPA correction)" ~ "Coronary artery disease",
        TRUE ~ originalname.outcome
    )
)
unique(res_single_combined$originalname.outcome)

level_order <- c(
    # Dementia phenotypes
    "Alzheimers", "Dementia", "Parkinson's disease",
    "Frontotemporal dementia", "Dementia with Lewy bodies",
    "Cerebrospinal fluid amyloid beta 42 levels",
    "Cerebrospinal fluid p-tau levels",
    # Autoimmunity
    "Type 1 diabetes", "Coronary artery disease",
    "Rheumatoid arthritis",
    # Bad health
    "Major depressive disorder", "Educational attainment (years of education)",
    "LDL cholesterol", "Cognitive performance"
)

res_single_combined <- res_single_combined %>% mutate(
    phenotype_group = case_when(
        originalname.outcome == "Alzheimers" |
            originalname.outcome == "Dementia" |
            originalname.outcome == "Parkinson's disease" |
            originalname.outcome == "Frontotemporal dementia" |
            originalname.outcome == "Dementia with Lewy bodies" |
            originalname.outcome == "Cerebrospinal fluid amyloid beta 42 levels" |
            originalname.outcome == "Cerebrospinal fluid p-tau levels" ~ "C",
        # Autoimmunity
        originalname.outcome == "Type 1 diabetes" |
            originalname.outcome == "Coronary artery disease" |
            originalname.outcome == "Rheumatoid arthritis" ~ "B",
        # Bad health
        TRUE ~ "A"
    )
)


res_single_combined <- res_single_combined %>%
    filter(originalname.outcome != "Frontotemporal dementia")
# %>%
# group_by(originalname.outcome)  %>%
# summarise(n())


res_single_combined %>% ggplot(
    aes(gene.exposure, originalname.outcome,
        # factor(originalname.outcome, levels = level_order),
        fill = b
    )
) +
    geom_tile(color = "black") +
    scale_fill_gradientn(
        #values = c(1, .7, .68, .6, .58, 0), #including FTD
        values = c(1, .68, .66, .55, .5, 0), # TODO FIX COLORS !!!!
        colours = c(
            "#EC7200", "#EAC959",
            "white",
            "#a0bbdb", "#2C7BB6", "#2C7BB6"
        )
    ) +
    geom_text(aes(label = stars), color = "black", size = 2) +
    ggforce::facet_col(vars(phenotype_group),
        space = "free",
        scales = "free_y"
    ) +
    guides(fill = guide_colourbar(title = "Beta(MR)")) +
    scale_x_discrete(guide = guide_axis(angle = -45)) +
    theme_minimal() +
    theme(
        strip.text.x = element_blank(),
        panel.background = element_rect(fill = "grey", colour = "grey"),
        panel.grid.major = element_blank(),
        axis.title = element_blank()
    ) +
    ggsave(
        file = "./data/99_results/heatmap_MR_AD_related_phenotypes.png",
        units = "cm",
        width = 27,
        height = 15
    )
 
> res_single_combined %>% filter(originalname.outcome == "Cognitive performance")

res_single_combined %>%
    group_by(gene.exposure) %>%
    summarize(n_complication = sum(p < 0.05)) %>%
    ggplot(
        data = ., aes(
            x = reorder(gene.exposure, -n_complication),
            y = n_complication
        )
    ) +
    geom_bar(
        stat = "identity",
        fill = "#2C7BB6"
    ) +
    coord_flip() +
    labs(
        title = "Number of complications associated with AD-associated proteins",
        y = "Number of complications",
        x = "Protein name"
    ) +
    theme_minimal() +
    theme(
        plot.title = element_text(size = 14, face = "bold"),
        legend.position = "none"
    ) +
    ggsave(file = "./data/99_results/barplot_MR_AD_related_phenotypes.png")

# Write data ---------------------------------------------------------------
## Extract table of studied phenotypes.
ao <- available_outcomes()

out_ID <- out_data %>% distinct(id.outcome)

ao <- ao %>% filter(id %in% out_ID$id.outcome)

# write_excel_csv(ao, file = "./data/99_results/phenotype_table.csv")
