# Collect outcome in one df per snp_list
# ----------------------------------------------
#
# Format data
#
suppressWarnings(library("TwoSampleMR"))


snplist <- fread(
    paste0("data/", protein_name, "_ldmat.txt.snplist"),
    header = F
) %>% pull(V1)

ad_file <- "../../data/alzheimers_data/out_data_AD_WIDE.csv"
ad_out_data <- read_outcome_data(
    snps = snplist,
    filename = ad_file,
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
ad_out_data <- ad_out_data %>% mutate(
    outcome = "Alzheimers",
    originalname.outcome = "Alzheimers",
    data_type = "cc",
    samplesize.outcome = 476107
)

demens_out_data <- extract_outcome_data(
    snps = snplist,
    outcomes = c(
        "ebi-a-GCST90001390", # Dementia with Lewy bodies
        "ebi-a-GCST90018894" # Parkinson
    )
) %>% mutate(
    data_type = "cc"
)

## Finngen data
out_finngen_dementia <- read_outcome_data(
    snps = snplist,
    filename = "../../data/studied_phenotypes/finngen_R10_F5_DEMENTIA",
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
)
out_finngen_dementia <- out_finngen_dementia %>% mutate(
    outcome = "Dementia",
    originalname.outcome = "Dementia",
    data_type = "cc", ,
    samplesize.outcome = 407717
)


## CSF
### Cerebrospinal fluid amyloid beta 42 levels
CSF_amyloid <- read.table("../../data/studied_phenotypes/CSF_adding_rsids_2_sumstat/data_cerebrospinal_fluid_biomarkers/GCST90129599_buildGRCh38.harm.tsv",
    header = TRUE,
    sep = "\t"
)
CSF_amyloid <- CSF_amyloid %>% dplyr::rename(
    chr = CHR,
    pos = BP,
    effect_allele = A1,
    other_allele = A2,
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
    snps = snplist,
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
)
out_CSF_amyloid <- out_CSF_amyloid %>% mutate(
    outcome = "Cerebrospinal fluid amyloid beta 42 levels",
    originalname.outcome = "Cerebrospinal fluid amyloid beta 42 levels",
    data_type = "quant",
    samplesize.outcome = 8074
)

### Cerebrospinal fluid p-tau levels

CSF_ptau <- read.table("../../data/studied_phenotypes/CSF_adding_rsids_2_sumstat/data_cerebrospinal_fluid_biomarkers/GCST90129600_buildGRCh38.harm.tsv",
    header = TRUE,
    sep = "\t"
)
CSF_ptau <- CSF_ptau %>% dplyr::rename(
    chr = CHR,
    pos = BP,
    effect_allele = A1,
    other_allele = A2,
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
    snps = snplist,
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
)

out_CSF_ptau <- out_CSF_ptau %>% mutate(
    outcome = "Cerebrospinal fluid p-tau levels",
    originalname.outcome = "Cerebrospinal fluid p-tau levels",
    data_type = "quant",
    samplesize.outcome = 7798
)

# Wrangle data ------------------------------------------------------------
## Combine finngen + CSF w. other outcome phenotypes
out_finngen <- full_join(out_finngen_dementia, ad_out_data)

out_CSF <- full_join(out_CSF_amyloid, out_CSF_ptau)

out_finngen_CSF <- full_join(out_finngen, out_CSF) %>%
    mutate(
        pos.outcome = as.character(pos.outcome),
        chr.outcome = as.character(chr.outcome)
    )

out_data_2 <- full_join(demens_out_data, out_finngen_CSF,
    by = join_by(
        SNP, beta.outcome, se.outcome, pval.outcome, data_type,
        eaf.outcome, effect_allele.outcome, other_allele.outcome,
        outcome, id.outcome, originalname.outcome, mr_keep.outcome,
        data_source.outcome, samplesize.outcome,
        pos == pos.outcome, chr == chr.outcome
    )
)

# out_data_2 %>%
#     group_by(originalname.outcome) %>%
#     summarise(n())

# length(snplist)

# Add traits ID
out_data_2 <- out_data_2 %>% mutate(
    trait_short = case_when(
        originalname.outcome == "Alzheimers" ~ "AD",
        originalname.outcome == "Cerebrospinal fluid amyloid beta 42 levels" ~ "Abeta",
        originalname.outcome == "Cerebrospinal fluid p-tau levels" ~ "ptau",
        originalname.outcome == "Dementia" ~ "DE",
        originalname.outcome == "Dementia with Lewy bodies" ~ "Lewy",
        originalname.outcome == "Parkinson's disease" ~ "PD"
    )
)

# out_data_2 %>%
#     group_by(samplesize.outcome) %>%
#     summarise(n())

write.table(out_data_2,
    file = paste0("./data/", protein_name, "_combined_outcome.csv"),
    sep = "\t", row.names = FALSE, quote = FALSE
)
