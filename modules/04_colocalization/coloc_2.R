#
#
# Read file for all outcome per "snp_list"

out_file <- paste0("./data/", protein_name, "_combined_outcome.csv")
out_data_full <- fread(out_file)

# Assign allele names to variables
allele1_col_name <- "effect_allele.outcome"
allele2_col_name <- "other_allele.outcome"
effect_col_name <- "beta.outcome"
rsid_col_name <- "SNP"
se_col_name <- "se.outcome"
# type_data <- "cc" # "cc or "quant"
# #trait_name <- ""

print(paste0("Preparing data for: ", protein_name))

# ---------------------------------------------------------------------------------
# p_fn <- "../../data/01_decode_data/coding_region/9574_11_BIN1_BIN1_region.csv"
# fread(p_fn, data.table = F) -> p_data

trait_list <- out_data_full %>%
    distinct(originalname.outcome) %>%
    # filter(originalname.outcome == "Alzheimers") %>%
    # filter(originalname.outcome == "Cerebrospinal fluid p-tau levels" |
    #     originalname.outcome == "Cerebrospinal fluid amyloid beta 42 levels") %>%
    pull(.)



for (i in trait_list) {
    # protein file should be read agian!
    fread(p_fn, data.table = F) -> p_data
    pqtl_snp <- significant_snps %>%
        filter(gene.exposure == protein_name) %>% pull(SNP)
    p_data %>% filter(rsids %in% pqtl_snp) %>% pull(Pos) -> pqtl_pos

    print(i)
    out_data <- out_data_full %>% filter(originalname.outcome == i)

    out_data <- out_data %>% rename(
        Allele1 = !!rlang::sym(allele1_col_name),
        Allele2 = !!rlang::sym(allele2_col_name),
        Effect = !!rlang::sym(effect_col_name),
        SNP = !!rlang::sym(rsid_col_name),
        SE = !!rlang::sym(se_col_name)
    )

    snplist <- fread(
        paste0("data/", protein_name, "_ldmat.txt.snplist"),
        header = F
    ) %>% pull(V1)

    out_data <- out_data %>% filter(SNP %in% snplist)


    p_data <- p_data %>% filter(rsids %in% out_data[, (SNP)])

    # remove duplicates
    out_data <- out_data %>%
        distinct(SNP, .keep_all = TRUE)
    p_data <- p_data %>%
        distinct(rsids, .keep_all = TRUE)

    # order
    ordering_vector <- match(out_data[["SNP"]], p_data$rsids)
    # Sort 'out_data' based on the ordering vector
    out_data <- out_data[order(ordering_vector), ]

    # check all equal
    print("check all equal: ")
    print(all(p_data$rsids == out_data[, (SNP)]))

    # Make alleles upper case
    out_data <- out_data %>%
        mutate(
            Allele1 := toupper(Allele1),
            Allele2 := toupper(Allele2)
        )

    # Filter out rows where alleles are not the same
    out_data <- out_data %>%
        filter((Allele1 == p_data$effectAllele &
            Allele2 == p_data$otherAllele) |
            (Allele2 == p_data$effectAllele &
                Allele1 == p_data$otherAllele))

    p_data <- p_data %>% filter(rsids %in% out_data[, (SNP)])

    print("check all equal rsids: ")
    print(all(p_data$rsids == out_data[, (SNP)]))


    # Flip alleles in 'out_data' if necessary
    out_data <- out_data %>%
        mutate(
            flip = ifelse((Allele1 == p_data$effectAllele & Allele2 == p_data$otherAllele), FALSE, TRUE),
            Effect = ifelse(flip, -Effect, Effect),
            temp_Allele1 = Allele1,
            temp_Allele2 = Allele2,
            Allele1 = ifelse(flip, temp_Allele2, Allele1),
            Allele2 = ifelse(flip, temp_Allele1, Allele2)
        ) %>%
        select(-temp_Allele1, -temp_Allele2, -flip)


    print("check all alleles equal: ")
    print(all(p_data$effectAllele == out_data$Allele1 & p_data$otherAllele == out_data$Allele2))

    # subset to +-  250 kb
    p_data <- p_data %>% filter(abs(Pos - pqtl_pos) < 250000)
    out_data <- out_data %>% filter(SNP %in% p_data[,"rsids"])
    #out_data <- out_data %>% filter(abs(pos - pqtl_pos) < 250000)
    
    print("check all equal rsids: ")
    print(all(p_data$rsids == out_data[, (SNP)]))


    # Load the LD mat
    ldmat <- fread(
        paste0("data/", protein_name, "_ldmat.txt.ld"),
        header = F,
        data.table = F
    )

    colnames(ldmat) <- snplist
    rownames(ldmat) <- snplist

    # select rsids
    ldmat <- ldmat[out_data$SNP, out_data$SNP]


    # Format for coloc -----------------------------------

    # Create df for protein
    D1 <- p_data %>%
        mutate(
            beta = as.numeric(Beta),
            varbeta = as.numeric(SE)^2,
            snp = rsids,
            position = Pos,
            MAF = ImpMAF
        )

    D1 <- list(
        beta = D1$beta,
        varbeta = D1$varbeta,
        snp = D1$snp,
        position = D1$position,
        type = "quant",
        N = mean(p_data$N),
        MAF = D1$MAF,
        sdY = 1,
        LD = as.matrix(ldmat)
    )

    # Create df for Outcome
    D2 <- out_data %>%
        mutate(
            beta = as.numeric(Effect),
            varbeta = as.numeric(SE)^2,
            snp = SNP,
            position = D1$position,
            MAF = eaf.outcome
        )

    if (out_data$data_type[1] == "cc") {
        D2 <- list(
            beta = D2$beta,
            varbeta = D2$varbeta,
            snp = D2$snp,
            position = D2$position,
            type = out_data$data_type[1],
            N = out_data$samplesize.outcome[1],
            LD = as.matrix(ldmat)
        )
    } else {
        D2 <- list(
            beta = D2$beta,
            varbeta = D2$varbeta,
            snp = D2$snp,
            position = D2$position,
            type = out_data$data_type[1],
            N = out_data$samplesize.outcome[1],
            LD = as.matrix(ldmat),
            MAF = D2$MAF,
            sdY = 1
        )
    }


    # ----------------------------------------------
    #
    # Coloc
    #
    #print("Doing COLOC")

    #my.res <- coloc.abf(
    #    dataset1 = D1,
    #    dataset2 = D2
    #)

    # Save  COLOC results ------------------------------------------

    # # Extract H4
    #H4 <- my.res$summary[["PP.H4.abf"]]
    #coloc_res <- my.res$summary
    #coloc_res["method"] <- "coloc"
    #coloc_res["prot_name"] <- protein_name
    #coloc_res["outcome"] <- out_data$originalname.outcome[1]
    #coloc_res_df <- data.frame(as.list(coloc_res))

    #write.table(coloc_res_df,
    #    file = paste0("./results/", out_data$trait_short[1], "_", protein_name, "_coloc.csv"),
    #    row.names = FALSE, sep = "\t", quote = FALSE
    #)

    # ----------------------------------------------
    #
    #     # SuSie
    #

    print(paste0("Running susie for: ", out_data$trait_short[1], "_", protein_name))
    susie.fn <- paste0("./results/", out_data$trait_short[1], "_", protein_name, "_susie_res.rds")

    S1 <- tryCatch(
         {
             runsusie(D1, maxit = 5000, repeat_until_convergence = FALSE) 
         },
         error = function(e) {
             message("Error in runsusie for D1: ", e$message)
             NULL
         }
     )

    S2 <- tryCatch(
         {
             runsusie(D2, maxit = 5000, repeat_until_convergence = FALSE)
         },
         error = function(e) {
             message("Error in runsusie for D2: ", e$message)
             NULL
         }
     )

    if (!is.null(S1) && !is.null(S2)) {
         susie.res <- coloc.susie(S1, S2)

        if(is.null(susie.res$nsnps)) {
             # Save results
            saveRDS(list(susie.res, D1, D2), file = susie.fn)
        } else {
            print("at least one dataset has no credible sets, nothing to colocalise")
        }        
     } else {
         msg <- paste(protein_name, "failed to run susie.")
         message(msg)
     }
}
# sensitivity(my.res, "H4 > 0.9")
