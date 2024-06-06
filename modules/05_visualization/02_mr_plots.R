#!/usr/bin/env Rscript


# Plots snps with p < bonferroni & F_stat >= 10

# Clear workspace ---------------------------------------------------------
rm(list = ls())

# Load libraries ----------------------------------------------------------
suppressPackageStartupMessages(library("tidyverse"))
library(ggforestplot)
library(ggplot2)
library(TwoSampleMR)
library(patchwork)


# Load data ---------------------------------------------------------------

res_single <- read.csv(
    # file = "data/02_analysis/mr_results_single.csv",
    file = "data/02_analysis/01_mr/MR_results_single_PAV.csv",
    header = TRUE, sep = "\t"
)


# Wrangle data ------------------------------------------------------------

# Bonferroni correction, (0.05/ n proteins)
n_gene <- res_single %>%
    distinct(gene.exposure) %>%
    dim()
bonferroni <- 0.05 / n_gene[1]



# filter significant genes
# res_single_sig <- res_single %>%
#     filter(p < bonferroni &
#         F_stat >= 10) 

# write.table(res_single_sig,
#     file = "data/02_analysis/01_mr/MR_results_single_significant_PAV.csv",
#     row.names = FALSE, sep = "\t", quote = FALSE
# )


# Plot --------------------------------------------------------------------

## Volcano plot
# bonferroni <- 0.005
# Mark up- and down regulated
res_single <- res_single %>% mutate(
    association = case_when(
        b > 0 & p < bonferroni & F_stat >= 10 ~ "UP",
        b < 0 & p < bonferroni & F_stat >= 10 ~ "DOWN",
        TRUE ~ "NEUTRAL"
    ),
    significant_genes = case_when(
        p < bonferroni &
            F_stat >= 10 ~ gene.exposure,
    ),
    minusLog10_pval = -log10(p)
)

# Mark the SNPs on chr19 to be excluded
# res_single  %>%  filter( chr.exposure == "chr19")
# excluded_snps <- res_single %>%
#     filter(p < bonferroni &
#         F_stat >= 10 &
#         chr.exposure == "chr19") %>%
#     pull(SNP)


res_single %>%
    group_by(association) %>%
    summarise(n())

### VOLCANO PLOT ###
# theme_set(theme_minimal() + theme(text = element_text(family = "Helvetica")))
# Estimate ORs and corresponding CIs by Wald ratio and the delta method
res_single$myFacet <- res_single$minusLog10_pval < 10


p_volcano <- ggplot(data = res_single, aes(
    x = b, y = minusLog10_pval,
    color = association, label = significant_genes
)) +
    geom_vline(xintercept = c(0), col = "#484848", linetype = "dashed") +
    geom_point(size = 2) +
    scale_color_manual(
        name = "Association",
        values = c(
            "UP" = "#EC7200",
            "DOWN" = "#32A29B",
            # "EXCLUDE" = "#474554",
            "NEUTRAL" = "#A3A1A8"
        ),
        labels = c(
            "DOWN" = "Inverse",
            "UP" = "Positive",
            # "EXCLUDE" = "Chr19",
            "NEUTRAL" = "Neutral"
        )
    ) +
    ggrepel::geom_text_repel(
        data = . %>%
            mutate(label = case_when(
                chr.exposure != "chr19" & p < bonferroni ~ significant_genes,
                chr.exposure == "chr19" & p < bonferroni ~ paste0(significant_genes, "*"),
                TRUE ~ ""
            )),
        aes(label = label),
        min.segment.length = 0,
        max.overlaps = Inf,
        # data = subset(res_single, chr.exposure != "chr19"),
        box.padding = unit(0.35, "lines"),
        show.legend = FALSE
    ) +
    geom_hline(
        data = subset(res_single, myFacet == "TRUE"),
        aes(yintercept = c(-log10(bonferroni))),
        col = "gray", linetype = "dashed"
    ) +
    facet_wrap(. ~ myFacet,
        nrow = 2,
        scales = "free_y"
    ) +
    ylab("-Log10Pval") +
    xlab(expression(beta)) +
    theme_minimal() +
    theme(
        plot.title = element_text(size = 14, face = "bold"),
        strip.text.x = element_blank(),
        panel.background = element_rect(colour = "grey"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
    )
# guides(fill = guide_legend(override.aes = aes(label = ""))) +
p_volcano + ggsave(
    file = "./data/99_results/volcano_MR.png",
    units = "cm",
    width = 27,
    height = 15
)


## Forest plot
## The 95% CI was based on the delta method standard error
res_single$myFacet <- res_single$b < -0.678392

p_forest <- res_single %>%
    filter(association != "NEUTRAL") %>%
    arrange(desc(b)) %>%
    forestplot(
        df = .,
        estimate = b,
        logodds = TRUE,
        se = se.delta,
        name = gene.exposure,
        colour = association
    ) +
    scale_color_manual(
        values = c(
            "DOWN" = "#32A29B",
            "UP" = "#EC7200"
            # "EXCLUDE" = "#474554"
        )
    ) +
    labs(
        x = "OR for MR effect size for Proteomics || Alzheimer's disease"
    ) +
    theme_minimal() +
    theme(
        legend.position = "none",
        axis.title.y = element_blank()
    )

p_forest + ggsave(
    file = "./data/99_results/forrest_MR.png"
    # width = 20,
    # height = 16,
    # units = "cm"
)

# TO DO create patchwork of both of them.
p_volcano + p_forest +
    plot_layout(
        widths = c(2, 1),
        # heights = unit(c(5, 1),
        # c('cm', 'null'))
    ) +
    ggsave(
        file = "./data/99_results/MR_combined.png",
        units = "cm",
        width = 40,
        height = 15
    )


# Write data ---------------------------------------------------------------
