#
# Metal GWAS meta-analysis
#
# --------------------------------------------------
#
#  env
#


PATH=${PATH}:/home/lag/master_thesis/modules/00_meta_analyses/generic-metal


# --------------------------------------------------
#
#  lift over 37 -> 38 Kunkle
#

#...


# --------------------------------------------------
#
#  Kunkle + Finngen hg38
#
# config AD_kunkle_finngen.txt
# run w/ SAMPLESIZE anf STDERR


metal AD_kunkle_finngen.txt 

metal AD_kunkle_finngen_wide.txt 

# --------------------------------------------------
#
#  Post analyses filter 
#


# input=../../data/metal_ms_science_finngen_hg37_1.tsv.gz
# Rscript postAnalysisFiltering.R $input


# input=../../data/metal_ms_science_finngen_hg37_noMHC_1.tsv.gz
# Rscript postAnalysisFiltering.R $input


# input=../../data/metal_ms_science_finngen_hg37_samplesize1.tsv.gz
# Rscript postAnalysisFiltering.R $input


# input=../../data/metal_ms_science_finngen_hg37_noMHC_samplesize1.tsv.gz
# Rscript postAnalysisFiltering.R $input



