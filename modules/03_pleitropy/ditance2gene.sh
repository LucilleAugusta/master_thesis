#!/bin/bash
# Pre-req: clean files from "./modules/03_pleitropy/run.sh" if files have data
#           touch count_cis_gene.txt
# run: ./modules/03_pleitropy/ditance2gene.sh 

# Set working directory --------------------------------------
cd /opt/projects/0002_ZS/00021_MastersThesis_LAG

# Path to the file with rs id, gene, and distance information
#rs_file="data/01_decode_data/associated_SNPs/pleiotropy/rs5841.txt"
input_dir="data/01_decode_data/associated_SNPs/pleiotropy_PAV"
count_file="data/01_decode_data/associated_SNPs/pleiotropy_PAV/count_cis_gene.txt"


for FILE in "$input_dir"/rs*
  do 
    
    Rscript ./modules/03_pleitropy/01_distance2gene_PAV.R "$FILE"
    
    # Do some kind of count for each file. 
    # try to add "my hit" from data/02_analysis/01_mr/mr_significant_meta_WIDE.csv

    #Count cis_gene_start
    snps_found=$(cat "$FILE" | wc -l)
    mr_hit=$(cat "$FILE" | awk '{print $10}' | tail -1 ) #only printing one time 
    start_count=$(cat "$FILE" | awk '{print $13}' | grep "TRUE" | wc -l)
    end_count=$(cat "$FILE" | awk '{print $14}' | grep "TRUE" | wc -l)
    
    echo "$FILE $mr_hit $snps_found $start_count $end_count" >> "data/01_decode_data/associated_SNPs/pleiotropy_PAV/count_cis_gene.txt"
  done


# Loop over the list of genes
# for gene in "${genes[@]}"
# do
#   echo "Processing ${gene}"
#   cat "$file" | grep -w "$gene" | grep protein_coding | \
#   awk -v ch="$chromosome" -v pos="$position" \
#   '$1 == ch &&  && $3 <= pos+1000000 && $3 >= pos-1000000 {print $2}'
# done


# start pos: cat data/resources/pc_gencodeGenes_filtered.csv | grep -w "NME4"| awk '{print $5}'
# end pos: cat data/resources/pc_gencodeGenes_filtered.csv | grep -w "NME4"| awk '{print $6}'
