#!/bin/bash
## run by: ./modules/03_pleitropy/run_clump.sh
# rm data/01_decode_data/associated_SNPs/pleiotropy/rs*

# Define your list of rsids
rsids=( "rs3730025"
"rs1859788"
"rs144092339"
"rs3745997"
"rs55667375"
"rs9299404"
"rs4984667"
"rs2075803"
"rs11078596"
"rs11248061"
"rs5848"
"rs61822669"
"rs2523708"
"rs651279"
"rs258341"
"rs429358"
"rs150037800"
"rs35103166" )

directory="data/01_decode_data/associated_SNPs"
cd ${directory}
# Loop over the list
for rs in "${rsids[@]}"
do
  echo "Processing ${rs}"
  touch ./pleiotropy_PAV/"$rs".txt
  grep -w -H ${rs} *.txt | awk '$8 < 5e-08' >> ./pleiotropy_PAV/"$rs".txt
  #grep -w -H ${rs} * | awk '$8 < 5e-08' | wc >> ./pleiotropy/"$rs".txt
done
