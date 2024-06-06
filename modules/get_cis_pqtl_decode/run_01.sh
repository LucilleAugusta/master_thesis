#
# ./....
#
# 0. get gen bed file
# loop step 1-3
# 1. download summary stat
# 2. extract generegion plus/minus 10 MB + variants with p < 5e-6
# 3. clean intermediate
#
#
# ------------------------------------
#
# 0. extract generegion plus/minus 10 MB
#

PATH=/home/gua/bin/bin:${PATH} 

gunzip -c gencode.v44.annotation.gtf.gz \
    | gtf2bed - \
    | grep -w gene \
    > ../../data/resources/gencodeGenes.bed

grep "protein_coding" ../../data/resources/gencodeGenes.bedgencodeGenes.bed \
     ../../data/resources/gencodeGenes.bed > pc_gencodeGenes.bed



# ------------------------------------
#
# 1. extract generegion plus/minus 10 MB
#


gene=name

wget 'pcsk9' -O $outfolder




# ------------------------------------
#
# 2. extract generegion plus/minus 10 MB
#

# sumstat=...summary_stat.gz

# 10 MB
locussize=10000000

