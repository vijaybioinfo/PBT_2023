#!/usr/bin/R

#########################################
# Shannon and Simpson Diversity indexes #
#########################################

# ---
# Author: Kevin Meza Landeros
# Date: 2022-01-15
# ---


VDJTOOLS="java -jar /home/kmlanderos/bin/vdjtools-1.2.1/vdjtools-1.2.1.jar"

OUT_DIR="/home/kmlanderos/tmp_large/pbtumor-all-Batch2/results/figures/tcr/vdjtools/outs/noJgene_removeMultipleBchain"
mkdir -p $OUT_DIR
# After liminating clonotypes with multiple beta chains
$VDJTOOLS CalcDiversityStats -m /home/kmlanderos/tmp_large/pbtumor-all-Batch2/results/figures/tcr/vdjtools/PBT_metadata_cd8_V2.tsv $OUT_DIR/PBT_cd8
$VDJTOOLS CalcDiversityStats -m /home/kmlanderos/tmp_large/pbtumor-all-Batch2/results/figures/tcr/vdjtools/PBT_metadata_cd4_V2.tsv $OUT_DIR/PBT_cd4
$VDJTOOLS CalcDiversityStats -m /home/kmlanderos/tmp_large/pbtumor-all-Batch2/results/figures/tcr/vdjtools/LUNG_metadata_cd8_V2.tsv $OUT_DIR/LUNG_cd8
$VDJTOOLS CalcDiversityStats -m /home/kmlanderos/tmp_large/pbtumor-all-Batch2/results/figures/tcr/vdjtools/LUNG_metadata_cd4_V2.tsv $OUT_DIR/LUNG_cd4
