#!/bin/bash
# 09_gene_annotation.sh
# Genes overlapping the haplotype blocks involved in LRLD events (genes located in the
# interval between the two blocks of an event are NOT included).
# Software: BEDTools (intersectBed)

set -euo pipefail
source "$(dirname "$0")/00_config.sh"
cd "${BLOCK_DIR}"

for i in {1..36}; do
  f=${DATASET}_LRLD_4_SNPs_chr_${i}_blks.bed
  [ -s "$f" ] || continue          # skip chromosomes without LRLD events
  intersectBed -wo -a ${GENE_BED} -b $f > ${DATASET}_LRLD_4_SNPs_chr_${i}_genes.txt
done

cat ${DATASET}_LRLD_4_SNPs_chr_*_genes.txt > ${DATASET}_LRLD_genes_autosomes.txt
echo "Unique genes overlapping LRLD blocks:"
cut -f4 ${DATASET}_LRLD_genes_autosomes.txt | sort -u | wc -l
