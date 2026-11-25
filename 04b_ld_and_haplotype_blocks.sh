#!/bin/bash
# 04b_ld_and_haplotype_blocks.sh
# (1) Pairwise LD (r2) per autosome, using the SNP subset from 04a.
#     All pairs with r2 >= 0.099 within 3 Mb are reported; the maximum distance searched is 3 Mb.
# (2) Haplotype blocks per autosome (Gabriel et al. method in PLINK), using all QC'd SNPs.
# Software: PLINK v1.9

set -euo pipefail
source "$(dirname "$0")/00_config.sh"
mkdir -p "${LD_DIR}" "${BLOCK_DIR}"

for i in {1..36}; do
  # Pairwise LD
  ${PLINK} ${OPTS} --file ${GENOTYPES}/${DATASET}_QC_IBS \
    --extract ${LD_DIR}/snp_for_LD_${DATASET}_chr_${i}.txt \
    --r2 --ld-window 1000000 --ld-window-kb 3000 --ld-window-r2 0.099 \
    --out ${LD_DIR}/${DATASET}_LD_chr_${i}

  # Haplotype blocks
  ${PLINK} ${OPTS} --file ${GENOTYPES}/${DATASET}_QC_IBS --chr ${i} \
    --blocks no-pheno-req no-small-max-span --blocks-max-kb 3000 \
    --out ${BLOCK_DIR}/${DATASET}_chr_${i}_haplotypes
done
