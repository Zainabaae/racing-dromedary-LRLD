#!/bin/bash
# 01_qc_plink.sh
# Quality control of autosomal SNP genotypes of Kuwaiti and Omani racing dromedaries.
# Software: PLINK v1.9, R
#
# Output:
#   ${DATASET}_QC_IBS      -> QC'd dataset (38 samples, 3,742,192 SNPs) used for LD, blocks and relatedness
#   ${DATASET}_QC_IBS_LD   -> LD-pruned subset (493,441 SNPs) used for PCA and homozygosity only

set -euo pipefail
source "$(dirname "$0")/00_config.sh"
cd "${GENOTYPES}"

# 1. Keep Kuwaiti and Omani racing camels
${PLINK} ${OPTS} --file ${RAW_PREFIX} --keep ${KEEP_LIST} --recode --out ${DATASET}_autosomal

# 2. Sample and SNP filters (PLINK applies --mind before the variant filters)
#    --mind 0.05 : remove samples with call rate < 95%
#    --geno 0.1  : remove SNPs with > 10% missing genotypes
#    --maf 0.05  : remove SNPs with minor allele frequency < 0.05
#    --hwe 1e-6  : remove SNPs deviating from Hardy-Weinberg equilibrium (P < 1e-6)
${PLINK} ${OPTS} --file ${DATASET}_autosomal \
  --geno 0.1 --maf 0.05 --hwe 0.000001 --mind 0.05 \
  --recode --out ${DATASET}_QC

# 3. Pairwise identity-by-state (IBS); remove one sample of each pair with IBS (DST) > 0.95
${PLINK} ${OPTS} --file ${DATASET}_QC --genome --out ${DATASET}_QC_genome
awk '{print $1,$3,$12}' ${DATASET}_QC_genome.genome > ${DATASET}_QC_genome_DST.txt   # FID1 FID2 DST

Rscript - "${DATASET}" <<'EOF'
ds  <- commandArgs(trailingOnly = TRUE)[1]
tmp <- read.table(paste0(ds, "_QC_genome_DST.txt"), header = TRUE, stringsAsFactors = FALSE)
cat("Mean pairwise IBS:", mean(tmp$DST), "\n")
exc <- tmp[which(tmp$DST > 0.95), ]
# FID equals IID in these files
out <- unique(data.frame(exc$FID1, exc$FID1))
write.table(out, "ibs_exclude_sample.txt", row.names = FALSE, col.names = FALSE,
            sep = "\t", quote = FALSE)
EOF

${PLINK} ${OPTS} --file ${DATASET}_QC --remove ibs_exclude_sample.txt \
  --recode --out ${DATASET}_QC_IBS

# 4. LD pruning (for PCA and homozygosity only): window 50 SNPs, step 5, r2 0.5
${PLINK} ${OPTS} --file ${DATASET}_QC_IBS --indep-pairwise 50 5 0.5 --out ${DATASET}_QC_IBS_LD
${PLINK} ${OPTS} --file ${DATASET}_QC_IBS --extract ${DATASET}_QC_IBS_LD.prune.in \
  --recode --out ${DATASET}_QC_IBS_LD
