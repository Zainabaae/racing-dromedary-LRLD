#!/bin/bash
# 00_config.sh
# Paths and software used by all shell scripts. Edit these lines only.
# Usage: the other shell scripts run `source scripts/00_config.sh` from the repository root.

PROJECT_DIR="/path/to/project"                 # top-level working directory
PLINK="plink"                                  # PLINK v1.9 executable
GENOTYPES="${PROJECT_DIR}/genotypes"           # PLINK .ped/.map files
LD_DIR="${PROJECT_DIR}/LD"                     # pairwise LD output (.ld)
BLOCK_DIR="${PROJECT_DIR}/LD_blocks"           # haplotype blocks and LRLD output
REF_FAI="${PROJECT_DIR}/reference/CamDro3.fa.fai"     # reference genome index (chromosome lengths)
GENE_BED="${PROJECT_DIR}/reference/CamDro3_genes.bed" # gene annotation (chr, start, end, gene)

# Input genotype prefix (all samples, autosomal SNPs) and list of racing camels to keep
RAW_PREFIX="${GENOTYPES}/drom_all_samples_autosomal"
KEEP_LIST="${GENOTYPES}/keep_racing_omn_kw.txt"

# Prefix used for the analysed racing dataset (Kuwait + Oman)
DATASET="drom_race_omn_kw"
OPTS="--allow-no-sex --autosome-num 36"        # dromedary: 36 autosomes
