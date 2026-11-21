# 04a_sample_snps_for_LD.R
# Randomly samples 30% of the QC'd SNPs on each autosome; pairwise LD (04b) is computed on
# these SNPs only, to keep the number of SNP pairs computationally manageable.
# NOTE: the original analysis was run without a fixed random seed. The SNP lists actually
# used are provided in data/snp_lists/ so that the analysis can be reproduced exactly;
# set.seed() is included here for future re-runs.

library(GenABEL)

project_dir <- "/path/to/project"
setwd(file.path(project_dir, "genotypes"))

convert.snp.ped("drom_race_omn_kw_QC_IBS.ped", "drom_race_omn_kw_QC_IBS.map",
                out = "HDdata_full.raw", mapHasHeaderLine = FALSE)
data <- load.gwaa.data(phe = "racing_project_pheno.txt", gen = "HDdata_full.raw", force = TRUE, sort = TRUE)

out_dir <- file.path(project_dir, "LD")
dir.create(out_dir, showWarnings = FALSE)
set.seed(2026)

for (i in 1:36) {
  snps <- data@gtdata@snpnames[data@gtdata@chromosome == i]
  pick <- sample(snps, round(0.30 * length(snps)), replace = FALSE)
  write.table(pick, file.path(out_dir, paste0("snp_for_LD_drom_race_omn_kw_chr_", i, ".txt")),
              row.names = FALSE, col.names = FALSE, quote = FALSE)
}
