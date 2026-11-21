# 02_relatedness_king.R
# Classifies pairs of racing dromedaries by KING kinship coefficient (phi).
# Input : <dataset>.relatedness2 (vcftools --relatedness2 output; KING-robust estimator)
#         [Add here the exact command used to produce this file, e.g.
#          vcftools --vcf <QC'd, LD-pruned VCF> --relatedness2 --out <dataset>]
# Output: lists of first-, second-, third-degree and unrelated pairs
# Thresholds follow Manichaikul et al. (2010).

project_dir <- "/path/to/project"
setwd(file.path(project_dir, "relatedness"))

rel   <- read.delim("drom_race_omn_kw_QC_IBS_LD.relatedness2", header = TRUE, stringsAsFactors = FALSE)
rel.f <- rel[, c(1, 2, 7)]                    # INDV1, INDV2, RELATEDNESS_PHI
rel.f <- rel.f[rel.f[, 1] != rel.f[, 2], ]    # drop self-comparisons
phi   <- rel.f[, 3]

cat("Mean phi:", mean(phi), " SD:", sd(phi), "\n")

classes <- list(
  duplicate_or_MZ = rel.f[phi > 0.354, ],
  first_degree    = rel.f[phi >= 0.177  & phi <= 0.354, ],
  second_degree   = rel.f[phi >= 0.0884 & phi <  0.177, ],
  third_degree    = rel.f[phi >= 0.0442 & phi <  0.0884, ],
  unrelated       = rel.f[phi <  0.0442, ])

for (k in names(classes)) {
  cat(k, ":", nrow(classes[[k]]), "pairs\n")
  write.table(classes[[k]], file = paste0(k, "_pairs_race_omn_kw.txt"),
              row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
}
