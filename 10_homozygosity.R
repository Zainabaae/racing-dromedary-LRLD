# 10_homozygosity.R
# Individual homozygosity and inbreeding coefficient (F) per racing group, using LD-pruned SNPs.
# Packages: GenABEL
# Pop codes: 15 = racing Kuwait, 20 = racing Oman (Sharqiya), 23 = racing Oman (Batinah)

library(GenABEL)

project_dir <- "/path/to/project"
setwd(file.path(project_dir, "genotypes"))

convert.snp.ped("drom_race_omn_kw_QC_IBS_LD.ped", "drom_race_omn_kw_QC_IBS_LD.map",
                out = "HDdata.raw", mapHasHeaderLine = FALSE)
data <- load.gwaa.data(phe = "racing_project_pheno.txt", gen = "HDdata.raw", force = TRUE, sort = TRUE)

H     <- hom(data)
hz    <- data.frame(id = rownames(H), H, stringsAsFactors = FALSE)
pheno <- read.table("racing_project_pheno.txt", header = TRUE, stringsAsFactors = FALSE)
hz    <- merge(hz, pheno[, c("id", "pop")], by = "id")
labels <- c("15" = "race.kw", "20" = "race.omn.sharqiya", "23" = "race.omn.batinah")
hz$group <- labels[as.character(hz$pop)]

dir.create(file.path(project_dir, "homozygosity"), showWarnings = FALSE)
setwd(file.path(project_dir, "homozygosity"))
write.table(hz, "homozygosity_race_omn_kw.txt", row.names = FALSE, sep = "\t", quote = FALSE)

# Summary per group (Wilcoxon signed-rank test of F against 0)
summ <- do.call(rbind, lapply(split(hz, hz$group), function(t)
  data.frame(group = t$group[1], mean_Hom = mean(t$Hom), sd_Hom = sd(t$Hom),
             mean_F = mean(t$F), F_p_value = wilcox.test(t$F)$p.value)))
write.table(summ, "homozygosity_summary_race_omn_kw.txt", row.names = FALSE, sep = "\t", quote = FALSE)

# Pairwise Mann-Whitney tests between groups
g <- unique(hz$group); pw <- data.frame()
for (i in seq_along(g)) for (j in seq_along(g)) if (i < j) {
  a <- hz[hz$group == g[i], ]; b <- hz[hz$group == g[j], ]
  pw <- rbind(pw, data.frame(group1 = g[i], group2 = g[j],
                             Hom_p_value = wilcox.test(a$Hom, b$Hom)$p.value,
                             F_p_value   = wilcox.test(a$F,   b$F)$p.value))
}
write.table(pw, "MW_test_hom_F_race_omn_kw.txt", row.names = FALSE, sep = "\t", quote = FALSE)
