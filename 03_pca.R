# 03_pca.R
# Principal component analysis of the 38 racing dromedaries using LD-pruned SNPs.
# Packages: GenABEL (data import), base R stats::prcomp (PCA)
# Input : <dataset>_QC_IBS_LD.ped/.map, phenotype file with columns id and pop
#         (pop codes: 15 = racing Kuwait, 20 = racing Oman (Sharqiya), 23 = racing Oman (Batinah))
# Output: PCA plots (PC1 vs PC2, PC1 vs PC3) and PC coordinates

library(GenABEL)

project_dir <- "/path/to/project"
setwd(file.path(project_dir, "genotypes"))

convert.snp.ped("drom_race_omn_kw_QC_IBS_LD.ped", "drom_race_omn_kw_QC_IBS_LD.map",
                out = "HDdata.raw", mapHasHeaderLine = FALSE)
data <- load.gwaa.data(phe = "racing_project_pheno.txt", gen = "HDdata.raw", force = TRUE, sort = TRUE)

# Genotypes coded 0/1/2; missing genotypes replaced by the SNP mean; SNPs centred and scaled
geno <- as.numeric(data@gtdata)
col.avg <- colMeans(geno, na.rm = TRUE)
for (i in seq_along(col.avg)) geno[is.na(geno[, i]), i] <- col.avg[i]
geno <- geno[, apply(geno, 2, var) > 0]          # drop monomorphic SNPs (required for scaling)
pca  <- prcomp(geno, center = TRUE, scale. = TRUE)
pve  <- round(100 * summary(pca)$importance[2, 1:3], 1)
print(pve)

pc    <- data.frame(id = rownames(pca$x), pca$x[, 1:3], stringsAsFactors = FALSE)
pheno <- read.table("racing_project_pheno.txt", header = TRUE, stringsAsFactors = FALSE)
pc    <- merge(pc, pheno[, c("id", "pop")], by = "id")

groups <- list("Racing Kuwait"            = list(pop = 15, col = "blue"),
               "Racing Oman (Batinah)"    = list(pop = 23, col = "red"),
               "Racing Oman (Sharqiya)"   = list(pop = 20, col = "darkgreen"))

dir.create(file.path(project_dir, "PCA"), showWarnings = FALSE)
setwd(file.path(project_dir, "PCA"))

plot_pcs <- function(a, b, fname) {
  tiff(fname, width = 150, height = 150, res = 300, units = "mm", compression = "lzw")
  par(cex = 0.7)
  plot(pc[[paste0("PC", a)]], pc[[paste0("PC", b)]], type = "n",
       xlab = paste0("PC", a, " (", pve[a], "%)"), ylab = paste0("PC", b, " (", pve[b], "%)"))
  for (g in names(groups)) {
    s <- pc[pc$pop == groups[[g]]$pop, ]
    points(s[[paste0("PC", a)]], s[[paste0("PC", b)]], pch = 19, col = groups[[g]]$col)
  }
  abline(h = 0, v = 0, lty = 3)
  legend("topleft", pch = 19, col = sapply(groups, `[[`, "col"), legend = names(groups))
  dev.off()
}
plot_pcs(1, 2, "PCA_race_omn_kw_PC1_PC2.tiff")
plot_pcs(1, 3, "PCA_race_omn_kw_PC1_PC3.tiff")

write.table(pc, "PCA_race_omn_kw_coordinates.txt", row.names = FALSE, sep = "\t", quote = FALSE)
