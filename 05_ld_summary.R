# 05_ld_summary.R
# Mean and SD of r2 per autosome. Note: only SNP pairs reported by PLINK (r2 >= 0.099,
# within 3 Mb) contribute, so values describe the retained pairs, not all possible pairs.

project_dir <- "/path/to/project"
setwd(file.path(project_dir, "LD"))

res <- data.frame()
for (i in 1:36) {
  ld  <- read.table(paste0("drom_race_omn_kw_LD_chr_", i, ".ld"), header = TRUE, stringsAsFactors = FALSE)
  res <- rbind(res, data.frame(chr = i, n_pairs = nrow(ld), mean_r2 = mean(ld$R2), sd_r2 = sd(ld$R2)))
}
write.table(res, "LD_summary_race_omn_kw.txt", row.names = FALSE, sep = "\t", quote = FALSE)
cat("Mean of per-chromosome mean r2:", mean(res$mean_r2), "\n")
cat("Mean of per-chromosome SD r2  :", mean(res$sd_r2), "\n")
