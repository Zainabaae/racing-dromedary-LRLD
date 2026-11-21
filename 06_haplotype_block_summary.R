# 06_haplotype_block_summary.R
# Number, SNP content and total size of haplotype blocks (> 4 SNPs) per autosome,
# and their coverage of the autosomal genome.
# Block sizes are computed in bp from BP1/BP2 (the PLINK "KB" column is in kilobases).

project_dir <- "/path/to/project"
setwd(file.path(project_dir, "LD_blocks"))

fai <- read.table(file.path(project_dir, "reference", "CamDro3.fa.fai"), stringsAsFactors = FALSE)
chr_len <- setNames(fai$V2, fai$V1)

res <- data.frame()
for (i in 1:36) {
  blk <- read.table(paste0("drom_race_omn_kw_chr_", i, "_haplotypes.blocks.det"),
                    header = TRUE, stringsAsFactors = FALSE)
  blk <- blk[blk$NSNPS > 4, ]
  size <- sum(blk$BP2 - blk$BP1 + 1)
  res <- rbind(res, data.frame(chr = i, n_blocks = nrow(blk), n_snps = sum(blk$NSNPS),
                               total_size_bp = size, chr_length_bp = chr_len[as.character(i)],
                               percent = 100 * size / chr_len[as.character(i)]))
}
write.table(res, "haplotype_blocks_summary_race_omn_kw.txt", row.names = FALSE, sep = "\t", quote = FALSE)

aut <- sum(chr_len[as.character(1:36)])
cat("Blocks:", sum(res$n_blocks), " SNPs:", sum(res$n_snps), "\n")
cat("Total block size (bp):", sum(res$total_size_bp),
    " = ", round(100 * sum(res$total_size_bp) / aut, 2), "% of", aut, "bp autosomes\n")
