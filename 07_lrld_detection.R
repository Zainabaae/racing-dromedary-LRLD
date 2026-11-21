# 07_lrld_detection.R
# Detection of long-range LD (LRLD) events.
# An LRLD event is a pair of distinct haplotype blocks (each with > 4 SNPs) that are
#   (i)  separated by >= 1 Mb (start of block B minus end of block A),
#   (ii) linked by >= 2 SNP pairs with r2 >= 0.6, and
#   (iii) involve >= 2 distinct SNPs in each block.
# Input : LD/drom_race_omn_kw_LD_chr_N.ld, LD_blocks/drom_race_omn_kw_chr_N_haplotypes.blocks.det
# Output: per chromosome
#   detail_LRLD_4_SNPs_chr_N.txt                          SNP pairs linking block pairs (criteria i-ii)
#   drom_race_omn_kw_unique_detailed_LRLD_4_SNPs_chr_N.txt one row per LRLD event (criteria i-iii)
#   drom_race_omn_kw_LRLD_4_SNPs_chr_N_blks.bed           blocks involved in LRLD events (for gene annotation)

library(GenomicRanges)

project_dir <- "/path/to/project"
ld_dir      <- file.path(project_dir, "LD")
blk_dir     <- file.path(project_dir, "LD_blocks")

## Step 1: keep SNP pairs with r2 >= 0.6
for (a in 1:36) {
  r    <- read.table(file.path(ld_dir, paste0("drom_race_omn_kw_LD_chr_", a, ".ld")),
                     header = TRUE, stringsAsFactors = FALSE)
  high <- r[r$R2 >= 0.6, ]
  write.table(high, file.path(ld_dir, paste0("high_r_drom_race_omn_kw_chr_", a, ".txt")),
              row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
}

## Step 2: assign SNPs to haplotype blocks; keep block pairs >= 1 Mb apart with >= 2 SNP pairs
for (a in 1:36) {
  ld <- read.table(file.path(ld_dir, paste0("high_r_drom_race_omn_kw_chr_", a, ".txt")),
                   header = TRUE, stringsAsFactors = FALSE)
  blocks <- read.table(file.path(blk_dir, paste0("drom_race_omn_kw_chr_", a, "_haplotypes.blocks.det")),
                       header = TRUE, stringsAsFactors = FALSE)
  blocks <- blocks[blocks$NSNPS > 4, ]
  blocks <- GRanges(blocks$CHR, IRanges(blocks$BP1, blocks$BP2))

  ld_tmp <- lapply(1:2, function(i) {
    if (i == 1) g <- unique(GRanges(ld$CHR_A, IRanges(ld$BP_A, ld$BP_A), BLOCK_A = NA, BLOCK_A_BP1 = NA, BLOCK_A_BP2 = NA))
    if (i == 2) g <- unique(GRanges(ld$CHR_B, IRanges(ld$BP_B, ld$BP_B), BLOCK_A = NA, BLOCK_A_BP1 = NA, BLOCK_A_BP2 = NA))
    ov <- findOverlaps(g, blocks)
    g[queryHits(ov)]$BLOCK_A     <- subjectHits(ov)
    g[queryHits(ov)]$BLOCK_A_BP1 <- start(blocks[subjectHits(ov)])
    g[queryHits(ov)]$BLOCK_A_BP2 <- end(blocks[subjectHits(ov)])
    if (i == 2) names(mcols(g)) <- c("BLOCK_B", "BLOCK_B_BP1", "BLOCK_B_BP2")
    g <- data.frame(g)
    g[, c(1:2, 6:8)]
  })

  ld <- merge(ld, ld_tmp[[1]], by.x = c("CHR_A", "BP_A"), by.y = c("seqnames", "start"))
  ld <- merge(ld, ld_tmp[[2]], by.x = c("CHR_B", "BP_B"), by.y = c("seqnames", "start"))
  ld <- ld[order(ld$BP_A), ]

  f.ld <- ld[which(ld$BLOCK_A != ld$BLOCK_B), ]               # SNPs in different blocks
  f.ld$block_AB <- paste0("Block_", f.ld$BLOCK_A, "_", f.ld$BLOCK_B)
  f.ld$distance <- f.ld$BLOCK_B_BP1 - f.ld$BLOCK_A_BP2        # inter-block distance

  LR    <- f.ld[f.ld$distance >= 1e6, ]
  freq  <- table(LR$block_AB)
  keep  <- names(freq)[freq > 1]                              # >= 2 SNP pairs
  f.LRLD <- f.ld[f.ld$block_AB %in% keep, ]

  write.table(f.LRLD, file.path(blk_dir, paste0("detail_LRLD_4_SNPs_chr_", a, ".txt")),
              row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
}

## Step 3: keep block pairs with >= 2 distinct SNPs in each block; write events and blocks
for (a in 1:36) {
  tmp <- read.delim(file.path(blk_dir, paste0("detail_LRLD_4_SNPs_chr_", a, ".txt")),
                    header = TRUE, stringsAsFactors = FALSE)
  nA <- tapply(tmp$SNP_A, tmp$block_AB, function(x) length(unique(x)))
  nB <- tapply(tmp$SNP_B, tmp$block_AB, function(x) length(unique(x)))
  keep <- names(nA)[nA >= 2 & nB[names(nA)] >= 2]

  d.LRLD <- tmp[tmp$block_AB %in% keep, ]
  events <- unique(d.LRLD[, c("BLOCK_A", "BLOCK_A_BP1", "BLOCK_A_BP2",
                              "BLOCK_B", "BLOCK_B_BP1", "BLOCK_B_BP2", "block_AB")])
  write.table(events, file.path(blk_dir, paste0("drom_race_omn_kw_unique_detailed_LRLD_4_SNPs_chr_", a, ".txt")),
              row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)

  if (nrow(events) > 0) {
    blks <- unique(rbind(
      data.frame(chr = a, start = events$BLOCK_A_BP1, stop = events$BLOCK_A_BP2, blk = paste0("Block_", events$BLOCK_A)),
      data.frame(chr = a, start = events$BLOCK_B_BP1, stop = events$BLOCK_B_BP2, blk = paste0("Block_", events$BLOCK_B))))
    write.table(blks, file.path(blk_dir, paste0("drom_race_omn_kw_LRLD_4_SNPs_chr_", a, "_blks.bed")),
                row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
  }
}
