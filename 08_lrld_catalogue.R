# 08_lrld_catalogue.R
# Builds the LRLD event catalogue (Supplementary Table S5) from the output of 07_lrld_detection.R
# (detail_LRLD_4_SNPs_chr_*.txt), applying the same final filter (>= 2 distinct SNPs per block),
# and calculates genomic coverage of LRLD blocks (each block counted once).
#
# Output (in the same folder):
#   lrld_snp_pairs.tsv        one row per supporting SNP pair (r2 >= 0.6)
#   lrld_event_catalogue.tsv  one row per LRLD event (block pair)

project_dir <- "/path/to/project"
setwd(file.path(project_dir, "LD_blocks"))

pairs_all  <- list()
events_all <- list()

for (a in 1:36) {
  f <- paste("detail_LRLD_4_SNPs_chr_", a, ".txt", sep = "")
  if (!file.exists(f)) { cat("missing:", f, "\n"); next }
  tmp <- read.delim(f, header = TRUE, stringsAsFactors = FALSE)
  if (nrow(tmp) == 0) next

  # same final filter as the detection script: >= 2 distinct SNPs in each block
  nA <- tapply(tmp$SNP_A, tmp$block_AB, function(x) length(unique(x)))
  nB <- tapply(tmp$SNP_B, tmp$block_AB, function(x) length(unique(x)))
  keep <- names(nA)[nA >= 2 & nB[names(nA)] >= 2]
  d <- tmp[tmp$block_AB %in% keep, ]
  if (nrow(d) == 0) next

  d$event_id <- paste("chr", a, "_", d$block_AB, sep = "")
  d$snp_distance_bp <- abs(d$BP_B - d$BP_A)

  pairs_all[[a]] <- data.frame(
    event_id = d$event_id, chr = a,
    SNP_A = d$SNP_A, BP_A = d$BP_A, SNP_B = d$SNP_B, BP_B = d$BP_B,
    snp_distance_bp = d$snp_distance_bp, r2 = d$R2, stringsAsFactors = FALSE)

  ev <- unique(d[, c("event_id", "BLOCK_A_BP1", "BLOCK_A_BP2", "BLOCK_B_BP1", "BLOCK_B_BP2")])
  ev$chr <- a
  ev$blockA_length_bp         <- ev$BLOCK_A_BP2 - ev$BLOCK_A_BP1 + 1
  ev$blockB_length_bp         <- ev$BLOCK_B_BP2 - ev$BLOCK_B_BP1 + 1
  ev$inter_block_distance_bp  <- ev$BLOCK_B_BP1 - ev$BLOCK_A_BP2
  ev$outer_interval_length_bp <- ev$BLOCK_B_BP2 - ev$BLOCK_A_BP1 + 1
  ev$n_supporting_pairs <- as.vector(table(d$event_id)[ev$event_id])
  ev$n_snps_blockA <- as.vector(nA[sub(paste("chr", a, "_", sep = ""), "", ev$event_id)])
  ev$n_snps_blockB <- as.vector(nB[sub(paste("chr", a, "_", sep = ""), "", ev$event_id)])
  ev$snps_blockA <- as.vector(tapply(d$SNP_A, d$event_id,
                       function(x) paste(sort(unique(x)), collapse = ","))[ev$event_id])
  ev$snps_blockB <- as.vector(tapply(d$SNP_B, d$event_id,
                       function(x) paste(sort(unique(x)), collapse = ","))[ev$event_id])
  ev$max_r2  <- as.vector(tapply(d$R2, d$event_id, max)[ev$event_id])
  ev$mean_r2 <- as.vector(tapply(d$R2, d$event_id, mean)[ev$event_id])
  events_all[[a]] <- ev
}

pairs  <- do.call(rbind, pairs_all)
events <- do.call(rbind, events_all)
events <- events[, c("event_id", "chr",
                     "BLOCK_A_BP1", "BLOCK_A_BP2", "BLOCK_B_BP1", "BLOCK_B_BP2",
                     "blockA_length_bp", "blockB_length_bp", "inter_block_distance_bp",
                     "outer_interval_length_bp", "n_snps_blockA", "n_snps_blockB",
                     "n_supporting_pairs", "max_r2", "mean_r2",
                     "snps_blockA", "snps_blockB")]
names(events)[3:6] <- c("blockA_start", "blockA_end", "blockB_start", "blockB_end")
events <- events[order(events$chr, events$blockA_start), ]
pairs  <- pairs[order(pairs$chr, pairs$BP_A, pairs$BP_B), ]

write.table(pairs,  "lrld_snp_pairs.tsv",       sep = "\t", quote = FALSE, row.names = FALSE)
write.table(events, "lrld_event_catalogue.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

cat("LRLD events:", nrow(events), "\n")
cat("Supporting SNP pairs:", nrow(pairs), "\n")
cat("Distinct blocks involved:",
    length(unique(c(paste(events$chr, events$blockA_start), paste(events$chr, events$blockB_start)))), "\n")
cat("Largest outer interval (bp):", max(events$outer_interval_length_bp), "\n")
cat("Smallest outer interval (bp):", min(events$outer_interval_length_bp), "\n")

# Genomic coverage of LRLD blocks (union of block intervals; each block counted once)
library(GenomicRanges)
gr  <- reduce(GRanges(c(events$chr, events$chr),
                      IRanges(c(events$blockA_start, events$blockB_start),
                              c(events$blockA_end,   events$blockB_end))))
fai <- read.table(file.path(project_dir, "reference", "CamDro3.fa.fai"), stringsAsFactors = FALSE)
aut <- sum(as.numeric(fai$V2[fai$V1 %in% as.character(1:36)]))
cat("LRLD block coverage (bp):", sum(width(gr)), " = ", round(100 * sum(width(gr)) / aut, 2), "% of autosomes\n")
