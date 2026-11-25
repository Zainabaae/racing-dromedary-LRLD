# racing-dromedary-LRLD

Analysis scripts for the genome-wide characterization of long-range linkage disequilibrium (LRLD) in racing dromedary camels (*Camelus dromedarius*) from Kuwait and Oman.

> Accompanies: Alfoudari A, Mohammad Z, Al Abri M, Bahbahani H. *Genome-wide characterization of long-range linkage disequilibrium in racing dromedary camels from Kuwait and Oman.* Journal of Animal Science (under review).
>
> Contact: Hussain Bahbahani (hussain.bahbahani@ku.edu.kw), Department of Biological Sciences, Kuwait University.

## Data

- **Sequence reads:** European Nucleotide Archive (ENA). The newly sequenced Omani racing camels are under study **PRJEB123824**; run accessions of all analysed animals are listed in Supplementary Table S1 of the article.
- **Reference genome:** CamDro3 (36 autosomes).
- **LRLD event catalogue:** the 4,444 LRLD events and the SNPs supporting them are provided in Supplementary Table S5 of the article.

Read processing (fastp), alignment (BWA-MEM), duplicate removal (Picard) and variant calling and filtering (GATK, bcftools) are described in the article; this repository starts from the filtered autosomal SNP genotypes in PLINK format.

## Software

| Software | Version | Used for |
|---|---|---|
| PLINK | 1.9 | QC, IBS, LD pruning, pairwise LD, haplotype blocks |
| VCFtools | 0.1.[x] | KING kinship (`--relatedness2`) |
| R | 4.4 | R scripts |
| R packages | GenABEL [x], GenomicRanges [x] | genotype import, PCA, interval operations |
| BEDTools | [x] | gene annotation (`intersectBed`) |

## Input files

| File | Description |
|---|---|
| `<prefix>.ped/.map` | filtered autosomal SNP genotypes of all sequenced animals (PLINK format) |
| `keep_racing_omn_kw.txt` | FID and IID of the 41 Kuwaiti and Omani racing camels entering QC |
| `racing_project_pheno.txt` | phenotype file for GenABEL: columns `id`, `sex`, `pop` (15 = racing Kuwait; 20 and 23 = racing Oman) |
| `CamDro3.fa.fai` | reference genome index (chromosome lengths) |
| `CamDro3_genes.bed` | gene annotation: chromosome, start, end, gene name |

## Usage

Set the paths in `00_config.sh` (shell scripts) and `project_dir` at the top of each R script, then run the scripts in order.

| Script | Step | Main output |
|---|---|---|
| `00_config.sh` | Paths and software settings | – |
| `01_qc_plink.sh` | Sample call rate (`--mind 0.05`), SNP filters (`--geno 0.1`, `--maf 0.05`, `--hwe 1e-6`), IBS filter (> 0.95), LD pruning (`--indep-pairwise 50 5 0.5`) | 41 → 38 samples; 3,742,192 SNPs (493,441 after pruning) |
| `02_relatedness_king.R` | KING kinship classes (LD-pruned SNPs); each pair counted once | lists of related pairs |
| `03_pca.R` | PCA on LD-pruned SNPs (missing genotypes replaced by the SNP mean; SNPs centred and scaled) | PCA plot (PC1 vs PC2) and coordinates |
| `04a_sample_snps_for_LD.R` | Random 30% of SNPs per autosome for pairwise LD | SNP lists |
| `04b_ld_and_haplotype_blocks.sh` | Pairwise LD (r² ≥ 0.099, ≤ 3 Mb) on the SNP subset; haplotype blocks (`--blocks no-pheno-req no-small-max-span --blocks-max-kb 3000`) on all SNPs | `.ld`, `.blocks.det` |
| `07_lrld_detection.R` | LRLD detection (see definition below) | LRLD events and block BED files per chromosome |
| `08_lrld_catalogue.R` | LRLD event catalogue and coverage of LRLD blocks | 4,444 events; 2,198 blocks; 61.57 Mb (3.18% of autosomes) |
| `09_gene_annotation.sh` | Genes overlapping LRLD blocks (genes between the two blocks of an event are not included) | gene lists |

## LRLD definition

An LRLD event is a pair of distinct haplotype blocks, each containing more than four SNPs, that are

1. separated by at least 1 Mb (from the end of the first block to the start of the second),
2. linked by at least two SNP pairs with r² ≥ 0.6, and
3. supported by at least two distinct SNPs in each block.

Because pairwise LD was calculated within a 3 Mb window, LRLD between blocks more than 3 Mb apart could not be detected.

## Notes

- Pairwise LD was calculated on a random 30% of SNPs per autosome to keep computation manageable; haplotype blocks were estimated from all QC'd SNPs.
- Sizes are in base pairs. The PLINK `.blocks.det` column `KB` is in kilobases and is not used for size calculations.
- Genomic coverage is calculated from the union of non-overlapping intervals, so blocks involved in several LRLD events are counted once.

## Licence

MIT (see `LICENSE`).
