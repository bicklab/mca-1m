#!/usr/bin/env Rscript
# Label-permutation enrichment test ("Quantitative Enrichment Analysis").
# Reviewer #2 asked specifically for this code and for clarification of the null.
#
# Null model: gene coordinates and mCA region coordinates are held FIXED; only the
# oncogene / tumour-suppressor LABELS are shuffled across protein-coding genes.
# Because positions never move, gene-dense loci (LILR cluster, MHC) contribute
# identically under the null and cannot by themselves generate enrichment.
#
# Statistic: number of recurrently GAINED regions containing >= 1 proto-oncogene
#            number of recurrently LOST  regions containing >= 1 tumour suppressor
# Each region contributes at most one gene per class (presence/absence), so
# gene-dense regions are not over-counted.
#
# MHC (chr6:28,510,120-33,480,577, GRCh38) is excluded.
#
# Reproduces the pooled and per-cohort p-values in the response to Reviewer #2:
#   gains / proto-oncogenes  pooled 0.003  (UKB 0.008, AoU 0.07, BioVU 0.09, TOPMed 0.01)
#   losses / tumour suppress pooled 0.001  (UKB 0.003, AoU 0.04, BioVU 0.02, TOPMed 0.06)

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
})

set.seed(20250101)
N_PERM <- 10000L

MHC <- GRanges("chr6", IRanges(28510120, 33480577))  # GRCh38

# --- inputs -----------------------------------------------------------------
# genes.tsv : chrom, start, end, symbol, class in {oncogene, tsg, both, neither}
#             (OncoKB annotation restricted to protein-coding genes)
# msar_regions.tsv : output of msar_define.R at the selected threshold, plus
#             a `direction` column in {gain, loss, cnloh} and `cohort`
#             ("pooled" for regions defined across all four cohorts)
genes   <- fread("data/oncokb_protein_coding_genes.tsv")
regions <- fread("results/msar_regions_selected.tsv")

gr_genes <- GRanges(genes$chrom, IRanges(genes$start, genes$end))
genes <- genes[!overlapsAny(gr_genes, MHC)]          # drop MHC
gr_genes <- GRanges(genes$chrom, IRanges(genes$start, genes$end))

# --- statistic --------------------------------------------------------------
# labels: character vector aligned to `genes`, values "oncogene"/"tsg"/"both"/"neither"
count_regions_with_class <- function(labels, regions_sub, target) {
  if (nrow(regions_sub) == 0L) return(0L)
  gr_reg <- GRanges(regions_sub$chrom,
                    IRanges(regions_sub$start_hg38, regions_sub$end_hg38))
  hit_genes <- labels %in% c(target, "both")
  if (!any(hit_genes)) return(0L)
  # presence/absence per region: at most one gene per class per region
  sum(countOverlaps(gr_reg, gr_genes[hit_genes]) > 0L)
}

run_test <- function(regions_sub, direction, target) {
  rs <- regions_sub[regions_sub$direction == direction, ]
  observed <- count_regions_with_class(genes$class, rs, target)
  null <- vapply(seq_len(N_PERM), function(i) {
    count_regions_with_class(sample(genes$class), rs, target)
  }, integer(1))
  # empirical p = fraction of permutations with enrichment >= observed
  # (+1/+1 correction so p is never exactly 0)
  p <- (sum(null >= observed) + 1) / (N_PERM + 1)
  data.table(direction = direction, gene_class = target,
             n_regions = nrow(rs), observed = observed,
             null_mean = mean(null), p = p)
}

out <- rbindlist(lapply(unique(regions$cohort), function(ch) {
  rs <- regions[cohort == ch]
  rbind(
    cbind(cohort = ch, run_test(rs, "gain", "oncogene")),
    cbind(cohort = ch, run_test(rs, "loss", "tsg"))
  )
}))

fwrite(out, "results/msar_permutation_enrichment.tsv", sep = "\t")
print(out)
