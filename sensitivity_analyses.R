#!/usr/bin/env Rscript
# Sensitivity analyses requested in Reviewer #1 comment 1(c) and reported in the
# response letter. Three checks:
#
#  A. Opposite-arm negative control. For each significant cis locus, test the same
#     variant against the SAME mCA type on the OPPOSITE arm of the same chromosome.
#     A true cis effect should be null there. Observed P range: 0.23-0.9.
#
#  B. 18q= phenotype refinement. Restrict the 18q= phenotype to events that overlap
#     (or nearly overlap) rs1496805 at PIK3C3. If the signal is real the effect
#     should strengthen; if it is a false positive it should attenuate.
#     Observed: beta 2.48 -> 4.57.
#
#  C. 18q= excluding 1KG-EUR-like participants. rs1496805 is essentially
#     monomorphic in EUR (AF 0.001), so EUR contributes no signal; removing it
#     reduces scope for stratification confounding.
#     Observed: beta = 3.17, SE = 0.50, P = 2e-10.

suppressPackageStartupMessages({
  library(data.table)
})

leads <- fread("results/cis_lead_variants.tsv")

# ---- A. opposite-arm control ------------------------------------------------
opposite_arm <- function(label) {
  if (grepl("p", label)) sub("p", "q", label) else sub("q", "p", label)
}

ctrl <- rbindlist(lapply(seq_len(nrow(leads)), function(i) {
  lab <- leads$mca_label[i]
  opp <- opposite_arm(lab)
  f <- sprintf("results/cis/META_%s_1.tbl", opp)
  if (!file.exists(f)) return(NULL)
  d <- fread(f)
  hit <- d[MarkerName == leads$rsid[i]]
  if (nrow(hit) == 0L) return(NULL)
  data.table(mca_label = lab, control_pheno = opp, rsid = leads$rsid[i],
             beta = hit$Effect, se = hit$StdErr, p = hit$`P-value`)
}))
fwrite(ctrl, "results/cis_opposite_arm_control.tsv", sep = "\t")

# ---- B. 18q= refined to events overlapping rs1496805 ------------------------
# Refine the phenotype, then re-run run_regenie_cis.sh with the refined pheno file.
PIK3C3_POS <- 41315963L
PAD <- 1e6L   # "nearly overlapping"

calls <- fread("data/mca_calls_pooled.tsv")
e18q <- calls[chrom == "chr18" & arm == "q" & type == "="]
e18q[, overlaps_pik3c3 := start <= (PIK3C3_POS + PAD) & end >= (PIK3C3_POS - PAD)]
cat(sprintf("18q= events: %d total, %d overlapping rs1496805 +/- 1Mb\n",
            nrow(e18q), sum(e18q$overlaps_pik3c3)))

pheno <- fread("pheno/chr18-cnloh-q.txt")
pheno[, mca_refined := as.integer(IID %in% e18q[overlaps_pik3c3 == TRUE]$sample_id)]
# controls stay the same; carriers of non-overlapping 18q= are dropped, not
# recoded as controls
pheno <- pheno[!(IID %in% e18q[overlaps_pik3c3 == FALSE]$sample_id)]
fwrite(pheno, "pheno/chr18-cnloh-q_refined.txt", sep = "\t")

# ---- C. 18q= excluding 1KG-EUR-like -----------------------------------------
keep <- fread("data/ancestry_assignments.tsv")[ancestry != "EUR", .(FID, IID)]
fwrite(keep, "keep_18q_nonEUR.id", sep = "\t", col.names = FALSE)

cat("Re-run run_regenie_cis.sh with:\n",
    "  --phenoColList mca_refined  (analysis B)\n",
    "  --keep keep_18q_nonEUR.id   (analysis C)\n")
