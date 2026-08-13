# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .R
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.19.1
#   kernelspec:
#     display_name: R
#     language: R
#     name: ir
# ---

# %% [markdown]
# # INTERPRET — what the TRS ordering actually means
#
# **Read-only.** Runs after `03_scavenge/scavenge_trs.R`. Answers the two
# questions that stand between the TRS table and any claim you could publish.
#
# **1. What are these clusters?** The reference labels them `Cluster1`…`Cluster31`
# and ships no cluster→cell-type map (it is in Satpathy\*/Granja\* et al. 2019
# Fig 1 / GEO GSE129785). "TRS is highest in Cluster2" is not a biological
# statement. But `colData$Group` records the sample of origin, and in this
# dataset CD34+ sorted samples are separable from BMMC and PBMC — so the
# cluster × Group composition tells you which clusters are progenitor-derived
# without needing anyone's annotation table. That is the exact axis a clonal-
# hematopoiesis hypothesis turns on.
#
# **2. Which loci actually drive it?** g-chromVAR weights a peak by the PP of
# variants inside it. A high-PP variant a few hundred bp outside every peak
# contributes nothing. In this run JAK2's lead sits 7 bp outside a peak and
# TERT's 6.7 kb outside, so the TRS may be carried by the p < 1e-6 loci rather
# than by the genome-wide significant ones. If so that is reportable — but only
# if stated, so measure it rather than assume either way.

# %%
# =============================================================================
# 0. Setup
# =============================================================================
.here <- function(f) {
  cands <- c(file.path(getwd(), f),
             file.path(getwd(), "00_setup", f),
             file.path(getwd(), "..", "00_setup", f),
             file.path(path.expand("~/repos/aou-yp/finemap_scavenge/00_setup"), f))
  hit <- cands[file.exists(cands)]
  if (!length(hit)) stop("cannot find ", f, " — setwd() into finemap_scavenge/ first")
  normalizePath(hit[1])
}
source(.here("config_load.R"))

suppressPackageStartupMessages({
  library(data.table); library(GenomicRanges); library(SummarizedExperiment)
})

log_header(paste0("INTERPRET — cluster identity and locus attribution [", GWAS_NAME, "]"))

need <- file.path(SC_DIR, paste0(c("scavenge_trs_", "trs_by_celltype_"),
                                 paste0("ext", CFG$scavenge$peak_extend_bp[1]), ".tsv"))
for (f in need)
  assert_that(file.exists(f),
              paste0("Missing ", f, " — run 03_scavenge/scavenge_trs.R first."))

PRIM  <- paste0("ext", CFG$scavenge$peak_extend_bp[1])
trs   <- fread(file.path(SC_DIR, paste0("scavenge_trs_", PRIM, ".tsv")))
by_ct <- fread(file.path(SC_DIR, paste0("trs_by_celltype_", PRIM, ".tsv")))
SE    <- readRDS(proj_path(CFG$paths$reference_local))
cd    <- as.data.frame(colData(SE))
if (is.null(rownames(cd))) rownames(cd) <- colnames(SE)

# %%
# =============================================================================
# 1. Cluster identity, inferred from sample of origin
# =============================================================================
assert_that("Group" %in% names(cd),
            paste0("colData has no `Group` column; available: ",
                   paste(names(cd), collapse = ", ")))
trs[, group := cd$Group[match(cell_id, rownames(cd))]]
log_kv("cells matched to a Group", sprintf("%d of %d", sum(!is.na(trs$group)), nrow(trs)))

cat("\n=== Group levels present (sample of origin) ===\n")
print(trs[, .N, by = group][order(-N)])

# This reference's `Group` is far more informative than a CD34 proxy: alongside
# the bulk samples (PBMC_Rep*, Bone_Marrow_Rep*, CD34_Progenitors_Rep*) it
# carries FACS-SORTED populations — HSC, MPP, LMPP, CMP, GMP, MEP, CLP, and the
# mature counterparts. A cluster made mostly of cells from the "HSC" sample is an
# HSC cluster. That is ground truth, not inference, so use it directly.
BULK_RX    <- "^(PBMC|Bone_Marrow|CD34_Progenitors)"
HSPC_SORTS <- c("HSC", "MPP", "LMPP", "CMP", "GMP", "MEP", "CLP")
STEM_SORTS <- c("HSC", "MPP", "LMPP")     # the stem/multipotent end specifically

trs[, base   := sub("_Rep[0-9]+$", "", as.character(group))]
trs[, sorted := !grepl(BULK_RX, group)]   # bulk samples contain everything, so
                                          # they cannot vote on cluster identity
trs[, is_hspc := base %in% HSPC_SORTS | grepl("^CD34_Progenitors", group)]
trs[, is_stem := base %in% STEM_SORTS]

cat("\n=== Sorted populations available to name clusters ===\n")
print(trs[sorted == TRUE, .N, by = base][order(-N)])
log_kv("sorted-population cells", sprintf("%d of %d (%.1f%%)",
       sum(trs$sorted), nrow(trs), 100 * mean(trs$sorted)))

# Cluster identity = the modal sorted population among that cluster's sorted
# cells, with purity. Low purity or few sorted cells means "unresolved", not a
# confident label.
ident <- trs[sorted == TRUE, .N, by = .(celltype, base)][order(celltype, -N)]
ident <- ident[, .(putative_identity = base[1], n_sorted = sum(N),
                   purity_pct = round(100 * N[1] / sum(N), 1),
                   runner_up = if (.N > 1) base[2] else NA_character_), by = celltype]
ident[n_sorted < 20 | purity_pct < 40,
      putative_identity := paste0("unresolved (", putative_identity, "?)")]

frac <- trs[, .(pct_HSPC = round(100 * mean(is_hspc), 1),
                pct_stem_HSC_MPP_LMPP = round(100 * mean(is_stem), 1),
                n_cells = .N), by = celltype]

out <- Reduce(function(a, b) merge(a, b, by = "celltype", all.x = TRUE),
              list(by_ct[, .(celltype, mean_TRS, median_TRS, pct_sig, TRS_vs_overall)],
                   ident, frac))
setorder(out, -mean_TRS)
cat("\n=== TRS ranking with cluster identity from sorted populations ===\n")
print(out)
fwrite(out, file.path(SC_DIR, "trs_by_celltype_annotated.tsv"), sep = "\t")
fwrite(dcast(trs, celltype ~ base, value.var = "cell_id", fun.aggregate = length),
       file.path(SC_DIR, "celltype_by_group.tsv"), sep = "\t")

# Two correlations, because they ask different questions. %HSPC lumps every
# progenitor together; CH should act at the stem/multipotent end specifically,
# so a signal that is real for CH should track HSC/MPP/LMPP more tightly than it
# tracks "any CD34+ progenitor".
for (v in c("pct_HSPC", "pct_stem_HSC_MPP_LMPP")) {
  ct <- suppressWarnings(cor.test(out$mean_TRS, out[[v]], method = "spearman",
                                  exact = FALSE))
  log_kv(paste0("Spearman(mean TRS, ", v, ")"),
         sprintf("rho = %.3f, p = %.3g", ct$estimate, ct$p.value))
}

top <- out[seq_len(min(6, .N))]
log_kv("top TRS clusters", paste0(top$celltype, " = ", top$putative_identity,
                                  " (", top$purity_pct, "% pure, TRS ",
                                  sprintf("%.2f", top$mean_TRS), ")", collapse = "; "))
cat("\n", strwrap(paste0(
  "Read the putative_identity column, not %HSPC alone: CD34+ spans HSC through ",
  "MEP and CLP, and clonal hematopoiesis should enrich at the stem/multipotent ",
  "end rather than across all progenitors. If the top clusters resolve to ",
  "HSC/MPP/LMPP the result is the expected one; if they resolve to a committed ",
  "progenitor or a mature population, that needs explaining before it is ",
  "reported as an HSPC signal. Clusters marked unresolved had too few sorted ",
  "cells to name."), width = 92), "\n", sep = "\n")

# %%
# =============================================================================
# 1b. TRS by SORTED POPULATION — the test that is not confounded
# =============================================================================
# Cluster-level correlation is a weak instrument: n = 31, and every cluster mixes
# populations (the top cluster is only 44.6% pure). The sorted samples are ground
# truth at the cell level, so compare TRS between them directly — hundreds to
# thousands of cells per population instead of 31 points.
lineage <- function(b) {
  fcase(b %in% c("HSC", "MPP", "LMPP"),                 "1_stem",
        b %in% c("CMP", "MEP"),                          "2_early_myeloid_prog",
        b %in% c("GMP", "CLP"),                          "3_committed_prog",
        b %in% c("Monocytes", "Dendritic_Cells", "BM_pDC"), "4_mature_myeloid",
        default =                                        "5_mature_lymphoid")
}
srt <- trs[sorted == TRUE]
srt[, lin := lineage(base)]

by_pop <- srt[, .(n_cells = .N, mean_TRS = round(mean(TRS), 3),
                  median_TRS = round(median(TRS), 3),
                  pct_sig = round(100 * mean(sig, na.rm = TRUE), 1)),
              by = .(lin, base)][order(-mean_TRS)]
cat("\n=== TRS by sorted population (ground truth, per-cell) ===\n"); print(by_pop)
fwrite(by_pop, file.path(SC_DIR, "trs_by_sorted_population.tsv"), sep = "\t")

by_lin <- srt[, .(n_cells = .N, mean_TRS = round(mean(TRS), 3),
                  pct_sig = round(100 * mean(sig, na.rm = TRUE), 1)),
              by = lin][order(lin)]
cat("\n=== TRS by lineage stage ===\n"); print(by_lin)

# Is the stem/early-progenitor compartment actually higher than everything else?
stem_early <- srt[lin %in% c("1_stem", "2_early_myeloid_prog"), TRS]
rest       <- srt[!lin %in% c("1_stem", "2_early_myeloid_prog"), TRS]
w <- wilcox.test(stem_early, rest, alternative = "greater")
log_kv("TRS: HSC/MPP/LMPP/CMP/MEP vs all other sorted cells",
       sprintf("median %.3f vs %.3f, Wilcoxon p = %.3g (n = %d vs %d)",
               median(stem_early), median(rest), w$p.value,
               length(stem_early), length(rest)))

# The specificity question: does the signal survive myeloid commitment? If GMP
# sits with the mature cells rather than with CMP, the enrichment is not
# "progenitors" but specifically the multipotent/common-myeloid compartment.
for (pair in list(c("CMP", "GMP"), c("HSC", "GMP"), c("CMP", "Monocytes"))) {
  a <- srt[base == pair[1], TRS]; b <- srt[base == pair[2], TRS]
  if (length(a) > 5 && length(b) > 5)
    log_kv(paste0("TRS ", pair[1], " vs ", pair[2]),
           sprintf("median %.3f vs %.3f, Wilcoxon p = %.3g",
                   median(a), median(b),
                   wilcox.test(a, b, alternative = "greater")$p.value))
}
log_msg("  > A large CMP-vs-GMP gap means the signal closes on myeloid ",
        "commitment — report it as HSC/CMP-specific, not as 'progenitor ",
        "enrichment', which the %HSPC correlation does not support.")

# %%
# =============================================================================
# 2. Which loci carry the landed weight?
# =============================================================================
src <- CFG$scavenge$trait_source; if (is.null(src)) src <- "primary"
pip_f  <- file.path(FM_DIR, if (identical(src, "primary")) "finemap_pip.tsv"
                            else paste0("finemap_pip_", src, ".tsv"))
peak_f <- file.path(FM_DIR, if (identical(src, "primary")) "peaks.tsv"
                            else paste0("peaks_", src, ".tsv"))
bed_f  <- file.path(SC_DIR, "ch_finemap.PP.hg19.bed")
if (!file.exists(bed_f)) bed_f <- file.path(SC_DIR, "ch_finemap.PP.bed")

bed <- fread(bed_f, col.names = c("chrom", "start", "end", "name", "pp"))
fm  <- fread(pip_f); pk <- fread(peak_f)

# MUST honour scavenge.peak_extend_bp, or this table describes a different
# configuration than the TRS it is explaining. Same helper as the preflight and
# Stage C so all three agree by construction.
ext <- CFG$scavenge$peak_extend_bp[1]
match_gr <- rowRanges(SE)
if (!is.null(ext) && is.finite(ext) && ext > 0) {
  match_gr <- GenomicRanges::resize(match_gr,
                                    width = GenomicRanges::width(match_gr) + 2 * ext,
                                    fix = "center")
}
log_kv("peak extension applied to attribution",
       if (!is.null(ext) && ext > 0) paste0("+/-", ext, " bp") else "none (0)")

ov <- findOverlaps(GRanges(bed$chrom, IRanges(bed$end, bed$end)), match_gr)
bed[, in_peak := FALSE]
if (length(ov)) bed[unique(queryHits(ov)), in_peak := TRUE]
bed[, locus := fm$locus_id[match(name, fm$variant_id)]]

att <- bed[, .(variants = .N, variants_in_peak = sum(in_peak),
               pp_total = round(sum(pp), 4),
               pp_landed = round(sum(pp[in_peak]), 4)), by = locus]
att[, lead_p := pk$lead_p[match(locus, pk$locus_id)]]
att[, genome_wide_sig := lead_p < CFG$peaks$p_threshold]
att[, share_of_landed := round(100 * pp_landed / sum(pp_landed), 1)]
setorder(att, -pp_landed)
cat("\n=== Landed PP attributed to each locus ===\n"); print(att)
fwrite(att, file.path(SC_DIR, "trs_locus_attribution.tsv"), sep = "\t")

gws_share <- 100 * att[genome_wide_sig == TRUE, sum(pp_landed)] / att[, sum(pp_landed)]
log_kv("loci contributing any weight", att[pp_landed > 0, .N])
log_kv("share of landed PP from genome-wide significant loci",
       sprintf("%.1f%%", gws_share))
log_kv(sprintf("share from sub-threshold (p >= %g) loci", CFG$peaks$p_threshold),
       sprintf("%.1f%%", 100 - gws_share))
if (gws_share < 50)
  log_msg("  > ⚠ most of the TRS signal comes from loci that are NOT ",
          "genome-wide significant. That is a legitimate result for a seeding ",
          "analysis, but it must be stated explicitly — the TRS is not a ",
          "readout of the four p < 5e-8 loci.")

cat("\n=== Loci contributing NOTHING (all variants outside every peak) ===\n")
print(att[pp_landed == 0][order(-pp_total)])

log_msg("")
log_msg("**Interpretation complete.** No pipeline output modified.")
for (f in c("trs_by_celltype_annotated.tsv", "celltype_by_compartment.tsv",
            "trs_locus_attribution.tsv"))
  push_to_bucket(file.path(SC_DIR, f), file.path("results/scavenge", f))

# %%
# =============================================================================
# 3. Which locus drives which cell type?
# =============================================================================
# TRS is a per-cell summary over all weighted peaks, so it cannot say whether the
# MEP/HSC signal and the CD8 T signal come from the same loci. They need not:
# g-chromVAR sums PP-weighted accessibility, and two loci whose peaks are open in
# different lineages produce one blended score. Going back to the peaks separates
# them. Accessibility is binarised (scATAC counts are near-binary at 500 bp) and
# averaged within each sorted population.
pop_all  <- sub("_Rep[0-9]+$", "", as.character(cd$Group))
srt_cols <- which(!grepl(BULK_RX, as.character(cd$Group)))

hit_idx <- unique(subjectHits(ov))
if (length(hit_idx)) {
  bin <- assay(SE, "counts")[hit_idx, srt_cols, drop = FALSE] > 0
  grp <- split(seq_along(srt_cols), pop_all[srt_cols])
  acc <- vapply(grp, function(ix) Matrix::rowMeans(bin[, ix, drop = FALSE]),
                numeric(length(hit_idx)))
  rownames(acc) <- as.character(rowRanges(SE)[hit_idx])

  # Which population is each weighted peak most open in, and how specific is it?
  peak_tab <- data.table(
    peak       = rownames(acc),
    top_pop    = colnames(acc)[max.col(acc, ties.method = "first")],
    top_acc    = round(apply(acc, 1, max), 4),
    specificity= round(apply(acc, 1, max) / pmax(apply(acc, 1, mean), 1e-9), 2))
  # attribute each peak to the locus whose variant landed in it, and its weight
  v2p <- data.table(v = queryHits(ov), p = subjectHits(ov))
  v2p[, `:=`(locus = bed$locus[v], pp = bed$pp[v])]
  pk_w <- v2p[, .(pp_in_peak = sum(pp), locus = locus[which.max(pp)]), by = p]
  peak_tab[, idx := hit_idx]
  peak_tab <- merge(peak_tab, pk_w, by.x = "idx", by.y = "p", all.x = TRUE)
  setorder(peak_tab, -pp_in_peak)
  cat("\n=== Weighted peaks: where each is most accessible ===\n")
  print(peak_tab[, .(locus, peak, pp_in_peak, top_pop, top_acc, specificity)])
  fwrite(peak_tab, file.path(SC_DIR, "peak_accessibility_by_population.tsv"), sep = "\t")

  # PP-weighted vote per locus: the population its peaks are most open in.
  loc_tab <- peak_tab[!is.na(locus), .(pp = sum(pp_in_peak), n_peaks = .N),
                      by = .(locus, top_pop)][order(locus, -pp)]
  loc_top <- loc_tab[, .(drives = top_pop[1], pp = round(pp[1], 3),
                         n_peaks = sum(n_peaks)), by = locus][order(-pp)]
  cat("\n=== Which population does each locus point at? ===\n"); print(loc_top)
  fwrite(loc_top, file.path(SC_DIR, "locus_to_population.tsv"), sep = "\t")

  cd8 <- loc_top[grepl("CD8", drives)]
  if (nrow(cd8))
    log_kv("loci whose peaks are most open in CD8 T cells",
           paste0(cd8$locus, " (PP ", cd8$pp, ")", collapse = "; "))
  log_kv("loci pointing at HSPC populations",
         paste0(loc_top[drives %in% c(HSPC_SORTS), locus], collapse = ", "))
  cat("\n", strwrap(paste0(
    "If a small number of loci account for the CD8 T signal while a different ",
    "set drives MEP/HSC/CMP, the second tier is a distinct, locus-specific ",
    "signal rather than a smear of the progenitor one — worth reporting ",
    "separately. If the same peaks are top-accessible in both, the CD8 tier is ",
    "more likely shared regulatory architecture or graph propagation from ",
    "neighbouring cells."), width = 92), "\n", sep = "\n")
  push_to_bucket(file.path(SC_DIR, "locus_to_population.tsv"),
                 "results/scavenge/locus_to_population.tsv")
}
