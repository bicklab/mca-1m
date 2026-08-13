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
# # STAGE C — SCAVENGE TRS for one GWAS, at every peak-extension setting
#
# Reads `$FS_GWAS`, writes to `results/<gwas>/scavenge/`.
#
# **Both extensions always run.** `ext0` is the primary — g-chromVAR's own peak
# assignment, unmodified. `ext250` tolerates a lead falling a few bp outside an
# arbitrary peak boundary. On the JAK2 GWAS the lead missed by 7 bp, so ext0 gave
# a p = 1.3e-49 locus zero weight while ext250 made it the top peak — and the
# cell-type ordering was identical either way. Running both per trait makes that
# robustness check standing output rather than a manual re-run. Accessibility, GC
# bias and background peaks are never altered; extension only changes which peak
# a variant is assigned to.
#
# **Requires `03_scavenge/build_ref_cache.R`.** GC bias, background peaks, the LSI
# embedding and the mutual-kNN graph are trait-independent and come from the
# cache, so every driver is scored on an identical cell graph and TRS rankings
# are comparable across traits.

# %%
# =============================================================================
# 0. Setup
# =============================================================================
.here <- function(f) {
  cands <- c(file.path(getwd(), f), file.path(getwd(), "00_setup", f),
             file.path(getwd(), "..", "00_setup", f),
             file.path(path.expand("~/repos/aou-yp/finemap_scavenge/00_setup"), f))
  hit <- cands[file.exists(cands)]
  if (!length(hit)) stop("cannot find ", f, " — setwd() into finemap_scavenge/ first")
  normalizePath(hit[1])
}
source(.here("config_load.R"))
source(.here("detect_build.R"))

suppressPackageStartupMessages({
  library(data.table); library(Matrix); library(ggplot2)
  library(SummarizedExperiment); library(GenomicRanges); library(GenomeInfoDb)
  library(chromVAR); library(gchromVAR); library(SCAVENGE)
})
SCV        <- CFG$scavenge
EXTENSIONS <- CFG$scavenge$peak_extend_bp
PRIMARY    <- EXTENSIONS[1]

log_header(paste0("STAGE C — SCAVENGE [", GWAS_NAME, "]"))
log_versions(c("SCAVENGE", "gchromVAR", "chromVAR", "SummarizedExperiment"))
log_kv("peak extensions", paste0(paste0(EXTENSIONS, " bp"), collapse = ", "))
log_kv("primary", paste0(PRIMARY, " bp"))

# %%
# =============================================================================
# 1. Reference cache + SE
# =============================================================================
CACHE <- file.path(CACHE_DIR, "scavenge_ref_cache.rds")
assert_that(file.exists(CACHE), paste0(
  "Reference cache missing at ", CACHE,
  "\n  Run once per app:  Rscript 03_scavenge/build_ref_cache.R"))
ref <- readRDS(CACHE)
log_kv("cache", paste0("built ", ref$built_at, " | ", format(ref$n_peaks, big.mark = ","),
                       " peaks x ", format(ref$n_cells, big.mark = ","), " cells"))

SE <- readRDS(proj_path(CFG$paths$reference_local))
if (!"counts" %in% assayNames(SE)) assayNames(SE)[1] <- "counts"
SE <- chromVAR::filterPeaks(SE, non_overlapping = TRUE)
assert_that(nrow(SE) == ref$n_peaks && ncol(SE) == ref$n_cells, paste0(
  "Reference does not match the cache (", nrow(SE), "x", ncol(SE), " vs ",
  ref$n_peaks, "x", ref$n_cells, "). Rebuild: Rscript 03_scavenge/build_ref_cache.R"))
rowData(SE)$bias <- ref$bias

bd <- detect_genome_build(rowRanges(SE), expected = SCV$reference_build,
                          allow_unverified = FALSE)
ref_build <- bd$build
log_kv("reference build", paste0(ref_build, " (via ", bd$how, " — ASSERTED)"))

# %%
# =============================================================================
# 2. Trait BED, build-matched
# =============================================================================
bed_candidates <- c(hg19 = file.path(SC_DIR, "ch_finemap.PP.hg19.bed"),
                    hg38 = file.path(SC_DIR, "ch_finemap.PP.bed"))
trait_bed <- unname(bed_candidates[ref_build])
assert_that(!is.na(trait_bed) && file.exists(trait_bed), paste0(
  "No ", ref_build, " trait BED for ", GWAS_NAME, " at ", trait_bed,
  "\n  Run: FS_GWAS=", GWAS_NAME, " Rscript 02_scavenge_prep/make_trait_bed.R"))
bed_build <- tolower(readLines(paste0(trait_bed, ".build"))[1])
if (bed_build == "grch38") bed_build <- "hg38"
assert_that(identical(bed_build, ref_build),
            paste0("BUILD MISMATCH — BED is ", bed_build, ", reference ", ref_build))
log_kv("trait BED", paste0(basename(trait_bed), " (", bed_build, " — ASSERTED)"))
log_kv("trait variants", nrow(fread(trait_bed, header = FALSE)))

cd_all <- as.data.frame(colData(SE))
if (is.null(rownames(cd_all))) rownames(cd_all) <- colnames(SE)
n_distinct <- function(x) length(unique(x[!is.na(x)]))
ct_col <- SCV$celltype_col
if (is.null(ct_col)) {
  known <- c("BioClassification", "Clusters", "ClusterName", "Group", "celltype",
             "cell_type", "CellType", "Annotation")
  hit <- known[known %in% names(cd_all)]
  cand <- names(cd_all)[vapply(cd_all, function(x)
    (is.character(x) || is.factor(x)) && n_distinct(x) >= 2 && n_distinct(x) <= 60,
    logical(1))]
  ct_col <- if (length(hit)) hit[1] else cand[which.max(vapply(cand, function(k)
    n_distinct(cd_all[[k]]), numeric(1)))]
}
assert_that(!is.na(ct_col) && ct_col %in% names(cd_all),
            paste0("No cell-type column; set scavenge.celltype_col. Have: ",
                   paste(names(cd_all), collapse = ", ")))
log_kv("cell-type column", paste0(ct_col, " (", n_distinct(cd_all[[ct_col]]), " populations)"))

# %%
# =============================================================================
# 3. One SCAVENGE run per extension
# =============================================================================
# For locus attribution: which loci actually carry the landed weight. This is the
# difference between "the TRS reflects your genome-wide significant loci" and "it
# reflects the suggestive ones", and it changes per extension, so it belongs here.
src <- CFG$scavenge$trait_source %||% "primary"
fm_f <- file.path(FM_DIR, if (identical(src, "primary")) "finemap_pip.tsv"
                          else paste0("finemap_pip_", src, ".tsv"))
pk_f <- file.path(FM_DIR, if (identical(src, "primary")) "peaks.tsv"
                          else paste0("peaks_", src, ".tsv"))
FM <- if (file.exists(fm_f)) fread(fm_f) else NULL
PK <- if (file.exists(pk_f)) fread(pk_f) else NULL
BEDT <- fread(trait_bed, col.names = c("chrom", "start", "end", "name", "pp"))

attribute_loci <- function(ext, tag) {
  if (is.null(FM) || is.null(PK)) return(invisible(NULL))
  ov <- findOverlaps(GRanges(BEDT$chrom, IRanges(BEDT$end, BEDT$end)),
                     match_ranges(rowRanges(SE), ext))
  b <- copy(BEDT)
  b[, in_peak := FALSE]
  if (length(ov)) b[unique(queryHits(ov)), in_peak := TRUE]
  b[, locus := FM$locus_id[match(name, FM$variant_id)]]
  att <- b[, .(variants = .N, variants_in_peak = sum(in_peak),
               pp_total = round(sum(pp), 4),
               pp_landed = round(sum(pp[in_peak]), 4)), by = locus]
  att[, lead_p := PK$lead_p[match(locus, PK$locus_id)]]
  att[, genome_wide_sig := lead_p < CFG$peaks$p_threshold]
  att[, share_of_landed := round(100 * pp_landed / max(sum(pp_landed), 1e-12), 1)]
  setorder(att, -pp_landed)
  fwrite(att, file.path(SC_DIR, paste0("locus_attribution_", tag, ".tsv")), sep = "\t")
  gws <- 100 * att[genome_wide_sig == TRUE, sum(pp_landed)] / max(att[, sum(pp_landed)], 1e-12)
  log_kv(paste0("share of landed PP from genome-wide significant loci (", tag, ")"),
         sprintf("%.1f%%", gws))
  zero <- att[pp_landed == 0 & genome_wide_sig == TRUE]
  if (nrow(zero))
    log_kv(paste0("genome-wide significant loci contributing NOTHING (", tag, ")"),
           paste0(zero$locus, " (PP ", zero$pp_total, ")", collapse = "; "))
  invisible(att)
}

match_ranges <- function(gr, ext) {
  if (!is.finite(ext) || ext <= 0) return(gr)
  GenomicRanges::resize(gr, width = GenomicRanges::width(gr) + 2 * ext, fix = "center")
}

run_scavenge <- function(ext) {
  tag <- paste0("ext", ext)
  log_msg(""); log_kv("=== EXTENSION", paste0(ext, " bp ==="))

  trait <- gchromVAR::importBedScore(match_ranges(rowRanges(SE), ext), trait_bed, colidx = 5)
  wts <- as.numeric(SummarizedExperiment::assay(trait))
  n_w <- sum(wts > 0)
  log_kv("peaks carrying trait weight", n_w)
  log_kv("PP landed", sprintf("%.4f", sum(wts)))
  attribute_loci(ext, tag)
  if (n_w == 0) {
    log_msg("  > zero weighted peaks at ", tag, " — skipping this extension.")
    return(NULL)
  }
  if (n_w < 10)
    log_msg("  > only ", n_w, " weighted peak(s); TRS will be dominated by a few ",
            "peaks. Treat as exploratory.")

  dev <- gchromVAR::computeWeightedDeviations(SE, trait, background_peaks = ref$bg)
  z <- as.numeric(SummarizedExperiment::assays(dev)[["z"]]); z[is.na(z)] <- 0
  log_kv("z-score range", sprintf("%.3f to %.3f", min(z), max(z)))

  seed_idx <- as.logical(SCAVENGE::seedindex(z, percent_cut = SCV$seed_percent))
  scale_factor <- SCAVENGE::cal_scalefactor(z_score = z, percent_cut = SCV$scalefactor_percent)
  log_kv("seed cells", sprintf("%d of %d (%.2f%%)", sum(seed_idx), length(seed_idx),
                               100 * mean(seed_idx)))
  assert_that(sum(seed_idx) > 0, "seedindex() returned no seed cells.")

  rw_args <- names(formals(SCAVENGE::randomWalk_sparse))
  np <- if ("queryCells" %in% rw_args)
    SCAVENGE::randomWalk_sparse(intM = ref$knn, queryCells = rownames(ref$knn)[seed_idx],
                                gamma = SCV$rw_gamma)
  else SCAVENGE::randomWalk_sparse(intM = ref$knn, seed_idx = which(seed_idx),
                                   gamma = SCV$rw_gamma)

  keep <- !(Matrix::rowSums(ref$knn) == 0 | np == 0)
  log_kv("degenerate cells dropped", sum(!keep))
  knn_k <- ref$knn[keep, keep]; np_k <- np[keep]
  TRS <- SCAVENGE::max_min_scale(
           SCAVENGE::capOutlierQuantile(x = np_k, q_ceiling = SCV$cap_outlier_quantile)
         ) * scale_factor
  assert_that(sd(TRS) > 0, "TRS is constant — propagation produced no differentiation.")
  log_kv("TRS range", sprintf("%.4f to %.4f (mean %.4f)", min(TRS), max(TRS), mean(TRS)))

  sig_args <- names(formals(SCAVENGE::get_sigcell_simple))
  call_args <- list(knn_sparse_mat = knn_k, seed_idx = seed_idx[keep],
                    topseed_npscore = np_k,
                    permutation_times = SCV$permutation_times,
                    true_cell_significance = SCV$cell_sig_alpha,
                    rda_output = FALSE, rw_gamma = SCV$rw_gamma)
  sig_res <- do.call(SCAVENGE::get_sigcell_simple, call_args[names(call_args) %in% sig_args])
  sig_col <- grep("true_cell", names(sig_res), value = TRUE)[1]
  sig <- if (!is.na(sig_col)) as.logical(sig_res[[sig_col]]) else rep(NA, length(TRS))
  log_kv("significant cells", sprintf("%d of %d (%.1f%%)", sum(sig, na.rm = TRUE),
                                      length(sig), 100 * mean(sig, na.rm = TRUE)))

  cd <- cd_all[keep, , drop = FALSE]
  trs <- data.table(cell_id = rownames(cd), celltype = as.character(cd[[ct_col]]),
                    z_score = z[keep], seed = seed_idx[keep], np_score = np_k,
                    TRS = TRS, sig = sig, extension_bp = ext)
  if ("Group" %in% names(cd)) trs[, group := as.character(cd$Group)]
  fwrite(trs, file.path(SC_DIR, paste0("scavenge_trs_", tag, ".tsv")), sep = "\t")

  by_ct <- trs[, .(n_cells = .N, n_seed = sum(seed), mean_TRS = mean(TRS),
                   median_TRS = median(TRS), sd_TRS = sd(TRS), mean_z = mean(z_score),
                   n_sig = sum(sig, na.rm = TRUE),
                   pct_sig = 100 * mean(sig, na.rm = TRUE)), by = celltype]
  by_ct[, TRS_vs_overall := mean_TRS / mean(trs$TRS)]
  setorder(by_ct, -mean_TRS)
  fwrite(by_ct, file.path(SC_DIR, paste0("trs_by_celltype_", tag, ".tsv")), sep = "\t")

  list(ext = ext, tag = tag, n_weighted_peaks = n_w, pp_landed = sum(wts),
       trs = trs, by_ct = by_ct, cd = cd)
}

RES <- lapply(EXTENSIONS, run_scavenge)
names(RES) <- paste0("ext", EXTENSIONS)
RES <- RES[!vapply(RES, is.null, logical(1))]
assert_that(length(RES) > 0, "No extension produced any weighted peaks.")

# %%
# =============================================================================
# 4. Do the extensions agree? (the standing robustness check)
# =============================================================================
cmp <- rbindlist(lapply(RES, function(r)
  r$by_ct[, .(extension_bp = r$ext, celltype, mean_TRS, pct_sig)]))
fwrite(cmp, file.path(SC_DIR, "trs_by_celltype_all_extensions.tsv"), sep = "\t")

if (length(RES) > 1) {
  w <- dcast(cmp, celltype ~ extension_bp, value.var = "mean_TRS")
  num <- setdiff(names(w), "celltype")
  rho <- cor(w[[num[1]]], w[[num[2]]], method = "spearman", use = "complete.obs")
  log_kv("agreement between extensions (Spearman of cell-type mean TRS)",
         sprintf("rho = %.4f", rho))
  log_kv("top-5 populations per extension",
         paste0(vapply(RES, function(r)
           paste0(r$tag, ": ", paste(r$by_ct$celltype[seq_len(min(5, nrow(r$by_ct)))],
                                     collapse = ">")), character(1)), collapse = "  |  "))
  if (rho < 0.9)
    log_msg("  > the two extensions disagree (rho < 0.9). Peak-assignment tolerance ",
            "is changing the answer, not just its magnitude — inspect before ",
            "reporting either.")
  cat("\n=== mean TRS by cell type, per extension ===\n"); print(w)
}

# %%
# =============================================================================
# 5. Plots — primary extension
# =============================================================================
P <- RES[[paste0("ext", PRIMARY)]]
if (is.null(P)) P <- RES[[1]]
umap_cols <- names(P$cd)[grepl("^(UMAP|umap)", names(P$cd))]

if (length(umap_cols) >= 2) {
  pd <- data.frame(x = P$cd[[umap_cols[1]]], y = P$cd[[umap_cols[2]]], TRS = P$trs$TRS)
  p1 <- ggplot(pd[order(pd$TRS), ], aes(x, y, colour = TRS)) +
    geom_point(size = 0.25, alpha = 0.85) +
    scale_colour_viridis_c(option = "magma") +
    labs(title = paste0("SCAVENGE TRS — ", GWAS_NAME),
         subtitle = paste0(format(nrow(pd), big.mark = ","), " cells | ",
                           P$n_weighted_peaks, " weighted peaks | ", P$tag,
                           " | ", ref_build),
         x = umap_cols[1], y = umap_cols[2]) +
    theme_classic(base_size = 11)
  ggsave(file.path(PLOT_DIR, paste0("scavenge_umap_TRS_", P$tag, ".png")), p1,
         width = 7, height = 6, dpi = 200)
}

ord <- P$by_ct[order(-median_TRS), celltype]
pv <- copy(P$trs)[, celltype := factor(celltype, levels = ord)]
p2 <- ggplot(pv, aes(celltype, TRS, fill = celltype)) +
  geom_violin(scale = "width", linewidth = 0.2, colour = "grey30") +
  geom_boxplot(width = 0.14, outlier.size = 0.2, fill = "white", linewidth = 0.25) +
  labs(title = paste0("TRS by population — ", GWAS_NAME),
       subtitle = paste0(P$tag, " (primary)"), x = NULL, y = "TRS") +
  theme_classic(base_size = 11) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(PLOT_DIR, paste0("scavenge_TRS_by_celltype_", P$tag, ".png")), p2,
       width = max(7, 0.32 * length(ord)), height = 5.5, dpi = 200)

if (length(RES) > 1) {
  w <- dcast(cmp, celltype ~ extension_bp, value.var = "mean_TRS")
  num <- setdiff(names(w), "celltype")
  p3 <- ggplot(w, aes(.data[[num[1]]], .data[[num[2]]])) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey60") +
    geom_point(size = 2, alpha = 0.8) +
    labs(title = paste0("Extension sensitivity — ", GWAS_NAME),
         subtitle = "each point is a cell population; on the line = unaffected",
         x = paste0("mean TRS (ext", num[1], ")"),
         y = paste0("mean TRS (ext", num[2], ")")) +
    theme_classic(base_size = 11)
  ggsave(file.path(PLOT_DIR, "scavenge_extension_sensitivity.png"), p3,
         width = 5.5, height = 5.5, dpi = 200)
}
log_kv("plots", PLOT_DIR)

for (f in list.files(SC_DIR, pattern = "[.]tsv$", full.names = TRUE))
  push_to_bucket(f, file.path("scavenge", basename(f)))
for (f in list.files(PLOT_DIR, pattern = "^scavenge_.*png$", full.names = TRUE))
  push_to_bucket(f, file.path("plots", basename(f)))

log_msg("")
log_msg("**STAGE C complete for ", GWAS_NAME, "** (extensions: ",
        paste(EXTENSIONS, collapse = ", "), " bp). Next: `04_report/make_report.R`.")
