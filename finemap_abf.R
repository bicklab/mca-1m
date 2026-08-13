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
# # STAGE A — LD-free fine-mapping of the CH GWAS meta-analysis (Wakefield/Maller aBF)
#
# Input: METAL fixed-effects (`SCHEME STDERR`) `.tbl`, GRCh38.
# Output: per-locus PIPs, 95%/99% credible sets, sentinels, QC panels.
#
# **Why aBF and nothing else** — aBF needs only marginal BETA and SE, so it is the
# only principled single-causal fine-mapper for a multi-cohort meta-analysis where
# no single LD panel is correct. SuSiE/FINEMAP/CARMA are deliberately absent
# (Appendix B). *No LD matrix is read anywhere in this script.*
#
# **Egress:** this stage reads variant-level summary statistics only — no
# participant-level data, so nothing here is subject to small-cell suppression.
# Everything written is an aggregate over the meta-analysis.
#
# **Memory.** This meta-analysis is ~58 M variants / 3.6 GiB on disk, which is
# roughly 7–8 GB once loaded (the `MarkerName` strings dominate). Peak RSS during
# coordinate parsing is around 12 GB, so give the app **≥ 32 GB RAM** for Stage A.
# The load path reads only the six columns it uses and renames in place rather
# than copying; the log records the in-memory size at two checkpoints so you can
# see where you actually landed.
#
# Run top-to-bottom in the Workbench R (`ir`) app after
# `00_setup/stage_inputs.sh sumstats` has staged the METAL `.tbl`.

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
source(.here("config_load.R"))          # -> CFG, PROJ_ROOT, DATA_DIR, RESULTS_DIR, set.seed()
source(.here("cis_arms.R"))             # -> load_centromeres(), assign_arm(), cis_restrict()
source(.here("load_sumstats.R"))        # -> load_sumstats(), lambda_gc()

suppressPackageStartupMessages({
  library(data.table)
})

# FM_DIR / PLOT_DIR / SUMSTATS come from config_load.R and are namespaced by
# $FS_GWAS, so one script serves every CHIP driver.
log_header(paste0("STAGE A — aBF fine-mapping [", GWAS_NAME, "]"))
log_versions(c("data.table", "locuszoomr", "ggplot2", "yaml"))
log_kv("seed", CFG$project$seed)
log_kv("GWAS build", CFG$project$gwas_build)

OMEGA     <- CFG$finemap$omega
WINDOW_BP <- CFG$finemap$window_kb   * 1000
SPAN_BP   <- CFG$peaks$peak_flank_kb * 1000
P_THR     <- CFG$peaks$p_threshold
MIN_PTS   <- CFG$peaks$min_variants_per_peak

# %%
# =============================================================================
# 1-2. Load, QC, and define the analysis universe
# =============================================================================
# The whole of this — column mapping, coordinate parsing, dedupe, lambda_GC, the
# optional cis-arm restriction, and the cohort-support filter — lives in
# 00_setup/load_sumstats.R so that Stage Q's QQ plot and this stage's fine-mapping
# describe byte-for-byte the same set of variants. A QQ plot drawn over a
# different universe than the one that gets fine-mapped only looks like a check.
#
# CIS is NULL unless this trait has cis metadata (cis_manifest.tsv, or a name like
# `<set>_chr9_cnloh_p_meta_1`), in which case the analysis is restricted to that
# chromosome arm and the log records how much the restriction actually removed.
cen <- if (!is.null(CIS))
  load_centromeres(CFG$project$gwas_build, CFG$cis$cytoband_path, log = log_msg) else NULL
if (!is.null(CIS))
  log_kv("cis trait", sprintf("chr%s %s %s-arm  (metadata from %s)",
                              CIS$chrom, CIS$event, CIS$arm, CIS$source))

.ld    <- load_sumstats(SUMSTATS, CFG, cis = CIS, cen = cen)
dt     <- .ld$dt
lambda <- .ld$stats$lambda
log_kv("in-memory size of the analysis table", format(object.size(dt), units = "auto"))


# %%
# =============================================================================
# 3-6. One analysis arm: peak calling -> aBF -> credible sets -> outputs
# =============================================================================
# Wrapped in a function because we run it twice on the SAME loaded data. The
# 3.6 GiB read dominates the runtime; peak calling and fine-mapping are seconds,
# so a second arm is nearly free.
#
#   primary arm  (p < 5e-8)  -> the reportable fine-mapping. Untouched.
#   scavenge arm (p < 1e-6)  -> seeds SCAVENGE only. At 5e-8 the trait BED held
#                               18 variants hitting 2 of 571,400 reference peaks,
#                               against ~1.7 expected by chance -- too thin for
#                               network propagation to mean anything. Suffixed
#                               outputs so the two never collide.
#
# A campaign over many cis regions is mostly null by construction, so "no variant
# reaches the threshold" has to be an outcome rather than a crash — otherwise the
# null regions read as run failures and bury the ones with signal. Controlled by
# peaks.allow_no_loci, which stays false for a single genome-wide GWAS where an
# empty result really is a mistake worth halting on.
ALLOW_NO_LOCI <- isTRUE(CFG$peaks$allow_no_loci)

# Empty-but-well-formed outputs, so every downstream reader sees the same schema
# whether or not the trait had a locus.
EMPTY_TABLES <- list(
  peaks = c("locus_id", "sentinel_id", "chr", "lead_pos", "lead_p",
            "n_sig_in_span", "near_neighbour"),
  finemap_pip = c("locus_id", "variant_id", "chr", "pos", "ea", "nea", "beta", "se",
                  "pval", "log_aBF", "aBF", "pip", "in_cs95", "in_cs99",
                  "direction", "n_studies"),
  credible_sets = c("locus_id", "chr", "sentinel_id", "lead_pos", "lead_p",
                    "max_pip_variant", "max_pip_pos", "max_pip", "max_pip_p",
                    "sentinel_is_max_pip", "cs95_size", "cs99_size", "cs95_sum_pip",
                    "cs99_sum_pip", "n_variants_in_window", "sum_pip",
                    "n_pass_ppmin", "pp_pass_ppmin"),
  sentinels = c("locus_id", "sentinel_id", "chr", "pos", "pval",
                "max_pip_variant", "max_pip", "cs95_size"))

write_no_loci <- function(suffix, why) {
  log_kv("no loci", why)
  for (nm in names(EMPTY_TABLES)) {
    f <- file.path(FM_DIR, paste0(nm, suffix, ".tsv"))
    writeLines(paste(EMPTY_TABLES[[nm]], collapse = "\t"), f)
    push_to_bucket(f, file.path("results/finemap", basename(f)))
  }
  writeLines(paste0("# no variant reached p < ", format(p_thr_last, scientific = TRUE)),
             file.path(FM_DIR, paste0("NO_LOCI", suffix, ".txt")))
  invisible(NULL)
}
p_thr_last <- NA_real_    # set by run_arm so write_no_loci can name the threshold

run_arm <- function(p_thr, suffix = "") {
  p_thr_last <<- p_thr
  log_msg("")
  log_kv("=== ARM", paste0(if (nzchar(suffix)) sub("^_", "", suffix) else "primary",
                           " (p < ", format(p_thr, scientific = TRUE), ") ==="))
# %%
  # =============================================================================
  # 3. Peak calling -> provisional sentinels
  # =============================================================================
  # locuszoomr::quick_peak(): greedy — take the smallest-p variant, drop everything
  # on the same chromosome within +/- `span`, repeat. `span` is a HALF-width, so
  # peak_flank_kb maps to span directly. See docs/verified_api_notes.md.
  #
  # Only p < P_THR rows can ever be peaks, so we hand it just those (identical
  # result, orders of magnitude less memory). Note quick_peak() takes a
  # *data.frame* — a data.table would mis-resolve its internal column lookups.
  #
  # The `rsid` column is not decoration. quick_peak() forwards to locuszoomr's
  # internal detect_cols(), which ALWAYS autodetects a SNP-id column even though
  # quick_peak exposes no argument for one — it needs exactly one column matching
  # /rs?id/i (or /snp/i) and otherwise dies with "unable to autodetect SNP id
  # column". None of chr/pos/pval match that, so a three-column frame fails.
  sig_df <- as.data.frame(dt[pval < p_thr, .(rsid = variant_id, chr, pos, pval)])
  if (nrow(sig_df) == 0) {
    if (ALLOW_NO_LOCI)
      return(write_no_loci(suffix, paste0("no variant reaches p < ",
                                          format(p_thr, scientific = TRUE),
                                          " (min p = ", signif(min(dt$pval), 3), ")")))
    assert_that(FALSE,
                paste0("No variant reaches p < ", p_thr,
                       " — nothing to fine-map. Loosen peaks.p_threshold to explore, ",
                       "or set peaks.allow_no_loci: true if a null region is expected."))
  }

  # Dependency-free equivalent, used if locuszoomr is unavailable (restricted egress).
  # Difference from upstream: min_points is applied to EVERY peak including the
  # first (upstream skips the check for peak 1), and a region failing the check is
  # consumed rather than re-scanned. Both paths are harmonised by the filter below.
  quick_peak_fallback <- function(df, p_cutoff, span, min_points) {
    idx <- order(df$pval)
    idx <- idx[df$pval[idx] < p_cutoff]
    pks <- integer(0)
    while (length(idx)) {
      top  <- idx[1]
      near <- which(df$chr[idx] == df$chr[top] & abs(df$pos[idx] - df$pos[top]) <= span)
      if (length(near) >= min_points) pks <- c(pks, top)
      idx <- idx[-near]
    }
    pks
  }

  # Prefer the published implementation, but never let it stop the run: its column
  # autodetection has already changed once between versions, and the fallback is a
  # faithful reimplementation of the same greedy rule.
  peak_idx <- NULL; peak_caller <- NULL
  if (requireNamespace("locuszoomr", quietly = TRUE)) {
    peak_idx <- tryCatch(
      locuszoomr::quick_peak(sig_df, npeaks = NA, p_cutoff = p_thr,
                             span = SPAN_BP, min_points = MIN_PTS,
                             chrom = "chr", pos = "pos", p = "pval"),
      error = function(e) {
        log_msg("  > locuszoomr::quick_peak failed (", conditionMessage(e),
                ") — using the built-in equivalent instead.")
        NULL
      })
    if (!is.null(peak_idx))
      peak_caller <- paste0("locuszoomr::quick_peak ",
                            as.character(utils::packageVersion("locuszoomr")))
  }
  if (is.null(peak_idx)) {
    peak_idx <- quick_peak_fallback(sig_df, p_thr, SPAN_BP, MIN_PTS)
    peak_caller <- if (requireNamespace("locuszoomr", quietly = TRUE))
      "built-in quick_peak_fallback (locuszoomr call errored)" else
      "built-in quick_peak_fallback (locuszoomr not installed)"
  }
  log_kv("peak caller", peak_caller)

  peaks <- data.table(sig_df[peak_idx, ])
  peaks[, rsid := NULL]                    # only existed to satisfy detect_cols()
  setnames(peaks, c("chr", "pos", "pval"), c("chr", "lead_pos", "lead_p"))

  # Re-apply min_variants_per_peak uniformly (upstream exempts the first peak).
  peaks[, n_sig_in_span := vapply(seq_len(.N), function(i)
    sum(sig_df$chr == chr[i] & abs(sig_df$pos - lead_pos[i]) <= SPAN_BP), integer(1))]
  dropped <- peaks[n_sig_in_span < MIN_PTS]
  if (nrow(dropped))
    log_kv("peaks dropped for < min_variants_per_peak", nrow(dropped))
  peaks <- peaks[n_sig_in_span >= MIN_PTS]
  if (nrow(peaks) == 0) {
    if (ALLOW_NO_LOCI)
      return(write_no_loci(suffix, paste0(
        "every candidate peak failed min_variants_per_peak = ", MIN_PTS,
        " (", nrow(sig_df), " variant(s) below the threshold, none with a ",
        "supporting neighbour within ", SPAN_BP / 1000, " kb)")))
    assert_that(FALSE, "Every candidate peak failed min_variants_per_peak.")
  }

  setorder(peaks, chr, lead_pos)
  peaks[, chr_n := suppressWarnings(as.integer(chr))]
  peaks[chr == "X", chr_n := 23L][chr == "Y", chr_n := 24L][chr == "MT", chr_n := 25L]
  setorder(peaks, chr_n, lead_pos)
  peaks[, locus_id := sprintf("L%03d_chr%s_%d", .I, chr, lead_pos)]

  # Attach the lead's identity/effect by exact coordinate match. Multi-allelic sites
  # survive the earlier dedupe (same chr:pos, different alleles), so collapse to the
  # strongest record per coordinate first — otherwise the merge duplicates the locus.
  # Join only the lead coordinates rather than sorting all 58 M rows: dt[order(...)]
  # would copy the entire table just to look up a few dozen positions.
  lead_lookup <- dt[peaks[, .(chr, pos = lead_pos)], on = .(chr, pos), nomatch = 0L,
                    .(chr, pos, pval, sentinel_id = variant_id,
                      sentinel_ea = ea, sentinel_nea = nea)]
  setorder(lead_lookup, chr, pos, pval)
  lead_lookup <- unique(lead_lookup, by = c("chr", "pos"))[, pval := NULL][]
  peaks <- merge(peaks, lead_lookup,
                 by.x = c("chr", "lead_pos"), by.y = c("chr", "pos"),
                 all.x = TRUE, sort = FALSE)
  setorder(peaks, chr_n, lead_pos)
  assert_that(!any(is.na(peaks$sentinel_id)),
              "A peak lead could not be matched back to a variant record.")

  # Broad regions (HLA) can surface as several leads < 2*span apart — flag, don't merge.
  peaks[, near_neighbour := FALSE]
  for (c_ in unique(peaks$chr)) {
    p_ <- peaks[chr == c_, lead_pos]
    if (length(p_) > 1) {
      d <- as.matrix(dist(p_)); diag(d) <- Inf
      peaks[chr == c_, near_neighbour := apply(d, 1, min) < 2 * SPAN_BP]
    }
  }
  if (any(peaks$near_neighbour))
    log_msg("  > ⚠ ", sum(peaks$near_neighbour), " lead(s) sit < ", 2 * SPAN_BP / 1000,
            " kb from another lead (broad-region artefact, e.g. HLA) — see peaks.tsv.")

  log_kv("N loci", nrow(peaks))
  fwrite(peaks[, .(locus_id, sentinel_id, chr, lead_pos, lead_p, n_sig_in_span,
                   near_neighbour)],
         file.path(FM_DIR, paste0("peaks", suffix, ".tsv")), sep = "\t")
  print(peaks[, .(locus_id, sentinel_id, chr, lead_pos, lead_p, n_sig_in_span)])

# %%
  # =============================================================================
  # 4. Approximate Bayes factors -> PIPs
  # =============================================================================
  # Wakefield / Maller et al. (2012), with V = SE^2 and W = omega:
  #     aBF = sqrt( V / (V+W) ) * exp( W * BETA^2 / (2 * V * (V+W)) )
  # Computed on the log scale — the exponent routinely exceeds 700 at a strong
  # lead, which overflows the direct form to Inf and destroys the PIPs.
  log_abf <- function(beta, se, omega) {
    V <- se^2
    0.5 * log(V / (V + omega)) + omega * beta^2 / (2 * V * (V + omega))
  }

  # Per-locus softmax: PIP_i = aBF_i / sum_j aBF_j, done as exp(l - max) / sum(...).
  setkey(dt, chr, pos)                     # binary search on chr, then scan in-window
  finemap_one <- function(i, omega, peaks) {
    lead <- peaks[i]
    w <- dt[.(lead$chr), nomatch = 0L][
      pos >= lead$lead_pos - WINDOW_BP & pos <= lead$lead_pos + WINDOW_BP]
    if (!nrow(w)) return(NULL)
    la  <- log_abf(w$beta, w$se, omega)
    pip <- exp(la - max(la)); pip <- pip / sum(pip)
    out <- data.table(locus_id = lead$locus_id, variant_id = w$variant_id,
                      chr = w$chr, pos = w$pos, ea = w$ea, nea = w$nea,
                      beta = w$beta, se = w$se, pval = w$pval,
                      log_aBF = la, aBF = exp(pmin(la, 700)), pip = pip)
    # METAL's Direction is one character per contributing cohort ('+'/'-'/'?').
    # A credible-set variant carried by one cohort out of five is a very different
    # claim from one carried by all five, and aBF cannot see the difference —
    # it only sees BETA and SE. Computed here, on the small window, rather than
    # over all 58 M rows.
    if ("direction" %in% names(w)) {
      out[, direction := w$direction]
      out[, n_plus    := nchar(gsub("[^+]", "", w$direction))]
      out[, n_minus   := nchar(gsub("[^-]", "", w$direction))]
      out[, n_missing := nchar(gsub("[^?]", "", w$direction))]
      out[, n_studies := n_plus + n_minus]
      # Every contributing cohort agreeing on the sign is the cheapest
      # heterogeneity statement there is, and it is what a supplemental table
      # gets asked for. aBF cannot see it — it reads only BETA and SE.
      out[, dir_concordant := (n_plus == 0L | n_minus == 0L) & n_studies > 0L]
    }
    # Optional passthroughs (effect-allele frequency, per-variant weight) when the
    # METAL run carried them. Attached here, on the small window, rather than
    # being propagated through every intermediate.
    for (opt in intersect(c("freq", "freq_se", "n"), names(w)))
      set(out, j = opt, value = w[[opt]])
    out
  }
# %%
  # =============================================================================
  # 5. Credible sets — computed at every omega, so the prior's influence is visible
  # =============================================================================
  # Smallest set, by descending PIP, whose cumulative PIP first reaches the level.
  mark_cs <- function(pip, level) {
    o  <- order(-pip)
    k  <- which(cumsum(pip[o]) >= level)[1]
    if (is.na(k)) k <- length(pip)          # cannot happen (sum == 1) but be safe
    out <- logical(length(pip)); out[o[seq_len(k)]] <- TRUE; out
  }
  lv <- CFG$finemap$cs_levels

  # omega = 0.04 asserts a prior SD of 0.2 on the effect — a convention from
  # COMMON-variant case-control GWAS. This meta-analysis has rare large-effect
  # leads (|beta| up to 2.6), which that prior calls implausible and shrinks hard,
  # leaving their aBF below the aBF mass a 500 kb window accrues for free. Running
  # a second omega does not change which answer is primary; it documents whether a
  # credible set is driven by the data or by the prior. Loci whose CS is stable
  # across omega are data-driven; loci that collapse from ~12,000 to a handful were
  # prior-limited. See docs/verified_api_notes.md and diagnose_universe.R.
  finemap_at_omega <- function(omega, peaks) {
    f <- rbindlist(lapply(seq_len(nrow(peaks)), finemap_one, omega = omega, peaks = peaks))
    f[, in_cs95 := mark_cs(pip, lv[1]), by = locus_id]
    f[, in_cs99 := mark_cs(pip, lv[2]), by = locus_id]
    setorder(f, locus_id, -pip)
    f[]
  }

  summarise_cs <- function(f, peaks) {
    s <- f[, {
      top <- .SD[which.max(pip)]
      .(max_pip_variant = top$variant_id, max_pip_chr = top$chr, max_pip_pos = top$pos,
        max_pip = top$pip, max_pip_p = top$pval,
        cs95_size = sum(in_cs95), cs99_size = sum(in_cs99),
        cs95_sum_pip = sum(pip[in_cs95]), cs99_sum_pip = sum(pip[in_cs99]),
        n_variants_in_window = .N, sum_pip = sum(pip),
        n_pass_ppmin = sum(in_cs95 & pip > CFG$finemap$pp_min_export),
        pp_pass_ppmin = sum(pip[in_cs95 & pip > CFG$finemap$pp_min_export]))
    }, by = locus_id]
    s <- merge(peaks[, .(locus_id, sentinel_id, chr, lead_pos, lead_p)], s,
               by = "locus_id", sort = FALSE)
    s[, sentinel_is_max_pip := sentinel_id == max_pip_variant][]
  }

  OMEGAS <- unique(c(OMEGA, as.numeric(CFG$finemap$omega_sensitivity)))
  OMEGAS <- OMEGAS[is.finite(OMEGAS)]
  log_kv("omega (primary)", OMEGA)
  log_kv("omega (sensitivity)", paste(setdiff(OMEGAS, OMEGA), collapse = ", "))

  FM <- lapply(OMEGAS, finemap_at_omega, peaks = peaks); names(FM) <- sprintf("%g", OMEGAS)
  CS <- lapply(FM, summarise_cs, peaks = peaks)

  fm <- FM[[sprintf("%g", OMEGA)]]          # primary result, used from here on
  cs <- CS[[sprintf("%g", OMEGA)]]

  log_kv("variants in fine-mapping windows (locus-variant pairs)",
         format(nrow(fm), big.mark = ","))
  log_kv("unique variants fine-mapped", format(uniqueN(fm$variant_id), big.mark = ","))
  if (any(!cs$sentinel_is_max_pip))
    log_msg("  > note: ", sum(!cs$sentinel_is_max_pip), " of ", nrow(cs),
            " loci have a max-PIP variant different from the min-P sentinel ",
            "(expected — aBF weights SE as well as effect size).")

  log_kv("median CS95 size", median(cs$cs95_size))
  log_kv("median CS99 size", median(cs$cs99_size))
  log_kv("loci with a single-variant CS95", sum(cs$cs95_size == 1))

  # --- side-by-side across omega ------------------------------------------------
  sens <- rbindlist(lapply(names(CS), function(k)
    CS[[k]][, .(omega = as.numeric(k), locus_id, chr, lead_pos, lead_p,
                max_pip_variant, max_pip, cs95_size, cs99_size,
                n_pass_ppmin, pp_pass_ppmin)]))
  setorder(sens, locus_id, omega)
  fwrite(sens, file.path(FM_DIR, paste0("credible_sets_omega_sensitivity", suffix, ".tsv")), sep = "\t")
  cat("\n=== Credible sets across omega ===\n")
  print(dcast(sens, locus_id + chr + lead_pos ~ omega,
              value.var = c("max_pip", "cs95_size", "pp_pass_ppmin")))

  prior_limited <- sens[, .(stable = uniqueN(max_pip_variant) == 1,
                            cs_ratio = max(cs95_size) / pmax(min(cs95_size), 1)),
                        by = locus_id][cs_ratio > 10]
  if (nrow(prior_limited))
    log_msg("  > ⚠ prior-limited loci (CS95 size changes >10x across omega): ",
            paste(prior_limited$locus_id, collapse = ", "),
            ". Their credible sets at omega = ", OMEGA,
            " reflect the prior more than the data — report them as unresolved, ",
            "or adopt the larger omega with that stated explicitly.")

# %%
  # =============================================================================
  # 6. Outputs
  # =============================================================================
  out_cols <- c("locus_id", "variant_id", "chr", "pos", "ea", "nea",
                intersect(c("freq", "freq_se", "n"), names(fm)),
                "beta", "se", "pval", "log_aBF", "aBF", "pip", "in_cs95", "in_cs99",
                intersect(c("direction", "n_plus", "n_minus", "n_missing",
                            "n_studies", "dir_concordant"), names(fm)))
  fwrite(fm[, ..out_cols], file.path(FM_DIR, paste0("finemap_pip", suffix, ".tsv")), sep = "\t")

  # Sensitivity runs get their own files. Stage B reads finemap_pip.tsv, so the
  # primary omega remains the one that feeds SCAVENGE unless you repoint it.
  for (k in setdiff(names(FM), sprintf("%g", OMEGA))) {
    f <- file.path(FM_DIR, sprintf("finemap_pip_omega%s%s.tsv", k, suffix))
    fwrite(FM[[k]][, ..out_cols], f, sep = "\t")
    log_kv(paste0("sensitivity output (omega = ", k, ")"), basename(f))
  }

  # Cohort support across the 95% credible sets — thin support is not an error, but
  # it changes how much weight a credible set deserves downstream.
  if ("n_studies" %in% names(fm)) {
    supp <- fm[in_cs95 == TRUE, .N, by = n_studies][order(n_studies)]
    log_kv("CS95 variants by contributing cohorts",
           paste0(supp$n_studies, " cohort(s): ", supp$N, collapse = "; "))
    thin <- fm[in_cs95 == TRUE & n_studies <= 1]
    if (nrow(thin))
      log_msg("  > ⚠ ", nrow(thin), " CS95 variant(s) are supported by a single cohort",
              " (max PIP ", sprintf("%.3f", max(thin$pip)), "). aBF sees only BETA/SE,",
              " so it cannot down-weight these — check them before Stage C.")
  }

  fwrite(cs[, .(locus_id, chr, sentinel_id, lead_pos, lead_p,
                max_pip_variant, max_pip_pos, max_pip, max_pip_p,
                sentinel_is_max_pip, cs95_size, cs99_size,
                cs95_sum_pip, cs99_sum_pip, n_variants_in_window, sum_pip,
                n_pass_ppmin, pp_pass_ppmin)],
         file.path(FM_DIR, paste0("credible_sets", suffix, ".tsv")), sep = "\t")

  fwrite(cs[, .(locus_id, sentinel_id, chr, pos = lead_pos, pval = lead_p,
                max_pip_variant, max_pip, cs95_size)],
         file.path(FM_DIR, paste0("sentinels", suffix, ".tsv")), sep = "\t")

  # --- optional VEP -------------------------------------------------------------
  # The image is very unlikely to ship the vep binary + its offline cache, and
  # egress blocks fetching them. If absent we still export a VCF of CS95 variants:
  # chr/pos/alleles only, no participant data, so it is safe to annotate outside.
  vep_bin <- Sys.which("vep")[[1]]
  cs95_vars <- unique(fm[in_cs95 == TRUE, .(chr, pos, variant_id, ea, nea)])
  vcf_path <- file.path(FM_DIR, paste0("credible_sets_cs95", suffix, ".vcf"))
  writeLines(c("##fileformat=VCFv4.2",
               "##note=CS95 variants from aBF fine-mapping; REF=Allele2(nea), ALT=Allele1(ea).",
               "##note=METAL does not record which allele is reference — VEP will",
               "##note=report REF mismatches for any variant where this guess is wrong.",
               "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"), vcf_path)
  fwrite(cs95_vars[order(chr, pos), .(chr, pos, variant_id, nea, ea,
                                      QUAL = ".", FILTER = ".", INFO = ".")],
         vcf_path, sep = "\t", col.names = FALSE, append = TRUE)

  run_vep <- CFG$finemap$run_vep
  if (identical(run_vep, "auto") && nzchar(vep_bin) || isTRUE(run_vep)) {
    vep_out <- file.path(FM_DIR, paste0("finemap_vep", suffix, ".tsv"))
    st <- system2(vep_bin, c("-i", shQuote(vcf_path), "-o", shQuote(vep_out),
                             "--tab", "--cache", "--offline", "--assembly", "GRCh38",
                             "--symbol", "--canonical", "--force_overwrite"))
    log_kv("VEP", if (st == 0) paste("ran ->", basename(vep_out))
                  else paste("FAILED (exit", st, ") — annotate", basename(vcf_path), "externally"))
  } else {
    log_kv("VEP", paste0("skipped — `vep` not on PATH. Annotate ",
                         basename(vcf_path), " with the Ensembl web VEP (GRCh38) ",
                         "outside the workbench; it contains no participant data."))
  }

  for (f in c(paste0("peaks", suffix, ".tsv"), paste0("finemap_pip", suffix, ".tsv"), paste0("credible_sets", suffix, ".tsv"), paste0("sentinels", suffix, ".tsv"),
              paste0("credible_sets_cs95", suffix, ".vcf"), paste0("credible_sets_omega_sensitivity", suffix, ".tsv")))
    push_to_bucket(file.path(FM_DIR, f), file.path("results/finemap", f))
  log_kv("outputs", paste0(FM_DIR, " (synced to ", gcs_out("results/finemap"), ")"))

  invisible(list(peaks = peaks, fm = fm, cs = cs, sens = sens))
}

# %%
# =============================================================================
# Run the arms. Verification and QC panels below operate on the PRIMARY arm.
# =============================================================================
primary <- run_arm(P_THR, "")

P_THR_SCAV <- CFG$peaks$p_threshold_scavenge
if (!is.null(P_THR_SCAV) && is.finite(P_THR_SCAV) && P_THR_SCAV != P_THR) {
  scav <- run_arm(P_THR_SCAV, "_scavenge")
  if (!is.null(scav))
    log_msg("  > SCAVENGE arm written with the `_scavenge` suffix. Point Stage B at it ",
            "with scavenge.trait_source: scavenge in config.yaml.")
}

# Everything below verifies and plots the PRIMARY arm. With no primary locus
# there is nothing to verify — the empty tables are already written, so stop here
# with a SUCCESS status rather than failing on an absent object.
#
# quit() rather than wrapping the rest of the file in `else { ... }`: this script
# is also opened as a notebook via jupytext, and `# %%` cell markers inside a
# brace block produce cells that are not independently runnable.
if (is.null(primary)) {
  log_msg("")
  log_msg("**STAGE A complete — no locus at p < ", format(P_THR, scientific = TRUE),
          ".** Empty result tables written; SCAVENGE will be skipped for this trait.")
  if (!interactive()) quit(status = 0, save = "no")
  stop("No locus at p < ", format(P_THR, scientific = TRUE),
       " — nothing below this point applies. (Interactive session: stopping here.)",
       call. = FALSE)
}

peaks <- primary$peaks; fm <- primary$fm; cs <- primary$cs

# %%
# =============================================================================
# 7. VERIFICATION — the checks the brief asks for, run as hard assertions
# =============================================================================
log_header("STAGE A — verification")

# (a) PIPs sum to 1 within every locus.
sums <- fm[, .(s = sum(pip)), by = locus_id]
worst <- sums[which.max(abs(s - 1))]
assert_that(all(abs(sums$s - 1) < 1e-8),
            paste0("PIPs do not sum to 1 in locus ", worst$locus_id, " (", worst$s, ")"))
log_kv("(a) sum(PIP) == 1 in all loci", sprintf("PASS (max |dev| = %.2e)",
                                                max(abs(sums$s - 1))))

# (b) Hand-recompute the top-PIP variant of the strongest locus, scalar and
#     on the NATURAL scale, from BETA/SE only — independent of the vectorised path.
best_locus <- cs[which.min(lead_p), locus_id]
sub <- fm[locus_id == best_locus]
hand_abf_scalar <- function(b, s, w) {
  V <- s * s
  sqrt(V / (V + w)) * exp(w * b * b / (2 * V * (V + w)))
}
hand <- vapply(seq_len(nrow(sub)),
               function(i) hand_abf_scalar(sub$beta[i], sub$se[i], OMEGA), numeric(1))
if (any(!is.finite(hand))) {
  log_kv("(b) hand-check", paste0("natural-scale aBF overflowed to Inf in ",
                                  sum(!is.finite(hand)), " variant(s) — this is exactly ",
                                  "why the pipeline works in logs; comparing on the log scale instead"))
  hand_log <- vapply(seq_len(nrow(sub)), function(i) {
    V <- sub$se[i]^2
    0.5 * log(V / (V + OMEGA)) + OMEGA * sub$beta[i]^2 / (2 * V * (V + OMEGA))
  }, numeric(1))
  hand_pip <- exp(hand_log - max(hand_log)); hand_pip <- hand_pip / sum(hand_pip)
} else {
  hand_pip <- hand / sum(hand)
}
top_i <- which.max(hand_pip)
cat("Locus ", best_locus, " — hand-computed top variant:\n", sep = "")
cat("  variant   : ", sub$variant_id[top_i], "\n",
    "  beta / se : ", signif(sub$beta[top_i], 6), " / ", signif(sub$se[top_i], 6), "\n",
    "  PIP (hand): ", signif(hand_pip[top_i], 10), "\n",
    "  PIP (pipe): ", signif(sub$pip[top_i], 10), "\n", sep = "")
assert_that(identical(sub$variant_id[top_i], sub$variant_id[which.max(sub$pip)]),
            "Hand-computed and vectorised top-PIP variants disagree.")
assert_that(isTRUE(all.equal(hand_pip, sub$pip, tolerance = 1e-9)),
            "Hand-computed and vectorised PIP vectors disagree beyond 1e-9.")
log_kv("(b) hand recomputation of top-PIP variant",
       paste0("PASS (", best_locus, ", ", sub$variant_id[top_i],
              ", PIP = ", signif(sub$pip[top_i], 4), ")"))

# (c) A locus with one strong signal should concentrate PIP on its lead.
log_kv("(c) max PIP at the strongest locus", sprintf("%.4f (CS95 size %d)",
       max(sub$pip), cs[locus_id == best_locus, cs95_size]))
if (max(sub$pip) < 0.1)
  log_msg("  > ⚠ lead PIP < 0.1 even at the strongest locus — either the window is ",
          "dense with near-equivalent variants (normal for a common-variant locus) ",
          "or omega is mis-specified. Inspect the QC panel before trusting the CS.")

# (d) Windows overlap when two leads are < 2*window apart, so a variant can carry
#     a PIP in more than one locus. Reported, not silently collapsed.
dupes <- fm[, .N, by = variant_id][N > 1]
log_kv("variants appearing in >1 locus window", nrow(dupes))

# %%
# =============================================================================
# 8. QC panels — P and PIP against position, per top locus
# =============================================================================
# Base graphics on purpose: no ggplot2/patchwork dependency, and these are
# throwaway diagnostics rather than figures for the paper.
top_loci <- cs[order(lead_p)][seq_len(min(CFG$finemap$qc_plot_top_n, nrow(cs))), locus_id]

# Prefer real LocusZoom figures (gene track underneath), fall back to the
# self-contained base-R panels. locuszoomr needs an EnsDb for the gene track, and
# app rebuilds have a habit of losing large annotation packages, so a run must
# not fail for want of a plot.
ENSDB <- c("EnsDb.Hsapiens.v86", "EnsDb.Hsapiens.v79")
ensdb <- ENSDB[vapply(ENSDB, requireNamespace, logical(1), quietly = TRUE)][1]
use_lz <- requireNamespace("locuszoomr", quietly = TRUE) && !is.na(ensdb)
log_kv("locus plots", if (use_lz) paste0("locuszoomr + ", ensdb, " (gene track)")
                      else "base-R panels (locuszoomr/EnsDb unavailable)")

plot_locus_lz <- function(d, lead, lid) {
  df <- as.data.frame(d[, .(rsid = variant_id, chr = as.integer(sub("X", "23", chr)),
                            pos, pval, pip)])
  df$chr[is.na(df$chr)] <- 23L
  loc <- locuszoomr::locus(data = df, seqname = df$chr[1],
                           xrange = c(min(df$pos), max(df$pos)),
                           ens_db = ensdb, chrom = "chr", pos = "pos",
                           p = "pval", labs = "rsid")
  png(file.path(PLOT_DIR, paste0("locuszoom_", lid, ".png")),
      width = 1200, height = 950, res = 130)
  locuszoomr::locus_plot(loc, labels = "index",
                         main = paste0(GWAS_NAME, " — ", lid))
  dev.off()
}

plot_locus_base <- function(d, lead, lid) {
  png(file.path(PLOT_DIR, paste0("locuszoom_", lid, ".png")),
      width = 1100, height = 800, res = 130)
  op <- par(mfrow = c(2, 1), mar = c(4, 4.5, 2.5, 1))
  plot(d$pos / 1e6, -log10(d$pval), pch = 19, cex = 0.45, col = "grey60",
       xlab = "", ylab = expression(-log[10](P)),
       main = paste0(GWAS_NAME, " — ", lid, "  (", nrow(d),
                     " variants, CS95 = ", lead$cs95_size, ")"))
  points(d[in_cs95 == TRUE, pos] / 1e6, -log10(d[in_cs95 == TRUE, pval]),
         pch = 19, cex = 0.6, col = "#1f77b4")
  abline(h = -log10(P_THR), lty = 3, col = "grey40")
  points(lead$lead_pos / 1e6, -log10(lead$lead_p), pch = 23, cex = 1.3,
         bg = "#d62728", col = "black")
  plot(d$pos / 1e6, d$pip, pch = 19, cex = 0.45, col = "grey60",
       xlab = paste0("chr", lead$chr, " position (Mb)"), ylab = "PIP (aBF)")
  points(d[in_cs95 == TRUE, pos] / 1e6, d[in_cs95 == TRUE, pip],
         pch = 19, cex = 0.6, col = "#1f77b4")
  points(lead$max_pip_pos / 1e6, lead$max_pip, pch = 23, cex = 1.3,
         bg = "#2ca02c", col = "black")
  legend("topright", bty = "n", cex = 0.75,
         legend = c("in 95% CS", "min-P sentinel", "max-PIP variant"),
         pch = c(19, 23, 23), col = c("#1f77b4", "black", "black"),
         pt.bg = c(NA, "#d62728", "#2ca02c"))
  par(op); dev.off()
}

for (lid in top_loci) {
  d <- fm[locus_id == lid][order(pos)]
  lead <- cs[locus_id == lid]
  ok <- FALSE
  if (use_lz) ok <- tryCatch({ plot_locus_lz(d, lead, lid); TRUE },
                             error = function(e) {
                               log_msg("  > locuszoomr failed on ", lid, " (",
                                       conditionMessage(e), ") — base-R panel instead")
                               FALSE })
  if (!ok) plot_locus_base(d, lead, lid)
}
log_kv("QC panels", paste0(length(top_loci), " PNGs in ", PLOT_DIR))
for (lid in top_loci)
  push_to_bucket(file.path(PLOT_DIR, paste0("locuszoom_", lid, ".png")),
                 file.path("results/plots", paste0("locuszoom_", lid, ".png")))

log_msg("")
log_msg("**STAGE A complete.** Next: `02_scavenge_prep/make_trait_bed.R`.")
