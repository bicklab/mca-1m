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
# # STAGE Q — QQ plot and genomic-control lambda, restricted to the cis arm
#
# One trait per run (`FS_GWAS=<name>`), exactly like every other stage. Writes
# `results/<name>/qc/qq_lambda.tsv` (a single row) and
# `results/<name>/plots/qq_<name>.png`. `05_cis_mca/make_supplement.R` stacks the
# per-trait rows into the lambda table.
#
# **What "restricted to the cis arm" means, and why it is checked rather than
# assumed.** These summary statistics are supposed to be cis to the mCA — only
# variants near the event, which in practice means the same chromosome arm. If
# that is true, the restriction removes nothing and the run says so. If it is not
# true, lambda is being computed over a mixture of a region under selection and a
# region that is not, and it means nothing. The `frac_on_named_arm` column is the
# number that settles it: 1.0 = the file was already cis.
#
# **Two lambdas per stratum, deliberately.** `lambda` comes from (BETA/SE)^2,
# which is exact under METAL's `SCHEME STDERR`; `lambda_p` comes from the P
# column. They should agree to ~3 decimals. When they do not, the P column was
# not derived from Effect/StdErr — worth knowing before you quote either.
#
# **Egress:** variant-level summary statistics only. No participant data.

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
source(.here("config_load.R"))     # -> CFG, FM_DIR, PLOT_DIR, SUMSTATS, CIS, ...
source(.here("cis_arms.R"))
source(.here("load_sumstats.R"))

suppressPackageStartupMessages(library(data.table))

QC_DIR <- file.path(GWAS_OUT, "qc")
dir.create(QC_DIR, recursive = TRUE, showWarnings = FALSE)

log_header(paste0("STAGE Q — QQ + lambda_GC [", GWAS_NAME, "]"))
log_versions(c("data.table", "qqman", "yaml"))
log_kv("GWAS build", CFG$project$gwas_build)

have_qqman <- requireNamespace("qqman", quietly = TRUE)
if (!have_qqman)
  log_msg("  > ⚠ qqman is not installed — falling back to a base-R QQ panel. ",
          "Install it with `Rscript 00_setup/install_packages.R`.")

# %%
# =============================================================================
# 1. Load, restrict to the arm, measure
# =============================================================================
cen <- load_centromeres(CFG$project$gwas_build,
                        cytoband_path = CFG$cis$cytoband_path,
                        log = log_msg)

if (is.null(CIS)) {
  log_msg("  > No cis metadata for '", GWAS_NAME, "' — running genome-wide. ",
          "Add a row to cis_manifest.tsv (05_cis_mca/run_cis_mca.sh writes it) ",
          "or name the trait like `<set>_chr9_cnloh_p_meta_1` to get the arm ",
          "restriction.")
} else {
  log_kv("cis trait", sprintf("chr%s %s %s-arm  [set=%s, meta=%s]",
                              CIS$chrom, CIS$event, CIS$arm, CIS$set, CIS$meta_tag))
}

ld <- load_sumstats(SUMSTATS, CFG, cis = CIS, cen = cen, keep_pvals = TRUE)
dt <- ld$dt; st <- ld$stats

# %%
# =============================================================================
# 2. QQ plot
# =============================================================================
# Two panels: the whole cis arm, and the analysis universe that Stage A will
# actually fine-map (arm ∩ tested in >= qc.min_cohorts cohorts). Showing only the
# second would hide the single-cohort variants that deflate lambda; showing only
# the first would not describe what gets fine-mapped.
MAX_POINTS <- as.numeric(CFG$qc$qq_max_points %||% 3e6)

# qqman::qq() plots every point and derives the expected quantiles from the
# vector it is handed, so thinning the input would move the expected axis and
# distort the line. Above MAX_POINTS we therefore draw the QQ ourselves, with the
# expected quantiles computed from the FULL n and only the plotting thinned.
qq_manual <- function(p, main, keep_top = 5e4, n_thin = 2e5) {
  p <- sort(p[is.finite(p) & p > 0 & p <= 1])
  n <- length(p)
  obs <- -log10(p)
  exp <- -log10(ppoints(n))                      # from the full n, not the sample
  idx <- unique(c(seq_len(min(keep_top, n)),
                  round(seq(min(keep_top, n) + 1, n, length.out = min(n_thin, n)))))
  idx <- idx[idx >= 1 & idx <= n]
  plot(exp[idx], obs[idx], pch = 19, cex = 0.4, col = "grey30",
       xlab = expression(Expected~~-log[10](italic(p))),
       ylab = expression(Observed~~-log[10](italic(p))), main = main)
  abline(0, 1, col = "red")
  invisible(length(idx))
}

draw_qq <- function(p, main) {
  p <- p[is.finite(p) & p > 0 & p <= 1]
  if (!length(p)) { plot.new(); title(paste0(main, "\n(no variants)")); return("empty") }
  if (have_qqman && length(p) <= MAX_POINTS) {
    qqman::qq(p, main = main, pch = 19, cex = 0.4, col = "grey30")
    "qqman::qq"
  } else {
    qq_manual(p, main)
    if (have_qqman) "base-R (thinned; over qc.qq_max_points)" else "base-R (qqman absent)"
  }
}

png_path <- file.path(PLOT_DIR, paste0("qq_", GWAS_NAME, ".png"))
png(png_path, width = 1500, height = 780, res = 140)
op <- par(mfrow = c(1, 2), mar = c(4.5, 4.5, 4, 1.5))

drawn_a <- draw_qq(st$p_arm, sprintf("%s\ncis arm: all variants (n = %s)",
                                     GWAS_NAME, format(length(st$p_arm), big.mark = ",")))
mtext(bquote(lambda[GC] == .(sprintf("%.3f", st$lambda_arm))), side = 3, line = -1.4,
      adj = 0.03, cex = 0.85)

drawn_b <- draw_qq(dt$pval, sprintf("analysis universe: >= %g cohort(s) (n = %s)",
                                    st$min_cohorts, format(nrow(dt), big.mark = ",")))
mtext(bquote(lambda[GC] == .(sprintf("%.3f", st$lambda))), side = 3, line = -1.4,
      adj = 0.03, cex = 0.85)
abline(h = -log10(CFG$peaks$p_threshold), lty = 3, col = "steelblue")

par(op); dev.off()
log_kv("QQ plot", paste0(png_path, "  [", drawn_a, " | ", drawn_b, "]"))
push_to_bucket(png_path, file.path("results/plots", basename(png_path)))

# %%
# =============================================================================
# 3. The one-row summary
# =============================================================================
# Every count that went into the lambdas is carried alongside them, so the table
# is auditable without reopening the .tbl.
num <- function(x, d = NA_real_) if (is.null(x) || !length(x)) d else as.numeric(x)
chr_str <- function(x, d = NA_character_) if (is.null(x) || !length(x)) d else as.character(x)

row <- data.table(
  gwas             = GWAS_NAME,
  set              = chr_str(CIS$set),
  chrom            = chr_str(CIS$chrom),
  event            = chr_str(CIS$event),
  arm              = chr_str(CIS$arm),
  meta_tag         = chr_str(CIS$meta_tag),
  arm_basis        = chr_str(st$arm_basis),
  file             = st$file,
  n_rows_file      = num(st$n_rows_file),
  n_after_qc       = num(st$n_after_qc),
  n_off_chrom      = num(st$n_off_chrom, 0),
  n_p_arm          = num(st$n_p_arm),
  n_acen           = num(st$n_acen),
  n_q_arm          = num(st$n_q_arm),
  frac_on_named_arm = num(st$frac_on_named_arm),
  n_removed_by_cis = num(st$n_dropped, 0),
  n_on_arm         = num(st$n_on_arm),
  max_cohorts      = num(st$max_cohorts),
  min_cohorts      = num(st$min_cohorts),
  n_removed_by_cohort_filter = num(st$n_dropped_universe, 0),
  n_analysed       = num(st$n_analysed),
  lambda_file      = num(st$lambda_file),
  lambda_file_p    = num(st$lambda_file_p),
  lambda_arm       = num(st$lambda_arm),
  lambda_arm_p     = num(st$lambda_arm_p),
  lambda           = num(st$lambda),          # headline: arm ∩ analysis universe
  lambda_p         = num(st$lambda_p),
  min_p            = num(st$min_p),
  n_gws            = num(st$n_gws),
  p_threshold      = num(CFG$peaks$p_threshold),
  qq_png           = basename(png_path)
)

out <- file.path(QC_DIR, "qq_lambda.tsv")
fwrite(row, out, sep = "\t")
push_to_bucket(out, file.path("results/qc", "qq_lambda.tsv"))

if (!is.null(st$cohort_support))
  fwrite(st$cohort_support, file.path(QC_DIR, "cohort_support.tsv"), sep = "\t")

log_kv("lambda_GC (headline, cis arm ∩ analysis universe)", sprintf("%.4f", st$lambda))
log_kv("summary row", out)

# %%
# =============================================================================
# 4. Verification
# =============================================================================
# lambda is the median chi-square over its null expectation. Recompute it from a
# completely different route — the empirical median p-value — as a check that the
# vector being plotted is the vector being summarised.
log_header("STAGE Q — verification")
med_p  <- median(dt$pval)
lam_chk <- qchisq(med_p, df = 1, lower.tail = FALSE) / qchisq(0.5, df = 1)
assert_that(abs(lam_chk - st$lambda_p) < 1e-6,
            paste0("lambda from the median p (", lam_chk,
                   ") disagrees with lambda_p (", st$lambda_p, ")"))
log_kv("(a) lambda_p reproduced from the median p-value",
       sprintf("PASS (median p = %.4g -> lambda = %.4f)", med_p, lam_chk))

if (!is.null(CIS) && !startsWith(st$arm_basis, "whole_chromosome")) {
  a <- assign_arm(dt$chr, dt$pos, cen)
  assert_that(all(a == CIS$arm, na.rm = TRUE) && !anyNA(a),
              paste0("After the cis restriction, ", sum(a != CIS$arm | is.na(a)),
                     " variant(s) are not on chr", CIS$chrom, CIS$arm))
  log_kv("(b) every analysed variant is on the named arm",
         paste0("PASS (", format(nrow(dt), big.mark = ","), " variants on chr",
                CIS$chrom, CIS$arm, ")"))
} else {
  log_kv("(b) arm check", paste0("n/a — basis is '", st$arm_basis, "'"))
}

log_msg("")
log_msg("**STAGE Q complete.** Next: `01_finemap/finemap_abf.R`.")
