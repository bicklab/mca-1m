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
# # STAGE B — fine-mapping output → SCAVENGE trait BED (+ build reconciliation)
#
# g-chromVAR ingests a BED whose **column 5 is the posterior probability**. We keep
# variants with `pip > 0.001` that sit inside a 95% credible set (the paper's
# "PP > 0.1% within the 95% CS" retention rule).
#
# **Build handling.** Our GWAS is GRCh38; the SCAVENGE hematopoiesis reference is
# hg19 (`docs/verified_api_notes.md` §6). Rather than guessing which BED Stage C
# will need, this script writes **both** and stamps each with a sidecar
# `.build` file. Stage C detects the reference's real build from `rowRanges(SE)`
# and *asserts* that it matches the sidecar of the BED it picks. Builds are never
# mixed silently, and neither stage has to trust the other's assumption.

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

# GenomicRanges/rtracklayer are needed only for the liftOver block, so they are
# loaded there — the BED itself must still be producible on a box where the
# Bioconductor stack has not been installed yet.
suppressPackageStartupMessages(library(data.table))

log_header(paste0("STAGE B — trait BED for SCAVENGE [", GWAS_NAME, "]"))
log_versions(c("data.table", "GenomicRanges", "rtracklayer", "IRanges"))

# %%
# =============================================================================
# 1. Select the variants that go downstream
# =============================================================================
# Stage A writes two arms. The reportable fine-mapping is always the primary one;
# which arm SEEDS SCAVENGE is a separate choice, because the trait BED needs
# genomic reach rather than genome-wide significance. See scavenge.trait_source.
src <- CFG$scavenge$trait_source
if (is.null(src)) src <- "primary"
pip_path <- file.path(FM_DIR, if (identical(src, "primary")) "finemap_pip.tsv"
                              else paste0("finemap_pip_", src, ".tsv"))
assert_that(file.exists(pip_path), paste0(
  "Missing ", pip_path, "\n  scavenge.trait_source = '", src, "'.",
  "\n  Run 01_finemap/finemap_abf.R (the '", src, "' arm needs ",
  "peaks.p_threshold_scavenge set), or switch trait_source to 'primary'."))
fm <- fread(pip_path)
log_kv("trait source arm", paste0(src, "  (", basename(pip_path), ")"))
log_kv("fine-mapped locus-variant pairs read", format(nrow(fm), big.mark = ","))

PP_MIN <- CFG$finemap$pp_min_export
sel <- fm[in_cs95 == TRUE & pip > PP_MIN]
log_kv("retained (in CS95 & PP > 0.1%)", format(nrow(sel), big.mark = ","))

# Overlapping +/-250 kb windows can hand one variant a PIP in two loci. g-chromVAR
# weights a peak by the PP of the variants inside it, so a duplicated row would
# double-count. Keep the strongest PP per physical variant.
setorder(sel, chr, pos, -pip)
n_pre <- nrow(sel)
sel <- sel[!duplicated(paste(chr, pos, ea, nea))]
log_kv("collapsed duplicates across overlapping windows", n_pre - nrow(sel))
log_kv("unique variants in trait BED", format(nrow(sel), big.mark = ","))
log_kv("sum of PP across trait BED", sprintf("%.3f", sum(sel$pip)))

assert_that(nrow(sel) > 0,
            "No variant survived the CS95 / PP>0.001 filter — nothing to feed SCAVENGE.")

# %%
# =============================================================================
# 2. Write the GRCh38 BED (PP in column 5)
# =============================================================================
# UCSC "chr" naming throughout: the reference peaks, the BSgenome, and the UCSC
# chain file all use it, so adopting it here avoids a seqlevels mismatch that
# would silently yield zero overlaps in importBedScore().
ucsc_chr <- function(x) {
  x <- sub("^chr", "", toupper(as.character(x)))
  x[x %in% c("MT", "M")] <- "M"
  paste0("chr", x)
}

bed38 <- data.table(
  chrom      = ucsc_chr(sel$chr),
  chromStart = as.integer(sel$pos) - 1L,   # BED is 0-based half-open
  chromEnd   = as.integer(sel$pos),
  name       = sel$variant_id,
  score      = sel$pip                     # <- column 5 = posterior probability
)
setorder(bed38, chrom, chromStart)

bed38_path <- file.path(SC_DIR, "ch_finemap.PP.bed")
fwrite(bed38, bed38_path, sep = "\t", col.names = FALSE, scipen = 999)
writeLines("GRCh38", paste0(bed38_path, ".build"))

# Canonical sort via bedtools when available (the spec asks for it); the R sort
# above is already deterministic, so this is belt-and-braces.
if (nzchar(Sys.which("bedtools"))) {
  tmp <- tempfile(fileext = ".bed")
  st <- system2("bedtools", c("sort", "-i", shQuote(bed38_path)), stdout = tmp)
  if (st == 0 && file.size(tmp) > 0) file.copy(tmp, bed38_path, overwrite = TRUE)
  log_kv("bedtools sort", if (st == 0) "applied" else paste("failed (exit", st, ")"))
} else {
  log_kv("bedtools sort", "bedtools not on PATH — used the R-side sort")
}
log_kv("GRCh38 BED", bed38_path)

# %%
# =============================================================================
# 3. liftOver GRCh38 -> hg19
# =============================================================================
# Produced unconditionally (when the chain is staged) so Stage C can pick whichever
# build the reference actually turns out to be, without re-running this stage.
chain_path <- proj_path(CFG$paths$chain_local)
bed19_path <- file.path(SC_DIR, "ch_finemap.PP.hg19.bed")

have_bioc <- requireNamespace("rtracklayer", quietly = TRUE) &&
             requireNamespace("GenomicRanges", quietly = TRUE)

if (!file.exists(chain_path) || !have_bioc) {
  why <- if (!have_bioc) "rtracklayer/GenomicRanges not installed"
         else paste0("chain file not staged at `", chain_path, "`")
  log_msg("- **liftOver**: SKIPPED — ", why, ".")
  log_msg("  > Run `00_setup/install_packages.R`, then `stage_inputs.sh push` ",
          "(needs internet) and `pull`.")
  log_msg("  > ⚠ If the Stage C reference is hg19, Stage C will HALT without the hg19 BED.")
} else {
  suppressPackageStartupMessages({
    library(GenomicRanges); library(rtracklayer)
  })
  chain <- import.chain(chain_path)

  gr38 <- GRanges(seqnames = bed38$chrom,
                  ranges   = IRanges(start = bed38$chromEnd, end = bed38$chromEnd),
                  name     = bed38$name,
                  score    = bed38$score)

  lifted  <- liftOver(gr38, chain)
  n_hits  <- lengths(lifted)

  n_unmapped   <- sum(n_hits == 0)
  n_multimap   <- sum(n_hits > 1)
  keep_multi   <- !isTRUE(CFG$liftover$drop_multimapped)
  keep_idx     <- if (keep_multi) which(n_hits >= 1) else which(n_hits == 1)

  gr19 <- unlist(lifted[keep_idx])
  # A variant that lands on a different chromosome is a mapping artefact, not a
  # coordinate update — drop it rather than carry a wrong locus into SCAVENGE.
  src_chrom <- rep(as.character(seqnames(gr38))[keep_idx], n_hits[keep_idx])
  wrong_chr <- as.character(seqnames(gr19)) != src_chrom
  n_wrongchr <- sum(wrong_chr)
  gr19 <- gr19[!wrong_chr]

  bed19 <- data.table(
    chrom      = as.character(seqnames(gr19)),
    chromStart = start(gr19) - 1L,
    chromEnd   = start(gr19),
    name       = gr19$name,
    score      = gr19$score
  )
  setorder(bed19, chrom, chromStart)
  fwrite(bed19, bed19_path, sep = "\t", col.names = FALSE, scipen = 999)
  writeLines("hg19", paste0(bed19_path, ".build"))

  drop_frac <- 1 - nrow(bed19) / nrow(bed38)
  log_kv("liftOver chain", basename(chain_path))
  log_kv("input variants (GRCh38)", nrow(bed38))
  log_kv("unmapped (dropped)", n_unmapped)
  log_kv("multi-mapped", paste0(n_multimap, if (keep_multi) " (kept)" else " (dropped)"))
  log_kv("mapped to a different chromosome (dropped)", n_wrongchr)
  log_kv("output variants (hg19)", nrow(bed19))
  log_kv("fraction lost", sprintf("%.3f%%", 100 * drop_frac))
  log_kv("PP retained after liftOver",
         sprintf("%.3f of %.3f", sum(bed19$score), sum(bed38$score)))
  if (drop_frac > CFG$liftover$max_drop_frac_warn)
    log_msg("  > ⚠ liftOver lost more than ",
            100 * CFG$liftover$max_drop_frac_warn, "% of credible-set variants. ",
            "Check whether any dropped variant carried a large PP — a lost causal ",
            "candidate would blunt the TRS signal.")

  lost <- bed38[!name %in% bed19$name][order(-score)]
  if (nrow(lost)) {
    fwrite(lost, file.path(SC_DIR, "liftover_dropped.tsv"), sep = "\t")
    log_kv("largest PP dropped by liftOver", sprintf("%.4f (%s)",
                                                     lost$score[1], lost$name[1]))
  }
  log_kv("hg19 BED", bed19_path)
}

# %%
# =============================================================================
# 4. Checks + sync
# =============================================================================
log_header("STAGE B — verification")

check_bed <- function(path) {
  if (!file.exists(path)) return(invisible(NULL))
  b <- fread(path, header = FALSE)
  nm <- basename(path)
  assert_that(ncol(b) >= 5, paste0(nm, ": fewer than 5 columns — g-chromVAR reads PP from col 5."))
  assert_that(is.numeric(b$V5), paste0(nm, ": column 5 is not numeric."))
  assert_that(all(b$V5 > 0 & b$V5 <= 1), paste0(nm, ": column 5 outside (0, 1]."))
  assert_that(all(b$V3 - b$V2 == 1L), paste0(nm, ": not all intervals are 1 bp."))
  assert_that(all(grepl("^chr", b$V1)), paste0(nm, ": seqnames are not UCSC-style."))
  log_kv(paste0("check ", nm),
         sprintf("PASS — %d rows, PP in (%.2e, %.4f], build=%s",
                 nrow(b), min(b$V5), max(b$V5),
                 readLines(paste0(path, ".build"))[1]))
}
check_bed(bed38_path)
check_bed(bed19_path)

for (f in list.files(SC_DIR, full.names = TRUE))
  push_to_bucket(f, file.path("results/scavenge", basename(f)))

log_msg("")
if (file.exists(bed19_path)) {
  log_msg("**STAGE B complete.** Both builds are on disk with `.build` sidecars; ",
          "Stage C picks the one matching the reference and asserts the match.")
} else {
  log_msg("**STAGE B complete — GRCh38 only.** The hg19 BED was not produced, so ",
          "Stage C will halt if the reference turns out to be hg19 (it is expected to be). ",
          "Stage the chain file and re-run this script before Stage C.")
}
