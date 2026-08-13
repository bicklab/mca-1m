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
# # BUILD REFERENCE CACHE — the trait-independent half of SCAVENGE, computed once
#
# Most of Stage C's cost has nothing to do with the trait. GC bias, background
# peaks, TF-IDF, the LSI SVD over 571,400 × 63,882, and the mutual-kNN graph are
# functions of the **reference alone**. Recomputing them for every CHIP driver ×
# every peak-extension setting would dominate the wall clock — roughly 10 drivers
# × 2 extensions × the full pipeline.
#
# Computing them once buys two things:
#   * time — each trait then only needs importBedScore → deviations → seeds →
#     random walk → permutations;
#   * comparability — every driver is scored on the *identical* cell graph, so
#     TRS rankings can be compared across traits without wondering whether the
#     embedding shifted.
#
# The cache records the parameters it was built with and is rebuilt automatically
# if any of them change (`00_setup/config_load.R` seeds the RNG, so the LSI and
# kNN are deterministic given those parameters).
#
# **Run once per app**, before any trait: `Rscript 03_scavenge/build_ref_cache.R`.
# ~30–60 min and the memory high-water mark of the whole pipeline (≥64 GB).

# %%
Sys.setenv(FS_NO_GWAS = "1")          # trait-independent: no results/<name>/ dir
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
  library(Matrix); library(SummarizedExperiment); library(GenomicRanges)
  library(chromVAR); library(SCAVENGE)
})
SCV <- CFG$scavenge
CACHE <- file.path(CACHE_DIR, "scavenge_ref_cache.rds")

log_header("BUILD REFERENCE CACHE")
log_versions(c("SCAVENGE", "chromVAR", "SummarizedExperiment", "Matrix", "irlba"))

# Anything here changing invalidates the cache — these are exactly the inputs the
# cached objects depend on.
key <- list(reference = basename(CFG$paths$reference_local),
            build     = SCV$reference_build,
            n_bg      = SCV$n_background_peaks,
            lsi_dims  = SCV$lsi_dims,
            knn_k     = SCV$knn_k,
            seed      = CFG$project$seed)
if (file.exists(CACHE)) {
  old <- readRDS(CACHE)
  if (identical(old$key, key)) {
    log_kv("cache", "up to date — nothing to do")
    cat("\nCache is current:\n"); str(key)
    quit(status = 0, save = "no")
  }
  log_kv("cache", "parameters changed — rebuilding")
  cat("\nwas:\n"); str(old$key); cat("now:\n"); str(key)
}

# %%
ref_path <- proj_path(CFG$paths$reference_local)
assert_that(file.exists(ref_path),
            paste0("Reference not found at ", ref_path,
                   "\n  Run: 00_setup/stage_inputs.sh fetch"))
SE <- readRDS(ref_path)
if (!"counts" %in% assayNames(SE)) assayNames(SE)[1] <- "counts"
log_kv("reference", paste0(basename(ref_path), " — ",
                           format(nrow(SE), big.mark = ","), " peaks x ",
                           format(ncol(SE), big.mark = ","), " cells"))

bd <- detect_genome_build(rowRanges(SE), expected = SCV$reference_build,
                          allow_unverified = FALSE)
log_kv("reference build", paste0(bd$build, " (via ", bd$how, ")"))

BSG <- if (bd$build == "hg19") "BSgenome.Hsapiens.UCSC.hg19" else "BSgenome.Hsapiens.UCSC.hg38"
assert_that(requireNamespace(BSG, quietly = TRUE),
            paste0(BSG, " is not installed — run 00_setup/install_packages.R."))

n0 <- nrow(SE)
SE <- chromVAR::filterPeaks(SE, non_overlapping = TRUE)
log_kv("peaks after filterPeaks", paste0(format(nrow(SE), big.mark = ","),
                                         " (dropped ", n0 - nrow(SE), ")"))

# %%
log_kv("step", "addGCBias")
SE <- chromVAR::addGCBias(SE, genome = getExportedValue(BSG, BSG))
log_kv("step", paste0("getBackgroundPeaks (niterations = ", SCV$n_background_peaks, ")"))
bg <- chromVAR::getBackgroundPeaks(SE, niterations = SCV$n_background_peaks)

log_kv("step", "tfidf -> do_lsi (the expensive one)")
tf  <- SCAVENGE::tfidf(bmat = assays(SE)[["counts"]], mat_binary = TRUE,
                       TF = TRUE, log_TF = TRUE)
lsi <- SCAVENGE::do_lsi(mat = tf, dims = SCV$lsi_dims)
rm(tf); invisible(gc())

log_kv("step", paste0("getmutualknn (k = ", SCV$knn_k, ")"))
knn <- SCAVENGE::getmutualknn(lsimat = lsi, num_k = SCV$knn_k)
log_kv("mutual-kNN graph", sprintf("%d x %d, %s edges", nrow(knn), ncol(knn),
                                   format(length(knn@x), big.mark = ",")))

# %%
# Store only what a trait run cannot recompute cheaply. The counts matrix and
# rowRanges are NOT cached — Stage C reloads the SE anyway for computeWeightedDeviations,
# and duplicating a multi-GB matrix on disk would buy nothing.
saveRDS(list(key = key,
             bias = rowData(SE)$bias,          # GC bias, re-attached by Stage C
             kept_peaks = rownames(SE) %||% seq_len(nrow(SE)),
             n_peaks = nrow(SE), n_cells = ncol(SE),
             bg = bg, lsi = lsi, knn = knn,
             build = bd$build, built_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
        CACHE, compress = FALSE)     # speed over size; this is scratch, not an artefact

log_kv("cache written", paste0(CACHE, " (",
       format(structure(file.size(CACHE), class = "object_size"), units = "auto"), ")"))
log_msg("")
log_msg("**Reference cache built.** Every trait now reuses this graph, so TRS ",
        "rankings are comparable across CHIP drivers. Next: `./run_pipeline.sh`.")
