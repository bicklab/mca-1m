#!/usr/bin/env Rscript
# Minimum shared altered region (MSAR) definition.
# Reviewer #2 (code availability) noted this was missing from the repository.
#
# Input : pooled aut-mCA call set, one row per event, with columns
#         sample_id, cohort, chrom, start, end, arm ("p"/"q"/"whole"), type ("+","-","=")
# Output: msar_regions.tsv  -> one row per (mCA label x threshold)
#         Supplementary Table 5 is the subset at the selected threshold per mCA.
#
# Definition: for a given mCA label (chrom+arm+type) and threshold f, the MSAR is the
# largest interval covered by at least f * N of the events with that label. It is
# computed from a coverage step function over event breakpoints, so it does not
# assume the events are nested.

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
calls_file <- if (length(args) >= 1) args[1] else "data/mca_calls_pooled.tsv"
out_file   <- if (length(args) >= 2) args[2] else "results/msar_regions.tsv"

THRESHOLDS <- c(0.50, 0.75, 0.80, 0.90)
MIN_EVENTS <- 25L   # same minimum used for the cis-GWAS phenotype inclusion

calls <- fread(calls_file)
stopifnot(all(c("chrom", "start", "end", "arm", "type") %in% names(calls)))

calls[, mca_label := paste0(chrom, ifelse(arm == "whole", "", arm), type)]

# Coverage step function over breakpoints; returns the widest run of positions
# covered by >= n_needed events.
msar_one <- function(start, end, n_needed) {
  bp <- sort(unique(c(start, end)))
  if (length(bp) < 2L) return(c(NA_real_, NA_real_))
  # coverage of each inter-breakpoint segment
  mid <- (head(bp, -1L) + tail(bp, -1L)) / 2
  cov <- vapply(mid, function(m) sum(start <= m & end >= m), integer(1))
  keep <- cov >= n_needed
  if (!any(keep)) return(c(NA_real_, NA_real_))
  # longest contiguous run of qualifying segments
  r <- rle(keep)
  ends <- cumsum(r$lengths)
  starts <- ends - r$lengths + 1L
  ok <- which(r$values)
  widths <- (bp[ends[ok] + 1L] - bp[starts[ok]])
  best <- ok[which.max(widths)]
  c(bp[starts[best]], bp[ends[best] + 1L])
}

res <- rbindlist(lapply(THRESHOLDS, function(f) {
  calls[, {
    n <- .N
    if (n < MIN_EVENTS) {
      NULL
    } else {
      iv <- msar_one(start, end, ceiling(f * n))
      list(threshold = f, n_events = n,
           start_hg38 = iv[1], end_hg38 = iv[2],
           length_mb = (iv[2] - iv[1]) / 1e6)
    }
  }, by = .(mca_label, chrom)]
}))

fwrite(res[order(chrom, mca_label, -threshold)], out_file, sep = "\t")
cat("wrote", out_file, "-", nrow(res), "rows\n")
