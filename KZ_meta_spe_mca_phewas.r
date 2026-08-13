library(survival)
library(tidyr)
library(dplyr)
library(haven)
library(reshape2)
library(tidyverse)
library(data.table)
library(lubridate)
library(ggplot2)
library(broom)
library(scales)
library(stringr)
library(metafor)

library(data.table)

folder <- "phewas_biovu250k_mca"

files <- list.files(folder, pattern = "_biovu250k\\.csv$", full.names = TRUE)

for (f in files) {
  dt <- fread(f)
  
  setnames(dt, old = names(dt)[1], new = "disease")
  
  fwrite(dt, f)
}

folder_path <- "phewas_biovu250k_mca"

file_list <- list.files(folder_path, pattern = "\\.csv$", full.names = TRUE)

data_list <- list()

for (file in file_list) {
  file_name <- gsub("\\.csv$", "", basename(file))
  data_list[[file_name]] <- read.csv(file)
}

names(data_list)

chromosome_types <- unique(gsub("(_AoU|_biovu250k|_ukb)$", "", names(data_list)))
chromosome_types

chr14q.CN.LOH_BioVU <- read.csv("phewas_biovu250k_mca/chr14q.CN.LOH_biovu250k.csv")
chr14q.CN.LOH_BioVU

chr14q.CN.LOH_aou <- read.csv("phewas_biovu250k_mca/chr14q.CN.LOH_AoU.csv")
chr14q.CN.LOH_aou

chr14q.CN.LOH_ukb <- read.csv("phewas_biovu250k_mca/chr14q.CN.LOH_ukb.csv")
chr14q.CN.LOH_ukb

combined_list <- list()

for (chrom_type in chromosome_types) {
  aou_data <- data_list[[paste0(chrom_type, "_AoU")]]
  biovu_data <- data_list[[paste0(chrom_type, "_biovu250k")]]
  ukb_data <- data_list[[paste0(chrom_type, "_ukb")]]
  
  colnames(ukb_data) <- paste0(colnames(ukb_data), "_ukb")
  combined_data <- merge(aou_data, biovu_data, by = "disease", all = TRUE, suffixes = c("_aou", "_biovu"))
  combined_data <- merge(combined_data, ukb_data, by.x = "disease",by.y = "disease_ukb", all = TRUE)
  
  combined_list[[paste0(chrom_type, "_combined")]] <- combined_data
}

names(combined_list)

ls(combined_list[["chr11p.CN.LOH_combined"]])

for (name in names(combined_list)) {
  cat("Processing:", name, "\n")
  
  combined_data <- combined_list[[name]]
  
  if ("total_samples_aou" %in% colnames(combined_data)) {
    combined_data$total_samples_aou[is.na(combined_data$total_samples_aou)] <- 0
  }
  
  if ("total_samples_biovu" %in% colnames(combined_data)) {
    combined_data$total_samples_biovu[is.na(combined_data$total_samples_biovu)] <- 0
  }
  
  if ("total_samples_ukb" %in% colnames(combined_data)) {
    combined_data$total_samples_ukb[is.na(combined_data$total_samples_ukb)] <- 0
  }
  
  combined_list[[name]] <- combined_data
}

str(combined_list)

names(combined_list)

combined_list[["chr11p.CN.LOH_combined"]]

combined_data <- combined_list[["chr20q.Loss_combined"]]

nrow(combined_data)

library(data.table)
library(metafor)

all_meta_results <- list()
meta_qc_summary <- list()

for (combined_name in names(combined_list)) {
  
  cat("Processing:", combined_name, "\n")
  flush.console()
  
  combined_data <- as.data.table(combined_list[[combined_name]])
  meta_results <- vector("list", nrow(combined_data))
  
  qc_counter <- list(
    total_rows = nrow(combined_data),
    skipped_total_n_lt50 = 0L,
    skipped_no_valid_cohort = 0L,
    skipped_meta_failed = 0L,
    kept_single = 0L,
    kept_meta = 0L
  )
  
  pb <- txtProgressBar(min = 0, max = nrow(combined_data), style = 3)
  
  for (i in seq_len(nrow(combined_data))) {
    setTxtProgressBar(pb, i)
    
    disease_i <- combined_data$disease[i]
    
    # 1) 提取三队列原始结果
    HRs <- c(
      combined_data$HR_aou[i],
      combined_data$HR_biovu[i],
      combined_data$HR_ukb[i]
    )
    
    CI_lower <- c(
      combined_data$CI_lower_aou[i],
      combined_data$CI_lower_biovu[i],
      combined_data$CI_lower_ukb[i]
    )
    
    CI_upper <- c(
      combined_data$CI_upper_aou[i],
      combined_data$CI_upper_biovu[i],
      combined_data$CI_upper_ukb[i]
    )
    
    ns <- c(
      combined_data$N_aou[i],
      combined_data$N_biovu[i],
      combined_data$N_ukb[i]
    )
    
    n_events <- c(
      combined_data$N_event_aou[i],
      combined_data$N_event_biovu[i],
      combined_data$N_event_ukb[i]
    )
    
    ps <- c(
      combined_data$p_value_aou[i],
      combined_data$p_value_biovu[i],
      combined_data$p_value_ukb[i]
    )
    
    logps <- c(
      combined_data$log_p_value_aou[i],
      combined_data$log_p_value_biovu[i],
      combined_data$log_p_value_ukb[i]
    )
    
    cohorts <- c("aou", "biovu", "ukb")
    
    # 强制 numeric
    HRs <- as.numeric(HRs)
    CI_lower <- as.numeric(CI_lower)
    CI_upper <- as.numeric(CI_upper)
    ns <- as.numeric(ns)
    n_events <- as.numeric(n_events)
    ps <- as.numeric(ps)
    logps <- as.numeric(logps)
    
    total_samples_combined <- sum(ns, na.rm = TRUE)
    if (is.na(total_samples_combined) || total_samples_combined < 50) {
      qc_counter$skipped_total_n_lt50 <- qc_counter$skipped_total_n_lt50 + 1L
      next
    }
    
    # 2) 先算 SE(log(HR))
    ses <- (log(CI_upper) - log(CI_lower)) / (2 * 1.96)
    betas <- log(HRs)
    
    # 3) robust QC 过滤
    valid_idx <- which(
      !is.na(HRs) &
      !is.na(CI_lower) &
      !is.na(CI_upper) &
      !is.na(ns) &
      !is.na(n_events) &
      is.finite(HRs) &
      is.finite(CI_lower) &
      is.finite(CI_upper) &
      is.finite(ses) &
      is.finite(betas) &
      HRs > 0.01 & HRs < 50 &
      CI_lower > 0.01 &
      CI_upper > 0.01 &
      CI_lower < CI_upper &
      ses > 0 &
      n_events >= 5
    )
    
    if (length(valid_idx) == 0) {
      qc_counter$skipped_no_valid_cohort <- qc_counter$skipped_no_valid_cohort + 1L
      next
    }
    
    HRs_valid <- HRs[valid_idx]
    betas_valid <- betas[valid_idx]
    ses_valid <- pmax(ses[valid_idx], 1e-4)
    ps_valid <- ps[valid_idx]
    ns_valid <- ns[valid_idx]
    logps_valid <- logps[valid_idx]
    cohorts_valid <- cohorts[valid_idx]
    
    # 4) single cohort 或 FE meta
    if (length(valid_idx) == 1) {
      meta_beta <- betas_valid
      meta_se <- ses_valid
      meta_p <- ps_valid
      meta_model <- "single"
      qc_counter$kept_single <- qc_counter$kept_single + 1L
      
    } else {
      meta_res <- tryCatch(
        rma(yi = betas_valid, sei = ses_valid, method = "FE"),
        error = function(e) NULL
      )
      
      if (is.null(meta_res)) {
        qc_counter$skipped_meta_failed <- qc_counter$skipped_meta_failed + 1L
        next
      }
      
      meta_beta <- as.numeric(meta_res$beta)
      meta_se <- as.numeric(meta_res$se)
      meta_p <- as.numeric(meta_res$pval)
      meta_model <- meta_res$method
      qc_counter$kept_meta <- qc_counter$kept_meta + 1L
    }
    
    # 5) 转回 HR 和 95%CI
    meta_HR <- exp(meta_beta)
    meta_CI_lower <- exp(meta_beta - 1.96 * meta_se)
    meta_CI_upper <- exp(meta_beta + 1.96 * meta_se)
    meta_logp <- ifelse(!is.na(meta_p) && meta_p > 0, -log10(meta_p), NA_real_)
    
    # 6) 保存
    meta_results[[i]] <- data.frame(
      phenotype = disease_i,
      
      HR_aou = combined_data$HR_aou[i],
      CI_lower_aou = combined_data$CI_lower_aou[i],
      CI_upper_aou = combined_data$CI_upper_aou[i],
      p_value_aou = combined_data$p_value_aou[i],
      log_p_value_aou = combined_data$log_p_value_aou[i],
      total_samples_aou = combined_data$N_aou[i],
      N_event_aou = combined_data$N_event_aou[i],
      
      HR_biovu = combined_data$HR_biovu[i],
      CI_lower_biovu = combined_data$CI_lower_biovu[i],
      CI_upper_biovu = combined_data$CI_upper_biovu[i],
      p_value_biovu = combined_data$p_value_biovu[i],
      log_p_value_biovu = combined_data$log_p_value_biovu[i],
      total_samples_biovu = combined_data$N_biovu[i],
      N_event_biovu = combined_data$N_event_biovu[i],
      
      HR_ukb = combined_data$HR_ukb[i],
      CI_lower_ukb = combined_data$CI_lower_ukb[i],
      CI_upper_ukb = combined_data$CI_upper_ukb[i],
      p_value_ukb = combined_data$p_value_ukb[i],
      log_p_value_ukb = combined_data$log_p_value_ukb[i],
      total_samples_ukb = combined_data$N_ukb[i],
      N_event_ukb = combined_data$N_event_ukb[i],
      
      total_samples_combined = total_samples_combined,
      valid_cohorts = paste(cohorts_valid, collapse = ";"),
      n_valid_cohorts = length(valid_idx),
      
      meta_beta = meta_beta,
      meta_se = meta_se,
      meta_p = meta_p,
      meta_logp = meta_logp,
      meta_HR = meta_HR,
      meta_CI_lower = meta_CI_lower,
      meta_CI_upper = meta_CI_upper,
      meta_model = meta_model
    )
  }
  
  close(pb)
  
  # 7) 安全拼接，避免空结果写出 NA 文件
  meta_results_df <- rbindlist(meta_results[!sapply(meta_results, is.null)], fill = TRUE)
  
  if (nrow(meta_results_df) == 0) {
    cat("No valid results for:", combined_name, "\n\n")
    flush.console()
    
    meta_qc_summary[[combined_name]] <- data.frame(
      combined_name = combined_name,
      total_rows = qc_counter$total_rows,
      skipped_total_n_lt50 = qc_counter$skipped_total_n_lt50,
      skipped_no_valid_cohort = qc_counter$skipped_no_valid_cohort,
      skipped_meta_failed = qc_counter$skipped_meta_failed,
      kept_single = qc_counter$kept_single,
      kept_meta = qc_counter$kept_meta,
      final_rows = 0
    )
    next
  }
  
  setorder(meta_results_df, -meta_logp)
  
  all_meta_results[[combined_name]] <- meta_results_df
  
  output_file <- paste0("meta_FE_", combined_name, ".csv")
  fwrite(meta_results_df, output_file)
  
  meta_qc_summary[[combined_name]] <- data.frame(
    combined_name = combined_name,
    total_rows = qc_counter$total_rows,
    skipped_total_n_lt50 = qc_counter$skipped_total_n_lt50,
    skipped_no_valid_cohort = qc_counter$skipped_no_valid_cohort,
    skipped_meta_failed = qc_counter$skipped_meta_failed,
    kept_single = qc_counter$kept_single,
    kept_meta = qc_counter$kept_meta,
    final_rows = nrow(meta_results_df)
  )
  
  cat("Finished:", combined_name, " | rows kept =", nrow(meta_results_df), "\n\n")
  flush.console()
}

qc_summary_df <- rbindlist(meta_qc_summary, fill = TRUE)
fwrite(qc_summary_df, "meta_qc_summary.csv")

cat("All meta analysis completed.\n")

test <- fread("meta_FE_chr20q.Loss_combined.csv")
test

library(data.table)
library(stringr)

meta_files <- list.files(pattern = "^meta_FE.*\\.csv$")

combined_meta <- rbindlist(lapply(meta_files, function(file) {
  
  df <- fread(file)
  mca_type <- file |>
    str_remove("^meta_FE_") |>
    str_remove("\\.csv$") |>
    str_remove("_combined$")
  
  df[, mca_type := mca_type]
  
  return(df)
}), fill = TRUE)

head(combined_meta)

log10(0.05/2000)

combined_meta_sig <- filter(combined_meta, combined_meta$meta_logp >= 4.6)
dim(combined_meta_sig)

combined_meta_sig

phecodeX <- read.csv("phecodeX_ICD_CM_map_flat.csv")

phecodeX <- phecodeX[,c(5,7)]
phecodeX <- distinct(phecodeX)
head(phecodeX)

combined_meta_sig2 <- merge(combined_meta_sig, phecodeX , by.x =  "phenotype", by.y = "phecode_string" , all.x = T)

table(combined_meta_sig2$category)

write.csv(combined_meta_sig2, file = "meta_spe_mCA_BUAT.csv")

combined_meta_sig_cv <- filter(combined_meta_sig2, combined_meta_sig2$category == "Cardiovascular")

combined_meta_sig_cv

combined_meta_sig_ge <- filter(combined_meta_sig2, combined_meta_sig2$category == "Genitourinary")

combined_meta_sig_ge

combined_meta_sig_inf <- filter(combined_meta_sig2, combined_meta_sig2$category == "Infections")
combined_meta_sig_inf


