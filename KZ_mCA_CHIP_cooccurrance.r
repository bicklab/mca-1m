library(broom)
library(scales)
library(stringr)
library(survival)
library(tidyr)
library(dplyr)
library(haven)
library(reshape2)
library(tidyverse)
library(data.table)
library(lubridate)
library(ggplot2) 
library(openxlsx)

data_4cohorts <- readRDS("mca_all_BUAT_KZ_04012026.rds")

dim(data_4cohorts)
head(data_4cohorts)

data_4cohorts$mca <- 1

data_chip_BUAT <- readRDS("CHIP_BUAT_KZ_0714.rds")

data_chip_BUAT

kmt2d <- filter(data_chip_BUAT, data_chip_BUAT$Gene.refGene == "KMT2D")

table(data_chip_BUAT$cohort)

chip_biovu <- fread("chip_vars_DP_15p_AD2_2p_bs.tsv")
head(chip_biovu)
dim(chip_biovu)

ls(chip_biovu)

chip_biovu[, minAD := apply(.SD, 1, function(x) {
  x <- x[x > 0]
  if(length(x) == 0) NA else min(x)
}), .SDcols = c("AD_1","AD_2","AD_3","AD_4")]

chip_biovu$cohort <- "BioVU"
chip_biovu$chip <- 1

chip_biovu <- chip_biovu[,c("Gene.refGene","Chr","Start","End","Ref","Alt","Func.refGene","ExonicFunc.refGene",
                     "transcriptOI","NonsynOI","minAD", "DP", "AF","GRID","cohort", "chip")]

#add ancestry
pcs <- fread("all_agd_covariates.txt")
pcs <- pcs[,c(1,16)]

head(pcs)

pcs$ancestry <- dplyr::recode(
  pcs$unsupervised_ancestry_cluster_relabel,
  "Admixed_(majority_ancestry_<_0.5)" = "MID",
  "k1_(EAS)" = "EAS",
  "k2_(SAS)" = "SAS",
  "k3_(AFR)" = "AFR",
  "k4_(AMR)" = "AMR",
  "k5_(EUR)" = "EUR"
)
table(pcs$ancestry)

pcs <- pcs[,c(1,3)]

chip_biovu <- merge(chip_biovu, pcs, by = "GRID", all.x = T)

names(chip_biovu)[1] <- "person_id"

data_chip_UAT <- filter(data_chip_BUAT, data_chip_BUAT$cohort != "BioVU")
dim(data_chip_UAT)
table(data_chip_UAT$ancestry)

ls(chip_biovu)
ls(data_chip_UAT)

data_chip_BUAT <- rbind(data_chip_UAT,chip_biovu)
dim(data_chip_BUAT)

saveRDS(data_chip_BUAT, file = "data_chip_BUAT_20260406.rds")

data_chip_BUAT <- readRDS("data_chip_BUAT_20260406.rds")

gene_chip <- c("DNMT3A","TET2","ASXL1","TP53","JAK2","MPL","CBL",
               "SRSF2","SF3B1")

chip_long <- data_chip_BUAT %>%
  filter(Gene.refGene %in% gene_chip) %>%
  select(person_id, cohort, Gene.refGene) %>%
  mutate(CHIP = 1)

names(data_4cohorts)[1] <- "person_id"

mca_long <- data_4cohorts %>%
  filter(!is.na(chr_type2), type != "Undetermined", chrom != "chrX") %>%
  mutate(mCA_type = paste0(chr_type2, ":", type)) %>%
  select(person_id, cohort, mCA_type) %>%
  mutate(mCA = 1)

pheno_df <- data_4cohorts %>%
  select(person_id, age, sex = computed_gender, cohort) %>%
  distinct()

all_combos <- expand.grid(gene = gene_chip, mca = unique(mca_long$mCA_type), stringsAsFactors = FALSE)

head(chip_long)

dim(pheno_df)
pheno_df <- pheno_df %>% distinct(person_id, .keep_all = TRUE)
dim(pheno_df)

gene_mca_pairs <- chip_long %>%
  inner_join(mca_long, by = "person_id") %>%
  count(Gene.refGene, mCA_type) %>%
  filter(n >= 5)

gene_mca_pairs

library(logistf)

logistf_results <- list()

for (i in seq_len(nrow(gene_mca_pairs))) {
  gene <- gene_mca_pairs$Gene.refGene[i]
  mca  <- gene_mca_pairs$mCA_type[i]

  chip_df <- chip_long %>%
    filter(Gene.refGene == gene) %>%
    select(person_id, CHIP = CHIP)

  mca_df <- mca_long %>%
    filter(mCA_type == mca) %>%
    select(person_id, mCA = mCA)

  df <- pheno_df %>%
    left_join(chip_df, by = "person_id") %>%
    left_join(mca_df, by = "person_id") %>%
    mutate(
      CHIP = replace_na(CHIP, 0),
      mCA  = replace_na(mCA, 0)
    )

  # 共同拥有该 gene 和 mCA 的个体数
  N <- df %>% filter(CHIP == 1 & mCA == 1) %>% nrow()

  # Firth logistic 回归
  fit <- logistf::logistf(CHIP ~ mCA + age + sex + cohort, data = df)

  # OR 和置信区间
  OR <- exp(fit$coefficients["mCA"])
  ci <- exp(confint(fit)["mCA", ])

  # 原始 p 值和格式化 p 值（防止显示为0）
  p_val_raw <- fit$prob["mCA"]
  p_val_fmt <- format.pval(p_val_raw, digits = 10, eps = 1e-300)

  # 保存结果
  logistf_results[[paste(gene, mca, sep = "_")]] <- data.frame(
    gene = gene,
    mCA_type = mca,
    N = N,
    OR = OR,
    CI_lower = ci[1],
    CI_upper = ci[2],
    p_value = p_val_raw,
    p_value_formatted = p_val_fmt
  )
}

logistf_summary_df <- bind_rows(logistf_results)

logistf_df <- bind_rows(logistf_results)
logistf_df$logp <- -log10(logistf_df$p_value)
logistf_df$adjusted_p <- p.adjust(logistf_df$p_value, method = "bonferroni")
logistf_df

final_plot_df <- logistf_df %>%
  filter(adjusted_p < 0.05) %>%
  separate(mCA_type, into = c("chr_full", "alteration"), sep = ":", remove = FALSE) %>%
  mutate(
    chr = gsub("(chr\\d+)[pq]$", "\\1", chr_full),     
    arm = gsub("chr\\d+", "", chr_full),               
    alteration = factor(alteration, levels = c("Gain", "Loss", "CN-LOH"))
  )

final_plot_df <- filter(final_plot_df, final_plot_df$chr != "chrX")

final_plot_df$pattern <- NA
final_plot_df$pattern <- ifelse(final_plot_df$arm == "q", "dot", final_plot_df$pattern)
final_plot_df$pattern <- ifelse(final_plot_df$arm == "p", "stripe", final_plot_df$pattern)

library(openxlsx)
write.xlsx(final_plot_df, file = "mca_chip_coocurrence_BUAT_0403.xlsx")

final_plot_df <- read.xlsx("mca_chip_coocurrence_BUAT_0403.xlsx")

final_plot_df2 <- filter(final_plot_df, final_plot_df$N >15)

library(ggpattern)

final_plot_df2

gene_chip <- c("DNMT3A", "TET2", "ASXL1", "TP53", "JAK2", "MPL", "CBL", "SRSF2", "SF3B1")
alteration_colors <- c("Gain" = "#56B4E9", "Loss" = "#E69F00", "CN-LOH" = "#009E73")

final_plot_df2$gene <- factor(final_plot_df2$gene, levels = rev(gene_chip))

p <- ggplot(final_plot_df2, aes(x = chr, y = gene, fill = alteration)) +
  geom_tile(aes(alpha = log(OR)), color = "white", size = 0.5) +
  scale_fill_manual(values = alteration_colors) +
  scale_alpha_continuous(range = c(0.3, 1)) +
  geom_tile_pattern(aes(pattern = pattern),
                    fill = NA, color = "black", size = 0.2,
                    pattern_fill = NA, pattern_density = 0.1,
                    pattern_angle = 45) +
  scale_x_discrete(
    limits = paste0("chr", 1:22),
    drop = FALSE
  ) +
  theme_minimal(base_size = 13) +
  labs(
    title = "Adjusted Co-occurrence of CHIP and mCA",
    subtitle = "Firth logistic regression (adjusted for age, sex, cohort)",
    x = "Chromosome", y = "CHIP Gene", fill = "mCA Type", alpha = "log(OR)"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p

library(sjPlot)
save_plot("Fig1C_cooccurrence_BUAT_0407.svg", fig = p, width=20, height=14)

library(dplyr)
library(tidyr)

logistf_results <- list()

cohorts <- unique(pheno_df$cohort)

for (coh in cohorts) {
  
  cat("Processing cohort:", coh, "\n")
  
  # 👉 当前 cohort 数据
  chip_sub <- chip_long %>% filter(cohort == coh)
  mca_sub  <- mca_long  %>% filter(cohort == coh)
  pheno_sub <- pheno_df %>% filter(cohort == coh)
  
  # 👉 在当前 cohort 内筛选 gene-mCA pairs
  gene_mca_pairs <- chip_sub %>%
    inner_join(mca_sub, by = "person_id") %>%
    count(Gene.refGene, mCA_type) %>%
    filter(n >= 2)
  
  for (i in seq_len(nrow(gene_mca_pairs))) {
    
    gene <- gene_mca_pairs$Gene.refGene[i]
    mca  <- gene_mca_pairs$mCA_type[i]
    
    chip_df <- chip_sub %>%
      filter(Gene.refGene == gene) %>%
      select(person_id, CHIP)
    
    mca_df <- mca_sub %>%
      filter(mCA_type == mca) %>%
      select(person_id, mCA)
    
    df <- pheno_sub %>%
      left_join(chip_df, by = "person_id") %>%
      left_join(mca_df, by = "person_id") %>%
      mutate(
        CHIP = replace_na(CHIP, 0),
        mCA  = replace_na(mCA, 0)
      )
    
    # 👉 N (co-occurrence)
    N <- df %>% filter(CHIP == 1 & mCA == 1) %>% nrow()
    
    if (N < 5) next
    
    # 👉 Firth logistic（⚠️ 不再需要 + cohort）
    fit <- tryCatch(
      logistf::logistf(CHIP ~ mCA + age + sex, data = df),
      error = function(e) NULL
    )
    
    if (is.null(fit)) next
    
    # 👉 提取结果
    OR <- exp(fit$coefficients["mCA"])
    ci <- exp(confint(fit)["mCA", ])
    
    p_val_raw <- fit$prob["mCA"]
    p_val_fmt <- format.pval(p_val_raw, digits = 10, eps = 1e-300)
    
    logistf_results[[paste(coh, gene, mca, sep = "_")]] <- data.frame(
      cohort = coh,
      gene = gene,
      mCA_type = mca,
      N = N,
      OR = OR,
      CI_lower = ci[1],
      CI_upper = ci[2],
      p_value = p_val_raw,
      p_value_formatted = p_val_fmt
    )
  }
}

logistf_summary_df <- bind_rows(logistf_results)

logistf_summary_df

table(logistf_summary_df$cohort)

logistf_summary_df

matched_df <- logistf_summary_df %>%
  semi_join(final_plot_df2, by = c("gene", "mCA_type"))

matched_df

library(dplyr)

single_cohort_pairs <- matched_df %>%
  group_by(gene, mCA_type) %>%
  filter(n_distinct(cohort) == 1) %>%
  ungroup()

single_cohort_pairs

library(dplyr)

matched_df_filtered <- matched_df %>%
  anti_join(
    single_cohort_pairs %>% distinct(gene, mCA_type),
    by = c("gene", "mCA_type")
  )

write.csv(matched_df_filtered, file = "mCA_CHIP_coocur_bychort.csv")

sig_mca <- logistf_summary_df %>%
  group_by(mCA_type) %>%
  filter(
    n_distinct(cohort) == 2 &
    all(p_value < 0.1)
  ) %>%
  ungroup()

table(logistf_summary_df$cohort)

sig_mca

summary(data_4cohorts$cf)

names(data_4cohorts)[1] <- "person_id"

data_chip_BUAT_5 <- filter(data_chip_BUAT, data_chip_BUAT$AF >= 0.05)

data_mca_BUAT_5 <- filter(data_4cohorts, data_4cohorts$cf >= 0.05)

chip_long <- data_chip_BUAT_5 %>%
  filter(Gene.refGene %in% gene_chip) %>%
  select(person_id, cohort, Gene.refGene) %>%
  mutate(CHIP = 1)

mca_long <- data_mca_BUAT_5 %>%
  filter(!is.na(chr_type2), type != "Undetermined", chrom != "chrX") %>%
  mutate(mCA_type = paste0(chr_type2, ":", type)) %>%
  select(person_id, mCA_type) %>%
  mutate(mCA = 1)

pheno_df <- data_4cohorts %>%
  select(person_id, age, sex = computed_gender, cohort) %>%
  distinct()

all_combos <- expand.grid(gene = gene_chip, mca = unique(mca_long$mCA_type), stringsAsFactors = FALSE)

dim(pheno_df)
pheno_df <- pheno_df %>% distinct(person_id, .keep_all = TRUE)
dim(pheno_df)

gene_mca_pairs <- chip_long %>%
  inner_join(mca_long, by = "person_id") %>%
  count(Gene.refGene, mCA_type) %>%
  filter(n >= 5)

logistf_results <- list()

for (i in seq_len(nrow(gene_mca_pairs))) {
  gene <- gene_mca_pairs$Gene.refGene[i]
  mca  <- gene_mca_pairs$mCA_type[i]

  chip_df <- chip_long %>%
    filter(Gene.refGene == gene) %>%
    select(person_id, CHIP = CHIP)

  mca_df <- mca_long %>%
    filter(mCA_type == mca) %>%
    select(person_id, mCA = mCA)

  df <- pheno_df %>%
    left_join(chip_df, by = "person_id") %>%
    left_join(mca_df, by = "person_id") %>%
    mutate(
      CHIP = replace_na(CHIP, 0),
      mCA  = replace_na(mCA, 0)
    )

  N <- df %>% filter(CHIP == 1 & mCA == 1) %>% nrow()

  fit <- logistf::logistf(CHIP ~ mCA + age + sex + cohort, data = df)

  OR <- exp(fit$coefficients["mCA"])
  ci <- exp(confint(fit)["mCA", ])

  p_val_raw <- fit$prob["mCA"]
  p_val_fmt <- format.pval(p_val_raw, digits = 10, eps = 1e-300)

  logistf_results[[paste(gene, mca, sep = "_")]] <- data.frame(
    gene = gene,
    mCA_type = mca,
    N = N,
    OR = OR,
    CI_lower = ci[1],
    CI_upper = ci[2],
    p_value = p_val_raw,
    p_value_formatted = p_val_fmt
  )
}

logistf_summary_df <- bind_rows(logistf_results)

logistf_df <- bind_rows(logistf_results)
logistf_df$logp <- -log10(logistf_df$p_value)
logistf_df$adjusted_p <- p.adjust(logistf_df$p_value, method = "bonferroni")
logistf_df

final_plot_df <- logistf_df %>%
  filter(adjusted_p < 0.05) %>%
  separate(mCA_type, into = c("chr_full", "alteration"), sep = ":", remove = FALSE) %>%
  mutate(
    chr = gsub("(chr\\d+)[pq]$", "\\1", chr_full),     
    arm = gsub("chr\\d+", "", chr_full),               
    alteration = factor(alteration, levels = c("Gain", "Loss", "CN-LOH"))
  )

final_plot_df <- filter(final_plot_df, final_plot_df$chr != "chrX")

final_plot_df$pattern <- NA
final_plot_df$pattern <- ifelse(final_plot_df$arm == "q", "dot", final_plot_df$pattern)
final_plot_df$pattern <- ifelse(final_plot_df$arm == "p", "stripe", final_plot_df$pattern)

final_plot_df2 <- filter(final_plot_df, final_plot_df$N >15)
final_plot_df2

gene_chip <- c("DNMT3A", "TET2", "ASXL1", "TP53", "JAK2", "MPL", "CBL", "SRSF2", "SF3B1")

final_plot_df2$gene <- factor(final_plot_df2$gene, levels = rev(gene_chip))

p <- ggplot(final_plot_df2, aes(x = chr, y = gene, fill = alteration)) +
  geom_tile(aes(alpha = log(OR)), color = "white", size = 0.5) +
  scale_fill_manual(values = alteration_colors) +
  scale_alpha_continuous(range = c(0.3, 1)) +
  geom_tile_pattern(aes(pattern = pattern),
                    fill = NA, color = "black", size = 0.2,
                    pattern_fill = NA, pattern_density = 0.1,
                    pattern_angle = 45) +
  scale_x_discrete(
    limits = paste0("chr", 1:22),
    drop = FALSE
  ) +
  theme_minimal(base_size = 13) +
  labs(
    title = "Adjusted Co-occurrence of CHIP and mCA",
    subtitle = "Firth logistic regression (adjusted for age, sex, cohort)",
    x = "Chromosome", y = "CHIP Gene", fill = "mCA Type", alpha = "log(OR)"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p

save_plot("Fig1C_cooccurrence_BUAT_AF5.svg", fig = p, width=24, height=12)


