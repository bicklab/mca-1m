library(broom)
library(scales)
library(stringr)
library(arrow)
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

mca_biovu250k <- fread("mocha_calls_AGD_250k_dedup_filtered.tsv")

data_all_pheno_cox7 <- readRDS("BioVU_PhenoData_250k_0922.rds") 
dim(data_all_pheno_cox7)

table(data_all_pheno_cox7$unsupervised_ancestry_cluster_relabel)

195000+110000

dim(data_all_pheno_cox7)

mapping_biovu <- data_all_pheno_cox7[,c("GRID", "person_id")]

mapping_biovu <- saveRDS(mapping_biovu, "mapping_biovu.rds")

head(mca_biovu250k)
dim(mca_biovu250k)

system("gsutil cp gs://bicklab-main-storage/Data/BioVU/ICA/sample_maps/AGD_March_2025_SampleIDMap.csv .")

sample_mapping <- read_csv('AGD_March_2025_SampleIDMap.csv')
 
mca_biovu250k <- mca_biovu250k %>% inner_join(sample_mapping, by=c('sample_id'='SAMPLE_ID')) %>%
    select(GRID, computed_gender:cf)

head(mca_biovu250k)

ls(mca_biovu250k)

biovu_info <- data_all_pheno_cox7[, .(GRID, age, unsupervised_ancestry_cluster_relabel, gender)]

data_mca_biovu250k <- merge(
  mca_biovu250k,
  biovu_info,
  by = "GRID", 
  all = F
)

dim(data_mca_biovu250k)

data_mca_auto <- filter(data_mca_biovu250k, data_mca_biovu250k$chrom != "chrX")
dim(data_mca_auto)

data_mca_auto %>% pull(GRID) %>% unique() %>% length()

mca_biovu_sex <- data_mca_biovu250k[,c("sample_name","gender")]

mca_biovu_sex <- distinct(mca_biovu_sex)
table(mca_biovu_sex$gender)

data_all_pheno_cox7$mca_auto <- NA
data_all_pheno_cox7$mca_auto <- ifelse(data_all_pheno_cox7$GRID %in% data_mca_auto$GRID, 1, 0)
table(data_all_pheno_cox7$mca_auto)

2144/(243207+2144)

#add CF
summary(data_mca_auto$cf)

cf10_id <- filter(data_mca_auto, data_mca_auto$cf > 0.1)

dim(cf10_id)
dim(data_mca_auto)

data_all_pheno_cox7$auto_mca_cf_l <- NA
data_all_pheno_cox7$auto_mca_cf_l <- ifelse(data_all_pheno_cox7$mca_auto ==0 , 0, data_all_pheno_cox7$auto_mca_cf_l)
data_all_pheno_cox7$auto_mca_cf_l <- ifelse(data_all_pheno_cox7$GRID %in% cf10_id$GRID, 1, data_all_pheno_cox7$auto_mca_cf_l)

data_all_pheno_cox7$auto_mca_cf_m <- NA
data_all_pheno_cox7$auto_mca_cf_m <- ifelse(data_all_pheno_cox7$mca_auto ==0 , 0, data_all_pheno_cox7$auto_mca_cf_m)
data_all_pheno_cox7$auto_mca_cf_m <- ifelse(is.na(data_all_pheno_cox7$auto_mca_cf_l), 1, data_all_pheno_cox7$auto_mca_cf_m)

table(data_all_pheno_cox7$auto_mca_cf_l, data_all_pheno_cox7$mca_auto, useNA = "always")

table(data_all_pheno_cox7$auto_mca_cf_l, data_all_pheno_cox7$auto_mca_cf_m, useNA = "always")

#High risk mCAs for CLL include loss of chr6/chr11/chr13/chr17, gain of chr12, or CNLOH of chr13
data_mca_auto$mac_high<-NA
data_mca_auto$mac_high<-ifelse(data_mca_auto$type=="Loss" & data_mca_auto$chrom =="chr6",1, data_mca_auto$mac_high)
data_mca_auto$mac_high<-ifelse(data_mca_auto$type=="Loss" & data_mca_auto$chrom =="chr11",1,data_mca_auto$mac_high)
data_mca_auto$mac_high<-ifelse(data_mca_auto$type=="Loss" & data_mca_auto$chrom =="chr13",1,data_mca_auto$mac_high)
data_mca_auto$mac_high<-ifelse(data_mca_auto$type=="Loss" & data_mca_auto$chrom =="chr17",1,data_mca_auto$mac_high)
data_mca_auto$mac_high<-ifelse(data_mca_auto$type=="Gain" & data_mca_auto$chrom =="chr12",1,data_mca_auto$mac_high)
data_mca_auto$mac_high<-ifelse(data_mca_auto$type=="CN-LOH" & data_mca_auto$chrom =="chr13",1,data_mca_auto$mac_high)
table(data_mca_auto$mac_high) #1305

id_mcahigh<-filter(data_mca_auto,data_mca_auto$mac_high==1)
id_mcahigh<-id_mcahigh[!duplicated(id_mcahigh$GRID),] #423
id_mcahigh<-id_mcahigh[,c("GRID","mac_high")]

data_all_pheno_cox7$auto_mca_highrisk <- NA
data_all_pheno_cox7$auto_mca_highrisk <- ifelse(data_all_pheno_cox7$mca_auto ==0 , 0, data_all_pheno_cox7$auto_mca_highrisk)
data_all_pheno_cox7$auto_mca_highrisk <- ifelse(data_all_pheno_cox7$GRID %in% id_mcahigh$GRID, 1, data_all_pheno_cox7$auto_mca_highrisk)

table(data_all_pheno_cox7$auto_mca_highrisk)

data_mca_mlox <- filter(data_mca_biovu250k, data_mca_biovu250k$chrom == "chrX" & data_mca_biovu250k$type == "Loss")
dim(data_mca_mlox)

data_all_pheno_cox7$mlox <- NA
data_all_pheno_cox7$mlox <- ifelse(data_all_pheno_cox7$GRID %in% data_mca_mlox$GRID, 1, 0)
table(data_all_pheno_cox7$mlox)

saveRDS(data_all_pheno_cox7,file = "data_all_pheno_mca_biovu.rds")

data_mca <- data_all_pheno_cox7[,c("GRID","mca_auto","auto_mca_highrisk","mlox")]

write.csv(data_mca, file = "id_mca_biovu250k.csv")

data_all_pheno_cox7 <- readRDS("data_all_pheno_mca_biovu.rds")

#age 40
data_all_pheno_cox7_40 <- filter(data_all_pheno_cox7, data_all_pheno_cox7$age >= 40)

table(data_all_pheno_cox7_40$mca_auto)

library(survival)
library(dplyr)
library(rlang)
library(broom)

run_cox_model_each_var <- function(data,
                                   event_col,
                                   event_time_col,
                                   baseline_time_col,
                                   last_record_time_col,
                                   covariates,
                                   variables_of_interest) {
  
  event_sym <- sym(event_col)
  event_time_sym <- sym(event_time_col)
  baseline_sym <- sym(baseline_time_col)
  last_sym <- sym(last_record_time_col)
  
  df <- data %>%
    mutate(
      followup_time = as.numeric(pmin(!!event_time_sym, !!last_sym, na.rm = TRUE) - !!baseline_sym),
      status = !!event_sym
    ) %>%
    filter(!is.na(followup_time), followup_time > 0)
  
  total_event <- sum(df$status == 1, na.rm = TRUE)

  surv_obj <- Surv(time = df$followup_time, event = df$status)

  results_list <- list()

  for (var in variables_of_interest) {
    var_sym <- sym(var)

    N <- df %>% filter(!!var_sym == 1) %>% nrow()
    N_event <- df %>% filter(!!var_sym == 1, status == 1) %>% nrow()
    
    event_rate <- if (total_event > 0) N_event / total_event else NA_real_

    formula_str <- paste("surv_obj ~", paste(c(covariates, var), collapse = " + "))
    cox_formula <- as.formula(formula_str)

    model <- coxph(cox_formula, data = df)

    result <- broom::tidy(model, exponentiate = TRUE, conf.int = TRUE) %>%
      filter(term == var) %>%
      mutate(
        variable = var,
        N = N,
        N_event = N_event,
        event_rate = event_rate,
        HR = estimate,
        CI_lower = conf.low,
        CI_upper = conf.high,
        p_value = p.value
      ) %>%
      select(variable, N, N_event, event_rate, HR, CI_lower, CI_upper, p_value)
    
    results_list[[var]] <- result
  }

  final_result <- bind_rows(results_list)

  return(final_result)
}

run_cox_model_2times <- function(data,
                                  event_col,
                                  event_time_col,
                                  baseline_time_col,
                                  last_record_time_col,
                                  age_blooddraw,
                                  covariates,
                                  variables_of_interest,
                                  return_df = FALSE,
                                  check_ph = FALSE) {   # 新增：check_ph 开关

  event_sym <- sym(event_col)
  event_time_sym <- sym(event_time_col)
  baseline_sym <- sym(baseline_time_col)
  last_sym <- sym(last_record_time_col)
  age_sym <- sym(age_blooddraw)

  df <- data %>%
    mutate(
      age_event = as.numeric((!!event_time_sym - !!baseline_sym) / 365.25) + !!age_sym,
      age_last  = as.numeric((!!last_sym       - !!baseline_sym) / 365.25) + !!age_sym,
      time1     = !!age_sym,
      time2     = pmin(age_event, age_last, na.rm = TRUE),
      status    = !!event_sym
    ) %>%
    filter(!is.na(time1), !is.na(time2), time2 > time1)

  if (return_df) return(df) 

  surv_obj <- Surv(time = df$time1, time2 = df$time2, event = df$status)

  total_event <- sum(df$status == 1, na.rm = TRUE)

  results_list <- list()

  for (var in variables_of_interest) {
    var_sym <- sym(var)

    N <- df %>% filter(!!var_sym == 1) %>% nrow()
    N_event <- df %>% filter(!!var_sym == 1, status == 1) %>% nrow()
    event_rate <- if (total_event > 0) N_event / total_event else NA_real_

    formula_str <- paste("surv_obj ~", paste(c(covariates, var), collapse = " + "))
    cox_formula <- as.formula(formula_str)

    model <- coxph(cox_formula, data = df)

    # 默认结果
    result <- broom::tidy(model, exponentiate = TRUE, conf.int = TRUE) %>%
      filter(term == var) %>%
      mutate(
        variable   = var,
        N          = N,
        N_event    = N_event,
        event_rate = event_rate,
        HR         = estimate,
        CI_lower   = conf.low,
        CI_upper   = conf.high,
        p_value    = p.value
      ) %>%
      select(variable, N, N_event, event_rate, HR, CI_lower, CI_upper, p_value)

    # 如果要检查 PH 假定，在这里加一列 PH 检验 p 值
    if (check_ph) {
      ph_test <- cox.zph(model)
      # ph_test$table 是一个矩阵，行名是协变量名和 "GLOBAL"
      # 取当前这个 var 对应行的 p 值（最后一列）
      if (var %in% rownames(ph_test$table)) {
        ph_p <- ph_test$table[var, "p"]
      } else {
        ph_p <- NA_real_
      }
      result <- result %>%
        mutate(PH_p = ph_p)
    }

    results_list[[var]] <- result
  }

  final_result <- bind_rows(results_list)
  return(final_result)
}

data_all_pheno_cox7_40$auto_mca_lowrisk <- NA
data_all_pheno_cox7_40$auto_mca_lowrisk <- ifelse(is.na(data_all_pheno_cox7_40$auto_mca_highrisk), 1, 0)
table(data_all_pheno_cox7_40$auto_mca_lowrisk)
data_all_pheno_cox7_40$auto_mca_lowrisk <- ifelse(data_all_pheno_cox7_40$auto_mca_highrisk == 1 & is.na(data_all_pheno_cox7_40$auto_mca_highrisk) == F, NA, data_all_pheno_cox7_40$auto_mca_lowrisk)

data_all_pheno_cox7_40$SHIFTED_SAMPLE_DATE <- as.Date(data_all_pheno_cox7_40$SHIFTED_SAMPLE_DATE)
data_all_pheno_cox7_40$observation_period_end_date <- as.Date(data_all_pheno_cox7_40$observation_period_end_date)

class(data_all_pheno_cox7_40$time_CA_121.21 )

data_all_pheno_cox7_40$time_CA_121.21 <- as.Date(data_all_pheno_cox7_40$time_CA_121.21)

cox_result_time <- run_cox_model_each_var(
  data = data_all_pheno_cox7_40,
  event_col = "CA_121.21",
  event_time_col = "time_CA_121.21",
  baseline_time_col = "SHIFTED_SAMPLE_DATE",
  last_record_time_col = "observation_period_end_date",
  covariates = c("age","gender","PC1_SUM","PC2_SUM","PC3_SUM","PC4_SUM","PC5_SUM"),
  variables_of_interest = c("mca_auto","auto_mca_highrisk","auto_mca_lowrisk","auto_mca_cf_l","auto_mca_cf_m")
)

cox_result_time

cox_result_time <- run_cox_model_2times(
  data = data_all_pheno_cox7_40,
  event_col = "CA_121.21",
  event_time_col = "time_CA_121.21",
  baseline_time_col = "SHIFTED_SAMPLE_DATE",
  last_record_time_col = "observation_period_end_date",
  age_blooddraw = "age",
  covariates = c("gender","PC1_SUM","PC2_SUM","PC3_SUM","PC4_SUM","PC5_SUM"),
  variables_of_interest = c("mca_auto","auto_mca_highrisk","auto_mca_lowrisk","auto_mca_cf_l","auto_mca_cf_m")
)

cox_result_time

write.csv(cox_result_time, file = "cox_mCA_CLL_biovu250k.csv")

cox_result_time <- run_cox_model_2times(
  data = data_all_pheno_cox7_40,
  event_col = "CA_121.21",
  event_time_col = "time_CA_121.21",
  baseline_time_col = "SHIFTED_SAMPLE_DATE",
  last_record_time_col = "observation_period_end_date",
  age_blooddraw = "age",
  covariates = c("gender","PC1_SUM","PC2_SUM","PC3_SUM","PC4_SUM","PC5_SUM"),
  variables_of_interest = c("mca_auto","auto_mca_highrisk")
)

cox_result_time

cox_result_time <- run_cox_model_each_var(
  data = data_all_pheno_cox7_40,
  event_col = "CA_121.21",
  event_time_col = "time_CA_121.21",
  baseline_time_col = "SHIFTED_SAMPLE_DATE",
  last_record_time_col = "observation_period_end_date",
  covariates = c("age","gender","PC1_SUM","PC2_SUM","PC3_SUM","PC4_SUM","PC5_SUM"),
  variables_of_interest = c("mca_auto","auto_mca_highrisk","auto_mca_lowrisk","auto_mca_cf_l","auto_mca_cf_m")
)

cox_result_time

cox_result_time <- run_cox_model_each_var(
  data = data_all_pheno_cox7_40,
  event_col = "CA_121.21",
  event_time_col = "time_CA_121.21",
  baseline_time_col = "SHIFTED_SAMPLE_DATE",
  last_record_time_col = "observation_period_end_date",
  covariates = c("age","gender","PC1_SUM","PC2_SUM","PC3_SUM","PC4_SUM","PC5_SUM"),
  variables_of_interest = "auto_mca_lowrisk"
)

cox_result_time

cox_result_time <- run_cox_model_each_var(
  data = data_all_pheno_cox7_40,
  event_col = "CA_121.21",
  event_time_col = "time_CA_121.21",
  baseline_time_col = "SHIFTED_SAMPLE_DATE",
  last_record_time_col = "observation_period_end_date",
  covariates = c("age","gender","PC1_SUM","PC2_SUM","PC3_SUM","PC4_SUM","PC5_SUM"),
  variables_of_interest = c("auto_mca_cf_l","auto_mca_cf_m")
)

cox_result_time

table(data_mca_auto$q_arm,data_mca_auto$p_arm)

data_mca_auto$p_count<-NA
data_mca_auto$q_count<-NA
data_mca_auto$p_count<-ifelse(data_mca_auto$p_arm=="Y"|data_mca_auto$p_arm=="T"|data_mca_auto$p_arm=="C",1,0)
data_mca_auto$q_count<-ifelse(data_mca_auto$q_arm=="Y"|data_mca_auto$q_arm=="T"|data_mca_auto$q_arm=="C",1,0)
table(data_mca_auto$p_count,data_mca_auto$q_count)

data_mca_auto$arm<-NA
data_mca_auto$arm<-ifelse(data_mca_auto$p_count ==1 & data_mca_auto$q_count ==0 ,"p",data_mca_auto$arm)
data_mca_auto$arm<-ifelse(data_mca_auto$q_count ==1 & data_mca_auto$p_count ==0 ,"q",data_mca_auto$arm)
table(data_mca_auto$arm)
data_mca_auto$chr_type <- paste(data_mca_auto$chrom,data_mca_auto$arm,sep = "")
data_mca_auto <- data_mca_auto %>% mutate(chr_type = str_replace(chr_type, "NA$", ""))

table(data_mca_auto$chr_type)

data_mca_auto$chr15gain<-NA
data_mca_auto$chr15gain<-ifelse(data_mca_auto$chr_type=="chr15" & data_mca_auto$type=="Gain",1,2)

data_mca_auto$chr_type2 <- data_mca_auto$chr_type
data_mca_auto$chr_type2 <- ifelse(data_mca_auto$chr_type2 == "chr15q", "chr15", data_mca_auto$chr_type2)

saveRDS(data_mca_auto, file = "data_mca_biovu250k.rds")

data_mca_auto <- readRDS("data_mca_biovu250k.rds")

head(data_mca_auto)

data_mca_auto$mca_type <- paste(data_mca_auto$chr_type2, data_mca_auto$type, sep = ".")

data_mca_auto_subtype <- filter(data_mca_auto, data_mca_auto$type != "Undetermined")

data_mca_auto_subtype[, mca_type := gsub("-", ".", mca_type)]

library(data.table)

# 确保是 data.table
setDT(data_mca_auto_subtype)

# 先去重（防止同一个GRID同一个mca_type重复）
dt <- unique(data_mca_auto_subtype[, .(GRID, mca_type)])

# 生成1/0矩阵
wide <- dcast(
  dt,
  GRID ~ mca_type,
  fun.aggregate = length,
  value.var = "mca_type"
)

# 把计数 >1 的情况压成1（保险）
cols <- setdiff(names(wide), "GRID")
wide[, (cols) := lapply(.SD, function(x) as.integer(x > 0)), .SDcols = cols]

wide

saveRDS(wide, file = "data_mca_BioVU250k_subtype.rds")

cols <- c(
  "chr11p.CN.LOH","chr11q.CN.LOH","chr12.Gain","chr13q.CN.LOH",
  "chr13q.Loss","chr14q.CN.LOH","chr15.CN.LOH","chr15.Gain",
  "chr17q.CN.LOH","chr1p.CN.LOH","chr1q.CN.LOH","chr20q.Loss",
  "chr22q.CN.LOH","chr3q.CN.LOH","chr4q.CN.LOH","chr6p.CN.LOH",
  "chr9p.CN.LOH","chr9q.CN.LOH"
)

wide[, lapply(.SD, sum, na.rm = TRUE), .SDcols = cols]



x <- load("mca_all_BUAT_KZ_06102025.Rdata")
x

ls(data_4cohorts)
table(data_4cohorts$chr_type2)

data_filtered <- data_mca_auto %>% filter(type != "Undetermined")

data_filtered <- filter(data_filtered, data_filtered$age >= 40)

table(data_filtered$computed_gender)

table(data_filtered$computed_gender,data_filtered$gender)

data_summary <- NA
data_summary <- data_filtered %>%
  mutate(sex_numeric = ifelse(gender == "male", 1, 0)) %>%
  group_by(type, chr_type2) %>%
  summarise(
    mean_age = mean(age, na.rm = TRUE),
    sem_age = sd(age, na.rm = TRUE) / sqrt(n()),
    mean_fraction_male = mean(sex_numeric, na.rm = TRUE),
    sem_fraction_male = sd(sex_numeric, na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = 'drop'
  )
data_summary<-filter(data_summary,data_summary$mean_fraction_male>0 & data_summary$mean_fraction_male<1)
data_summary<-filter(data_summary,data_summary$n>=10)

data_summary <- data_summary %>%
  mutate(show_label = ifelse(chr_type2 %in% c("chr15","chr20q","chr10q"), TRUE, FALSE))

data_summary<-filter(data_summary,data_summary$chr_type2 != "chrX" & data_summary$chr_type2 != "chrXp" & data_summary$chr_type2 != "chrXq")

data_summary <- data_summary %>% arrange(mean_fraction_male)
data_summary

p <- ggplot(data_summary, aes(x = mean_age, y = mean_fraction_male, color = type, size = n)) +
  geom_point() +
  geom_segment(aes(x = mean_age - sem_age, xend = mean_age + sem_age,
                   y = mean_fraction_male, yend = mean_fraction_male), linewidth = 0.5) +
  
  geom_segment(aes(x = mean_age, xend = mean_age,
                   y = mean_fraction_male - sem_fraction_male,
                   yend = mean_fraction_male + sem_fraction_male), linewidth = 0.5) +

  geom_text(data = filter(data_summary, show_label),
            aes(label = chr_type2),
            hjust = 1.1, vjust = 1.1, size = 3.5, show.legend = FALSE) +

  scale_color_manual(values = c("Gain" = "#56B4E9", "Loss" = "#E69F00", "CN-LOH" = "#009E73")) +
  scale_size_continuous(range = c(2, 8)) +

  labs(
    title = "mCA types with age and sex (BioVU)",
    x = "Mean age (years)",
    y = "Fraction male",
    size = "Sample Size",
    color = "mCA Type"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "right")
p

library(sjPlot)

save_plot("age_sex_mca_biovu.svg", fig = p, width=18, height=15)

data_mca_auto

data_mca_auto$ancestry <- dplyr::recode(
  data_mca_auto$unsupervised_ancestry_cluster_relabel,
  "Admixed_(majority_ancestry_<_0.5)" = "Admixed",
  "k1_(EAS)" = "EAS",
  "k2_(SAS)" = "SAS",
  "k3_(AFR)" = "AFR",
  "k4_(AMR)" = "AMR",
  "k5_(EUR)" = "EUR"
)
table(data_mca_auto$ancestry)

table(data_mca_auto$chr_type2)

counts <- data_mca_auto %>%
  group_by(ancestry, chr_type2,type) %>%
  summarise(count = n()) %>%
  ungroup()

total_counts <- data_mca_auto %>%
  group_by(ancestry) %>%
  summarise(total = n()) %>%
  ungroup()

counts <- counts %>%
  left_join(total_counts, by = "ancestry") %>%
  mutate(proportion = (count / total) * 100) 

combined_data <- filter(counts,counts$ancestry=="EUR" | counts$ancestry=="AFR" )
combined_data <- filter(combined_data, combined_data$type!= "Undetermined" )

ls(combined_data)
eur_data <- combined_data %>%
  dplyr::filter(ancestry == "EUR") %>%
  dplyr::select(chr_type2,type, eur_proportion = proportion, count_eur = count)

afr_data <- combined_data %>%
  dplyr::filter(ancestry == "AFR") %>%
  dplyr::select(chr_type2,type, afr_proportion = proportion, count_afr = count)

combined_data <- merge(eur_data, afr_data, by = c("chr_type2","type"))
combined_data <- filter(combined_data,combined_data$count_eur >= 5 | combined_data$count_afr >= 5)
combined_data <- filter(combined_data,combined_data$count_eur >= 2 & combined_data$count_afr >= 2)

combined_data$count <- combined_data$count_eur + combined_data$count_afr

combined_data <- combined_data %>% arrange(eur_proportion)
combined_data

p <- ggplot(combined_data, aes(x = eur_proportion, y = afr_proportion, color = type, label = chr_type2, size = count)) +
  geom_point() +  
  geom_text(vjust = -0.5, hjust = 0.5) +  
  scale_x_continuous(name = "EUR Proportion (%)", limits = c(0, 4)) +
  scale_y_continuous(name = "AFR Proportion (%)", limits = c(0, 4)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey") + 
  scale_color_manual(values = c("Gain" = "#56B4E9", "Loss" = "#E69F00", "CN-LOH" = "#009E73")) + 
  scale_size_continuous(range = c(1, 5), name = "Count") + 
  theme_minimal() +
  labs(title = "Chromosome Type Proportion in EUR vs AFR",
       subtitle = "Colored by Chromosome Event Type")
p

p <- ggplot(combined_data, aes(x = eur_proportion, y = afr_proportion, color = type, label = chr_type2, size = count)) +
  geom_point() +  
  #geom_text(vjust = -0.5, hjust = 0.5) +  
  scale_x_continuous(name = "EUR Proportion (%)", limits = c(0, 4)) +
  scale_y_continuous(name = "AFR Proportion (%)", limits = c(0, 4)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey") + 
  scale_color_manual(values = c("Gain" = "#56B4E9", "Loss" = "#E69F00", "CN-LOH" = "#009E73")) + 
  scale_size_continuous(range = c(1, 5), name = "Count") + 
  theme_minimal() +
  labs(title = "Chromosome Type Proportion in EUR vs AFR",
       subtitle = "Colored by Chromosome Event Type")
p

save_plot("ancestry_mca_BioVU250k.svg", fig = p, width=18, height=15)

x <- load("mca_all_BUAT_KZ_06102025.Rdata")
x

data_uat <- filter(data_4cohorts, data_4cohorts$cohort != "BioVU")

ls(data_uat)

names(data_mca_auto)[1] <- "sample_id"
data_mca_auto$cohort <- "BioVU"
ls(data_mca_auto)

common_cols <- intersect(names(data_mca_auto), names(data_uat))

data_4cohorts <- data.table::rbindlist(list(
  data_mca_auto[, ..common_cols],
  data_uat[, ..common_cols]
), use.names = TRUE)

saveRDS(data_4cohorts, file = "mca_all_BUAT_KZ_04012026.rds")

data_4cohorts <- readRDS("mca_all_BUAT_KZ_04012026.rds")

data_filtered <- data_4cohorts %>% filter(type != "Undetermined")
table(data_filtered$cohort, data_filtered$ancestry, useNA = "always")

summary(data_filtered$cf)

data_filtered <- filter(data_4cohorts, data_4cohorts$chrom!= "X" & data_4cohorts$chrom!= "chrX")

library(dplyr)

unique_counts_by_type <- data_filtered %>%
  group_by(cohort, type) %>%
  summarise(unique_n = n_distinct(sample_id), .groups = "drop")

unique_counts_overall <- data_filtered %>%
  group_by(cohort) %>%
  summarise(unique_n = n_distinct(sample_id)) %>%
  mutate(type = "All")

unique_counts_combined <- bind_rows(unique_counts_overall, unique_counts_by_type) %>%
  arrange(cohort, type)

print(unique_counts_combined)

355/2144

cf_stats_by_type <- data_filtered %>%
  group_by(cohort, type) %>%
  summarise(
    median_cf = median(cf, na.rm = TRUE),
    iqr_lower = quantile(cf, 0.25, na.rm = TRUE),
    iqr_upper = quantile(cf, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

cf_stats_overall <- data_filtered %>%
  group_by(cohort) %>%
  summarise(
    median_cf = median(cf, na.rm = TRUE),
    iqr_lower = quantile(cf, 0.25, na.rm = TRUE),
    iqr_upper = quantile(cf, 0.75, na.rm = TRUE)
  ) %>%
  mutate(type = "All")

cf_stats_combined <- bind_rows(cf_stats_overall, cf_stats_by_type) %>%
  arrange(cohort, type)

cf_stats_combined

data_filtered_10 <- filter(data_filtered, data_filtered$cf >= 0.1)

unique_counts_by_type <- data_filtered_10 %>%
  group_by(cohort, type) %>%
  summarise(unique_n = n_distinct(sample_id), .groups = "drop")

unique_counts_overall <- data_filtered_10 %>%
  group_by(cohort) %>%
  summarise(unique_n = n_distinct(sample_id)) %>%
  mutate(type = "All")

unique_counts_combined <- bind_rows(unique_counts_overall, unique_counts_by_type) %>%
  arrange(cohort, type)

print(unique_counts_combined)

1295/414809
4253/473067
530/104426
1891/245351

cf_stats_by_type <- data_filtered_10 %>%
  group_by(cohort, type) %>%
  summarise(
    median_cf = median(cf, na.rm = TRUE),
    iqr_lower = quantile(cf, 0.25, na.rm = TRUE),
    iqr_upper = quantile(cf, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

cf_stats_overall <- data_filtered_10 %>%
  group_by(cohort) %>%
  summarise(
    median_cf = median(cf, na.rm = TRUE),
    iqr_lower = quantile(cf, 0.25, na.rm = TRUE),
    iqr_upper = quantile(cf, 0.75, na.rm = TRUE)
  ) %>%
  mutate(type = "All")

cf_stats_combined <- bind_rows(cf_stats_overall, cf_stats_by_type) %>%
  arrange(cohort, type)

cf_stats_combined

pairwise.wilcox.test(data_filtered$age, data_filtered$type, p.adjust.method = "bonferroni")

library(rstatix)

pairwise_t_test_result <- data_filtered %>%
  pairwise_t_test(age ~ type, p.adjust.method = "bonferroni") %>%
  select(group1, group2, p, p.adj, p.adj.signif)

pairwise_t_test_result

data_filtered <- data_4cohorts %>% filter(type != "Undetermined")
table(data_filtered$cohort, data_filtered$ancestry, useNA = "always")

data_summary<-NA
data_summary <- data_filtered %>%
  mutate(sex_numeric = ifelse(computed_gender == "M", 1, 0)) %>%
  group_by(type, chr_type2) %>%
  summarise(
    mean_age = mean(age, na.rm = TRUE),
    sem_age = sd(age, na.rm = TRUE) / sqrt(n()),
    mean_fraction_male = mean(sex_numeric, na.rm = TRUE),
    sem_fraction_male = sd(sex_numeric, na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = 'drop'
  )
data_summary<-filter(data_summary,data_summary$mean_fraction_male>0 & data_summary$mean_fraction_male<1)
data_summary<-filter(data_summary,data_summary$n>=150)

data_summary <- data_summary %>%
  mutate(show_label = ifelse(chr_type2 %in% c("chr15","chr20q","chr10q"), TRUE, FALSE))

data_summary<-filter(data_summary,data_summary$chr_type2 != "chrX" & data_summary$chr_type2 != "chrXp" & data_summary$chr_type2 != "chrXq")

data_summary <- data_summary %>% arrange(mean_fraction_male)
data_summary

p <- ggplot(data_summary, aes(x = mean_age, y = mean_fraction_male, color = type, size = n)) +
  geom_point() +
  geom_segment(aes(x = mean_age - sem_age, xend = mean_age + sem_age,
                   y = mean_fraction_male, yend = mean_fraction_male), linewidth = 0.5) +
  
  geom_segment(aes(x = mean_age, xend = mean_age,
                   y = mean_fraction_male - sem_fraction_male,
                   yend = mean_fraction_male + sem_fraction_male), linewidth = 0.5) +

  geom_text(data = filter(data_summary, show_label),
            aes(label = chr_type2),
            hjust = 1.1, vjust = 1.1, size = 3.5, show.legend = FALSE) +

  scale_color_manual(values = c("Gain" = "#56B4E9", "Loss" = "#E69F00", "CN-LOH" = "#009E73")) +
  scale_size_continuous(range = c(2, 8)) +
  scale_x_continuous(limits = c(55, 75))+
  labs(
    title = "mCA types with age and sex (BUAT all)",
    x = "Mean age (years)",
    y = "Fraction male",
    size = "Sample Size",
    color = "mCA Type"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "right")
p

save_plot("age_sex_mca_BUAT.svg", fig = p, width=18, height=15)

id_chrX <- filter(data_filtered, data_filtered$chrom == "X" | data_filtered$chrom == "chrX")

dim(id_chrX)

data_filtered2 <- data_filtered %>%
  filter(!sample_id %in% id_chrX$sample_id)

data_summary<-NA
data_summary <- data_filtered2 %>%
  mutate(sex_numeric = ifelse(computed_gender == "M", 1, 0)) %>%
  group_by(type, chr_type2) %>%
  summarise(
    mean_age = mean(age, na.rm = TRUE),
    sem_age = sd(age, na.rm = TRUE) / sqrt(n()),
    mean_fraction_male = mean(sex_numeric, na.rm = TRUE),
    sem_fraction_male = sd(sex_numeric, na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = 'drop'
  )
data_summary<-filter(data_summary,data_summary$mean_fraction_male>0 & data_summary$mean_fraction_male<1)
data_summary<-filter(data_summary,data_summary$n>=110)

data_summary <- data_summary %>%
  mutate(show_label = ifelse(chr_type2 %in% c("chr15","chr20q","chr10q"), TRUE, FALSE))

data_summary<-filter(data_summary,data_summary$chr_type2 != "chrX" & data_summary$chr_type2 != "chrXp" & data_summary$chr_type2 != "chrXq")

data_summary <- data_summary %>% arrange(mean_age)
data_summary

p <- ggplot(data_summary, aes(x = mean_age, y = mean_fraction_male, color = type, size = n)) +
  geom_point() +
  geom_segment(aes(x = mean_age - sem_age, xend = mean_age + sem_age,
                   y = mean_fraction_male, yend = mean_fraction_male), linewidth = 0.5) +
  
  geom_segment(aes(x = mean_age, xend = mean_age,
                   y = mean_fraction_male - sem_fraction_male,
                   yend = mean_fraction_male + sem_fraction_male), linewidth = 0.5) +

  geom_text(data = filter(data_summary, show_label),
            aes(label = chr_type2),
            hjust = 1.1, vjust = 1.1, size = 3.5, show.legend = FALSE) +

  scale_color_manual(values = c("Gain" = "#56B4E9", "Loss" = "#E69F00", "CN-LOH" = "#009E73")) +
  scale_size_continuous(range = c(2, 8)) +
  scale_x_continuous(limits = c(55, 75))+
  labs(
    title = "mCA types with age and sex (BUAT all)",
    x = "Mean age (years)",
    y = "Fraction male",
    size = "Sample Size",
    color = "mCA Type"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "right")
p

save_plot("age_sex_mca_BUAT_nochrX.svg", fig = p, width=18, height=15)

counts <- data_filtered %>%
  group_by(ancestry, chr_type2,type) %>%
  summarise(count = n()) %>%
  ungroup()

total_counts <- data_filtered %>%
  group_by(ancestry) %>%
  summarise(total = n()) %>%
  ungroup()

counts <- counts %>%
  left_join(total_counts, by = "ancestry") %>%
  mutate(proportion = (count / total) * 100) 

combined_data <- filter(counts,counts$ancestry=="EUR" | counts$ancestry=="AFR" )
combined_data <- filter(combined_data, combined_data$type!= "Undetermined" )
combined_data <- filter(combined_data, combined_data$chr_type2!= "chrX" & combined_data$chr_type2!= "chrXp" & combined_data$chr_type2!= "chrXq")

ls(combined_data)
eur_data <- combined_data %>%
  dplyr::filter(ancestry == "EUR") %>%
  dplyr::select(chr_type2,type, eur_proportion = proportion, count_eur = count)

afr_data <- combined_data %>%
  dplyr::filter(ancestry == "AFR") %>%
  dplyr::select(chr_type2,type, afr_proportion = proportion, count_afr = count)

combined_data <- merge(eur_data, afr_data, by = c("chr_type2","type"))
combined_data <- filter(combined_data,combined_data$count_eur >= 20 | combined_data$count_afr >= 20)
combined_data$count <- combined_data$count_eur + combined_data$count_afr

eur_total <- data_filtered %>%
  filter(ancestry == "EUR") %>%
  summarise(total_eur = n())

afr_total <- data_filtered %>%
  filter(ancestry == "AFR") %>%
  summarise(total_afr = n())

combined_data <- combined_data %>%
  mutate(
    total_eur = eur_total$total_eur,
    total_afr = afr_total$total_afr
  )

library(binom)

ci_eur <- binom::binom.confint(
  combined_data$count_eur,
  combined_data$total_eur,
  method = "wilson"
)

ci_afr <- binom::binom.confint(
  combined_data$count_afr,
  combined_data$total_afr,
  method = "wilson"
)

combined_data <- combined_data %>%
  mutate(
    eur_lower = ci_eur$lower * 100,
    eur_upper = ci_eur$upper * 100,
    afr_lower = ci_afr$lower * 100,
    afr_upper = ci_afr$upper * 100
  )

combined_data <- combined_data %>% arrange(eur_proportion)
combined_data

compare_proportions <- function(eur_proportion, count_eur,
                                afr_proportion, count_afr) {
  total_eur <- round(count_eur / eur_proportion)
  total_afr <- round(count_afr / afr_proportion)

  eur_fail <- total_eur - count_eur
  afr_fail <- total_afr - count_afr

  test <- prop.test(c(count_eur, count_afr), c(total_eur, total_afr))

  list(
    total_eur = total_eur,
    total_afr = total_afr,
    eur_rate = paste0(round(100 * eur_proportion, 2), "%"),
    afr_rate = paste0(round(100 * afr_proportion, 2), "%"),
    p_value = signif(test$p.value, 3),
    estimate_diff = round(diff(test$estimate), 4),
    conf_int = round(test$conf.int, 4),
    significant = ifelse(test$p.value < 0.05, "Yes", "No"),
    method = test$method
  )
}

compare_proportions(
  eur_proportion = 0.013951855,
  count_eur = 437,
  afr_proportion = 0.029893648,
  count_afr = 104
)

-log10(0.05/3311)

p <- ggplot(combined_data, aes(x = eur_proportion, y = afr_proportion, color = type, label = chr_type2, size = count)) +
  geom_point() +  
  geom_text(vjust = -0.5, hjust = 0.5) +  
  scale_x_continuous(name = "EUR Proportion (%)", limits = c(0, 4)) +
  scale_y_continuous(name = "AFR Proportion (%)", limits = c(0, 4)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey") + 
  scale_color_manual(values = c("Gain" = "#56B4E9", "Loss" = "#E69F00", "CN-LOH" = "#009E73")) + 
  scale_size_continuous(range = c(1, 5), name = "Count") + 
  theme_minimal() +
  labs(title = "Chromosome Type Proportion in EUR vs AFR",
       subtitle = "Colored by Chromosome Event Type")
p

p <- ggplot(combined_data, aes(x = eur_proportion, y = afr_proportion, color = type, label = chr_type2, size = count)) +
  geom_point() +  
  #geom_text(vjust = -0.5, hjust = 0.5) +  
  scale_x_continuous(name = "EUR Proportion (%)", limits = c(0, 4)) +
  scale_y_continuous(name = "AFR Proportion (%)", limits = c(0, 4)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey") + 
  scale_color_manual(values = c("Gain" = "#56B4E9", "Loss" = "#E69F00", "CN-LOH" = "#009E73")) + 
  scale_size_continuous(range = c(1, 5), name = "Count") + 
  theme_minimal() +
  labs(title = "Chromosome Type Proportion in EUR vs AFR",
       subtitle = "Colored by Chromosome Event Type")
p

p <- ggplot(combined_data, aes(x = eur_proportion, y = afr_proportion, color = type, label = chr_type2, size = count)) +
  geom_point() +
  
  geom_errorbar(aes(ymin = afr_lower, ymax = afr_upper),
                width = 0,
                alpha = 0.4,        # 透明
                linewidth = 1) +   # 变细
  
  geom_errorbarh(aes(xmin = eur_lower, xmax = eur_upper),
                 height = 0,
                 alpha = 0.4,
                 linewidth = 1) +
  
  scale_x_continuous(name = "EUR Proportion (%)", limits = c(0, 4.5)) +
  scale_y_continuous(name = "AFR Proportion (%)", limits = c(0, 4.5)) +
  
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey") + 
  scale_color_manual(values = c("Gain" = "#56B4E9", "Loss" = "#E69F00", "CN-LOH" = "#009E73")) + 
  scale_size_continuous(range = c(1, 5), name = "Count") + 
  
  theme_minimal() +
  labs(title = "Chromosome Type Proportion in EUR vs AFR",
       subtitle = "Colored by Chromosome Event Type")

p

library(sjPlot)
save_plot("ancestry_mca_BUAT_CI.svg", fig = p, width=18, height=15)

library(dplyr)

suptable_overall <- combined_data %>%
  select(
    chr_type2,
    type,
    eur_proportion,
    count_eur,
    afr_proportion,
    count_afr,
    eur_lower,
    eur_upper,
    afr_lower,
    afr_upper
  ) %>%
  mutate(
    cohort = "Overall"
  )

counts_cohort <- data_filtered %>%
  group_by(cohort, ancestry, chr_type2, type) %>%
  summarise(count = n(), .groups = "drop")

total_counts_cohort <- data_filtered %>%
  group_by(cohort, ancestry) %>%
  summarise(total = n(), .groups = "drop")

counts_cohort <- counts_cohort %>%
  left_join(total_counts_cohort, by = c("cohort","ancestry")) %>%
  mutate(proportion = count / total)

eur_data <- counts_cohort %>%
  filter(ancestry == "EUR") %>%
  select(cohort, chr_type2, type,
         eur_proportion = proportion,
         count_eur = count,
         total_eur = total)

afr_data <- counts_cohort %>%
  filter(ancestry == "AFR") %>%
  select(cohort, chr_type2, type,
         afr_proportion = proportion,
         count_afr = count,
         total_afr = total)

combined_cohort <- merge(eur_data, afr_data,
                         by = c("cohort","chr_type2","type"))

library(binom)

ci_eur <- binom::binom.confint(
  combined_cohort$count_eur,
  combined_cohort$total_eur,
  method = "wilson"
)

ci_afr <- binom::binom.confint(
  combined_cohort$count_afr,
  combined_cohort$total_afr,
  method = "wilson"
)

combined_cohort <- combined_cohort %>%
  mutate(
    eur_lower = ci_eur$lower * 100,
    eur_upper = ci_eur$upper * 100,
    afr_lower = ci_afr$lower * 100,
    afr_upper = ci_afr$upper * 100,
    
    eur_proportion = eur_proportion * 100,
    afr_proportion = afr_proportion * 100
  )

suptable_cohort <- combined_cohort %>%
  select(
    cohort,
    chr_type2,
    type,
    eur_proportion,
    count_eur,
    afr_proportion,
    count_afr,
    eur_lower,
    eur_upper,
    afr_lower,
    afr_upper
  )

final_suptable <- bind_rows(suptable_overall, suptable_cohort)

final_suptable <- final_suptable %>%
  arrange(chr_type2, type, cohort)

final_suptable <- filter(final_suptable, final_suptable$type != "Undetermined")
final_suptable <- filter(final_suptable, final_suptable$cohort != "UKB")

head(final_suptable)

write.csv(final_suptable,
          "mCA_proportions_anecstry_BUAT.csv",
          row.names = FALSE)

library(dplyr)
library(ggplot2)

biovu_data <- final_suptable %>%
  filter(cohort == "BioVU") %>%
  filter(type != "Undetermined") %>%
  filter(!chr_type2 %in% c("chrX", "chrXp", "chrXq")) %>%
  mutate(count = count_eur + count_afr)

p_biovu <- ggplot(
  biovu_data,
  aes(
    x = eur_proportion,
    y = afr_proportion,
    color = type,
    label = chr_type2,
    size = count
  )
) +
  geom_point() +
  geom_errorbar(
    aes(ymin = afr_lower, ymax = afr_upper),
    width = 0,
    alpha = 0.2,
    linewidth = 0.8
  ) +
  geom_errorbarh(
    aes(xmin = eur_lower, xmax = eur_upper),
    height = 0,
    alpha = 0.2,
    linewidth = 0.8
  ) +
  scale_x_continuous(name = "EUR Proportion (%)", limits = c(0, 4.5)) +
  scale_y_continuous(name = "AFR Proportion (%)", limits = c(0, 4.5)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey") +
  scale_color_manual(values = c(
    "Gain" = "#56B4E9",
    "Loss" = "#E69F00",
    "CN-LOH" = "#009E73"
  )) +
  scale_size_continuous(range = c(1, 5), name = "Count") +
  theme_minimal() +
  labs(
    title = "Chromosome Type Proportion in EUR vs AFR (BioVU)",
    subtitle = "Colored by Chromosome Event Type"
  )

p_biovu

library(sjPlot)
save_plot("ancestry_mca_Biovu_CI.svg", fig = p_biovu, width=18, height=15)

library(dplyr)
library(ggplot2)

topmed_data <- final_suptable %>%
  filter(cohort == "TOPMed") %>%
  filter(type != "Undetermined") %>%
  filter(!chr_type2 %in% c("chrX", "chrXp", "chrXq")) %>%
  mutate(count = count_eur + count_afr)

p_topmed <- ggplot(
  topmed_data,
  aes(
    x = eur_proportion,
    y = afr_proportion,
    color = type,
    label = chr_type2,
    size = count
  )
) +
  geom_point() +
  geom_errorbar(
    aes(ymin = afr_lower, ymax = afr_upper),
    width = 0,
    alpha = 0.2,
    linewidth = 0.8
  ) +
  geom_errorbarh(
    aes(xmin = eur_lower, xmax = eur_upper),
    height = 0,
    alpha = 0.2,
    linewidth = 0.8
  ) +
  scale_x_continuous(name = "EUR Proportion (%)", limits = c(0, 5)) +
  scale_y_continuous(name = "AFR Proportion (%)", limits = c(0, 5)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey") +
  scale_color_manual(values = c(
    "Gain" = "#56B4E9",
    "Loss" = "#E69F00",
    "CN-LOH" = "#009E73"
  )) +
  scale_size_continuous(range = c(1, 5), name = "Count") +
  theme_minimal() +
  labs(
    title = "Chromosome Type Proportion in EUR vs AFR (TopMed)",
    subtitle = "Colored by Chromosome Event Type"
  )

p_topmed

library(sjPlot)
save_plot("ancestry_mca_TOPMed_CI.svg", fig = p_topmed, width=18, height=15)

data_4cohorts <- readRDS("mca_all_BUAT_KZ_04012026.rds")

data_auto_BUAT <- filter(data_4cohorts, data_4cohorts$chrom != "X" & data_4cohorts$chrom != "chrX")
head(data_auto_BUAT)

dim(data_auto_BUAT)

data_auto_BUAT <- data_auto_BUAT %>%
  add_count(sample_id, cohort, name = "Count")

table(data_auto_BUAT$cohort, data_auto_BUAT$Count)

library(dplyr)
library(ggplot2)

plot_data <- data_auto_BUAT %>%
  mutate(count_group = case_when(
    Count == 1 ~ "1",
    Count == 2 ~ "2",
    Count == 3 ~ "3",
    Count == 4 ~ "4",
    Count > 4 ~ ">4",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(count_group)) %>%
  mutate(count_group = factor(count_group, levels = c("1", "2", "3", "4", ">4")))

count_summary <- plot_data %>%
  group_by(count_group, cohort) %>%
  summarise(n = n(), .groups = "drop")

p <- ggplot(count_summary, aes(x = count_group, y = n, fill = cohort)) +
  geom_bar(stat = "identity") +
  labs(
    x = "Duplicate Count (1, 2, 3, 4, >4)",
    y = "Number of Records",
    title = "Sample Duplication Frequency by Cohort (Stacked)",
    fill = "Cohort"
  ) +
  theme_minimal(base_size = 14)
p

group_totals <- count_summary %>%
  group_by(count_group) %>%
  summarise(
    total_n = sum(n),
    label = paste0("N=", total_n, "\n", sprintf("%.1f%%", 100 * total_n / sum(n)))
  )
group_totals

save_plot("mca_count_BUAT.svg", fig = p, width=18, height=13)

x <- load("data_BUAT_all.Rdata")
x

dim(data_BUAT_all)
head(data_BUAT_all)

data_4cohorts <- readRDS("mca_all_BUAT_KZ_04012026.rds")

head(data_4cohorts)

table(data_all_pheno_cox7_40$gender)

data_mca_biovu <- data_all_pheno_cox7_40[,c("GRID","unsupervised_ancestry_cluster_relabel","gender","age","mca_auto")]

data_mca_biovu$smoking <- "Unknown"
data_mca_biovu$cohort <- "BioVU250k"
names(data_mca_biovu)[1] <- "sample_id"

data_mca_biovu$ancestry <- dplyr::recode(
  data_mca_biovu$unsupervised_ancestry_cluster_relabel,
  "Admixed_(majority_ancestry_<_0.5)" = "Admixed",
  "k1_(EAS)" = "EAS",
  "k2_(SAS)" = "SAS",
  "k3_(AFR)" = "AFR",
  "k4_(AMR)" = "AMR", 
  "k5_(EUR)" = "EUR"
)
table(data_mca_biovu$ancestry)

data_mca_biovu$sex <- dplyr::recode(
  data_mca_biovu$gender,
  "female" = "Female",
  "male" = "Male"
)
table(data_mca_biovu$sex)

data_mca_biovu <- data_mca_biovu[, c("sample_id", "ancestry", "sex", "age","smoking","mca_auto","cohort")]

data_UAT_all <- filter(data_BUAT_all, data_BUAT_all$cohort != "BioVU")

data_BUAT_all2 <- rbind(data_UAT_all,data_mca_biovu )

dim(data_BUAT_all2)
head(data_BUAT_all2)

summary(data_4cohorts$c)

# Optional: cf > 0.1
data_4cohorts_high <- filter(data_4cohorts, data_4cohorts$cf >= 0.05)

data_BUAT_all2$mca_high <- NA
data_BUAT_all2$mca_high <- ifelse(data_BUAT_all2$sample_id %in% data_4cohorts_high$sample_id, 1, 0)

table(data_BUAT_all2$mca_high, data_BUAT_all2$mca_auto)

data_BUAT_all3 <- filter(data_BUAT_all2, data_BUAT_all2$cohort == "AoU")
dim(data_BUAT_all3)

library(mgcv)
data_BUAT_anc <- data_BUAT_all3 %>% filter(ancestry %in% c("EUR", "AFR", "EAS"))

gam_model <- gam(mca_high ~ s(age, bs = "cs") + ancestry, data = data_BUAT_anc, family = binomial)

library(mgcv)
data_BUAT_anc <- data_BUAT_all2 %>% filter(ancestry %in% c("EUR", "AFR", "EAS"))

gam_model <- gam(mca_auto ~ s(age, bs = "cs") + ancestry, data = data_BUAT_anc, family = binomial)

age_points <- expand.grid(
  age = c(30, 40, 50, 60, 70, 80), 
  ancestry = factor(c("EUR", "AFR", "EAS")) 
)

predictions <- predict(gam_model, newdata = age_points, type = "response", se.fit = TRUE)
age_points$fit <- predictions$fit  
age_points$se <- predictions$se.fit 

age_points$lower <- age_points$fit - 1.96 * age_points$se
age_points$upper <- age_points$fit + 1.96 * age_points$se


p <- ggplot(age_points, aes(x = age, y = fit, color = ancestry)) +
  geom_pointrange(aes(ymin = lower, ymax = upper), size = 0.8, shape = 15, fill = "white") + 
  #geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.6, size = 0.6) +
  geom_line(size = 1) +
  scale_color_manual(
    values = c("EUR" = "cornflowerblue", "AFR" = "brown1", "EAS" = "#009E73")
  ) +
  labs(
    title = "Age vs mCA high CF by Ancestry Groups",
    x = "Age",
    y = "mCA Status Probability",
    color = "Ancestry Group"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 15),  
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.position = "right",  
    legend.title = element_text(size = 14), 
    legend.text = element_text(size = 12)  
  )
p

p <- ggplot(age_points, aes(x = age, y = fit, color = ancestry)) +
  geom_pointrange(aes(ymin = lower, ymax = upper), size = 0.8, shape = 15, fill = "white") + 
  geom_line(size = 1) +
  scale_color_manual(
    values = c("EUR" = "cornflowerblue", "AFR" = "brown1", "EAS" = "#009E73")
  ) +
  labs(
    title = "Age vs mCA Status by Ancestry Groups",
    x = "Age",
    y = "mCA Status Probability",
    color = "Ancestry Group"
  ) +
  coord_cartesian(ylim = c(0, 0.1)) + 
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 15),  
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.position = "right",  
    legend.title = element_text(size = 14), 
    legend.text = element_text(size = 12)  
  )

p

library(sjPlot)
save_plot("mca_pre_age_anc_BUAT_4080.svg", fig = p, width=18, height=15)

p <- ggplot(age_points, aes(x = age, y = fit, color = ancestry)) +
  geom_pointrange(aes(ymin = lower, ymax = upper), size = 0.8, shape = 15, fill = "white") + 
  #geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.6, size = 0.6) +
  geom_line(size = 1) +
  scale_color_manual(
    values = c("EUR" = "cornflowerblue", "AFR" = "brown1", "EAS" = "#009E73")
  ) +
  labs(
    title = "Age vs mCA Status by Ancestry Groups",
    x = "Age",
    y = "mCA Status Probability",
    color = "Ancestry Group"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 15),  
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.position = "right",  
    legend.title = element_text(size = 14), 
    legend.text = element_text(size = 12)  
  )
p

library(sjPlot)
save_plot("mca_pre_age_anc_BUAT.svg", fig = p, width=18, height=15)

data_filtered <- data_BUAT_all3 %>%
  filter(!is.na(age), ancestry %in% c("EUR", "AFR", "EAS")) %>%
  mutate(
    age_bin = cut(
      age,
      breaks = c(-Inf, 40, 50, 60, 70, 80, Inf),
      labels = c("<40", "40–50", "50–60", "60–70", "70–80", "80+"),
      right = FALSE
    )
  )

age_mca_summary <- data_filtered %>%
  group_by(ancestry, age_bin) %>%
  summarise(
    n = n(),
    mca_count = sum(mca_high == 1, na.rm = TRUE),
    prevalence = mca_count / n,
    .groups = "drop"
  )

p <- ggplot(age_mca_summary, aes(x = age_bin, y = prevalence, group = ancestry, color = ancestry)) +
  geom_line(aes(group = ancestry), size = 1.2) +
  geom_point(size = 2) +
  labs(
    x = "Age Group",
    y = "mCA Prevalence",
    title = "mCA Prevalence by Age and Ancestry",
    color = "Ancestry"
  ) +
  theme_minimal(base_size = 14)
p

library(sjPlot)
save_plot("mca_pre_age_anc_highCF.svg", fig = p, width=18, height=15)

dim(data_all_pheno_cox7)

summary(data_all_pheno_cox7$age)
table(data_all_pheno_cox7$gender)

smoking <- fread("agd_ever_smoker_df.tsv")
data_all_pheno_cox7 <- merge(data_all_pheno_cox7, smoking, by = "person_id", all.x = T)
table(data_all_pheno_cox7$ever_smoker, useNA = "always")

data_all_pheno_cox7 <- data_all_pheno_cox7 %>%
  mutate(
    ever_smoker = ifelse(is.na(ever_smoker), "FALSE", ever_smoker)
  )

table(data_all_pheno_cox7$ever_smoker, useNA = "always")

120124/(120124+71916)

245351*0.374

data_all_pheno_mca <- filter(data_all_pheno_cox7, data_all_pheno_cox7$mca_auto == 0)
dim(data_all_pheno_mca)

summary(data_all_pheno_mca$age)
table(data_all_pheno_mca$gender)
sd(data_all_pheno_mca$age)
table(data_all_pheno_mca$ever_smoker, useNA = "always")

3089/243207

table(data_all_pheno_mca$ever_smoker)

119324/(71229+119324)

90960/243207

saveRDS(data_all_pheno_cox7_40, file = "data_phewas_mca_biovu_40.rds")

system("gsutil cp data_phewas_mca_biovu_40.rds gs://bicklab-main-storage/Users/Kun_Zhao/", intern = TRUE)

disease_cols <- colnames(data_all_pheno_cox7_40)[8:3584]
disease_cols

table(data_all_pheno_cox7_40$unsupervised_ancestry_cluster_relabel)

data_all_pheno_cox7_40_eur <- filter(data_all_pheno_cox7_40, data_all_pheno_cox7_40$unsupervised_ancestry_cluster_relabel == "k5_(EUR)")

dim(data_all_pheno_cox7_40_eur)

saveRDS(data_all_pheno_cox7_40_eur, file = "data_phewas_mca_biovu_40_eur.rds")

system("gsutil cp data_phewas_mca_biovu_40_eur.rds gs://bicklab-main-storage/Users/Kun_Zhao/", intern = TRUE)

data_all_pheno_cox7_40_afr <- filter(data_all_pheno_cox7_40, data_all_pheno_cox7_40$unsupervised_ancestry_cluster_relabel == "k3_(AFR)")

dim(data_all_pheno_cox7_40_afr)

saveRDS(data_all_pheno_cox7_40_afr, file = "data_phewas_mca_biovu_40_afr.rds")

system("gsutil cp data_phewas_mca_biovu_40_afr.rds gs://bicklab-main-storage/Users/Kun_Zhao/", intern = TRUE)

mca_subtype <- readRDS("data_mca_BioVU250k_subtype.rds")

dim(mca_subtype)

table(mca_subtype$chr20q.Loss, useNA = "always")
table(mca_subtype$chr15.Gain, useNA = "always")

data_all_pheno_cox8_40 <- merge(data_all_pheno_cox7_40, mca_subtype, by = "GRID", all.x = T)

mca_cols <- grep("^chr", names(data_all_pheno_cox8_40), value = TRUE)

data_all_pheno_cox8_40[, (mca_cols) := lapply(.SD, function(x) {
  is_na <- is.na(x)      
  x[x == 0] <- NA        
  x[is_na] <- 0         
  
  return(x)
}), .SDcols = mca_cols]

table(data_all_pheno_cox8_40$chr1.CN.LOH, useNA = "ifany")

saveRDS(data_all_pheno_cox8_40, file = "data_phewas_mca_spe_biovu_40.rds")

data_all_pheno_cox8_40 <- readRDS("data_phewas_mca_spe_biovu_40.rds")

table(data_all_pheno_cox8_40$mlox)

system("gsutil cp data_phewas_mca_spe_biovu_40.rds gs://bicklab-main-storage/Users/Kun_Zhao/", intern = TRUE)

system("gsutil ls gs://fc-secure-343d99f4-8f3b-46da-ac8c-0f9f75cc97b3/submissions/5c3c1932-82c8-42eb-9e6c-e478cc6cc799/**/phewas_*.csv", intern = TRUE)

system("gsutil -m cp 'gs://fc-secure-343d99f4-8f3b-46da-ac8c-0f9f75cc97b3/submissions/**/phewas*.csv' .")

system("gsutil cp gs://fc-secure-343d99f4-8f3b-46da-ac8c-0f9f75cc97b3/submissions/01dd334f-ebae-49e5-bf59-d9f64bb86c76/run_phewas_cox/41414e5e-da6e-4051-950a-68d6e2807d1d/call-Phewas/phewas_mca_biovu250k_age_AFR.csv .")

data_all_pheno_cox8_40_mlox <- filter(data_all_pheno_cox8_40, data_all_pheno_cox8_40$gender == "female")

table(data_all_pheno_cox8_40_mlox$mlox)

dim(data_all_pheno_cox8_40_mlox)

saveRDS(data_all_pheno_cox8_40_mlox, file = "data_all_pheno_cox8_40_mlox.rds")

system("gsutil cp data_all_pheno_cox8_40_mlox.rds gs://bicklab-main-storage/Users/Kun_Zhao/", intern = TRUE)

library(data.table)

# EUR
f1 <- "phewas_mca_biovu250k/phewas_mca_biovu250k_age_EUR.csv"
dt1 <- fread(f1)
dt1[, term := "mca_eur"]
fwrite(dt1, f1)

# AFR
f2 <- "phewas_mca_biovu250k/phewas_mca_biovu250k_age_AFR.csv"
dt2 <- fread(f2)
dt2[, term := "mca_afr"]
fwrite(dt2, f2)

files <- list.files("phewas_mca_biovu250k", pattern = "\\.csv$", full.names = TRUE)

lapply(files[1:5], function(f) names(fread(f)))

phecodex <- fread("phecodeX_ICD_CM_map_flat.csv")
head(phecodex)

phecodex <- phecodex[,c("phecode_string","phecode")]
distinct <- distinct(phecodex)

df_all <- rbindlist(lapply(files, function(f) {
  dt <- fread(f)
  
  setnames(dt, "disease", "Disease", skip_absent = TRUE)
  
  dt[, Disease := as.character(Disease)]
  phecodex[, phecode := as.character(phecode)]
  
  dt <- merge(
    dt,
    phecodex[, .(phecode, phecode_string)],
    by.x = "Disease",
    by.y = "phecode",
    all.x = TRUE
  )
  
  dt[, Disease := phecode_string]
  dt[, phecode_string := NULL]
  
  return(dt)
}), fill = TRUE)

str(df_all)

df_all <- unique(df_all, by = c("Disease", "term"))

library(data.table)

setDT(df_all)

# 1️⃣ 类型转换（很重要）
df_all[, HR := as.numeric(HR)]
df_all[, CI_lower := as.numeric(CI_lower)]
df_all[, CI_upper := as.numeric(CI_upper)]

# 2️⃣ 计算 beta（= log(HR)）
df_all[, estimate := log(HR)]

# 3️⃣ 计算标准误（基于 log CI）
df_all[, std.error := (log(CI_upper) - log(CI_lower)) / (2 * 1.96)]

# 4️⃣ 计算 z statistic
df_all[, statistic := estimate / std.error]

# 5️⃣ rename 其他列
setnames(df_all,
         old = c("p_value", "N"),
         new = c("p.value", "total_samples"),
         skip_absent = TRUE)

# 6️⃣ 过滤不合法值（强烈建议）
df_all <- df_all[
  is.finite(estimate) &
  is.finite(std.error) &
  std.error > 0
]

# 7️⃣ 最终数据
df_final <- df_all[, .(
  Disease,
  term,
  estimate,        # ✔ 现在是 log(HR)
  std.error,       # ✔ 对应 log-scale
  statistic,
  p.value,
  total_samples,
  log_p_value
)]

df_final <- df_all[, .(
  Disease,
  term,
  estimate,
  std.error,
  statistic,
  p.value,
  total_samples,
  log_p_value
)]

summary(df_final$estimate)
summary(df_final$std.error)
summary(df_final$statistic)

df_final

df_final <- df_final %>% arrange(desc(log_p_value))
df_final

setDT(df_final)

# 按 term 分组
terms <- unique(df_final$term)

# 循环写文件
for (t in terms) {
  dt_sub <- df_final[term == t]
  
  fwrite(
    dt_sub,
    file = file.path("phewas_biovu250k_mca", paste0(t, "_Ntotal_biovu250k.csv"))
  )
}

files <- list.files("phewas_mca_biovu250k", pattern = "\\.csv$", full.names = TRUE)

phecodex <- fread("phecodeX_ICD_CM_map_flat.csv")
head(phecodex)

df_all <- rbindlist(lapply(files, function(f) {
  dt <- fread(f)
  
  setnames(dt, "disease", "Disease", skip_absent = TRUE)
  
  dt[, Disease := as.character(Disease)]
  phecodex[, phecode := as.character(phecode)]
  
  dt <- merge(
    dt,
    phecodex[, .(phecode, phecode_string)],
    by.x = "Disease",
    by.y = "phecode",
    all.x = TRUE
  )
  
  dt[, Disease := phecode_string]
  dt[, phecode_string := NULL]
  
  return(dt)
}), fill = TRUE)

df_all <- unique(df_all, by = c("Disease", "term"))
df_all

setDT(df_all)

# 按 term 分组
terms <- unique(df_all$term)

# 循环写文件
for (t in terms) {
  dt_sub <- df_all[term == t]
  
  fwrite(
    dt_sub,
    file = file.path("phewas_biovu250k_mca", paste0(t, "_biovu250k.csv"))
  )
}

test <- fread("phewas_biovu250k_mca/chr11p.CN.LOH_biovu250k.csv")
test

folder <- "phewas_biovu250k_mca"

files <- list.files(folder, pattern = "^results_.*\\.csv$", full.names = TRUE)

for (f in files) {
  fname <- basename(f)
  
  # 去掉 "results_" 和 ".csv"
  core <- sub("^results_", "", fname)
  core <- sub("\\.csv$", "", core)
  
  # 新名字
  new_name <- paste0(core, "_AoU.csv")
  
  # 重命名
  file.rename(f, file.path(folder, new_name))
}

files <- system('gsutil -m ls "gs://fc-secure-343d99f4-8f3b-46da-ac8c-0f9f75cc97b3/submissions/**/phewas_ukb_*.csv"', intern = TRUE)

head(files)
length(files)

writeLines(files, "file_list.txt")

system("gsutil -m cp -I phewas_biovu250k_mca/ < file_list.txt")

folder <- "phewas_biovu250k_mca"

files <- list.files(folder, pattern = "^phewas_ukb_.*\\.csv$", full.names = TRUE)

for (f in files) {
  fname <- basename(f)
  
  # 去掉 "results_" 和 ".csv"
  core <- sub("^phewas_ukb_.", "", fname)
  core <- sub("\\.csv$", "", core)
  
  # 新名字
  new_name <- paste0(core, "_ukb.csv")
  
  # 重命名
  file.rename(f, file.path(folder, new_name))
}

files <- list.files(folder, pattern = "^hr.*_ukb\\.csv$", full.names = TRUE)

for (f in files) {
  fname <- basename(f)
  
  # 新文件名：前面加 c
  new_name <- paste0("c", fname)
  
  file.rename(f, file.path(folder, new_name))
}

library(logistf)
library(dplyr)
library(broom)

data_all_pheno_cox8_40 <- readRDS("data_phewas_mca_spe_biovu_40.rds")

result_phewas_cox <- read.csv("meta_spe_mCA_BUAT.csv")

result_phewas_cox

table(result_phewas_cox$mca_type)

dim(result_phewas_cox)
head(result_phewas_cox)

phecodex <- fread("phecodeX_ICD_CM_map_flat.csv")

phecodex <- phecodex[,c(4,5)]
phecodex <- distinct(phecodex)

head(phecodex)

result_phewas_cox <- merge(result_phewas_cox, phecodex, by.x = "phenotype", by.y = "phecode_string", all = F)
dim(result_phewas_cox)

phenotype_sig <- result_phewas_cox %>% distinct(phecode)

table(result_phewas_cox$mca_type)

target <- filter(result_phewas_cox, result_phewas_cox$mca_type == "chr9p.CN.LOH")

target

phenotype_sig <- target$phecode
length(phenotype_sig)

phenotype_sig

phenotype_sig <- phenotype_sig[phenotype_sig %in% colnames(data_all_pheno_cox8_40)]
length(phenotype_sig)

target_mca_types <- "chr9p.CN.LOH"

#-PheWAS mCA by Firth Logistic regression-----
disease_cols <- phenotype_sig

firth_results_list <- list()

for (i in seq_along(target_mca_types)) {
  mca_type <- target_mca_types[i]
  
  cat("Processing MCA Type:", mca_type, "（Loop", i, "of", length(target_mca_types), "）\n")
  
  firth_results <- list()
  
  for (j in seq_along(disease_cols)) {
    disease_col <- disease_cols[j]
    
    cat(" Processing Disease:", disease_col, "（Sub-loop", j, "of", length(disease_cols), "）\n")
    
    unique_values <- unique(data_all_pheno_cox8_40[[disease_col]])
    if (!all(unique_values %in% c(0, 1))) {
      warning(paste("Column", disease_col, "contains non-binary values. Skipping."))
      next
    }
    
    data_logistf <- data_all_pheno_cox8_40 %>%
      filter(!is.na(!!sym(disease_col)), !is.na(!!sym(mca_type)))  # Remove NA
    
    # Change mCA types to numeric
    if (!is.numeric(data_logistf[[mca_type]])) {
      data_logistf[[mca_type]] <- as.numeric(data_logistf[[mca_type]] != 0)
    }
    
    total_cases <- sum(data_logistf[[disease_col]] == 1)
    
    # Model
    formula <- as.formula(paste0("`", disease_col, "` ~ ", mca_type,
                             "+ age+gender+PC1_SUM+PC2_SUM+PC3_SUM+PC4_SUM+PC5_SUM"))
    # Firth logistic
    fit_result <- tryCatch({
      model <- logistf::logistf(
        formula,
        data = data_logistf,
        control = logistf.control(
          maxit = 15,
          lconv = 1e-04,
          gconv = 1e-04
        )
      )
        
      ci_index <- which(names(model$coefficients) == mca_type)

      tibble::tibble(
        term = mca_type,
        estimate = model$coefficients[ci_index],
        std_error = model$se[ci_index],
        OR = exp(model$coefficients[ci_index]),
        conf_low = exp(model$ci.lower[ci_index]),
        conf_high = exp(model$ci.upper[ci_index]),
        p_value = model$prob[ci_index],
        log_p_value = -log10(model$prob[ci_index]),
        total_cases = total_cases
      )
    }, error = function(e) {
      warning(paste("Error fitting model for", disease_col, "-", e$message))
      return(NULL)
    })
    
    if (!is.null(fit_result)) {
      firth_results[[disease_col]] <- fit_result
    }
  }
  
  if (length(firth_results) > 0) {
    combined_result <- bind_rows(firth_results, .id = "Disease")
    
    # Save
    output_file <- paste0("results_spe_mCA_phewas_firth_log/", mca_type, "_FirthLogis_BioVU250k.csv")
    write.csv(combined_result, output_file, row.names = FALSE)
    
    firth_results_list[[mca_type]] <- combined_result
  } else {
    warning(paste("No valid results for MCA Type:", mca_type))
  }
}

cat("Firth logistic regression completed. Results saved in 'results_spe_mCA_phewas_firth_log' directory.\n")

combined_result

x <- load("data_all_pheno_cox_BioVU.Rdata")
x

dim(data_all_pheno_cox)

table(data_all_pheno_cox$mca_status)
table(data_all_pheno_cox7$mca_auto)

data_all_pheno_cox$mca_status

mca_array <- data_all_pheno_cox[,c("person_id","mca_status")]
names(mca_array)[2] <- "mca_array"

data_compare <- merge(
  data_all_pheno_cox7,
  mca_array,
  by.x = "GRID",
  by.y = "person_id",
  all = FALSE
)

dim(data_compare$mca_auto)

data_compare40 <- filter(data_compare, data_compare$age >= 40)

table(data_compare40$mca_auto, data_compare40$mca_array)

library(epiR)

tab <- table(
  test  = factor(data_compare40$mca_array, levels = c(1, 0)),
  truth = factor(data_compare40$mca_auto,  levels = c(1, 0))
)

epiR::epi.tests(tab)

df_long <- data_compare %>%
  tidyr::pivot_longer(
    cols = c(mca_auto, mca_array),
    names_to = "method",
    values_to = "mCA"
  )

lm <- lme4::glmer(
  mCA ~ method + age + gender + (1 | person_id),
  family = binomial,
  data = df_long
)

summary(lm)


