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

install.packages("vcfR")

library(vcfR)

vcf <- read.vcfR("tert.vcf")

gt <- extract.gt(vcf, element = "GT")

dim(gt)

snp_row <- gt[1, ]

tert <- data.frame(
  id = colnames(gt),           
  tert = as.character(snp_row)
)

head(tert)

table(tert$tert, useNA = "always")

tert$id <- sub("_.*", "", tert$id)

head(tert)

names(tert)[1] <- "ID_VUMC"

x <- load("UKB_data_tet2_specific_0320.Rdata")
x

dim(data_all_pheno_cox2)
table(data_all_pheno_cox2$mca_status,useNA = "always")
summary(data_all_pheno_cox2$baseline_age)

tert$ID_VUMC <- as.character(tert$ID_VUMC)
data_all_pheno_cox2$ID_VUMC <- as.character(data_all_pheno_cox2$ID_VUMC)

data_all_tert <- merge(data_all_pheno_cox2, tert, by = "ID_VUMC", all = F)
dim(data_all_tert)

table(data_all_tert$tert)

library(mgcv)

gam_model <- gam(mca_status ~ s(baseline_age, bs = "cs") + tert, data = data_all_tert, family = binomial)

age_points <- expand.grid(
  baseline_age = c(40, 50, 60, 70, 80), 
  tert = factor(c("0/0", "0/1", "1/1"), levels = c("0/0", "0/1", "1/1")) 
)

predictions <- predict(gam_model, newdata = age_points, type = "response", se.fit = TRUE)
age_points$fit <- predictions$fit  
age_points$se <- predictions$se.fit 

age_points$lower <- age_points$fit - 1.96 * age_points$se
age_points$upper <- age_points$fit + 1.96 * age_points$se

print(age_points)

p <- ggplot(age_points, aes(x = baseline_age, y = fit, color = tert)) +
  geom_pointrange(aes(ymin = lower, ymax = upper), size = 0.8, shape = 15, fill = "white") + 
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.6, size = 0.6) +
  geom_line(size = 1) +
  scale_color_manual(
    values = c("0/0" = "#009E73", "0/1" = "#56B4E9", "1/1" = "#E69F00")
  ) +
  labs(
    title = "Age vs mCA Status by TERT Genotypes (UKB)",
    x = "Age",
    y = "mCA Status Probability",
    color = "TERT Genotypes"
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

system("gsutil cp gs://bicklab-main-storage/Users/Yash_Pershad/telo_vgwas_pheno_cov_v2.tsv .")

mltl <- fread("telo_vgwas_pheno_cov_v2.tsv")

dim(mltl)
head(mltl)

summary(mltl$telo_rint)

summary(mltl$telo_rint_resid_sq)

a <- load("mCA_clean_BUA.Rdata")
a

mCA_clean_BUA$chrom <- gsub("chr", "", mCA_clean_BUA$chrom)
table(mCA_clean_BUA$chrom,useNA="always")

mCA_clean_BUA$p_count<-NA
mCA_clean_BUA$q_count<-NA
mCA_clean_BUA$p_count<-ifelse(mCA_clean_BUA$p_arm=="Y"|mCA_clean_BUA$p_arm=="T"|mCA_clean_BUA$p_arm=="C",1,0)
mCA_clean_BUA$q_count<-ifelse(mCA_clean_BUA$q_arm=="Y"|mCA_clean_BUA$q_arm=="T"|mCA_clean_BUA$q_arm=="C",1,0)
table(mCA_clean_BUA$p_count,mCA_clean_BUA$q_count)
mCA_clean_BUA$arm<-NA
mCA_clean_BUA$arm<-ifelse(mCA_clean_BUA$p_count ==1 & mCA_clean_BUA$q_count ==0 ,"p",mCA_clean_BUA$arm)
mCA_clean_BUA$arm<-ifelse(mCA_clean_BUA$q_count ==1 & mCA_clean_BUA$p_count ==0 ,"q",mCA_clean_BUA$arm)
table(mCA_clean_BUA$arm)
mCA_clean_BUA$chr_type <- paste(mCA_clean_BUA$chrom,mCA_clean_BUA$arm,sep = "")
mCA_clean_BUA <- mCA_clean_BUA %>% mutate(chr_type = str_replace(chr_type, "NA$", ""))


table(mCA_clean_BUA$chr_type, useNA="always")

table(mCA_clean_BUA$chr_type, useNA="always")
mCA_clean_BUA$type <- ifelse(is.na(mCA_clean_BUA$type), "NA", mCA_clean_BUA$type)
table(mCA_clean_BUA$type, useNA="always")
mCA_clean_BUA$mca_type <- paste(mCA_clean_BUA$chr_type,mCA_clean_BUA$type,sep="_")
table(mCA_clean_BUA$mca_type, useNA="always")

mCA_clean_BUA$chr_type2 <- mCA_clean_BUA$chr_type
mCA_clean_BUA$chr_type2 <- ifelse(mCA_clean_BUA$chr_type2=="chr15q","chr15",mCA_clean_BUA$chr_type2)
mCA_clean_BUA$mca_type2 <- paste(mCA_clean_BUA$chr_type2,mCA_clean_BUA$type,sep="_")
table(mCA_clean_BUA$mca_type2)

mCA_clean_BUA$mca_type <- factor(mCA_clean_BUA$mca_type, levels = unique(mCA_clean_BUA$mca_type))
mCA_clean_BUA$mca_type <- relevel(mCA_clean_BUA$mca_type, ref = "NA_NA")

mCA_clean_BUA$mca_type2 <- factor(mCA_clean_BUA$mca_type2, levels = unique(mCA_clean_BUA$mca_type2))
mCA_clean_BUA$mca_type2 <- relevel(mCA_clean_BUA$mca_type2, ref = "NA_NA")

table(mCA_clean_BUA$type, useNA = "always")

mCA_clean_BUA <- filter(mCA_clean_BUA, is.na(mCA_clean_BUA$type) | mCA_clean_BUA$type != "Undetermined")
table(mCA_clean_BUA$type,useNA = "always")

mCA_clean_BUA <- filter(mCA_clean_BUA, is.na(mCA_clean_BUA$chrom) |mCA_clean_BUA$chrom != "X")
table(mCA_clean_BUA$chrom,useNA = "always")
ls(mCA_clean_BUA)

save(mCA_clean_BUA, file = "mCA_clean_BUA.Rdata")

x <- load("mCA_clean_BUA.Rdata")
x

mCA_clean_ukb <- filter(mCA_clean_BUA, mCA_clean_BUA$cohort == "UKB")

dim(mCA_clean_ukb)
head(mCA_clean_ukb)

colnames(mCA_clean_ukb)

mCA_clean_ukb <- mCA_clean_ukb[,c(1,5,10,24,26)]

mltl$IID <- as.character(mltl$IID)

mltl_ukb <- merge(mltl, mCA_clean_ukb, by.x = "IID", by.y = "person_id", all = F)
dim(mltl_ukb)

colnames(mltl_ukb)

head(mltl_ukb)

mltl_ukb$age_squ <- mltl_ukb$baseline_age * mltl_ukb$baseline_age

model <- glm(mca_auto ~ telo_rint + baseline_age + age_squ +genetic_sex + smoking + PC1+PC2+PC3+PCD4+PC5, data = mltl_ukb, family = binomial)
summary(model)

mca_types <- unique(mltl_ukb$mca_type2)
mca_types <- mca_types[mca_types != "NA_NA"]
mca_types

table(mltl_ukb$mca_type2)

results_summary <- data.frame(
  mca_type = character(),
  Estimate = numeric(),
  Std_Error = numeric(),
  z_value = numeric(),
  p_value = numeric(),
  N = integer(),
  stringsAsFactors = FALSE
)

mltl_ukb$mca_type2 <- as.character(mltl_ukb$mca_type2)

mca_types <- unique(mltl_ukb$mca_type2)
mca_types <- mca_types[mca_types != "NA_NA"]

for (mca in mca_types) {
  temp_data <- mltl_ukb %>%
    filter(mca_type2 == "NA_NA" | mca_type2 == !!mca) %>%
    mutate(mca_binary = as.integer(mca_type2 == !!mca))
  
  N_mca <- sum(temp_data$mca_binary, na.rm = TRUE)
  
  print(paste("Analyzing mca_type:", mca, "N =", N_mca))
  
  if (nrow(temp_data) > 0 && N_mca > 0) {
    model <- glm(mca_binary ~ telo_rint + baseline_age + age_squ + genetic_sex + 
                   smoking + PC1 + PC2 + PC3 + PCD4 + PC5, 
                 data = temp_data, family = binomial)
    model_summary <- summary(model)
    
    estimate <- model_summary$coefficients["telo_rint", "Estimate"]
    std_error <- model_summary$coefficients["telo_rint", "Std. Error"]
    z_value <- model_summary$coefficients["telo_rint", "z value"]
    p_value <- model_summary$coefficients["telo_rint", "Pr(>|z|)"]
    
    results_summary <- rbind(results_summary, 
                             data.frame(mca_type = mca, 
                                        Estimate = estimate, 
                                        Std_Error = std_error, 
                                        z_value = z_value, 
                                        p_value = p_value,
                                        N = N_mca))
  } else {
    warning(paste("No sufficient data for mca_type:", mca))
  }
}

results_summary <- results_summary %>%
  arrange(p_value)

print(results_summary)

results_summary

write.csv(results_summary, "Measureed_ltl_mca_ukb.csv", row.names = FALSE)

library(dplyr)
library(logistf)

mca_types2 <- mca_types[157:178]
mca_types <- mca_types2

length(mca_types)

results_summary <- data.frame(
  mca_type = character(),
  Estimate = numeric(),
  Std_Error = numeric(),
  OR = numeric(),
  CI_Lower = numeric(),
  CI_Upper = numeric(),
  p_value = numeric(),
  N = integer(),
  stringsAsFactors = FALSE
)

mltl_ukb$mca_type2 <- as.character(mltl_ukb$mca_type2)

for (mca in mca_types) {
  temp_data <- mltl_ukb %>%
    filter(mca_type2 == "NA_NA" | mca_type2 == !!mca) %>%
    mutate(mca_binary = as.integer(mca_type2 == !!mca))
  
  N_mca <- sum(temp_data$mca_binary, na.rm = TRUE)
  
  print(paste("Analyzing mca_type:", mca, "N =", N_mca))
  
  if (nrow(temp_data) > 0 && N_mca > 0) {
    fit <- logistf(
      formula = mca_binary ~ telo_rint + baseline_age + age_squ + genetic_sex + 
        smoking + PC1 + PC2 + PC3 + PCD4 + PC5, 
      data = temp_data
    )
    
    if ("telo_rint" %in% names(fit$coefficients)) {
    estimate <- fit$coefficients["telo_rint"]
    std_error <- if (!is.null(dimnames(fit$var))) sqrt(fit$var["telo_rint", "telo_rint"]) else NA
    OR <- exp(estimate)
    CI <- tryCatch(exp(confint(fit)["telo_rint", ]), error = function(e) c(NA, NA))
    p_value <- fit$prob["telo_rint"]
        } else {
    estimate <- std_error <- OR <- p_value <- CI <- NA
    warning(paste("Variable telo_rint not found in model for mca_type:", mca))
        }
    
    results_summary <- rbind(results_summary, 
                             data.frame(mca_type = mca, 
                                        Estimate = estimate, 
                                        Std_Error = std_error, 
                                        OR = OR,
                                        CI_Lower = CI[1],
                                        CI_Upper = CI[2],
                                        p_value = p_value,
                                        N = N_mca))
  } else {
    warning(paste("No sufficient data for mca_type:", mca))
  }
}

results_summary <- results_summary %>%
  arrange(p_value)

print(results_summary)

results_summary <- results_summary %>%
  arrange(p_value)
results_summary

write.csv(results_summary, "firth_mltl_mca_results_ukb_3.csv", row.names = FALSE)

mltl_result1 <- read.csv("firth_mltl_mca_results_ukb_1.csv")
mltl_result2 <- read.csv("firth_mltl_mca_results_ukb_2.csv")

mltl_result_all <- rbind(mltl_result1, mltl_result2, results_summary)

mltl_result_all$logP <- -log10(mltl_result_all$p_value)

mltl_result_all <- mltl_result_all %>% arrange(p_value)
mltl_result_all

- log10(0.05/178)

write.csv(mltl_result_all, "firth_mltl_mca_results_ukb_all.csv", row.names = FALSE)
