##-----Mediation Analysis------##
#Load packages
library(tidyr)
library(dplyr)
library(tidyverse)
library(data.table)
library(lubridate)
library(broom)
library(survival)
#UKB
load("Users_Kun_Zhao_UKB_mCA_clean2.Rdata")
cll_prs_ukb <- fread('Data_UKBioBank_prs_cll_prs_yp09302024.tsv') 
UKB_mCA_clean <- merge(UKB_mCA_clean,cll_prs_ukb, by.x = "ID_VUMC", by.y = "IID",all.x = T)
summary(UKB_mCA_clean$cll_prs)
UKB_mCA_clean <- filter(UKB_mCA_clean, is.na(UKB_mCA_clean$cll_prs)==F)
#Merge AoU
load("data_all_pheno_cox_PhecodeX.Rdata")
data_ancestry <- fread('ancestry_preds.tsv') 
data_all_pheno_cox<- merge(data_all_pheno_cox,data_ancestry,by.x = "person_id", by.y = "research_id",all.x = T)
table(data_all_pheno_cox$ancestry_pred,useNA = "always")
data_eur <- filter(data_all_pheno_cox , data_all_pheno_cox$ancestry_pred == "eur")

CLL_prs_aou <- fread('cll_pgs_v7_yp09182024.tsv') 
data_eur <- merge(data_eur,CLL_prs_aou,by.x = "person_id", by.y = "IID", all.x=T)
summary(data_eur$cll_prs)

data_all_pheno_cox2 <- filter(data_eur, is.na(data_eur$cll_prs)==F)
data_all_pheno_cox2 <- filter(data_all_pheno_cox2,data_all_pheno_cox2$age>=40 & data_all_pheno_cox2$age<=90)
summary(data_all_pheno_cox2$age)

data_all_pheno_cox2$mca_auto <- data_all_pheno_cox2$mca_status
data_all_pheno_cox2$mca_auto <- ifelse(data_all_pheno_cox2$chrom == "chrX" & data_all_pheno_cox2$mca_status == "1" , 0 ,data_all_pheno_cox2$mca_auto)
table(data_all_pheno_cox2$mca_auto)
table(data_all_pheno_cox2$mca_status)

AoU_mCA_clean <- data_all_pheno_cox2[,c("person_id", "age", "age_squ", "sex_at_birth", "smoking", "cll_prs", "ltl_prs",
                                        "cll_262_coh", "ltl_262_coh", "mca_auto", "cf_cat", "cf_max","chrom", "mca_highrisk",
                                        "p_arm","q_arm","type","cll")]

UKB_mCA_clean2 <- UKB_mCA_clean[,c("ID_VUMC", "baseline_age", "age2", "genetic_sex","eversmoked_0","cll_prs", "ltl_prs",
                                   "cll_262_coh", "ltl_262_coh", "mca_status", "cf_cat", "cf_max","chrom", "mca_highrisk",
                                   "p_arm","q_arm","type","cll")]


names(UKB_mCA_clean2)[1] <- "person_id"
names(UKB_mCA_clean2)[2] <- "age"
names(UKB_mCA_clean2)[3] <- "age_squ"
names(UKB_mCA_clean2)[4] <- "sex_at_birth"
names(UKB_mCA_clean2)[5] <- "smoking"
names(UKB_mCA_clean2)[10] <- "mca_auto"

UKB_mCA_clean2$cf_cat <- ifelse(is.na(UKB_mCA_clean2$cf_cat),0,UKB_mCA_clean2$cf_cat)
UKB_mCA_clean2$mca_highrisk <- ifelse(is.na(UKB_mCA_clean2$mca_highrisk),0,UKB_mCA_clean2$mca_highrisk)
table(UKB_mCA_clean2$cf_cat, useNA = "always")
table(UKB_mCA_clean2$mca_highrisk, useNA = "always")

UKB_mCA_carriers <- filter(UKB_mCA_clean2, UKB_mCA_clean2$mca_auto ==1)
summary(UKB_mCA_carriers$cf_max)

UKB_mCA_clean2$cohort <- "UKB"
AoU_mCA_clean$cohort <- "AoU"
mCA_clean_A_U <- rbind(UKB_mCA_clean2, AoU_mCA_clean)

#Mediation
install.packages("mediation")
library(mediation)

model.m <- glm(mca_auto ~ cll_prs + age + age_squ+ sex_at_birth+ smoking, 
               data = mCA_clean_A_U, family = binomial)
summary(model.m)

model.y <- glm(cll ~ mca_auto + cll_prs + age + age_squ+ sex_at_birth+ smoking, 
               data = mCA_clean_A_U, family = binomial)
summary(model.y)

med.out <- mediate(model.m, model.y, treat = "cll_prs", mediator = "mca_auto", robustSE = TRUE, sims = 100)
summary(med.out)

#High risk 
table(mCA_clean_A_U$mca_highrisk)
mCA_clean_A_U <- filter(mCA_clean_A_U,mCA_clean_A_U$mca_highrisk==0 | mCA_clean_A_U$mca_highrisk==1)
mCA_clean_A_U <- mCA_clean_A_U %>%
  mutate(mca_highrisk = ifelse(mca_highrisk == 2, NA, mca_highrisk)) 
table(mCA_clean_A_U$mca_highrisk)

mCA_clean_A_U$mca_highrisk <- ifelse(mCA_clean_A_U$mca_highrisk ==1, 0, mCA_clean_A_U$mca_highrisk)
mCA_clean_A_U$mca_highrisk <- ifelse(mCA_clean_A_U$mca_highrisk ==2, 1, mCA_clean_A_U$mca_highrisk)
table(mCA_clean_A_U$mca_highrisk)

model.m <- glm(mca_highrisk ~ cll_prs + age + age_squ+ sex_at_birth+ smoking, 
               data = mCA_clean_A_U, family = binomial)
summary(model.m)

model.y <- glm(cll ~ mca_highrisk + cll_prs + age + age_squ+ sex_at_birth+ smoking, 
               data = mCA_clean_A_U, family = binomial)
summary(model.y)

med.out <- mediate(model.m, model.y, treat = "cll_prs", mediator = "mca_highrisk", robustSE = TRUE, sims = 100)
summary(med.out)
