##-----Mediation Analysis------##
#Load packages
library(tidyr)
library(dplyr)
library(tidyverse)
library(data.table)
library(lubridate)
library(broom)
library(survival)
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
