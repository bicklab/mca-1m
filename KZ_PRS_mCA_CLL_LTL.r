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

library(tidyr)
library(dplyr)
library(tidyverse)
library(data.table)
library(lubridate)
library(broom)
library(survival)
library(mgcv)
library(ggplot2)

x <- load("data_all_BUA_CLL_LTL.Rdata")
x
dim(data_all_BUA)

table(data_all_BUA$cohort)

library(data.table)

setDT(data_all_BUA)

data_all_BUA[, .(
  min_ltl = min(ltl_prs, na.rm = TRUE),
  q1_ltl  = quantile(ltl_prs, 0.25, na.rm = TRUE),
  median_ltl = median(ltl_prs, na.rm = TRUE),
  mean_ltl = mean(ltl_prs, na.rm = TRUE),
  q3_ltl  = quantile(ltl_prs, 0.75, na.rm = TRUE),
  max_ltl = max(ltl_prs, na.rm = TRUE),
  
  min_cll = min(cll_prs, na.rm = TRUE),
  q1_cll  = quantile(cll_prs, 0.25, na.rm = TRUE),
  median_cll = median(cll_prs, na.rm = TRUE),
  mean_cll = mean(cll_prs, na.rm = TRUE),
  q3_cll  = quantile(cll_prs, 0.75, na.rm = TRUE),
  max_cll = max(cll_prs, na.rm = TRUE)
), by = cohort]

system("gsutil cp gs://bicklab-main-storage/Users/Kun_Zhao/all_3_cohorts_-1_to_1_normalized.tsv .")

data_all_BUA <- fread("all_3_cohorts_-1_to_1_normalized.tsv")
head(data_all_BUA)

library(data.table)

setDT(data_all_BUA)

data_all_BUA[, .(
  min_ltl = min(ltl_prs, na.rm = TRUE),
  q1_ltl  = quantile(ltl_prs, 0.25, na.rm = TRUE),
  median_ltl = median(ltl_prs, na.rm = TRUE),
  mean_ltl = mean(ltl_prs, na.rm = TRUE),
  q3_ltl  = quantile(ltl_prs, 0.75, na.rm = TRUE),
  max_ltl = max(ltl_prs, na.rm = TRUE),
  
  min_cll = min(cll_prs, na.rm = TRUE),
  q1_cll  = quantile(cll_prs, 0.25, na.rm = TRUE),
  median_cll = median(cll_prs, na.rm = TRUE),
  mean_cll = mean(cll_prs, na.rm = TRUE),
  q3_cll  = quantile(cll_prs, 0.75, na.rm = TRUE),
  max_cll = max(cll_prs, na.rm = TRUE)
), by = cohort]

head(data_all_BUA)

table(data_all_BUA$cohort)

data_all_UA <- filter(data_all_BUA, data_all_BUA$cohort != "bioVU")

ls(data_all_UA)

summary(data_all_UA$cll_prs)

data_all_pheno_cox7 <- readRDS("data_all_pheno_mca_biovu.rds")

table(data_all_pheno_cox7$mca_auto)
table(data_all_pheno_cox7$auto_mca_highrisk)

smoking <- fread("agd_ever_smoker_df.tsv")
data_all_pheno_cox7 <- merge(data_all_pheno_cox7, smoking, by = "person_id", all.x = T)
table(data_all_pheno_cox7$ever_smoker, useNA = "always")

data_all_pheno_cox7 <- data_all_pheno_cox7 %>%
  mutate(
    ever_smoker = ifelse(is.na(ever_smoker), "FALSE", ever_smoker)
  )

table(data_all_pheno_cox7$ever_smoker, useNA = "always")

table(data_all_pheno_cox7$unsupervised_ancestry_cluster_relabel)

data_biovu <- data_all_pheno_cox7[,c("GRID", "age", "age2", "mca_auto", "auto_mca_highrisk", "gender",
                                    "ever_smoker","CA_121.21","unsupervised_ancestry_cluster_relabel")]

data_biovu_eur <- filter(data_biovu, data_biovu$unsupervised_ancestry_cluster_relabel == "k5_(EUR)")

dim(data_biovu_eur)

cll_prs <- fread("concat_25_cll_prs_result.tsv")

summary(cll_prs$CLL_PRS_z)
summary(cll_prs$PRS)

names(cll_prs)[2] <- "cll_prs"

ltl_prs <- fread("concat_8_telo_prs_result.tsv")
names(ltl_prs)[2] <- "ltl_prs"

data_biovu_eur <- merge(data_biovu_eur, cll_prs, by.x = "GRID", by.y = "#IID", all = F)

data_biovu_eur <- merge(data_biovu_eur, ltl_prs, by.x = "GRID", by.y = "#IID", all = F)

dim(data_biovu_eur)

data_all_UA <- data_all_UA %>%
  select(-cf_cat)

data_all_UA <- data_all_UA %>%
  select(-cll_262_coh, -ltl_262_coh)

colnames(data_all_UA)
colnames(data_biovu_eur)

setDT(data_biovu_eur)

data_biovu_eur[, `:=`(
  person_id = GRID,
  age_squ = age2,
  sex_at_birth = gender,
  smoking = ever_smoker,
  mca_highrisk = auto_mca_highrisk,
  cll = CA_121.21,
  ltl_prs = telo_PRS_z,
  cll_prs = CLL_PRS_z,
  cohort = "BioVU"
)]

data_biovu_eur <- data_biovu_eur[, .(
  person_id,
  age,
  age_squ,
  sex_at_birth,
  smoking,
  mca_auto,
  mca_highrisk,
  cll,
  ltl_prs,
  cll_prs,
  cohort
)]

data_biovu_eur

colnames(data_all_UA)
colnames(data_biovu_eur)

data_all_BUA_eur <- rbind(data_all_UA, data_biovu_eur)

dim(data_all_BUA_eur)

data_all_BUA_eur$ltl_262_coh <- 2
data_all_BUA_eur$ltl_262_coh[data_all_BUA_eur$ltl_prs >= quantile(data_all_BUA_eur$ltl_prs,0.8)] <- 3
data_all_BUA_eur$ltl_262_coh[data_all_BUA_eur$ltl_prs < quantile(data_all_BUA_eur$ltl_prs,0.2)] <- 1
data_all_BUA_eur$ltl_262_coh<-as.factor(data_all_BUA_eur$ltl_262_coh)
table(data_all_BUA_eur$ltl_262_coh)

data_all_BUA_eur$cll_262_coh <- 2
data_all_BUA_eur$cll_262_coh[data_all_BUA_eur$cll_prs >= quantile(data_all_BUA_eur$cll_prs,0.8)] <- 3
data_all_BUA_eur$cll_262_coh[data_all_BUA_eur$cll_prs < quantile(data_all_BUA_eur$cll_prs,0.2)] <- 1
data_all_BUA_eur$cll_262_coh<-as.factor(data_all_BUA_eur$cll_262_coh)
table(data_all_BUA_eur$cll_262_coh)

data_all_BUA_eur <- data_all_BUA_eur %>%
  group_by(cohort) %>%
  mutate(
    cll_262 = case_when(
      cll_prs >= quantile(cll_prs, 0.8, na.rm = TRUE) ~ 3,
      cll_prs <  quantile(cll_prs, 0.2, na.rm = TRUE) ~ 1,
      TRUE ~ 2
    )
  ) %>%
  ungroup()

data_all_BUA_eur$cll_262 <- as.factor(data_all_BUA_eur$cll_262)

data_all_BUA_eur <- data_all_BUA_eur %>%
  group_by(cohort) %>%
  mutate(
    ltl_262 = case_when(
      ltl_prs >= quantile(ltl_prs, 0.8, na.rm = TRUE) ~ 3,
      ltl_prs <  quantile(ltl_prs, 0.2, na.rm = TRUE) ~ 1,
      TRUE ~ 2
    )
  ) %>%
  ungroup()

data_all_BUA_eur$ltl_262 <- as.factor(data_all_BUA_eur$ltl_262)

data_all_BUA_eur <- data_all_BUA_eur %>%
  group_by(cohort) %>%
  mutate(
    cll_333 = case_when(
      cll_prs >= quantile(cll_prs, 2/3, na.rm = TRUE) ~ 3,
      cll_prs <  quantile(cll_prs, 1/3, na.rm = TRUE) ~ 1,
      TRUE ~ 2
    )
  ) %>%
  ungroup()

data_all_BUA_eur$cll_333 <- as.factor(data_all_BUA_eur$cll_333)

data_all_BUA_eur <- data_all_BUA_eur %>%
  group_by(cohort) %>%
  mutate(
    ltl_333 = case_when(
      ltl_prs >= quantile(ltl_prs, 2/3, na.rm = TRUE) ~ 3,
      ltl_prs <  quantile(ltl_prs, 1/3, na.rm = TRUE) ~ 1,
      TRUE ~ 2
    )
  ) %>%
  ungroup()

data_all_BUA_eur$ltl_333 <- as.factor(data_all_BUA_eur$ltl_333)

table(data_all_BUA_eur$cll_262, data_all_BUA_eur$cll_262_coh)
table(data_all_BUA_eur$cll_262)

table(data_all_BUA_eur$ltl_262)

saveRDS(data_all_BUA_eur, file = "data_all_BUA_CLL_LTL_eur_2026.RDS")

data_all_BUA_eur <- readRDS("data_all_BUA_CLL_LTL_eur_2026.RDS")

data_all_B_eur <- filter(data_all_BUA_eur, data_all_BUA_eur$cohort == "BioVU")
data_all_U_eur <- filter(data_all_BUA_eur, data_all_BUA_eur$cohort == "AoUv8")
data_all_A_eur <- filter(data_all_BUA_eur, data_all_BUA_eur$cohort == "BioVU")

gam_model <- gam(mca_auto ~ s(age, bs = "cs") + cll_262, data = data_all_B_eur, family = binomial)

age_points <- expand.grid(
  age = c(40, 50, 60, 70, 80), 
  cll_262 = factor(c("1", "2", "3"), levels = c("1", "2", "3")) 
)

predictions <- predict(gam_model, newdata = age_points, type = "response", se.fit = TRUE)
age_points$fit <- predictions$fit  
age_points$se <- predictions$se.fit 

age_points$lower <- age_points$fit - 1.96 * age_points$se
age_points$upper <- age_points$fit + 1.96 * age_points$se

print(age_points)

p <- ggplot(age_points, aes(x = age, y = fit, color = cll_262)) +
  geom_pointrange(aes(ymin = lower, ymax = upper), size = 0.8, shape = 15, fill = "white") + 
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.6, size = 0.6) +
  geom_line(size = 1) +
  scale_color_manual(
    values = c("1" = "#009E73", "2" = "#56B4E9", "3" = "#E69F00"),
    labels = c("1" = "Low", "2" = "Intermediate", "3" = "High")
  ) +
  labs(
    title = "Age vs mCA Status by CLL-PRS Groups",
    x = "Age",
    y = "mCA Status Probability",
    color = "CLL-PRS Group"
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

gam_model <- gam(mca_auto ~ s(age, bs = "cs") + cll_262, data = data_all_BUA_eur, family = binomial)

age_points <- expand.grid(
  age = c(40, 50, 60, 70, 80), 
  cll_262 = factor(c("1", "2", "3"), levels = c("1", "2", "3")) 
)

predictions <- predict(gam_model, newdata = age_points, type = "response", se.fit = TRUE)
age_points$fit <- predictions$fit  
age_points$se <- predictions$se.fit 

age_points$lower <- age_points$fit - 1.96 * age_points$se
age_points$upper <- age_points$fit + 1.96 * age_points$se

print(age_points)

p <- ggplot(age_points, aes(x = age, y = fit, color = cll_262)) +
  geom_pointrange(aes(ymin = lower, ymax = upper), size = 0.8, shape = 15, fill = "white") + 
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.6, size = 0.6) +
  geom_line(size = 1) +
  scale_color_manual(
    values = c("1" = "#009E73", "2" = "#56B4E9", "3" = "#E69F00"),
    labels = c("1" = "Low", "2" = "Intermediate", "3" = "High")
  ) +
  labs(
    title = "Age vs mCA Status by CLL-PRS Groups",
    x = "Age",
    y = "mCA Status Probability",
    color = "CLL-PRS Group"
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
save_plot("plot_CLL_PRS_mCA_BUAT.svg", fig = p, width=20, height=15)

gam_model <- gam(mca_auto ~ s(age, bs = "cs") + ltl_262, data = data_all_BUA_eur, family = binomial)

age_points <- expand.grid(
  age = c(40, 50, 60, 70, 80), 
  ltl_262 = factor(c("1", "2", "3"), levels = c("1", "2", "3")) 
)

predictions <- predict(gam_model, newdata = age_points, type = "response", se.fit = TRUE)
age_points$fit <- predictions$fit  
age_points$se <- predictions$se.fit 

age_points$lower <- age_points$fit - 1.96 * age_points$se
age_points$upper <- age_points$fit + 1.96 * age_points$se

print(age_points)

p<- ggplot(age_points, aes(x = age, y = fit, color = ltl_262)) +
  geom_pointrange(aes(ymin = lower, ymax = upper), size = 0.8, shape = 15, fill = "white") + 
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.6, size = 0.6) +
  geom_line(size = 1) +
  scale_color_manual(
    values = c("1" = "#009E73", "2" = "#56B4E9", "3" = "#E69F00"),
    labels = c("1" = "Low", "2" = "Intermediate", "3" = "High")
  ) +
  labs(
    title = "Age vs mCA Status by LTL-PRS Groups",
    x = "Age",
    y = "mCA Status Probability",
    color = "LTL-PRS Group"
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
save_plot("plot_LTL_PRS_mCA_BUAT.svg", fig = p, width=20, height=15)

gam_model <- gam(mca_auto ~ s(age, bs = "cs") + ltl_262, data = data_all_B_eur, family = binomial)

age_points <- expand.grid(
  age = c(40, 50, 60, 70, 80), 
  ltl_262 = factor(c("1", "2", "3"), levels = c("1", "2", "3")) 
)

predictions <- predict(gam_model, newdata = age_points, type = "response", se.fit = TRUE)
age_points$fit <- predictions$fit  
age_points$se <- predictions$se.fit 

age_points$lower <- age_points$fit - 1.96 * age_points$se
age_points$upper <- age_points$fit + 1.96 * age_points$se

print(age_points)

p<- ggplot(age_points, aes(x = age, y = fit, color = ltl_262)) +
  geom_pointrange(aes(ymin = lower, ymax = upper), size = 0.8, shape = 15, fill = "white") + 
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.6, size = 0.6) +
  geom_line(size = 1) +
  scale_color_manual(
    values = c("1" = "#009E73", "2" = "#56B4E9", "3" = "#E69F00"),
    labels = c("1" = "Low", "2" = "Intermediate", "3" = "High")
  ) +
  labs(
    title = "Age vs mCA Status by LTL-PRS Groups",
    x = "Age",
    y = "mCA Status Probability",
    color = "LTL-PRS Group"
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



#High Risk
table(data_all_BUA_eur$mca_highrisk,data_all_BUA_eur$mca_auto, useNA = "always")

mCA_clean_BUA_high <- filter(data_all_BUA_eur,data_all_BUA_eur$mca_highrisk==0 | data_all_BUA_eur$mca_highrisk==1)
mCA_clean_BUA_high <- mCA_clean_BUA_high %>%
  mutate(mca_highrisk = ifelse(mca_highrisk == 2, NA, mca_highrisk)) 
table(mCA_clean_BUA_high$mca_highrisk)

mCA_clean_BUA_high$mca_highrisk <- ifelse(mCA_clean_BUA_high$mca_highrisk ==1, 0, mCA_clean_BUA_high$mca_highrisk)
mCA_clean_BUA_high$mca_highrisk <- ifelse(mCA_clean_BUA_high$mca_highrisk ==3, 1, mCA_clean_BUA_high$mca_highrisk)
table(mCA_clean_BUA_high$mca_highrisk)

gam_model <- gam(mca_highrisk ~ s(age, bs = "cs") + cll_262, data = mCA_clean_BUA_high, family = binomial)

age_points <- expand.grid(
  age = c(40, 50, 60, 70, 80), 
  cll_262 = factor(c("1", "2", "3"), levels = c("1", "2", "3")) 
)
predictions <- predict(gam_model, newdata = age_points, type = "response", se.fit = TRUE)
age_points$fit <- predictions$fit  
age_points$se <- predictions$se.fit 

age_points$lower <- age_points$fit - 1.96 * age_points$se
age_points$upper <- age_points$fit + 1.96 * age_points$se

print(age_points)

p <- ggplot(age_points, aes(x = age, y = fit, color = cll_262)) +
  geom_pointrange(aes(ymin = lower, ymax = upper), size = 0.8, shape = 15, fill = "white") + 
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 1, size = 1) +
  geom_line(size = 1) +
  scale_color_manual(
    values = c("1" = "#009E73", "2" = "#56B4E9", "3" = "#E69F00"),
    labels = c("1" = "Low", "2" = "Intermediate", "3" = "High")
  ) +
  labs(
    title = "Age vs high-risk mCA Status by CLL-PRS Groups",
    x = "Age",
    y = "high-risk mCA Probability",
    color = "CLL-PRS Group"
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

save_plot("plot_CLL_PRS_highriskmCA_BUA.svg", fig = p, width=20, height=15)

gam_model <- gam(mca_highrisk ~ s(age, bs = "cs") + ltl_262, data = mCA_clean_BUA_high, family = binomial)

age_points <- expand.grid(
  age = c(40, 50, 60, 70, 80), 
  ltl_262 = factor(c("1", "2", "3"), levels = c("1", "2", "3")) 
)
predictions <- predict(gam_model, newdata = age_points, type = "response", se.fit = TRUE)
age_points$fit <- predictions$fit  
age_points$se <- predictions$se.fit 

age_points$lower <- age_points$fit - 1.96 * age_points$se
age_points$upper <- age_points$fit + 1.96 * age_points$se

print(age_points)

p <- ggplot(age_points, aes(x = age, y = fit, color = ltl_262)) +
  geom_pointrange(aes(ymin = lower, ymax = upper), size = 0.8, shape = 15, fill = "white") + 
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 1, size = 1) +
  geom_line(size = 1) +
  scale_color_manual(
    values = c("1" = "#009E73", "2" = "#56B4E9", "3" = "#E69F00"),
    labels = c("1" = "Low", "2" = "Intermediate", "3" = "High")
  ) +
  labs(
    title = "Age vs high-risk mCA Status by LTL-PRS Groups",
    x = "Age",
    y = "high-risk mCA Probability",
    color = "LTL-PRS Group"
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

save_plot("plot_LTL_PRS_highriskmCA_BUA.svg", fig = p, width=20, height=15)

dim(data_all_BUA_eur)

table(data_all_BUA_eur$cll_262_coh)

library(data.table)

setDT(data_all_BUA_eur)

data_all_BUA_eur[, mca_cll := 
  fifelse(mca_auto == 0 & cll_262 %in% c(1, 2), 0,
  fifelse(mca_auto == 0 & cll_262 == 3, 1,
  fifelse(mca_auto == 1 & cll_262 %in% c(1, 2), 2,
  fifelse(mca_auto == 1 & cll_262 == 3, 3, NA_integer_))))
]

table(data_all_BUA_eur$mca_cll)
table(data_all_BUA_eur$cohort)

head(data_all_B_eur)

data_all_BUA_eur[, mca_cll := factor(mca_cll, levels = 0:3)]

#biovu
data_all_B_eur <- filter(data_all_BUA_eur, data_all_BUA_eur$cohort == "BioVU")

data_all_pheno_cox7 <- readRDS("data_all_pheno_mca_biovu.rds")

#age 40
data_all_pheno_cox7_40 <- filter(data_all_pheno_cox7, data_all_pheno_cox7$age >= 40)

data_all_pheno_cox7_40_eur <- merge(data_all_B_eur, data_all_pheno_cox7_40, by.x = "person_id", by.y = "GRID", all = F)

dim(data_all_pheno_cox7_40_eur)

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

    formula_str <- paste("surv_obj ~", paste(c(covariates, var), collapse = " + "))
    cox_formula <- as.formula(formula_str)

    model <- coxph(cox_formula, data = df)

    tidy_res <- broom::tidy(model, exponentiate = TRUE, conf.int = TRUE)

    # 👉 抓所有该变量的level
    var_rows <- tidy_res %>%
      filter(grepl(paste0("^", var), term))

    if (nrow(var_rows) == 0) next

    # 👉 每个 level 单独算 N
    level_results <- lapply(1:nrow(var_rows), function(i) {

      term_name <- var_rows$term[i]

      # 提取 level（例如 mca_cll3 → 3）
      level_value <- gsub(var, "", term_name)

      df_sub <- df %>% filter(!!var_sym == as.numeric(level_value))

      N <- nrow(df_sub)
      N_event <- sum(df_sub$status == 1, na.rm = TRUE)

      event_rate <- if (total_event > 0) N_event / total_event else NA_real_

      data.frame(
        variable = var,
        level = level_value,
        N = N,
        N_event = N_event,
        event_rate = event_rate,
        HR = var_rows$estimate[i],
        CI_lower = var_rows$conf.low[i],
        CI_upper = var_rows$conf.high[i],
        p_value = var_rows$p.value[i]
      )
    })

    results_list[[var]] <- bind_rows(level_results)
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

table(data_all_pheno_cox7_40_eur$mca_cll)

cox_result_time <- run_cox_model_each_var(
  data = data_all_pheno_cox7_40_eur,
  event_col = "CA_121.21",
  event_time_col = "time_CA_121.21",
  baseline_time_col = "SHIFTED_SAMPLE_DATE",
  last_record_time_col = "observation_period_end_date",
  covariates = c("age.x","gender","PC1_SUM","PC2_SUM","PC3_SUM","PC4_SUM","PC5_SUM"),
  variables_of_interest = c("mca_cll")
)

cox_result_time

run_mediation_batch <- function(df,
                                protein_list,
                                chip_var,
                                outcome_var,
                                covariates){

  results <- list()

  # -------- Total effect model (c path) --------
  formula_total <- as.formula(
    paste(outcome_var, "~", chip_var, "+", paste(covariates, collapse = "+"))
  )
  model_total <- lm(formula_total, data = df)

  beta_total <- coef(summary(model_total))[chip_var, "Estimate"]
  p_total    <- coef(summary(model_total))[chip_var, "Pr(>|t|)"]

  for (prot in protein_list){

    # -------- a path --------
    formula_a <- as.formula(
      paste(prot, "~", chip_var, "+", paste(covariates, collapse = "+"))
    )
    model_a <- lm(formula_a, data = df)

    if (!(chip_var %in% rownames(coef(summary(model_a))))) next

    beta_a <- coef(summary(model_a))[chip_var, "Estimate"]
    se_a   <- coef(summary(model_a))[chip_var, "Std. Error"]

    # -------- b path + direct effect --------
    formula_b <- as.formula(
      paste(outcome_var, "~", prot, "+", chip_var, "+", paste(covariates, collapse = "+"))
    )
    model_b <- lm(formula_b, data = df)

    if (!(prot %in% rownames(coef(summary(model_b))))) next

    beta_b <- coef(summary(model_b))[prot, "Estimate"]
    se_b   <- coef(summary(model_b))[prot, "Std. Error"]

    beta_direct <- coef(summary(model_b))[chip_var, "Estimate"]
    p_direct    <- coef(summary(model_b))[chip_var, "Pr(>|t|)"]

    # -------- indirect effect --------
    indirect <- beta_a * beta_b
    se_indirect <- sqrt((beta_b^2 * se_a^2) + (beta_a^2 * se_b^2))
    z <- indirect / se_indirect
    p_med <- 2 * (1 - pnorm(abs(z)))

    # -------- proportion mediated --------
    prop_mediated <- indirect / beta_total
    percent_mediated <- prop_mediated * 100

    results[[prot]] <- data.frame(
      protein = prot,

      beta_total = beta_total,
      p_total    = p_total,

      beta_direct = beta_direct,
      p_direct    = p_direct,

      beta_a = beta_a,
      beta_b = beta_b,

      indirect = indirect,
      se_indirect = se_indirect,
      p_med = p_med,

      prop_mediated = prop_mediated,
      percent_mediated = percent_mediated
    )
  }

  final <- dplyr::bind_rows(results)
  final$FDR <- p.adjust(final$p_med, method = "fdr")

  return(final)
}

data_all_BUA <- fread("all_3_cohorts_z_normalized.tsv")
dim(data_all_BUA)

colnames(data_all_BUA)
dim(data_all_BUA)
table(data_all_BUA$cohort)

data_all_BUA_4090 <- filter(data_all_BUA, data_all_BUA$age >=40 & data_all_BUA$age <=90)

data_all_BUA_4090 <- data_all_BUA_4090 %>%
  mutate(
    smoking = case_when(
      smoking %in% c(TRUE, "TRUE", "Yes") ~ "Yes",
      smoking %in% c(FALSE, "FALSE", "No") ~ "No",
      TRUE ~ "Unknown"
    )
  )

data_all_BUA_4090 <- data_all_BUA_4090 %>%
  mutate(
    sex_at_birth = case_when(
      tolower(sex_at_birth) == "female" ~ "Female",
      tolower(sex_at_birth) == "male" ~ "Male",
      TRUE ~ NA_character_
    )
  )

data_all_B_4090 <- filter(data_all_BUA_4090, data_all_BUA_4090$cohort == "BioVU")
data_all_U_4090 <- filter(data_all_BUA_4090, data_all_BUA_4090$cohort == "UKB")
data_all_A_4090 <- filter(data_all_BUA_4090, data_all_BUA_4090$cohort == "AoUv8")

mediation_results <- run_mediation_batch(
  df = data_all_B_4090,
  protein_list = c("mca_auto", "mca_highrisk"),
  chip_var = "cll_prs_z",
  outcome_var = "cll",
  covariates = c("age","age_squ","sex_at_birth", "smoking")
)
mediation_results

mediation_results

library(mediation)

model.m <- glm(mca_auto ~ cll_prs_z + age + age_squ+ sex_at_birth+ smoking, 
               data = data_all_BUA_4090, family = binomial)
summary(model.m)

model.y <- glm(cll ~ mca_auto + cll_prs_z + age + age_squ+ sex_at_birth+ smoking, 
               data = data_all_BUA_4090, family = binomial)
summary(model.y)

med.out <- mediate(model.m, model.y, treat = "cll_prs_z", mediator = "mca_auto", robustSE = TRUE, sims = 100)
summary(med.out)


