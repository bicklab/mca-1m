#-PheWAS meta-analysis-------
#-2.1.AoU
system("gsutil cp gs://fc-secure-cb192ac6-30ba-46b9-92ee-896a6e36c63e/kun/mCH_PheWAS_results_cox_phecodeX.csv .", intern = TRUE)
system("gsutil cp gs://fc-secure-cb192ac6-30ba-46b9-92ee-896a6e36c63e/kun/lCH_PheWAS_results_cox_phecodeX.csv .", intern = TRUE)

PheWAS_mCA_Aou <- read.csv("mCA_PheWAS_results_cox_phecodeX.csv")
PheWAS_lCH_Aou <- read.csv("lCH_PheWAS_results_cox_phecodeX.csv")
PheWAS_mCH_Aou <- read.csv("mCH_PheWAS_results_cox_phecodeX.csv")
#-2.2.UKB
PheWAS_mCA_UKB <- read.csv("PheWAS_mCA_ukb.csv")
PheWAS_lCH_UKB <- read.csv("PheWAS_lCH_ukb.csv")
PheWAS_mCH_UKB <- read.csv("PheWAS_mCH_ukb.csv")
#-2.3.BioVU
PheWAS_mCA_BioVU <- read.csv("PheWAS_mCA_cox_BioVU.csv")
PheWAS_mCA_BioVU <- filter(PheWAS_mCA_BioVU, is.na(PheWAS_mCA_BioVU$log_p_value) ==F)
#--Meta-----
install.packages("metafor")
library(metafor)
#-mCH
combined_data <- merge(PheWAS_mCH_UKB, PheWAS_mCH_Aou, by = "disease",all=T, suffixes = c("_c1", "_c2"))
combined_data$total_samples_c1 <- ifelse(is.na(combined_data$total_samples_c1),0,combined_data$total_samples_c1)
combined_data$total_samples_c2 <- ifelse(is.na(combined_data$total_samples_c2),0,combined_data$total_samples_c2)
combined_data$total_sample <- combined_data$total_samples_c1 + combined_data$total_samples_c2
combined_data_mCH <- combined_data
meta_results <- list()
for (i in 1:nrow(combined_data)) {
  betas <- c(combined_data$estimate_c1[i], combined_data$estimate_c2[i])
  ses <- c(combined_data$std.error_c1[i], combined_data$std.error_c2[i])
  ns <- c(combined_data$total_samples_c1[i], combined_data$total_samples_c2[i])
  ps <- c(combined_data$p.value_c1[i], combined_data$p.value_c2[i])
  non_na_indices <- which(!is.na(betas) & !is.na(ses))
  
  if (length(non_na_indices) == 1) {
    meta_beta <- betas[non_na_indices]
    meta_se <- ses[non_na_indices]
    meta_p <- ps[non_na_indices]
    
  } else if (length(non_na_indices) > 1) {
    meta_res <- rma(yi = betas[non_na_indices], sei = ses[non_na_indices], method = "FE")
    meta_beta <- meta_res$beta
    meta_se <- meta_res$se
    meta_p <- meta_res$pval
    
  } else {
    next
  }
  
  meta_results[[i]] <- data.frame(
    phenotype = combined_data$disease[i],
    meta_beta = meta_beta,
    meta_se = meta_se,
    meta_p = meta_p,
    meta_logp = -log10(meta_p)
  )
}
meta_results_mCH <- do.call(rbind, meta_results)
write.csv(meta_results_mCH, "meta_results_mCH.csv", row.names = FALSE)
meta_results_mCH <- read.csv("meta_results_mCH.csv")

combined_data$total_sample <- combined_data$total_samples_c1 + combined_data$total_samples_c2
#-lCH
combined_data <- merge(PheWAS_lCH_UKB, PheWAS_lCH_Aou, by = "disease",all=T, suffixes = c("_c1", "_c2"))
combined_data$total_samples_c1 <- ifelse(is.na(combined_data$total_samples_c1),0,combined_data$total_samples_c1)
combined_data$total_samples_c2 <- ifelse(is.na(combined_data$total_samples_c2),0,combined_data$total_samples_c2)
combined_data$total_sample <- combined_data$total_samples_c1 + combined_data$total_samples_c2
combined_data_lCH <- combined_data
meta_results <- list()
for (i in 1001:1692) {
  betas <- c(combined_data$estimate_c1[i], combined_data$estimate_c2[i])
  ses <- c(combined_data$std.error_c1[i], combined_data$std.error_c2[i])
  ns <- c(combined_data$total_samples_c1[i], combined_data$total_samples_c2[i])
  ps <- c(combined_data$p.value_c1[i], combined_data$p.value_c2[i])
  non_na_indices <- which(!is.na(betas) & !is.na(ses))
  
  if (length(non_na_indices) == 1) {
    meta_beta <- betas[non_na_indices]
    meta_se <- ses[non_na_indices]
    meta_p <- ps[non_na_indices]  
    
  } else if (length(non_na_indices) > 1) {
    meta_res <- rma(yi = betas[non_na_indices], sei = ses[non_na_indices], method = "FE")
    meta_beta <- meta_res$beta
    meta_se <- meta_res$se
    meta_p <- meta_res$pval
    
  } else {
    next
  }
  
  meta_results[[i]] <- data.frame(
    phenotype = combined_data$disease[i],
    meta_beta = meta_beta,
    meta_se = meta_se,
    meta_p = meta_p,
    meta_logp = -log10(meta_p)
  )
}
meta_results_lCH <- do.call(rbind, meta_results)
write.csv(meta_results_lCH, "meta_results_lCH.csv", row.names = FALSE)
meta_results_lCH <- read.csv("meta_results_lCH.csv")
#-mCA
combined_data <- merge(PheWAS_mCA_UKB, PheWAS_mCA_Aou, by = "disease",all=T, suffixes = c("_c1", "_c2"))
combined_data <- merge(combined_data, PheWAS_mCA_BioVU, by = "disease",all=T)
meta_results <- list()
for (i in 1:nrow(combined_data)) {
  # Extrave beta, se, n and p
  betas <- c(combined_data$estimate_c1[i], combined_data$estimate_c2[i], combined_data$estimate[i])
  ses <- c(combined_data$std.error_c1[i], combined_data$std.error_c2[i], combined_data$std.error[i])
  ns <- c(combined_data$total_samples_c1[i], combined_data$total_samples_c2[i], combined_data$total_samples[i])
  ps <- c(combined_data$p.value_c1[i], combined_data$p.value_c2[i], combined_data$p.value[i])
  # check beta and se (non- NA)
  non_na_indices <- which(!is.na(betas) & !is.na(ses))
  
  if (length(non_na_indices) == 1) {
    # If only one cohort has the phecode , keep it
    meta_beta <- betas[non_na_indices]
    meta_se <- ses[non_na_indices]
    meta_p <- ps[non_na_indices] 
    
  } else if (length(non_na_indices) > 1) {
    # If 2 or 3 cohorts have the phecode，use metafor to do Meta
    meta_res <- rma(yi = betas[non_na_indices], sei = ses[non_na_indices], method = "FE")
    
    # extract Meta results
    meta_beta <- meta_res$beta
    meta_se <- meta_res$se
    meta_p <- meta_res$pval
    
  } else {
    # Skip if no results
    next
  }
  
  # Save
  meta_results[[i]] <- data.frame(
    phenotype = combined_data$disease[i],
    meta_beta = meta_beta,
    meta_se = meta_se,
    meta_p = meta_p,
    meta_logp = -log10(meta_p)
  )
}
meta_results_mCA <- do.call(rbind, meta_results)
write.csv(meta_results_mCA, "meta_results_mCA.csv", row.names = FALSE)

system("gsutil cp meta_results_mCA.csv gs://fc-secure-cb192ac6-30ba-46b9-92ee-896a6e36c63e/kun/")
system("gsutil cp meta_results_lCH.csv gs://fc-secure-cb192ac6-30ba-46b9-92ee-896a6e36c63e/kun/")
system("gsutil cp meta_results_mCH.csv gs://fc-secure-cb192ac6-30ba-46b9-92ee-896a6e36c63e/kun/")

#--MANHATTAN plot (meta-mCA) -----
phecodeX <- read.csv("phecodeX_ICD_CM_map_flat.csv")
phecodeX2 <- phecodeX[,c(4,5,7)]
phecodeX2 <- distinct(phecodeX2)
meta_results_mCA<- read.csv("meta_results_mCA.csv")
results2 <- merge(meta_results_mCA,phecodeX2,by.x = "phenotype", by.y = "phecode_string",all.x = T)
results2 <- filter(results2,is.na(results2$category) == F)
names(results2)[1]<-"outcome"
names(results2)[2]<-"estimate"
names(results2)[5]<-"log_p_value"

results2 <- results2 %>%
  ungroup() %>%  
  dplyr::mutate(outcome_index = row_number(),   
                shape = ifelse(estimate > 0, 24, 25))  
top_10_outcomes <- results2 %>%
  arrange(desc(log_p_value)) %>%
  slice_head(n = 20) %>%
  pull(outcome) 

results2 <- results2 %>%
  mutate(top_10_label = ifelse(outcome %in% top_10_outcomes, outcome, NA))  

threshold_pvalue_0_05 <- 4.3  

results2 <- results2 %>%
  arrange(category, outcome_index) %>%
  dplyr::mutate(outcome_index = row_number())  

results2_clean <- results2 %>%
  filter(!is.na(log_p_value), !is.na(outcome_index))

table(results2_clean$category)
results2_clean$category<- as.factor(results2_clean$category)

results2_clean$top_10_label <- ifelse(is.na(results2_clean$top_10_label), "", results2_clean$top_10_label)

results2_clean$top_10_label <- ifelse(results2_clean$phecode == "CA_125", "Other malignant neoplasms", results2_clean$top_10_label)

results2_clean$category <- factor(results2_clean$category, levels = unique(results2_clean$category))
class(results2_clean$outcome_index)

# Define breaks
split_data <- split(results2_clean, results2_clean$category)
mid_values <- lapply(split_data, function(df) mean(df$outcome_index))
breaks <- unlist(mid_values)

# Define labels
labels <- as.character(unique(results2_clean$category))

ggplot(results2_clean, aes(x = outcome_index, y = log_p_value)) +
  geom_point(aes(color = category, shape = as.factor(shape), fill = category), size = 3) +
  scale_shape_manual(values = c(24, 25)) +
  geom_hline(yintercept = threshold_pvalue_0_05, linetype = "dashed", color = "red") +
  geom_text(aes(label = top_10_label), size = 2, vjust = -1) +  
  theme_minimal() +
  labs(x = "Category", y = "-log10(p.value)", title = "Manhattan Plot of meta PheWAS of mCA") +
  scale_color_viridis_d(name = "Category", option = "plasma") +  
  scale_fill_viridis_d(name = "Category", option = "plasma") +  
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_continuous(breaks = breaks,
                     labels = labels)

#--MANHATTAN plot (meta-m_CH) -----
phecodeX <- read.csv("phecodeX_ICD_CM_map_flat.csv")
phecodeX2 <- phecodeX[,c(4,5,7)]
phecodeX2 <- distinct(phecodeX2)
meta_results_mCH<- read.csv("meta_results_mCH.csv")
results2 <- merge(meta_results_mCH,phecodeX2,by.x = "phenotype", by.y = "phecode_string",all.x = T)
results2 <- filter(results2,is.na(results2$category) == F)
names(results2)[1]<-"outcome"
names(results2)[2]<-"estimate"
names(results2)[5]<-"log_p_value"

results2 <- results2 %>%
  ungroup() %>% 
  dplyr::mutate(outcome_index = row_number(),  
                shape = ifelse(estimate > 0, 24, 25))  

top_10_outcomes <- results2 %>%
  arrange(desc(log_p_value)) %>%
  slice_head(n = 20) %>%
  pull(outcome) 

results2 <- results2 %>%
  mutate(top_10_label = ifelse(outcome %in% top_10_outcomes, outcome, NA))  


threshold_pvalue_0_05 <- 4.3 

results2 <- results2 %>%
  arrange(category, outcome_index) %>%
  dplyr::mutate(outcome_index = row_number()) 

results2_clean <- results2 %>%
  filter(!is.na(log_p_value), !is.na(outcome_index))

table(results2_clean$category)
results2_clean$category<- as.factor(results2_clean$category)

results2_clean$top_10_label <- ifelse(is.na(results2_clean$top_10_label), "", results2_clean$top_10_label)

results2_clean$top_10_label <- ifelse(results2_clean$phecode == "CA_125", "Other malignant neoplasms", results2_clean$top_10_label)


results2_clean$category <- factor(results2_clean$category, levels = unique(results2_clean$category))
class(results2_clean$outcome_index)

# Define breaks
split_data <- split(results2_clean, results2_clean$category)
mid_values <- lapply(split_data, function(df) mean(df$outcome_index))
breaks <- unlist(mid_values)

# Define labels
labels <- as.character(unique(results2_clean$category))

ggplot(results2_clean, aes(x = outcome_index, y = log_p_value)) +
  geom_point(aes(color = category, shape = as.factor(shape), fill = category), size = 3) +
  scale_shape_manual(values = c(24, 25)) +
  geom_hline(yintercept = threshold_pvalue_0_05, linetype = "dashed", color = "red") +
  geom_text(aes(label = top_10_label), size = 2, vjust = -1) +  
  theme_minimal() +
  labs(x = "Category", y = "-log10(p.value)", title = "Manhattan Plot of meta PheWAS of m_CH") +
  scale_color_viridis_d(name = "Category", option = "plasma") +  
  scale_fill_viridis_d(name = "Category", option = "plasma") +  
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_continuous(breaks = breaks,
                     labels = labels)