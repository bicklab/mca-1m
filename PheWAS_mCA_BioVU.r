library(survival)
library(tidyr)
library(dplyr)
library(haven)

library(data.table)
library(lubridate)
library(ggplot2)
library(broom)
library(scales)

system("gsutil cp gs://bicklab-main-storage/Users/Yash_Pershad/mega_phecodes.txt .", intern = TRUE)

bioVU_all <- fread('mega_phecodes.txt') #590,3810

head(bioVU_all)

system("gsutil cp gs://fc-secure-343d99f4-8f3b-46da-ac8c-0f9f75cc97b3/notebooks/biovu_mca_calls..txt .", intern = TRUE)

system("gsutil cp gs://fc-secure-343d99f4-8f3b-46da-ac8c-0f9f75cc97b3/notebooks/BioVU_mCA_calls.txt .", intern = TRUE)

data_biovu1 <- fread('biovu_mca_calls..txt') 
data_biovu2 <- fread('BioVU_mCA_calls.txt') 

dim(data_biovu1)
dim(data_biovu2)

head(data_biovu1)
mca_biovu_auto <- filter(data_biovu1, data_biovu1$chrom != "X")
mca_biovu_auto %>% pull(sample_id) %>% unique() %>% length()

15100+803+3068+10462

head(data_biovu2)

summary(data_biovu2$YAgeAtDNA)

ls(data_biovu2)

table_ancestry <- table(data_biovu2$Race)
prop.table(table_ancestry) * 100

data_biovu <- data_biovu1[,c(1,3:22,48)]
data_biovu$p_count<-NA
data_biovu$q_count<-NA
data_biovu$p_count<-ifelse(data_biovu$p_arm=="Y"|data_biovu$p_arm=="T"|data_biovu$p_arm=="C",1,0)
data_biovu$q_count<-ifelse(data_biovu$q_arm=="Y"|data_biovu$q_arm=="T"|data_biovu$q_arm=="C",1,0)
table(data_biovu$p_count,data_biovu$q_count)


data_biovu$arm<-NA
data_biovu$arm<-ifelse(data_biovu$p_count ==1 & data_biovu$q_count ==0 ,"p",data_biovu$arm)
data_biovu$arm<-ifelse(data_biovu$q_count ==1 & data_biovu$p_count ==0 ,"q",data_biovu$arm)
table(data_biovu$arm)
data_biovu$chr_type <- paste("chr",data_biovu$chrom,data_biovu$arm,sep = "")
data_biovu <- data_biovu %>% mutate(chr_type = str_replace(chr_type, "NA$", ""))
names(data_biovu)[23] <- "age"
names(data_biovu)[9] <- "n_sites"
names(data_biovu)[10] <- "n_hets"
names(data_biovu)[19] <- "baf_conc"
ls(data_biovu)
data_biovu$cohort<-"BioVU"

dim(data_biovu2)
ls(data_biovu2)

data_biovu2 <- data_biovu2[,c(1,10:20)]
data_biovu$mca_status <- 1
BioVU_mCA <- merge(data_biovu2,data_biovu,by.y = "sample_id",by.x = "ourSid",all.x = T)

dim(BioVU_mCA)

table(BioVU_mCA$computed_gender,useNA="always")
table(BioVU_mCA$mca_status,useNA="always")

system("gsutil cp gs://fc-secure-343d99f4-8f3b-46da-ac8c-0f9f75cc97b3/notebooks/phecodeX_ICD_CM_map_flat.csv .", intern = TRUE)

#Load Phecode "phecode_icd10.csv"
phecodex <- read.csv("phecodeX_ICD_CM_map_flat.csv")
phecodex <- phecodex[,c("phecode","phecode_string","category")]
phecodex <- distinct(phecodex)
bioVU_all2 <- merge(bioVU_all,phecodex, by="phecode",all.x = T)
table(bioVU_all2$category,useNA = "always")
test<-filter(bioVU_all2,is.na(bioVU_all2$category))
table(test$phecode)
bioVU_all2$phecode_string <- ifelse(is.na(bioVU_all2$phecode_string),bioVU_all2$phecode,bioVU_all2$phecode_string)
bioVU_all2 <- bioVU_all2 %>%
  mutate(value = 1)

#Generate new dataset by each phenotype and start time (for COX regression)
data_pheno_cox <-  bioVU_all2[,c("GRID","phecode_string","first_phecode_date")]
names(data_pheno_cox)[1]<-"person_id"
names(data_pheno_cox)[2]<-"Phenotype"
names(data_pheno_cox)[3]<-"condition_start_datetime"

data_pheno_cox <- distinct(data_pheno_cox) #5903587
dim(data_pheno_cox)

# keep the earliest condition_start_datetime
data_pheno_cox_filtered <- data_pheno_cox %>%
  group_by(person_id, Phenotype) %>% 
  slice_min(condition_start_datetime) %>% 
  ungroup()

# Make sure the 'person_id' + 'Phenotype' is unique
data_pheno_cox_filtered <- data_pheno_cox_filtered %>%
  mutate(unique_id = rowid(person_id, Phenotype))
data_pheno_cox_filtered <- data_pheno_cox_filtered %>%
  filter(!is.na(Phenotype) & Phenotype != "")

# Convert into a Wide-format  (N=82874)
data_wide <- data_pheno_cox_filtered %>%
  mutate(value = 1) %>% 
  pivot_wider(names_from = Phenotype, values_from = value, values_fill = list(value = 0), id_cols = person_id) 

data_wide #82974 × 3514

# Create time columns as "time_'Disease'"
data_time <- data_pheno_cox_filtered %>%
  pivot_wider(names_from = Phenotype, values_from = condition_start_datetime, 
              names_prefix = "time_", values_fill = list(condition_start_datetime = NA), id_cols = person_id)

data_time #82974 × 3514

# Conbime data
data_final <- data_wide %>%
  left_join(data_time, by = "person_id")

data_final #82974 × 7027

save(data_final, file = "data_final.Rdata")

system("gsutil cp data_final.Rdata gs://bicklab-main-storage/Users/Kun_Zhao")

a <- load("data_final.Rdata")
a

ls(BioVU_mCA)

dim(BioVU_mCA)

BioVU_mCA$mca_status <- ifelse(is.na(BioVU_mCA$mca_status),0,BioVU_mCA$mca_status)
table(BioVU_mCA$mca_status)

summary(BioVU_mCA$YAgeAtDNA)

BioVU_mCA$age2 <- BioVU_mCA$YAgeAtDNA * BioVU_mCA$YAgeAtDNA

data_all_pheno_cox <- merge(data_final,BioVU_mCA,by.x ="person_id", by.y = "ourSid",all = F)

data_all_pheno_cox #83036 × 7066

save(data_all_pheno_cox,file = "data_all_pheno_cox_BioVU.Rdata")

system("gsutil cp data_all_pheno_cox_BioVU.Rdata gs://bicklab-main-storage/Users/Kun_Zhao")

system("gsutil cp gs://bicklab-main-storage/Users/Kun_Zhao/data_all_pheno_cox_BioVU.Rdata .")

x <- load("data_all_pheno_cox_BioVU.Rdata")
x

data_all_pheno_cox2 <- filter(data_all_pheno_cox,data_all_pheno_cox$YAgeAtDNA >=40 & data_all_pheno_cox$YAgeAtDNA <=90)

dim(data_all_pheno_cox2) #51202 × 7064

data_all_pheno_cox2$lastRecordDate <- as.Date(data_all_pheno_cox2$lastRecordDate, format = "%m/%d/%y")
data_all_pheno_cox2$DNADate <- as.Date(data_all_pheno_cox2$DNADate, format = "%m/%d/%y")

data_all_pheno_cox2$survial_cll <- NA
data_all_pheno_cox2$survial_cll <- as.numeric(difftime(data_all_pheno_cox2$`time_Chronic lymphoid leukemia`, data_all_pheno_cox2$DNADate, units ="days"))/365.25
summary(data_all_pheno_cox2$survial_cll)

data_unique <- data_all_pheno_cox2 %>%
  distinct(person_id, .keep_all = TRUE)

dim(data_unique)

table(data_unique$mca_status)

data_unique$mca_auto <- NA
data_unique$mca_auto <- ifelse(data_unique$person_id %in% data_biovu_auto$sample_id, 1, 0)

table(data_unique$mca_auto)

table(data_unique$mca_auto, useNA = "always")

#-PheWAS mCA by COX regression-----
cox_results <- list()
disease_cols <- colnames(data_all_pheno_cox2)[2:3514] #2:3514

for (disease_col in disease_cols) {
  time_col <- paste0("time_", disease_col)
  
  # Print "time_col"
  print(time_col)
  
  # Check time_col
  if (any(time_col %in% colnames(data_all_pheno_cox2))) {
    
    data_surv <- data_all_pheno_cox2 %>%
      mutate(
        surv_time = case_when(
          !!sym(disease_col) == 1 ~ as.numeric(difftime(!!sym(time_col), DNADate, units = "days")) / 365.25,
          !!sym(disease_col) == 0 ~ as.numeric(difftime(lastRecordDate, DNADate, units = "days")) / 365.25
        )
      ) %>%
      filter(surv_time > 0)
    
    total_samples <- sum(data_surv[[disease_col]] == 1)
    mca_status_1_outcome <- sum(data_surv$mca_status == 1 & data_surv[[disease_col]] == 1)
    
    surv_obj <- Surv(time = data_surv$surv_time, event = data_surv[[disease_col]])
    
    cox_model <- coxph(surv_obj ~ mca_status + YAgeAtDNA + age2 + SmokingStatus + computed_gender + Ethnicity, data = data_surv)
    
    cox_result <- tidy(cox_model) %>%
      filter(term == "mca_status") %>%
      mutate(
        total_samples = total_samples,
        mca_status_1_outcome = mca_status_1_outcome,
        log_p_value = -log10(p.value)  # -log10(p.value)
      )
    
    # Save
    cox_results[[disease_col]] <- cox_result
  } else {
    warning(paste("The", time_col, "Not exist!"))
  }
}
results_df <- bind_rows(cox_results, .id = "disease")
write.csv(results_df, "PheWAS_mCA_cox_BioVU.csv", row.names = FALSE)
system("gsutil cp PheWAS_mCA_cox_BioVU.csv gs://bicklab-main-storage/Users/Kun_Zhao")

results_df

summary(results_df$log_p_value)

results_df2 <- filter(results_df,results_df$log_p_value> 20)
results_df2

write.csv(results_df, "PheWAS_mCA_cox_BioVU2.csv", row.names = FALSE)
system("gsutil cp PheWAS_mCA_cox_BioVU2.csv gs://bicklab-main-storage/Users/Kun_Zhao")

data_biovu_auto <- filter(data_biovu1, data_biovu1$chrom != "X")
dim(data_biovu_auto)

x <- load("data_pheno_BioVU_4090.Rdata")
x

data_all_pheno_cox3$mca_auto <- NA
data_all_pheno_cox3$mca_auto <- ifelse(data_all_pheno_cox3$person_id %in% data_biovu_auto$sample_id, 1, 0)
table(data_all_pheno_cox3$mca_auto, useNA = "always")

dim(data_all_pheno_cox3)
data_all_pheno_cox3 <- data_all_pheno_cox3 %>% distinct(person_id, .keep_all = TRUE)
dim(data_all_pheno_cox3)

data_all_pheno_cox3$lastRecordDate <- as.Date(data_all_pheno_cox3$lastRecordDate, format = "%m/%d/%y")
data_all_pheno_cox3$DNADate <- as.Date(data_all_pheno_cox3$DNADate, format = "%m/%d/%y")

#-PheWAS mCA by COX regression-----
cox_results <- list()
disease_cols <- colnames(data_all_pheno_cox3)[2:3514] #2:3514

for (i in seq_along(disease_cols)) {
  disease_col <- disease_cols[i]
  time_col <- paste0("time_", disease_col)

  cat("Processing:", time_col, "（Loop", i, "of", length(disease_cols), "）\n")
    
  if (!any(time_col %in% colnames(data_all_pheno_cox3))) {
    warning(paste("Time column", time_col, "does not exist! Skipping", disease_col))
    next
  }
  
  unique_values <- unique(data_all_pheno_cox3[[disease_col]])
  if (!all(unique_values %in% c(0, 1))) {
    warning(paste("Column", disease_col, "contains non-binary values. Skipping this column."))
    next
  }

  data_surv <- data_all_pheno_cox3 %>%
    mutate(
      surv_time = case_when(
        !!sym(disease_col) == 1 ~ as.numeric(difftime(!!sym(time_col), DNADate, units = "days")) / 365.25,
        !!sym(disease_col) == 0 ~ as.numeric(difftime(lastRecordDate, DNADate, units = "days")) / 365.25
      )
    ) %>%
    filter(surv_time > 0, !is.infinite(surv_time), !is.na(surv_time)) 

  total_samples <- sum(data_surv[[disease_col]] == 1)
  if (total_samples < 10) {  # Filter the sample size
    warning(paste("Not enough events for", disease_col, "- Skipping."))
    next
  }

  surv_obj <- Surv(time = data_surv$surv_time, event = data_surv[[disease_col]])

  cox_result <- tryCatch({
    cox_model <- coxph(surv_obj ~ mca_auto + YAgeAtDNA + age2 + SmokingStatus + computed_gender.y + ancestry, data = data_surv)

    tidy(cox_model) %>%
      filter(term == "mca_auto") %>%
      mutate(
        total_samples = total_samples,
        log_p_value = -log10(p.value)
      )
  }, error = function(e) {
    warning(paste("Cox model failed for", disease_col, "with error:", e$message))
    return(NULL)
  })

  if (!is.null(cox_result)) {
    cox_results[[disease_col]] <- cox_result
  }
}

results_df <- bind_rows(cox_results, .id = "disease")
results_df <- results_df %>% arrange(desc(log_p_value))
results_df

write.csv(results_df, "PheWAS_automCA_cox_BioVU_0519.csv", row.names = FALSE)
system("gsutil cp PheWAS_automCA_cox_BioVU_0519.csv gs://bicklab-main-storage/Users/Kun_Zhao")

table(data_all_pheno_cox3$ancestry)

data_eur_pheno <- filter(data_all_pheno_cox3, data_all_pheno_cox3$ancestry == "EUR")
dim(data_eur_pheno)

#-PheWAS mCA by COX regression-----
cox_results <- list()
disease_cols <- colnames(data_all_pheno_cox3)[2:3514] #2:3514

for (i in seq_along(disease_cols)) {
  disease_col <- disease_cols[i]
  time_col <- paste0("time_", disease_col)

  cat("Processing:", time_col, "（Loop", i, "of", length(disease_cols), "）\n")
    
  if (!any(time_col %in% colnames(data_eur_pheno))) {
    warning(paste("Time column", time_col, "does not exist! Skipping", disease_col))
    next
  }
  
  unique_values <- unique(data_eur_pheno[[disease_col]])
  if (!all(unique_values %in% c(0, 1))) {
    warning(paste("Column", disease_col, "contains non-binary values. Skipping this column."))
    next
  }

  data_surv <- data_eur_pheno %>%
    mutate(
      surv_time = case_when(
        !!sym(disease_col) == 1 ~ as.numeric(difftime(!!sym(time_col), DNADate, units = "days")) / 365.25,
        !!sym(disease_col) == 0 ~ as.numeric(difftime(lastRecordDate, DNADate, units = "days")) / 365.25
      )
    ) %>%
    filter(surv_time > 0, !is.infinite(surv_time), !is.na(surv_time)) 

  total_samples <- sum(data_surv[[disease_col]] == 1)
  if (total_samples < 10) {  # Filter the sample size
    warning(paste("Not enough events for", disease_col, "- Skipping."))
    next
  }

  surv_obj <- Surv(time = data_surv$surv_time, event = data_surv[[disease_col]])

  cox_result <- tryCatch({
    cox_model <- coxph(surv_obj ~ mca_auto + YAgeAtDNA + age2 + SmokingStatus + computed_gender.y, data = data_surv)

    tidy(cox_model) %>%
      filter(term == "mca_auto") %>%
      mutate(
        total_samples = total_samples,
        log_p_value = -log10(p.value)
      )
  }, error = function(e) {
    warning(paste("Cox model failed for", disease_col, "with error:", e$message))
    return(NULL)
  })

  if (!is.null(cox_result)) {
    cox_results[[disease_col]] <- cox_result
  }
}

results_df <- bind_rows(cox_results, .id = "disease")
results_df <- results_df %>% arrange(desc(log_p_value))
results_df

write.csv(results_df, "PheWAS_automCA_cox_BioVU_EUR.csv", row.names = FALSE)
system("gsutil cp PheWAS_automCA_cox_BioVU_EUR.csv gs://bicklab-main-storage/Users/Kun_Zhao")

data_afr_pheno <- filter(data_all_pheno_cox3, data_all_pheno_cox3$ancestry == "AFR")
dim(data_afr_pheno)

#-PheWAS mCA by COX regression-----
cox_results <- list()
disease_cols <- colnames(data_all_pheno_cox3)[2:3514] #2:3514

for (i in seq_along(disease_cols)) {
  disease_col <- disease_cols[i]
  time_col <- paste0("time_", disease_col)

  cat("Processing:", time_col, "（Loop", i, "of", length(disease_cols), "）\n")
    
  if (!any(time_col %in% colnames(data_afr_pheno))) {
    warning(paste("Time column", time_col, "does not exist! Skipping", disease_col))
    next
  }
  
  unique_values <- unique(data_afr_pheno[[disease_col]])
  if (!all(unique_values %in% c(0, 1))) {
    warning(paste("Column", disease_col, "contains non-binary values. Skipping this column."))
    next
  }

  data_surv <- data_afr_pheno %>%
    mutate(
      surv_time = case_when(
        !!sym(disease_col) == 1 ~ as.numeric(difftime(!!sym(time_col), DNADate, units = "days")) / 365.25,
        !!sym(disease_col) == 0 ~ as.numeric(difftime(lastRecordDate, DNADate, units = "days")) / 365.25
      )
    ) %>%
    filter(surv_time > 0, !is.infinite(surv_time), !is.na(surv_time)) 

  total_samples <- sum(data_surv[[disease_col]] == 1)
  if (total_samples < 10) {  # Filter the sample size
    warning(paste("Not enough events for", disease_col, "- Skipping."))
    next
  }

  surv_obj <- Surv(time = data_surv$surv_time, event = data_surv[[disease_col]])

  cox_result <- tryCatch({
    cox_model <- coxph(surv_obj ~ mca_auto + YAgeAtDNA + age2 + SmokingStatus + computed_gender.y, data = data_surv)

    tidy(cox_model) %>%
      filter(term == "mca_auto") %>%
      mutate(
        total_samples = total_samples,
        log_p_value = -log10(p.value)
      )
  }, error = function(e) {
    warning(paste("Cox model failed for", disease_col, "with error:", e$message))
    return(NULL)
  })

  if (!is.null(cox_result)) {
    cox_results[[disease_col]] <- cox_result
  }
}

results_df <- bind_rows(cox_results, .id = "disease")
results_df <- results_df %>% arrange(desc(log_p_value))
results_df

-log10(0.05/1952)

write.csv(results_df, "PheWAS_automCA_cox_BioVU_AFR.csv", row.names = FALSE)
system("gsutil cp PheWAS_automCA_cox_BioVU_AFR.csv gs://bicklab-main-storage/Users/Kun_Zhao")

table(data_biovu2$Race)

race <- data_biovu2[,c("ourSid","Race")]

data_all_pheno_cox3 <- merge(data_all_pheno_cox3, race, by.x = "person_id", by.y ="ourSid", all.x = T)

table(data_all_pheno_cox3$mca_auto)

table(data_all_pheno_cox3$Race,data_all_pheno_cox3$ancestry)

data_eas_pheno <- filter(data_all_pheno_cox3, data_all_pheno_cox3$Race == "A")
dim(data_eas_pheno)

#-PheWAS mCA by COX regression-----
cox_results <- list()
disease_cols <- colnames(data_all_pheno_cox3)[2:3514] #2:3514

for (i in seq_along(disease_cols)) {
  disease_col <- disease_cols[i]
  time_col <- paste0("time_", disease_col)

  cat("Processing:", time_col, "（Loop", i, "of", length(disease_cols), "）\n")
    
  if (!any(time_col %in% colnames(data_eas_pheno))) {
    warning(paste("Time column", time_col, "does not exist! Skipping", disease_col))
    next
  }
  
  unique_values <- unique(data_eas_pheno[[disease_col]])
  if (!all(unique_values %in% c(0, 1))) {
    warning(paste("Column", disease_col, "contains non-binary values. Skipping this column."))
    next
  }

  data_surv <- data_eas_pheno %>%
    mutate(
      surv_time = case_when(
        !!sym(disease_col) == 1 ~ as.numeric(difftime(!!sym(time_col), DNADate, units = "days")) / 365.25,
        !!sym(disease_col) == 0 ~ as.numeric(difftime(lastRecordDate, DNADate, units = "days")) / 365.25
      )
    ) %>%
    filter(surv_time > 0, !is.infinite(surv_time), !is.na(surv_time)) 

  total_samples <- sum(data_surv[[disease_col]] == 1)
  if (total_samples < 5) {  # Filter the sample size
    warning(paste("Not enough events for", disease_col, "- Skipping."))
    next
  }

  surv_obj <- Surv(time = data_surv$surv_time, event = data_surv[[disease_col]])

  cox_result <- tryCatch({
    cox_model <- coxph(surv_obj ~ mca_auto + YAgeAtDNA + age2 + SmokingStatus + computed_gender.y, data = data_surv)

    tidy(cox_model) %>%
      filter(term == "mca_auto") %>%
      mutate(
        total_samples = total_samples,
        log_p_value = -log10(p.value)
      )
  }, error = function(e) {
    warning(paste("Cox model failed for", disease_col, "with error:", e$message))
    return(NULL)
  })

  if (!is.null(cox_result)) {
    cox_results[[disease_col]] <- cox_result
  }
}

results_df <- bind_rows(cox_results, .id = "disease")
results_df <- results_df %>% arrange(desc(log_p_value))
results_df

write.csv(results_df, "PheWAS_automCA_cox_BioVU_EAS.csv", row.names = FALSE)
system("gsutil cp PheWAS_automCA_cox_BioVU_EAS.csv gs://bicklab-main-storage/Users/Kun_Zhao")

###CHIP PheWAS
chip_calls <- read.csv("../../chip_calls_bioVU.csv")

chip_calls

chip_calls<-chip_calls[!duplicated(chip_calls$GRID),]
dim(chip_calls)

jak2 <- filter(chip_calls, chip_calls$Gene.refGene == "JAK2")

jak2$test <- NA
jak2$test <- ifelse(jak2$GRID %in% data_final$person_id, 1, 0)
table(jak2$test)

a <- load("data_final.Rdata")
a

data_final

data_biovu2

data_biovu2$chip <- NA
data_biovu2$chip <- ifelse(data_biovu2$ourSid %in% chip_calls$GRID,1,0)
table(data_biovu2$chip)

data_biovu2$age2 <- data_biovu2$YAgeAtDNA * data_biovu2$YAgeAtDNA

data_all_pheno_cox <- merge(data_biovu2,data_final,by.x ="ourSid", by.y = "person_id",all.x = T)

dim(data_all_pheno_cox)

data_all_pheno_cox$lastRecordDate <- as.Date(data_all_pheno_cox$lastRecordDate, format = "%m/%d/%y")
data_all_pheno_cox$DNADate <- as.Date(data_all_pheno_cox$DNADate, format = "%m/%d/%y")

dim(data_biovu2)

head(data_all_pheno_cox)

cox_results <- list()

#-PheWAS CHIP by COX regression-----
disease_cols <- colnames(data_all_pheno_cox)[3481:3535] #23:3535
disease_cols

for (disease_col in disease_cols) {
  time_col <- paste0("time_", disease_col)
  
  # Print "time_col"
  print(time_col)
  
  # Check time_col
  if (any(time_col %in% colnames(data_all_pheno_cox))) {
    
    data_surv <- data_all_pheno_cox %>%
      mutate(
        surv_time = case_when(
          !!sym(disease_col) == 1 ~ as.numeric(difftime(!!sym(time_col), DNADate, units = "days")) / 365.25,
          !!sym(disease_col) == 0 ~ as.numeric(difftime(lastRecordDate, DNADate, units = "days")) / 365.25
        )
      ) %>%
      filter(surv_time > 0)
    
    total_samples <- sum(data_surv[[disease_col]] == 1)
    chip_status_1_outcome <- sum(data_surv$chip == 1 & data_surv[[disease_col]] == 1)
    
    surv_obj <- Surv(time = data_surv$surv_time, event = data_surv[[disease_col]])
    
    cox_model <- coxph(surv_obj ~ chip + YAgeAtDNA + age2 + SmokingStatus + computed_gender + Ethnicity, data = data_surv)
    
    cox_result <- tidy(cox_model) %>%
      filter(term == "chip") %>%
      mutate(
        total_samples = total_samples,
        chip_status_1_outcome = chip_status_1_outcome,
        log_p_value = -log10(p.value)  # -log10(p.value)
      )
    
    # Save
    cox_results[[disease_col]] <- cox_result
  } else {
    warning(paste("The", time_col, "Not exist!"))
  }
}
results_df <- bind_rows(cox_results, .id = "disease")

results_df

results_df2<- distinct(results_df)
results_df2

write.csv(results_df, "PheWAS_CHIP_cox_BioVU.csv", row.names = FALSE)

results_df <- read.csv("PheWAS_CHIP_cox_BioVU.csv")

results2_df <- filter(results_df,is.na(results_df$estimate) == F)
dim(results2_df)

results2_df <- results2_df %>%
  mutate(log_p_value = -log10(p.value),  
         outcome_index = row_number(),  
         shape = ifelse(estimate > 0, 24, 25)) 


top_10_outcomes <- results2_df %>%
  arrange(desc(log_p_value)) %>%
  slice_head(n = 30) %>%
  pull(disease) 


results2_df

phecode <- read.csv("phecodeX_ICD_CM_map_flat.csv")
phecode <- phecode[,c("phecode_string","category")]
phecode <- distinct(phecode)
phecode

results2_df <- merge(results2_df,phecode,by.x = "disease", by.y = "phecode_string" ,all.x = T)
table(results2_df$category, useNA = "always")

results2_df <- results2_df %>%
  mutate(top_10_label = ifelse(disease %in% top_10_outcomes, disease, NA))
threshold_pvalue_0_05 <- 4.3  

results2_df <- results2_df %>%
  arrange(category, outcome_index) %>%
  mutate(outcome_index = row_number()) 

results2_df <- filter(results2_df,is.na(results2_df$category) == F)


table(results2_df$top_10_label)

results2_df$top_10_label <- ifelse(nchar(results2_df$top_10_label)>30,NA,results2_df$top_10_label)

table(results2_df$top_10_label)

summary(results2_df$log_p_value)
results2_df$log_p_value<- ifelse(results2_df$log_p_value >40, 40, results2_df$log_p_value)
summary(results2_df$log_p_value)

library(ggrepel)
ggplot(results2_df, aes(x = outcome_index, y = log_p_value)) +
  geom_point(aes(color = category, shape = as.factor(shape), fill = category), size = 3) +
  scale_shape_manual(values = c(24, 25)) +
  geom_hline(yintercept = threshold_pvalue_0_05, linetype = "dashed", color = "red") +
  geom_text_repel(aes(label = top_10_label), size = 2.5, max.overlaps = 20, min.segment.length = 0) +  # ggrepel set positions of labels
  theme_minimal() +
  labs(x = "Category", y = "-log10(p.value)", title = "Manhattan Plot of CHIP PheWAS in BioVU") +
  scale_color_viridis_d(name = "Category", option = "plasma") +  # viridis, magma, plasma, inferno, cividis
  scale_fill_viridis_d(name = "Category", option = "plasma") +  # Same (fill colour)
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_continuous(breaks = results2_df %>% group_by(category) %>% summarize(mid = mean(outcome_index)) %>% pull(mid),
                     labels = results2_df %>% group_by(category) %>% summarize(cat = first(category)) %>% pull(cat))


library(ggrepel)
p <- ggplot(results2_df, aes(x = outcome_index, y = log_p_value)) +
  geom_point(aes(color = category, shape = as.factor(shape), fill = category), size = 1) +
  scale_shape_manual(values = c(24, 25)) +
  geom_hline(yintercept = threshold_pvalue_0_05, linetype = "dashed", color = "red") +
  geom_text_repel(aes(label = top_10_label), size = 1, max.overlaps = 20, min.segment.length = 0) +  # ggrepel set positions of labels
  theme_minimal() +
  labs(x = "", y = "-log10(p.value)", title = "Manhattan Plot of CHIP PheWAS in BioVU") +
  scale_color_viridis_d(name = "Category", option = "plasma") +  # viridis, magma, plasma, inferno, cividis
  scale_fill_viridis_d(name = "Category", option = "plasma") +  # Same (fill colour)
  theme(legend.position = "none",  # Remove legend
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_continuous(breaks = results2_df %>% group_by(category) %>% summarize(mid = mean(outcome_index)) %>% pull(mid),
                     labels = results2_df %>% group_by(category) %>% summarize(cat = first(category)) %>% pull(cat))
p
# Adjust the plot size and save it as a file
ggsave("manhattan_plot.png", plot = p, width = 8, height = 2, units = "in")

library(ggplot2)
library(ggrepel)
library(dplyr)

p <- ggplot(results2_df, aes(x = outcome_index, y = log_p_value)) +
  geom_point(aes(color = category, shape = as.factor(shape), fill = category), size = 1) +  # Set size and fill
  scale_shape_manual(values = c(24, 25), labels = c("Risk", "Protective"), guide = guide_legend(override.aes = list(fill = c("black", "black")))) +  # Use solid shapes
  geom_hline(yintercept = threshold_pvalue_0_05, linetype = "dashed", color = "red") +
  geom_text_repel(aes(label = top_10_label), size = 1, max.overlaps = 20, min.segment.length = 0) +  # Small text for labels
  theme_minimal() +
  labs(x = "", y = "-log10(p.value)", title = "Manhattan Plot of CHIP PheWAS in BioVU") +
  scale_color_viridis_d(name = "Category", option = "plasma", guide = "none") +  # Remove category legend
  scale_fill_viridis_d(name = "Category", option = "plasma", guide = "none") +   # Remove fill legend
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 3),  # Smaller x-axis text size
    axis.text.y = element_text(size = 3),  # Smaller y-axis text size
    axis.title.y = element_text(size = 4),  # Smaller y-axis title
    plot.title = element_text(hjust = 0.5, size = 6),  # Smaller title size
    legend.position = c(0.95, 0.95),  # Place legend inside the plot (top-right corner)
    legend.background = element_rect(fill = alpha('white', 0.5)),  # Set transparent background for legend
    legend.key.size = unit(0.2, "cm"),  # Adjust the size of the legend keys
    legend.title = element_text(size = 3),  # Smaller legend title (Effect Direction)
    legend.text = element_text(size = 3)  # Smaller legend text (Risk/Protective)
  ) +
  scale_x_continuous(breaks = results2_df %>% group_by(category) %>% summarize(mid = mean(outcome_index)) %>% pull(mid),
                     labels = results2_df %>% group_by(category) %>% summarize(cat = first(category)) %>% pull(cat)) +
  guides(shape = guide_legend(title = "Effect Direction", override.aes = list(size = 1)))  # Show shape legend with custom title

p

ggsave("manhattan_plot.png", plot = p, width = 6, height = 1.5, units = "in")

x <- load("data_all_pheno_cox_BioVU.Rdata")
x

data_all_pheno_cox2 <- filter(data_all_pheno_cox,data_all_pheno_cox$YAgeAtDNA >=40 & data_all_pheno_cox$YAgeAtDNA <=90)

dim(data_all_pheno_cox2) #51202 × 7064

data_all_pheno_cox2$lastRecordDate <- as.Date(data_all_pheno_cox2$lastRecordDate, format = "%m/%d/%y")
data_all_pheno_cox2$DNADate <- as.Date(data_all_pheno_cox2$DNADate, format = "%m/%d/%y")

#-PheWAS mCA by COX regression-----
cox_results <- list()
disease_cols <- colnames(data_all_pheno_cox2)[2:3514] #2:3514

table(data_all_pheno_cox2$`Acute kidney failure`,useNA = "always")

biovu_mca <- read.csv("BioVU_mca_type_data.csv")

data_all_pheno_cox3 <- merge(data_all_pheno_cox2,biovu_mca, by.x = "person_id", by.y = "sample_id", all.x = T)

dim(data_all_pheno_cox3)

target_mca_types <- colnames(data_all_pheno_cox3)[7077:7082]
target_mca_types

data_all_pheno_cox3[, 7065:7082] <- lapply(data_all_pheno_cox3[, 7065:7082], function(column) {
  ifelse(is.na(column), 0, column)
})

table(data_all_pheno_cox3$chr11p.CN.LOH,useNA = "always")

biovu_AFR <- fread("20200515_biallelic_mega_recalled.chr1-22.grid.AA.filt1.fam",header = FALSE)
biovu_EUR <- fread("20200518_biallelic_mega_recalled.chr1-22.grid.EU.filt1.fam",header = FALSE)

data_all_pheno_cox3$ancestry <- NA

data_all_pheno_cox3 <- data_all_pheno_cox3 %>% mutate(ancestry = ifelse(person_id %in% biovu_AFR$V2, "AFR", ancestry))
data_all_pheno_cox3 <- data_all_pheno_cox3 %>% mutate(ancestry = ifelse(person_id %in% biovu_EUR$V2, "EUR", ancestry))
table(data_all_pheno_cox3$ancestry,useNA = "always")

data_all_pheno_cox3$ancestry <- ifelse(is.na(data_all_pheno_cox3$ancestry),"Others",data_all_pheno_cox3$ancestry)

#-PheWAS mCA by COX regression-----
cox_results <- list()
disease_cols <- colnames(data_all_pheno_cox2)[2:3514] #2:3514

summary(data_all_pheno_cox3$age2)
table(data_all_pheno_cox3$computed_gender, useNA="always")

data_gender <- data_biovu2[,c("ourSid", "computed_gender")]

data_all_pheno_cox3 <- merge(data_all_pheno_cox3,data_gender, by.x = "person_id", by.y = "ourSid", all.x = T)

table(data_all_pheno_cox3$computed_gender.y, useNA = "always")

save(data_all_pheno_cox3, file = "data_pheno_BioVU_4090.Rdata")

library(survival)
library(dplyr)
library(broom)

cox_results_list <- list()

for (i in seq_along(target_mca_types)) {
  mca_type <- target_mca_types[i]
  
  cat("Processing MCA Type:", mca_type, "（Loop", i, "of", length(target_mca_types), "）\n")
  
  cox_results <- list() 
  
  for (j in seq_along(disease_cols)) {
    disease_col <- disease_cols[j]
    time_col <- paste0("time_", disease_col)
    
    cat("  Processing Disease:", disease_col, "（Sub-loop", j, "of", length(disease_cols), "）\n")
    
    unique_values <- unique(data_all_pheno_cox3[[disease_col]])
    if (!all(unique_values %in% c(0, 1))) {
      warning(paste("Column", disease_col, "contains non-binary values. Skipping this column."))
      next
    }
    
    # 构建生存时间
    data_surv <- data_all_pheno_cox3 %>%
      mutate(
        surv_time = case_when(
          !!sym(disease_col) == 1 ~ as.numeric(difftime(!!sym(time_col), DNADate, units = "days")) / 365.25,
          !!sym(disease_col) == 0 ~ as.numeric(difftime(lastRecordDate, DNADate, units = "days")) / 365.25
        )
      ) %>%
      filter(surv_time > 0)
    
    # 确保 mca_type 是数值型
    if (!is.numeric(data_surv[[mca_type]])) {
      data_surv[[mca_type]] <- as.numeric(data_surv[[mca_type]] != 0)  # 转换为 0/1
    }
    
    total_samples <- sum(data_surv[[disease_col]] == 1)
    
    surv_obj <- Surv(time = data_surv$surv_time, event = data_surv[[disease_col]])
    
    # 构建 Cox 模型，捕获错误
    cox_result <- tryCatch({
      formula <- as.formula(paste("surv_obj ~", mca_type, "+ YAgeAtDNA + age2 + SmokingStatus + computed_gender.y + ancestry"))
      cox_model <- coxph(formula, data = data_surv)
      
      # 提取结果
      tidy(cox_model) %>%
        filter(term == mca_type) %>%
        mutate(
          total_samples = total_samples,
          log_p_value = -log10(p.value)
        )
    }, error = function(e) {
      warning(paste("Error in Cox model for disease:", disease_col, "and mca_type:", mca_type, "-", e$message))
      return(NULL)
    })
    
    if (!is.null(cox_result)) {
      cox_results[[disease_col]] <- cox_result
    }
  }
  
  # 合并所有 disease_col 的结果
  if (length(cox_results) > 0) {
    combined_result <- bind_rows(cox_results, .id = "Disease")
    
    # 保存为 CSV
    output_file <- paste0("results/", mca_type, "_BioVU.csv")  # 保存到 results 目录
    dir.create("results", showWarnings = FALSE)
    write.csv(combined_result, output_file, row.names = FALSE)
    
    # 存储结果到列表
    cox_results_list[[mca_type]] <- combined_result
  } else {
    warning(paste("No valid results for MCA Type:", mca_type))
  }
}

cat("Cox regression analysis completed. Results saved as individual CSV files in 'results' directory.\n")

biovu_AFR <- fread("20200515_biallelic_mega_recalled.chr1-22.grid.AA.filt1.fam",header = FALSE)
biovu_EUR <- fread("20200518_biallelic_mega_recalled.chr1-22.grid.EU.filt1.fam",header = FALSE)

data_all_pheno_cox3$ancestry <- NA

data_all_pheno_cox3 <- data_all_pheno_cox3 %>% mutate(ancestry = ifelse(person_id %in% biovu_AFR$V2, "AFR", ancestry))
data_all_pheno_cox3 <- data_all_pheno_cox3 %>% mutate(ancestry = ifelse(person_id %in% biovu_EUR$V2, "EUR", ancestry))
table(data_all_pheno_cox3$ancestry,useNA = "always")

data_all_pheno_cox3$ancestry <- ifelse(is.na(data_all_pheno_cox3$ancestry),"Others",data_all_pheno_cox3$ancestry)

#-PheWAS mCA by COX regression-----
cox_results <- list()
disease_cols <- colnames(data_all_pheno_cox2)[2:3514] #2:3514

summary(data_all_pheno_cox3$age2)
table(data_all_pheno_cox3$computed_gender, useNA="always")

x<- load("data_pheno_BioVU_4090.Rdata")
x

table(data_all_pheno_cox3$computed_gender.y, useNA = "always")

table(data_all_pheno_cox3$ancestry,useNA = "always")

data_all_pheno_cox3[, 7065:7082] <- lapply(data_all_pheno_cox3[, 7065:7082], function(column) {
  ifelse(is.na(column), 0, column)
})

library(readxl)
result_phewas_cox <- read_excel("meta_sepcific_mca_phewas_sig_P46.xlsx")

dim(result_phewas_cox)
head(result_phewas_cox)

table(result_phewas_cox$mca_type)

target <- filter(result_phewas_cox, result_phewas_cox$mca_type == "chr9p.CN.LOH_combined")

phenotype_sig <- target %>% distinct(phenotype)

phenotype_sig <- phenotype_sig$phenotype

length(phenotype_sig)

phenotype_sig <- phenotype_sig[phenotype_sig %in% colnames(data_all_pheno_cox3)]
length(phenotype_sig)

target_mca_types <- colnames(data_all_pheno_cox3)[7065:7082] #7065:7082
target_mca_types

target_mca_types <- "chr9p.CN.LOH"

#-PheWAS mCA by Firth Logistic regression-----
disease_cols <- phenotype_sig
disease_cols

library(logistf)
library(dplyr)
library(broom)

#dir.create("results_spe_mCA_phewas_firth_log")

firth_results_list <- list()

for (i in seq_along(target_mca_types)) {
  mca_type <- target_mca_types[i]
  
  cat("Processing MCA Type:", mca_type, "（Loop", i, "of", length(target_mca_types), "）\n")
  
  firth_results <- list()
  
  for (j in seq_along(disease_cols)) {
    disease_col <- disease_cols[j]
    
    cat("  Processing Disease:", disease_col, "（Sub-loop", j, "of", length(disease_cols), "）\n")
    
    unique_values <- unique(data_all_pheno_cox3[[disease_col]])
    if (!all(unique_values %in% c(0, 1))) {
      warning(paste("Column", disease_col, "contains non-binary values. Skipping."))
      next
    }
    
    data_logistf <- data_all_pheno_cox3 %>%
      filter(!is.na(!!sym(disease_col)), !is.na(!!sym(mca_type)))  # Remove NA
    
    # Change mCA types to numeric
    if (!is.numeric(data_logistf[[mca_type]])) {
      data_logistf[[mca_type]] <- as.numeric(data_logistf[[mca_type]] != 0)
    }
    
    total_cases <- sum(data_logistf[[disease_col]] == 1)
    
    # Model
    formula <- as.formula(paste0("`", disease_col, "` ~ ", mca_type,
                             " + YAgeAtDNA + age2 + SmokingStatus + computed_gender.y + ancestry"))
    # Firth logistic
    fit_result <- tryCatch({
      model <- logistf::logistf(formula, data = data_logistf)
  
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
    output_file <- paste0("results_spe_mCA_phewas_firth_log/", mca_type, "_FirthLogis_BioVU.csv")
    write.csv(combined_result, output_file, row.names = FALSE)
    
    firth_results_list[[mca_type]] <- combined_result
  } else {
    warning(paste("No valid results for MCA Type:", mca_type))
  }
}

cat("Firth logistic regression completed. Results saved in 'results_spe_mCA_phewas_firth_log' directory.\n")

ltl_prs <- fread("ltl_prs_scores_biovu.tsv")

dim(ltl_prs)
ls(ltl_prs)

head(data_all_pheno_cox3)

ltl_data <- merge(ltl_prs,data_all_pheno_cox3, by.x = "IID", by.y = "person_id", all.x = T)

summary(ltl_data$TOTAL_PRS)

model <- lm(TOTAL_PRS ~ `Heart failure`+ YAgeAtDNA + age2 + SmokingStatus + computed_gender.y + ancestry, data = ltl_data)
summary(model)
