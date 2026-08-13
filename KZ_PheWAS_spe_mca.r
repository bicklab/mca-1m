library(logistf)
library(dplyr)
library(broom)

library(broom)
library(scales)
library(arrow)
library(survival)
library(tidyr)
library(dplyr)
library(haven)
library(data.table)
library(lubridate)
library(ggplot2)

data_all_pheno_cox3 <- readRDS("data_all_pheno_cox3_UKB.rds")
dim(data_all_pheno_cox3)

target_mca_types <- colnames(data_all_pheno_cox3)[4842:4859]
target_mca_types

result_phewas_cox <- read.csv("meta_spe_mCA_BUAT.csv")

dim(result_phewas_cox)
head(result_phewas_cox)

phenotype_sig <- result_phewas_cox %>% distinct(phenotype)

table(result_phewas_cox$mca_type)

target <- filter(result_phewas_cox, result_phewas_cox$mca_type == "chr9p.CN.LOH")

phenotype_sig <- target$phenotype
length(phenotype_sig)
phenotype_sig

phenotype_sig <- phenotype_sig[phenotype_sig %in% colnames(data_all_pheno_cox3)]
length(phenotype_sig)

#-PheWAS mCA by Firth Logistic regression-----
disease_cols <- phenotype_sig
disease_cols

target_mca_types <- "chr9p.CN.LOH"

firth_results_list <- list()

for (i in seq_along(target_mca_types)) {
  mca_type <- target_mca_types[i]
  
  cat("Processing MCA Type:", mca_type, "（Loop", i, "of", length(target_mca_types), "）\n")
  
  firth_results <- list()
  
  for (j in seq_along(disease_cols)) {
    disease_col <- disease_cols[j]
    
    cat(" Processing Disease:", disease_col, "（Sub-loop", j, "of", length(disease_cols), "）\n")
    
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
                             "+ baseline_age + age2 + smoking_0 + genetic_sex + PC1 + PC2 + PC3 + PCD4 + PC5"))
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
    output_file <- paste0("results_spe_mCA_phewas_firth_log/", mca_type, "_FirthLogis_UKB.csv")
    write.csv(combined_result, output_file, row.names = FALSE)
    
    firth_results_list[[mca_type]] <- combined_result
  } else {
    warning(paste("No valid results for MCA Type:", mca_type))
  }
}

cat("Firth logistic regression completed. Results saved in 'results_spe_mCA_phewas_firth_log' directory.\n")

combined_result


