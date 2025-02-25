#Combine 4 datasets of mCA
#-1. Create & combine datasets
#---1.1 mCA-----
ukb_mca <- fread('mCA_calls_ukb.tsv') #17865
data_covari <- fread('demographic_data_Aug15_participant.tsv')
data_age <- data_covari[,c("ID_VUMC","baseline_age")]
ukb_mca <- ukb_mca[,c(3:24)]
ukb_mca <- merge(ukb_mca,data_age,by="ID_VUMC",all.x = T)
ukb_mca$p_count<-NA
ukb_mca$q_count<-NA
ukb_mca$p_count<-ifelse(ukb_mca$p_arm=="Y"|ukb_mca$p_arm=="T"|ukb_mca$p_arm=="C",1,0)
ukb_mca$q_count<-ifelse(ukb_mca$q_arm=="Y"|ukb_mca$q_arm=="T"|ukb_mca$q_arm=="C",1,0)
table(ukb_mca$p_count,ukb_mca$q_count)
ukb_mca$arm<-NA
ukb_mca$arm<-ifelse(ukb_mca$p_count ==1 & ukb_mca$q_count ==0 ,"p",ukb_mca$arm)
ukb_mca$arm<-ifelse(ukb_mca$q_count ==1 & ukb_mca$p_count ==0 ,"q",ukb_mca$arm)
table(ukb_mca$arm)
ukb_mca$chr_type <- paste(ukb_mca$chrom,ukb_mca$arm,sep = "")
ukb_mca <- ukb_mca %>% mutate(chr_type = str_replace(chr_type, "NA$", ""))
names(ukb_mca)[23] <- "age"
names(ukb_mca)[1] <- "sample_id"
ukb_mca$chr_type2<-ukb_mca$chr_type
ukb_mca$chr_type2<- ifelse(ukb_mca$chr_type2=="chr15q","chr15",ukb_mca$chr_type2)
ukb_mca$cohort <- "UKB"
ukb_mca$ancestry <- "EUR"
ls(ukb_mca)
#-----Combine 4 cohorts-----
data_4cohorts <- merge(data_TOP_AoU_bioVU,ukb_mca,all = T)
table(data_4cohorts$cohort) #51,094
data_4cohorts$mca <- 1
save(data_4cohorts,file = "data_4cohorts.Rdata")
#-----Plot: Age & male plot-----
data_filtered <- data_4cohorts %>% filter(type != "Undetermined")
data_filtered <- filter(data_filtered,age>40 & age <90)
table(data_filtered$cohort)
table(data_filtered$type)
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
data_summary<-filter(data_summary,data_summary$n>=100)


data_summary <- data_summary %>%
  mutate(show_label = ifelse(chr_type2 %in% c("chr15","chr20q","chr12"), TRUE, FALSE))

data_summary<-filter(data_summary,data_summary$chr_type2 != "chrX" & data_summary$chr_type2 != "chrXp" & data_summary$chr_type2 != "chrXq")

ggplot(data_summary, aes(x = mean_age, y = mean_fraction_male, color = type, alpha = n)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_fraction_male - sem_fraction_male, ymax = mean_fraction_male + sem_fraction_male), width = 0.2) +
  geom_errorbarh(aes(xmin = mean_age - sem_age, xmax = mean_age + sem_age), height = 0.2) +
  geom_text(data = filter(data_summary, show_label), aes(label = chr_type2), hjust = 1.1, vjust = 1.1) +
  scale_color_manual(values = c("Gain" = "#56B4E9", "Loss" = "#E69F00", "CN-LOH" = "#009E73")) +
  scale_alpha_continuous(range = c(0.5, 1)) +
  labs(title = "Carriers of different mCA types with different age and sex distributions",
       x = "Mean age (years)",
       y = "Fraction male",
       alpha = "Sample Size") +
  theme_minimal() +
  theme(legend.position = "right",
        legend.title = element_blank())
#-----Plot: Gender -----
# Gender calculate
data_filtered <- filter(data_4cohorts,age>40 & age<90)
table(data_filtered$computed_gender, useNA="always")
counts <- data_filtered %>%
  group_by(computed_gender, chr_type2,type) %>%
  summarise(count = n()) %>%
  ungroup()

# 计算每种gender的总计数
total_counts <- data_filtered %>%
  group_by(computed_gender) %>%
  summarise(total = n()) %>%
  ungroup()

# 将总计数合并到计数数据框中
counts <- counts %>%
  left_join(total_counts, by = "computed_gender") %>%
  mutate(proportion = (count / total) * 100)  # 计算比例，转化为百分比

# 筛选male和female的比例并合并到一个数据框中
combined_data <- filter(counts,counts$computed_gender=="F" | counts$computed_gender=="M" )
combined_data <- filter(combined_data, combined_data$type!= "Undetermined" )
combined_data <- filter(combined_data, combined_data$chr_type2!= "chrX" & combined_data$chr_type2!= "chrXp" & combined_data$chr_type2!= "chrXq")

male_data <- combined_data %>%
  filter(computed_gender == "M") %>%
  select(chr_type2,type, male_proportion = proportion, count_male = count)

female_data <- combined_data %>%
  filter(computed_gender == "F") %>%
  select(chr_type2,type, female_proportion = proportion, count_female = count)

combined_data <- merge(female_data, male_data, by = c("chr_type2","type"))
combined_data <- filter(combined_data,combined_data$count_male >= 20 | combined_data$count_female >= 20)

combined_data$count <- combined_data$count_male + combined_data$count_female

ggplot(combined_data, aes(x = male_proportion, y = female_proportion, color = type, label = chr_type2, size = count)) +
  geom_point() +  
  geom_text(vjust = -0.5, hjust = 0.5) +  # chr_type
  scale_x_continuous(name = "Male Proportion (%)", limits = c(0, 3)) +
  scale_y_continuous(name = "Female Proportion (%)", limits = c(0, 3)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey") + 
  scale_color_manual(values = c("Gain" = "#56B4E9", "Loss" = "#E69F00", "CN-LOH" = "#009E73")) +  
  scale_size_continuous(range = c(1, 5), name = "Count") +  # Point Size
  theme_minimal() +
  labs(title = "Chromosome Type Proportion in Male vs Female",
       subtitle = "Colored by Chromosome Event Type")

#Combine 4 datasets of CHIP
#---1.2 CHIP----
#-----Combine 4 cohort-----
#ukb
ukb_chip <- fread('chip_calls_ukb.txt') #17026
ukb_chip <- ukb_chip[,1:14]
ukb_chip2 <- merge(ukb_chip,vumc_id,by.x = "ID",by.y = "ID_Broad",all.x = T)
names(ukb_chip2)[15]<-"person_id"
ukb_chip2$cohort <- "UKB"
ukb_chip2 <- ukb_chip2[,2:16]
ukb_chip <-ukb_chip2
#AoU
chip <- fread('CHIP_calls_batch2_nomyeloidcancers_06162023.txt') 
chip2 <- fread('CHIP_calls_nomyeloidcancers_12082022.txt') 
chip_carriers <- merge(chip,chip2,all = T)
aou_chip <- distinct(chip_carriers) #10464
aou_chip <- aou_chip[,c(1:8,10,28,29,38,39,44)]
aou_chip$cohort <- "AoU"
ls(aou_chip)
#BioVU
biovu_chip <- read.csv(file = "chip_calls_bioVU.csv")
biovu_chip<-biovu_chip[,c(1:7,9,27,28,36,39,40)]
names(biovu_chip)[11]<-"person_id"
biovu_chip$cohort <- "BioVU"
biovu_chip$minAD <- NA
ls(biovu_chip)
str(biovu_chip)
biovu_chip$AF <- as.numeric(biovu_chip$AF)
#Merge
aou_chip$Start <- as.integer(aou_chip$Start)
aou_chip$End <- as.integer(aou_chip$End)
aou_chip$AF <- as.numeric(aou_chip$AF)
data_3cohorts_chip <- merge(ukb_chip,aou_chip,all = T)
str(data_3cohorts_chip)
data_3cohorts_chip$person_id <- as.character(data_3cohorts_chip$person_id)
data_3cohorts_chip <- merge(data_3cohorts_chip,biovu_chip,all = T) #33240
data_3cohorts_chip$chip <- 1
data_3cohorts_chip$person_id.x <- ifelse(is.na(data_3cohorts_chip$person_id.x),data_3cohorts_chip$person_id.y, data_3cohorts_chip$person_id.x)
data_3cohorts_chip$cohort.x <- ifelse(is.na(data_3cohorts_chip$cohort.x),data_3cohorts_chip$cohort.y, data_3cohorts_chip$cohort.x)
data_3cohorts_chip <- data_3cohorts_chip[,c(1:15,18)]
names(data_3cohorts_chip)[14] <- "person_id"
names(data_3cohorts_chip)[15] <- "cohort"
save(data_3cohorts_chip,file = "data_3cohorts_chip.Rdata")

#TOPMED
topmed_chip <- fread('TOPMed_CHIPcalls.tsv') #6158
topmed_chip <- topmed_chip[,c("ALT","CHROM","ExonicFunc","Gene","POS","REF","Sample","VAF")]
topmed_chip$DP<-NA
topmed_chip$End<-NA
topmed_chip$Func.refGene<-NA
topmed_chip$minAD<-NA
topmed_chip$NonsynOI<-NA
topmed_chip$transcriptOI<-NA
names(topmed_chip)[1]<-"Alt"
names(topmed_chip)[2]<-"Chr"
names(topmed_chip)[3]<-"ExonicFunc.refGene"
names(topmed_chip)[4]<-"Gene.refGene"
names(topmed_chip)[5]<-"Start"
names(topmed_chip)[6]<-"Ref"
names(topmed_chip)[7]<-"person_id"
names(topmed_chip)[8]<-"AF"
topmed_chip$cohort <- "TOPMed"
topmed_chip$chip <- 1
data_4cohorts_chip <- merge(data_3cohorts_chip,topmed_chip,all = T) #33240
data_4cohorts_chip$chip <-1
data_4cohorts_chip$person_id.x <- ifelse(is.na(data_4cohorts_chip$person_id.x),data_4cohorts_chip$person_id.y, data_4cohorts_chip$person_id.x)
data_4cohorts_chip$cohort.x <- ifelse(is.na(data_4cohorts_chip$cohort.x),data_4cohorts_chip$cohort.y, data_4cohorts_chip$cohort.x)
data_4cohorts_chip <- data_4cohorts_chip[,c(1:15,20)]
names(data_4cohorts_chip)[14] <- "person_id"
names(data_4cohorts_chip)[15] <- "cohort"
save(data_4cohorts_chip,file = "data_4cohorts_chip.Rdata")


#-2. Combine mCA and CHIP------
ls(data_4cohorts)
ls(data_4cohorts_chip)
names(data_4cohorts)[1] <- "person_id"
mCA_4c <- data_4cohorts[,c("person_id","cohort","type","chrom","arm","chr_type","chr_type2","mca")]
CHIP_4c <- data_4cohorts_chip[,c("person_id","chip","cohort","Gene.refGene","Chr","AF")]
combined_data <- inner_join(mCA_4c, CHIP_4c, by = c("person_id","cohort"), relationship = "many-to-many")
combined_data <- filter(combined_data, combined_data$chrom != "chrX" & combined_data$chrom != "X")
combined_data <- filter(combined_data, combined_data$type != "Undetermined")
table(combined_data$cohort)
#-calculate OR
combined_data$chr_type3 <- paste(combined_data$chr_type2, combined_data$type, sep = ":")
contingency_table <- table(combined_data$Gene.refGene, combined_data$chr_type3)
table(combined_data$Gene.refGene)
#install.packages("epitools")
library(epitools)

gene_chip <- c("ASXL1","CBL","CREBBP","DNMT3A","EP300","GNB1","JAK2","KRAS","MPL","PPM1D","SF3B1","SRSF2","TET2","TP53","U2AF1")
# 初始化存储 OR、P 值和 Co-occurrence 结果
or_results <- list()

# 循环遍历每个基因和 mCA 类型的组合
for (Gene.refGene in gene_chip) {
  for (chr_type3 in colnames(contingency_table)) {
    # 提取 2x2 表
    sub_table <- matrix(c(
      contingency_table[Gene.refGene, chr_type3],             # 共现的个体数
      sum(contingency_table[Gene.refGene,]) - contingency_table[Gene.refGene, chr_type3],  # 基因存在但mCA不存在
      sum(contingency_table[, chr_type3]) - contingency_table[Gene.refGene, chr_type3],  # mCA存在但基因不存在
      sum(contingency_table) - sum(contingency_table[Gene.refGene,]) - sum(contingency_table[, chr_type3]) + contingency_table[Gene.refGene, chr_type3]
    ), nrow = 2)
    
    # 检查是否有零单元格
    if (all(sub_table > 0)) {
      # 计算 OR 和 Fisher P 值
      or <- oddsratio(sub_table)
      p_value <- fisher.test(sub_table)$p.value
      co_occurrence_count <- contingency_table[Gene.refGene, chr_type3]
      
      # 保存 OR、p 值和 Co-occurrence 个体数
      or_results[[paste(Gene.refGene, chr_type3, sep = "_")]] <- list(or = or$measure[2, 1], p_value = p_value, co_occurrence = co_occurrence_count)
    }
  }
}

# Convert the results to dataframe
final_results <- do.call(rbind, lapply(names(or_results), function(name) {
  data.frame(
    gene = unlist(strsplit(name, "_"))[1],
    mca_type = unlist(strsplit(name, "_"))[2],
    or = or_results[[name]]$or,
    p_value = or_results[[name]]$p_value,
    co_occurrence = or_results[[name]]$co_occurrence  # 添加 Co-occurrence 个体数
  )
}))

# Excluded OR=NA and Bonferroni adjusted
final_results <- final_results[!is.na(final_results$or), ]
final_results$adjusted_p_value <- p.adjust(final_results$p_value, method = "bonferroni")
final_results2 <- filter(final_results,final_results$adjusted_p_value<=0.05)
final_results2 <- final_results2 %>%
  separate(mca_type, 
           into = c("chromosome", "arm_alteration"), 
           sep = "(?<=\\d)(?=[pq]:|:)",  
           remove = FALSE) %>% 
  separate(arm_alteration, 
           into = c("arm", "alteration_type"), 
           sep = ":", 
           fill = "right", 
           remove = TRUE) 

#Draw plot
final_results2$pattern <- ifelse(final_results2$arm == "p", "stripe", "dot")
alteration_colors <- c("Gain" = "#56B4E9", "Loss" = "#E69F00", "CN-LOH" = "#009E73")

library(ggplot2)
library(dplyr)
#devtools::install_github("coolbutuseless/ggpattern")
library(ggpattern)

chromosome_order <- paste0("chr", 1:22)

ggplot(final_results2, aes(x = chromosome, y = gene, fill = alteration_type)) +
  geom_tile(aes(alpha = log(or)), color = "white", size = 0.5) + 
  scale_fill_manual(values = alteration_colors) + 
  scale_alpha_continuous(range = c(0.2, 1)) +  
  geom_tile_pattern(aes(pattern = pattern), 
                    fill = NA, color = "black", size = 0.2, 
                    pattern_fill = NA, pattern_density = 0.1, 
                    pattern_angle = 45) + 
  scale_x_discrete(limits = chromosome_order) +  
  theme_minimal() +
  labs(title = "Co-occurrence of Mutations by Chromosome and Gene",
       x = "mCA by Chromosome", y = "CHIP Genes", fill = "Alteration Type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

results_show <- final_results2[,c(1,3:6,8,9)]

