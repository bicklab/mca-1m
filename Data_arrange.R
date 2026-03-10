#Allofus mCA paper#
#-Library Packages-------
library(survival)
library(tidyr)
library(dplyr)
library(tidyverse)
library(data.table)
library(lubridate)
install.packages("survminer")
library(survminer)
install.packages("doBy")
library(doBy)
library(haven)
install.packages("Publish")
library(Publish)
install.packages("epiDisplay")
library(epiDisplay)
library(reshape2)
install.packages("skimr")
library(skimr)
install.packages("DataExplorer")
library(DataExplorer)
install.packages("Gmisc")
library(Gmisc)
library(ggplot2)
install.packages("psych")
library(psych)
install.packages("forestplot")
library(forestplot)
#Regression
library(broom)
install.packages("kableExtra")
library(kableExtra)
library(scales)
library(stringr)
library(ggrepel)
#----1.Load data#-------
#mCAs#
system("gsutil cp gs://fc-secure-cb192ac6-30ba-46b9-92ee-896a6e36c63e/yp_mca_calls_021924.txt .", intern = TRUE)
mca_df <- fread('yp_mca_calls_021924.txt') #19799
#data_condition <- fread('data_condition.txt') 
head(mca_df, 10)
mca_df$mca<-1
mca_df$cf_cat<-NA
mca_df$cf_cat[0 < mca_df$cf & mca_df$cf < 0.1] <- 1
mca_df$cf_cat[mca_df$cf >= 0.1] <- 2
table(mca_df$type,mca_df$cf_cat,useNA = "always")
#High risk mCAs for CLL include loss of chr6/chr11/chr13/chr17, gain of chr12, or CNLOH of chr13
mca_df$mac_high<-NA
mca_df$mac_high<-ifelse(mca_df$type=="Loss" & mca_df$chrom =="chr6",1,mca_df$mac_high)
mca_df$mac_high<-ifelse(mca_df$type=="Loss" & mca_df$chrom =="chr11",1,mca_df$mac_high)
mca_df$mac_high<-ifelse(mca_df$type=="Loss" & mca_df$chrom =="chr13",1,mca_df$mac_high)
mca_df$mac_high<-ifelse(mca_df$type=="Loss" & mca_df$chrom =="chr17",1,mca_df$mac_high)
mca_df$mac_high<-ifelse(mca_df$type=="Gain" & mca_df$chrom =="chr12",1,mca_df$mac_high)
mca_df$mac_high<-ifelse(mca_df$type=="CN-LOH" & mca_df$chrom =="chr13",1,mca_df$mac_high)
table(mca_df$mac_high) #1305
id_mcahigh<-filter(mca_df,mca_df$mac_high==1)
id_mcahigh<-id_mcahigh[!duplicated(id_mcahigh$sample_id),] #997
id_mcahigh<-id_mcahigh[,c("sample_id","mac_high")]
#id_unique
mca_df2 <- mca_df %>% group_by(sample_id) %>% mutate(Count = n()) %>% ungroup()
mca_cfmax <- mca_df2 %>% group_by(sample_id) %>% summarize(cf_max = max(cf), .groups = 'drop')
data_cfmax <- mca_df2 %>% left_join(mca_cfmax, by = "sample_id")
mca_idunique<- data_cfmax %>% group_by(sample_id) %>% slice(1) %>% ungroup()
mca_idunique$cf_cat<-NA
mca_idunique$cf_cat[0 < mca_idunique$cf_max & mca_idunique$cf_max < 0.1] <- 1
mca_idunique$cf_cat[mca_idunique$cf_max >= 0.1] <- 2 #14678
#Revision
data_cfmax_unique <- data_cfmax %>%
  group_by(sample_id) %>%
  slice_max(cf, with_ties = FALSE) %>%
  ungroup()
ls(data_cfmax_unique)
cols_to_replace <- c("baf_conc", "bdev", "bdev_se", "beg_GRCh38", "cf", "cf_max", "chrom", 
                     "computed_gender", "Count", "dnmt3a", "end_GRCh38", "length", 
                     "lod_baf_conc", "lod_baf_phase", "lod_lrr_baf", "n_flips", 
                     "n_hets", "n_sites", "n50_hets", "p_arm", "q_arm", "rel_cov", 
                     "rel_cov_se", "SV", "tet2", "type")

data_cfmax_subset <- data_cfmax_unique %>%  select(sample_id, all_of(cols_to_replace))
names(data_cfmax_subset)[1] <- "person_id"
data_all_replaced <- data_all %>%
  select(-all_of(cols_to_replace)) %>%
  left_join(data_cfmax_subset, by = "person_id")
summary(data_all_replaced$cf_max)
summary(data_all$cf_max)
data_all_replaced<-data_all
save(data_all,file = "data_all.Rdata")
#merge mca_high_risk
mca_idunique<-merge(mca_idunique,id_mcahigh,by="sample_id",all.x = T)
mca_idunique<-subset(mca_idunique,select = -mac_high.x)
names(mca_idunique)[names(mca_idunique) == "mac_high.y"] <-"mca_highrisk"
mca_idunique$mca_highrisk<-ifelse(is.na(mca_idunique$mca_highrisk),2,mca_idunique$mca_highrisk) #997 (6.79%)
save(mca_idunique, file="mca_idunique.Rdata")
#TET2 & DNMT3A
table(mca_df$type)
mca_df$tet2<-NA
mca_df$tet2<-ifelse(mca_df$chrom=="chr4" & mca_df$q_arm =="Y" & mca_df$type == "LOSS",1,mca_df$tet2)
mca_df$tet2<-ifelse(mca_df$chrom=="chr4" & mca_df$q_arm =="Y" & mca_df$type == "CN-LOH",1,mca_df$tet2) 
mca_tet2<-filter(mca_df,mca_df$tet2==1) #133
mca_df$dnmt3a<-NA
mca_df$dnmt3a<-ifelse(mca_df$chrom=="chr2" & mca_df$p_arm =="Y" & mca_df$type == "LOSS",1,mca_df$dnmt3a)
mca_df$dnmt3a<-ifelse(mca_df$chrom=="chr2" & mca_df$p_arm =="Y" & mca_df$type == "CN-LOH",1,mca_df$dnmt3a) #124
mca_dnmt3a<-filter(mca_df,mca_df$dnmt3a==1) #124
mca_idunique$tet2<-ifelse(mca_idunique$sample_id %in% mca_tet2$sample_id,1,2) #118
mca_idunique$dnmt3a<-ifelse(mca_idunique$sample_id %in% mca_dnmt3a$sample_id,1,2) #105
#add to data_all
data_all$tet2<-data_all$mca_status
data_all$dnmt3a<-data_all$mca_status
data_all$tet2<-ifelse(data_all$person_id %in% mca_tet2$sample_id,2,data_all$tet2) #116
data_all$dnmt3a<-ifelse(data_all$person_id %in% mca_dnmt3a$sample_id,2,data_all$dnmt3a) #103
table(data_all$dnmt3a)
data_all$dnmt3a<-as.factor(data_all$dnmt3a)
data_all$tet2<-as.factor(data_all$tet2)
save(data_all,file = "data_all.Rdata")
#All of us
head(dataset_86645058_survey_df, 10)
head(dataset_86645058_person_df, 10)
head(dataset_86645058_measurement_df, 10)
head(dataset_86645058_condition_df, 10)

save(dataset_86645058_survey_df,dataset_86645058_person_df,dataset_86645058_measurement_df,dataset_86645058_condition_df, file="genotyping.Rdata")
write.table(dataset_86645058_survey_df, "genotyping_survey.txt", sep="\t", row.names=F, col.names=T, quote=F)
write.table(dataset_86645058_person_df, "genotyping_person.txt", sep="\t", row.names=F, col.names=T, quote=F)
write.table(dataset_86645058_measurement_df, "genotyping_measurement.txt", sep="\t", row.names=F, col.names=T, quote=F)
write.table(dataset_86645058_condition_df, "genotyping_condition.txt", sep="\t", row.names=F, col.names=T, quote=F)
dataset_86645058_measurement_df <- fread('genotyping_measurement.txt')
dataset_86645058_survey_df <- fread('genotyping_survey.txt')
dataset_86645058_condition_df <- fread('genotyping_condition.txt')
#Ancestry data
system("gsutil -u $GOOGLE_PROJECT cp gs://fc-aou-datasets-controlled/v7/wgs/short_read/snpindel/aux/ancestry/ancestry_preds.tsv .", intern=T)
data_ancestry <- fread('ancestry_preds.tsv') 
#----2.Reshape datasets-----
table(dataset_86645058_survey_df$question)
table(dataset_86645058_person_df$race)
table(dataset_86645058_measurement_df$standard_concept_name)
data_geno<-NA
data_geno<-merge(dataset_86645058_person_df,mca_idunique,by.x = "person_id",by.y ="sample_id",all.x = T)
table(data_geno$mca)
data_geno$mca_status<-0
data_geno$mca_status<-ifelse(is.na(data_geno$mca),0,1)
table(data_geno$mca_status) #14632
#Uni question
ls(dataset_86645058_survey_df)
survey_ques<-dataset_86645058_survey_df %>% select(question) %>% distinct() %>% arrange(question) %>% mutate(number = row_number())
ls(dataset_86645058_condition_df)
condi_ques<-dataset_86645058_condition_df %>% dplyr::select(standard_concept_name) %>% distinct() %>% arrange(standard_concept_name) %>% mutate(number = dplyr::row_number())
head(dataset_86645058_measurement_df)
measure_ques <- dataset_86645058_measurement_df %>%
  dplyr::select(standard_concept_name) %>%
  distinct() %>%
  arrange(standard_concept_name) %>%
  mutate(number = dplyr::row_number())

lymohocyte_data <- dataset_86645058_measurement_df %>% filter(standard_concept_name == "Lymphocytes/100 leukocytes in Blood")
head(lymohocyte_data)
summary(lymohocyte_data$value_as_number)
#filter question number
survey_selected_no<-c(677,472,632,154,463,475,482,486,648:650)
survey_selected<-filter(survey_ques,survey_ques$number %in% survey_selected_no)
data_survey_selected <- filter(dataset_86645058_survey_df,dataset_86645058_survey_df$question %in% survey_selected$question)
#reshape to wide
data_survey_selected<-subset(data_survey_selected,select = -survey)
#keep the mini date
data_survey_selected$date <- as.POSIXct(data_survey_selected$survey_datetime, format ="%Y-%m-%d %H:%M:%S", tz = "UTC")
min_dates <- data_survey_selected %>% group_by(person_id) %>%
  summarize(min_date = min(survey_datetime)) %>%
  ungroup()
data_survey_selected <- data_survey_selected %>% left_join(min_dates, by ="person_id")
data_survey_selected <- subset(data_survey_selected,select = -survey_datetime)
data_survey_selected <- subset(data_survey_selected,select = -date)
data_survey_id<-NA
data_survey_id<- data_survey_selected %>% pivot_wider(names_from = question, values_from = answer)
#merge survey to person
data_geno2<-merge(data_geno,data_survey_id,by="person_id",all.x = T)
ls(data_geno2)
names(data_geno2)[36]<-"famhist_cancer"
names(data_geno2)[37]<-"income"
names(data_geno2)[38]<-"smoking"
names(data_geno2)[39]<-"employment"
names(data_geno2)[40]<-"education"
names(data_geno2)[41]<-"overall_health"
names(data_geno2)[42]<-"drink"
names(data_geno2)[43]<-"sex_bio"
names(data_geno2)[44]<-"overall_mental"
names(data_geno2)[45]<-"overall_phy"
names(data_geno2)[46]<-"famhist_heart_blood"
#add age
data_geno2$date_of_birth <- ymd_hms(data_geno2$date_of_birth, tz ="UTC")
data_geno2$min_date <- ymd_hms(data_geno2$min_date, tz ="UTC")
data_geno2$age <-as.numeric(difftime(data_geno2$min_date, data_geno2$date_of_birth, units ="days"))/365.25
#Adjust column
data_geno2$smoking<-gsub("100 Cigs Lifetime: ","",data_geno2$smoking)
data_geno2$smoking<-ifelse(data_geno2$smoking=="NA"|data_geno2$smoking=="PMI: Prefer Not To Answer"|data_geno2$smoking=="PMI: Skip"|data_geno2$smoking=="NULL",NA,data_geno2$smoking)
data_geno2$sex_at_birth<-ifelse(data_geno2$sex_at_birth=="I prefer not to answer"|data_geno2$sex_at_birth=="No matching concept"|data_geno2$sex_at_birth=="None"|data_geno2$sex_at_birth=="PMI: Skip",NA,data_geno2$sex_at_birth)
data_geno2$race<-ifelse(data_geno2$race=="I prefer not to answer"|data_geno2$race=="PMI: Skip",NA,data_geno2$race)
data_geno2$ethnicity<-ifelse(data_geno2$ethnicity=="PMI: Prefer Not To Answer"|data_geno2$ethnicity=="PMI: Skip"|data_geno2$ethnicity=="No matching concept",NA,data_geno2$ethnicity)
data_geno2$income<-gsub("Annual Income: ","",data_geno2$income)
data_geno2$income<-ifelse(data_geno2$income=="NA"|data_geno2$income=="PMI: Prefer Not To Answer"|data_geno2$income=="PMI: Skip",NA,data_geno2$income)
data_geno2$employment<-gsub("Employment Status: ","",data_geno2$employment)
data_geno2$employment<-ifelse(data_geno2$employment=="Employed For Wages" |data_geno2$employment=="Homemaker"|data_geno2$employment=="Out Of Work Less Than One"|data_geno2$employment=="Out Of Work One Or More"|
                                data_geno2$employment=="Retired" |data_geno2$employment=="Self Employed"|data_geno2$employment=="Student"|data_geno2$employment=="Unable To Work",data_geno2$employment,NA)
data_geno2$education<-gsub("Highest Grade: ","",data_geno2$education)
data_geno2$education<-ifelse(data_geno2$education=="NA"|data_geno2$education=="PMI: Prefer Not To Answer"|data_geno2$education=="PMI: Skip",NA,data_geno2$education)
data_geno2$drink<-gsub("6 or More Drinks Occurrence: ","",data_geno2$drink)
data_geno2$drink<-ifelse(data_geno2$drink=="NA"|data_geno2$drink=="PMI: Prefer Not To Answer"|data_geno2$drink=="PMI: Skip",NA,data_geno2$drink)
data_geno2$overall_health<-gsub("General Health: ","",data_geno2$overall_health)
data_geno2$overall_health<-ifelse(data_geno2$overall_health=="NA"|data_geno2$overall_health=="PMI: Prefer Not To Answer"|data_geno2$overall_health=="PMI: Skip"|data_geno2$overall_health=="NULL",NA,data_geno2$overall_health)
data_geno2$overall_mental<-gsub("General Mental Health: ","",data_geno2$overall_mental)
data_geno2$overall_mental<-ifelse(data_geno2$overall_mental=="NA"|data_geno2$overall_mental=="PMI: Prefer Not To Answer"|data_geno2$overall_mental=="PMI: Skip"|data_geno2$overall_mental=="NULL",NA,data_geno2$overall_mental)
data_geno2$overall_phy<-gsub("General Physical Health: ","",data_geno2$overall_phy)
data_geno2$overall_phy<-ifelse(data_geno2$overall_phy=="NA"|data_geno2$overall_phy=="PMI: Prefer Not To Answer"|data_geno2$overall_phy=="PMI: Skip"|data_geno2$overall_phy=="NULL",NA,data_geno2$overall_phy)
table(data_geno2$drink)

data_clean<-data_geno2
data_clean$cf_cat<-ifelse(is.na(data_clean$cf_cat)==T,0,data_clean$cf_cat)
data_clean$mca_highrisk<-ifelse(is.na(data_clean$mca_highrisk)==T,0,data_clean$mca_highrisk)

data_clean$income2<-"10-50k"
data_clean$income2<-ifelse(data_clean$income == "less 10k","less 10k",data_clean$income2)
data_clean$income2<-ifelse(data_clean$income == "50k 75k"|data_clean$income == "75k 100k","50-100k",data_clean$income2)
data_clean$income2<-ifelse(data_clean$income == "100k 150k"|data_clean$income == "150k 200k","100-200k",data_clean$income2)
data_clean$income2<-ifelse(data_clean$income == "more 200k","more 200k",data_clean$income2)
table(data_clean$income2,data_clean$income)
data_clean$education2<-data_clean$education
data_clean$education2<-ifelse(data_clean$education == "Five Through Eight"|data_clean$education == "Never Attended"|
                                data_clean$education == "Nine Through Eleven"|data_clean$education == "One Through Four","Less than twelve",data_clean$education2)
table(data_clean$education2,data_clean$education)

#survey_datetime
result <- data %>%
  group_by(person_id) %>%
  summarize(max_survey_datetime = max(survey_datetime)) %>%
  ungroup()
#Table 1
getTABLE1 <- function(table_name, group_name){
  library("Gmisc")
  library("magrittr")
  getT1Stat <- function(varname,digits = 2,group_var,show = 1){
    getDescriptionStatsBy(x = data_clean[,varname],by = group_var, header_count = TRUE,
                          digits = digits, add_total_col = TRUE, show_all_values = FALSE,
                          default_ref = show, hrzl_prop = FALSE, useNA = "no", statistics = TRUE,
                          statistics.sig_lim = 10^-4, html = TRUE)
  }
  ### Traits
  table <- list()
  table[["Sex"]] <- getT1Stat("sex_at_birth",2,group_name)
  table[["Age at blood draw"]] <- getT1Stat("age",2,group_name)
  table[["Race"]] <- getT1Stat("race",2,group_name)
  table[["Ethnicity"]] <- getT1Stat("ethnicity",2,group_name)
  table[["Income"]] <- getT1Stat("income2",2,group_name)
  table[["Employment"]] <- getT1Stat("employment",2,group_name)
  table[["Education"]] <- getT1Stat("education2",2,group_name)
  table[["Smoking (100packs)"]] <- getT1Stat("smoking",2,group_name)
  table[["Drinking Frequency"]] <- getT1Stat("drink",2,group_name)
  table[["General Health"]] <- getT1Stat("overall_health",2,group_name)
  table[["General Mental Health"]] <- getT1Stat("overall_mental",2,group_name)
  table[["General Physical Health"]] <- getT1Stat("overall_phy",2,group_name)
  ### Combine the results, generate html
  mergeDesc(table,htmlTable_args = list(caption = paste("Table1.Basic Characteristics by",table_name,sep = " "))) %>% htmlTable::addHtmlTableStyle(css.rgroup = "")
}
## Create Table1
getTABLE1(group_name = data_clean$mca_status,table_name = "mCA")
getTABLE1(group_name = data_clean$cf_cat,table_name = "cell fraction (10%)")
getTABLE1(group_name = data_clean$mca_highrisk,table_name = "mCA (high risk)")

#Plot_Count
table(data_clean$Count)
data_clean$count10<-ifelse(data_clean$Count>10,11,data_clean$Count)
g <- ggplot(data_clean, aes(x=count10))
p<-g + geom_bar(fill="darkred")+geom_text(stat='count', aes(label=..count..), size=6, vjust=-0.4) + theme_bw()+ labs(x="Number of mCA Mutations", y="Count", title="Number of Mutations per Individual") +theme(plot.title = element_text(hjust = 0.5))+ theme(text = element_text(size = 25))
p

#Create age_group
data_clean$agegroup<-NA
data_clean$agegroup[data_clean$age<30] <-"1"
data_clean$agegroup[data_clean$age>=30 & data_clean$age<40] <-"2"
data_clean$agegroup[data_clean$age>=40 & data_clean$age<50] <-"3"
data_clean$agegroup[data_clean$age>=50 & data_clean$age<60] <-"4"
data_clean$agegroup[data_clean$age>=60 & data_clean$age<70] <-"5"
data_clean$agegroup[data_clean$age>=70] <-"6"
table(data_clean$agegroup,data_clean$mca_status)

system("gsutil cp data_clean.Rdata gs://fc-secure-cb192ac6-30ba-46b9-92ee-896a6e36c63e/kun/")
system("gsutil cp data_geno2.Rdata gs://fc-secure-cb192ac6-30ba-46b9-92ee-896a6e36c63e/kun/")
system("gsutil cp mca_idunique.Rdata gs://fc-secure-cb192ac6-30ba-46b9-92ee-896a6e36c63e/kun/")
system("gsutil cp genotyping.Rdata gs://fc-secure-cb192ac6-30ba-46b9-92ee-896a6e36c63e/kun/")
write.table(dataset_86645058_condition_df, "data_condition.txt", sep="\t", row.names=F, col.names=T, quote=F)
system("gsutil cp data_condition.txt gs://fc-secure-cb192ac6-30ba-46b9-92ee-896a6e36c63e/kun/")
system("gsutil cp data_condition.Rdata gs://fc-secure-cb192ac6-30ba-46b9-92ee-896a6e36c63e/kun/")

# Plot 2 by Race
filtered_race <- data_clean %>% filter(race %in% c("White", "Asian", "Black or African American"))
filtered_race <- filtered_race %>%
  mutate(agegroup = factor(agegroup, levels = 1:6, labels = c("<=30", "30-40", "40-50", "50-60", "60-70", ">70")))

rate_data <- filtered_race %>% group_by(agegroup, race) %>% summarize(mca_rate = mean(mca_status)) %>% ungroup()


ggplot(rate_data, aes(x = agegroup, y = mca_rate, color = race, group = race)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "Proportion of individuals with detectable mCA by age group",
    x = "Age group",
    y = "Proportion of individuals with mCA"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


describeBy(data_clean$age,data_clean$race)

#Plot 2.1 Race+ CF group
filtered_race$cf_cat<-ifelse(is.na(filtered_race$cf_cat),0,filtered_race$cf_cat)
filtered_race <- filtered_race %>%
  mutate(cf_cat = factor(cf_cat, levels = c(1, 2), labels = c("small", "large (>=10%)")))

df_summary <- filtered_race %>%
  group_by(agegroup, race, cf_cat) %>%
  summarize(mCA_rate = n()) %>%
  mutate(mCA_rate = mCA_rate / sum(mCA_rate)) %>%
  ungroup()
df_summary <- df_summary %>%
  filter(!is.na(cf_cat))

ggplot(df_summary, aes(x = agegroup, y = mCA_rate, color = race, linetype = cf_cat, group = interaction(race, cf_cat))) +
  geom_line(size = 1) + scale_y_continuous(labels = percent) +
  labs(
    title = "Proportion of individuals with detectable mCA by age group",
    x = "Age group",
    y = "Proportion of individuals with mCA",
    color = "Race",
    linetype = "Cell fraction"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right"
  )

#The Latest Survey Time
dataset_86645058_survey_df$survey_datetime <- as.POSIXct(dataset_86645058_survey_df$survey_datetime, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
max_dates <- data_condition[,c("person_id","max_date")]
data_clean<- merge(data_clean,max_dates,by="person_id",all.x = T)
save(data_clean, file="data_clean.Rdata")
#Adjusted covariates
data_condition2$age_squ <- data_condition2$age * data_condition2$age

#----3.Define outcomes----
system("gsutil cp gs://fc-secure-cb192ac6-30ba-46b9-92ee-896a6e36c63e/pershy1/tls_Icd10cmHumanReadableMap_US1000124_20240301.tsv .", intern = TRUE)
icd10code <- fread('tls_Icd10cmHumanReadableMap_US1000124_20240301.tsv') 
ls(icd10code)
icd10code2<-icd10code[,c("referencedComponentId","mapTarget")]
names(icd10code2)[2]<-"icd10"
icd10code <- filter(icd10code,icd10code$mapTarget!="")
icd10code <- icd10code %>% mutate(referencedComponentId = as.integer(referencedComponentId))
icd10code2 <- icd10code2 %>% mutate(referencedComponentId = as.integer(referencedComponentId))
icd10code2 <- icd10code2 %>% distinct(referencedComponentId, .keep_all = TRUE)
data_condition10<-merge(dataset_86645058_condition_df,icd10code2,by.x = "standard_concept_code", by.y = "referencedComponentId",all.x = T)
data_condition10$icd10<-ifelse(is.na(data_condition10$source_vocabulary) == F & data_condition10$source_vocabulary == "ICD10CM", data_condition10$source_concept_code,data_condition10$icd10)
test<-filter(data_condition10,is.na(data_condition10$icd10))
head(data_condition10)
save(data_condition10,file="data_condition10.Rdata")
system("gsutil cp data_condition10.Rdata gs://fc-secure-cb192ac6-30ba-46b9-92ee-896a6e36c63e/kun/")

#heart failure
condi_heartfail_icd10 <- data_condition10 %>%
  filter(grepl("I50", icd10, ignore.case = TRUE) | grepl("P29.0", icd10, ignore.case = TRUE) |grepl("I11.0", icd10, ignore.case = TRUE)) %>%
  arrange(person_id, condition_start_datetime) %>%
  distinct(person_id, .keep_all = TRUE) %>% 
  mutate(heart_failure = 1) %>%  
  select(person_id, condition_start_datetime, heart_failure, standard_concept_name)
names(condi_heartfail_icd10)[4]<-"type_HF"
names(condi_heartfail_icd10)[2]<-"start_date_HF"
data_condition2<- merge(data_clean,condi_heartfail_icd10,by="person_id",all.x = T)
#----3.Analysis (Single Disease)------
#------3.1. CLL chronic lymphocytic leukemia------
condi_cll <- NA
condi_cll <- data_condition10 %>%
  filter(grepl("C91.1", icd10, ignore.case = TRUE)) %>%
  arrange(person_id, condition_start_datetime) %>%
  distinct(person_id, .keep_all = TRUE) %>% 
  mutate(cll = 1) %>%  
  select(person_id, condition_start_datetime, cll, standard_concept_name)
names(condi_cll)[4]<-"type_CLL"
names(condi_cll)[2]<-"start_date_cll"

#CLL survival 
data_condition2<- merge(data_condition2,condi_cll,by="person_id",all.x = T)
data_condition2$start_date_cll <- ymd_hms(data_condition2$start_date_cll, tz ="UTC")
data_condition2$survial_cll <- NA
data_condition2$survial_cll <-as.numeric(difftime(data_condition2$start_date_cll, data_condition2$min_date, units ="days"))/365.25
summary(data_condition2$survial_cll)
data_condition2$survial_cll <- ifelse(is.na(data_condition2$survial_cll),data_condition2$survial_all,data_condition2$survial_cll)
summary(data_condition2$survial_cll)
data_condition2$cll<- ifelse(is.na(data_condition2$cll),0,data_condition2$cll)

#K-M plot
df_filtered <- data_condition2 %>% filter(survial_cll > 0)
surv_obj <- Surv(time = df_filtered$survial_cll, event = df_filtered$cll == 1)
fit <- survfit(surv_obj ~ mca_status, data = df_filtered)
ggsurv<-ggsurvplot(
  fit, 
  data = df_filtered,
  fun = "event",
  pval = TRUE,
  risk.table = TRUE,
  conf.int = TRUE,
  xlab = "Time (Years)",
  ylab = "Probability of CLL",
  legend.title = "mCA",
  legend.labs = c("Non-Carriers", "mCA Carriers"),
  palette = c("#2E9FDF", "#E7B800"),
  ggtheme = theme_minimal()
)
ggsurv

plot_mca_kmcat(data_condition2, "survial_cll", "cll", "cf_cat", "CLL")
plot_mca_km_cll(data_condition2, "survial_cll", "cll", "mca_highrisk", "CLL")

#-*Cox Analysis formula----------
analyze_cox_model <- function(data, time_var, event_var, group_var, covariates) {
  df_filtered <- data %>% filter(.data[[time_var]] > 0)
  df_filtered <- df_filtered %>% drop_na(all_of(covariates))
  surv_obj <- Surv(time = df_filtered[[time_var]], event = df_filtered[[event_var]] == 1)
  covariate_formula <- paste(covariates, collapse = " + ")
  formula <- as.formula(paste("surv_obj ~", group_var, "+", covariate_formula))
  cox_fit <- coxph(formula, data = df_filtered)
  summary_cox <- summary(cox_fit)
  print(summary_cox)
  results <- summary_cox$coefficients
  colnames(results) <- c("Coefficient", "Exp(Coefficient)", "Standard Error", "z value", "p value")
  print(results)
}
data_condition$mca_highrisk<-factor(data_condition, levels = c(0, 1, 2), labels = c("non_carriers", "high_risk", "low_risk"))
data_condition$mca_highrisk<-as.factor(data_condition$mca_highrisk)
data_condition$cf_cat<-as.factor(data_condition$cf_cat)
analyze_cox_model(data_condition, "survial_cll", "cll", "mca_status", c("age", "age_squ", "sex_at_birth", "smoking", "race"))
analyze_cox_model(data_condition, "survial_cll", "cll", "mca_highrisk", c("age", "age_squ", "sex_at_birth", "smoking", "race"))
analyze_cox_model(data_condition, "survial_cll", "cll", "cf_cat", c("age", "age_squ", "sex_at_birth", "smoking", "race"))
#CHIP
analyze_cox_model(data_all_wgs, "survial_hf", "heart_failure", "chip", c("age", "age_squ", "sex_at_birth", "smoking", "race"))

plot_cox_model <- function(data, time_var, event_var, group_vars, covariates) {
  df_filtered <- data %>% filter(.data[[time_var]] > 0)
  df_filtered <- df_filtered %>% drop_na(all_of(covariates))
  surv_obj <- Surv(time = df_filtered[[time_var]], event = df_filtered[[event_var]] == 1)
  
  results_list <- list()
  row_labels <- list()
  
  for (group_var in group_vars) {
    df_filtered[[group_var]] <- factor(df_filtered[[group_var]], levels = unique(df_filtered[[group_var]]))
    covariate_formula <- paste(covariates, collapse = " + ")
    formula <- as.formula(paste("surv_obj ~", group_var, "+", covariate_formula))
    cox_fit <- coxph(formula, data = df_filtered)
    summary_cox <- summary(cox_fit)
    print(summary_cox)
    coef_table <- summary_cox$coefficients
    conf_int <- summary_cox$conf.int
    results <- cbind(coef_table, conf_int)
    results <- results[grep(group_var, rownames(results)), , drop = FALSE]
    print(results)
    results_list[[group_var]] <- results
    row_labels[[group_var]] <- rownames(results)
  }

  final_results <- do.call(rbind, results_list)
  

  tabletext <- rbind(
    c("Traits", "HR", "95% CI", "P value"),
    cbind(
      Variable = unlist(row_labels),
      HR = sprintf("%.2f", final_results[, "exp(coef)"]),
      `95% CI` = paste0(sprintf("%.2f", final_results[, "lower .95"]), " - ", sprintf("%.2f", final_results[, "upper .95"])),
      `p value` = format(final_results[, "Pr(>|z|)"], scientific = F)
    )
  )
  

  mean_values <- c(NA, final_results[, "exp(coef)"])
  lower_values <- c(NA, final_results[, "lower .95"])
  upper_values <- c(NA, final_results[, "upper .95"])
  
  forestplot(labeltext = tabletext, 
             mean = mean_values, 
             lower = lower_values, 
             upper = upper_values,
             zero = 1, 
             boxsize = 0.2,
             lineheight = unit(8, "mm"),
             col = fpColors(box = "royalblue", line = "darkblue", summary = "royalblue"),
             xlab = "Hazard Ratio",
             title = "Forest Plot for Cox Regression Results")
}
plot_cox_model(data_condition2, "survial_cll", "cll", c("mca_status", "cf_cat", "mca_highrisk"), c("age", "age_squ", "sex_at_birth", "smoking", "race"))
plot_cox_model(data_condition2, "survial_hf", "heart_failure", c("mca_status", "cf_cat"), c("age", "age_squ", "sex_at_birth", "smoking", "race"))

#----4.Plot (chr15 Gain)-------
data_age <- data_all[,c("person_id","age","race")]
data_cfmax_unique <- merge(data_cfmax_unique,data_age,by.x = "sample_id",by.y = "person_id", all.x = T)
ggplot(data_cfmax_unique, aes(x = age)) +
  geom_histogram(binwidth = 1, fill = "skyblue", color = "black", alpha = 0.7) +
  theme_minimal() +
  labs(title = "Age Distribution in mCA carriers",
       x = "Age",
       y = "Frequency")


data_filtered <- merge(data_cfmax,data_age,by.x = "sample_id",by.y = "person_id", all.x = T)
ls(data_filtered)
table(data_filtered$q_arm,data_filtered$p_arm)
data_filtered$p_count<-NA
data_filtered$q_count<-NA
data_filtered$p_count<-ifelse(data_filtered$p_arm=="Y"|data_filtered$p_arm=="T"|data_filtered$p_arm=="C",1,0)
data_filtered$q_count<-ifelse(data_filtered$q_arm=="Y"|data_filtered$q_arm=="T"|data_filtered$q_arm=="C",1,0)
table(data_filtered$p_count,data_filtered$q_count)
data_filtered$arm<-NA
data_filtered$arm<-ifelse(data_filtered$p_count ==1 & data_filtered$q_count ==0 ,"p",data_filtered$arm)
data_filtered$arm<-ifelse(data_filtered$q_count ==1 & data_filtered$p_count ==0 ,"q",data_filtered$arm)
table(data_filtered$arm)
data_filtered$chr_type <- paste(data_filtered$chrom,data_filtered$arm,sep = "")
data_filtered <- data_filtered %>% mutate(chr_type = str_replace(chr_type, "NA$", ""))
data_filtered$chr15gain<-NA
data_filtered$chr15gain<-ifelse(data_filtered$chr_type=="chr15" & data_filtered$type=="Gain",1,2)
#save(data_filtered,file = "data_filtered.Rdata")
data_filtered<- merge(data_filtered,data_ancestry,by.x = "sample_id", by.y = "research_id",all.x = T)
table(data_filtered$ancestry_pred,useNA = "always")

data_filtered <- data_filtered %>% filter(type != "Undetermined")
table(data_filtered$type)
data_summary<-NA
data_summary <- data_filtered %>%
  mutate(sex_numeric = ifelse(computed_gender == "M", 1, 0)) %>%
  group_by(type, chr_type) %>%
  summarise(
    mean_age = mean(age, na.rm = TRUE),
    sem_age = sd(age, na.rm = TRUE) / sqrt(n()),
    mean_fraction_male = mean(sex_numeric, na.rm = TRUE),
    sem_fraction_male = sd(sex_numeric, na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = 'drop'
  )
data_summary<-filter(data_summary,data_summary$mean_fraction_male>0 & data_summary$mean_fraction_male<1)
data_summary<-filter(data_summary,data_summary$n>=50)

data_summary <- data_summary %>%
  mutate(show_label = ifelse(chr_type %in% c("chr15", "chr2","chr15q","chr15p"), TRUE, FALSE))

ggplot(data_summary, aes(x = mean_age, y = mean_fraction_male, color = type, alpha = n)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_fraction_male - sem_fraction_male, ymax = mean_fraction_male + sem_fraction_male), width = 0.2) +
  geom_errorbarh(aes(xmin = mean_age - sem_age, xmax = mean_age + sem_age), height = 0.2) +
  geom_text(data = filter(data_summary, show_label), aes(label = chr_type), hjust = 1.1, vjust = 1.1) +
  scale_color_manual(values = c("Gain" = "blue", "CN-LOH" = "lightgreen", "Loss" = "red")) +
  scale_alpha_continuous(range = c(0.5, 1)) +
  labs(title = "Carriers of different mCA types have different age and sex distributions",
       x = "Mean age (years)",
       y = "Fraction male",
       alpha = "Sample Size") +
  theme_minimal() +
  theme(legend.position = "right",
        legend.title = element_blank())
#Race plot
data_summary2<-NA
data_summary2 <- data_filtered %>%
  mutate(sex_numeric = ifelse(ancestry_pred == "afr", 1, 0)) %>%
  group_by(type, chr_type) %>%
  summarise(
    mean_age = mean(age, na.rm = TRUE),
    sem_age = sd(age, na.rm = TRUE) / sqrt(n()),
    mean_fraction_eur = mean(sex_numeric, na.rm = TRUE),
    sem_fraction_eur = sd(sex_numeric, na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = 'drop'
  )
data_summary2<-filter(data_summary2,data_summary2$mean_fraction_eur>0 & data_summary2$mean_fraction_eur<1)
data_summary2<-filter(data_summary2,data_summary2$n>=50)
data_summary2 <- data_summary2 %>%
  mutate(show_label = ifelse(chr_type %in% c("chr15", "chr2","chr15q"), TRUE, FALSE))

ggplot(data_summary2, aes(x = mean_age, y = mean_fraction_eur, color = type, alpha = n)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_fraction_eur - sem_fraction_eur, ymax = mean_fraction_eur + sem_fraction_eur), width = 0.2) +
  geom_errorbarh(aes(xmin = mean_age - sem_age, xmax = mean_age + sem_age), height = 0.2) +
  geom_text(data = filter(data_summary2, show_label), aes(label = chr_type), hjust = 1.1, vjust = 1.1) +
  scale_color_manual(values = c("Gain" = "blue", "CN-LOH" = "lightgreen", "Loss" = "red")) +
  scale_alpha_continuous(range = c(0.5, 1)) +
  labs(title = "Carriers of different mCA types have different age and ancestry distributions",
       x = "Mean age (years)",
       y = "Fraction EUR",
       alpha = "Sample Size") +
  theme_minimal() +
  theme(legend.position = "right",
        legend.title = element_blank())
#Association between chr15gain and smoking
load("data_filtered")
ls(data_filtered)
data_chr15g <- data_filtered[,c("sample_id","chr15gain")]
data_chr15g_unique <- data_chr15g %>%
  group_by(sample_id) %>%
  summarise(chr15gain = ifelse(any(chr15gain == 1), 1, 0), .groups = 'drop')
table(data_chr15g_unique$chr15gain) #123
data_all <- merge(data_all,data_chr15g_unique,by.x = "person_id",by.y = "sample_id",all.x = T)
save(data_all,file = "data_all.Rdata")
table(data_all$chr15gain)
data_all$chr15gain<-ifelse(data_all$chr15gain=="0",2,data_all$chr15gain)
data_all$chr15gain<-ifelse(is.na(data_all$chr15gain),0,data_all$chr15gain)
data_all$chr15gain<-as.factor(data_all$chr15gain)

#chr15gain vs mCA(non-chr15gain)
analysis_15gain<-filter(data_all,data_all$chr15gain != 0)
analysis_15gain <- analysis_15gain %>%
  mutate(chr15gain_ref0 = relevel(chr15gain, ref = "2"))

model <- glm(chr15gain_ref0 ~ smoking + age + age_squ + sex_at_birth + race, data = analysis_15gain, family = binomial)
summary(model)
install.packages("gtsummary")
library(gtsummary)
tbl_regression(model, exponentiate = TRUE)
#chr15gain vs non-mCA
analysis_15gain<-filter(data_all,data_all$chr15gain != 2)
model <- glm(chr15gain ~ smoking + age + age_squ + sex_at_birth + race, data = analysis_15gain, family = binomial)
summary(model)
tbl_regression(model, exponentiate = TRUE)
# chr 15 : 20,500,000 to 31,000,000
data_15gain<-filter(data_filtered,data_filtered$chr15gain == "1")

#Age Plot
data_chr15 <- filter(data_filtered,data_filtered$chr_type == "chr15" & data_filtered$type == "Gain")
ggplot(data_chr15, aes(x = age)) +
  geom_histogram(binwidth = 1, fill = "skyblue", color = "black", alpha = 0.7) +
  theme_minimal() +
  labs(title = "Age Distribution in data_chr15",
       x = "Age",
       y = "Frequency")

data_chr15q <- filter(data_filtered,data_filtered$chr_type == "chr15q" & data_filtered$type == "Gain")
ggplot(data_chr15q, aes(x = age)) +
  geom_histogram(binwidth = 1, fill = "skyblue", color = "black", alpha = 0.7) +
  theme_minimal() +
  labs(title = "Age Distribution in data_chr15q",
       x = "Age",
       y = "Frequency")
summary(data_chr15q$age)

#Combine chr15 and chr15q
data_filtered$chr_type2<-data_filtered$chr_type
data_filtered$chr_type2<- ifelse(data_filtered$chr_type2=="chr15q","chr15",data_filtered$chr_type2)



#-TOPMed------
data_topmed <- read.csv(file = "long.mocha.filtered.csv" , sep = "\t")
data_topmed_age <- fread('topmed_tcell_cov.txt') 
data_topmed_ances <- fread('Workspaces_bicklab-tcr-quant_Hannahs_Folder_all_ancestries.txt') 
data_topmed_ances <- data_topmed_ances[,c("sample","ancestry")]
data_topmed_all <- merge(data_topmed,data_topmed_age, by.x = "sample_id",by.y = "sample",all.x = T)
summary(data_topmed_all$AgeAtBloodDraw)
data_topmed_all <- merge(data_topmed_all,data_topmed_ances, by.x = "sample_id",by.y = "sample",all.x = T)
table(data_topmed_all$ancestry)
save(data_topmed_all, file = "data_topmed_all.Rdata")

#------Combine TOPMed and AoU------
data_AoU <- data_filtered[,c(1:22,28,30:33,35)]
names(data_AoU)[28]<-"ancestry"
data_AoU <- data_AoU %>% mutate(ancestry = str_to_upper(ancestry))

data_topmed_all$computed_gender <-data_topmed_all$sex

data_topmed_all$p_count<-NA
data_topmed_all$q_count<-NA
data_topmed_all$p_count<-ifelse(data_topmed_all$p_arm=="Y"|data_topmed_all$p_arm=="T"|data_topmed_all$p_arm=="C",1,0)
data_topmed_all$q_count<-ifelse(data_topmed_all$q_arm=="Y"|data_topmed_all$q_arm=="T"|data_topmed_all$q_arm=="C",1,0)
table(data_topmed_all$p_count,data_topmed_all$q_count)
data_topmed_all$arm<-NA
data_topmed_all$arm<-ifelse(data_topmed_all$p_count ==1 & data_topmed_all$q_count ==0 ,"p",data_topmed_all$arm)
data_topmed_all$arm<-ifelse(data_topmed_all$q_count ==1 & data_topmed_all$p_count ==0 ,"q",data_topmed_all$arm)
table(data_topmed_all$arm)
data_topmed_all$chr_type <- paste(data_topmed_all$chrom,data_topmed_all$arm,sep = "")
data_topmed_all <- data_topmed_all %>% mutate(chr_type = str_replace(chr_type, "NA$", ""))
names(data_topmed_all)[24] <- "age"
data_topmed_all<- data_topmed_all[,c(1:22,24:29)]
ls(data_AoU)
data_topmed_all$cohort<-"TOPMed"
data_AoU$cohort<-"AoU"
data_TOP_AoU <- merge(data_topmed_all,data_AoU,all = T)
save(data_TOP_AoU,file = "data_TOP_AoU.Rdata")
table(data_TOP_AoU$cohort)
#------TOPMed only--------
table(data_topmed_all$type)
data_filtered <- data_topmed_all %>% filter(type != "Undetermined")
table(data_filtered$type)
data_summary<-NA
data_summary <- data_filtered %>%
  mutate(sex_numeric = ifelse(computed_gender == "M", 1, 0)) %>%
  group_by(type, chr_type) %>%
  summarise(
    mean_age = mean(age, na.rm = TRUE),
    sem_age = sd(age, na.rm = TRUE) / sqrt(n()),
    mean_fraction_male = mean(sex_numeric, na.rm = TRUE),
    sem_fraction_male = sd(sex_numeric, na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = 'drop'
  )
data_summary<-filter(data_summary,data_summary$mean_fraction_male>0 & data_summary$mean_fraction_male<1)
data_summary<-filter(data_summary,data_summary$n>=20)

data_summary <- data_summary %>%
  mutate(show_label = ifelse(chr_type %in% c("chr15", "chr2","chr15q","chr15p"), TRUE, FALSE))

ggplot(data_summary, aes(x = mean_age, y = mean_fraction_male, color = type, alpha = n)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_fraction_male - sem_fraction_male, ymax = mean_fraction_male + sem_fraction_male), width = 0.2) +
  geom_errorbarh(aes(xmin = mean_age - sem_age, xmax = mean_age + sem_age), height = 0.2) +
  geom_text(data = filter(data_summary, show_label), aes(label = chr_type), hjust = 1.1, vjust = 1.1) +
  scale_color_manual(values = c("Gain" = "blue", "CN-LOH" = "lightgreen", "Loss" = "red")) +
  scale_alpha_continuous(range = c(0.5, 1)) +
  labs(title = "Carriers of different mCA types have different age and sex distributions",
       x = "Mean age (years)",
       y = "Fraction male",
       alpha = "Sample Size") +
  theme_minimal() +
  theme(legend.position = "right",
        legend.title = element_blank())
#------TOPMed & AoU--------
data_TOP_AoU$chr_type2<-data_TOP_AoU$chr_type
data_TOP_AoU$chr_type2<- ifelse(data_TOP_AoU$chr_type2=="chr15q","chr15",data_TOP_AoU$chr_type2)

data_filtered <- data_TOP_AoU %>% filter(type != "Undetermined")
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
data_summary<-filter(data_summary,data_summary$n>=50)

data_summary <- data_summary %>%
  mutate(show_label = ifelse(chr_type2 %in% c("chr15", "chr2","chr15q","chr15p"), TRUE, FALSE))

data_summary<-filter(data_summary,data_summary$chr_type2 != "chrX" & data_summary$chr_type2 != "chrXp" & data_summary$chr_type2 != "chrXq")

ggplot(data_summary, aes(x = mean_age, y = mean_fraction_male, color = type, alpha = n)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_fraction_male - sem_fraction_male, ymax = mean_fraction_male + sem_fraction_male), width = 0.2) +
  geom_errorbarh(aes(xmin = mean_age - sem_age, xmax = mean_age + sem_age), height = 0.2) +
  geom_text(data = filter(data_summary, show_label), aes(label = chr_type2), hjust = 1.1, vjust = 1.1) +
  scale_color_manual(values = c("Gain" = "blue", "CN-LOH" = "lightgreen", "Loss" = "red")) +
  scale_alpha_continuous(range = c(0.5, 1)) +
  labs(title = "Carriers of different mCA types have different age and sex distributions",
       x = "Mean age (years)",
       y = "Fraction male",
       alpha = "Sample Size") +
  theme_minimal() +
  theme(legend.position = "right",
        legend.title = element_blank())

#--------Ancestry plot
# 计算每种ancestry中的每种chr_type的计数
ls(data_TOP_AoU)
counts <- data_TOP_AoU %>%
  group_by(ancestry, chr_type2,type) %>%
  summarise(count = n()) %>%
  ungroup()

# 计算每种ancestry的总计数
total_counts <- data_TOP_AoU %>%
  group_by(ancestry) %>%
  summarise(total = n()) %>%
  ungroup()

# 将总计数合并到计数数据框中
counts <- counts %>%
  left_join(total_counts, by = "ancestry") %>%
  mutate(proportion = (count / total) * 100)  # 计算比例，转化为百分比

# 筛选EUR和AFR的比例并合并到一个数据框中
combined_data <- filter(counts,counts$ancestry=="EUR" | counts$ancestry=="AFR" )
combined_data <- filter(combined_data, combined_data$type!= "Undetermined" )
combined_data <- filter(combined_data, combined_data$chr_type2!= "chrX" & combined_data$chr_type2!= "chrXp" & combined_data$chr_type2!= "chrXq")

eur_data <- combined_data %>%
  filter(ancestry == "EUR") %>%
  select(chr_type2,type, eur_proportion = proportion, count_eur = count)

afr_data <- combined_data %>%
  filter(ancestry == "AFR") %>%
  select(chr_type2,type, afr_proportion = proportion, count_afr = count)

combined_data <- merge(eur_data, afr_data, by = c("chr_type2","type"))
combined_data <- filter(combined_data,combined_data$count_eur >= 20 | combined_data$count_afr >= 20)
combined_data$count <- combined_data$count_eur + combined_data$count_afr


ggplot(combined_data, aes(x = eur_proportion, y = afr_proportion, color = type, label = chr_type2, size = count)) +
  geom_point() +  # 绘制散点，大小根据count变量
  geom_text(vjust = -0.5, hjust = 0.5) +  # 在每个点旁边显示chr_type
  scale_x_continuous(name = "EUR Proportion (%)", limits = c(0, 2)) +
  scale_y_continuous(name = "AFR Proportion (%)", limits = c(0, 2)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey") + 
  scale_color_manual(values = c("Gain" = "#56B4E9", "Loss" = "#E69F00", "CN-LOH" = "#009E73")) +  # 定义颜色
  scale_size_continuous(range = c(1, 5), name = "Count") +  # 调整点的大小范围
  theme_minimal() +
  labs(title = "Chromosome Type Proportion in EUR vs AFR",
       subtitle = "Colored by Chromosome Event Type")

#--------Gender plot
# Gender calculate
table(data_TOP_AoU$computed_gender, useNA="always")
counts <- data_TOP_AoU %>%
  group_by(computed_gender, chr_type2,type) %>%
  summarise(count = n()) %>%
  ungroup()

# 计算每种ancestry的总计数
total_counts <- data_TOP_AoU %>%
  group_by(computed_gender) %>%
  summarise(total = n()) %>%
  ungroup()

# 将总计数合并到计数数据框中
counts <- counts %>%
  left_join(total_counts, by = "computed_gender") %>%
  mutate(proportion = (count / total) * 100)  # 计算比例，转化为百分比

# 筛选EUR和AFR的比例并合并到一个数据框中
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
  geom_point() +  # 绘制散点，大小根据count变量
  geom_text(vjust = -0.5, hjust = 0.5) +  # 在每个点旁边显示chr_type
  scale_x_continuous(name = "Male Proportion (%)", limits = c(0, 2)) +
  scale_y_continuous(name = "Female Proportion (%)", limits = c(0, 2)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey") + 
  scale_color_manual(values = c("Gain" = "#56B4E9", "Loss" = "#E69F00", "CN-LOH" = "#009E73")) +  # 定义颜色
  scale_size_continuous(range = c(1, 5), name = "Count") +  # 调整点的大小范围
  theme_minimal() +
  labs(title = "Chromosome Type Proportion in Male vs Female",
       subtitle = "Colored by Chromosome Event Type")

#-BioVU------
data_biovu <- fread('biovu_mca_calls..txt') 
data_biovu <- data_biovu[,c(1:22,48)]
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
data_biovu$ancestry<-NA
#ancestry of biovu
biovu_AFR <- fread("20200515_biallelic_mega_recalled.chr1-22.grid.AA.filt1.fam",header = FALSE)
biovu_EUR <- fread("20200518_biallelic_mega_recalled.chr1-22.grid.EU.filt1.fam",header = FALSE)
data_biovu <- data_biovu %>% mutate(ancestry = ifelse(sample_id %in% biovu_AFR$V2, "AFR", ancestry))
data_biovu <- data_biovu %>% mutate(ancestry = ifelse(sample_id %in% biovu_EUR$V2, "EUR", ancestry))
table(data_biovu$ancestry,useNA = "always")
#------BioVU only--------
table(data_biovu$type)
data_filtered <- data_biovu %>% filter(type != "Undetermined")
table(data_filtered$type)
data_summary<-NA
data_summary <- data_filtered %>%
  mutate(sex_numeric = ifelse(computed_gender == "M", 1, 0)) %>%
  group_by(type, chr_type) %>%
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
data_summary<-filter(data_summary,data_summary$chr_type != "chrX" & data_summary$chr_type != "chrXp" & data_summary$chr_type != "chrXq")

data_summary <- data_summary %>%
  mutate(show_label = ifelse(chr_type %in% c("chr15", "chr2","chr15q","chr15p"), TRUE, FALSE))

ggplot(data_summary, aes(x = mean_age, y = mean_fraction_male, color = type, alpha = n)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_fraction_male - sem_fraction_male, ymax = mean_fraction_male + sem_fraction_male), width = 0.2) +
  geom_errorbarh(aes(xmin = mean_age - sem_age, xmax = mean_age + sem_age), height = 0.2) +
  geom_text(data = filter(data_summary, show_label), aes(label = chr_type), hjust = 1.1, vjust = 1.1) +
  scale_color_manual(values = c("Gain" = "blue", "CN-LOH" = "lightgreen", "Loss" = "red")) +
  scale_alpha_continuous(range = c(0.5, 1)) +
  labs(title = "Carriers of different mCA types have different age and sex distributions",
       x = "Mean age (years)",
       y = "Fraction male",
       alpha = "Sample Size") +
  theme_minimal() +
  theme(legend.position = "right",
        legend.title = element_blank())
#------Combine BioVU, TOPMed and AoU------
data_biovu$cohort<-"BioVU"
data_biovu$ancestry<-NA
data_biovu$chr_type2<-data_biovu$chr_type
data_biovu$chr_type2<- ifelse(data_biovu$chr_type2=="chr15q","chr15",data_biovu$chr_type2)

data_TOP_AoU_bioVU <- merge(data_TOP_AoU,data_biovu,all = T)
save(data_TOP_AoU_bioVU,file = "data_TOP_AoU_bioVU.Rdata")
table(data_TOP_AoU_bioVU$cohort)
system("gsutil cp data_TOP_AoU_bioVU.Rdata gs://fc-secure-cb192ac6-30ba-46b9-92ee-896a6e36c63e/kun/")

#------BioVU, TOPMed & AoU--------
data_filtered <- data_TOP_AoU_bioVU %>% filter(type != "Undetermined")
data_filtered <- filter(data_filtered,age>60 & age <90)
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
data_summary<-filter(data_summary,data_summary$n>=60)


data_summary <- data_summary %>%
  mutate(show_label = ifelse(chr_type2 %in% c("chr15","chr20q","chr12"), TRUE, FALSE))

data_summary<-filter(data_summary,data_summary$chr_type2 != "chrX" & data_summary$chr_type2 != "chrXp" & data_summary$chr_type2 != "chrXq")

ggplot(data_summary, aes(x = mean_age, y = mean_fraction_male, color = type, alpha = n)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_fraction_male - sem_fraction_male, ymax = mean_fraction_male + sem_fraction_male), width = 0.2) +
  geom_errorbarh(aes(xmin = mean_age - sem_age, xmax = mean_age + sem_age), height = 0.2) +
  geom_text(data = filter(data_summary, show_label), aes(label = chr_type2), hjust = 1.1, vjust = 1.1) +
  scale_color_manual(values = c("Gain" = "blue", "CN-LOH" = "lightgreen", "Loss" = "red")) +
  scale_alpha_continuous(range = c(0.5, 1)) +
  labs(title = "Carriers of different mCA types have different age and sex distributions",
       x = "Mean age (years)",
       y = "Fraction male",
       alpha = "Sample Size") +
  theme_minimal() +
  theme(legend.position = "right",
        legend.title = element_blank())

#--------Ancestry plot
# 计算每种ancestry中的每种chr_type的计数
ls(data_TOP_AoU_bioVU)
data_filtered <- filter(data_TOP_AoU_bioVU,age>40 & age<90)
counts <- data_filtered %>%
  group_by(ancestry, chr_type2,type) %>%
  summarise(count = n()) %>%
  ungroup()

# 计算每种ancestry的总计数
total_counts <- data_filtered %>%
  group_by(ancestry) %>%
  summarise(total = n()) %>%
  ungroup()

# 将总计数合并到计数数据框中
counts <- counts %>%
  left_join(total_counts, by = "ancestry") %>%
  mutate(proportion = (count / total) * 100)  # 计算比例，转化为百分比

# 筛选EUR和AFR的比例并合并到一个数据框中
combined_data <- filter(counts,counts$ancestry=="EUR" | counts$ancestry=="AFR" )
combined_data <- filter(combined_data, combined_data$type!= "Undetermined" )
combined_data <- filter(combined_data, combined_data$chr_type2!= "chrX" & combined_data$chr_type2!= "chrXp" & combined_data$chr_type2!= "chrXq")

eur_data <- combined_data %>%
  filter(ancestry == "EUR") %>%
  select(chr_type2,type, eur_proportion = proportion, count_eur = count)

afr_data <- combined_data %>%
  filter(ancestry == "AFR") %>%
  select(chr_type2,type, afr_proportion = proportion, count_afr = count)

combined_data <- merge(eur_data, afr_data, by = c("chr_type2","type"))
combined_data <- filter(combined_data,combined_data$count_eur >= 20 | combined_data$count_afr >= 20)
combined_data$count <- combined_data$count_eur + combined_data$count_afr


ggplot(combined_data, aes(x = eur_proportion, y = afr_proportion, color = type, label = chr_type2, size = count)) +
  geom_point() +  # 绘制散点，大小根据count变量
  geom_text(vjust = -0.5, hjust = 0.5) +  # 在每个点旁边显示chr_type
  scale_x_continuous(name = "EUR Proportion (%)", limits = c(0, 2.2)) +
  scale_y_continuous(name = "AFR Proportion (%)", limits = c(0, 2.2)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey") + 
  scale_color_manual(values = c("Gain" = "#56B4E9", "Loss" = "#E69F00", "CN-LOH" = "#009E73")) +  # 定义颜色
  scale_size_continuous(range = c(1, 5), name = "Count") +  # 调整点的大小范围
  theme_minimal() +
  labs(title = "Chromosome Type Proportion in EUR vs AFR",
       subtitle = "Colored by Chromosome Event Type")

#--------Gender plot
# Gender calculate
data_filtered <- filter(data_TOP_AoU_bioVU,age>40 & age<90)
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

# 筛选EUR和AFR的比例并合并到一个数据框中
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
  geom_point() +  # 绘制散点，大小根据count变量
  geom_text(vjust = -0.5, hjust = 0.5) +  # 在每个点旁边显示chr_type
  scale_x_continuous(name = "Male Proportion (%)", limits = c(0, 2)) +
  scale_y_continuous(name = "Female Proportion (%)", limits = c(0, 2)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey") + 
  scale_color_manual(values = c("Gain" = "#56B4E9", "Loss" = "#E69F00", "CN-LOH" = "#009E73")) +  # 定义颜色
  scale_size_continuous(range = c(1, 5), name = "Count") +  # 调整点的大小范围
  theme_minimal() +
  labs(title = "Chromosome Type Proportion in Male vs Female",
       subtitle = "Colored by Chromosome Event Type")

#----5.Overlapping-------
overlap_cnloh<-filter(data_filtered,data_filtered$Count>1 & data_filtered$type=="CN-LOH")
overlap_cnloh<- overlap_cnloh %>%
  group_by(sample_id, chrom) %>%
  filter(n() > 1) %>%
  ungroup()
#function of selecting overlapping
library(purrr)
contains_other_interval <- function(beg1, end1, beg2, end2) {
  return((beg1 < beg2 & end1 > beg2) | (beg1 < end2 & end1 > end2))
}
#Run function in CN-LOH
overlap_cnloh_filtered <- NA
overlap_cnloh_filtered <- overlap_cnloh %>%
  group_by(sample_id, chrom) %>%
  filter(
    any(
      map2_lgl(seq_along(beg_GRCh38), beg_GRCh38, function(i, beg) {
        any(map2_lgl(beg_GRCh38[-i], end_GRCh38[-i], ~ contains_other_interval(beg, end_GRCh38[i], .x, .y)))
      })
    )
  ) %>%
  ungroup()
#Gain/ Loss
overlap_gl<-filter(data_filtered,data_filtered$Count>1 & data_filtered$type!="CN-LOH")
overlap_gl<-filter(overlap_gl,overlap_gl$Count>1 & overlap_gl$type!="Undetermined")
overlap_gl<- overlap_gl %>%
  group_by(sample_id, chrom) %>%
  filter(n() > 1) %>%
  ungroup()
#Run function in G/L
overlap_gl_filtered <- NA
overlap_gl_filtered <- overlap_gl %>%
  group_by(sample_id, chrom) %>%
  filter(
    any(
      map2_lgl(seq_along(beg_GRCh38), beg_GRCh38, function(i, beg) {
        any(map2_lgl(beg_GRCh38[-i], end_GRCh38[-i], ~ contains_other_interval(beg, end_GRCh38[i], .x, .y)))
      })
    )
  ) %>%
  ungroup()
#Calculate overlapping rate
overlap_length <- function(beg1, end1, beg2, end2) {
  overlap_beg <- max(beg1, beg2)
  overlap_end <- min(end1, end2)
  if (overlap_beg <= overlap_end) {
    return(overlap_end - overlap_beg)
  } else {
    return(0)
  }
}
overlap_gl_filtered2<-NA
overlap_gl_filtered2 <- overlap_gl_filtered %>%
  group_by(sample_id, chrom) %>%
  mutate(
    overlap_rate = map_dbl(seq_along(beg_GRCh38), function(i) {
      overlaps <- map2_dbl(beg_GRCh38[-i], end_GRCh38[-i], ~ overlap_length(beg_GRCh38[i], end_GRCh38[i], .x, .y))
      max_overlap <- max(overlaps, na.rm = TRUE)
      if (max_overlap > 0) {
        return(max_overlap / (end_GRCh38[i] - beg_GRCh38[i]))
      } else {
        return(0)
      }
    })
  ) %>%
  ungroup()

overlap_gl_filtered3<-filter(overlap_gl_filtered2,overlap_gl_filtered2$overlap_rate>0)
overlap_gl_filtered3$overlap_rate<-ifelse(overlap_gl_filtered3$overlap_rate==1,NA,overlap_gl_filtered3$overlap_rate)
overlap_gl_filtered20<-filter(overlap_gl_filtered3,overlap_gl_filtered3$overlap_rate>0.2)
describeBy(overlap_gl_filtered20$overlap_rate)
table(overlap_gl_filtered20$type,overlap_gl_filtered20$chrom)

#----7.Gain15 and Diseases
table(data_mca$chr15gain)
data_mca <- filter(data_all,data_all$mca_status==1)
data_mca <- data_mca %>% filter(chr15gain %in% c(1, 2))
data_mca <- data_mca %>% mutate(chr15gain = factor(chr15gain, levels = c(2, 1)))
#Cox analysis
analyze_cox_model(data_all, "survial_cll", "cll", "chr15gain", c("age", "age_squ", "sex_at_birth", "smoking", "race"))
analyze_cox_model(data_mca, "survial_cad", "cad", "chr15gain", c("age", "age_squ", "sex_at_birth", "smoking", "race"))
analyze_cox_model(data_mca, "survial_hf", "heart_failure", "chr15gain", c("age", "age_squ", "sex_at_birth", "smoking", "race"))
analyze_cox_model(data_mca, "survial_ai", "AI", "chr15gain", c("age", "age_squ", "sex_at_birth", "smoking", "race"))
analyze_cox_model(data_mca, "survial_ckd", "ckd", "chr15gain", c("age", "age_squ", "sex_at_birth", "smoking", "race"))
#Linear analysis
lm_model <- lm(ALC_count ~ chr15gain + age + age_squ + sex_at_birth + smoking + race, data = data_mca)
summary(lm_model)
lm_model <- lm(WBC_count ~ chr15gain + age + age_squ + sex_at_birth + smoking + race, data = data_mca)
summary(lm_model)
#Pair (Case control)
install.packages("MatchIt")
library(MatchIt)
ls(data_mca)
table(data_mca$chr15gain)
matchit_model<-NA
matchit_model <- matchit(chr15gain ~ age + computed_gender, 
                         data = data_mca, 
                         method = "nearest", 
                         distance = "robust_mahalanobis",
                         ratio = 1,
                         caliper = c(age = 5 / sd(data_mca$age)),
                         exact = ~ computed_gender)
matched_data <- match.data(matchit_model)
class(matched_data$chr15gain)
table(matched_data$computed_gender,matched_data$chr15gain)
describeBy(matched_data$age,matched_data$chr15gain)
print(matched_data)
table(matched_data$chr15gain)

matched_data$strata <- matched_data$subclass
clogit_model <- clogit(cad ~ chr15gain + age + race + smoking+ strata(strata), data = matched_data)
summary(clogit_model)
table(matched_data$race)

#----6.mCA PheWAS-------
#------6.1. Create dataset------
phecode <- read.csv("phecodeX_ICD_WHO_map_flat.csv")
table(phecode$vocabulary_id) #All ICD10
#Keep one decimal place of ICD_10 code of data_condition10 
process_icd10 <- function(x) {
  # 将输入转换为字符格式
  x <- as.character(x)
  
  # 使用正则表达式提取数字部分（包括小数点）
  number <- sub("^[A-Za-z]*", "", x)
  
  # 保留小数点后一位，使用正则截断
  formatted_number <- sub("(\\d+\\.\\d).*", "\\1", number)
  
  # 提取字母部分
  prefix <- sub("[0-9.]+$", "", x)
  
  # 拼接字母部分和处理后的数字部分
  return(paste0(prefix, formatted_number))
}
data_condition10$icd10_processed <- sapply(data_condition10$icd10, process_icd10)
save(data_condition10,file = "data_condition.Rdata")
#megre
#1:Many or many:1 
data_phenocode<-merge(data_condition10,phecode,by.x = "icd10_processed",by.y="icd",all.x = T, allow.cartesian = TRUE)
head(data_phenocode)
test<-filter(data_phenocode,is.na(data_phenocode$icd10)==F & is.na(data_phenocode$phecode)==T)
miss_con<-data.frame(table(test$icd10_processed))
#Deal with the ICD_process according to the UNmerged missing code
process_icd10_based_on_miss_con <- function(icd10_value, miss_con_var1) {
  if (icd10_value %in% miss_con_var1) {
    # 如果有小数点，删除小数点及其后的内容
    if (grepl("\\.", icd10_value)) {
      return(sub("\\..*", "", icd10_value))
    } else {
      # 如果没有小数点，则添加 ".0"
      return(paste0(icd10_value, ".0"))
    }
  } else {
    # 如果不在 miss_con 中，则保持不变
    return(icd10_value)
  }
}
data_condition10$icd10_processed <- sapply(data_condition10$icd10_processed, process_icd10_based_on_miss_con, miss_con_var1 = miss_con$Var1)
#merge again
data_phenocode2<-merge(data_condition10,phecode,by.x = "icd10_processed",by.y="icd",all.x = T,, allow.cartesian = TRUE)
test2<-filter(data_phenocode2,is.na(data_phenocode2$icd10)==F & is.na(data_phenocode2$phecode_string)==T)
miss_con2<-data.frame(table(test2$icd10_processed))
head(data_phenocode2)
phenotype <- data.frame(table(data_phenocode2$phecode_string))
phenotype_cat <- data.frame(table(data_phenocode2$category))
data_phenocode2$phecode_string <- ifelse(is.na(data_phenocode2$phecode_string),data_phenocode2$icd10,data_phenocode2$phecode_string)
save(data_phenocode2,file = "data_phenocode2_phenoX.Rdata")
#Generate new dataset by each phenotype
data_phenocode3<-filter(data_phenocode2,is.na(data_phenocode2$phecode_string)==F) #8018046
str(data_phenocode3)
data_phenocode3 <- data_phenocode3 %>%
  mutate(value = 1)
# data_phenocode3 <- data_phenocode3 %>%
#   mutate(Phenotype = gsub("[^[:alnum:]_]", "_", Phenotype)) 
data_phenocode3 <- data_phenocode3 %>%
  mutate(value = as.numeric(value))

data_pheno <- data_phenocode3[,c("person_id","phecode_string","value")]
names(data_pheno)[2]<-"Phenotype"
data_pheno <- data_pheno %>% distinct() 
data_pheno <- data_pheno %>%
  filter(!is.na(Phenotype) & Phenotype != "")

data_wide <- pivot_wider(
  data_pheno, 
  names_from = Phenotype, 
  values_from = value, 
  values_fill = list(value = 0) 
)

data_wide <- data_wide %>%
  group_by(person_id) %>% 
  summarise(across(everything(), max))
save(data_wide,file = "data_Pheotype_wide.Rdata")

data_all_pheno<-merge(data_all,data_wide,by="person_id",all.x = T)
data_all_pheno <- data_all_pheno %>%
  mutate(across(82:259, ~replace_na(.x, 0)))
save(data_all_pheno,file = "data_all_Pheotype.Rdata")

#Generate new dataset by each phenotype and start time (for COX regression)
data_pheno_cox <-  data_phenocode3[,c("person_id","phecode_string","condition_start_datetime")]
names(data_pheno_cox)[2]<-"Phenotype"
data_pheno_cox <- distinct(data_pheno_cox) #4749524
# 按 person_id 和 Phenotype 分组，保留 condition_start_datetime 最早的一次
data_pheno_cox_filtered <- data_pheno_cox %>%
  group_by(person_id, Phenotype) %>%  # 按 person_id 和 Phenotype 分组
  slice_min(condition_start_datetime) %>%  # 保留每组中最早的 condition_start_datetime
  ungroup()
# 1. 使用 rowid() 创建唯一键 'person_id' + 'Phenotype'，以确保每个 id 和疾病唯一
data_pheno_cox_filtered <- data_pheno_cox_filtered %>%
  mutate(unique_id = rowid(person_id, Phenotype))
data_pheno_cox_filtered <- data_pheno_cox_filtered %>%
  filter(!is.na(Phenotype) & Phenotype != "")
# 2. 创建宽格式数据，Phenotype 作为列名，标记 1
data_wide <- data_pheno_cox_filtered %>%
  mutate(value = 1) %>%  # 创建标记列，用于表示患有某种疾病
  pivot_wider(names_from = Phenotype, values_from = value, values_fill = list(value = 0), id_cols = person_id)  # 使用 person_id 作为 id
# 3. 创建对应的时间列，列名格式为 time_疾病
data_time <- data_pheno_cox_filtered %>%
  pivot_wider(names_from = Phenotype, values_from = condition_start_datetime, 
              names_prefix = "time_", values_fill = list(condition_start_datetime = NA), id_cols = person_id)  # 使用 person_id 作为 id
# 4. 将标记数据和时间数据合并
data_final <- data_wide %>%
  left_join(data_time, by = "person_id") #90489，id is distinct
data_all_pheno_cox<-merge(data_all,data_final,by="person_id",all.x = T)
data_all_pheno_cox <- data_all_pheno_cox %>%
  mutate_at(vars(82:462), ~ replace(., is.na(.),0))

save(data_final,file = "data_wide_for_cox_with_date.Rdata")
save(data_all_pheno_cox,file = "data_all_pheno_cox_PhecodeX.Rdata")
system("gsutil cp data_all_pheno_cox_PhecodeX.Rdata gs://fc-secure-cb192ac6-30ba-46b9-92ee-896a6e36c63e/kun/")
#------6.2. PheWAS mCA (logistics)------
library(broom)
#results <- data.frame()
for (i in 82:259) {
  outcome_var <- names(data_all_pheno)[i]  # 获取结局变量的列名
  
  # 计算总样本量和 mca_status 为 1 时的结局个数
  total_samples <- sum(data_all_pheno[[outcome_var]] == 1, na.rm = TRUE)
  mca_status_1_outcome <- sum(data_all_pheno[[outcome_var]] == 1 & data_all_pheno$mca_status == 1, na.rm = TRUE)  # mca_status 为 1 时的结局个数
  
  # 定义公式，使用暴露变量和调整变量
  formula <- as.formula(paste(outcome_var, "~ mca_status + age + age_squ + race + sex_at_birth + smoking"))
  
  # 运行Logistic回归模型
  model <- glm(formula, data = data_all_pheno, family = binomial)
  
  # 提取模型结果，并转换为数据框
  model_summary <- tidy(model) %>%
    filter(term == "mca_status") %>%
    mutate(outcome = outcome_var,  # 添加结局变量的列名
           total_samples = total_samples,  # 添加总样本量
           mca_status_1_outcome = mca_status_1_outcome)  # 添加 mca_status 为 1 时的结局个数
  
  # 将结果合并到总结果数据框中
  results <- bind_rows(results, model_summary)
}
write.csv(results, "PheWAS_results.csv", row.names = FALSE)
#------6.3. Manhattan plot (logistics)------
library(ggrepel)
system("gsutil cp gs://fc-secure-cb192ac6-30ba-46b9-92ee-896a6e36c63e/hpoisner/backup_files/phecode_definitions1.2.csv .")
phecode_definitions <- read.csv(file = "phecode_definitions1.2.csv")
phecode_definitions2 <- phecode_definitions %>%
  mutate(phenotype = gsub("[^[:alnum:]_]", "_", phenotype)) 
phecode_definitions2 <- phecode_definitions2[,c("phenotype","category")]
results2 <- merge(results,phecode_definitions2,by.x = "outcome", by.y = "phenotype",all.x = T)
results2$category <- ifelse(is.na(results2$category),"neoplasms",results2$category) #One condition didnt merged which is neoplasms
# 1. 计算 -log10(p.value) 和创建 outcome_index 列
results2 <- results2 %>%
  mutate(log_p_value = -log10(p.value),  # 计算 -log10(p.value)
         outcome_index = row_number(),   # 创建 outcome_index 列
         shape = ifelse(estimate > 0, 24, 25))  # 根据 estimate 设置三角形方向

# 2. 找到 -log10(p.value) 最高的前10个结局变量
top_10_outcomes <- results2 %>%
  arrange(desc(log_p_value)) %>%
  slice_head(n = 10) %>%
  pull(outcome)  # 提取前10个结局变量的名称

# 3. 标记前10个结局变量
results2 <- results2 %>%
  mutate(top_10_label = ifelse(outcome %in% top_10_outcomes, outcome, NA))  # 只有前10个结局变量会显示标签

# 4. 设置显著性阈值
threshold_pvalue_0_05 <- -log10(0.05)  # p-value < 0.05 的阈值

results2 <- results2 %>%
  arrange(category, outcome_index) %>%
  mutate(outcome_index = row_number())  # 重新创建 outcome_index，以确保它在每个 category 中是连续的

# 5. 绘制曼哈顿图，添加显著性阈值线并只显示前10个结局变量的图例
ggplot(results2, aes(x = outcome_index, y = log_p_value)) +
  geom_point(aes(color = category, shape = as.factor(shape), fill = category), size = 3) +  # 填充和边框颜色一致
  scale_shape_manual(values = c(24, 25)) +  # 设置实心三角形的 shape 值
  geom_hline(yintercept = threshold_pvalue_0_05, linetype = "dashed", color = "red") +  # 添加 p-value < 0.05 的阈值线
  geom_text_repel(aes(label = top_10_label), size = 3, max.overlaps = 15, min.segment.length = 0) +  # 使用 ggrepel 自动调整标签位置
  theme_minimal() +
  labs(x = "Category", y = "-log10(p.value)", title = "Manhattan Plot with Categories Colored") +
  scale_color_brewer(palette = "Set3", name = "Category") +  # 使用颜色调色板
  scale_fill_brewer(palette = "Set3", name = "Category") +  # 设置填充颜色与边框颜色一致
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 45, hjust = 1)) +  # 调整x轴标签的角度
  scale_x_continuous(breaks = results2 %>% group_by(category) %>% summarize(mid = mean(outcome_index)) %>% pull(mid),
                     labels = results2 %>% group_by(category) %>% summarize(cat = first(category)) %>% pull(cat))  # 设置 x 轴标签为 category 名称
#------6.4. PheWAS mCA (COX regression)------
cox_results <- list()
summary(data_all_pheno_cox$age)
disease_cols <- colnames(data_all_pheno_cox)[82:100] #82:462
#time_cols <- colnames(data_all_pheno_cox)[260:437]
for (disease_col in disease_cols) {
  time_col <- paste0("time_", disease_col)
  if (time_col %in% colnames(data_all_pheno_cox)) {
    data_surv <- data_all_pheno_cox %>%
      mutate(
        surv_time = case_when(
          !!sym(disease_col) == 1 ~ as.numeric(difftime(!!sym(time_col), min_date, units = "days")) / 365.25,  
          !!sym(disease_col) == 0 ~ as.numeric(difftime(max_date, min_date, units = "days")) / 365.25  
        )
      ) %>%
      filter(surv_time > 0)  
    
    total_samples <- sum(data_surv[[disease_col]] == 1)
    mca_status_1_outcome <- sum(data_surv$mca_status == 1 & data_surv[[disease_col]] == 1)
    

    surv_obj <- Surv(time = data_surv$surv_time, event = data_surv[[disease_col]])

    cox_model <- coxph(surv_obj ~ mca_status + age + age_squ+race+smoking+sex_at_birth, data = data_surv)

    cox_result <- tidy(cox_model) %>%
      filter(term == "mca_status") %>% 
      mutate(total_samples = total_samples,  
             mca_status_1_outcome = mca_status_1_outcome) 
    
    cox_results[[disease_col]] <- cox_result
  } else {
    warning(paste("Time:", time_col, "Not exist"))
  }
}
results_df <- bind_rows(cox_results, .id = "disease")
write.csv(results_df, "PheWAS_results_cox.csv", row.names = FALSE)
#----7.Myeloid CH and Lymphoid CH-----
#------7.1. Define mCH and lCH-----
mca_df$mCH<-NA
mca_df$mCH<-ifelse(mca_df$type=="Loss" & mca_df$chrom =="chr12" & mca_df$q_arm !="N",1,mca_df$mCH)
mca_df$mCH<-ifelse(mca_df$type=="Loss" & mca_df$chrom =="chr20" & mca_df$q_arm !="N",1,mca_df$mCH)
mca_df$mCH<-ifelse(mca_df$type=="Loss" & mca_df$chrom =="chr5" & mca_df$q_arm !="N",1,mca_df$mCH)
mca_df$mCH<-ifelse(mca_df$type=="Gain" & mca_df$chrom =="chr1" & mca_df$q_arm !="N",1,mca_df$mCH)
mca_df$mCH<-ifelse(mca_df$type=="Gain" & mca_df$chrom =="chr9" & mca_df$p_arm !="N",1,mca_df$mCH)
mca_df$mCH<-ifelse(mca_df$type=="CN-LOH" & mca_df$chrom =="chr22" & mca_df$q_arm !="N",1,mca_df$mCH)
mca_df$mCH<-ifelse(mca_df$type=="CN-LOH" & mca_df$chrom =="chr9" & mca_df$p_arm !="N",1,mca_df$mCH)
mca_df$mCH<-ifelse(mca_df$type=="CN-LOH" & mca_df$chrom =="chr14" & mca_df$q_arm !="N",1,mca_df$mCH)
mca_df$mCH<-ifelse(mca_df$type=="Gain" & mca_df$chrom =="chr8" & mca_df$q_arm =="Y" & mca_df$p_arm =="Y",1,mca_df$mCH)
mca_df$mCH<-ifelse(mca_df$type=="Gain" & mca_df$chrom =="chr8" & mca_df$q_arm =="T" & mca_df$p_arm =="T",1,mca_df$mCH)
mca_df$mCH<-ifelse(mca_df$type=="Gain" & mca_df$chrom =="chr8" & mca_df$q_arm =="Y" & mca_df$p_arm =="T",1,mca_df$mCH)
mca_df$mCH<-ifelse(mca_df$type=="Gain" & mca_df$chrom =="chr8" & mca_df$q_arm =="T" & mca_df$p_arm =="Y",1,mca_df$mCH)
table(mca_df$mCH)

mca_df$lCH<-NA
mca_df$lCH <- ifelse(mca_df$type == "Loss" & mca_df$chrom == "chr10" & mca_df$p_arm != "N", 1, mca_df$lCH)
mca_df$lCH <- ifelse(mca_df$type == "Loss" & mca_df$chrom == "chr10" & mca_df$q_arm != "N", 1, mca_df$lCH)
mca_df$lCH <- ifelse(mca_df$type == "Loss" & mca_df$chrom == "chr11" & mca_df$q_arm != "N", 1, mca_df$lCH)
mca_df$lCH <- ifelse(mca_df$type == "Loss" & mca_df$chrom == "chr13" & mca_df$q_arm != "N", 1, mca_df$lCH)
mca_df$lCH <- ifelse(mca_df$type == "Loss" & mca_df$chrom == "chr14" & mca_df$q_arm != "N", 1, mca_df$lCH)
mca_df$lCH <- ifelse(mca_df$type == "Loss" & mca_df$chrom == "chr15" & mca_df$q_arm != "N", 1, mca_df$lCH)
mca_df$lCH <- ifelse(mca_df$type == "Loss" & mca_df$chrom == "chr17" & mca_df$p_arm != "N", 1, mca_df$lCH)
mca_df$lCH <- ifelse(mca_df$type == "Loss" & mca_df$chrom == "chr1" & mca_df$p_arm != "N", 1, mca_df$lCH)
mca_df$lCH <- ifelse(mca_df$type == "Loss" & mca_df$chrom == "chr1" & mca_df$q_arm != "N", 1, mca_df$lCH)
mca_df$lCH <- ifelse(mca_df$type == "Loss" & mca_df$chrom == "chr22" & mca_df$q_arm != "N", 1, mca_df$lCH)
mca_df$lCH <- ifelse(mca_df$type == "Loss" & mca_df$chrom == "chr6" & mca_df$q_arm != "N", 1, mca_df$lCH)
mca_df$lCH <- ifelse(mca_df$type == "Loss" & mca_df$chrom == "chr7" & mca_df$q_arm != "N", 1, mca_df$lCH)
mca_df$lCH <- ifelse(mca_df$type == "Loss" & mca_df$chrom == "chr8" & mca_df$p_arm != "N", 1, mca_df$lCH)
mca_df$lCH <- ifelse(mca_df$type == "Gain" & mca_df$chrom == "chr12" & mca_df$q_arm != "N", 1, mca_df$lCH)
mca_df$lCH <- ifelse(mca_df$type == "Gain" & mca_df$chrom == "chr15" & mca_df$q_arm != "N", 1, mca_df$lCH)
mca_df$lCH <- ifelse(mca_df$type == "Gain" & mca_df$chrom == "chr17" & mca_df$q_arm != "N", 1, mca_df$lCH)
mca_df$lCH <- ifelse(mca_df$type == "Gain" & mca_df$chrom == "chr22" & mca_df$q_arm != "N", 1, mca_df$lCH)
mca_df$lCH <- ifelse(mca_df$type == "Gain" & mca_df$chrom == "chr2" & mca_df$p_arm != "N", 1, mca_df$lCH)
mca_df$lCH <- ifelse(mca_df$type == "Gain" & mca_df$chrom == "chr3" & mca_df$q_arm != "N", 1, mca_df$lCH)
mca_df$lCH <- ifelse(mca_df$type == "Gain" & mca_df$chrom == "chr8" & mca_df$q_arm != "N", 1, mca_df$lCH)
mca_df$lCH <- ifelse(mca_df$type == "Gain" & mca_df$chrom == "chr9" & mca_df$q_arm != "N", 1, mca_df$lCH)
mca_df$lCH <- ifelse(mca_df$type == "CN-LOH" & mca_df$chrom == "chr16" & mca_df$p_arm != "N", 1, mca_df$lCH)
mca_df$lCH <- ifelse(mca_df$type == "CN-LOH" & mca_df$chrom == "chr1" & mca_df$q_arm != "N", 1, mca_df$lCH)
mca_df$lCH <- ifelse(mca_df$type == "CN-LOH" & mca_df$chrom == "chr7" & mca_df$q_arm != "N", 1, mca_df$lCH)
mca_df$lCH <- ifelse(mca_df$type == "CN-LOH" & mca_df$chrom == "chr13" & mca_df$q_arm != "N", 1, mca_df$lCH)
mca_df$lCH <- ifelse(mca_df$type == "CN-LOH" & mca_df$chrom == "chr12" & mca_df$q_arm != "N", 1, mca_df$lCH)
mca_df$lCH <- ifelse(mca_df$type == "CN-LOH" & mca_df$chrom == "chr9" & mca_df$q_arm != "N", 1, mca_df$lCH)
mca_df$lCH <- ifelse(mca_df$type == "Gain" & mca_df$chrom == "chr12" & mca_df$q_arm =="Y" & mca_df$p_arm =="Y",1,mca_df$lCH)
mca_df$lCH <- ifelse(mca_df$type == "Gain" & mca_df$chrom == "chr12" & mca_df$q_arm =="T" & mca_df$p_arm =="Y",1,mca_df$lCH)
mca_df$lCH <- ifelse(mca_df$type == "Gain" & mca_df$chrom == "chr12" & mca_df$q_arm =="Y" & mca_df$p_arm =="T",1,mca_df$lCH)
mca_df$lCH <- ifelse(mca_df$type == "Gain" & mca_df$chrom == "chr12" & mca_df$q_arm =="T" & mca_df$p_arm =="T",1,mca_df$lCH)
mca_df$lCH <- ifelse(mca_df$type == "Gain" & mca_df$chrom == "chr18" & mca_df$q_arm =="Y" & mca_df$p_arm =="Y",1,mca_df$lCH)
mca_df$lCH <- ifelse(mca_df$type == "Gain" & mca_df$chrom == "chr18" & mca_df$q_arm =="T" & mca_df$p_arm =="Y",1,mca_df$lCH)
mca_df$lCH <- ifelse(mca_df$type == "Gain" & mca_df$chrom == "chr18" & mca_df$q_arm =="Y" & mca_df$p_arm =="T",1,mca_df$lCH)
mca_df$lCH <- ifelse(mca_df$type == "Gain" & mca_df$chrom == "chr18" & mca_df$q_arm =="T" & mca_df$p_arm =="T",1,mca_df$lCH)
mca_df$lCH <- ifelse(mca_df$type == "Gain" & mca_df$chrom == "chr19" & mca_df$q_arm =="Y" & mca_df$p_arm =="Y",1,mca_df$lCH)
mca_df$lCH <- ifelse(mca_df$type == "Gain" & mca_df$chrom == "chr19" & mca_df$q_arm =="T" & mca_df$p_arm =="Y",1,mca_df$lCH)
mca_df$lCH <- ifelse(mca_df$type == "Gain" & mca_df$chrom == "chr19" & mca_df$q_arm =="Y" & mca_df$p_arm =="T",1,mca_df$lCH)
mca_df$lCH <- ifelse(mca_df$type == "Gain" & mca_df$chrom == "chr19" & mca_df$q_arm =="T" & mca_df$p_arm =="T",1,mca_df$lCH)
table(mca_df$lCH)

mca_df$aCH <- NA
mca_df$aCH <- ifelse(mca_df$type == "Gain" & mca_df$chrom == "chr21" & mca_df$q_arm != "N", 1, mca_df$aCH)
mca_df$aCH <- ifelse(mca_df$type == "CN-LOH" & mca_df$chrom == "chr11" & mca_df$q_arm != "N", 1, mca_df$aCH)
mca_df$aCH <- ifelse(mca_df$type == "CN-LOH" & mca_df$chrom == "chr16" & mca_df$q_arm != "N", 1, mca_df$aCH)
mca_df$aCH <- ifelse(mca_df$type == "CN-LOH" & mca_df$chrom == "chr1" & mca_df$p_arm != "N", 1, mca_df$aCH)
mca_df$aCH <- ifelse(mca_df$type == "CN-LOH" & mca_df$chrom == "chr17" & mca_df$p_arm != "N", 1, mca_df$aCH)

mca_df$lCH <- ifelse(is.na(mca_df$aCH) ==F & mca_df$aCH==1,NA,mca_df$lCH)
mca_df$mCH <- ifelse(is.na(mca_df$aCH) ==F & mca_df$aCH==1,NA,mca_df$mCH)
mca_df$mCH <- ifelse(is.na(mca_df$lCH) ==F & is.na(mca_df$mCH) ==F  & mca_df$lCH==1 & mca_df$mCH==1,NA,mca_df$mCH)
mca_df$lCH <- ifelse(is.na(mca_df$lCH) ==F & is.na(mca_df$mCH) ==F  & mca_df$lCH==1 & mca_df$mCH==1,NA,mca_df$lCH)

table(mca_df$mCH)
table(mca_df$mCH,mca_df$lCH)
test<-filter(mca_df,mca_df$mCH==1 & mca_df$lCH==1)
#------7.2. PheWAS of mCH------
mch_person_ids <- mca_df %>%
  filter(mCH == 1) %>%
  pull(sample_id)
lch_person_ids <- mca_df %>%
  filter(lCH == 1) %>%
  pull(sample_id)
chip_person_ids <- chip_carriers$person_id

data_all_pheno_cox$lCH<-NA
data_all_pheno_cox <- data_all_pheno_cox %>%
  mutate(lCH = ifelse(person_id %in% lch_person_ids, 1, 0))
data_all_pheno_cox$lCH <- ifelse(data_all_pheno_cox$lCH == 0 & data_all_pheno_cox$mca_status==1,NA,data_all_pheno_cox$lCH)
table(data_all_pheno_cox$lCH) #2260


data_all_wgs_cox <- data_all_wgs_cox_PhecodeX
data_all_wgs_cox$mCH<-NA
data_all_wgs_cox <- data_all_wgs_cox %>%
  mutate(mCH = ifelse(person_id %in% chip_person_ids, 1, data_all_wgs_cox$mCH))
data_all_wgs_cox <- data_all_wgs_cox %>%
  mutate(mCH = ifelse(person_id %in% mch_person_ids, 1, data_all_wgs_cox$mCH))
data_all_wgs_cox$mCH <- ifelse(is.na(data_all_wgs_cox$mCH),0,data_all_wgs_cox$mCH)
data_all_wgs_cox$mCH <- ifelse(data_all_wgs_cox$mCH == 0 & data_all_wgs_cox$mca_status==1,NA,data_all_wgs_cox$mCH)
table(data_all_wgs_cox$mCH) #10213
save(data_all_wgs_cox,file = "data_all_wgs_cox_PhecodeX.Rdata")
save(data_all_pheno_cox,file = "data_all_pheno_cox_PhecodeX.Rdata")
system("gsutil cp data_all_wgs_cox_PhecodeX.Rdata gs://fc-secure-cb192ac6-30ba-46b9-92ee-896a6e36c63e/kun/")
system("gsutil cp data_all_pheno_cox_PhecodeX.Rdata gs://fc-secure-cb192ac6-30ba-46b9-92ee-896a6e36c63e/kun/")


#mCH PheWAS in data_all_wgs_cox
cox_results <- list()
disease_cols <- colnames(data_all_wgs_cox)[82:259]
for (disease_col in disease_cols) {
  # 动态找到对应的时间列（格式为 time_疾病）
  time_col <- paste0("time_", disease_col)
  
  # 确保时间列存在
  if (time_col %in% colnames(data_all_wgs_cox)) {
    
    # 计算生存时间：对于结局为1的记录，生存时间为疾病发生的时间减去min_date；对于结局为0的记录，生存时间为max_date减去min_date
    data_surv <- data_all_wgs_cox %>%
      mutate(
        surv_time = case_when(
          !!sym(disease_col) == 1 ~ as.numeric(difftime(!!sym(time_col), min_date, units = "days")) / 365.25,  # 生存时间转换为年
          !!sym(disease_col) == 0 ~ as.numeric(difftime(max_date, min_date, units = "days")) / 365.25  # 未患病的生存时间也转换为年
        )
      ) %>%
      filter(surv_time > 0)  # 舍去生存时间 <= 0 的行
    
    # 计算样本量和 mCH == 1 且结局发生的个数
    total_samples <- sum(data_surv[[disease_col]] == 1)
    mca_outcome <- sum(is.na(data_surv$mCH)==F & data_surv$mCH == 1 & data_surv[[disease_col]] == 1)
    
    # 创建 Surv 对象
    surv_obj <- Surv(time = data_surv$surv_time, event = data_surv[[disease_col]])
    
    # 执行 Cox 回归，调整变量为 mCH, age, age_squ, race, smoking, sex_at_birth
    cox_model <- coxph(surv_obj ~ mCH + age + age_squ + race + smoking + sex_at_birth, data = data_surv)
    
    # 提取回归结果，并添加样本量和 mCH == 1 时的结局个数
    cox_result <- tidy(cox_model) %>%
      filter(term == "mCH") %>%  # 只保留 mCH 的结果
      mutate(total_samples = total_samples,  # 添加总样本量
             mca_outcome = mca_outcome)  # 添加 mCH == 1 时的结局个数
    
    # 保存结果到列表中
    cox_results[[disease_col]] <- cox_result
  } else {
    warning(paste("时间列", time_col, "不存在，跳过该疾病"))
  }
}
# 将所有结果合并成一个数据框
results_df <- bind_rows(cox_results, .id = "disease")
write.csv(results_df,file = "PheWAS_m_CH_cox.csv")

#MANHATTAN plot-----
phecodeX <- read.csv("phecodeX_ICD_WHO_map_flat.csv")
phecodeX2 <- phecodeX[,c(4,5,7)]
phecodeX2 <- distinct(phecodeX2)
results_df <- read.csv("mCH_PheWAS_results_cox_phecodeX.csv")
results2 <- merge(results_df,phecodeX2,by.x = "disease", by.y = "phecode_string",all.x = T)
results2 <- filter(results2,is.na(results2$estimate) == F)
names(results2)[1]<-"outcome"
# 1. 计算 -log10(p.value) 和创建 outcome_index 列
results2 <- results2 %>%
  mutate(log_p_value = -log10(p.value),  # 计算 -log10(p.value)
         outcome_index = row_number(),   # 创建 outcome_index 列
         shape = ifelse(estimate > 0, 24, 25))  # 根据 estimate 设置三角形方向

# 2. 找到 -log10(p.value) 最高的前10个结局变量
top_10_outcomes <- results2 %>%
  arrange(desc(log_p_value)) %>%
  slice_head(n = 10) %>%
  pull(outcome)  # 提取前10个结局变量的名称

# 3. 标记前10个结局变量
results2 <- filter(results2,results2$total_samples>=20)

results2 <- results2 %>%
  mutate(top_10_label = ifelse(outcome %in% top_10_outcomes, outcome, NA))  # 只有前10个结局变量会显示标签

# 4. 设置显著性阈值
threshold_pvalue_0_05 <- 4.3  # p-value < 0.05 的阈值

results2 <- results2 %>%
  arrange(category, outcome_index) %>%
  mutate(outcome_index = row_number())  # 重新创建 outcome_index，以确保它在每个 category 中是连续的

# 过滤掉缺失值和超出范围的值
results2_clean <- results2 %>%
  filter(!is.na(log_p_value), !is.na(outcome_index))

table(results2_clean$category)

results2_clean <- filter(results2_clean,is.na(results2_clean$category) ==F)
results2_clean$category<- as.factor(results2_clean$category)
# 5. 绘制曼哈顿图，添加显著性阈值线并只显示前10个结局变量的图例
num_categories <- n_distinct(results2_clean$category)
ggplot(results2_clean, aes(x = outcome_index, y = log_p_value)) +
  geom_point(aes(color = category, shape = as.factor(shape), fill = category), size = 3) +
  scale_shape_manual(values = c(24, 25)) +
  geom_hline(yintercept = threshold_pvalue_0_05, linetype = "dashed", color = "red") +
  geom_text_repel(aes(label = top_10_label), size = 2.5, max.overlaps = 20, min.segment.length = 0) +  # 使用 ggrepel 自动调整标签位置
  theme_minimal() +
  labs(x = "Category", y = "-log10(p.value)", title = "Manhattan Plot of m_CH PheWAS in All of Us (40-90)") +
  scale_color_viridis_d(name = "Category", option = "plasma") +  # 可以选择 viridis, magma, plasma, inferno, cividis
  scale_fill_viridis_d(name = "Category", option = "plasma") +  # 使用相同的 option 保持一致
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_continuous(breaks = results2_clean %>% group_by(category) %>% summarize(mid = mean(outcome_index)) %>% pull(mid),
                     labels = results2_clean %>% group_by(category) %>% summarize(cat = first(category)) %>% pull(cat))

#------7.3. PheWAS of lCH------
#lCH PheWAS in data_all_pheno_cox
cox_results <- list()
disease_cols <- colnames(data_all_wgs_cox)[82:259]
for (disease_col in disease_cols) {
  # 动态找到对应的时间列（格式为 time_疾病）
  time_col <- paste0("time_", disease_col)
  
  # 确保时间列存在
  if (time_col %in% colnames(data_all_pheno_cox)) {
    
    # 计算生存时间：对于结局为1的记录，生存时间为疾病发生的时间减去min_date；对于结局为0的记录，生存时间为max_date减去min_date
    data_surv <- data_all_pheno_cox %>%
      mutate(
        surv_time = case_when(
          !!sym(disease_col) == 1 ~ as.numeric(difftime(!!sym(time_col), min_date, units = "days")) / 365.25,  # 生存时间转换为年
          !!sym(disease_col) == 0 ~ as.numeric(difftime(max_date, min_date, units = "days")) / 365.25  # 未患病的生存时间也转换为年
        )
      ) %>%
      filter(surv_time > 0)  # 舍去生存时间 <= 0 的行
    
    # 计算样本量和 lch == 1 且结局发生的个数
    total_samples <- sum(data_surv[[disease_col]] == 1)
    chip_1_outcome <- sum( is.na(data_surv$lCH)==F & data_surv$lCH == 1 & data_surv[[disease_col]] == 1)
    
    # 创建 Surv 对象
    surv_obj <- Surv(time = data_surv$surv_time, event = data_surv[[disease_col]])
    
    # 执行 Cox 回归，调整变量为 lCH, age, age_squ, race, smoking, sex_at_birth
    cox_model <- coxph(surv_obj ~ lCH + age + age_squ + race + smoking + sex_at_birth, data = data_surv)
    
    # 提取回归结果，并添加样本量和 lCH == 1 时的结局个数
    cox_result <- tidy(cox_model) %>%
      filter(term == "lCH") %>%  # 只保留 lCH 的结果
      mutate(total_samples = total_samples,  # 添加总样本量
             chip_1_outcome = chip_1_outcome)  # 添加 lCH == 1 时的结局个数
    
    # 保存结果到列表中
    cox_results[[disease_col]] <- cox_result
  } else {
    warning(paste("时间列", time_col, "不存在，跳过该疾病"))
  }
}
# 将所有结果合并成一个数据框
results_lCH <- bind_rows(cox_results, .id = "disease")
write.csv(results_lCH,file = "PheWAS_l_CH_cox.csv")
#MANHATTAN plot------
results_lCH <-  read.csv("lCH_PheWAS_results_cox_phecodeX.csv")
results2 <- merge(results_lCH,phecodeX2,by.x = "disease", by.y = "phecode_string",all.x = T)
results2 <- filter(results2,is.na(results2$estimate) == F)
names(results2)[1]<-"outcome"
# 1. 计算 -log10(p.value) 和创建 outcome_index 列
results2 <- results2 %>%
  mutate(log_p_value = -log10(p.value),  # 计算 -log10(p.value)
         outcome_index = row_number(),   # 创建 outcome_index 列
         shape = ifelse(estimate > 0, 24, 25))  # 根据 estimate 设置三角形方向

# 2. 找到 -log10(p.value) 最高的前10个结局变量
top_10_outcomes <- results2 %>%
  arrange(desc(log_p_value)) %>%
  slice_head(n = 10) %>%
  pull(outcome)  # 提取前10个结局变量的名称

# 3. 标记前10个结局变量
results2 <- results2 %>%
  mutate(top_10_label = ifelse(outcome %in% top_10_outcomes, outcome, NA))  # 只有前10个结局变量会显示标签

# 4. 设置显著性阈值
threshold_pvalue_0_05 <- 4.3  # p-value < 0.05 的阈值

results2 <- results2 %>%
  arrange(category, outcome_index) %>%
  mutate(outcome_index = row_number())  # 重新创建 outcome_index，以确保它在每个 category 中是连续的

results2<- filter(results2,results2$total_samples>50)
table(results2$category)
results2_clean <- filter(results2,is.na(results2$category) ==F)
# 5. 绘制曼哈顿图，添加显著性阈值线并只显示前10个结局变量的图例
ggplot(results2_clean, aes(x = outcome_index, y = log_p_value)) +
  geom_point(aes(color = category, shape = as.factor(shape), fill = category), size = 3) +
  scale_shape_manual(values = c(24, 25)) +
  geom_hline(yintercept = threshold_pvalue_0_05, linetype = "dashed", color = "red") +
  geom_text_repel(aes(label = top_10_label), size = 2.5, max.overlaps = 20, min.segment.length = 0) +  # 使用 ggrepel 自动调整标签位置
  theme_minimal() +
  labs(x = "Category", y = "-log10(p.value)", title = "Manhattan Plot of l_CH PheWAS in All of Us (40-90)") +
  scale_color_viridis_d(name = "Category", option = "plasma") +  # 可以选择 viridis, magma, plasma, inferno, cividis
  scale_fill_viridis_d(name = "Category", option = "plasma") +  # 使用相同的 option 保持一致
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_continuous(breaks = results2_clean %>% group_by(category) %>% summarize(mid = mean(outcome_index)) %>% pull(mid),
                     labels = results2_clean %>% group_by(category) %>% summarize(cat = first(category)) %>% pull(cat))
#------9.4. Draw Scatter Plot of mCH and lCH (without A-CH)--------
results_mCH <- read.csv("PheWAS_m_CH_cox.csv")
phewas_mCH <- results_mCH %>%
  rename_with(~ paste0("mCH_", .), -disease) 
results_lCH <- read.csv("PheWAS_l_CH_cox.csv")
phewas_lCH <- results_lCH %>%
  rename_with(~ paste0("lCH_", .), -disease) 
# 基于 outcome 列合并两个数据集
phewas_combined <- phewas_mCH %>%
  inner_join(phewas_lCH, by = "disease")
names(phewas_combined)[2]<-"outcome"
phewas_filtered <- phewas_combined %>%
  filter(mCH_p.value <= 0.05 | lCH_p.value <= 0.05)
phewas_filtered <- phewas_filtered %>%
  filter(mCH_total_samples >=20 & lCH_total_samples >=20 & lCH_chip_1_outcome>0)
#95%CI
phewas_filtered <- phewas_filtered %>%
  mutate(
    mCH_ci_lower = mCH_estimate - 1.96 * mCH_std.error,  # mCA 95% CI 下限
    mCH_ci_upper = mCH_estimate + 1.96 * mCH_std.error,  # mCA 95% CI 上限
    lCH_ci_lower = lCH_estimate - 1.96 * lCH_std.error,  # CHIP 95% CI 下限
    lCH_ci_upper = lCH_estimate + 1.96 * lCH_std.error   # CHIP 95% CI 上限
  )
max_range <- max(abs(phewas_filtered$mCH_estimate), abs(phewas_filtered$lCH_estimate))

#
phewas_filtered <- phewas_filtered %>%
  mutate(
    log_mCH_total_samples = log10(mCH_total_samples),  # 对 mca_total_samples 取对数（+1避免 log(0)）
    log_lCH_total_samples = log10(lCH_total_samples),  # 对 chip_total_samples 取对数（+1避免 log(0)）
    ellipse_a = rescale(log_mCH_total_samples, to = c(0.05, 0.3)),  # 缩放横向直径到合理范围
    ellipse_b = rescale(log_lCH_total_samples, to = c(0.05, 0.3))  # 缩放纵向直径到合理范围
  )

# 2. 创建颜色分类和对应的 p 值用于颜色深浅
phewas_filtered <- phewas_filtered %>%
  mutate(
    color_group = case_when(
      mCH_p.value <= 0.05 & lCH_p.value <= 0.05 ~ "Both Significant",
      mCH_p.value <= 0.05 & lCH_p.value > 0.05 ~ "m_CH Significant",
      mCH_p.value > 0.05 & lCH_p.value <= 0.05 ~ "l_CH Significant"
    ),
    # 计算透明度，p 值越小，透明度越大（取值范围在0.3到1之间）
    alpha_value = case_when(
      color_group == "Both Significant" ~ rescale(1 - (mCH_p.value + lCH_p.value) / 2, to = c(0.3, 1)),
      color_group == "m_CH Significant" ~ rescale(1 - mCH_p.value, to = c(0.3, 1)),
      color_group == "l_CH Significant" ~ rescale(1 - lCH_p.value, to = c(0.3, 1))
    )
  )

# 3. 筛选出 mca_p.value 和 chip_p.value 都 <= 0.05 的疾病
significant_points2 <- phewas_filtered %>% filter(phewas_filtered$color_group == "Both Significant")
#  绘制散点图，确保坐标轴在0处相交，刻度一致，并为每个点分配不同的颜色
#install.packages("scales")
#install.packages("ggforce")
library(ggforce)
library(scales)
library(ggrepel)

ggplot(phewas_filtered, aes(x = mCH_estimate, y = lCH_estimate, fill = color_group)) +
  geom_ellipse(aes(x0 = mCH_estimate, y0 = lCH_estimate, a = ellipse_a, b = ellipse_b, angle = 0, alpha = alpha_value), 
               color = NA) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  geom_text_repel(data = significant_points2, aes(label = outcome), size = 3, 
                  max.overlaps = 15, box.padding = 0.5, point.padding = 0.5, 
                  segment.color = "gray", segment.size = 0.5) +
  scale_fill_manual(values = c("m_CH Significant" = "#E69F00", 
                               "l_CH Significant" = "#56B4E9", 
                               "Both Significant" = "#009E73")) +
  theme_classic() +
  labs(
    title = "Scatter Plot of Significant mCH and lCH PheWAS Estimates",
    x = "mCH Estimate",
    y = "lCH Estimate",
    fill = "Significance",
    alpha = "p-value Transparency"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.position = "right",
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  ) +
  coord_fixed(ratio = 1, xlim = c(-max_range, max_range), ylim = c(-max_range, max_range))




#----8. CLL PRS----
system("gsutil cp gs://fc-secure-cb192ac6-30ba-46b9-92ee-896a6e36c63e/pershy1/cll_pgs/cll_pgs_v7_yp09182024.tsv .", intern = TRUE)
CLL_prs_aou <- fread('cll_pgs_v7_yp09182024.tsv') 
data_all_pheno_cox <- merge(data_all_pheno_cox,CLL_prs_aou,by.x = "person_id", by.y = "IID", all.x=T)
summary(data_all_pheno_cox$cll_prs)
save(data_all_pheno_cox,file = "data_all_pheno_cox_PhecodeX.Rdata")

data_all_prs <- filter(data_all_pheno_cox,is.na(data_all_pheno_cox$cll_prs)==F)
data_all_prs$cll_262 <- 2
data_all_prs$cll_262[data_all_prs$cll_prs >= quantile(data_all_prs$cll_prs,0.8)] <- 3
data_all_prs$cll_262[data_all_prs$cll_prs < quantile(data_all_prs$cll_prs,0.2)] <- 1
data_all_prs$cll_262<-as.factor(data_all_prs$cll_262)
table(data_all_prs$cll_262)
#CLL_PRS and mCA (mediator) and CLL
install.packages("mediation")
library(mediation)
model.m <- glm(mca_status ~ cll_prs + age + age_squ+ sex_at_birth+ smoking+ race, 
               data = data_all_pheno_cox, family = binomial)
summary(model.m)

model.y <- glm(cll ~ mca_status + cll_prs + age + age_squ+ sex_at_birth+ smoking+ race, 
               data = data_all_pheno_cox, family = binomial)
summary(model.y)

model.y <- glm(cll ~ mca_status + cll_prs + age + age_squ+ sex_at_birth+ smoking+ race + mca_status*cll_prs, 
               data = data_all_pheno_cox, family = binomial)
summary(model.y)

med.out <- mediate(model.m, model.y, treat = "cll_prs", mediator = "mca_status", robustSE = TRUE, sims = 100)
summary(med.out)

#CLL_PRS and mCA highrisk (mediator) and CLL
install.packages("mediation")
library(mediation)
table(data_all_pheno_cox$mca_highrisk)

data_all_pheno_cox <- filter(data_all_pheno_cox,data_all_pheno_cox$mca_highrisk!=2)
model.m <- glm(mca_highrisk ~ cll_prs + age + age_squ+ sex_at_birth+ smoking+ race, 
               data = data_all_pheno_cox, family = binomial)
summary(model.m)

model.y <- glm(cll ~ mca_highrisk + cll_prs + age + age_squ+ sex_at_birth+ smoking+ race, 
               data = data_all_pheno_cox, family = binomial)
summary(model.y)

model.inter <- glm(cll ~ mca_highrisk + cll_prs + age + age_squ+ sex_at_birth+ smoking+ race + mca_highrisk*cll_prs, 
                   data = data_all_pheno_cox_highrisk, family = binomial)
summary(model.inter)

med.out <- mediate(model.m, model.y, treat = "cll_prs", mediator = "mca_highrisk", sims = 100)
summary(med.out)


#-BioVU PheWAS (jupyter)
data_biovu2 <- fread('BioVU_mCA_calls.txt') 
system("gsutil cp gs://fc-secure-cb192ac6-30ba-46b9-92ee-896a6e36c63e/yp_mca_calls_021924.txt .", intern = TRUE)
mca_df <- fread('yp_mca_calls_021924.txt') #19799

table(data_biovu$Autosomes)
data_biovu_all <- merge(data_biovu2,data_biovu,by.x = "ourSid", by.y = "sample_id",all.x = T)
table(data_biovu_all$Autosomes,data_biovu_all$cohort)
