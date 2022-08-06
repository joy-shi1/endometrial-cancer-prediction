# DEVELOPMENT AND VALIDATION OF A RISK PREDICTION MODEL FOR ENDOMETRIAL CANCER
#
# PROGRAMMER: Joy Shi
# DATE UPDATED: 1/11/2022
#
# PURPOSE: to develop and validate a risk prediction model for endometrial
#	   cancer
#
# COHORT FOLLOW UP PERIOD: NHS participant age 45 to 85
#
# CASE DEFINITION: all endometrial cancer cases confirmed by medical records, 
#                  death records or tumor registry
#
# INPUT FILES:
#   - NHS I and NHS II questionnaire data (merged and initially cleaned in SAS)
#     	/udd/resjo/E2C2_prediction/E2C2_Validation_2a_ExportedData_NHSI.csv
#     	/udd/resjo/E2C2_prediction/E2C2_Validation_2b_ExportedData_NHSII.csv
#   - NHS I genetic data for specific genetic variants related to endometrial cancer
#     	/udd/resjo/E2C2_prediction/HRC/
#   - Reference data required to develop the risk prediction model:
#     	- Pooled data from case-control studies in Epidemiology of Endometrial Cancer
#         Consortium (E2C2)
# 	    - Baseline rates of endometrial cancer (from SEER)
#       - Baseline rates of competing risks (from CDC WONDER, SEER and BRFSS)
#       - List of predictors and SNPs to be considered in the prediction model
#       - Reference information on SNPs
#	- Cleaned NHANES dataset, used as a reference dataset to estimate the
#   underlying distribution of risk factors. Files can be found in:
#	    /udd/resjo/E2C2_prediction/Input/
# - PLCO questionnaire data:
#     /proj/nhgens/nhgen0z/MODELING/ENDOMETRIAL_CANCER/VALIDATION_PLCO/plco182042916.sas7bdat
#
# EXCLUSIONS AT BASELINE (AGE 45-48):
#  - History of cancer (except non-melona skin cancer)
#  - History of hysterectomy
#  - Restricted to postmenopausal white women
#
# EXCLUSIONS DURING FOLLOW-UP: Follow up ends if participant
#  - Develops the outcome of interest (endometrial cancer)
#  - Dies
#  - Is lost to follow-up
#  - Undergoes a hysterectomy
#  - Develops any cancer (other than non-melanoma skin cancer) 
#
# PREDICTORS: Variables considered for inclusion in the model included:
#  - Age
#  - Smoking status
#  - BMI
#  - Parity
#  - Age at menarche
#  - Hormone therapy use
#  - Oral contraceptive use
#  - Age at first birth
#  - Diabetes
#  - Hypertension
#
# STATISTICAL ANALYSIS: An absolute risk prediction model will be developed
# using the iCare package in R. The model requires three inputs:
#  - A model for the log relative risk parameters for the predictors included
#    in the model; group LASSO will be used on pooled data from E2C2
#    to perform variable selection and regularization
#  - Marginal age-specific incidence rates for endometrial cancer and its
#    competing risks
#  - A reference dataset (NHANES) to estimate the risk factor distribution for
#    the underlying population
# The model will be validated in three cohorts: NHS I, NHS II and PLCO. We
# will develop a clinical-only model (which includes only data from questionnaires)
# and a clinical plus genetic model (which additionally includes previously
# identified SNPs related to endometrial cancer risk). We will assess discrimination
# using AUC and calibration using expected-to-observed ratios. 

# ---------------------------------------------------------------
# -------------------------- (1) Setup --------------------------
# ---------------------------------------------------------------

.libPaths('~/R/x86_64-pc-linux-gnu-library/') 
setwd('/udd/resjo/E2C2_prediction/')

library('withr')
with_makevars(c(PKG_CFLAGS = "-std=c11"), install.packages("rlang"), assignment = "+=")

library('openxlsx')
library('dplyr')
library('purrr')
library('tidyr')
library('stringr')
library('forcats')
library('fastDummies')
library('vcfR')
library('iCARE')
library('haven')
library('gglasso')
library('RColorBrewer')
library('cowplot')
library('broom')

start.age <- 45

# ---------------------------------------------------------------
# ------------ (2) Importing model development inputs -----------
# ---------------------------------------------------------------

# E2C2 Data
e2c2.all <- read.csv('Input/E2C2_Development_1_Clean_E2C2_All.csv')
e2c2.nomiss <- read.csv('Input/E2C2_Development_1_Clean_E2C2_NoMiss.csv')
e2c2.exclusions <- read.csv('Input/E2C2_Development_1_Clean_E2C2_Exclusions.csv')

# Baseline rates of endometrial cancer
base.endo.all <- read.csv('Input/E2C2_Development_2_Endo_All.csv')
base.endo.current <- read.csv('Input/E2C2_Development_2_Endo_Current.csv')
base.endo.nhsi <- read.csv('Input/E2C2_Development_2_Endo_NHSI.csv')
base.endo.nhsii <- read.csv('Input/E2C2_Development_2_Endo_NHSII.csv')
base.endo.plco <- read.csv('Input/E2C2_Development_2_Endo_PLCO.csv')

# Baseline rates of competing risks
base.competing.mortality <- read.csv('Input/E2C2_Development_3_Competing_Mortality.csv')
base.competing.cancer <- read.csv('Input/E2C2_Development_3_Competing_Cancer.csv')
base.competing.current <- read.csv('Input/E2C2_Development_3_Competing_Current.csv')
base.competing.nhsi <- read.csv('Input/E2C2_Development_3_Competing_NHSI.csv')
base.competing.nhsii <- read.csv('Input/E2C2_Development_3_Competing_NHSII.csv')
base.competing.plco <- read.csv('Input/E2C2_Development_3_Competing_PLCO.csv')

# List of predictors
predictors <- read.csv('Input/E2C2_Development_4_Predictors.csv') %>% pull()
snps <- read.csv('Input/E2C2_Development_4_SNPs.csv')

# NHANES dataset for underlying risk factor distribution
nhanes.current <- read.csv('Input/E2C2_Development_5_NHANES_Current.csv')
nhanes.nhsi <- read.csv('Input/E2C2_Development_5_NHANES_NHSI.csv')
nhanes.nhsii <- read.csv('Input/E2C2_Development_5_NHANES_NHSII.csv')
nhanes.plco <- read.csv('Input/E2C2_Development_5_NHANES_PLCO.csv')

# Reference dataset for SNPs
snp.ref <- read.csv('Input/SNP_Reference.csv') %>%
  mutate(SNP38=paste("X", Chr, ".", Position38, sep="")) %>%
  mutate(SNP37=paste("X", Chr, ".", Position37, sep=""))

# ---------------------------------------------------------------
# ---------- (3) NHS I: Importing/cleaning genetic data ---------
# ---------------------------------------------------------------

# --- a. Importing genetic data ---

affy.vcf <- lapply(Sys.glob("./*/*_AffymetrixData.vcf"), read.vcfR)
humancore.vcf <- lapply(Sys.glob("./*/*_HumanCoreExData2.vcf"), read.vcfR)
illumina.vcf <- lapply(Sys.glob("./*/*_IlluminaHumanHapData.vcf"), read.vcfR)
omni.vcf <- lapply(Sys.glob("./*/*_OmniExpressData.vcf"), read.vcfR)
onco.vcf <- lapply(Sys.glob("./*/*_OncoArrayData.vcf"), read.vcfR)


# --- b. Converting VCF to data frames

affy <- do.call(cbind, lapply(1:11, function(i){data.frame(t(extract.gt(affy.vcf[[i]], return.alleles=T, IDtoRowNames=T)))})) 
humancore <- do.call(cbind, lapply(1:12, function(i){data.frame(t(extract.gt(humancore.vcf[[i]], return.alleles=T, IDtoRowNames=T)))})) 
illumina <- do.call(cbind, lapply(1:12, function(i){data.frame(t(extract.gt(illumina.vcf[[i]], return.alleles=T, IDtoRowNames=T)))})) 
omni <- do.call(cbind, lapply(1:12, function(i){data.frame(t(extract.gt(omni.vcf[[i]], return.alleles=T, IDtoRowNames=T)))})) 
onco <- do.call(cbind, lapply(1:12, function(i){data.frame(t(extract.gt(onco.vcf[[i]], return.alleles=T, IDtoRowNames=T)))})) 


# --- c. Combining datasets ---

nhs.genetic <- cbind(id=rownames(affy), affy) %>%
  bind_rows(cbind(id=rownames(humancore), humancore)) %>%
  bind_rows(cbind(id=rownames(illumina), illumina)) %>%
  bind_rows(cbind(id=rownames(omni), omni)) %>%
  bind_rows(cbind(id=rownames(onco), onco)) %>%
  # Remove duplicate IDs
  mutate(dup=duplicated(id)) %>%
  filter(dup==F) %>%
  mutate(id=as.numeric(id)) %>%
  mutate(gwas=1)
  
  
# --- d. Converting genotypes to number of reference alleles ---

for (i in grep("X", colnames(nhs.genetic), value=T)){
  nhs.genetic[,i] <- str_count(nhs.genetic[,i], 
    snp.ref[snp.ref$SNP37==i, 
    "Alternate.Allele"])
  colnames(nhs.genetic)[colnames(nhs.genetic)==i] <- snp.ref[snp.ref$SNP37==i, "rsID"]
  }


# ---------------------------------------------------------------
# ------ (4) NHS I: Importing/cleaning questionnaire data -------
# ---------------------------------------------------------------

# --- a. Importing data ---
nhsi.import <- read.csv("Input/E2C2_Validation_2a_ExportedData_NHSI.csv")

# --- b. Determining exclusions ---
nhsi.all <- nhsi.import %>%
  arrange(id, cancer_dxmonth) %>%
  distinct(id, .keep_all=T) %>%
  left_join(nhs.genetic, by="id") %>%
  
  # Calculating year of cancer diagnosis/death
  mutate(endo.year=endo_dxmonth/12+1900) %>%
  mutate(death.year=(deadmonth/12)+1900) %>%
  mutate(cancer.year=(cancer_dxmonth/12+1900)) %>%
  mutate(cancer.year=ifelse((is.na(cancer.year)|endo.year<cancer.year) & !is.na(endo.year), 
    endo.year, cancer.year)) %>%
  
  # Determining exclusions at baseline
  mutate(in.all=1) %>%
  mutate(ex.cancer=ifelse(!is.na(cancer.year) & cancer.year<=1976, cancer.year, 9999)) %>%
  mutate(ex.hyst=ifelse(!is.na(hyst76) & hyst76==1, 1976, 9999)) %>%
  mutate(ex.date=pmap(list(ex.cancer, ex.hyst), min)) %>%
  mutate(ex.cancer=ifelse(ex.cancer!=9999 & ex.cancer==ex.date, 1, 0)) %>%
  mutate(ex.hyst=ifelse(ex.hyst!=9999 & ex.hyst==ex.date, 1, 0)) %>%
  mutate(ex.hyst=ifelse(ex.cancer==1, 0, ex.hyst)) %>%
  rowwise() %>%
  mutate(menopause.age = median(c(!!! rlang::syms(grep('menoage', names(.), value=TRUE))), na.rm=T)) %>%
  ungroup() %>%
  mutate(ex.white=ifelse(!is.na(race) & race!=1 & ex.hyst==0 & ex.cancer==0, 1, 0)) %>%
  mutate(ex=ifelse(ex.hyst==1|ex.cancer==1|ex.white==1, 1, 0)) %>%
  mutate(include=1-ex)
  
nhsi.inclusion <- nhsi.all %>% 
  summarize_at(vars(in.all, ex, ex.cancer, ex.hyst, ex.white, include), list(~sum(., na.rm=T)))

# --- c. Converting data from wide to long ---
nhsi.long <- nhsi.all %>%
  
  # Removing people excluded at baseline
  filter(include==1) %>%
  
  # Converting data from wide to long
  pivot_longer(cols=c(starts_with("bmi"), starts_with("hyst"), starts_with("menstat"),
                      starts_with("anyhrtever"), starts_with("ephrtever"), 
		      starts_with("ephrtdur"), starts_with("ehrtever"), starts_with("ehrtdur"), 
		      starts_with("smoking"), starts_with("packyears"), starts_with("ocever"), 
		      starts_with("ocdur"), starts_with("parity"), starts_with("hypertension"),
                      starts_with("alcoholgr"), starts_with("completed")),
               names_to=c(".value", "year"), names_sep="(?<=[A-Za-z])(?=[0-9])") %>%
 
   # Creating variable for time (year/age)
  mutate(year=ifelse(as.numeric(year)<20, as.numeric(year)+2000, as.numeric(year)+1900)) %>%
  mutate(age=year-(yobf+1900+mobf/12)) %>%
 
   # Determining timing of events or meeting exclusion criteria over follow-up
  mutate(endo.ever=ifelse(!is.na(endo.year), 1, 0)) %>%
  mutate(hyst=ifelse(!is.na(hyst) & hyst==1, year, 9999)) %>%
  group_by(id) %>%
  mutate(hyst.year=min(hyst)) %>%
  ungroup() %>%
  mutate(include=ifelse(!is.na(age) & !is.na(menopause.age) & age>=menopause.age & 
    !is.na(completed) & age>start.age, year, 9999)) %>%
  mutate(exclude=ifelse(!is.na(completed), year+2, 0000),
         exclude=ifelse(age>85, year, exclude),
         exclude=ifelse(!is.na(death.year) & death.year<exclude, death.year, exclude),
         exclude=ifelse(hyst.year<exclude, hyst.year, exclude),
         exclude=ifelse(!is.na(cancer.year) & cancer.year<exclude, cancer.year, exclude),
         exclude=ifelse(exclude>=2014, 2014, exclude)) %>%
  group_by(id) %>%
 
  # Calculating years for start and end of follow-up
  mutate(start.year=min(include)) %>%
  mutate(end.year=max(exclude)) %>%
  ungroup() %>%
  mutate(start.year=ifelse(start.year>=end.year, 9999, start.year)) %>%
  group_by(id) %>%
  
  # Filling in missing values for parity, ocever, ocdur and smoking (carry last observation forward)
  fill(parity, ocever, ocdur, smoking) %>%
  ungroup() %>%
  
  # Cleaning/recoding covariates
  mutate(education=3) %>%
  mutate(education1=0) %>%
  mutate(education2=0) %>%
  mutate(education3=1) %>%
  mutate(ehrtdur=ehrtdur/12) %>%
  mutate(ephrtdur=ephrtdur/12) %>%
  mutate(ocdur=ocdur/12) %>%
  mutate(diabetes=ifelse(!is.na(db_dxyear) & db_dxyear<=year, 1, 0)) %>%
  mutate(hypertension=ifelse(is.na(hypertension) & !is.na(completed), 0, hypertension),
         hypertension=ifelse(hypertension==2, 0, hypertension)) %>%
  mutate(ehrtever=ifelse(anyhrtever==0, 0, ehrtever)) %>%
  mutate(ephrtever=ifelse(anyhrtever==0, 0, ephrtever)) %>%  
  mutate(anyhrtever=ifelse(is.na(anyhrtever) & lag(anyhrtever)==1 & id==lag(id), 1, anyhrtever)) %>%
  mutate(ehrtever=ifelse(is.na(ehrtever) & lag(ehrtever)==1 & id==lag(id), 1, ehrtever)) %>%
  mutate(ephrtever=ifelse(is.na(ephrtever) & lag(ephrtever)==1 & id==lag(id), 1, ephrtever)) %>%
  arrange(id, -year) %>%
  mutate(anyhrtever=ifelse(is.na(anyhrtever) & lag(anyhrtever)==0 & id==lag(id), 0, anyhrtever)) %>%
  arrange(id, -year) %>%
  mutate(ehrtever=ifelse(is.na(ehrtever) & lag(ehrtever)==0 & id==lag(id), 0, ehrtever)) %>%
  mutate(ephrtever=ifelse(is.na(ephrtever) & lag(ephrtever)==0 & id==lag(id), 0, ephrtever)) %>%  
  arrange(id, year) %>%
  mutate(ocdur=ifelse(ocever==0, 0, ocdur)) %>%
  mutate(ehrtdur=ifelse(ehrtever==0, 0, ehrtdur)) %>%
  mutate(ephrtdur=ifelse(ephrtever==0, 0, ephrtdur)) %>%
  mutate(bmi.group=cut(bmi, breaks=c(-Inf, 18.5, 25, 30, 35, Inf), right=F, labels=c(1, 2, 3, 4, 5))) %>%
  mutate(parity.group=ifelse(parity>=4, 4, parity)) %>%
  mutate(menarchage.group=ifelse(menarchage==0, NA, ifelse(menarchage<=10, 10, 
    ifelse(menarchage>=16, 16, menarchage)))) %>%  
  mutate(menarchage2.group=ifelse(menarchage==0, NA, ifelse(menarchage<=9, 9, 
    ifelse(menarchage %in% c(10, 11), 10, ifelse(menarchage %in% c(12, 13), 12,
    ifelse(menarchage %in% c(14, 15), 14, ifelse(menarchage >=16, 16, NA))))))) %>%
  mutate(ephrtdur.group=ifelse(ephrtdur>0 & ephrtdur<=5, 1, ifelse(ephrtdur>5 & ephrtdur<=10, 
    2, ifelse(ephrtdur!=999 & ephrtdur!=0, 3, ephrtdur)))) %>%
  mutate(ehrtdur.group=ifelse(ehrtdur>0 & ehrtdur<=5, 1, ifelse(ehrtdur>5 & ehrtdur<=10, 
    2, ifelse(ehrtdur!=999 & ehrtdur!=0, 3, ehrtdur)))) %>%  
  mutate(ocdur.group=ifelse(ocdur>0 & ocdur<=5, 1, ifelse(ocdur>5 & ocdur<=10, 2, 
    ifelse(ocdur!=999 & ocdur!=0, 3, ocdur)))) %>%
  mutate(firstbirthage.group=cut(firstbirthage, breaks=c(-Inf, 20, 25, 30, 35, Inf), 
    labels=seq(1:5), right=F)) %>%
  mutate(firstbirthage.group=ifelse(!is.na(parity.group) & parity.group==0, 9, firstbirthage.group)) %>%
  
  # Creating dummy variables for categorical variables
  dummy_cols(., select_columns=c("smoking", "bmi.group", "parity.group", "menarchage.group", 
    "menarchage2.group", "ehrtdur.group", "ephrtdur.group", "ocdur.group", "firstbirthage.group"), ignore_na=T) 
  
names(nhsi.long) <- sub("_", "", names(nhsi.long))


# --- d. Creating interaction terms and finalizing dataset ---

nhsi.final <- nhsi.long %>%
  filter(year==start.year) %>%

  # Interaction between ever OC use and BMI
  mutate(ocever.bmigroup3=ifelse(bmi.group==3, ocever, 0)) %>%
  mutate(ocever.bmigroup4=ifelse(bmi.group==4, ocever, 0)) %>%
  mutate(ocever.bmigroup5=ifelse(bmi.group==5, ocever, 0)) %>%
  
  # Interaction between duration of OC use and BMI
  mutate(ocdur1.bmigroup3=ifelse(ocdur.group==1 & bmi.group==3, 1, 0)) %>%
  mutate(ocdur1.bmigroup4=ifelse(ocdur.group==1 & bmi.group==4, 1, 0)) %>%
  mutate(ocdur1.bmigroup5=ifelse(ocdur.group==1 & bmi.group==5, 1, 0)) %>%
  mutate(ocdur2.bmigroup3=ifelse(ocdur.group==2 & bmi.group==3, 1, 0)) %>%
  mutate(ocdur2.bmigroup4=ifelse(ocdur.group==2 & bmi.group==4, 1, 0)) %>%
  mutate(ocdur2.bmigroup5=ifelse(ocdur.group==2 & bmi.group==5, 1, 0)) %>%
  mutate(ocdur3.bmigroup3=ifelse(ocdur.group==3 & bmi.group==3, 1, 0)) %>%
  mutate(ocdur3.bmigroup4=ifelse(ocdur.group==3 & bmi.group==4, 1, 0)) %>%
  mutate(ocdur3.bmigroup5=ifelse(ocdur.group==3 & bmi.group==5, 1, 0)) %>%

  # Interaction between any HT use and BMI 
  mutate(anyhrtever.bmigroup3=ifelse(bmi.group==3, anyhrtever, 0)) %>%
  mutate(anyhrtever.bmigroup4=ifelse(bmi.group==4, anyhrtever, 0)) %>%
  mutate(anyhrtever.bmigroup5=ifelse(bmi.group==5, anyhrtever, 0)) %>%

  # Interaction between e-only HT use and BMI
  mutate(ehrtever.bmigroup3=ifelse(bmi.group==3 & ehrtever==1, 1, 0)) %>%
  mutate(ehrtever.bmigroup4=ifelse(bmi.group==4 & ehrtever==1, 1, 0)) %>%
  mutate(ehrtever.bmigroup5=ifelse(bmi.group==5 & ehrtever==1, 1, 0)) %>%

  # Interaction between EP HT use and BMI
  mutate(ephrtever.bmigroup3=ifelse(bmi.group==3 & ephrtever==1, 1, 0)) %>%
  mutate(ephrtever.bmigroup4=ifelse(bmi.group==4 & ephrtever==1, 1, 0)) %>%
  mutate(ephrtever.bmigroup5=ifelse(bmi.group==5 & ephrtever==1, 1, 0)) %>%

  # Interaction between duration of e-only HT use and BMI
  mutate(ehrtdur1.bmigroup3=ifelse(ehrtdur.group==1 & bmi.group==3, 1, 0)) %>%
  mutate(ehrtdur1.bmigroup4=ifelse(ehrtdur.group==1 & (bmi.group==4|bmi.group==5), 1, 0)) %>%
  mutate(ehrtdur2.bmigroup3=ifelse(ehrtdur.group==2 & bmi.group==3, 1, 0)) %>%
  mutate(ehrtdur2.bmigroup4=ifelse(ehrtdur.group==2 & (bmi.group==4|bmi.group==5), 1, 0)) %>%
  mutate(ehrtdur3.bmigroup3=ifelse(ehrtdur.group==3 & bmi.group==3, 1, 0)) %>%
  mutate(ehrtdur3.bmigroup4=ifelse(ehrtdur.group==3 & (bmi.group==4|bmi.group==5), 1, 0)) %>%

  # Interaction between duration of EP HT use and BMI
  mutate(ephrtdur1.bmigroup3=ifelse(ephrtdur.group==1 & bmi.group==3, 1, 0)) %>%
  mutate(ephrtdur1.bmigroup4=ifelse(ephrtdur.group==1 & (bmi.group==4|bmi.group==5), 1, 0)) %>%
  mutate(ephrtdur2.bmigroup3=ifelse(ephrtdur.group==2 & bmi.group==3, 1, 0)) %>%
  mutate(ephrtdur2.bmigroup4=ifelse(ephrtdur.group==2 & (bmi.group==4|bmi.group==5), 1, 0)) %>%
  mutate(ephrtdur3.bmigroup3=ifelse(ephrtdur.group==3 & bmi.group==3, 1, 0)) %>%
  mutate(ephrtdur3.bmigroup4=ifelse(ephrtdur.group==3 & (bmi.group==4|bmi.group==5), 1, 0)) %>%
 
  # Making sure variable names match up and finalizing time-related variables
  rename(diabetes1=diabetes) %>%
  rename(hypertension1=hypertension) %>%
  rename(ehrtever1=ehrtever) %>%
  rename(ephrtever1=ephrtever) %>%
  mutate(birth.year=1900+yobf+mobf/12) %>%
  mutate(study.entry.age=(start.year-birth.year)) %>%
  mutate(study.exit.age=ifelse(!is.na(endo.year) & end.year==endo.year, end.year-birth.year+10, end.year-birth.year)) %>%
  mutate(observed.outcome=endo.ever) %>%
  mutate(time.of.onset=ifelse(!is.na(endo.year), endo.year-birth.year-study.entry.age, Inf)) %>%
  mutate(observed.followup=study.exit.age-study.entry.age) %>%
  mutate(study.entry.age=floor(study.entry.age)) %>%
  mutate(study.exit.age=ceil(study.exit.age)) %>%  
  mutate(ephrtdur.group3=0) %>%
  mutate(gwas=ifelse(is.na(gwas), 0, 1)) %>%
  filter(time.of.onset>=0.5)
  

# ---------------------------------------------------------------
# ------ (5) NHS II: Importing/cleaning questionnaire data ------
# ---------------------------------------------------------------

# --- a. Importing data ---

nhsii.import <- read.csv("Input/E2C2_Validation_2b_ExportedData_NHSII.csv")


# --- b. Determining exclusions ---

nhsii.all <- nhsii.import %>%
  arrange(id, cancer_dxmonth) %>%
  distinct(id, .keep_all=T) %>%  

  # Calculating year of cancer diagnosis/death
  mutate(endo.year=endo_dxmonth/12+1900) %>%
  mutate(death.year=(deadmonth/12)+1900) %>%
  mutate(cancer.year=(cancer_dxmonth/12+1900)) %>%
  mutate(cancer.year=ifelse((is.na(cancer.year)|endo.year<cancer.year) & 
    !is.na(endo.year), endo.year, cancer.year)) %>%
  
  # Determining exclusions at baseline
  mutate(in.all=1) %>%
  mutate(ex.cancer=ifelse(!is.na(cancer.year) & cancer.year<=1989, cancer.year, 9999)) %>%
  mutate(ex.hyst=ifelse(!is.na(hyst89) & hyst89==1, 1989, 9999)) %>%
  mutate(ex.date=pmap(list(ex.cancer, ex.hyst), min)) %>%
  mutate(ex.cancer=ifelse(ex.cancer!=9999 & ex.cancer==ex.date, 1, 0)) %>%
  mutate(ex.hyst=ifelse(ex.hyst!=9999 & ex.hyst==ex.date, 1, 0)) %>%
  mutate(ex.hyst=ifelse(ex.cancer==1, 0, ex.hyst)) %>%
  rowwise() %>%
  mutate(menopause.age = median(c(!!! rlang::syms(grep('menoage', names(.), value=TRUE))), na.rm=T)) %>%
  ungroup() %>%
  mutate(ex.white=ifelse(!is.na(race) & race!=1 & ex.hyst==0 & ex.cancer==0, 1, 0)) %>%
  mutate(ex=ifelse(ex.hyst==1|ex.cancer==1|ex.white==1, 1, 0)) %>%
  mutate(include=1-ex)

nhsii.inclusion <- nhsii.all %>% summarize_at(vars(in.all, ex, ex.cancer, 
  ex.hyst, ex.white, include), list(~sum(., na.rm=T)))


# --- c. Converting Data from Wide to Long ---

nhsii.long <- nhsii.all %>%
  
  # Removing people excluded at baseline
  filter(include==1) %>%
  
  # Converting data from wide to long
  pivot_longer(cols=c(starts_with("bmi"), starts_with("hyst"), starts_with("menstat"),
    starts_with("anyhrtever"), starts_with("ephrtever"), starts_with("ephrtdur"),
    starts_with("ehrtever"), starts_with("ehrtdur"), starts_with("smoking"), 
    starts_with("packyears"), starts_with("ocever"), starts_with("ocdur"), 
    starts_with("parity"), starts_with("hypertension"), starts_with("alcoholgr"), 
    starts_with("retmo")), 
    names_to=c(".value", "year"), names_sep="(?<=[A-Za-z])(?=[0-9])") %>%
  
  # Creating variable for time (year/age)
  mutate(year=ifelse(as.numeric(year)<20, as.numeric(year)+2000, as.numeric(year)+1900)) %>%
  mutate(yob=birthday/12+1900) %>%
  mutate(age=year-yob) %>%
  
  # Determining timing of events or meeting exclusion criteria over follow-up
  mutate(endo.ever=ifelse(!is.na(endo.year), 1, 0)) %>%
  rename(completed=retmo) %>%
  mutate(hyst=ifelse(!is.na(hyst), year, 9999)) %>%
  group_by(id) %>%
  mutate(hyst.year=min(hyst)) %>%
  ungroup() %>%
  mutate(include=ifelse(!is.na(age) & !is.na(menopause.age) & age>=menopause.age 
    & !is.na(completed) & age>start.age, year, 9999)) %>%
  mutate(exclude=ifelse(!is.na(completed), year+2, 0000),
         exclude=ifelse(age>85, year, exclude),
         exclude=ifelse(!is.na(death.year) & death.year<exclude, death.year, exclude),
         exclude=ifelse(hyst.year<exclude, hyst.year, exclude),
         exclude=ifelse(!is.na(cancer.year) & cancer.year<exclude, cancer.year, exclude),
         exclude=ifelse(exclude>=2014, 2014, exclude)) %>%
  
  # Calculating years for start and end of follow-up
  group_by(id) %>%
  mutate(start.year=min(include)) %>%
  mutate(end.year=max(exclude)) %>%
  ungroup() %>%
  mutate(start.year=ifelse(start.year>=end.year, 9999, start.year)) %>%
  
  # Filling in missing values for parity, ocever and smoking (carry last observation forward)
  group_by(id) %>%
  fill(parity, ocever, smoking) %>%
  ungroup() %>%
  group_by(id) %>%
  fill(parity, ocever, smoking) %>%
  ungroup() %>%
  
  # Cleaning/recoding variables
  mutate(education=3) %>%
  mutate(education1=0) %>%
  mutate(education2=0) %>%
  mutate(education3=1) %>%
  mutate(ehrtdur=ehrtdur/12) %>%
  mutate(ephrtdur=ephrtdur/12) %>%
  mutate(ocdur=ocdur/12) %>%  
  mutate(famhx=famhx01) %>%
  mutate(diabetes=ifelse(!is.na(db_dxmonth) & (db_dxmonth/12+1900)<=year, 1, 0)) %>%
  mutate(hypertension=ifelse(is.na(hypertension) & !is.na(completed), 0, hypertension),
         hypertension=ifelse(hypertension==2, 0, hypertension)) %>%
  mutate(ehrtever=ifelse(anyhrtever==0, 0, ehrtever)) %>%
  mutate(ephrtever=ifelse(anyhrtever==0, 0, ephrtever)) %>%  
  mutate(anyhrtever=ifelse(is.na(anyhrtever) & lag(anyhrtever)==1 & id==lag(id), 1, anyhrtever)) %>%
  mutate(ehrtever=ifelse(is.na(ehrtever) & lag(ehrtever)==1 & id==lag(id), 1, ehrtever)) %>%
  mutate(ephrtever=ifelse(is.na(ephrtever) & lag(ephrtever)==1 & id==lag(id), 1, ephrtever)) %>%
  arrange(id, -year) %>%
  mutate(anyhrtever=ifelse(is.na(anyhrtever) & lag(anyhrtever)==0 & id==lag(id), 0, anyhrtever)) %>%
  arrange(id, -year) %>%
  mutate(ehrtever=ifelse(is.na(ehrtever) & lag(ehrtever)==0 & id==lag(id), 0, ehrtever)) %>%
  mutate(ephrtever=ifelse(is.na(ephrtever) & lag(ephrtever)==0 & id==lag(id), 0, ephrtever)) %>%  
  arrange(id, year) %>%
  mutate(ocdur=ifelse(ocever==0, 0, ocdur)) %>%
  mutate(ehrtdur=ifelse(ehrtever==0, 0, ehrtdur)) %>%
  mutate(ephrtdur=ifelse(ephrtever==0, 0, ephrtdur)) %>%
  mutate(bmi.group=cut(bmi, breaks=c(-Inf, 18.5, 25, 30, 35, Inf), 
    right=F, labels=c(1, 2, 3, 4, 5))) %>%
  mutate(parity.group=ifelse(parity>=4, 4, parity)) %>%
  mutate(menarchage.group=ifelse(menarchage==0, NA, ifelse(menarchage<=10, 10, 
    ifelse(menarchage>=16, 16, menarchage)))) %>%  
  mutate(menarchage2.group=ifelse(menarchage==0, NA, ifelse(menarchage<=9, 9, 
    ifelse(menarchage %in% c(10, 11), 10, ifelse(menarchage %in% c(12, 13), 12,
    ifelse(menarchage %in% c(14, 15), 14, ifelse(menarchage >=16, 16, NA))))))) %>%  
  mutate(ephrtdur.group=ifelse(ephrtdur>0 & ephrtdur<=5, 1, ifelse(ephrtdur>5 & ephrtdur<=10, 
    2, ifelse(ephrtdur!=999 & ephrtdur!=0, 3, ephrtdur)))) %>%
  mutate(ehrtdur.group=ifelse(ehrtdur>0 & ehrtdur<=5, 1, ifelse(ehrtdur>5 & 
    ehrtdur<=10, 2, ifelse(ehrtdur!=999 & ehrtdur!=0, 3, ehrtdur)))) %>%
  mutate(ocdur.group=ifelse(ocdur>0 & ocdur<=5, 1, ifelse(ocdur>5 & ocdur<=10, 
    2, ifelse(ocdur!=999 & ocdur!=0, 3, ocdur)))) %>%
  mutate(firstbirthage.group=cut(firstbirthage, breaks=c(-Inf, 20, 25, 30, 35, Inf), 
    labels=seq(1:5), right=F)) %>%
  mutate(firstbirthage.group=ifelse(!is.na(parity.group) & parity.group==0, 9, 
    firstbirthage.group)) %>%
  dummy_cols(., select_columns=c("smoking", "bmi.group", "parity.group", 
    "menarchage.group", "menarchage2.group", "ehrtdur.group", "ephrtdur.group", 
    "ocdur.group", "firstbirthage.group"), ignore_na=T) 

names(nhsii.long) <- sub("_", "", names(nhsii.long))


# --- d. Creating interaction terms and finalizing dataset ---

nhsii.final <- nhsii.long %>%
  filter(year==start.year) %>%

  # Interaction between ever OC use and BMI
  mutate(ocever.bmigroup3=ifelse(bmi.group==3, ocever, 0)) %>%
  mutate(ocever.bmigroup4=ifelse(bmi.group==4, ocever, 0)) %>%
  mutate(ocever.bmigroup5=ifelse(bmi.group==5, ocever, 0)) %>%
  
  # Interaction between duration of OC use and BMI
  mutate(ocdur1.bmigroup3=ifelse(ocdur.group==1 & bmi.group==3, 1, 0)) %>%
  mutate(ocdur1.bmigroup4=ifelse(ocdur.group==1 & bmi.group==4, 1, 0)) %>%
  mutate(ocdur1.bmigroup5=ifelse(ocdur.group==1 & bmi.group==5, 1, 0)) %>%
  mutate(ocdur2.bmigroup3=ifelse(ocdur.group==2 & bmi.group==3, 1, 0)) %>%
  mutate(ocdur2.bmigroup4=ifelse(ocdur.group==2 & bmi.group==4, 1, 0)) %>%
  mutate(ocdur2.bmigroup5=ifelse(ocdur.group==2 & bmi.group==5, 1, 0)) %>%
  mutate(ocdur3.bmigroup3=ifelse(ocdur.group==3 & bmi.group==3, 1, 0)) %>%
  mutate(ocdur3.bmigroup4=ifelse(ocdur.group==3 & bmi.group==4, 1, 0)) %>%
  mutate(ocdur3.bmigroup5=ifelse(ocdur.group==3 & bmi.group==5, 1, 0)) %>%
  
  # Interaction between any HT use and BMI 
  mutate(anyhrtever.bmigroup3=ifelse(bmi.group==3, anyhrtever, 0)) %>%
  mutate(anyhrtever.bmigroup4=ifelse(bmi.group==4, anyhrtever, 0)) %>%
  mutate(anyhrtever.bmigroup5=ifelse(bmi.group==5, anyhrtever, 0)) %>%
  
  # Interaction between e-only HT use and BMI
  mutate(ehrtever.bmigroup3=ifelse(bmi.group==3 & ehrtever==1, 1, 0)) %>%
  mutate(ehrtever.bmigroup4=ifelse(bmi.group==4 & ehrtever==1, 1, 0)) %>%
  mutate(ehrtever.bmigroup5=ifelse(bmi.group==5 & ehrtever==1, 1, 0)) %>%
  
  # Interaction between EP HT use and BMI
  mutate(ephrtever.bmigroup3=ifelse(bmi.group==3 & ephrtever==1, 1, 0)) %>%
  mutate(ephrtever.bmigroup4=ifelse(bmi.group==4 & ephrtever==1, 1, 0)) %>%
  mutate(ephrtever.bmigroup5=ifelse(bmi.group==5 & ephrtever==1, 1, 0)) %>%
  
  # Interaction between duration of e-only HT use and BMI
  mutate(ehrtdur1.bmigroup3=ifelse(ehrtdur.group==1 & bmi.group==3, 1, 0)) %>%
  mutate(ehrtdur1.bmigroup4=ifelse(ehrtdur.group==1 & (bmi.group==4|bmi.group==5), 1, 0)) %>%
  mutate(ehrtdur2.bmigroup3=ifelse(ehrtdur.group==2 & bmi.group==3, 1, 0)) %>%
  mutate(ehrtdur2.bmigroup4=ifelse(ehrtdur.group==2 & (bmi.group==4|bmi.group==5), 1, 0)) %>%
  mutate(ehrtdur3.bmigroup3=ifelse(ehrtdur.group==3 & bmi.group==3, 1, 0)) %>%
  mutate(ehrtdur3.bmigroup4=ifelse(ehrtdur.group==3 & (bmi.group==4|bmi.group==5), 1, 0)) %>%
  
  # Interaction between duration of EP HT use and BMI
  mutate(ephrtdur1.bmigroup3=ifelse(ephrtdur.group==1 & bmi.group==3, 1, 0)) %>%
  mutate(ephrtdur1.bmigroup4=ifelse(ephrtdur.group==1 & (bmi.group==4|bmi.group==5), 1, 0)) %>%
  mutate(ephrtdur2.bmigroup3=ifelse(ephrtdur.group==2 & bmi.group==3, 1, 0)) %>%
  mutate(ephrtdur2.bmigroup4=ifelse(ephrtdur.group==2 & (bmi.group==4|bmi.group==5), 1, 0)) %>%
  mutate(ephrtdur3.bmigroup3=ifelse(ephrtdur.group==3 & bmi.group==3, 1, 0)) %>%
  mutate(ephrtdur3.bmigroup4=ifelse(ephrtdur.group==3 & (bmi.group==4|bmi.group==5), 1, 0)) %>%
  
  # Making sure variable names match up and finalizing time-related variables
  rename(diabetes1=diabetes) %>%
  rename(hypertension1=hypertension) %>%
  rename(ehrtever1=ehrtever) %>%
  rename(ephrtever1=ephrtever) %>%
  mutate(birth.year=1900+birthday/12) %>%
  mutate(study.entry.age=(start.year-birth.year)) %>%
  mutate(study.exit.age=ifelse(!is.na(endo.year) & end.year==endo.year, end.year-birth.year+10, end.year-birth.year)) %>%
  mutate(observed.outcome=endo.ever) %>%
  mutate(time.of.onset=ifelse(!is.na(endo.year), endo.year-birth.year-study.entry.age, Inf)) %>%
  mutate(observed.followup=study.exit.age-study.entry.age) %>%
  mutate(study.entry.age=floor(study.entry.age)) %>%
  mutate(study.exit.age=ceil(study.exit.age)) %>%
  filter(time.of.onset>=0.5)


# ---------------------------------------------------------------
# ---------- (6) PLCO: Importing/cleaning genetic data ----------
# ---------------------------------------------------------------

# --- a. Importing data ---

plco.chr1 <- lapply(Sys.glob("Input/PLCO SNPs/*-chr1-*.vcf.gz"), read.vcfR)
plco.chr2 <- lapply(Sys.glob("Input/PLCO SNPs/*-chr2-*.vcf.gz"), read.vcfR)
plco.chr6 <- lapply(Sys.glob("Input/PLCO SNPs/*-chr6-*.vcf.gz"), read.vcfR)
plco.chr8 <- lapply(Sys.glob("Input/PLCO SNPs/*-chr8-*.vcf.gz"), read.vcfR)
plco.chr9 <- lapply(Sys.glob("Input/PLCO SNPs/*-chr9-*.vcf.gz"), read.vcfR)
plco.chr11 <- lapply(Sys.glob("Input/PLCO SNPs/*-chr11-*.vcf.gz"), read.vcfR)
plco.chr12 <- lapply(Sys.glob("Input/PLCO SNPs/*-chr12-*.vcf.gz"), read.vcfR)
plco.chr13 <- lapply(Sys.glob("Input/PLCO SNPs/*-chr13-*.vcf.gz"), read.vcfR)
plco.chr14 <- lapply(Sys.glob("Input/PLCO SNPs/*-chr14-*.vcf.gz"), read.vcfR)
plco.chr15 <- lapply(Sys.glob("Input/PLCO SNPs/*-chr15-*.vcf.gz"), read.vcfR)
plco.chr17 <- lapply(Sys.glob("Input/PLCO SNPs/*-chr17-*.vcf.gz"), read.vcfR)


# --- b. Converting VCF to data frames

vcf2df <- function(i, j){data.frame(t(extract.gt(i[[j]], return.alleles=T, IDtoRowNames=T)))}

plco.gwas1 <- lapply(ls(pattern="plco.chr"), function(i){
  x <- get(i)
  return <- do.call(bind_rows, lapply(1:length(x), function(j){
    data.frame(t(extract.gt(x[[j]], return.alleles=T, IDtoRowNames=T)))
  }))
  return$plco_id <- rownames(return)
  return(return)
})


# --- c. Combining datasets ---

plco.gwas2 <- Reduce(function(x,y,...) merge(x,y,all=T,...), plco.gwas1) %>%
  mutate(plco_id=substr(plco_id, 1, regexpr("\\_", plco_id)-1)) 

colnames(plco.gwas2)[-1] <- paste("X", substr(colnames(plco.gwas2)[-1], 4, 
  nchar(colnames(plco.gwas2)[-1])-4), sep="")


# d. Renaming variables
for (i in snp.ref$SNP38){
  plco.gwas2[,i] <- str_count(plco.gwas2[,i], snp.ref[snp.ref$SNP38==i, 
    "Alternate.Allele"])
  colnames(plco.gwas2)[colnames(plco.gwas2)==i] <- snp.ref[snp.ref$SNP38==i, "rsID"]
}

plco.genetic <- plco.gwas2 %>% mutate(genetic=1)


# ---------------------------------------------------------------
# ------- (7) PLCO: Importing/cleaning questionnaire data -------
# ---------------------------------------------------------------

# --- a. Importing data ---
plco.import <- read_sas("/proj/nhgens/nhgen0z/MODELING/ENDOMETRIAL_CANCER/VALIDATION_PLCO/plco182042916.sas7bdat")


# --- b. Determining exclusions at baseline ---

plco.all <- plco.import %>% 
  left_join(plco.genetic, by="plco_id") %>%
  mutate(ex.questionnaire=ifelse(bq_returned==0, 1, 0)) %>%
  mutate(ex.age=ifelse(bq_cohort_entryage>75|is.na(bq_cohort_entryage), 1, 0)) %>%
  mutate(ex.age=ifelse(ex.questionnaire==1, 0, ex.age)) %>%
  mutate(ex.cancer=ifelse(trial_ph_any==1|trial_ph_any==9, 1, 0)) %>%
  mutate(ex.cancer=ifelse(ex.questionnaire==1|ex.age==1, 0, ex.cancer)) %>%
  mutate(ex.hyst=ifelse(hyster_f==1|hyster_f==2|is.na(hyster_f), 1, 0)) %>%
  mutate(ex.hyst=ifelse(ex.questionnaire==1|ex.age==1|ex.cancer==1, 0, ex.hyst)) %>% 
  mutate(ex.race=ifelse(race7!=1|is.na(race7), 1, 0)) %>%
  mutate(ex.race=ifelse(ex.questionnaire==1|ex.age==1|ex.cancer==1|ex.hyst==1, 0, ex.race)) %>%  
  mutate(ex=ifelse(ex.questionnaire==1|ex.age==1|ex.race==1|ex.cancer==1|ex.hyst==1, 1, 0)) %>%
  mutate(include=1-ex) %>%
  mutate(all=1)

plco.inclusion <- plco.all %>% summarize_at(vars(all, ex, starts_with("ex"), include), sum)


# --- c. Cleaning variables, part 1 ---

plco.long <- plco.all %>%
  
  # Removing people excluded at baseline
  filter(include==1) %>%
  
  # Creaing variable for time (year/age)
  mutate(age=bq_cohort_entryage) %>%
  mutate(year=rndyear) %>%
  
  # Cleaning/renaming covariates 
  mutate(education=ifelse(educat %in% c(1,2,3), 1, ifelse(educat %in% c(4,5), 2, 
    ifelse(educat %in% c(6, 7), 3, NA)))) %>%
  mutate(smoking=ifelse(cig_stat==0, 1, ifelse(cig_stat==1, 3, 
    ifelse(cig_stat==2, 2, NA)))) %>%
  rename(packyears=pack_years) %>%
  mutate(bmi=bmi_curr) %>%
  mutate(parity=livec+stillb) %>%
  mutate(parity=ifelse(is.na(stillb), livec, parity)) %>%
  mutate(firstbirthage.group=ifelse(fchilda %in% c(1,2), 1, fchilda-1)) %>%
  mutate(firstbirthage.group=ifelse(firstbirthage.group %in% c(5, 6), 5, firstbirthage.group)) %>%
  mutate(firstbirthage=ifelse(fchilda==1, 14, ifelse(fchilda==2, 17.5, ifelse(fchilda==3, 22.5, 
    ifelse(fchilda==4, 27.5, ifelse(fchilda==5, 32.5, ifelse(fchilda==6, 37.5, 
    ifelse(fchilda==7, 42.5, NA)))))))) %>%
  mutate(menarchage.group=ifelse(fmenstr==1, 10, ifelse(fmenstr==2, 11, ifelse(fmenstr==3, 13,
    ifelse(fmenstr==4, 15, ifelse(fmenstr==5, 16, NA)))))) %>%
  mutate(menarchage2.group=ifelse(fmenstr==1, 9, ifelse(fmenstr==2, 10,
    ifelse(fmenstr==3, 12, ifelse(fmenstr==4, 14, ifelse(fmenstr==5, 16, NA)))))) %>%
  mutate(anyhrtever=sqxbq_hrt) %>%
  mutate(hrt.agestop1=ifelse(!is.na(sqxbq_hrt_agestop), sqxbq_hrt_agestop, sqxbq_hrt_ageswitch)) %>%
  mutate(hrt.agestop1=ifelse(!is.na(sqxbq_hrt_chng) & sqxbq_hrt_chng==1, age, hrt.agestop1)) %>%
  mutate(hrt.agestop1=ifelse(hrt.agestop1>age, age, hrt.agestop1)) %>%
  mutate(hrtdur1=hrt.agestop1-sqxbq_hrt_age) %>%
  mutate(hrtdur1=ifelse(hrtdur1==0, 0.5, hrtdur1)) %>%
  mutate(hrtdur1=ifelse(hrtdur1<0, 0, hrtdur1)) %>%
  mutate(hrtdur2=ifelse(sqxbq_hrt_age==age, 0.5, hrtdur1)) %>%
  mutate(hrtdur1=ifelse(sqxbq_hrt_age>age, 0, hrtdur1)) %>%
  mutate(hrt.agestop2=ifelse(sqxbq_hrt2curr==1, age, NA)) %>%
  mutate(hrtdur2=hrt.agestop2-sqxbq_hrt_ageswitch) %>%
  mutate(hrtdur2=ifelse(hrtdur2==0, 0.5, hrtdur2)) %>%
  mutate(hrtdur2=ifelse(hrtdur2<0, 0, hrtdur2)) %>%
  mutate(hrtdur2=ifelse(is.na(hrtdur2) & is.na(sqxbq_hrt_type2) & 
    is.na(sqxbq_hrt_ageswitch) & is.na(sqxbq_hrt2curr), 0, hrtdur2)) %>%
  mutate(hrtdur2=ifelse(sqxbq_hrt_ageswitch==age, 0.5, hrtdur2)) %>%
  mutate(hrtdur2=ifelse(sqxbq_hrt_ageswitch>age, 0, hrtdur2)) %>%
  mutate(ehrtdur=ifelse(sqxbq_hrt_type1==1 & (is.na(sqxbq_hrt_type2)|sqxbq_hrt_type2!=1), hrtdur1,
    ifelse(sqxbq_hrt_type1!=1 & !is.na(sqxbq_hrt_type2) & sqxbq_hrt_type2==1, hrtdur2,
    ifelse(sqxbq_hrt_type1==1 & !is.na(sqxbq_hrt_type2) & sqxbq_hrt_type2==1, hrtdur1+hrtdur2, NA)))) %>%
  mutate(ephrtdur=ifelse(sqxbq_hrt_type1==3 & (is.na(sqxbq_hrt_type2)|sqxbq_hrt_type2!=3), hrtdur1,
    ifelse(sqxbq_hrt_type1!=3 & !is.na(sqxbq_hrt_type2) & sqxbq_hrt_type2==3, hrtdur2,
    ifelse(sqxbq_hrt_type1==3 & !is.na(sqxbq_hrt_type2) & sqxbq_hrt_type2==3, hrtdur1+hrtdur2, NA)))) %>%
  mutate(ehrtever=ifelse(anyhrtever==0, 0, NA)) %>%
  mutate(ehrtever=ifelse(is.na(ehrtever) & !is.na(sqxbq_hrt_type1) & 
    sqxbq_hrt_type1!=1 & !is.na(sqxbq_hrt_type2) & sqxbq_hrt_type2!=1, 0, ehrtever)) %>%
  mutate(ehrtever=ifelse(is.na(ehrtever) & !is.na(sqxbq_hrt_type1) & 
    sqxbq_hrt_type1!=1 & is.na(sqxbq_hrt_ageswitch) & is.na(sqxbq_hrt_type2), 0, ehrtever)) %>%
  mutate(ehrtever=ifelse(is.na(ehrtever) & ehrtdur>0, 1, ehrtever)) %>%
  mutate(ehrtever=ifelse(is.na(ehrtever) & ehrtdur==0, 0, ehrtever)) %>%
  mutate(ehrtever=ifelse(is.na(ehrtever) & ((!is.na(sqxbq_hrt_type1) & 
    sqxbq_hrt_type1==1 & sqxbq_hrt_age<=age)|(!is.na(sqxbq_hrt_type2) & 
    sqxbq_hrt_type2==1 & sqxbq_hrt_ageswitch<=age)), 1, ehrtever)) %>%
  mutate(ehrtdur=ifelse(ehrtever==0, 0, ehrtdur)) %>%
  mutate(ephrtever=ifelse(anyhrtever==0, 0, NA)) %>%
  mutate(ephrtever=ifelse(is.na(ephrtever) & !is.na(sqxbq_hrt_type1) & 
    sqxbq_hrt_type1!=3 & !is.na(sqxbq_hrt_type2) & sqxbq_hrt_type2!=3, 0, ephrtever)) %>%
  mutate(ephrtever=ifelse(is.na(ephrtever) & !is.na(sqxbq_hrt_type1) & 
    sqxbq_hrt_type1!=3 & is.na(sqxbq_hrt_ageswitch) & is.na(sqxbq_hrt_type2), 0, ephrtever)) %>%
  mutate(ephrtever=ifelse(is.na(ephrtever) & ephrtdur>0, 1, ephrtever)) %>%
  mutate(ephrtever=ifelse(is.na(ephrtever) & ephrtdur==0, 0, ephrtever)) %>%
  mutate(ephrtever=ifelse(is.na(ephrtever) & ((!is.na(sqxbq_hrt_type1) & 
    sqxbq_hrt_type1==3 & sqxbq_hrt_age<=age)|(!is.na(sqxbq_hrt_type2) & 
    sqxbq_hrt_type2==3 & sqxbq_hrt_ageswitch<=age)), 1, ephrtever)) %>%
  mutate(ephrtdur=ifelse(ephrtever==0, 0, ephrtdur)) %>%
  mutate(ocever=bcontr_f) %>%
  mutate(ocdur.group=ifelse(bcontrt==0, 0, ifelse(bcontrt %in% c(3, 4, 5), 1, 
    ifelse(bcontrt==2, 2, ifelse(bcontrt==1, 3, NA))))) %>%
  mutate(diabetes1=diabetes_f) %>%
  mutate(hypertension1=hyperten_f) %>%
  mutate(bmi.group=cut(bmi, breaks=c(-Inf, 18.5, 25, 30, 35, Inf), right=F, 
    labels=c(1, 2, 3, 4, 5))) %>%
  mutate(parity.group=ifelse(parity>=4, 4, parity)) %>%
  mutate(ehrtdur.group=ifelse(ehrtdur==0, 0, ifelse(ehrtdur>0 & ehrtdur<=5, 1, 
    ifelse(ehrtdur>5 & ehrtdur<=10, 2, ifelse(ehrtdur>10, 3, NA))))) %>%
  mutate(ephrtdur.group=ifelse(ephrtdur==0, 0, ifelse(ephrtdur>0 & ephrtdur<=5, 
    1, ifelse(ephrtdur>5 & ephrtdur<=10, 2, ifelse(ephrtdur>10, 3, NA))))) %>%
  mutate(menarchage.group12=0) %>%
  mutate(menarchage.group14=0) %>%
  mutate(firstbirthage=ifelse(!is.na(parity.group) & parity.group==0, 0, firstbirthage)) %>%  
  mutate(firstbirthage.group=ifelse(!is.na(parity.group) & parity.group==0, 9, firstbirthage.group)) %>%  
  dummy_cols(., select_columns=c("education", "smoking", "bmi.group", 
    "parity.group", "menarchage.group", "menarchage2.group", "ehrtdur.group", 
    "ephrtdur.group", "ocdur.group", "firstbirthage.group"), ignore_na=T) 

names(plco.long) <- sub("_", "", names(plco.long))  


# --- d. Cleaning variables, part 2 ---

plco.long2 <- plco.long %>%
  
  # Interaction between ever OC use and BMI
  mutate(ocever.bmigroup3=ifelse(bmi.group==3, ocever, 0)) %>%
  mutate(ocever.bmigroup4=ifelse(bmi.group==4, ocever, 0)) %>%
  mutate(ocever.bmigroup5=ifelse(bmi.group==5, ocever, 0)) %>%

  # Interaction between duration of OC use and BMI
  mutate(ocdur1.bmigroup3=ifelse(ocdur.group==1 & bmi.group==3, 1, 0)) %>%
  mutate(ocdur1.bmigroup4=ifelse(ocdur.group==1 & bmi.group==4, 1, 0)) %>%
  mutate(ocdur1.bmigroup5=ifelse(ocdur.group==1 & bmi.group==5, 1, 0)) %>%
  mutate(ocdur2.bmigroup3=ifelse(ocdur.group==2 & bmi.group==3, 1, 0)) %>%
  mutate(ocdur2.bmigroup4=ifelse(ocdur.group==2 & bmi.group==4, 1, 0)) %>%
  mutate(ocdur2.bmigroup5=ifelse(ocdur.group==2 & bmi.group==5, 1, 0)) %>%
  mutate(ocdur3.bmigroup3=ifelse(ocdur.group==3 & bmi.group==3, 1, 0)) %>%
  mutate(ocdur3.bmigroup4=ifelse(ocdur.group==3 & bmi.group==4, 1, 0)) %>%
  mutate(ocdur3.bmigroup5=ifelse(ocdur.group==3 & bmi.group==5, 1, 0)) %>%
  
  # Interaction between any HT use and BMI 
  mutate(anyhrtever.bmigroup3=ifelse(bmi.group==3, anyhrtever, 0)) %>%
  mutate(anyhrtever.bmigroup4=ifelse(bmi.group==4, anyhrtever, 0)) %>%
  mutate(anyhrtever.bmigroup5=ifelse(bmi.group==5, anyhrtever, 0)) %>%
  
  # Interaction between e-only HT use and BMI
  mutate(ehrtever.bmigroup3=ifelse(bmi.group==3 & ehrtever==1, 1, 0)) %>%
  mutate(ehrtever.bmigroup4=ifelse(bmi.group==4 & ehrtever==1, 1, 0)) %>%
  mutate(ehrtever.bmigroup5=ifelse(bmi.group==5 & ehrtever==1, 1, 0)) %>%
  
  # Interaction between EP HT use and BMI
  mutate(ephrtever.bmigroup3=ifelse(bmi.group==3 & ephrtever==1, 1, 0)) %>%
  mutate(ephrtever.bmigroup4=ifelse(bmi.group==4 & ephrtever==1, 1, 0)) %>%
  mutate(ephrtever.bmigroup5=ifelse(bmi.group==5 & ephrtever==1, 1, 0)) %>%
  
  # Interaction between duration of e-only HT use and BMI
  mutate(ehrtdur1.bmigroup3=ifelse(ehrtdur.group==1 & bmi.group==3, 1, 0)) %>%
  mutate(ehrtdur1.bmigroup4=ifelse(ehrtdur.group==1 & (bmi.group==4|bmi.group==5), 1, 0)) %>%
  mutate(ehrtdur2.bmigroup3=ifelse(ehrtdur.group==2 & bmi.group==3, 1, 0)) %>%
  mutate(ehrtdur2.bmigroup4=ifelse(ehrtdur.group==2 & (bmi.group==4|bmi.group==5), 1, 0)) %>%
  mutate(ehrtdur3.bmigroup3=ifelse(ehrtdur.group==3 & bmi.group==3, 1, 0)) %>%
  mutate(ehrtdur3.bmigroup4=ifelse(ehrtdur.group==3 & (bmi.group==4|bmi.group==5), 1, 0)) %>%
  
  # Interaction between duration of EP HT use and BMI
  mutate(ephrtdur1.bmigroup3=ifelse(ephrtdur.group==1 & bmi.group==3, 1, 0)) %>%
  mutate(ephrtdur1.bmigroup4=ifelse(ephrtdur.group==1 & (bmi.group==4|bmi.group==5), 1, 0)) %>%
  mutate(ephrtdur2.bmigroup3=ifelse(ephrtdur.group==2 & bmi.group==3, 1, 0)) %>%
  mutate(ephrtdur2.bmigroup4=ifelse(ephrtdur.group==2 & (bmi.group==4|bmi.group==5), 1, 0)) %>%
  mutate(ephrtdur3.bmigroup3=ifelse(ephrtdur.group==3 & bmi.group==3, 1, 0)) %>%
  mutate(ephrtdur3.bmigroup4=ifelse(ephrtdur.group==3 & (bmi.group==4|bmi.group==5), 1, 0)) %>%
  
  # Cleaning remaining variables
  rename(ehrtever1=ehrtever) %>%
  rename(ephrtever1=ephrtever) %>%
  mutate(hyst.age=ifelse(sqxbqhystera==1, 37.5, ifelse(sqxbqhystera==2, 42.5, 
    ifelse(sqxbqhystera==3, 47.5, ifelse(sqxbqhystera==4, 52.5,
    ifelse(sqxbqhystera==5, 57.5, ifelse(sqxbqhystera==6, 65, 
    ifelse(sqxbqhystera==7, 75, ifelse(sqxbqhystera==8, 85, NA))))))))) %>%  
  mutate(hyst.age=ifelse(is.na(hyst.age), 99, hyst.age)) %>%
  mutate(study.entry.age=age) %>%
  mutate(exit.days=ifelse(fstcanexitdays-bqcohort_entrydays>=endoexitdays-bqcohort_entrydays, 
    endoexitdays-bqcohort_entrydays, fstcanexitdays-bqcohort_entrydays)) %>%
  mutate(observed.outcome=ifelse(endocstatus_cat==1 & endois_first_dx==1 & 
    exit.days==endoexitdays-bqcohort_entrydays, 1, 0)) %>%
  mutate(exit.days=ifelse(exit.days>=(hyst.age-study.entry.age)*365.25 & 
    observed.outcome==0, (hyst.age-study.entry.age)*365.25, exit.days)) %>%
  mutate(time.of.onset=ifelse(observed.outcome==1, exit.days/365.25, Inf)) %>%
  mutate(observed.followup=ifelse(observed.outcome==0, exit.days/365.25, 10)) %>%
  mutate(study.exit.age=ceil(study.entry.age+observed.followup))

plco.final <- plco.long2 %>% filter(observed.followup>0) %>% filter(time.of.onset>=0.5)

# ---------------------------------------------------------------
# -------- (8) Current: Creating dataset from NHANES data -------
# ---------------------------------------------------------------

current.final <-nhanes.current %>%
  mutate(anyhrtever=as.numeric(as.character(anyhrtever))) %>%
  mutate(ehrtever1=as.numeric(as.character(ehrtever1))) %>%
  mutate(ephrtever1=as.numeric(as.character(ephrtever1))) %>%
  mutate(ocever=as.numeric(as.character(ocever))) %>%
  mutate(diabetes1=as.numeric(as.character(diabetes1))) %>%
  mutate(hypertension1=as.numeric(as.character(hypertension1))) %>%
  mutate(study.entry.age=age) %>%
  mutate(observed.outcome=rbinom(nrow(nhanes.current), 1, p=0.4)) %>%
  mutate(random.time=runif(nrow(nhanes.current))) %>%
  mutate(time.of.onset=ifelse(observed.outcome==1, random.time*10, Inf)) %>%
  mutate(study.exit.age=age+20)


# ---------------------------------------------------------------
# ---------------------- (9) Group LASSO ------------------------
# ---------------------------------------------------------------


# --- a. Creating list of predictors for the model ---

# Covariate List
endo.cov.info <- list(
  list(name="education2", type="continuous"),
  list(name="education3", type="continuous"),
  list(name="smoking2", type="continuous"),
  list(name="smoking3", type="continuous"),
  list(name="bmi.group1", type="continuous"),
  list(name="bmi.group3", type="continuous"),
  list(name="bmi.group4", type="continuous"),
  list(name="bmi.group5", type="continuous"),
  list(name="parity.group1", type="continuous"),
  list(name="parity.group2", type="continuous"),
  list(name="parity.group3", type="continuous"),
  list(name="parity.group4", type="continuous"),
  list(name="menarchage2.group10", type="continuous"),
  list(name="menarchage2.group12", type="continuous"),
  list(name="menarchage2.group14", type="continuous"),
  list(name="menarchage2.group16", type="continuous"),
  list(name="anyhrtever", type="continuous"),
  list(name="ehrtever1", type="continuous"),  
  list(name="ehrtdur.group1", type="continuous"),  
  list(name="ehrtdur.group2", type="continuous"), 
  list(name="ehrtdur.group3", type="continuous"), 
  list(name="ephrtever1", type="continuous"),
  list(name="ephrtdur.group1", type="continuous"),   
  list(name="ephrtdur.group2", type="continuous"),  
  list(name="ephrtdur.group3", type="continuous"),    
  list(name="ocever", type="continuous"),
  list(name="ocdur.group1", type="continuous"),
  list(name="ocdur.group2", type="continuous"),
  list(name="ocdur.group3", type="continuous"),
  list(name="firstbirthage.group2", type="continuous"),
  list(name="firstbirthage.group3", type="continuous"),
  list(name="firstbirthage.group4", type="continuous"),
  list(name="firstbirthage.group5", type="continuous"),
  list(name="firstbirthage.group9", type="continuous"),
  list(name="diabetes1", type="continuous"),
  list(name="hypertension1", type="continuous"),
  list(name="ocever.bmigroup3", type="continuous"),
  list(name="ocever.bmigroup4", type="continuous"),
  list(name="ocever.bmigroup5", type="continuous"),
  list(name="ocdur1.bmigroup3", type="continuous"),
  list(name="ocdur1.bmigroup4", type="continuous"),
  list(name="ocdur1.bmigroup5", type="continuous"),
  list(name="ocdur2.bmigroup3", type="continuous"),
  list(name="ocdur2.bmigroup4", type="continuous"),
  list(name="ocdur2.bmigroup5", type="continuous"),
  list(name="ocdur3.bmigroup3", type="continuous"),
  list(name="ocdur3.bmigroup4", type="continuous"),
  list(name="ocdur3.bmigroup5", type="continuous"),
  list(name="anyhrtever.bmigroup3", type="continuous"),
  list(name="anyhrtever.bmigroup4", type="continuous"),
  list(name="anyhrtever.bmigroup5", type="continuous"),
  list(name="ehrtever.bmigroup3", type="continuous"),
  list(name="ehrtever.bmigroup4", type="continuous"),
  list(name="ehrtever.bmigroup5", type="continuous"),
  list(name="ehrtdur1.bmigroup3", type="continuous"),
  list(name="ehrtdur1.bmigroup4", type="continuous"),
  list(name="ehrtdur2.bmigroup3", type="continuous"),
  list(name="ehrtdur2.bmigroup4", type="continuous"),
  list(name="ehrtdur3.bmigroup3", type="continuous"),
  list(name="ehrtdur3.bmigroup4", type="continuous"),
  list(name="ephrtever.bmigroup3", type="continuous"),
  list(name="ephrtever.bmigroup4", type="continuous"),
  list(name="ephrtever.bmigroup5", type="continuous"),
  list(name="ephrtdur1.bmigroup3", type="continuous"),
  list(name="ephrtdur1.bmigroup4", type="continuous"),
  list(name="ephrtdur2.bmigroup3", type="continuous"),
  list(name="ephrtdur2.bmigroup4", type="continuous"),
  list(name="ephrtdur3.bmigroup3", type="continuous"),
  list(name="ephrtdur3.bmigroup4", type="continuous")
)

# Formula
endo.formula <- Y ~ education2 + education3 + smoking2 + smoking3 + 
  bmi.group1 + bmi.group3 + bmi.group4 + bmi.group5 + 
  parity.group1 + parity.group2 + parity.group3 + parity.group4 +
  menarchage2.group10 + menarchage2.group12 + menarchage2.group14 + menarchage2.group16 +
  anyhrtever + 
  ehrtever1 + ehrtdur.group1 + ehrtdur.group2 + ehrtdur.group3 + 
  ephrtever1 + ephrtdur.group1 + ephrtdur.group2 + ephrtdur.group3 + 
  ocever + ocdur.group1 + ocdur.group2 + ocdur.group3 + 
  firstbirthage.group2 + firstbirthage.group3 + firstbirthage.group4 + firstbirthage.group5 + firstbirthage.group9 +
  diabetes1 + hypertension1 + 
  ocever.bmigroup3 + ocever.bmigroup4 + ocever.bmigroup5 + 
  ocdur1.bmigroup3 + ocdur1.bmigroup4 + ocdur1.bmigroup5 +
  ocdur2.bmigroup3 + ocdur2.bmigroup4 + ocdur2.bmigroup5 +
  ocdur3.bmigroup3 + ocdur3.bmigroup4 + ocdur3.bmigroup5 +
  anyhrtever.bmigroup3 + anyhrtever.bmigroup4 + anyhrtever.bmigroup5 + 
  ehrtever.bmigroup3 + ehrtever.bmigroup4 + ehrtever.bmigroup5 + 
  ehrtdur1.bmigroup3 + ehrtdur1.bmigroup4 + 
  ehrtdur2.bmigroup3 + ehrtdur2.bmigroup4 + 
  ehrtdur3.bmigroup3 + ehrtdur3.bmigroup4 + 
  ephrtever.bmigroup3 + ephrtever.bmigroup4 + ephrtever.bmigroup5 + 
  ephrtdur1.bmigroup3 + ephrtdur1.bmigroup4 + 
  ephrtdur2.bmigroup3 + ephrtdur2.bmigroup4 + 
  ephrtdur3.bmigroup3 + ephrtdur3.bmigroup4

dataset.predictors <- all.vars(endo.formula)[-1]


# --- b. Group LASSO ---

data <- e2c2.nomiss %>% mutate(foldid=as.numeric(factor(site)))
d.train <- data[,c("casecon", "foldid", predictors)] %>% drop_na  
foldid <- d.train$foldid
nfolds <- length(unique(foldid))
y.train <- d.train %>% select(casecon) %>% mutate(casecon=ifelse(casecon==0, -1, 1)) %>% as.matrix()
x.train <- d.train %>% select(-casecon, -foldid) %>% as.matrix()
grid <- 10^seq(10, -5, length=200)
group.var <- str_replace(predictors, "999", "missing") %>%
  str_replace_all(., "([0-9])", "") %>%
  fct_inorder() %>%
  as.numeric()
penalty <- c(1-as.numeric(grepl("site|missing|999|age.group50|age.group55|age.group6|age.group7|age.group8", 
                                predictors)))*sqrt(as.numeric(ave(predictors, group.var, FUN=length)))
penalty.factor <- aggregate(penalty, by=list(group.var), mean) %>% pull(x)

# Cross-validation for different lambda values
if (1==0){
  cv.out <- cv.gglasso(x.train, y.train, group=group.var, nfolds=nfolds, foldid=foldid, 
                       lambda=grid, loss="logit", pf=penalty.factor)
  saveRDS(cv.out, file="Output/E2C2_Validation_3_LASSO_results.rds")
}
cv.out <- readRDS("Output/E2C2_Validation_3_LASSO_results.rds")

# Pulling out lambda with minimum values
# lasso.results <- coef(cv.out$gglasso.fit, s = cv.out$lambda.1se)
lasso.results <- coef(cv.out$gglasso.fit, s = cv.out$lambda.min)

# Creating list of RRs (With Interaction Terms)
lasso.results.filter <- lasso.results %>% 
  data.frame() %>%
  mutate(var=rownames(lasso.results)) %>%
  filter(var %in% dataset.predictors)
endo.log.RR <- lasso.results.filter$X1 
names(endo.log.RR) <- lasso.results.filter$var


# ---------------------------------------------------------------
# ------------ (10) Model development and validation ------------
# ---------------------------------------------------------------

# --- a. NHS I: Clinical model in whole dataset ---

# Developing model 
nhsi10.clinicalrisk <- computeAbsoluteRisk(
  model.formula=endo.formula,
  model.cov.info=endo.cov.info,
  model.log.RR=endo.log.RR,
  model.ref.dataset=nhanes.nhsi[,dataset.predictors],
  model.ref.dataset.weights=nhanes.nhsi$weights,
  apply.age.start=start.age,
  apply.age.interval.length=40,
  model.disease.incidence.rates=base.endo.nhsi,
  model.competing.incidence.rates=base.competing.nhsi,
  apply.cov.profile=nhsi.final[,all.vars(endo.formula)[-1]],
  n.imp=5, use.c.code=1, return.lp=T, return.refs.risk=T)

# Validating Model
nhsi.clinical.model <- list(
  model.formula=endo.formula,
  model.cov.info=endo.cov.info,
  model.log.RR=endo.log.RR,
  model.ref.dataset=nhanes.nhsi[,dataset.predictors],
  model.ref.dataset.weights=nhanes.nhsi$weights,
  apply.age.start=start.age,
  model.disease.incidence.rates=base.endo.nhsi,
  model.competing.incidence.rates=base.competing.nhsi,
  apply.cov.profile=nhsi.final[,all.vars(endo.formula)[-1]],
  n.imp=5, use.c.code=1, return.lp=T, return.refs.risk=T)

nhsi10.clinical <- ModelValidation(
  study.data=nhsi.final, 
  predicted.risk.interval=10, 
  iCARE.model.object=nhsi.clinical.model, 
  number.of.percentiles=10)


# --- b. NHS I: Clinical model among participants with genetic data ---

nhsi.final.cg <- nhsi.final %>% filter(gwas==1)

# Model Development
nhsi.genetic.model1 <- list(
  model.formula=endo.formula,
  model.cov.info=endo.cov.info,
  model.log.RR=endo.log.RR,
  model.ref.dataset=nhanes.nhsi[,dataset.predictors],
  model.ref.dataset.weights=nhanes.nhsi$weights,
  model.disease.incidence.rates=base.endo.nhsi,
  model.competing.incidence.rates=base.competing.nhsi,
  apply.cov.profile=nhsi.final.cg[,all.vars(endo.formula)[-1]],
  n.imp=5, use.c.code=1, return.lp=T, return.refs.risk=T)
  
# Model Validation
nhsi10.genetic1 <- ModelValidation(
  study.data=nhsi.final.cg, 
  predicted.risk.interval=10, 
  iCARE.model.object=nhsi.genetic.model1, 
  number.of.percentiles=10)

# --- c. NHS I: Clinical+genetic model among participants with genetic data ---

# Model development
nhsi.genetic.model2 <- list(
  model.formula=endo.formula,
  model.cov.info=endo.cov.info,
  model.snp.info=snps,
  model.log.RR=endo.log.RR,
  model.ref.dataset=nhanes.nhsi[,dataset.predictors],
  model.ref.dataset.weights=nhanes.nhsi$weights,
  apply.age.start=start.age,
  model.disease.incidence.rates=base.endo.nhsi,
  model.competing.incidence.rates=base.competing.nhsi,
  apply.cov.profile=nhsi.final.cg[,all.vars(endo.formula)[-1]],
  apply.snp.profile=nhsi.final.cg[,snps$snp.name],
  model.bin.fh.name=NA,
  n.imp=5, use.c.code=1, return.lp=T, return.refs.risk=T)

nhsi10.genetic2 <- ModelValidation(
  study.data=nhsi.final.cg, 
  predicted.risk.interval=10, 
  iCARE.model.object=nhsi.genetic.model2, 
  number.of.percentiles=10)


# --- d. NHS II: Clinical model in whole dataset ---

# Model development
nhsii10.clinicalrisk <- computeAbsoluteRisk(
  model.formula=endo.formula,
  model.cov.info=endo.cov.info,
  model.log.RR=endo.log.RR,
  model.ref.dataset=nhanes.nhsii[,dataset.predictors],
  model.ref.dataset.weights=nhanes.nhsii$weights,
  apply.age.start=45,
  apply.age.interval.length=40,
  model.disease.incidence.rates=base.endo.nhsii,
  model.competing.incidence.rates=base.competing.nhsii,
  apply.cov.profile=nhsii.final[,all.vars(endo.formula)[-1]],
  n.imp=5, use.c.code=1, return.lp=T, return.refs.risk=T)

# Model validation
nhsii.clinical.model <- list(
  model.formula=endo.formula,
  model.cov.info=endo.cov.info,
  model.log.RR=endo.log.RR,
  model.ref.dataset=nhanes.nhsii[,dataset.predictors],
  model.ref.dataset.weights=nhanes.nhsii$weights,
  model.disease.incidence.rates=base.endo.nhsii,
  model.competing.incidence.rates=base.competing.nhsii,
  apply.cov.profile=nhsii.final[,all.vars(endo.formula)[-1]],
  n.imp=5, use.c.code=1, return.lp=T, return.refs.risk=T)

nhsii10.clinical <- ModelValidation(
  study.data=nhsii.final, 
  predicted.risk.interval=10, 
  iCARE.model.object=nhsii.clinical.model, 
  number.of.percentiles=10)


# --- e. PLCO: Clinical model in whole dataset ---

# Model development
plco10.clinicalrisk <- computeAbsoluteRisk(
  model.formula=endo.formula,
  model.cov.info=endo.cov.info,
  model.log.RR=endo.log.RR,
  model.ref.dataset=nhanes.plco[,dataset.predictors],
  model.ref.dataset.weights=nhanes.plco$weights,
  apply.age.start=start.age,
  apply.age.interval.length=40,
  model.disease.incidence.rates=base.endo.plco,
  model.competing.incidence.rates=base.competing.plco,
  apply.cov.profile=plco.final[,all.vars(endo.formula)[-1]],
  n.imp=5, use.c.code=1, return.lp=T, return.refs.risk=T)

# Model validation
plco.clinical.model <- list(
  model.formula=endo.formula,
  model.cov.info=endo.cov.info,
  model.log.RR=endo.log.RR,
  model.ref.dataset=nhanes.plco[,dataset.predictors],
  model.ref.dataset.weights=nhanes.plco$weights,
  apply.age.start=start.age,
  apply.age.interval.length=10,
  model.disease.incidence.rates=base.endo.plco,
  model.competing.incidence.rates=base.competing.plco,
  apply.cov.profile=plco.final[,all.vars(endo.formula)[-1]],
  n.imp=5, use.c.code=1, return.lp=T, return.refs.risk=T)

plco10.clinical <- ModelValidation(
  study.data=plco.final, 
  predicted.risk.interval=10, 
  iCARE.model.object=plco.clinical.model, 
  number.of.percentiles=10)


# --- f. PLCO: Clinical model among participants with genetic data ---

plco.final.cg <- plco.final %>% filter(genetic==1)

# Model Development
plco.genetic.model1 <- list(
  model.formula=endo.formula,
  model.cov.info=endo.cov.info,
  model.log.RR=endo.log.RR,
  model.ref.dataset=nhanes.plco[,dataset.predictors],
  model.ref.dataset.weights=nhanes.plco$weights,
  model.disease.incidence.rates=base.endo.plco,
  model.competing.incidence.rates=base.competing.plco,
  apply.cov.profile=plco.final.cg[,all.vars(endo.formula)[-1]],
  n.imp=5, use.c.code=1, return.lp=T, return.refs.risk=T)

plco10.genetic1 <- ModelValidation(
  study.data=plco.final.cg, 
  predicted.risk.interval=10, 
  iCARE.model.object=plco.genetic.model1, 
  number.of.percentiles=10)


# --- g. PLCO: Clinical+genetic model among participants with genetic data ---

# Model development
plco.genetic.model2 <- list(
  model.formula=endo.formula,
  model.cov.info=endo.cov.info,
  model.snp.info=snps,
  model.log.RR=endo.log.RR,
  model.ref.dataset=nhanes.plco[,dataset.predictors],
  model.ref.dataset.weights=nhanes.plco$weights,
  apply.age.start=start.age,
  model.disease.incidence.rates=base.endo.plco,
  model.competing.incidence.rates=base.competing.plco,
  apply.cov.profile=plco.final.cg[,all.vars(endo.formula)[-1]],
  apply.snp.profile=plco.final.cg[,snps$snp.name],
  model.bin.fh.name=NA,
  n.imp=5, use.c.code=1, return.lp=T, return.refs.risk=T)

# Model validation
plco10.genetic2 <- ModelValidation(
  study.data=plco.final.cg, 
  predicted.risk.interval=10, 
  iCARE.model.object=plco.genetic.model2,
  number.of.percentiles=10)


# --- h. Current: Clinical model ---

# Model development
current.clinical.absmodel <- computeAbsoluteRisk(
  model.formula=endo.formula,
  model.cov.info=endo.cov.info,
  model.log.RR=endo.log.RR,
  model.ref.dataset=nhanes.current[,dataset.predictors],
  model.ref.dataset.weights=nhanes.current$weights,
  apply.age.start=start.age,
  apply.age.interval.length=40,
  model.disease.incidence.rates=base.endo.current,
  model.competing.incidence.rates=base.competing.current,
  apply.cov.profile=current.final[,all.vars(endo.formula)[-1]],
  n.imp=5, use.c.code=1, return.lp=T, return.refs.risk=T)

# Model validation
current.clinical.model <- list(
  model.formula=endo.formula,
  model.cov.info=endo.cov.info,
  model.log.RR=endo.log.RR,
  model.ref.dataset=nhanes.current[,dataset.predictors],
  model.ref.dataset.weights=nhanes.current$weights,
  apply.age.start=start.age,
  apply.age.interval.length=40,
  model.disease.incidence.rates=base.endo.current,
  model.competing.incidence.rates=base.competing.current,
  apply.cov.profile=current.final[,all.vars(endo.formula)[-1]],
  n.imp=5, use.c.code=1, return.lp=T, return.refs.risk=T)

current10.clinical <- ModelValidation(
  study.data=current.final, 
  predicted.risk.interval=10, 
  iCARE.model.object=current.clinical.model, 
  number.of.percentiles=10)


# ---------------------------------------------------------------
# ---------------------- (11) Figures 1-2 -----------------------
# ---------------------------------------------------------------

# --- a. Pulling out estimates from model validation output ---

t <- lapply(c("nhsi10.clinical", "nhsii10.clinical", 
  "plco10.clinical", "plco10.genetic1"), function(x){
    output <- get(x)
    return(
      cbind(output[["Category_Results"]],
      HL=rep(output[["Hosmer_Lemeshow_Results"]]$statistic, 10),
      HL.p=rep(output[["HL_pvalue"]], 10),
      GOF=rep(output[["RR_test_result"]]$statistic, 10),
      GOF.p=rep(output[["RR_test_pvalue"]], 10),                         
      Study=rep(x, 10)))}) %>%
  bind_rows() %>%
  data.frame() %>%
  # Converting from charatcer to numeric
  mutate(Observed_Absolute_Risk=as.numeric(as.character(Observed_Absolute_Risk))*100) %>%
  mutate(Predicted_Absolute_Risk=as.numeric(as.character(Predicted_Absolute_Risk))*100) %>%
  mutate(CI_Absolute_Risk_Lower=as.numeric(as.character(CI_Absolute_Risk_Lower))*100) %>%
  mutate(CI_Absolute_Risk_Upper=as.numeric(as.character(CI_Absolute_Risk_Upper))*100) %>%
  mutate(Observed_Relative_Risk=as.numeric(as.character(Observed_Relative_Risk))) %>%
  mutate(Predicted_Relative_Risk=as.numeric(as.character(Predicted_Relative_Risk))) %>%
  mutate(CI_Relative_Risk_Lower=as.numeric(as.character(CI_Relative_Risk_Lower))) %>%
  mutate(CI_Relative_Risk_Upper=as.numeric(as.character(CI_Relative_Risk_Upper))) %>%
  mutate(Label_Absolute=ifelse(HL.p<0.001, sprintf("HL-statistic = %.1f\np < 0.001", HL), 
    sprintf("HL-statistic = %.1f\np = %.3f", HL, HL.p))) %>%
  mutate(Label_RR=ifelse(GOF.p<0.001, sprintf("GOF \U03C7\U00B2-statistic = %.1f\np < 0.001", GOF), 
    sprintf("GOF X\U00B2-statistic = %.1f\np = %.3f", GOF, GOF.p))) %>%
  group_by(Study) %>%
  # Creating labels for plots
  mutate(Label_Absolute=ifelse(row_number()==1, Label_Absolute, " ")) %>%
  mutate(Label_RR=ifelse(row_number()==1, Label_RR, " ")) %>%
  ungroup()


# --- b. Figure 1 ---

fig3a1 <- t %>% 
  filter(Study=="nhsi10.clinical") %>%
  ggplot(aes(x=Predicted_Absolute_Risk, y=Observed_Absolute_Risk)) +
  geom_abline(intercept=0, slope=1, linetype="dashed", color="grey") +
  geom_errorbar(aes(ymin=CI_Absolute_Risk_Lower, ymax=CI_Absolute_Risk_Upper), width=0.025) + 
  geom_point() + 
  geom_text(aes(label=Label_Absolute), y=0, x=3.6, hjust=1, vjust=0) +
  scale_x_continuous(" \n Expected Absolute 10-Year Risk (%)", limits=c(0, 3.5)) +
  scale_y_continuous(" \nObserved Absolute 10-Year Risk (%) \n ", limits=c(0, 3.5)) + 
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill=NA), 
    strip.background=element_blank())

fig3a2 <- t %>% 
  filter(Study=="nhsi10.clinical") %>%
  ggplot(aes(x=Predicted_Relative_Risk, y=Observed_Relative_Risk)) +
  geom_abline(intercept=0, slope=1, linetype="dashed", color="grey") +
  geom_errorbar(aes(ymin=CI_Relative_Risk_Lower, ymax=CI_Relative_Risk_Upper), width=0.025) + 
  geom_point() + 
  scale_x_continuous(" \nExpected Relative 10-Year Risk", trans="log2", 
    limits=c(0.09, 4), breaks=c(2^(-4:4))) +
  scale_y_continuous(" \nObserved Relative 10-Year Risk \n ", trans="log2", 
    limits=c(0.09, 4), breaks=c(2^(-4:4))) + 
  geom_text(aes(label=Label_RR), y=log2(0.09), x=log2(4.1), hjust=1, vjust=0) +
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill=NA), 
    strip.background=element_blank())

fig3b1 <- t %>% 
  filter(Study=="nhsii10.clinical") %>%
  ggplot(aes(x=Predicted_Absolute_Risk, y=Observed_Absolute_Risk)) +
  geom_abline(intercept=0, slope=1, linetype="dashed", color="grey") +
  geom_errorbar(aes(ymin=CI_Absolute_Risk_Lower, ymax=CI_Absolute_Risk_Upper), width=0.025) + 
  geom_point() + 
  geom_text(aes(label=Label_Absolute), y=0, x=3.5, hjust=1, vjust=0) +
  scale_x_continuous(" \nExpected Absolute 10-Year Risk (%)", limits=c(0, 3.5)) +
  scale_y_continuous(" \nObserved Absolute 10-Year Risk (%) \n ", limits=c(0, 3.5)) + 
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill=NA), 
    strip.background=element_blank())

fig3b2 <- t %>% 
  filter(Study=="nhsii10.clinical") %>%
  ggplot(aes(x=Predicted_Relative_Risk, y=Observed_Relative_Risk)) +
  geom_abline(intercept=0, slope=1, linetype="dashed", color="grey") +
  geom_errorbar(aes(ymin=CI_Relative_Risk_Lower, ymax=CI_Relative_Risk_Upper), width=0.025) + 
  geom_point() + 
  scale_x_continuous(" \nExpected Relative 10-Year Risk", trans="log2", 
    limits=c(0.09, 4), breaks=c(2^(-4:4))) +
  scale_y_continuous(" \nObserved Relative 10-Year Risk \n ", trans="log2", 
    limits=c(0.09, 4), breaks=c(2^(-4:4))) + 
  geom_text(aes(label=Label_RR), y=log2(0.09), x=log2(4.1), hjust=1, vjust=0) +
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill=NA), 
    strip.background=element_blank())

fig3 <- plot_grid(fig3a1, fig3a2, fig3b1, fig3b2, 
                   nrow=2,
                   labels=c("A", " ", "B", " "),
                   label_size=11,
                   rel_widths=c(5.25/11,5.75/11))

ggsave("./Output/Figure_1.pdf", plot=fig3, device="pdf", width=8, height=7.5, units="in")  


# --- d. Figure 2 ---

fig4a1 <- t %>% 
  filter(Study=="plco10.clinical") %>%
  ggplot(aes(x=Predicted_Absolute_Risk, y=Observed_Absolute_Risk)) +
  geom_abline(intercept=0, slope=1, linetype="dashed", color="grey") +
  geom_errorbar(aes(ymin=CI_Absolute_Risk_Lower, ymax=CI_Absolute_Risk_Upper), width=0.025) + 
  geom_point() + 
  geom_text(aes(label=Label_Absolute), y=0, x=4, hjust=1, vjust=0) +
  scale_x_continuous(" \nExpected Absolute 10-Year Risk (%)", limits=c(0, 4)) +
  scale_y_continuous(" \nObserved Absolute 10-Year Risk (%) \n ", limits=c(0, 4)) + 
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill=NA), strip.background=element_blank())

fig4a2 <- t %>% 
  filter(Study=="plco10.clinical") %>%
  ggplot(aes(x=Predicted_Relative_Risk, y=Observed_Relative_Risk)) +
  geom_abline(intercept=0, slope=1, linetype="dashed", color="grey") +
  geom_errorbar(aes(ymin=CI_Relative_Risk_Lower, ymax=CI_Relative_Risk_Upper), width=0.025) + 
  geom_point() + 
  scale_x_continuous(" \nExpected Relative 10-Year Risk", trans="log2", limits=c(0.09, 4), breaks=c(2^(-4:4))) +
  scale_y_continuous(" \nObserved Relative 10-Year Risk \n ", trans="log2", limits=c(0.09, 4), breaks=c(2^(-4:4))) + 
  geom_text(aes(label=Label_RR), y=log2(0.09), x=log2(4.1), hjust=1, vjust=0) +
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill=NA), strip.background=element_blank())

fig4b1 <- t %>% 
  filter(Study=="plco10.genetic1") %>%
  ggplot(aes(x=Predicted_Absolute_Risk, y=Observed_Absolute_Risk)) +
  geom_abline(intercept=0, slope=1, linetype="dashed", color="grey") +
  geom_errorbar(aes(ymin=CI_Absolute_Risk_Lower, ymax=CI_Absolute_Risk_Upper), width=0.025) + 
  geom_point() + 
  geom_text(aes(label=Label_Absolute), y=0, x=4, hjust=1, vjust=0) +
  scale_x_continuous(" \nExpected Absolute 10-Year Risk (%)", limits=c(0, 4)) +
  scale_y_continuous(" \nObserved Absolute 10-Year Risk (%) \n ", limits=c(0, 4)) + 
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill=NA), strip.background=element_blank())

fig4b2 <- t %>% 
  filter(Study=="plco10.genetic1") %>%
  ggplot(aes(x=Predicted_Relative_Risk, y=Observed_Relative_Risk)) +
  geom_abline(intercept=0, slope=1, linetype="dashed", color="grey") +
  geom_errorbar(aes(ymin=CI_Relative_Risk_Lower, ymax=CI_Relative_Risk_Upper), width=0.025) + 
  geom_point() + 
  scale_x_continuous(" \nExpected Relative 10-Year Risk", trans="log2", limits=c(0.09, 4), breaks=c(2^(-4:4))) +
  scale_y_continuous(" \nObserved Relative 10-Year Risk \n ", trans="log2", limits=c(0.09, 4), breaks=c(2^(-4:4))) + 
  geom_text(aes(label=Label_RR), y=log2(0.09), x=log2(4.1), hjust=1, vjust=0) +
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill=NA), strip.background=element_blank())

fig4 <- plot_grid(fig4a1, fig4a2, fig4b1, fig4b2, 
                   nrow=2,
                   labels=c("A", " ", "B", " "),
                   label_size=11,
                   rel_widths=c(5.25/11,5.75/11))


ggsave("./Output/Figure_2.pdf", plot=fig4, device="pdf", width=8, height=7.5, units="in")  


# ---------------------------------------------------------------
# ------------------------ (12) Figure 3 ------------------------
# ---------------------------------------------------------------

# Calculating baseline 10-year risks 
t1 <- base.endo.current %>%
  data.frame() %>%
  mutate(Survival=1-rate.corrected) %>%
  mutate(Risk10=1-Survival^10)

# Calculating predicted risk for each individual in NHANES data
t2 <- data.frame(weight=current.final$weight, age=current.final$age, 
                 obs.risk=current10.clinical$Subject_Specific_Predicted_Absolute_Risk) %>%
  left_join(t1, c("age"="age.group")) %>%
  mutate(rr=obs.risk/Risk10)

# Calculating quantiles of risk
quantiles <- c(wtd.quantile(t2$rr, weights=t2$weight, 
                            probs=c(0.05, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99)), 1)

# Calculating baseline cumulative risk
t3 <- base.endo.current %>%
  data.frame() %>%
  mutate(survival=1-rate.corrected) %>%
  mutate(cum.survival=cumprod(survival)) %>%
  mutate(cum.risk=1-cum.survival) 

out <- t3$survival
for (i in 1:9){
  out <- out[1:(length(out)-1)] * t3$survival[(i+1):length(t3$survival)]
}

t3.final <- t3 %>% 
  mutate(survival10=c(out, rep(NA, 9))) %>%
  mutate(risk10=1-survival10)

# Legend
legend <- t3.final[["cum.risk"]] %*% t(quantiles) %>%
  data.frame() %>%
  mutate(Age=40:85) %>%
  pivot_longer(cols=c(starts_with("X"), starts_with("V")), 
               names_to="Percentile", values_to="Risk") %>%
  mutate(Percentile=as.character(as.numeric(gsub("\\D", "", Percentile)))) %>%
  mutate(Percentile=ifelse(Percentile=="8", "Average", Percentile)) %>%
  mutate(Percentile=as.factor(Percentile)) %>%
  mutate(Percentile=fct_relevel(Percentile, "5")) %>%
  ggplot(aes(x=Age, y=Risk, linetype=Percentile, color=Percentile)) +
  geom_line(size=0.75) +
  scale_x_continuous(" \nAge") +
  scale_y_continuous("Cumulative Absolute Risk\n ", limits=c(0, 0.255),
                     breaks=seq(0, 0.25, 0.05), labels=seq(0, 0.25, 0.05)) +
  scale_color_manual(values=c(brewer.pal(n = 7, name = "Blues"), "black")) +
  scale_linetype_manual(values=c(rep("solid", 7), "dotted")) + 
  theme_bw() +
  theme(text=element_text(size=9))
legend <- get_legend(legend)

# Plotting cumulative risks across various quantiles of risk
figs2a <- t3.final[["cum.risk"]] %*% t(quantiles) %>%
  data.frame() %>%
  mutate(Age=40:85) %>%
  pivot_longer(cols=c(starts_with("X"), starts_with("V")), 
               names_to="Percentile", values_to="Risk") %>%
  mutate(Percentile=as.character(as.numeric(gsub("\\D", "", Percentile)))) %>%
  mutate(Percentile=ifelse(Percentile=="8", "Average", Percentile)) %>%
  mutate(Percentile=as.factor(Percentile)) %>%
  mutate(Percentile=fct_relevel(Percentile, "5")) %>%
  ggplot(aes(x=Age, y=Risk, linetype=Percentile, color=Percentile)) +
  geom_line(size=0.75) +
  scale_x_continuous(" \nAge") +
  scale_y_continuous("Cumulative Absolute Risk (%)\n ", limits=c(0, 0.255),
                     breaks=seq(0, 0.25, 0.05), labels=seq(0, 0.25, 0.05)) +
  scale_color_manual(values=c(brewer.pal(n = 7, name = "Blues"), "black")) +
  scale_linetype_manual(values=c(rep("solid", 7), "dotted")) + 
  theme_bw() +
  theme(legend.position="none", text=element_text(size=9))

# Plotting 10-year absolute risks across various quantiles of risk
figs2b <- t3.final[["risk10"]] %*% t(quantiles) %>%
  data.frame() %>%
  mutate(Age=40:85) %>%
  filter(Age<=75) %>%
  pivot_longer(cols=c(starts_with("X"), starts_with("V")), 
               names_to="Percentile", values_to="Risk") %>%
  mutate(Percentile=as.character(as.numeric(gsub("\\D", "", Percentile)))) %>%
  mutate(Percentile=ifelse(Percentile=="8", "Average", Percentile)) %>%
  mutate(Percentile=as.factor(Percentile)) %>%
  mutate(Percentile=fct_relevel(Percentile, "5")) %>%
  ggplot(aes(x=Age, y=Risk, linetype=Percentile, color=Percentile)) +
  geom_line(size=0.75) +
  scale_x_continuous(" \nAge") +
  scale_y_continuous("Absolute 10-Year Risk\n ", limits=c(0, 0.1),
                     breaks=seq(0, 0.1, 0.02), labels=seq(0, 0.1, 0.02)) +
  scale_color_manual(values=c(brewer.pal(n = 7, name = "Blues"), "black")) +
  scale_linetype_manual(values=c(rep("solid", 7), "dotted")) + 
  theme_bw() +
  theme(legend.position="none", text=element_text(size=9))

# Bar graph of distribution of cumulative risk
t4 <- t2 %>%
  mutate(cum.risk=rr*t3[t3$age.group==85, "cum.risk"]) %>%
  mutate(risk.group=ifelse(cum.risk<=0.025, "0 to <2.5%",
    ifelse(cum.risk<=0.05, "2.5% to <5%",
    ifelse(cum.risk<=0.075, "5% to <7.5%",
    ifelse(cum.risk<=0.10, "7.5% to <10%",
    ifelse(cum.risk<=0.125, "10% to <12.5%",
    ifelse(cum.risk<=0.15, "12.5% to <15%",
    ifelse(cum.risk<=0.175, "15% to <17.5%",
    ifelse(cum.risk<=0.2, "17.5% to <20%", "≥20%"))))))))) %>%
  mutate(risk.group=factor(risk.group, levels=c(
    "0 to <2.5%",
    "2.5% to <5%",
    "5% to <7.5%",
    "7.5% to <10%",
    "10% to <12.5%",
    "12.5% to <15%",
    "15% to <17.5%",
    "17.5% to <20%", 
    "≥20%"))) %>%   
  count(risk.group, wt=weight) %>%
  mutate(percent=n/sum(n)*100) %>%
  mutate(label=sprintf("%.1f", percent))

figs2c <- ggplot(data=t4, aes(x=risk.group, y=percent)) +
  geom_bar(stat="identity", fill=brewer.pal(n = 7, name = "Blues")[7]) +
  geom_text(aes(label=label), vjust=-0.3, size=(5/14)*7) +
  scale_x_discrete(" \nCumulative Absolute Risk",
                   labels=c("0 to <2.5%",
                            "2.5% to <5%",
                            "5% to <7.5%",
                            "7.5% to <10%",
                            "10% to <12.5%",
                            "12.5% to <15%",
                            "15% to <17.5%",
                            "17.5% to <20%", 
                            expression(phantom(0)>="20%"))) +
  scale_y_continuous("Proportion (%)\n \n ", expand=c(0,0.1), limit=c(0, 40)) +
  theme_bw() +
  theme(legend.position="none", text=element_text(size=9), axis.text.x=element_text(angle=45, hjust=1),
        panel.grid.major.x=element_blank())

figs2 <- plot_grid(figs2a, legend, figs2c, NULL,
                   nrow=2,
                   labels=c("A", " ", "B", " "),
                   label_size=10, 
                   rel_widths=c(4/9, 1/9, 4/9, 1/9),
                   rel_heights=c(4.25/9, 4.75/9))

ggsave("./Output/Figure_3.pdf", plot=figs2, device="pdf", width=4, height=7, units="in", encoding=) 


# Calculating what percentile of women have RR>=3
rr.greater3 <- weighted.mean(t2$rr>3, w=t2$weight)

# Calculating how many women have predicted cumulative risk >=20%
cumrisk.greater20 <- weighted.mean(t2$rr>(0.20/t3[which(t3$age.group==85), "cum.risk"]), w=t2$weight)

# Calculating average cumulative abs risk
cumrisk.avg <- weighted.mean(t3[which(t3$age.group==85),]$cum.risk)

write.csv(cbind(rr.greater3, cumrisk.greater20, cumrisk.avg), "./Output/Manuscript_Text_Current_Risks.csv", row.names=F)

# ---------------------------------------------------------------
# ---------------------- (13) Tables 1-3 ------------------------
# ---------------------------------------------------------------

# --- a. Table 1 ---

t1 <- do.call(bind_rows, lapply(c(
  "nhsi", "nhsi10.clinical", "nhsi10.genetic1", "nhsi10.genetic2", 
  "nhsii", "nhsii10.clinical", 
  "plco", "plco10.clinical", "plco10.genetic1", "plco10.genetic2"), function(x){
    if(exists(x)){
      output <- get(x)
      results <- c(
        analysis=x,
        n=length(output[["Subject_Specific_Observed_Outcome"]]),
        events=sum(output[["Subject_Specific_Observed_Outcome"]]),
        auc=output[["AUC"]],
        auc.ll=output[["CI_AUC"]][1],
        auc.ul=output[["CI_AUC"]][2])
        return(results)
      } else{
        results <- c(analysis=x)
        return(results)
      }})) %>%
  data.frame() %>%
  mutate_at(c("n", "auc", "auc.ll", "auc.ul"), as.numeric) %>%
  mutate(auc=ifelse(!is.na(auc), sprintf("%.3f (%.3f, %.3f)", auc, auc.ll, auc.ul), "")) %>%
  mutate(n=trimws(sprintf("%s", format(n, big.mark=",")))) %>%
  mutate(n=ifelse(n=="NA", "", n)) %>%
  select(analysis, n, events, auc)

write.csv(t1, "./Output/Table_1.csv", row.names=F, na="")


# --- b. Table 2 ---

# Extracing validation results for NHS I
t.nhsi <- bind_cols(outcome=nhsi10.clinical$Subject_Specific_Observed_Outcome,
                    p_risk=nhsi10.clinical$Subject_Specific_Predicted_Absolute_Risk,
                    p_rs=nhsi10.clinical$Subject_Specific_Risk_Score) %>%
  mutate(decile=cut(p_rs, breaks=c(-Inf, quantile(nhsi10.clinical$Subject_Specific_Risk_Score, 
    probs=seq(0.1, 0.9, 0.1)), Inf), labels=c(1:10))) %>%
  mutate(n=1) %>%
  group_by(decile) %>%
  summarize_at(vars(n, outcome), list(n=sum)) %>%
  bind_cols(nhsi10.clinical$Category_Results) %>%
  select(-Categories) %>%
  mutate_if(is.character, as.numeric) %>%
  mutate(var_obsrisk=(Observed_Absolute_Risk)*(1-Observed_Absolute_Risk)/n_n) %>%
  mutate(var_log_exp_by_obs=var_obsrisk/(Observed_Absolute_Risk)^2) %>%
  mutate(EO=Predicted_Absolute_Risk/Observed_Absolute_Risk) %>%
  mutate(EO_LCI=exp(log(EO) - 1.96 * sqrt(var_log_exp_by_obs))) %>%
  mutate(EO_UCI=exp(log(EO) + 1.96 * sqrt(var_log_exp_by_obs))) %>%
  mutate(Study="NHS") %>%
  mutate(decile=as.numeric(as.character(decile))) %>%
  add_row(Study="NHS", decile=11, 
          n_n=length(nhsi10.clinical$Subject_Specific_Observed_Outcome), 
          outcome_n=sum(nhsi10.clinical$Subject_Specific_Observed_Outcome), 
          EO=nhsi10.clinical$Overall_Expected_to_Observed_Ratio, 
          EO_LCI=nhsi10.clinical$CI_Overall_Expected_to_Observed_Ratio[1], 
          EO_UCI=nhsi10.clinical$CI_Overall_Expected_to_Observed_Ratio[2])

# Extracting validation results for NHS II
t.nhsii <- bind_cols(outcome=nhsii10.clinical$Subject_Specific_Observed_Outcome,
                     p_risk=nhsii10.clinical$Subject_Specific_Predicted_Absolute_Risk,
                     p_rs=nhsii10.clinical$Subject_Specific_Risk_Score) %>%
  mutate(decile=cut(p_rs, breaks=c(-Inf, quantile(nhsii10.clinical$Subject_Specific_Risk_Score, 
    probs=seq(0.1, 0.9, 0.1)), Inf), labels=c(1:10))) %>%
  mutate(n=1) %>%
  group_by(decile) %>%
  summarize_at(vars(n, outcome), list(n=sum)) %>%
  bind_cols(nhsii10.clinical$Category_Results) %>%
  select(-Categories) %>%
  mutate_if(is.character, as.numeric) %>%
  mutate(var_obsrisk=(Observed_Absolute_Risk)*(1-Observed_Absolute_Risk)/n_n) %>%
  mutate(var_log_exp_by_obs=var_obsrisk/(Observed_Absolute_Risk)^2) %>%
  mutate(EO=Predicted_Absolute_Risk/Observed_Absolute_Risk) %>%
  mutate(EO_LCI=exp(log(EO) - 1.96 * sqrt(var_log_exp_by_obs))) %>%
  mutate(EO_UCI=exp(log(EO) + 1.96 * sqrt(var_log_exp_by_obs))) %>%
  mutate(Study="NHS II") %>%
  mutate(decile=as.numeric(as.character(decile))) %>%
  add_row(Study="NHS II", decile=11, 
          n_n=length(nhsii10.clinical$Subject_Specific_Observed_Outcome), 
          outcome_n=sum(nhsii10.clinical$Subject_Specific_Observed_Outcome), 
          EO=nhsii10.clinical$Overall_Expected_to_Observed_Ratio, 
          EO_LCI=nhsii10.clinical$CI_Overall_Expected_to_Observed_Ratio[1], 
          EO_UCI=nhsii10.clinical$CI_Overall_Expected_to_Observed_Ratio[2])  

# Combining results to create table
t2 <- bind_rows(t.nhsi, t.nhsii) %>%
  mutate(EO=ifelse(EO==Inf, "N/A", sprintf("%.2f (%.2f, %.2f)", EO, EO_LCI, EO_UCI))) %>%
  mutate(Observed.AbRisk=ifelse(decile!=11, sprintf("%.2f (%.2f, %.2f)", 
    Observed_Absolute_Risk*100, CI_Absolute_Risk_Lower*100, CI_Absolute_Risk_Upper*100), "")) %>%
  mutate(Expected.AbRisk=ifelse(decile!=11, sprintf("%.2f", Predicted_Absolute_Risk*100), "")) %>%
  mutate(Observed.RelRisk=ifelse(decile!=11, sprintf("%.2f (%.2f, %.2f)", 
    Observed_Relative_Risk, CI_Relative_Risk_Lower, CI_Relative_Risk_Upper), "")) %>%
  mutate(Expected.RelRisk=ifelse(decile!=11, sprintf("%.2f", Predicted_Relative_Risk), "")) %>%
  rename(N=n_n, Outcome=outcome_n, Decile=decile) %>%
  mutate(N=trimws(sprintf("%s", format(N, big.mark=",")))) %>%
  add_row(Study="NHS", Decile=0) %>%
  add_row(Study="NHS II", Decile=0) %>%
  select(Study, Decile, N, Outcome, Observed.RelRisk, Expected.RelRisk, 
    Observed.AbRisk, Expected.AbRisk, EO) %>%
  arrange(Study, Decile)

write.csv(t2, "./Output/Table_2.csv", row.names=F, na="")


# --- c. Table 3 ---

# Extracting validation results for clinical-only model in PLCO
t.plco.c <- bind_cols(outcome=plco10.clinical$Subject_Specific_Observed_Outcome,
                      p_risk=plco10.clinical$Subject_Specific_Predicted_Absolute_Risk,
                      p_rs=plco10.clinical$Subject_Specific_Risk_Score) %>%
  mutate(decile=cut(p_rs, breaks=c(-Inf, quantile(plco10.clinical$Subject_Specific_Risk_Score, 
    probs=seq(0.1, 0.9, 0.1)), Inf), labels=c(1:10))) %>%
  mutate(n=1) %>%
  group_by(decile) %>%
  summarize_at(vars(n, outcome), list(n=sum)) %>%
  bind_cols(plco10.clinical$Category_Results) %>%
  select(-Categories) %>%
  mutate_if(is.character, as.numeric) %>%
  mutate(var_obsrisk=(Observed_Absolute_Risk)*(1-Observed_Absolute_Risk)/n_n) %>%
  mutate(var_log_exp_by_obs=var_obsrisk/(Observed_Absolute_Risk)^2) %>%
  mutate(EO=Predicted_Absolute_Risk/Observed_Absolute_Risk) %>%
  mutate(EO_LCI=exp(log(EO) - 1.96 * sqrt(var_log_exp_by_obs))) %>%
  mutate(EO_UCI=exp(log(EO) + 1.96 * sqrt(var_log_exp_by_obs))) %>%
  mutate(Study="PLCO-C") %>%
  mutate(decile=as.numeric(as.character(decile))) %>%
  add_row(Study="PLCO-C", decile=11, 
          n_n=length(plco10.clinical$Subject_Specific_Observed_Outcome), 
          outcome_n=sum(plco10.clinical$Subject_Specific_Observed_Outcome), 
          EO=plco10.clinical$Overall_Expected_to_Observed_Ratio, 
          EO_LCI=plco10.clinical$CI_Overall_Expected_to_Observed_Ratio[1], 
          EO_UCI=plco10.clinical$CI_Overall_Expected_to_Observed_Ratio[2])    

# Extracting validation results for C+G model in PLCO
t.plco.g <- bind_cols(outcome=plco10.genetic1$Subject_Specific_Observed_Outcome,
                      p_risk=plco10.genetic1$Subject_Specific_Predicted_Absolute_Risk,
                      p_rs=plco10.genetic1$Subject_Specific_Risk_Score) %>%
  mutate(decile=cut(p_rs, breaks=c(-Inf, quantile(plco10.genetic1$Subject_Specific_Risk_Score, 
    probs=seq(0.1, 0.9, 0.1)), Inf), labels=c(1:10))) %>%
  mutate(n=1) %>%
  group_by(decile) %>%
  summarize_at(vars(n, outcome), list(n=sum)) %>%
  bind_cols(plco10.genetic1$Category_Results) %>%
  select(-Categories) %>%
  mutate_if(is.character, as.numeric) %>%
  mutate(var_obsrisk=(Observed_Absolute_Risk)*(1-Observed_Absolute_Risk)/n_n) %>%
  mutate(var_log_exp_by_obs=var_obsrisk/(Observed_Absolute_Risk)^2) %>%
  mutate(EO=Predicted_Absolute_Risk/Observed_Absolute_Risk) %>%
  mutate(EO_LCI=exp(log(EO) - 1.96 * sqrt(var_log_exp_by_obs))) %>%
  mutate(EO_UCI=exp(log(EO) + 1.96 * sqrt(var_log_exp_by_obs))) %>%
  mutate(Study="PLCO-CG") %>%
  mutate(decile=as.numeric(as.character(decile))) %>%
  add_row(Study="PLCO-CG", decile=11, 
          n_n=length(plco10.genetic1$Subject_Specific_Observed_Outcome), 
          outcome_n=sum(plco10.genetic1$Subject_Specific_Observed_Outcome), 
          EO=plco10.genetic1$Overall_Expected_to_Observed_Ratio, 
          EO_LCI=plco10.genetic1$CI_Overall_Expected_to_Observed_Ratio[1], 
          EO_UCI=plco10.genetic1$CI_Overall_Expected_to_Observed_Ratio[2]) 

# Combining results to create table
t3 <- bind_rows(t.plco.c, t.plco.g) %>%
  mutate(EO=ifelse(EO==Inf, "N/A", sprintf("%.2f (%.2f, %.2f)", EO, EO_LCI, EO_UCI))) %>%
  mutate(Observed.AbRisk=ifelse(decile!=11, sprintf("%.2f (%.2f, %.2f)", 
    Observed_Absolute_Risk*100, CI_Absolute_Risk_Lower*100, CI_Absolute_Risk_Upper*100), "")) %>%
  mutate(Expected.AbRisk=ifelse(decile!=11, sprintf("%.2f", Predicted_Absolute_Risk*100), "")) %>%
  mutate(Observed.RelRisk=ifelse(decile!=11, sprintf("%.2f (%.2f, %.2f)", 
    Observed_Relative_Risk, CI_Relative_Risk_Lower, CI_Relative_Risk_Upper), "")) %>%
  mutate(Expected.RelRisk=ifelse(decile!=11, sprintf("%.2f", Predicted_Relative_Risk), "")) %>%
  rename(N=n_n, Outcome=outcome_n, Decile=decile) %>%
  mutate(N=trimws(sprintf("%s", format(N, big.mark=",")))) %>%
  add_row(Study="PLCO-C", Decile=0) %>%
  add_row(Study="PLCO-CG", Decile=0) %>%
  select(Study, Decile, N, Outcome, Observed.RelRisk, Expected.RelRisk, 
    Observed.AbRisk, Expected.AbRisk, EO) %>%
  arrange(Study, Decile)

write.csv(t3, "./Output/Table_3.csv", row.names=F, na="")


# ---------------------------------------------------------------
# ----------------- (14) Supplementary Figures ------------------
# ---------------------------------------------------------------

# --- a. Figure S1 ---

# Determining exclusions post-baseline
nhsi.inclusion2 <- nhsi.long %>%
  mutate(baseline=ifelse(!is.na(year) & !is.na(start.year) & year==start.year, 1, 0)) %>%
  group_by(id) %>%
  mutate(included=max(baseline)) %>%
  mutate(excluded.year=max(exclude)) %>%
  distinct(id, .keep_all=T) %>%
  mutate(ex.cancer2=ifelse(!is.na(cancer.year) & !is.na(excluded.year) & 
    cancer.year==excluded.year & included==0, 1, 0)) %>%
  mutate(ex.hyst2=ifelse(!is.na(hyst.year) & !is.na(yobf) & 
    !is.na(menopause.age) & !is.na(excluded.year) & included==0 & 
    (hyst.year<=(yobf+1900+menopause.age+2)|(hyst.year-(yobf+1900))<=start.age|
    hyst.year==excluded.year), 1, 0)) %>%
  mutate(ex.hyst2=ifelse(ex.cancer2==1, 0, ex.hyst2)) %>%
  mutate(ex.death2=ifelse(!is.na(death.year) & !is.na(excluded.year) & 
    death.year==excluded.year & included==0, 1, 0)) %>%
  mutate(ex.death2=ifelse(ex.cancer2==1|ex.hyst2==1, 0, ex.death2)) %>%
  mutate(ex.questionnaire2=ifelse(!is.na(yobf) & !is.na(menopause.age) & 
    !is.na(excluded.year) & included==0 &  (1900+yobf+menopause.age+2)>=excluded.year, 1, 0)) %>%
  mutate(ex.questionnaire2=ifelse(ex.cancer2==1|ex.hyst2==1|ex.death2==1, 0, 
    ex.questionnaire2)) %>%
  mutate(ex.menopausal2=ifelse(is.na(menopause.age), 1, 0)) %>%
  mutate(ex.menopausal2=ifelse(ex.cancer2==1|ex.hyst2==1|ex.death2==1|
    ex.questionnaire2==1, 0, ex.menopausal2)) %>%
  mutate(ex.questionnaire2=ifelse(ex.cancer2==0 & ex.hyst2==0 & ex.death2==0 & 
    ex.questionnaire2==0 & included==0, 1, ex.questionnaire2)) %>%
  mutate(ex2=ifelse(ex.cancer2==1|ex.hyst2==1|ex.menopausal2==1|ex.death2==1|
    ex.questionnaire2==1, 1, 0)) %>%
  mutate(in2=1-ex2) %>%
  ungroup() %>%
  summarize_at(vars(ex2, ex.cancer2, ex.hyst2, ex.death2, ex.questionnaire2, 
    ex.menopausal2, in2), sum)

nhsi.inclusion3 <- data.frame(ex.cancer2=nhsi.inclusion2$in2-nrow(nhsi.final)) %>%
  mutate(ex2=ex.cancer2) %>%
  mutate(in2=-ex2)

# Combining pre and post exclusions
t1 <- nhsi.inclusion
t2 <- nhsi.inclusion2 %>% bind_rows(nhsi.inclusion3) %>% summarize_all(sum, na.rm=T)
sf1b <- bind_cols(t1, t2)

write.csv(sf1b, "./Output/Figure_S1b.csv", row.names=F, na="")


# --- b. Figure S2 ---

# Determining exclusions at baseline
nhsii.inclusion2 <- nhsii.long %>%
  mutate(baseline=ifelse(!is.na(year) & !is.na(start.year) & year==start.year, 1, 0)) %>%
  group_by(id) %>%
  mutate(included=max(baseline)) %>%
  mutate(excluded.year=max(exclude)) %>%
  distinct(id, .keep_all=T) %>%
  mutate(ex.cancer2=ifelse(!is.na(cancer.year) & !is.na(excluded.year) & cancer.year==excluded.year & included==0, 1, 0)) %>%
  mutate(ex.hyst2=ifelse(!is.na(hyst.year) & !is.na(yob) & !is.na(menopause.age) & !is.na(excluded.year) & included==0 & (hyst.year<=(floor(yob)+menopause.age+2)|(hyst.year-floor(yob))<=start.age|hyst.year==excluded.year), 1, 0)) %>%
  # mutate(ex.hyst2=ifelse(!is.na(hyst.year) & !is.na(exclude) & hyst.year==exclude & included==0, 1, 0)) %>%
  mutate(ex.hyst2=ifelse(ex.cancer2==1, 0, ex.hyst2)) %>%
  mutate(ex.death2=ifelse(!is.na(death.year) & !is.na(excluded.year) & death.year==excluded.year & included==0, 1, 0)) %>%
  mutate(ex.death2=ifelse(ex.cancer2==1|ex.hyst2==1, 0, ex.death2)) %>%
  mutate(ex.questionnaire2=ifelse(!is.na(yob) & !is.na(menopause.age) & !is.na(excluded.year) & included==0 &  (floor(yob)+menopause.age+2)>=excluded.year, 1, 0)) %>%
  mutate(ex.questionnaire2=ifelse(ex.cancer2==1|ex.hyst2==1|ex.death2==1, 0, ex.questionnaire2)) %>%
  mutate(ex.menopausal2=ifelse(is.na(menopause.age), 1, 0)) %>%
  mutate(ex.menopausal2=ifelse(ex.cancer2==1|ex.hyst2==1|ex.death2==1|ex.questionnaire2==1, 0, ex.menopausal2)) %>%
  mutate(ex.questionnaire2=ifelse(ex.cancer2==0 & ex.hyst2==0 & ex.death2==0 & ex.questionnaire2==0 & included==0, 1, ex.questionnaire2)) %>%
  mutate(ex2=ifelse(ex.cancer2==1|ex.hyst2==1|ex.menopausal2==1|ex.death2==1|ex.questionnaire2==1, 1, 0)) %>%
  mutate(in2=1-ex2) %>%
  ungroup() %>%
  summarize_at(vars(ex2, ex.cancer2, ex.hyst2, ex.death2, ex.questionnaire2, ex.menopausal2, in2), sum)

# Determining exclusions over follow-up
nhsii.inclusion3 <- data.frame(ex.cancer2=nhsii.inclusion2$in2-nrow(nhsii.final)) %>%
  mutate(ex2=ex.cancer2) %>%
  mutate(in2=-ex2)

# Combining pre- and post-baseline exclusions
t1 <- nhsii.inclusion
t2 <- nhsii.inclusion2 %>% bind_rows(nhsii.inclusion3) %>% summarize_all(sum, na.rm=T)
sf1c <- bind_cols(t1, t2)

write.csv(sf1c, "./Output/Figure_S1c.csv", row.names=F, na="")


# --- c. Figure S3 ---

# PLCO exclusions at baseline
plco.inclusion2 <- plco.long2 %>%
  mutate(ex=ifelse(observed.followup<=0|time.of.onset<0.5, 1, 0)) %>%
  mutate(ex.hyst=ifelse(hyst.age<=study.entry.age & ex==1, 1, 0)) %>%
  mutate(ex.cancer=ifelse(ex==1 & ex.hyst==0, 1, 0)) %>%
  mutate(ex.genetic=ifelse(ex==0 & is.na(genetic), 1, 0)) %>%
  summarize_at(vars(ex, ex.cancer, ex.hyst, ex.genetic), sum) %>%
  mutate(include=-ex)

# Summarizing exclusions
sf1d <- plco.inclusion %>% bind_rows(plco.inclusion2) %>% summarize_all(sum, na.rm=T) %>%
  mutate(genetic.include=include-ex.genetic)

write.csv(sf1d, "./Output/Figure_S1d.csv", row.names=F, na="")

# --- d. Figure S4 ---

# Counting number of participants in complete-case analysis
data <- e2c2.nomiss[,c("casecon", "site", predictors)] %>% drop_na

# Summarizing number of people that were excluded for each criteria
exclusions <- e2c2.exclusions %>% 
  arrange(-casecon) %>%
  mutate(missing.exclude=include-rev(table(data$casecon))) %>%
  mutate(final.include=rev(table(data$casecon)))

write.csv(exclusions, "./Output/Figure_S1a.csv", row.names=F, na="")

# --- e. Figures S5-6

  # --- i. Pulling out estimates from validation results ---
  
  t <- rbind(base.endo.nhsi[,2] %*% t(nhsi10.clinical$Category_Specific_Predicted_Relative_Risk),
             base.endo.nhsii[,2] %*% t(nhsii10.clinical$Category_Specific_Predicted_Relative_Risk),
             base.endo.plco[,2] %*% t(plco10.clinical$Category_Specific_Predicted_Relative_Risk),
             base.endo.plco[,2] %*% t(plco10.genetic1$Category_Specific_Predicted_Relative_Risk)) %>%
    data.frame() %>%
    mutate(Average=c(base.endo.nhsi[,2], base.endo.nhsii[,2], 
                     base.endo.plco[,2], base.endo.plco[,2])) %>%
    mutate(Age=rep(40:85, 4)) %>%
    mutate(Study=c(rep("NHS", length(40:85)), rep("NHSII", length(40:85)), 
                   rep("PLCO", length(40:85)), rep("PLCO-G", length(40:85)))) %>%
    # Converting data from wide to long
    pivot_longer(cols=c("Average", starts_with("X")), values_to="Risk", names_to="Decile") %>%
    mutate(Decile=ifelse(Decile!="Average", substr(Decile, 2, nchar(Decile)), Decile)) %>%
    filter(Age>=start.age) %>%  
    arrange(Study, Decile, Age) %>%
    # Calculating 10-year risk and cumulative risk
    mutate(Risk=Risk) %>%
    group_by(Study, Decile) %>%
    mutate(Survival=1-Risk) %>%
    mutate(CumSurvival=cumprod(Survival)) %>%
    mutate(RiskCum=1-CumSurvival) %>%
    mutate(Risk10=1-(Survival*lead(Survival,1)*lead(Survival,2)*lead(Survival,3)*
                       lead(Survival,4)*lead(Survival,5)*lead(Survival,6)*lead(Survival,7)*
                       lead(Survival,8)*lead(Survival,9))) %>%
    ungroup() %>%
    # Converting Decile variable into factor
    mutate(Decile=factor(Decile)) %>%
    mutate(Decile=fct_relevel(Decile, "10", after=9))

  # --- ii. Creating legend for plot ---
  
  legend <- t %>% 
    filter(Study=="NHS") %>%
    ggplot(aes(x=Age, y=RiskCum, linetype=Decile, color=Decile)) +
    geom_line(size=0.75) +
    scale_x_continuous("\nAge") +
    scale_y_continuous("Cumulative Absolute Risk\n") +
    scale_color_manual(values=c(brewer.pal(n = 10, name = "Spectral"), "black")) +
    scale_linetype_manual(values=c(rep("solid", 10), "dotted")) + 
    theme_bw() +
    theme(text=element_text(size=9), legend.key.height = unit(0.5, "cm"))
  legend <- get_legend(legend)

  
  # --- iii. Figure S5 ---
  
  fig1a <- t %>% 
    filter(Study=="NHS") %>%
    ggplot(aes(x=Age, y=RiskCum, linetype=Decile, color=Decile)) +
    geom_line(size=0.75) +
    scale_x_continuous(" \nAge") +
    scale_y_continuous("Cumulative Absolute Risk\n ", limits=c(0, 0.15)) +
    scale_color_manual(values=c(brewer.pal(n = 10, name = "Spectral"), "black")) +
    scale_linetype_manual(values=c(rep("solid", 10), "dotted")) + 
    theme_bw() +
    theme(legend.position="none", text=element_text(size=9))
  
  fig1b <- t %>% 
    filter(Study=="NHS" & Age<=75) %>%
    ggplot(aes(x=Age, y=Risk10, linetype=Decile, color=Decile)) +
    geom_line(size=0.75) +
    scale_x_continuous(" \nAge") +
    scale_y_continuous("Absolute 10-Year Risk\n ", limits=c(0, 0.053), 
                       labels=seq(0, 0.05, 0.01), breaks=seq(0, 0.05, 0.01)) +
    scale_color_manual(values=c(brewer.pal(n = 10, name = "Spectral"), "black")) +
    scale_linetype_manual(values=c(rep("solid", 10), "dotted")) + 
    theme_bw() +
    theme(legend.position="none", text=element_text(size=9))
  
  fig1c <- t %>% 
    filter(Study=="NHSII") %>%
    ggplot(aes(x=Age, y=RiskCum, linetype=Decile, color=Decile)) +
    geom_line(size=0.75) +
    scale_x_continuous(" \nAge") +
    scale_y_continuous("Cumulative Absolute Risk\n ", limits=c(0, 0.15)) +
    scale_color_manual(values=c(brewer.pal(n = 10, name = "Spectral"), "black")) +
    scale_linetype_manual(values=c(rep("solid", 10), "dotted")) + 
    theme_bw() +
    theme(legend.position="none", text=element_text(size=9))
  
  fig1d <- t %>% 
    filter(Study=="NHSII" & Age<=75) %>%
    ggplot(aes(x=Age, y=Risk10, linetype=Decile, color=Decile)) +
    geom_line(size=0.75) +
    scale_x_continuous(" \nAge") +
    scale_y_continuous("Absolute 10-Year Risk\n ", limits=c(0, 0.053), 
                       labels=seq(0, 0.05, 0.01), breaks=seq(0, 0.05, 0.01)) +
    scale_color_manual(values=c(brewer.pal(n = 10, name = "Spectral"), "black")) +
    scale_linetype_manual(values=c(rep("solid", 10), "dotted")) + 
    theme_bw() +
    theme(legend.position="none", text=element_text(size=9))
  
  plots <- plot_grid(fig1a, fig1b, fig1c, fig1d, 
                     nrow=2,
                     labels=c("A", " ", "B", " "),
                     label_size=10)
  
  fig1 <- plot_grid(plots, legend, nrow=1, rel_widths=c(8/9, 1/9))
  
  ggsave("./Output/Figure_S5.pdf", plot=fig1, device="pdf", width=7, height=6, units="in")  

  
  # --- iv. Figure S6 ---
  
  fig2a <- t %>% 
    filter(Study=="PLCO") %>%
    ggplot(aes(x=Age, y=RiskCum, linetype=Decile, color=Decile)) +
    geom_line(size=0.75) +
    scale_x_continuous(" \nAge") +
    scale_y_continuous("Cumulative Absolute Risk\n ", limits=c(0, 0.15)) +
    scale_color_manual(values=c(brewer.pal(n = 10, name = "Spectral"), "black")) +
    scale_linetype_manual(values=c(rep("solid", 10), "dotted")) + 
    theme_bw() +
    theme(legend.position="none", text=element_text(size=9))
  
  fig2b <- t %>% 
    filter(Study=="PLCO" & Age<=75) %>%
    ggplot(aes(x=Age, y=Risk10, linetype=Decile, color=Decile)) +
    geom_line(size=0.75) +
    scale_x_continuous(" \nAge") +
    scale_y_continuous("Absolute 10-Year Risk\n ", limits=c(0, 0.05)) +
    scale_color_manual(values=c(brewer.pal(n = 10, name = "Spectral"), "black")) +
    scale_linetype_manual(values=c(rep("solid", 10), "dotted")) + 
    theme_bw() +
    theme(legend.position="none", text=element_text(size=9))
  
  fig2c <- t %>% 
    filter(Study=="PLCO-G") %>%
    ggplot(aes(x=Age, y=RiskCum, linetype=Decile, color=Decile)) +
    geom_line(size=0.75) +
    scale_x_continuous(" \nAge") +
    scale_y_continuous("Cumulative Absolute Risk\n ", limits=c(0, 0.15)) +
    scale_color_manual(values=c(brewer.pal(n = 10, name = "Spectral"), "black")) +
    scale_linetype_manual(values=c(rep("solid", 10), "dotted")) + 
    theme_bw() +
    theme(legend.position="none", text=element_text(size=9))
  
  fig2d <- t %>% 
    filter(Study=="PLCO-G" & Age<=75) %>%
    ggplot(aes(x=Age, y=Risk10, linetype=Decile, color=Decile)) +
    geom_line(size=0.75) +
    scale_x_continuous(" \nAge") +
    scale_y_continuous("Absolute 10-Year Risk\n ", limits=c(0, 0.05)) +
    scale_color_manual(values=c(brewer.pal(n = 10, name = "Spectral"), "black")) +
    scale_linetype_manual(values=c(rep("solid", 10), "dotted")) + 
    theme_bw() +
    theme(legend.position="none", text=element_text(size=9))
  
  plots <- plot_grid(fig2a, fig2b, fig2c, fig2d, 
                     nrow=2,
                     labels=c("A", " ", "B", " "),
                     label_size=10)
  
  fig2 <- plot_grid(plots, legend, nrow=1, rel_widths=c(8/9, 1/9))
  
  ggsave("./Output/Figure_S6.pdf", plot=fig2, device="pdf", width=7, height=6, units="in")  


  # --- v. Saving values for range of cumulative risks ---
  
  t.save <- rbind(t[which(t$Decile==1 & t$Age==85),],
                  t[which(t$Decile==10 & t$Age==85),])
  write.csv(t.save, "./Output/Manuscript_Text_CumRisk_Range.csv", row.names=F)


# ---------------------------------------------------------------
# ----------------- (15) Supplementary Tables -------------------
# ---------------------------------------------------------------


# --- a. Table S1 ---
e2c2.id <- e2c2.nomiss[,c("id", predictors)] %>% drop_na %>% pull(id)

data <- e2c2.all %>% filter(id %in% e2c2.id)

# Summarizing age by study
t.age <- data %>% 
  group_by(site) %>%
  summarize_at("age", c(mean, sd)) %>%
  mutate(age=sprintf("%.1f \U00B1 %.1f", fn1, fn2)) %>%
  select(site, age)

# Calculating number of cases and controls by study
t.n <- data %>%
  group_by(site, casecon) %>%
  tally() %>%
  ungroup() %>%
  spread(casecon, n, sep="")

# Combining summary statistics
st1 <- full_join(t.age, t.n, by="site") %>%
  select(site, age, casecon1, casecon0) %>%
  mutate(site=ifelse(site==211, 011, site)) %>%
  arrange(site)

write.xlsx(st1, "./Output/Table_S1.xlsx", rowNames=F, overwrite=T)


# --- b. Table S2 ---

t <- e2c2.all %>%
  group_by(site) %>%
  
  # Determining if ever measured, stratified by study site
  summarize_at(
    c("education", "smoking", "packyears", "alcohol", "bmi", 
    "parity", "firstbirthage", "lastbirthage", "menarchage", 
    "anyhrtever", "ephrtever", "ehrtever",
    "ephrtdur", "ehrtdur", "ocever", "ocdur",
    "diabetes", "endom", "hypertension", "famhx"), 
    ~ifelse(sum(!is.na(.))>0, "\U2713", NA)) %>%
  
  # Transposing table
  gather("var", "check", -site) %>%
  spread(site, check, sep="") %>%
  mutate(total=rowSums(!is.na(dplyr::select(., -var)))) %>%
  rbind(data.frame(t(colSums(!is.na(.[,]))))) %>%
  mutate(var=ifelse(var=="20", "total", var)) %>%
  mutate(total=ifelse(var=="total", NA, total))    

# Ordering variables
t.var <- c(
  "demographic", "education", 
  "lifestyle", "smoking", "packyears", "alcohol", "bmi", 
  "repro", "parity", "firstbirthage", "lastbirthage", "menarchage", 
    "hrtever", "anyhrtever", "ephrtever", "ehrtever",
  "hrtdur", "ephrtdur", "ehrtdur", "ocever", "ocdur",
  "clinical", "diabetes", "endom", "hypertension", "famhx",
  "total")
t.var <- data.frame(cbind(var=t.var, order=seq(1:length(t.var))), stringsAsFactors=F)

# Renaming variables to each study name
st2 <- t %>%
  full_join(t.var, by="var") %>%
  mutate(order=as.numeric(order)) %>%
  arrange(order) %>%
  rename(EDGE=site101) %>%
  rename(FHCRC=site104) %>%
  rename(WISE=site105) %>%
  rename(HAW=site106) %>%
  rename(PCCS=site109) %>%
  rename(CECS=site111) %>%
  rename(USECCS=site202) %>%
  rename(AECPA=site204) %>%
  rename(BAWHS=site205) %>%
  rename('USC/LA'=site206) %>%
  rename(ANECS=site207) %>%
  rename(RPEDS=site210) %>%
  rename(SV=site303) %>%
  rename(WNYDS=site304) %>%
  rename('Milano 1'=site305) %>%
  rename('Milano 2'=site306) %>%
  rename(ECDBWGEI=site307) %>%
  rename(IMS=site405) %>%
  rename(Screenwide=site406) %>%
  select(var, everything())

write.xlsx(st2, './Output/Table_S2.xlsx', rowNames=F, na="", overwrite=T)


# --- c. Table S3 ---

e2c2.id <- e2c2.nomiss[,c("id", predictors)] %>% drop_na %>% pull(id)

# Table for participants included in the analysis
  
  # Set duration missing if never used
  data <- e2c2.all %>% 
    filter(id %in% e2c2.id) %>%
    mutate(ephrtdur1=ifelse(ephrtdur!=0, ephrtdur, NA)) %>%
    mutate(ehrtdur1=ifelse(ehrtdur!=0, ehrtdur, NA)) %>%
    mutate(ocdur1=ifelse(ocdur!=0, ocdur, NA)) %>%
    mutate(packyears1=ifelse(packyears!=0, packyears, NA))
  
  # Calculating mean/SD for continuous variables
  t.continuous1 <- data %>%
    group_by(casecon) %>%
    summarize_at(vars(age, bmi, firstbirthage, lastbirthage, menarchage), 
                 list(mean=mean, sd=sd), na.rm=T) %>%
    ungroup() %>%
    gather(var, value, -casecon) %>%
    separate(var, c("var", "measure")) %>%
    spread(measure, value) %>%
    mutate(summary=sprintf("%.1f \U00B1 %.1f", mean, sd)) %>%
    select(-mean, -sd) %>%
    spread(casecon, summary, sep="")
  
  # Calculating median/IQR for continuous variables
  t.continuous2 <- data %>%
    group_by(casecon) %>%
    summarize_at(vars(pregs, packyears1, parity, alcohol, ephrtdur1, ehrtdur1, ocdur1), 
                 list(median=median, ~quantile(., probs=0.25, na.rm=T), 
                      ~quantile(., probs=0.75, na.rm=T)), na.rm=T) %>%
    ungroup() %>%
    gather(var, value, -casecon) %>%
    separate(var, c("var", "measure"), sep="_") %>%
    spread(measure, value) %>%
    mutate(summary=sprintf("%.1f (%.1f, %.1f)", median, quantile..2, quantile..3)) %>%
    select(casecon, var, summary) %>%
    spread(casecon, summary, sep="")
  
  # Calculating frequency/percentage for categorical variables
  t.binary <- data %>%
    group_by(casecon) %>%
    summarize_at(vars(hispanic, 
                      education1, education2, education3, 
                      smoking1, smoking2, smoking3, 
                      grade1, grade2, grade3, grade4, 
                      anyhrtever, ephrtever, ehrtever, ocever, 
                      diabetes, endom, hypertension, famhx),
                 list(freq=sum, percent=mean), na.rm=T) %>%
    ungroup() %>%
    gather(var, value, -casecon) %>%
    separate(var, c("var", "measure")) %>%
    spread(measure, value) %>%
    mutate(summary=sprintf("%s (%.1f)", format(freq, big.mark=","), percent*100)) %>%
    select(-freq, -percent) %>%
    spread(casecon, summary, sep="")
  
  # Calculating n/% missing for each variable
  t.n <- data %>%
    group_by(casecon) %>%
    summarize_at(vars(id, age, bmi, firstbirthage, lastbirthage, menarchage,
                      pregs, packyears1, parity, alcohol, ephrtdur1, ehrtdur1, ocdur1,
                      hispanic, education, smoking, 
                      grade, anyhrtever, ephrtever, ehrtever, ocever, 
                      diabetes, endom, hypertension, famhx),
                 list(~sum(!is.na(.))), na.rm=T) %>%
    gather(var, n, id:famhx) %>%
    spread(casecon, n) %>%
    rename(n0="0") %>%
    rename(n1="1") %>%
    mutate(var=ifelse(var=="id", "n", var)) %>%
    mutate(missing0=100-n0/nrow(data[(data$casecon==0),])*100) %>%
    mutate(missing1=100-n1/nrow(data[(data$casecon==1),])*100) %>%
    mutate(missing0=ifelse(var=="ephrtdur1", 100-n0/sum(data[which(data$casecon==0),"ephrtever"]==1, na.rm=T)*100, missing0)) %>%
    mutate(missing1=ifelse(var=="ephrtdur1", 100-n1/sum(data[which(data$casecon==1),"ephrtever"]==1, na.rm=T)*100, missing1)) %>%
    mutate(missing0=ifelse(var=="ehrtdur1", 100-n0/sum(data[which(data$casecon==0),"ehrtever"]==1, na.rm=T)*100, missing0)) %>%
    mutate(missing1=ifelse(var=="ehrtdur1", 100-n1/sum(data[which(data$casecon==1),"ehrtever"]==1, na.rm=T)*100, missing1)) %>%
    mutate(missing0=ifelse(var=="ocdur1", 100-n0/sum(data[which(data$casecon==0),"ocever"]==1, na.rm=T)*100, missing0)) %>%
    mutate(missing1=ifelse(var=="ocdur1", 100-n1/sum(data[which(data$casecon==1),"ocever"]==1, na.rm=T)*100, missing1)) %>%
    mutate(missing0=ifelse(var=="packyears1", 100-n0/nrow(data[which(data$casecon==0 & (data$smoking==2|data$smoking==3)),])*100, missing0)) %>%
    mutate(missing1=ifelse(var=="packyears1", 100-n1/nrow(data[which(data$casecon==1 & (data$smoking==2|data$smoking==3)),])*100, missing1)) %>%
    mutate(n0=ifelse(var=="n", sprintf("n = %s", n0), sprintf("%s (%.1f)", format(n0, big.mark=","), missing0))) %>%
    mutate(n1=ifelse(var=="n", sprintf("n = %s", n1), sprintf("%s (%.1f)", format(n1, big.mark=","), missing1))) %>%
    select(-missing0, -missing1)

# Table for participants NOT included in the analysis
  
  # Set duration missing if never used
  data2 <- e2c2.all %>% 
    filter(!(id %in% e2c2.id)) %>%
    mutate(ephrtdur1=ifelse(ephrtdur!=0, ephrtdur, NA)) %>%
    mutate(ehrtdur1=ifelse(ehrtdur!=0, ehrtdur, NA)) %>%
    mutate(ocdur1=ifelse(ocdur!=0, ocdur, NA)) %>%
    mutate(packyears1=ifelse(packyears!=0, packyears, NA))
  
  # Calculating mean/SD for continuous variables
  t2.continuous1 <- data2 %>%
    group_by(casecon) %>%
    summarize_at(vars(age, bmi, firstbirthage, lastbirthage, menarchage), 
                 list(mean=mean, sd=sd), na.rm=T) %>%
    ungroup() %>%
    gather(var, value, -casecon) %>%
    separate(var, c("var", "measure")) %>%
    spread(measure, value) %>%
    mutate(summary=sprintf("%.1f \U00B1 %.1f", mean, sd)) %>%
    select(-mean, -sd) %>%
    spread(casecon, summary, sep="")
  
  # Calculating median/IQR for continuous variables
  t2.continuous2 <- data2 %>%
    group_by(casecon) %>%
    summarize_at(vars(pregs, packyears1, parity, alcohol, ephrtdur1, ehrtdur1, ocdur1), 
                 list(median=median, ~quantile(., probs=0.25, na.rm=T), 
                      ~quantile(., probs=0.75, na.rm=T)), na.rm=T) %>%
    ungroup() %>%
    gather(var, value, -casecon) %>%
    separate(var, c("var", "measure"), sep="_") %>%
    spread(measure, value) %>%
    mutate(summary=sprintf("%.1f (%.1f, %.1f)", median, quantile..2, quantile..3)) %>%
    select(casecon, var, summary) %>%
    spread(casecon, summary, sep="")
  
  # Calculating frequency/percentage for categorical variables
  t2.binary <- data2 %>%
    group_by(casecon) %>%
    summarize_at(vars(hispanic, 
                      education1, education2, education3, 
                      smoking1, smoking2, smoking3, 
                      grade1, grade2, grade3, grade4, 
                      anyhrtever, ephrtever, ehrtever, ocever, 
                      diabetes, endom, hypertension, famhx),
                 list(freq=sum, percent=mean), na.rm=T) %>%
    ungroup() %>%
    gather(var, value, -casecon) %>%
    separate(var, c("var", "measure")) %>%
    spread(measure, value) %>%
    mutate(summary=sprintf("%s (%.1f)", format(freq, big.mark=","), percent*100)) %>%
    select(-freq, -percent) %>%
    spread(casecon, summary, sep="")
  
  # Calculating n/% missing for each variable
  t2.n <- data2 %>%
    group_by(casecon) %>%
    summarize_at(vars(id, age, bmi, firstbirthage, lastbirthage, menarchage,
                      pregs, packyears1, parity, alcohol, ephrtdur1, ehrtdur1, ocdur1,
                      hispanic, education, smoking, 
                      grade, anyhrtever, ephrtever, ehrtever, ocever, 
                      diabetes, endom, hypertension, famhx),
                 list(~sum(!is.na(.))), na.rm=T) %>%
    gather(var, n, id:famhx) %>%
    spread(casecon, n) %>%
    rename(n0="0") %>%
    rename(n1="1") %>%
    mutate(var=ifelse(var=="id", "n", var)) %>%
    mutate(missing0=100-n0/nrow(data2[(data2$casecon==0),])*100) %>%
    mutate(missing1=100-n1/nrow(data2[(data2$casecon==1),])*100) %>%
    mutate(missing0=ifelse(var=="ephrtdur1", 100-n0/sum(data2[which(data2$casecon==0),"ephrtever"]==1, na.rm=T)*100, missing0)) %>%
    mutate(missing1=ifelse(var=="ephrtdur1", 100-n1/sum(data2[which(data2$casecon==1),"ephrtever"]==1, na.rm=T)*100, missing1)) %>%
    mutate(missing0=ifelse(var=="ehrtdur1", 100-n0/sum(data2[which(data2$casecon==0),"ehrtever"]==1, na.rm=T)*100, missing0)) %>%
    mutate(missing1=ifelse(var=="ehrtdur1", 100-n1/sum(data2[which(data2$casecon==1),"ehrtever"]==1, na.rm=T)*100, missing1)) %>%
    mutate(missing0=ifelse(var=="ocdur1", 100-n0/sum(data2[which(data2$casecon==0),"ocever"]==1, na.rm=T)*100, missing0)) %>%
    mutate(missing1=ifelse(var=="ocdur1", 100-n1/sum(data2[which(data2$casecon==1),"ocever"]==1, na.rm=T)*100, missing1)) %>%
    mutate(missing0=ifelse(var=="packyears1", 100-n0/nrow(data2[which(data2$casecon==0 & (data2$smoking==2|data2$smoking==3)),])*100, missing0)) %>%
    mutate(missing1=ifelse(var=="packyears1", 100-n1/nrow(data2[which(data2$casecon==1 & (data2$smoking==2|data2$smoking==3)),])*100, missing1)) %>%
    mutate(n0=ifelse(var=="n", sprintf("n = %s", n0), sprintf("%s (%.1f)", format(n0, big.mark=","), missing0))) %>%
    mutate(n1=ifelse(var=="n", sprintf("n = %s", n1), sprintf("%s (%.1f)", format(n1, big.mark=","), missing1))) %>%
    select(-missing0, -missing1)  
  
# Ordering rows of tables
t.var <- c(
  "n",
  "demographic",
  "age",
  "hispanic", 
  "education", "education1", "education2", "education3", 
  "lifestyle",
  "smoking", "smoking1", "smoking2", "smoking3", "packyears1",
  "alcohol",
  "bmi", 
  "repro",
  "parity", 
  "firstbirthage", "lastbirthage", "menarchage", 
  "hrtever", "anyhrtever", "ephrtever", "ehrtever",
  "hrtdur", "ephrtdur1", "ehrtdur1",
  "ocever", "ocdur1",
  "clinical",
  "diabetes",
  "endom", 
  "hypertension", 
  "famhx")
t.var <- data.frame(cbind(var=t.var, order=seq(1:length(t.var))), stringsAsFactors=F)

# Combining all summary statistics
st3a <- bind_rows(t.continuous1, t.continuous2, t.binary) %>%
  full_join(t.n, by="var") %>%
  full_join(t.var, by="var") %>%
  mutate(order=as.numeric(order)) %>%
  arrange(order) %>%
  filter(!is.na(order)) %>%
  mutate(n0=trimws(n0)) %>%
  mutate(casecon0=trimws(casecon0)) %>%
  mutate(n1=trimws(n1)) %>%
  mutate(casecon1=trimws(casecon1)) %>%
  select(var, n1, casecon1, n0, casecon0) %>%
  data.frame()

st3b <- bind_rows(t2.continuous1, t2.continuous2, t2.binary) %>%
  full_join(t2.n, by="var") %>%
  full_join(t.var, by="var") %>%
  mutate(order=as.numeric(order)) %>%
  arrange(order) %>%
  filter(!is.na(order)) %>%
  mutate(n0=trimws(n0)) %>%
  mutate(casecon0=trimws(casecon0)) %>%
  mutate(n1=trimws(n1)) %>%
  mutate(casecon1=trimws(casecon1)) %>%
  select(var, n1, casecon1, n0, casecon0) %>%
  rename(n1.missing=n1) %>%
  rename(casecon1.missing=casecon1) %>%
  rename(n0.missing=n0) %>%
  rename(casecon0.missing=casecon0) %>%
  data.frame()

st3 <- full_join(st3a, st3b, by="var")

write.xlsx(st3, "./Output/Table_S3.xlsx", na="", rowNames=F, overwrite=T)


# --- d. Table S4 ---

# NHS I: all eligible participants
data1.all <- nhsi.final

  # Summarizing continuous variables
  t1.continuous <- data1.all %>%
    summarize_at(vars(age, year, firstbirthage), 
                 list(median=median, min=min, max=max, mean=mean, sd=sd), na.rm=T) %>%
    pivot_longer(cols=everything(), names_to="var", values_to="value") %>%
    separate(var, c("var", "measure"), sep="_") %>%
    pivot_wider(id_cols=var, names_from="measure", values_from="value") %>%
    mutate(summary=sprintf("%.0f (%.0f, %.0f)", median, min, max)) %>%
    mutate(summary=ifelse(var=="firstbirthage", 
      sprintf("%.1f \U00B1 %.1f", mean, sd), summary)) %>%
    select(var, summary) %>%
    rename(summary1=summary)  

  # Summarizing categorical variables
  t1.binary <- data1.all %>%
    summarize_at(vars(
      education1, education2, education3, 
      smoking1, smoking2, smoking3, 
      bmi.group1, bmi.group2, bmi.group3, bmi.group4, bmi.group5,
      parity.group0, parity.group1, parity.group2, parity.group3, parity.group4,
      firstbirthage.group1, firstbirthage.group2, firstbirthage.group3, 
      firstbirthage.group4, firstbirthage.group5, firstbirthage.group9,
      menarchage2.group9, menarchage2.group10, menarchage2.group12, 
      menarchage2.group14, menarchage2.group16,
      anyhrtever, 
      ehrtever1, ehrtdur.group0, ehrtdur.group1, ehrtdur.group2, ehrtdur.group3,
      ephrtever1, ephrtdur.group0, ephrtdur.group1, ephrtdur.group2, ephrtdur.group3,
      ocever, ocdur.group0,  ocdur.group1,  ocdur.group2,  ocdur.group3,
      diabetes1, hypertension1),
      list(freq=sum, percent=mean), na.rm=T) %>%
  ungroup() %>%
  gather(var, value) %>%
  separate(var, c("var", "measure"), sep="_") %>%
  spread(measure, value) %>%
  mutate(summary=sprintf("%s (%.1f)", format(freq, big.mark=","), percent*100)) %>%
  mutate(summary=trimws(summary, which="both")) %>%
  select(var, summary) %>%
  rename(summary1=summary)  

  # Obtaining n and % missing for each variable
  t1.n <- data.frame(var="n", summary1=sprintf("%s", format(nrow(data1.all), big.mark=","))) %>% 
    mutate_all(as.character)


# NHS I: all eligible participants with genetic data
data1.g <- nhsi.final %>% filter(gwas==1)

  # Summarizing continuous variables
  t2.continuous <- data1.g %>%
    summarize_at(vars(age, year, firstbirthage), 
                 list(median=median, min=min, max=max, mean=mean, sd=sd), na.rm=T) %>%
    pivot_longer(cols=everything(), names_to="var", values_to="value") %>%
    separate(var, c("var", "measure"), sep="_") %>%
    pivot_wider(id_cols=var, names_from="measure", values_from="value") %>%
    mutate(summary=sprintf("%.0f (%.0f, %.0f)", median, min, max)) %>%
    mutate(summary=ifelse(var=="firstbirthage", 
      sprintf("%.1f \U00B1 %.1f", mean, sd), summary)) %>%
    select(var, summary) %>%
    rename(summary2=summary)  

  # Summarizing categorical variables
  t2.binary <- data1.g %>%
    summarize_at(vars(
      education1, education2, education3, 
      smoking1, smoking2, smoking3, 
      bmi.group1, bmi.group2, bmi.group3, bmi.group4, bmi.group5,
      parity.group0, parity.group1, parity.group2, parity.group3, parity.group4,
      firstbirthage.group1, firstbirthage.group2, firstbirthage.group3, 
      firstbirthage.group4, firstbirthage.group5, firstbirthage.group9,
      menarchage2.group9, menarchage2.group10, menarchage2.group12, 
      menarchage2.group14, menarchage2.group16,
      anyhrtever, 
      ehrtever1, ehrtdur.group0, ehrtdur.group1, ehrtdur.group2, ehrtdur.group3,
      ephrtever1, ephrtdur.group0, ephrtdur.group1, ephrtdur.group2, ephrtdur.group3,
      ocever, ocdur.group0,  ocdur.group1,  ocdur.group2,  ocdur.group3,
      diabetes1, hypertension1),
      list(freq=sum, percent=mean), na.rm=T) %>%
    ungroup() %>%
    gather(var, value) %>%
    separate(var, c("var", "measure"), sep="_") %>%
    spread(measure, value) %>%
    mutate(summary=sprintf("%s (%.1f)", format(freq, big.mark=","), percent*100)) %>%
    mutate(summary=trimws(summary, which="both")) %>%
    select(var, summary) %>%
    rename(summary2=summary)  

  # Obtaining n and % missing for each variable
  t2.n <- data.frame(var="n", summary2=sprintf("%s", 
    format(nrow(data1.g), big.mark=","))) %>% 
    mutate_all(as.character)


# NHS I: all eligible participants 
data2.all <- nhsii.final
  
  # Summarizing continuous variables
  t3.continuous <- data2.all %>%
    summarize_at(vars(age, year, firstbirthage), 
                 list(median=median, min=min, max=max, mean=mean, sd=sd), na.rm=T) %>%
    pivot_longer(cols=everything(), names_to="var", values_to="value") %>%
    separate(var, c("var", "measure"), sep="_") %>%
    pivot_wider(id_cols=var, names_from="measure", values_from="value") %>%
    mutate(summary=sprintf("%.0f (%.0f, %.0f)", median, min, max)) %>%
    mutate(summary=ifelse(var=="firstbirthage", sprintf("%.1f \U00B1 %.1f", mean, sd), summary)) %>%
    select(var, summary) %>%
    rename(summary3=summary)  
  
  # Summarizing categorical variables
  t3.binary <- data2.all %>%
    summarize_at(vars(
      education1, education2, education3, 
      smoking1, smoking2, smoking3, 
      bmi.group1, bmi.group2, bmi.group3, bmi.group4, bmi.group5,
      parity.group0, parity.group1, parity.group2, parity.group3, parity.group4,
      firstbirthage.group1, firstbirthage.group2, firstbirthage.group3, 
      firstbirthage.group4, firstbirthage.group5, firstbirthage.group9,
      menarchage2.group9, menarchage2.group10, menarchage2.group12, 
      menarchage2.group14, menarchage2.group16,
      anyhrtever, 
      ehrtever1, ehrtdur.group0, ehrtdur.group1, ehrtdur.group2, ehrtdur.group3,
      ephrtever1, ephrtdur.group0, ephrtdur.group1, ephrtdur.group2, ephrtdur.group3,
      ocever, ocdur.group0,  ocdur.group1,  ocdur.group2,  ocdur.group3,
      diabetes1, hypertension1),
      list(freq=sum, percent=mean), na.rm=T) %>%
    ungroup() %>%
    gather(var, value) %>%
    separate(var, c("var", "measure"), sep="_") %>%
    spread(measure, value) %>%
    mutate(summary=sprintf("%s (%.1f)", format(freq, big.mark=","), percent*100)) %>%
    mutate(summary=trimws(summary, which="both")) %>%
    select(var, summary) %>%
    rename(summary3=summary)  
  
  # Obtaining n and % missing for each variable
  t3.n <- data.frame(var="n", summary3=sprintf("%s", format(nrow(data2.all), big.mark=","))) %>% 
    mutate_all(as.character)

# PLCO: all eligible participants
data3.all <- plco.final

  # Summarizing continuous variables
  t5.continuous <- data3.all %>%
    summarize_at(vars(age, year), 
                 list(median=median, min=min, max=max), na.rm=T) %>%
    pivot_longer(cols=everything(), names_to="var", values_to="value") %>%
    separate(var, c("var", "measure"), sep="_") %>%
    pivot_wider(id_cols=var, names_from="measure", values_from="value") %>%
    mutate(summary=sprintf("%.0f (%.0f, %.0f)", median, min, max)) %>%
    select(var, summary) %>%
    rename(summary5=summary)  

  # Summarizing categorical variables
  t5.binary <- data3.all %>%
    summarize_at(vars(
      education1, education2, education3, 
      smoking1, smoking2, smoking3, 
      bmi.group1, bmi.group2, bmi.group3, bmi.group4, bmi.group5,
      parity.group0, parity.group1, parity.group2, parity.group3, parity.group4,
      firstbirthage.group1, firstbirthage.group2, firstbirthage.group3, 
      firstbirthage.group4, firstbirthage.group5, firstbirthage.group9,
      menarchage2.group9, menarchage2.group10, menarchage2.group12, 
      menarchage2.group14, menarchage2.group16,
      anyhrtever, 
      ehrtever1, ehrtdur.group0, ehrtdur.group1, ehrtdur.group2, ehrtdur.group3,
      ephrtever1, ephrtdur.group0, ephrtdur.group1, ephrtdur.group2, ephrtdur.group3,
      ocever, ocdur.group0,  ocdur.group1,  ocdur.group2,  ocdur.group3,
      diabetes1, hypertension1),
      list(freq=sum, percent=mean), na.rm=T) %>%
    ungroup() %>%
    gather(var, value) %>%
    separate(var, c("var", "measure"), sep="_") %>%
    spread(measure, value) %>%
    mutate(summary=sprintf("%s (%.1f)", format(freq, big.mark=","), percent*100)) %>%
    mutate(summary=trimws(summary, which="both")) %>%
    select(var, summary) %>%
    rename(summary5=summary)  

  # Obtaining n and % missing for each variable
  t5.n <- data.frame(var="n", summary5=sprintf("%s", format(nrow(data3.all), big.mark=","))) %>% 
    mutate_all(as.character)
  
  
# PLCO: all eligible participants with genetic data
data3.genetic <- plco.final %>% filter(genetic==1)
  
  # Summarizing continuous variables
  t6.continuous <- data3.genetic %>%
    summarize_at(vars(age, year), 
                 list(median=median, min=min, max=max), na.rm=T) %>%
    pivot_longer(cols=everything(), names_to="var", values_to="value") %>%
    separate(var, c("var", "measure"), sep="_") %>%
    pivot_wider(id_cols=var, names_from="measure", values_from="value") %>%
    mutate(summary=sprintf("%.0f (%.0f, %.0f)", median, min, max)) %>%
    select(var, summary) %>%
    rename(summary5=summary)  
  
  # Summarizing categorical variables
  t6.binary <- data3.genetic %>%
    summarize_at(vars(
      education1, education2, education3, 
      smoking1, smoking2, smoking3, 
      bmi.group1, bmi.group2, bmi.group3, bmi.group4, bmi.group5,
      parity.group0, parity.group1, parity.group2, parity.group3, parity.group4,
      firstbirthage.group1, firstbirthage.group2, firstbirthage.group3, 
      firstbirthage.group4, firstbirthage.group5, firstbirthage.group9,
      menarchage2.group9, menarchage2.group10, menarchage2.group12, 
      menarchage2.group14, menarchage2.group16,
      anyhrtever, 
      ehrtever1, ehrtdur.group0, ehrtdur.group1, ehrtdur.group2, ehrtdur.group3,
      ephrtever1, ephrtdur.group0, ephrtdur.group1, ephrtdur.group2, ephrtdur.group3,
      ocever, ocdur.group0,  ocdur.group1,  ocdur.group2,  ocdur.group3,
      diabetes1, hypertension1),
      list(freq=sum, percent=mean), na.rm=T) %>%
    ungroup() %>%
    gather(var, value) %>%
    separate(var, c("var", "measure"), sep="_") %>%
    spread(measure, value) %>%
    mutate(summary=sprintf("%s (%.1f)", format(freq, big.mark=","), percent*100)) %>%
    mutate(summary=trimws(summary, which="both")) %>%
    select(var, summary) %>%
    rename(summary5=summary)  
  
  # Obtaining n and % missing for each variable
  t6.n <- data.frame(var="n", summary5=sprintf("%s", format(nrow(data3.genetic), big.mark=","))) %>% 
    mutate_all(as.character)

# Combining all summary statistics       
t.var <- c(
  "n", "year", "age",
  "demographic",
  "education", "education1", "education2", "education3", 
  "lifestyle",
  "smoking", "smoking1", "smoking2", "smoking3",
  "bmi.group", "bmi.group1", "bmi.group2", "bmi.group3", "bmi.group4", "bmi.group5", 
  "repro",
  "parity.group", "parity.group0", "parity.group1", "parity.group2", 
  "parity.group3", "parity.group4",
  "firstbirthage.group", "firstbirthage.group1", "firstbirthage.group2", 
  "firstbirthage.group3", "firstbirthage.group4", "firstbirthage.group5",
  "firstbirthage.group9",
  "menarchage2.group", "menarchage2.group9", "menarchage2.group10", 
  "menarchage2.group12", "menarchage2.group14", "menarchage2.group16",
  "anyhrtever", 
  "ehrtever1", "ehrtdur.group", "ehrtdur.group0", "ehrtdur.group1", 
  "ehrtdur.group2", "ehrtdur.group3",
  "ephrtever1", "ephrtdur.group", "ephrtdur.group0", "ephrtdur.group1", 
  "ephrtdur.group2", "ephrtdur.group3",
  "ocever", "ocdur.group", "ocdur.group0", "ocdur.group1", "ocdur.group2", "ocdur.group3",
  "clinical",
  "diabetes1",
  "hypertension1")
t.var <- data.frame(cbind(var=t.var, order=seq(1:length(t.var))), stringsAsFactors=F)

st4 <- t.var %>%
  full_join(bind_rows(t1.continuous, t1.binary, t1.n), by="var") %>%
  full_join(bind_rows(t2.continuous, t2.binary, t2.n), by="var") %>%
  full_join(bind_rows(t3.continuous, t3.binary, t3.n), by="var") %>%
  full_join(bind_rows(t5.continuous, t5.binary, t5.n), by="var") %>%
  full_join(bind_rows(t6.continuous, t6.binary, t6.n), by="var") %>%
  data.frame %>%
  filter(!is.na(order)) %>%
  select(-order)

write.xlsx(st4, "./Output/Table_S4.xlsx", rowNames=F, na="", overwrite=T)


# --- e. Table S5 ---

nhsi.seeryear <- "1989-1993"
nhsii.seeryear <- "2003-2007"
plco.seeryear <- "1996-2000"
current.seeryear <- "2013-2017"

nhsi.cdcyear <- 1988
nhsii.cdcyear <- 2004
plco.cdcyear <- 1997
current.cdcyear <- 2017

nhsi.brfssyear <- 1988
nhsii.brfssyear <- c(2006, 2008)
plco.brfssyear <- c(1996, 1998) 
current.brfssyear <- c(2016, 2018)

# Pulling baseline data used for NHS I
t.nhsi <- base.endo.all %>% 
  filter(brfss.year %in% nhsi.brfssyear & seer.year==nhsi.seeryear) %>%
  group_by(age.group) %>%
  mutate(hystavg=mean(hyst)) %>%
  ungroup() %>%
  filter(brfss.year==nhsi.brfssyear[1]) %>%
  mutate(rate.corrected=sprintf("%.1f", rate/(1-hystavg))) %>%
  mutate(hystavg=sprintf("%.3f", hystavg)) %>%
  mutate(rate=sprintf("%.1f", rate)) %>%
  filter(age.group>=start.age) %>%
  select(age.group, rate, hystavg, rate.corrected) %>%
  full_join(base.competing.mortality, by=c("age.group"="Age")) %>%
  filter(CDC.Year==nhsi.cdcyear) %>%
  mutate(mortality.rate=sprintf("%.1f", Mortality.Rate*100000)) %>%
  select(age.group, rate, hystavg, rate.corrected, mortality.rate) %>%
  full_join(base.competing.cancer, by=c("age.group"="Age")) %>%
  filter(SEER.Year==nhsi.seeryear) %>%
  mutate(cancer.rate=sprintf("%.1f", Cancer.Rate*100000)) %>%
  select(-Cancer.Rate, -SEER.Year) %>%
  filter(!is.na(rate)) %>%
  add_row(rate=paste("SEER", nhsi.seeryear), hystavg=paste("BRFSS ", nhsi.brfssyear[1], "-", tail(nhsi.brfssyear,1), sep=""),
          mortality.rate=paste("CDC WONDER", nhsi.cdcyear), cancer.rate=paste("SEER", nhsi.seeryear)) %>%
  arrange(!is.na(age.group), age.group) %>%
  mutate(age.group=ifelse(!is.na(age.group), paste(age.group, "-", age.group+4, sep=""), NA)) %>%
  mutate(age.group=ifelse(is.na(age.group), "NHS", age.group))

# Pulling baseline data used for NHS II
t.nhsii <- base.endo.all %>% 
  filter(brfss.year %in% nhsii.brfssyear & seer.year==nhsii.seeryear) %>%
  group_by(age.group) %>%
  mutate(hystavg=mean(hyst)) %>%
  ungroup() %>%
  filter(brfss.year==nhsii.brfssyear[1]) %>%
  mutate(rate.corrected=sprintf("%.1f", rate/(1-hystavg))) %>%
  mutate(hystavg=sprintf("%.3f", hystavg)) %>%
  mutate(rate=sprintf("%.1f", rate)) %>%
  filter(age.group>=start.age) %>%
  select(age.group, rate, hystavg, rate.corrected) %>%
  full_join(base.competing.mortality, by=c("age.group"="Age")) %>%
  filter(CDC.Year==nhsii.cdcyear) %>%
  mutate(age.gr=cut(age.group, breaks=seq(40, 85, by=5), labels=seq(40, 80, by=5), right=F)) %>%
  group_by(age.gr) %>%
  mutate(mortality.rate=sprintf("%.1f", sum(Deaths)/sum(Population)*100000)) %>%
  ungroup() %>%
  full_join(base.competing.cancer, by=c("age.group"="Age")) %>%
  filter(SEER.Year==nhsii.seeryear) %>%
  mutate(cancer.rate=sprintf("%.1f", Cancer.Rate*100000)) %>%
  select(age.group, rate, hystavg, rate.corrected, mortality.rate, cancer.rate) %>%
  filter(!is.na(rate)) %>%
  add_row(rate=paste("SEER", nhsii.seeryear), hystavg=paste("BRFSS ", nhsii.brfssyear[1], "-", tail(nhsii.brfssyear,1), sep=""),
          mortality.rate=paste("CDC WONDER", nhsii.cdcyear), cancer.rate=paste("SEER", nhsii.seeryear)) %>%
  arrange(!is.na(age.group), age.group) %>%
  mutate(age.group=ifelse(!is.na(age.group), paste(age.group, "-", age.group+4, sep=""), NA)) %>%
  mutate(age.group=ifelse(is.na(age.group), "NHS II", age.group))

# Pulling baseline data used for PLCO
t.plco <- base.endo.all %>% 
  filter(brfss.year %in% plco.brfssyear & seer.year==plco.seeryear) %>%
  group_by(age.group) %>%
  mutate(hystavg=mean(hyst)) %>%
  ungroup() %>%
  filter(brfss.year==plco.brfssyear[1]) %>%
  mutate(rate.corrected=sprintf("%.1f", rate/(1-hystavg))) %>%
  mutate(hystavg=sprintf("%.3f", hystavg)) %>%
  mutate(rate=sprintf("%.1f", rate)) %>%
  filter(age.group>=start.age) %>%
  select(age.group, rate, hystavg, rate.corrected) %>%
  full_join(base.competing.mortality, by=c("age.group"="Age")) %>%
  filter(CDC.Year==plco.cdcyear) %>%
  mutate(mortality.rate=sprintf("%.1f", Mortality.Rate*100000)) %>%
  select(age.group, rate, hystavg, rate.corrected, mortality.rate) %>%
  full_join(base.competing.cancer, by=c("age.group"="Age")) %>%
  filter(SEER.Year==nhsi.seeryear) %>%
  mutate(cancer.rate=sprintf("%.1f", Cancer.Rate*100000)) %>%
  select(-Cancer.Rate, -SEER.Year) %>%
  filter(!is.na(rate)) %>%
  add_row(rate=paste("SEER", plco.seeryear), hystavg=paste("BRFSS ", plco.brfssyear[1], "-", tail(plco.brfssyear,1), sep=""),
          mortality.rate=paste("CDC WONDER", plco.cdcyear), cancer.rate=paste("SEER", plco.seeryear)) %>%
  arrange(!is.na(age.group), age.group) %>%
  mutate(age.group=ifelse(!is.na(age.group), paste(age.group, "-", age.group+4, sep=""), NA)) %>%
  mutate(age.group=ifelse(is.na(age.group), "PLCO", age.group))

# Combining data to create final table
st5 <- bind_rows(t.nhsi, t.nhsii, t.plco) %>% 
  select(age.group, hystavg, rate, rate.corrected, mortality.rate, cancer.rate)

write.csv(st5, "./Output/Table_S5.csv", row.names=F, na="")


# --- f. Table S6 ---

data1 <- nhanes.nhsi
data2 <- nhanes.nhsii

# Obtaining counts and %s for categorical variables
t1.binary <- data1 %>%
  summarize_at(vars(education1, education2, education3,
                    smoking1, smoking2, smoking3,
                    bmi.group1, bmi.group2, bmi.group3, bmi.group4, bmi.group5, 
                    # bmi.group6,
                    parity.group0, parity.group1, parity.group2, parity.group3, parity.group4, 
                    firstbirthage.group1, firstbirthage.group2, firstbirthage.group3, firstbirthage.group4, firstbirthage.group5, firstbirthage.group9,
                    menarchage2.group9, menarchage2.group10, menarchage2.group12, menarchage2.group14, menarchage2.group16,
                    anyhrtever,
                    ehrtever1, ehrtdur.group0, ehrtdur.group1, ehrtdur.group2, ehrtdur.group3,
                    ephrtever1, ephrtdur.group0, ephrtdur.group1, ephrtdur.group2, ephrtdur.group3,
                    ocever, ocdur.group0, ocdur.group1, ocdur.group2, ocdur.group3,
                    diabetes1, hypertension1),
               list(~wtd.mean(., weights=weights)), na.rm=T) %>%
  pivot_longer(cols=everything(), names_to="var", values_to="value") %>%
  mutate(summary1=sprintf("%.1f", value*100)) %>%
  select(var, summary1)

t2.binary <- data2 %>%
  summarize_at(vars(education1, education2, education3,
                    smoking1, smoking2, smoking3,
                    bmi.group1, bmi.group2, bmi.group3, bmi.group4, bmi.group5, 
                    # bmi.group6,
                    parity.group0, parity.group1, parity.group2, parity.group3, parity.group4, 
                    firstbirthage.group1, firstbirthage.group2, firstbirthage.group3, firstbirthage.group4, firstbirthage.group5, firstbirthage.group9,
                    menarchage2.group9, menarchage2.group10, menarchage2.group12, menarchage2.group14, menarchage2.group16,
                    anyhrtever,
                    ehrtever1, ehrtdur.group0, ehrtdur.group1, ehrtdur.group2, ehrtdur.group3,
                    ephrtever1, ephrtdur.group0, ephrtdur.group1, ephrtdur.group2, ephrtdur.group3,
                    ocever, ocdur.group0, ocdur.group1, ocdur.group2, ocdur.group3,
                    diabetes1, hypertension1),
               list(~wtd.mean(., weights=weights)), na.rm=T) %>%
  pivot_longer(cols=everything(), names_to="var", values_to="value") %>%
  mutate(summary2=sprintf("%.1f", value*100)) %>%
  select(var, summary2)

# Formating RR's obtained from group LASSO
t.rr <- endo.log.RR %>% data.frame() %>% rename(rr=".") %>% mutate(var=names(endo.log.RR))

# Combining Tables
t.var <- c("demographics",
           "education", "education1", "education2", "education3",
           "lifestyle",
           "smoking", "smoking1", "smoking2", "smoking3",
           "bmi", "bmi.group1", "bmi.group2", "bmi.group3", "bmi.group4", "bmi.group5", 
           "reproductive",
           "parity", "parity.group0", "parity.group1", "parity.group2", "parity.group3", "parity.group4",
           "firstbirthage.group", "firstbirthage.group1", "firstbirthage.group2", "firstbirthage.group3", "firstbirthage.group4", 
           "firstbirthage.group5", "firstbirthage.group9",
           "menarchage", "menarchage2.group9", "menarchage2.group10", "menarchage2.group12", 
           "menarchage2.group14", "menarchage2.group16",
           "anyhrtever",
           "ehrtever1", "ehrtdur", "ehrtdur.group0", "ehrtdur.group1", "ehrtdur.group2", "ehrtdur.group3",
           "ephrtever1", "ephrtdur", "ephrtdur.group0", "ephrtdur.group1", "ephrtdur.group2", "ephrtdur.group3",
           "ocever", "ocdur", "ocdur.group0", "ocdur.group1", "ocdur.group2", "ocdur.group3",
           "clinical",
           "diabetes1", "hypertension1",
           "oc.bmi", 
           "ocever.bmigroup3", "ocever.bmigroup4", "ocever.bmigroup5", 
           "ocdur1.bmigroup3", "ocdur1.bmigroup4", "ocdur1.bmigroup5", 
           "ocdur2.bmigroup3", "ocdur2.bmigroup4", "ocdur2.bmigroup5", 
           "ocdur3.bmigroup3", "ocdur3.bmigroup4", "ocdur3.bmigroup5", 
           "anyhrtever.bmi",
           "anyhrtever.bmigroup3", "anyhrtever.bmigroup4", "anyhrtever.bmigroup5", 
           "ehrt.bmi", 
           "ehrtever.bmigroup3", "ehrtever.bmigroup4", "ehrtever.bmigroup5", 
           "ehrtdur1.bmigroup3", "ehrtdur1.bmigroup4", 
           "ehrtdur2.bmigroup3", "ehrtdur2.bmigroup4", 
           "ehrtdur3.bmigroup3", "ehrtdur3.bmigroup4", 
           "ephrt.bmi",
           "ephrtever.bmigroup3", "ephrtever.bmigroup4", "ephrtever.bmigroup5",
           "ephrtdur1.bmigroup3", "ephrtdur1.bmigroup4", 
           "ephrtdur2.bmigroup3", "ephrtdur2.bmigroup4", 
           "ephrtdur3.bmigroup3", "ephrtdur3.bmigroup4"
) %>%
  data.frame() %>%
  mutate(var=as.character(.)) %>%
  mutate(order=row_number()) %>%
  select(var, order)

# Formatting final table
st6 <- t.var %>%
  left_join(bind_rows(t1.continuous, t1.binary), by="var") %>%
  left_join(bind_rows(t2.continuous, t2.binary), by="var") %>%
  left_join(t.rr, by="var") %>%
  mutate(rr=ifelse(!is.na(rr) & rr!="-", sprintf("%.2f", exp(as.numeric(as.character((rr))))), rr)) %>%
  mutate(rr=ifelse(grepl("[[:digit:]]", var)==T & is.na(rr), "(ref)", rr))

write.csv(st6, "./Output/Table_S6.csv", row.names=F, na="")
