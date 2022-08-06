################################################################################
#                                                                              #
# Endometrial Cancer Risk Prediction Model: Cleaning NHANES Data               #
#                                                                              #
# Author: Joy Shi                                                              #
# Last updated: 09/15/2021                                                     #
#                                                                              #
# Purpose of program: importing and cleaning BRFSS data to obtain prevalence   #
#                     of hysterectomy; combining with SEER data to obtain      #
#                     baseline rates of endometrial cancer                     #
#                                                                              #
################################################################################

# -------------------------------- (1) Set up ----------------------------------

setwd("D:/Dropbox/PhD Dissertation/Paper 3 - E2C2/Datasets")

if (!require('tidyverse')) install.packages('tidyverse'); library('tidyverse')
if (!require('foreign')) install.packages('foreign'); library('foreign')
if (!require('mice')) install.packages('mice'); library('mice')
if (!require('fastDummies')) install.packages('fastDummies'); library('fastDummies')

nhsi.nhanesyear <- "1999"
nhsii.nhanesyear <- "2007"
plco.nhanesyear <- "1999"
current.nhanesyear <- "2017"

start.age <- 45

years <- c(nhsi.nhanesyear, nhsii.nhanesyear, plco.nhanesyear, current.nhanesyear)


# -------------- (2) Importing, cleaning and imputing NHANES Data --------------

# -- a. 1999-2000 Data --
if (1999 %in% years){
  # Cleaning Data
  nhanes1999 <- lapply(Sys.glob("NHANES 1999-2000 Data/*.XPT"), read.xport) %>% 
    reduce(full_join, by="SEQN") %>%
    filter(RIAGENDR==2) %>%
    filter(RIDRETH1==3) %>%
    filter(RIDSTATR==2) %>%
    filter(RIDAGEYR>=start.age & RIDAGEYR<=75) %>%
    rename(age=RIDAGEYR) %>%
    mutate(education=ifelse(DMDEDUC2<=3, 1, ifelse(DMDEDUC2==4, 2, ifelse(DMDEDUC2==5, 3, NA)))) %>%
    mutate(education=factor(education)) %>%
    rename(weights=WTMEC2YR) %>%
    mutate(menarchage=ifelse(RHQ010!=999 & RHQ010!=0, RHQ010, NA)) %>%
    rename(parity=RHQ160) %>%
    mutate(parity=ifelse(!is.na(RHD130) & RHD130==2, 0, parity)) %>%
    mutate(parity=ifelse(parity==99|parity==77, NA, parity)) %>%
    rename(firstbirthage=RHQ180) %>%
    mutate(firstbirthage=ifelse(parity==0, 0, firstbirthage)) %>%
    mutate(firstbirthage=ifelse(firstbirthage==999, NA, firstbirthage)) %>%
    mutate(hyst=ifelse(RHD280!=9, 2-RHD280, NA)) %>%
    mutate(hyst=factor(hyst)) %>%
    mutate(ocever=ifelse(RHQ420!=9, 2-RHQ420, NA)) %>%
    mutate(ocever=factor(ocever)) %>%
    mutate(ocdur=ifelse(RHQ460Q!=99, RHQ460Q, NA)) %>%
    mutate(ocdur=ifelse(RHQ460U==1, ocdur/12, ocdur)) %>% 
    mutate(ocdur=ifelse(ocever==0, 0, ocdur)) %>%
    mutate(anyhrtever=ifelse(RHQ540!=9, 2-RHQ540, NA)) %>%
    mutate(anyhrtever=factor(anyhrtever)) %>%
    mutate(ehrtever=ifelse(RHQ554!=9, 2-RHQ554, NA)) %>%
    mutate(ehrtever=ifelse(anyhrtever==0, 0, ehrtever)) %>%
    mutate(ehrtever=factor(ehrtever)) %>%
    mutate(ehrtdur=ifelse(RHQ560Q!=99, RHQ560Q, NA)) %>%
    mutate(ehrtdur=ifelse(RHQ560U==1, ehrtdur/12, ehrtdur)) %>% 
    mutate(ehrtdur=ifelse(ehrtever==0, 0, ehrtdur)) %>%
    mutate(ephrtever=ifelse(RHQ570!=9, 2-RHQ570, NA)) %>%
    mutate(ephrtever=ifelse(anyhrtever==0, 0, ephrtever)) %>%
    mutate(ephrtever=factor(ephrtever)) %>%
    mutate(ephrtdur=ifelse(RHQ576Q!=99, RHQ576Q, NA)) %>%
    mutate(ephrtdur=ifelse(RHQ576U==1, ephrtdur/12, ephrtdur)) %>% 
    mutate(ephrtdur=ifelse(ephrtever==0, 0, ephrtdur)) %>%
    mutate(smoking=ifelse(SMQ020==1, 2, ifelse(SMQ020==2, 1, NA))) %>%
    mutate(smoking=ifelse(!is.na(SMQ040) & (SMQ040==1|SMQ040==2), 3, smoking)) %>%
    mutate(smoking=factor(smoking)) %>%
    rename(bmi=BMXBMI) %>%
    mutate(hypertension=ifelse(BPQ020!=9, 2-BPQ020, NA)) %>%
    mutate(hypertension=factor(hypertension)) %>%
    mutate(diabetes=ifelse(DIQ010==1, 1, ifelse(DIQ010==2|DIQ010==3, 0, NA))) %>%
    mutate(diabetes=factor(diabetes)) %>%
    select(weights, age, diabetes, hypertension, bmi, smoking, education, parity, ocever, ocdur, menarchage, 
           anyhrtever, ehrtever, ephrtever, ehrtdur, ephrtdur, firstbirthage, hyst) %>%
    mutate(row=row_number())
  
  # Imputing missing data
  nhanes1999.imputed <- nhanes1999 %>% 
    dplyr::select(-weights, -row) %>% 
    mice(., m=10, meth=c("", "logreg", "logreg", "pmm", "polyreg", "polyreg", "pmm", "logreg", "pmm", "pmm",
                         "logreg", "logreg", "logreg", "pmm", "pmm", "pmm", "logreg"), seed=1484) %>%
    complete(., action="long")
  
  # Adding weights
  nhanes1999.imputed$weights <- rep(nhanes1999$weights, 10)
}

# -- b. 2005-2006 Data --
if (2005 %in% years){
  # Cleaning Data
  nhanes2005 <- lapply(Sys.glob("NHANES 2005-2006 Data/*.XPT"), read.xport) %>% 
    reduce(full_join, by="SEQN") %>%
    filter(RIAGENDR==2) %>%
    filter(RIDRETH1==3) %>%
    filter(RIDSTATR==2) %>%
    filter(RIDAGEYR>=start.age & RIDAGEYR<=75) %>%
    rename(age=RIDAGEYR) %>%
    mutate(education=ifelse(DMDEDUC2<=3, 1, ifelse(DMDEDUC2==4, 2, ifelse(DMDEDUC2==5, 3, NA)))) %>%
    mutate(education=factor(education)) %>%
    rename(weights=WTMEC2YR) %>%
    mutate(menarchage=ifelse(RHQ010!=999 & RHQ010!=0, RHQ010, NA)) %>%
    rename(parity=RHQ160) %>%
    mutate(parity=ifelse(!is.na(RHQ131) & RHQ131==2, 0, parity)) %>%
    mutate(parity=ifelse(parity==99|parity==77, NA, parity)) %>%
    rename(firstbirthage=RHQ180) %>%
    mutate(firstbirthage=ifelse(parity==0, 0, firstbirthage)) %>%
    mutate(firstbirthage=ifelse(firstbirthage==999, NA, firstbirthage)) %>%
    mutate(hyst=ifelse(RHD280!=9, 2-RHD280, NA)) %>%
    mutate(hyst=factor(hyst)) %>%
    mutate(ocever=ifelse(RHQ420!=9, 2-RHQ420, NA)) %>%
    mutate(ocever=factor(ocever)) %>%
    mutate(ocdur=ifelse(RHQ460Q!=99, RHQ460Q, NA)) %>%
    mutate(ocdur=ifelse(RHQ460U==1, ocdur/12, ocdur)) %>% 
    mutate(ocdur=ifelse(ocever==0, 0, ocdur)) %>%
    mutate(anyhrtever=ifelse(RHQ540!=9, 2-RHQ540, NA)) %>%
    mutate(anyhrtever=factor(anyhrtever)) %>%
    mutate(ehrtever=ifelse(RHQ554!=9, 2-RHQ554, NA)) %>%
    mutate(ehrtever=ifelse(anyhrtever==0, 0, ehrtever)) %>%
    mutate(ehrtever=factor(ehrtever)) %>%
    mutate(ehrtdur=ifelse(RHQ560Q!=99, RHQ560Q, NA)) %>%
    mutate(ehrtdur=ifelse(RHQ560U==1, ehrtdur/12, ehrtdur)) %>% 
    mutate(ehrtdur=ifelse(ehrtever==0, 0, ehrtdur)) %>%
    mutate(ephrtever=ifelse(RHQ570!=9, 2-RHQ570, NA)) %>%
    mutate(ephrtever=ifelse(anyhrtever==0, 0, ephrtever)) %>%
    mutate(ephrtever=factor(ephrtever)) %>%
    mutate(ephrtdur=ifelse(RHQ576Q!=99, RHQ576Q, NA)) %>%
    mutate(ephrtdur=ifelse(RHQ576U==1, ephrtdur/12, ephrtdur)) %>% 
    mutate(ephrtdur=ifelse(ephrtever==0, 0, ephrtdur)) %>%
    mutate(smoking=ifelse(SMQ020==1, 2, ifelse(SMQ020==2, 1, NA))) %>%
    mutate(smoking=ifelse(!is.na(SMQ040) & (SMQ040==1|SMQ040==2), 3, smoking)) %>%
    mutate(smoking=factor(smoking)) %>%
    rename(bmi=BMXBMI) %>%
    mutate(hypertension=ifelse(BPQ020!=9, 2-BPQ020, NA)) %>%
    mutate(hypertension=factor(hypertension)) %>%
    mutate(diabetes=ifelse(DIQ010==1, 1, ifelse(DIQ010==2|DIQ010==3, 0, NA))) %>%
    mutate(diabetes=factor(diabetes)) %>%
    select(weights, age, diabetes, hypertension, bmi, smoking, education, parity, ocever, ocdur, menarchage, 
           anyhrtever, ehrtever, ehrtdur, ephrtever, ephrtdur, hyst, firstbirthage) %>%
    mutate(row=row_number())
  
  # Imputing missing data
  nhanes2005.imputed <- nhanes2005 %>% 
    dplyr::select(-weights, -row) %>% 
    mice(., m=10, meth=c("", "", "logreg", "pmm", "polyreg", "polyreg", "pmm", "logreg", "pmm", "pmm",
                         "logreg", "logreg", "pmm", "logreg", "pmm", "logreg", "pmm"), seed=1484) %>%
    complete(., action="long")
  
  nhanes2005.imputed$weights <- rep(nhanes2005$weights, 10)
}

# -- c. 2007-2008 Data --
if (2007 %in% years){
  # Cleaning data
  nhanes2007 <- lapply(Sys.glob("NHANES 2007-2008 Data/*.XPT"), read.xport) %>% 
    reduce(full_join, by="SEQN") %>%
    filter(RIAGENDR==2) %>%
    filter(RIDRETH1==3) %>%
    filter(RIDSTATR==2) %>%
    filter(RIDAGEYR>=start.age & RIDAGEYR<=75) %>%
    rename(age=RIDAGEYR) %>%
    mutate(education=ifelse(DMDEDUC2<=3, 1, ifelse(DMDEDUC2==4, 2, ifelse(DMDEDUC2==5, 3, NA)))) %>%
    mutate(education=factor(education)) %>%
    rename(weights=WTMEC2YR) %>%
    mutate(menarchage=ifelse(RHQ010!=999 & RHQ010!=0, RHQ010, NA)) %>%
    rename(parity=RHQ160) %>%
    mutate(parity=ifelse(!is.na(RHQ131) & RHQ131==2, 0, parity)) %>%
    rename(firstbirthage=RHD180) %>%
    mutate(firstbirthage=ifelse(parity==0, 0, firstbirthage)) %>%
    mutate(firstbirthage=ifelse(firstbirthage==999, NA, firstbirthage)) %>%
    mutate(hyst=ifelse(RHD280!=9, 2-RHD280, NA)) %>%
    mutate(hyst=factor(hyst)) %>%
    mutate(ocever=ifelse(RHQ420!=9, 2-RHQ420, NA)) %>%
    mutate(ocever=factor(ocever)) %>%
    mutate(ocdur=ifelse(RHQ460Q!=99, RHQ460Q, NA)) %>%
    mutate(ocdur=ifelse(RHQ460U==1, ocdur/12, ocdur)) %>% 
    mutate(ocdur=ifelse(ocever==0, 0, ocdur)) %>%
    mutate(anyhrtever=ifelse(RHQ540!=9, 2-RHQ540, NA)) %>%
    mutate(anyhrtever=factor(anyhrtever)) %>%
    mutate(ehrtever=ifelse(RHQ554!=9, 2-RHQ554, NA)) %>%
    mutate(ehrtever=ifelse(anyhrtever==0, 0, ehrtever)) %>%
    mutate(ehrtever=factor(ehrtever)) %>%
    mutate(ehrtdur=ifelse(RHQ560Q!=99, RHQ560Q, NA)) %>%
    mutate(ehrtdur=ifelse(RHQ560U==1, ehrtdur/12, ehrtdur)) %>% 
    mutate(ehrtdur=ifelse(ehrtever==0, 0, ehrtdur)) %>%
    mutate(ephrtever=ifelse(RHQ570!=9, 2-RHQ570, NA)) %>%
    mutate(ephrtever=ifelse(anyhrtever==0, 0, ephrtever)) %>%
    mutate(ephrtever=factor(ephrtever)) %>%
    mutate(ephrtdur=ifelse(RHQ576Q!=99, RHQ576Q, NA)) %>%
    mutate(ephrtdur=ifelse(RHQ576U==1, ephrtdur/12, ephrtdur)) %>% 
    mutate(ephrtdur=ifelse(ephrtever==0, 0, ephrtdur)) %>%
    mutate(smoking=ifelse(SMQ020==1, 2, ifelse(SMQ020==2, 1, NA))) %>%
    mutate(smoking=ifelse(!is.na(SMQ040) & (SMQ040==1|SMQ040==2), 3, smoking)) %>%
    mutate(smoking=factor(smoking)) %>%
    rename(bmi=BMXBMI) %>%
    mutate(hypertension=ifelse(BPQ020!=9, 2-BPQ020, NA)) %>%
    mutate(hypertension=factor(hypertension)) %>%
    mutate(diabetes=ifelse(DIQ010==1, 1, ifelse(DIQ010==2|DIQ010==3, 0, NA))) %>%
    mutate(diabetes=factor(diabetes)) %>%
    select(weights, age, hypertension, diabetes, bmi, education, smoking, ocever, ocdur, menarchage, hyst, parity, 
           anyhrtever, ehrtever, ephrtever, ephrtdur, ehrtdur, firstbirthage) %>%
    mutate(row=row_number())
  
  # Imputing missing data
  nhanes2007.imputed <- nhanes2007 %>% 
    dplyr::select(-weights, -row) %>% 
    mice(., m=10, meth=c("", "logreg", "logreg", "pmm", "polyreg", "polyreg", "polyreg", "pmm", "pmm", "polyreg", "pmm",
                         "polyreg", "polyreg", "polyreg", "pmm", "pmm", "pmm"), seed=10432) %>%
    complete(., action="long")
  
  nhanes2007.imputed$weights <- rep(nhanes2007$weights, 10)
}

# -- d. 2009-2010 Data --
if (2009 %in% years){
  # Cleaning data
  nhanes2009 <- lapply(Sys.glob("NHANES 2009-2010 Data/*.XPT"), read.xport) %>% 
    reduce(full_join, by="SEQN") %>%
    filter(RIAGENDR==2) %>%
    filter(RIDRETH1==3) %>%
    filter(RIDSTATR==2) %>%
    filter(RIDAGEYR>=start.age & RIDAGEYR<=75) %>%
    rename(age=RIDAGEYR) %>%
    mutate(education=ifelse(DMDEDUC2<=3, 1, ifelse(DMDEDUC2==4, 2, ifelse(DMDEDUC2==5, 3, NA)))) %>%
    mutate(education=factor(education)) %>%
    rename(weights=WTMEC2YR) %>%
    mutate(menarchage=ifelse(RHQ010!=999 & RHQ010!=0, RHQ010, NA)) %>%
    rename(parity=RHQ160) %>%
    mutate(parity=ifelse(!is.na(RHQ131) & RHQ131==2, 0, parity)) %>%
    rename(firstbirthage=RHD180) %>%
    mutate(firstbirthage=ifelse(parity==0, 0, firstbirthage)) %>%
    mutate(firstbirthage=ifelse(firstbirthage==999, NA, firstbirthage)) %>%
    mutate(hyst=ifelse(RHD280!=9, 2-RHD280, NA)) %>%
    mutate(hyst=factor(hyst)) %>%
    mutate(ocever=ifelse(RHQ420!=9, 2-RHQ420, NA)) %>%
    mutate(ocever=factor(ocever)) %>%
    mutate(ocdur=ifelse(RHQ460Q!=99, RHQ460Q, NA)) %>%
    mutate(ocdur=ifelse(RHQ460U==1, ocdur/12, ocdur)) %>% 
    mutate(ocdur=ifelse(ocever==0, 0, ocdur)) %>%
    mutate(anyhrtever=ifelse(RHQ540!=9, 2-RHQ540, NA)) %>%
    mutate(anyhrtever=factor(anyhrtever)) %>%
    mutate(ehrtever=ifelse(RHQ554!=9, 2-RHQ554, NA)) %>%
    mutate(ehrtever=ifelse(anyhrtever==0, 0, ehrtever)) %>%
    mutate(ehrtever=factor(ehrtever)) %>%
    mutate(ehrtdur=ifelse(RHQ560Q!=99, RHQ560Q, NA)) %>%
    mutate(ehrtdur=ifelse(RHQ560U==1, ehrtdur/12, ehrtdur)) %>% 
    mutate(ehrtdur=ifelse(ehrtever==0, 0, ehrtdur)) %>%
    mutate(ephrtever=ifelse(RHQ570!=9, 2-RHQ570, NA)) %>%
    mutate(ephrtever=ifelse(anyhrtever==0, 0, ephrtever)) %>%
    mutate(ephrtever=factor(ephrtever)) %>%
    mutate(ephrtdur=ifelse(RHQ576Q!=99, RHQ576Q, NA)) %>%
    mutate(ephrtdur=ifelse(RHQ576U==1, ephrtdur/12, ephrtdur)) %>% 
    mutate(ephrtdur=ifelse(ephrtever==0, 0, ephrtdur)) %>%
    mutate(smoking=ifelse(SMQ020==1, 2, ifelse(SMQ020==2, 1, NA))) %>%
    mutate(smoking=ifelse(!is.na(SMQ040) & (SMQ040==1|SMQ040==2), 3, smoking)) %>%
    mutate(smoking=factor(smoking)) %>%
    rename(bmi=BMXBMI) %>%
    mutate(hypertension=ifelse(BPQ020!=9, 2-BPQ020, NA)) %>%
    mutate(hypertension=factor(hypertension)) %>%
    mutate(diabetes=ifelse(DIQ010==1, 1, ifelse(DIQ010==2|DIQ010==3, 0, NA))) %>%
    mutate(diabetes=factor(diabetes)) %>%
    select(weights, age, hypertension, diabetes, bmi, smoking, education, ocever, ocdur, menarchage, parity, hyst,
           anyhrtever, ehrtever, ephrtever, ephrtdur, ehrtdur, firstbirthage) %>%
    mutate(row=row_number())
  
  # Imputing missing data
  nhanes2009.imputed <- nhanes2009 %>% 
    dplyr::select(-weights, -row) %>% 
    mice(., m=10, meth=c("", "logreg", "logreg", "pmm", "polyreg", "polyreg", "logreg", "pmm", "pmm", "pmm", "logreg",
                         "logreg", "logreg", "logreg", "pmm", "pmm", "pmm"), seed=10432) %>%
    complete(., action="long")
  
  # Adding weights
  nhanes2009.imputed$weights <- rep(nhanes2009$weights, 10)
}

# -- e. 2017-2018 Data --
if (2017 %in% years){
  # Cleaning data for 2011 first (needed for imputation)
  nhanes2011 <- lapply(Sys.glob("NHANES 2011-2012 Data/*.XPT"), read.xport) %>% 
    reduce(full_join, by="SEQN") %>%
    filter(RIAGENDR==2) %>%
    filter(RIDRETH1==3) %>%
    filter(RIDSTATR==2) %>%
    filter(RIDAGEYR>=start.age & RIDAGEYR<=75) %>%
    rename(age=RIDAGEYR) %>%
    mutate(education=ifelse(DMDEDUC2<=3, 1, ifelse(DMDEDUC2==4, 2, ifelse(DMDEDUC2==5, 3, NA)))) %>%
    mutate(education=factor(education)) %>%
    rename(weights=WTMEC2YR) %>%
    mutate(menarchage=ifelse(RHQ010!=999 & RHQ010!=0, RHQ010, NA)) %>%
    rename(parity=RHQ160) %>%
    mutate(parity=ifelse(!is.na(RHQ131) & RHQ131==2, 0, parity)) %>%
    rename(firstbirthage=RHD180) %>%
    mutate(firstbirthage=ifelse(parity==0, 0, firstbirthage)) %>%
    mutate(firstbirthage=ifelse(firstbirthage==999, NA, firstbirthage)) %>%
    mutate(hyst=ifelse(RHD280!=9, 2-RHD280, NA)) %>%
    mutate(hyst=factor(hyst)) %>%
    mutate(ocever=ifelse(RHQ420!=9, 2-RHQ420, NA)) %>%
    mutate(ocever=factor(ocever)) %>%
    mutate(ocdur=ifelse(RHQ460Q!=99, RHQ460Q, NA)) %>%
    mutate(ocdur=ifelse(RHQ460U==1, ocdur/12, ocdur)) %>%
    mutate(ocdur=ifelse(ocever==0, 0, ocdur)) %>%
    mutate(anyhrtever=ifelse(RHQ540!=9, 2-RHQ540, NA)) %>%
    mutate(anyhrtever=factor(anyhrtever)) %>%
    mutate(ehrtever=ifelse(RHQ554!=9, 2-RHQ554, NA)) %>%
    mutate(ehrtever=ifelse(anyhrtever==0, 0, ehrtever)) %>%
    mutate(ehrtever=factor(ehrtever)) %>%
    mutate(ehrtdur=ifelse(RHQ560Q!=99, RHQ560Q, NA)) %>%
    mutate(ehrtdur=ifelse(RHQ560U==1, ehrtdur/12, ehrtdur)) %>% 
    mutate(ehrtdur=ifelse(ehrtever==0, 0, ehrtdur)) %>%
    mutate(ephrtever=ifelse(RHQ570!=9, 2-RHQ570, NA)) %>%
    mutate(ephrtever=ifelse(anyhrtever==0, 0, ephrtever)) %>%
    mutate(ephrtever=factor(ephrtever)) %>%
    mutate(ephrtdur=ifelse(RHQ576Q!=99, RHQ576Q, NA)) %>%
    mutate(ephrtdur=ifelse(RHQ576U==1, ephrtdur/12, ephrtdur)) %>% 
    mutate(ephrtdur=ifelse(ephrtever==0, 0, ephrtdur)) %>%
    mutate(smoking=ifelse(SMQ020==1, 2, ifelse(SMQ020==2, 1, NA))) %>%
    mutate(smoking=ifelse(!is.na(SMQ040) & (SMQ040==1|SMQ040==2), 3, smoking)) %>%
    mutate(smoking=factor(smoking)) %>%
    rename(bmi=BMXBMI) %>%
    mutate(hypertension=ifelse(BPQ020!=9, 2-BPQ020, NA)) %>%
    mutate(hypertension=factor(hypertension)) %>%
    mutate(diabetes=ifelse(DIQ010==1, 1, ifelse(DIQ010==2|DIQ010==3, 0, NA))) %>%
    mutate(diabetes=factor(diabetes)) %>%
    mutate(year=2011) %>%
    select(year, weights, age, hypertension, diabetes, bmi, smoking, education, ocever, ocdur, menarchage, parity, hyst,
           anyhrtever, ehrtever, ephrtever, ephrtdur, ehrtdur, firstbirthage) %>%
    mutate(row=row_number())
  
  # Cleaning data for 2017-2018
  nhanes2017 <- lapply(Sys.glob("NHANES 2017-2018 Data/*.XPT"), read.xport) %>% 
    reduce(full_join, by="SEQN") %>%
    filter(RIAGENDR==2) %>%
    filter(RIDRETH1==3) %>%
    filter(RIDSTATR==2) %>%
    filter(RIDAGEYR>=start.age & RIDAGEYR<=75) %>%
    rename(age=RIDAGEYR) %>%
    mutate(education=ifelse(DMDEDUC2<=3, 1, ifelse(DMDEDUC2==4, 2, ifelse(DMDEDUC2==5, 3, NA)))) %>%
    mutate(education=factor(education)) %>%
    rename(weights=WTMEC2YR) %>%
    mutate(menarchage=ifelse(RHQ010!=999 & RHQ010!=0, RHQ010, NA)) %>%
    rename(parity=RHQ160) %>%
    mutate(parity=ifelse(!is.na(RHQ131) & RHQ131==2, 0, parity)) %>%
    rename(firstbirthage=RHD180) %>%
    mutate(firstbirthage=ifelse(parity==0, 0, firstbirthage)) %>%
    mutate(firstbirthage=ifelse(firstbirthage==999, NA, firstbirthage)) %>%
    mutate(hyst=ifelse(RHD280!=9, 2-RHD280, NA)) %>%
    mutate(hyst=factor(hyst)) %>%
    mutate(ocever=ifelse(RHQ420!=9, 2-RHQ420, NA)) %>%
    mutate(ocever=factor(ocever)) %>%
    # mutate(ocdur=ifelse(RHQ460Q!=99, RHQ460Q, NA)) %>%
    # mutate(ocdur=ifelse(RHQ460U==1, ocdur/12, ocdur)) %>% 
    mutate(ocdur=ifelse(ocever==0, 0, NA)) %>%
    mutate(anyhrtever=ifelse(RHQ540!=9, 2-RHQ540, NA)) %>%
    mutate(anyhrtever=factor(anyhrtever)) %>%
    mutate(ehrtever=ifelse(RHQ554!=9, 2-RHQ554, NA)) %>%
    mutate(ehrtever=ifelse(anyhrtever==0, 0, ehrtever)) %>%
    mutate(ehrtever=factor(ehrtever)) %>%
    mutate(ehrtdur=ifelse(RHQ560Q!=99, RHQ560Q, NA)) %>%
    mutate(ehrtdur=ifelse(RHQ560U==1, ehrtdur/12, ehrtdur)) %>% 
    mutate(ehrtdur=ifelse(ehrtever==0, 0, ehrtdur)) %>%
    mutate(ephrtever=ifelse(RHQ570!=9, 2-RHQ570, NA)) %>%
    mutate(ephrtever=ifelse(anyhrtever==0, 0, ephrtever)) %>%
    mutate(ephrtever=factor(ephrtever)) %>%
    mutate(ephrtdur=ifelse(RHQ576Q!=99, RHQ576Q, NA)) %>%
    mutate(ephrtdur=ifelse(RHQ576U==1, ephrtdur/12, ephrtdur)) %>% 
    mutate(ephrtdur=ifelse(ephrtever==0, 0, ephrtdur)) %>%
    mutate(smoking=ifelse(SMQ020==1, 2, ifelse(SMQ020==2, 1, NA))) %>%
    mutate(smoking=ifelse(!is.na(SMQ040) & (SMQ040==1|SMQ040==2), 3, smoking)) %>%
    mutate(smoking=factor(smoking)) %>%
    rename(bmi=BMXBMI) %>%
    mutate(hypertension=ifelse(BPQ020!=9, 2-BPQ020, NA)) %>%
    mutate(hypertension=factor(hypertension)) %>%
    mutate(diabetes=ifelse(DIQ010==1, 1, ifelse(DIQ010==2|DIQ010==3, 0, NA))) %>%
    mutate(diabetes=factor(diabetes)) %>%
    mutate(year=2017) %>%
    select(year, weights, age, hypertension, diabetes, bmi, smoking, education, ocever, ocdur, menarchage, parity, hyst,
           anyhrtever, ehrtever, ephrtever, ephrtdur, ehrtdur, firstbirthage) %>%
    mutate(row=row_number())
  
  # Imputing missing data
  nhanes2017.imputed <- nhanes2017 %>%
    # need firstbirthage data from 2011 cycle for imputation
    bind_rows(nhanes2011) %>%
    dplyr::select(-weights, -row) %>%
    select(year, age, hypertension, diabetes, smoking, education, bmi, hyst, ocever, parity, anyhrtever, menarchage, ehrtever, ephrtever, ephrtdur, ehrtdur, 
           firstbirthage, ocdur) %>%
    mice(., m=10, meth=c("", "", "", "",  "polyreg", "polyreg", "pmm", "logreg", "logreg", "pmm", "logreg", "pmm", "logreg", 
                         "logreg", "pmm", "pmm", "pmm", "pmm"), seed=10432) %>%
    complete(., action="long") %>%
    filter(year==2017)
  
  nhanes2017.imputed$weights <- rep(nhanes2017$weights, 10)
}


# ---------------------- (3) Restructuring NHANES Data -------------------------

for (i in ls(pattern=".imputed")){
  temp <- get(i) %>% 
    mutate(education1=ifelse(education==1, 1, 0)) %>%
    mutate(education2=ifelse(education==2, 1, 0)) %>%
    mutate(education3=ifelse(education==3, 1, 0)) %>%
    mutate(smoking1=ifelse(smoking==1, 1, 0)) %>%
    mutate(smoking2=ifelse(smoking==2, 1, 0)) %>%
    mutate(smoking3=ifelse(smoking==3, 1, 0)) %>%
    mutate(bmi.group=cut(bmi, breaks=c(-Inf, 18.5, 25, 30, 35,Inf), right=F, 
      labels=c(1, 2, 3, 4, 5))) %>%
    mutate(parity.group=ifelse(parity>=4, 4, parity)) %>%
    mutate(menarchage.group=ifelse(menarchage==0, NA, ifelse(menarchage<=10, 10, 
      ifelse(menarchage>=16, 16, menarchage)))) %>%  
    mutate(menarchage2.group=ifelse(menarchage==0, NA, ifelse(menarchage<=9, 9,
      ifelse(menarchage %in% c(10,11), 10, ifelse(menarchage %in% c(12,13), 12,
      ifelse(menarchage %in% c(14,15), 14, ifelse(menarchage >=16, 16, NA))))))) %>%                             
    mutate(ephrtdur.group=ifelse(ephrtdur>0 & ephrtdur<=5, 1, ifelse(ephrtdur>5 & ephrtdur<=10,
      2, ifelse(ephrtdur!=999 & ephrtdur!=0, 3, ephrtdur)))) %>%
    mutate(ehrtdur.group=ifelse(ehrtdur>0 & ehrtdur<=5, 1, ifelse(ehrtdur>5 & 
      ehrtdur<=10, 2, ifelse(ehrtdur!=999 & ehrtdur!=0, 3, ehrtdur)))) %>%
    mutate(firstbirthage.int=ifelse(firstbirthage==999, 0, firstbirthage)) %>%
    mutate(firstbirthage.sqint=(firstbirthage.int^2)) %>%
    mutate(firstbirthage.missing2=ifelse(firstbirthage==0, 1, 0)) %>%
    mutate(firstbirthage.group=cut(firstbirthage, breaks=c(-Inf, 20, 25, 30, 35, Inf), 
      labels=seq(1:5), right=F)) %>%
    mutate(firstbirthage.group=ifelse(!is.na(parity.group) & parity.group==0, 9, 
      firstbirthage.group)) %>%  
    mutate(ocdur.group=ifelse(ocdur>0 & ocdur<=5, 1, ifelse(ocdur>5 & ocdur<=10, 
      2, ifelse(ocdur!=999 & ocdur!=0, 3, ocdur)))) %>%
    # Interaction between ever OC use and BMI
    mutate(ocever.bmigroup1=ifelse(bmi.group==1, ocever, 0)) %>%
    mutate(ocever.bmigroup2=ifelse(bmi.group==2, ocever, 0)) %>%
    mutate(ocever.bmigroup3=ifelse(bmi.group==3, ocever, 0)) %>%
    mutate(ocever.bmigroup4=ifelse(bmi.group==4, ocever, 0)) %>%
    mutate(ocever.bmigroup5=ifelse(bmi.group==5, ocever, 0)) %>%
    mutate(ocever.bmigroup6=ifelse(bmi.group==6, ocever, 0)) %>%  
    # Interaction between duration of OC use and BMI
    mutate(ocdur1.bmigroup1=ifelse(ocdur.group==1 & bmi.group==1, 1, 0)) %>%
    mutate(ocdur1.bmigroup2=ifelse(ocdur.group==1 & bmi.group==2, 1, 0)) %>%
    mutate(ocdur1.bmigroup3=ifelse(ocdur.group==1 & bmi.group==3, 1, 0)) %>%
    mutate(ocdur1.bmigroup4=ifelse(ocdur.group==1 & bmi.group==4, 1, 0)) %>%
    mutate(ocdur1.bmigroup5=ifelse(ocdur.group==1 & bmi.group==5, 1, 0)) %>%
    mutate(ocdur1.bmigroup6=ifelse(ocdur.group==1 & bmi.group==6, 1, 0)) %>%
    mutate(ocdur2.bmigroup1=ifelse(ocdur.group==2 & bmi.group==1, 1, 0)) %>%
    mutate(ocdur2.bmigroup2=ifelse(ocdur.group==2 & bmi.group==2, 1, 0)) %>%
    mutate(ocdur2.bmigroup3=ifelse(ocdur.group==2 & bmi.group==3, 1, 0)) %>%
    mutate(ocdur2.bmigroup4=ifelse(ocdur.group==2 & bmi.group==4, 1, 0)) %>%
    mutate(ocdur2.bmigroup5=ifelse(ocdur.group==2 & bmi.group==5, 1, 0)) %>%
    mutate(ocdur2.bmigroup6=ifelse(ocdur.group==2 & bmi.group==6, 1, 0)) %>%
    mutate(ocdur3.bmigroup1=ifelse(ocdur.group==3 & bmi.group==1, 1, 0)) %>%
    mutate(ocdur3.bmigroup2=ifelse(ocdur.group==3 & bmi.group==2, 1, 0)) %>%
    mutate(ocdur3.bmigroup3=ifelse(ocdur.group==3 & bmi.group==3, 1, 0)) %>%
    mutate(ocdur3.bmigroup4=ifelse(ocdur.group==3 & bmi.group==4, 1, 0)) %>%
    mutate(ocdur3.bmigroup5=ifelse(ocdur.group==3 & bmi.group==5, 1, 0)) %>%
    mutate(ocdur3.bmigroup6=ifelse(ocdur.group==3 & bmi.group==6, 1, 0)) %>%  
    # Interaction between any HT use and BMI 
    mutate(anyhrtever.bmigroup1=ifelse(bmi.group==1, anyhrtever, 0)) %>%
    mutate(anyhrtever.bmigroup2=ifelse(bmi.group==2, anyhrtever, 0)) %>%
    mutate(anyhrtever.bmigroup3=ifelse(bmi.group==3, anyhrtever, 0)) %>%
    mutate(anyhrtever.bmigroup4=ifelse(bmi.group==4, anyhrtever, 0)) %>%
    mutate(anyhrtever.bmigroup5=ifelse(bmi.group==5, anyhrtever, 0)) %>%
    mutate(anyhrtever.bmigroup6=ifelse(bmi.group==6, anyhrtever, 0)) %>%
    # Interaction between e-only HT use and BMI
    mutate(ehrtever.bmigroup1=ifelse(bmi.group==1 & ehrtever==1, 1, 0)) %>%
    mutate(ehrtever.bmigroup2=ifelse(bmi.group==2 & ehrtever==1, 1, 0)) %>%
    mutate(ehrtever.bmigroup3=ifelse(bmi.group==3 & ehrtever==1, 1, 0)) %>%
    mutate(ehrtever.bmigroup4=ifelse(bmi.group==4 & ehrtever==1, 1, 0)) %>%
    mutate(ehrtever.bmigroup5=ifelse(bmi.group==5 & ehrtever==1, 1, 0)) %>%
    mutate(ehrtever.bmigroup6=ifelse(bmi.group==6 & ehrtever==1, 1, 0)) %>%
    # Interaction between EP HT use and BMI
    mutate(ephrtever.bmigroup1=ifelse(bmi.group==1 & ephrtever==1, 1, 0)) %>%
    mutate(ephrtever.bmigroup2=ifelse(bmi.group==2 & ephrtever==1, 1, 0)) %>%
    mutate(ephrtever.bmigroup3=ifelse(bmi.group==3 & ephrtever==1, 1, 0)) %>%
    mutate(ephrtever.bmigroup4=ifelse(bmi.group==4 & ephrtever==1, 1, 0)) %>%
    mutate(ephrtever.bmigroup5=ifelse(bmi.group==5 & ephrtever==1, 1, 0)) %>%
    mutate(ephrtever.bmigroup6=ifelse(bmi.group==6 & ephrtever==1, 1, 0)) %>%
    # Interaction between duration of e-only HT use and BMI
    mutate(ehrtdur1.bmigroup1=ifelse(ehrtdur.group==1 & bmi.group==1, 1, 0)) %>%
    mutate(ehrtdur1.bmigroup2=ifelse(ehrtdur.group==1 & bmi.group==2, 1, 0)) %>%
    mutate(ehrtdur1.bmigroup3=ifelse(ehrtdur.group==1 & bmi.group==3, 1, 0)) %>%
    mutate(ehrtdur1.bmigroup4=ifelse(ehrtdur.group==1 & (bmi.group==4|bmi.group==5), 1, 0)) %>%
    mutate(ehrtdur2.bmigroup1=ifelse(ehrtdur.group==2 & bmi.group==1, 1, 0)) %>%
    mutate(ehrtdur2.bmigroup2=ifelse(ehrtdur.group==2 & bmi.group==2, 1, 0)) %>%
    mutate(ehrtdur2.bmigroup3=ifelse(ehrtdur.group==2 & bmi.group==3, 1, 0)) %>%
    mutate(ehrtdur2.bmigroup4=ifelse(ehrtdur.group==2 & (bmi.group==4|bmi.group==5), 1, 0)) %>%
    mutate(ehrtdur3.bmigroup1=ifelse(ehrtdur.group==3 & bmi.group==1, 1, 0)) %>%
    mutate(ehrtdur3.bmigroup2=ifelse(ehrtdur.group==3 & bmi.group==2, 1, 0)) %>%
    mutate(ehrtdur3.bmigroup3=ifelse(ehrtdur.group==3 & bmi.group==3, 1, 0)) %>%
    mutate(ehrtdur3.bmigroup4=ifelse(ehrtdur.group==3 & (bmi.group==4|bmi.group==5), 1, 0)) %>%
    # Interaction between duration of EP HT use and BMI
    mutate(ephrtdur1.bmigroup1=ifelse(ephrtdur.group==1 & bmi.group==1, 1, 0)) %>%
    mutate(ephrtdur1.bmigroup2=ifelse(ephrtdur.group==1 & bmi.group==2, 1, 0)) %>%
    mutate(ephrtdur1.bmigroup3=ifelse(ephrtdur.group==1 & bmi.group==3, 1, 0)) %>%
    mutate(ephrtdur1.bmigroup4=ifelse(ephrtdur.group==1 & (bmi.group==4|bmi.group==5), 1, 0)) %>%
    mutate(ephrtdur2.bmigroup1=ifelse(ephrtdur.group==2 & bmi.group==1, 1, 0)) %>%
    mutate(ephrtdur2.bmigroup2=ifelse(ephrtdur.group==2 & bmi.group==2, 1, 0)) %>%
    mutate(ephrtdur2.bmigroup3=ifelse(ephrtdur.group==2 & bmi.group==3, 1, 0)) %>%
    mutate(ephrtdur2.bmigroup4=ifelse(ephrtdur.group==2 & (bmi.group==4|bmi.group==5), 1, 0)) %>%
    mutate(ephrtdur3.bmigroup1=ifelse(ephrtdur.group==3 & bmi.group==1, 1, 0)) %>%
    mutate(ephrtdur3.bmigroup2=ifelse(ephrtdur.group==3 & bmi.group==2, 1, 0)) %>%
    mutate(ephrtdur3.bmigroup3=ifelse(ephrtdur.group==3 & bmi.group==3, 1, 0)) %>%
    mutate(ephrtdur3.bmigroup4=ifelse(ephrtdur.group==3 & (bmi.group==4|bmi.group==5), 1, 0)) %>%
    rename(diabetes1=diabetes) %>%
    rename(hypertension1=hypertension) %>%
    rename(ehrtever1=ehrtever) %>%
    rename(ephrtever1=ephrtever) %>%
    mutate(ehrtdur.bmigroup=(as.numeric(as.character(bmi.group))-1)*ehrtdur.group) %>%
    dummy_cols(., select_columns=c("bmi.group", "parity.group", "menarchage.group", 
      "menarchage2.group", "ehrtdur.group", "ephrtdur.group", "ocdur.group", "firstbirthage.group"), ignore_na=T)  
  
  names(temp) <- sub("_", "", names(temp))
  assign(paste0(i, ".ref", sep=""), temp)
}


# --------------- (3) Assigning NHANES Data to Validation Cohorts --------------

nhsi.ref.dataset <- get(paste("nhanes", nhsi.nhanesyear, ".imputed.ref", sep=""))
nhsii.ref.dataset <- get(paste("nhanes", nhsii.nhanesyear, ".imputed.ref", sep=""))
plco.ref.dataset <- get(paste("nhanes", plco.nhanesyear, ".imputed.ref", sep=""))
current.ref.dataset <- get(paste("nhanes", current.nhanesyear, ".imputed.ref", sep=""))


# ---------------------------- (4) Exporting data ------------------------------

write.csv(nhsi.ref.dataset, "../Analysis/Output/E2C2_Development_5_NHANES_NHSI.csv", row.names=F)
write.csv(nhsii.ref.dataset, "../Analysis/Output/E2C2_Development_5_NHANES_NHSII.csv", row.names=F)
write.csv(plco.ref.dataset, "../Analysis/Output/E2C2_Development_5_NHANES_PLCO.csv", row.names=F)
write.csv(current.ref.dataset, "../Analysis/Output/E2C2_Development_5_NHANES_Current.csv", row.names=F)

