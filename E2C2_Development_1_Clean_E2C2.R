################################################################################
#                                                                              #
# Endometrial Cancer Risk Prediction Model: Cleaning E2C2 Data                 #
#                                                                              #
# Author: Joy Shi                                                              #
# Last updated: 11/06/2021                                                     #
#                                                                              #
# Purpose of program: importing and cleaning E2C2 data                         #
#                                                                              #
################################################################################


# -------------------------------- (1) Set up ----------------------------------

setwd("D:/Dropbox/PhD Dissertation/Paper 3 - E2C2/Datasets")

if (!require('tidyverse')) install.packages('tidyverse'); library('tidyverse')
if (!require('haven')) install.packages('haven'); library('haven')
if (!require('fastDummies')) install.packages('fastDummies'); library('fastDummies')

start.age <- 45


# --------------- (2) Importing and cleaning questionnaire data ----------------

# -- a. Importing and binding data --
data_import <- lapply(Sys.glob("core*_harmonized.sas7bdat"), read_sas)
data_bind <- bind_rows(data_import)

# -- b. Cleaning variables --
data_clean <- data_bind %>%
  filter(site!=107) %>%
  filter(site!=208) %>%
  filter(site!=211) %>%
  mutate(dobyear=ifelse(dobyear==8888, NA, dobyear)) %>%
  mutate(age=ifelse(dxage==888, NA, dxage))%>%
  mutate(dxyear=ifelse(dxyear==8888, NA, dxyear)) %>%
  mutate(intyear=ifelse(intyear==8888, NA, intyear)) %>%
  mutate(race=ifelse(race==8, NA, race)) %>%
  mutate(hispanic=ifelse(hispanic==8, NA, 2-hispanic)) %>%
  mutate(education=ifelse(education==888, NA, education)) %>%
  mutate(pathrev=ifelse(pathrev==8, NA, pathrev)) %>%
  mutate(grade=ifelse(grade==8, NA, grade)) %>%
  mutate(height=ifelse(height==888, NA, height)) %>%
  mutate(weight=ifelse(wgtref==888, NA, wgtref)) %>%
  mutate(bmi=weight/(height/100)^2) %>%
  mutate(bmi=ifelse(bmi>quantile(bmi, 0.999, na.rm=T), quantile(bmi, 0.999, na.rm=T), bmi)) %>%
  mutate(bmi=ifelse(bmi<quantile(bmi, 0.001, na.rm=T), quantile(bmi, 0.001, na.rm=T), bmi)) %>%
  mutate(smoking=ifelse(smoking==8, NA, smoking)) %>%
  mutate(packyears=ifelse(packyears==888, NA, packyears)) %>%
  mutate(alcohol=ifelse(alcohol==888, NA, alcohol)) %>%
  mutate(alcoholgr=ifelse(alcograms==8888, NA, alcograms)) %>%
  mutate(complpregs=ifelse(complpregs==888, NA, complpregs)) %>%
  mutate(pregs=ifelse(pregs==888, NA, pregs)) %>%
  mutate(parity=ifelse(parity==888, NA, parity)) %>%
  mutate(lastpregage=ifelse(lastpregage==888, NA, lastpregage)) %>%
  mutate(firstbirthage=ifelse(firstbirthage %in% c(888, 88, 0), NA, firstbirthage)) %>%
  mutate(lastbirthage=ifelse(lastbirthage==888, NA, lastbirthage)) %>%
  mutate(menarchage=ifelse(menarchage==88, NA, menarchage)) %>%
  mutate(menarchage=ifelse(menarchage>21, 21, menarchage)) %>%
  mutate(menstat=ifelse(menstat==8, NA, menstat)) %>%
  mutate(anyhrtever=ifelse(any_hrt_ever==8, NA, 2-any_hrt_ever)) %>%
  mutate(ephrtever=ifelse(combo==8, NA, 2-combo)) %>%
  mutate(ephrtdur=ifelse(combomos==888, NA, combomos/12)) %>%
  mutate(ehrtever=ifelse(estrogen==8, NA, 2-estrogen)) %>%
  mutate(ehrtdur=ifelse(estrogenmos==888, NA, estrogenmos/12)) %>%
  mutate(ocever=ifelse(ocever==8, NA, 2-ocever)) %>%
  mutate(ocdur=ifelse(ocmos==888, NA, ocmos/12)) %>%
  mutate(diabetes=ifelse(diabetes==8, NA, 2-diabetes)) %>%
  mutate(endom=ifelse(endom==8, NA, 2-endom)) %>%
  mutate(hypertension=ifelse(hypertension==8, NA, 2-hypertension)) %>%
  mutate(famhx=ifelse(endocancersis==1|endocancermom==1|endocancerdau==1, 1, 0)) %>%
  mutate(famhx=ifelse(is.na(famhx) & (!is.na(endocancersis)|!is.na(endocancermom)|!is.na(endocancerdau)), 0, famhx)) %>%
  #mutate(famhx=ifelse(is.na(famhx), 0, famhx)) %>%
  mutate(casecon=2-casecon) %>%
  mutate(age=ifelse(is.na(age), intyear-dobyear, age)) %>%
  select(-any_hrt_ever)

# -- c. Identifying exclusions --
exclusions <- data_clean %>%
  mutate(all=1) %>%
  mutate(age.exclude=ifelse(is.na(age)|age<start.age|age>85, 1, 0)) %>%
  mutate(menstat.exclude=ifelse(menstat!=3, 1, 0)) %>%
  mutate(menstat.exclude=ifelse(is.na(menstat.exclude), 0, menstat.exclude)) %>%
  mutate(menstat.exclude=ifelse(age.exclude==1, 0, menstat.exclude)) %>%
  mutate(race.exclude=ifelse(race!=1, 1, 0)) %>%
  mutate(race.exclude=ifelse(is.na(race.exclude), 0, race.exclude)) %>%
  mutate(race.exclude=ifelse(menstat.exclude==1|age.exclude==1, 0, race.exclude)) %>%
  mutate(exclude=ifelse(race.exclude==1|menstat.exclude==1|age.exclude==1, 1, 0)) %>%
  mutate(include=ifelse(race.exclude==0 & menstat.exclude==0 & age.exclude==0, 1, 0)) %>%
  group_by(casecon) %>%
  summarize_at(vars(all, exclude, age.exclude, menstat.exclude, race.exclude, include), list(~sum(., na.rm=T)))  

# -- d. Removing excluded IDs --
train_clinical <- data_clean %>% filter(race==1|is.na(race)) %>% filter(menstat==3|is.na(menstat)) %>% filter(!is.na(age) & age>=start.age & age<=85) %>%
  mutate(ocdur=ifelse(is.na(ocdur) & site %in% unique(data_clean[which(!is.na(data_clean$ocdur)),]$site) & ocever==0, 0, ocdur)) %>%
  mutate(ehrtdur=ifelse(is.na(ehrtdur) & site %in% unique(data_clean[which(!is.na(data_clean$ehrtdur)),]$site) & ehrtever==0, 0, ehrtdur)) %>%
  mutate(ehrtever=ifelse(is.na(ehrtever) & site %in% unique(data_clean[which(!is.na(data_clean$ehrtever)),]$site) & ehrtdur==0, 0, ehrtever)) %>%
  mutate(ephrtdur=ifelse(is.na(ephrtdur) & site %in% unique(data_clean[which(!is.na(data_clean$ephrtdur)),]$site) & ephrtever==0, 0, ephrtdur)) %>%
  mutate(ephrtever=ifelse(is.na(ephrtever) & site %in% unique(data_clean[which(!is.na(data_clean$ephrtever)),]$site) & ephrtdur==0, 0, ephrtever)) %>%
  mutate(firstbirthage=ifelse(site %in% unique(data_clean[which(!is.na(data_clean$firstbirthage)),]$site) & parity==0, 0, firstbirthage)) 

# -- e. Creating dataset with no missing values --
train_clinical_nomiss <- train_clinical %>%
  mutate(firstbirthage=ifelse(site %in% unique(train_clinical[which(!is.na(train_clinical$firstbirthage)),]$site), firstbirthage, 999)) %>%
  mutate(lastbirthage=ifelse(site %in% unique(train_clinical[which(!is.na(train_clinical$lastbirthage)),]$site), lastbirthage, 999)) %>%
  mutate(menarchage=ifelse(site %in% unique(train_clinical[which(!is.na(train_clinical$menarchage)),]$site), menarchage, 999)) %>%
  mutate(pregs=ifelse(site %in% unique(train_clinical[which(!is.na(train_clinical$pregs)),]$site), pregs, 999)) %>%
  mutate(packyears=ifelse(site %in% unique(train_clinical[which(!is.na(train_clinical$packyears)),]$site), packyears, 999)) %>%
  mutate(parity=ifelse(site %in% unique(train_clinical[which(!is.na(train_clinical$parity)),]$site), parity, 999)) %>%
  mutate(alcohol=ifelse(site %in% unique(train_clinical[which(!is.na(train_clinical$alcohol)),]$site), alcohol, 999)) %>%
  mutate(ephrtdur=ifelse(site %in% unique(train_clinical[which(!is.na(train_clinical$ephrtdur)),]$site), ephrtdur, 999)) %>%
  mutate(ehrtdur=ifelse(site %in% unique(train_clinical[which(!is.na(train_clinical$ehrtdur)),]$site), ehrtdur, 999)) %>%
  mutate(ocdur=ifelse(site %in% unique(train_clinical[which(!is.na(train_clinical$ocdur)),]$site), ocdur, 999)) %>%
  mutate(education=ifelse(site %in% unique(train_clinical[which(!is.na(train_clinical$education)),]$site), education, 999)) %>%
  mutate(smoking=ifelse(site %in% unique(train_clinical[which(!is.na(train_clinical$smoking)),]$site), smoking, 999)) %>%  
  mutate(anyhrtever=ifelse(site %in% unique(train_clinical[which(!is.na(train_clinical$anyhrtever)),]$site), anyhrtever, 999)) %>%
  mutate(ephrtever=ifelse(site %in% unique(train_clinical[which(!is.na(train_clinical$ephrtever)),]$site), ephrtever, 999)) %>%
  mutate(ehrtever=ifelse(site %in% unique(train_clinical[which(!is.na(train_clinical$ehrtever)),]$site), ehrtever, 999)) %>%
  mutate(ocever=ifelse(site %in% unique(train_clinical[which(!is.na(train_clinical$ocever)),]$site), ocever, 999)) %>%
  mutate(diabetes=ifelse(site %in% unique(train_clinical[which(!is.na(train_clinical$diabetes)),]$site), diabetes, 999)) %>%
  mutate(endom=ifelse(site %in% unique(train_clinical[which(!is.na(train_clinical$endom)),]$site), endom, 999)) %>%
  mutate(hypertension=ifelse(site %in% unique(train_clinical[which(!is.na(train_clinical$hypertension)),]$site), hypertension, 999)) %>%
  mutate(famhx=ifelse(site %in% unique(train_clinical[which(!is.na(train_clinical$famhx)),]$site), famhx, 999)) %>%
  mutate(agesq=age^2) %>%
  mutate(bmisq=bmi^2) %>%
  mutate(paritysq=parity^2) %>%
  mutate(menarchagesq=menarchage^2) %>%
  mutate(ephrtdur.int=ifelse(ephrtdur==999, 0, ephrtdur)) %>%
  mutate(ephrtdur.sqint=(ephrtdur.int^2)) %>%
  mutate(ephrtdur.missing=ifelse(ephrtdur==999, 1, 0)) %>%
  mutate(ephrtdur.group=ifelse(ephrtdur>0 & ephrtdur<=5, 1, ifelse(ephrtdur>5 & 
    ephrtdur<=10, 2, ifelse(ephrtdur!=999 & ephrtdur!=0, 3, ephrtdur)))) %>%
  mutate(ehrtdur.int=ifelse(ehrtdur==999, 0, ehrtdur)) %>%
  mutate(ehrtdur.sqint=(ehrtdur.int^2)) %>%
  mutate(ehrtdur.missing=ifelse(ehrtdur==999, 1, 0)) %>%
  mutate(ehrtdur.group=ifelse(ehrtdur>0 & ehrtdur<=5, 1, ifelse(ehrtdur>5 & 
    ehrtdur<=10, 2, ifelse(ehrtdur!=999 & ehrtdur!=0, 3, ehrtdur)))) %>%
  mutate(firstbirthage.int=ifelse(firstbirthage==999, 0, firstbirthage)) %>%
  mutate(firstbirthage.sqint=(firstbirthage.int^2)) %>%
  mutate(firstbirthage.missing1=ifelse(firstbirthage==999, 1, 0)) %>%
  mutate(firstbirthage.missing2=ifelse(firstbirthage==0, 1, 0)) %>%
  mutate(firstbirthage.group=as.numeric(as.character(cut(firstbirthage, 
    breaks=c(-Inf, 1, 20, 25, 30, 35, 40, 998, Inf), labels=c(9, 1:6, 999), right=F)))) %>%
  mutate(ocdur.int=ifelse(ocdur==999, 0, ocdur)) %>%
  mutate(ocdur.sqint=(ocdur.int^2)) %>%
  mutate(ocdur.missing=ifelse(ocdur==999, 1, 0)) %>%  
  mutate(ocdur.group=ifelse(ocdur>0 & ocdur<=5, 1, ifelse(ocdur>5 & ocdur<=10, 
    2, ifelse(ocdur!=999 & ocdur!=0, 3, ocdur)))) %>%
  mutate(age.group=cut(age, breaks=c(-Inf, seq(50, 85, by=5), Inf), 
    labels=c(49, seq(50, 85, by=5)), right=F)) %>%
  mutate(bmi.group=cut(bmi, breaks=c(-Inf, 18.5, 25, 30, 35, Inf), right=F, 
    labels=c(1, 2, 3, 4, 5))) %>%
  mutate(parity.group=ifelse(parity>=4 & parity!=999, 4, parity)) %>%
  mutate(menarchage.group=ifelse(menarchage<=10, 10, ifelse(menarchage==12.75, 12, 
    ifelse(menarchage>=16 & menarchage!=999, 16, menarchage)))) %>% 
  mutate(menarchage2.group=ifelse(menarchage<=9, 9, ifelse(menarchage %in% c(10, 11), 
    10, ifelse(menarchage %in% c(12, 12.75, 13), 12, ifelse(menarchage %in% c(14, 15),
    14, ifelse(menarchage >= 16, 16, NA)))))) %>%
  # Interaction between ever OC use and BMI
  mutate(ocever.bmi=ocever*bmi) %>%
  mutate(ocever.bmisq=ocever*bmisq) %>%
  mutate(ocever.bmigroup2=ifelse(bmi.group==2, ocever, 0)) %>%
  mutate(ocever.bmigroup3=ifelse(bmi.group==3, ocever, 0)) %>%
  mutate(ocever.bmigroup4=ifelse(bmi.group==4, ocever, 0)) %>%
  mutate(ocever.bmigroup5=ifelse(bmi.group==5, ocever, 0)) %>%
  mutate(ocever.bmigroup6=ifelse(bmi.group==6, ocever, 0)) %>%  
  # Interaction between duration of OC use and BMI
  mutate(ocdur1.bmigroup1=ifelse(ocdur.group==1 & bmi.group==1, 1, 0)) %>%
  mutate(ocdur1.bmigroup3=ifelse(ocdur.group==1 & bmi.group==3, 1, 0)) %>%
  mutate(ocdur1.bmigroup4=ifelse(ocdur.group==1 & bmi.group==4, 1, 0)) %>%
  mutate(ocdur1.bmigroup5=ifelse(ocdur.group==1 & bmi.group==5, 1, 0)) %>%
  mutate(ocdur1.bmigroup6=ifelse(ocdur.group==1 & bmi.group==6, 1, 0)) %>%
  mutate(ocdur2.bmigroup1=ifelse(ocdur.group==2 & bmi.group==1, 1, 0)) %>%
  mutate(ocdur2.bmigroup3=ifelse(ocdur.group==2 & bmi.group==3, 1, 0)) %>%
  mutate(ocdur2.bmigroup4=ifelse(ocdur.group==2 & bmi.group==4, 1, 0)) %>%
  mutate(ocdur2.bmigroup5=ifelse(ocdur.group==2 & bmi.group==5, 1, 0)) %>%
  mutate(ocdur2.bmigroup6=ifelse(ocdur.group==2 & bmi.group==6, 1, 0)) %>%
  mutate(ocdur3.bmigroup1=ifelse(ocdur.group==3 & bmi.group==1, 1, 0)) %>%
  mutate(ocdur3.bmigroup3=ifelse(ocdur.group==3 & bmi.group==3, 1, 0)) %>%
  mutate(ocdur3.bmigroup4=ifelse(ocdur.group==3 & bmi.group==4, 1, 0)) %>%
  mutate(ocdur3.bmigroup5=ifelse(ocdur.group==3 & bmi.group==5, 1, 0)) %>%
  mutate(ocdur3.bmigroup6=ifelse(ocdur.group==3 & bmi.group==6, 1, 0)) %>%  
  # Interaction between any HT use and BMI 
  mutate(anyhrtever.bmi=anyhrtever*bmi) %>%
  mutate(anyhrtever.bmisq=anyhrtever*bmisq) %>%
  mutate(anyhrtever.bmigroup2=ifelse(bmi.group==2, anyhrtever, 0)) %>%
  mutate(anyhrtever.bmigroup3=ifelse(bmi.group==3, anyhrtever, 0)) %>%
  mutate(anyhrtever.bmigroup4=ifelse(bmi.group==4, anyhrtever, 0)) %>%
  mutate(anyhrtever.bmigroup5=ifelse(bmi.group==5, anyhrtever, 0)) %>%
  mutate(anyhrtever.bmigroup6=ifelse(bmi.group==6, anyhrtever, 0)) %>%
  # Interaction between e-only HT use and BMI
  mutate(ehrtever.bmi=ifelse(ehrtever==1, bmi, 0)) %>%
  mutate(ehrtever.bmisq=ifelse(ehrtever==1, bmisq, 0)) %>%
  mutate(ehrtever.bmigroup2=ifelse(bmi.group==2 & ehrtever==1, 1, 0)) %>%
  mutate(ehrtever.bmigroup3=ifelse(bmi.group==3 & ehrtever==1, 1, 0)) %>%
  mutate(ehrtever.bmigroup4=ifelse(bmi.group==4 & ehrtever==1, 1, 0)) %>%
  mutate(ehrtever.bmigroup5=ifelse(bmi.group==5 & ehrtever==1, 1, 0)) %>%
  mutate(ehrtever.bmigroup6=ifelse(bmi.group==6 & ehrtever==1, 1, 0)) %>%
  # Interaction between EP HT use and BMI
  mutate(ephrtever.bmi=ifelse(ephrtever==1, bmi, 0)) %>%
  mutate(ephrtever.bmisq=ifelse(ephrtever==1, bmisq, 0)) %>%
  mutate(ephrtever.bmigroup2=ifelse(bmi.group==2 & ephrtever==1, 1, 0)) %>%
  mutate(ephrtever.bmigroup3=ifelse(bmi.group==3 & ephrtever==1, 1, 0)) %>%
  mutate(ephrtever.bmigroup4=ifelse(bmi.group==4 & ephrtever==1, 1, 0)) %>%
  mutate(ephrtever.bmigroup5=ifelse(bmi.group==5 & ephrtever==1, 1, 0)) %>%
  mutate(ephrtever.bmigroup6=ifelse(bmi.group==6 & ephrtever==1, 1, 0)) %>%
  # Interaction between duration of e-only HT use and BMI
  mutate(ehrtdur1.bmigroup1=ifelse(ehrtdur.group==1 & bmi.group==1, 1, 0)) %>%
  mutate(ehrtdur1.bmigroup3=ifelse(ehrtdur.group==1 & bmi.group==3, 1, 0)) %>%
  mutate(ehrtdur1.bmigroup4=ifelse(ehrtdur.group==1 & (bmi.group==4|bmi.group==5), 1, 0)) %>%
  mutate(ehrtdur2.bmigroup1=ifelse(ehrtdur.group==2 & bmi.group==1, 1, 0)) %>%
  mutate(ehrtdur2.bmigroup3=ifelse(ehrtdur.group==2 & bmi.group==3, 1, 0)) %>%
  mutate(ehrtdur2.bmigroup4=ifelse(ehrtdur.group==2 & (bmi.group==4|bmi.group==5), 1, 0)) %>%
  mutate(ehrtdur3.bmigroup1=ifelse(ehrtdur.group==3 & bmi.group==1, 1, 0)) %>%
  mutate(ehrtdur3.bmigroup3=ifelse(ehrtdur.group==3 & bmi.group==3, 1, 0)) %>%
  mutate(ehrtdur3.bmigroup4=ifelse(ehrtdur.group==3 & (bmi.group==4|bmi.group==5), 1, 0)) %>%
  # Interaction between duration of EP HT use and BMI
  mutate(ephrtdur1.bmigroup1=ifelse(ephrtdur.group==1 & bmi.group==1, 1, 0)) %>%
  mutate(ephrtdur1.bmigroup3=ifelse(ephrtdur.group==1 & bmi.group==3, 1, 0)) %>%
  mutate(ephrtdur1.bmigroup4=ifelse(ephrtdur.group==1 & (bmi.group==4|bmi.group==5), 1, 0)) %>%
  mutate(ephrtdur2.bmigroup1=ifelse(ephrtdur.group==2 & bmi.group==1, 1, 0)) %>%
  mutate(ephrtdur2.bmigroup3=ifelse(ephrtdur.group==2 & bmi.group==3, 1, 0)) %>%
  mutate(ephrtdur2.bmigroup4=ifelse(ephrtdur.group==2 & (bmi.group==4|bmi.group==5), 1, 0)) %>%
  mutate(ephrtdur3.bmigroup1=ifelse(ephrtdur.group==3 & bmi.group==1, 1, 0)) %>%
  mutate(ephrtdur3.bmigroup3=ifelse(ephrtdur.group==3 & bmi.group==3, 1, 0)) %>%
  mutate(ephrtdur3.bmigroup4=ifelse(ephrtdur.group==3 & (bmi.group==4|bmi.group==5), 1, 0)) %>%
  mutate(ehrtdur.bmigroup=(as.numeric(as.character(bmi.group))-1)*ehrtdur.group) %>%
  mutate(ehrtdur.bmigroup=ifelse(ehrtdur.bmigroup>=999, 0, ehrtdur.bmigroup)) %>%
  mutate(ehrtdur.bmigroup999=ifelse(ehrtdur.group==999, 1, 0)) %>%
  dummy_cols(., select_columns=c("race", "education", "smoking", "menstat", "grade", 
    "site", "age.group", "bmi.group", "parity.group", "menarchage.group", 
    "menarchage2.group", "firstbirthage.group", "famhx", "diabetes", 
    "hypertension", "ephrtever", "ehrtever", "ephrtdur.group", "ehrtdur.group", 
    "ocdur.group"), ignore_na=T)

# -- f.  Renaming dummy variables --
names(train_clinical_nomiss) <- sub("_", "", names(train_clinical_nomiss))

train_clinical <- train_clinical %>%
  dummy_cols(., select_columns=c("race", "education", "smoking", "menstat", "grade", "site"), ignore_na=T)
names(train_clinical) <- sub("_", "", names(train_clinical))


# ------------------------ (3) Exporting cleaned data --------------------------
write.csv(train_clinical, "../Analysis/Output/E2C2_Development_1_Clean_E2C2_All.csv", row.names=F)
write.csv(train_clinical_nomiss, "../Analysis/Output/E2C2_Development_1_Clean_E2C2_NoMiss.csv", row.names=F)
write.csv(exclusions, "../Analysis/Output/E2C2_Development_1_Clean_E2C2_Exclusions.csv", row.names=F)

