################################################################################
#                                                                              #
# Endometrial Cancer Risk Prediction Model: Obtaining Baseline Endo Rates      #
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
if (!require('haven')) install.packages('haven'); library('haven')

start.age <- 45

run.brfss <- 0
nhsi.brfssyear <- 1988
nhsi.seeryear <- "1989-1993"
nhsii.brfssyear <- c(2006, 2008)
nhsii.seeryear <- "2003-2007"
plco.brfssyear <- c(1996, 1998) 
plco.seeryear <- "1996-2000"
current.brfssyear <- c(2016, 2018)
current.seeryear <- "2013-2017"


# ------------------- (2) Importing and cleaning BRFSS data --------------------

if (1==run.brfss){
brfss <- lapply(Sys.glob("BRFSS Data/*.XPT"), function(x){
  dat <- read.xport(x) %>% 
    select(contains("SEX"), contains("AGE"), contains("FINALWT"), contains("RACE"), 
           contains("HADHYST"), contains("LLCPWT"), contains("STATE"))
  dat$filename <- tools::file_path_sans_ext(basename(x))
  dat
}) %>%
  bind_rows() %>% 
  mutate(year=as.numeric(str_sub(filename, -2, -1)),
         year=ifelse(year<20, 2000+year, 1900+year),
         race=ifelse(year<=2000, X_RACEG, ifelse(year<=2012, X_RACEGR2, X_RACEG21)),
         age=ifelse(year<=2012, AGE, X_AGE80),
         hyst=ifelse(year<=2000, HADHYST, HADHYST2),
         hyst=ifelse(hyst==2, 0, hyst),
         weight=ifelse(year<=2010, X_FINALWT, X_LLCPWT),
         sex=ifelse(year<=2016, SEX, SEX1)) %>%
  filter(sex==2 & race==1 & age>=18 & (hyst==1|hyst==0)) %>%
  filter(X_STATE %in% c(2, 6, 9, 13, 15, 16, 19, 21, 22, 25, 26, 34, 35, 36, 49, 53)) %>%
  select(year, age, weight, hyst) %>%
  arrange(year, age) %>%
  group_by(year, age) %>%
  mutate(hyst=weighted.mean(hyst, w=weight)) %>%
  mutate(weight=sum(weight)) %>%
  distinct()
}

# --------------------- (3) Exporting cleaned BRFSS data -----------------------
if (1==run.brfss){
write.table(brfss, "../Analysis/Output/E2C2_Development_2_Clean_BRFSS.txt", sep="\t")
}


# ------------------------ (4) Re-importing BRFSS data -------------------------

hyst.prev <- read.table("../Analysis/Output/E2C2_Development_2_Clean_BRFSS.txt") %>% 
  rename(brfss.year=year) %>%
  mutate(age.group=cut(age, breaks=c(seq(15, 85, by=5), Inf), labels=seq(15, 85, by=5), right=F)) %>%
  group_by(brfss.year, age.group) %>%
  mutate(hyst=weighted.mean(hyst, w=weight)) %>%
  ungroup() %>%
  distinct(brfss.year, age.group, .keep_all=T) %>%
  dplyr::select(brfss.year, age.group, hyst) %>%
  arrange(brfss.year, age.group) %>%
  mutate(age.group=as.numeric(as.character(age.group)))


# ------------------------- (5) Importing SEER data ----------------------------

seer.endo <- read.csv("SEER Endometrial Cancer Incidence Rates.csv") %>%
  rename(age.group=Age_group) %>%
  mutate(age.group=as.numeric(str_sub(age.group, 1, 2))) %>%
  filter(age.group>=15) %>%
  rename(seer.year=Year, rate=Rate)


# -------------------- (6) Combining BRFSS and SEER data -----------------------

endo.rates <- full_join(hyst.prev, seer.endo, by="age.group")

# -- a. Cleaning data for NHS I -- 
nhsi.disease.incidence.rates <- endo.rates %>%
  filter(brfss.year %in% nhsi.brfssyear & seer.year==nhsi.seeryear) %>%
  group_by(age.group) %>%
  mutate(hystavg=mean(hyst)) %>%
  ungroup() %>%
  filter(brfss.year==nhsi.brfssyear[1]) %>%
  mutate(rate.corrected=rate/(1-hystavg)/100000) %>%
  full_join(data.frame(age.group=seq(40, 85)), by="age.group") %>%
  arrange(age.group) %>%
  fill(rate.corrected) %>%
  dplyr::select(age.group, rate.corrected) %>%
  filter(age.group>=40 & age.group<=85) %>%
  as.matrix()

# -- b. Cleaning data for NHS II -- 
nhsii.disease.incidence.rates <- endo.rates %>%
  filter(brfss.year %in% nhsii.brfssyear & seer.year==nhsii.seeryear) %>%
  group_by(age.group) %>%
  mutate(hystavg=mean(hyst)) %>%
  ungroup() %>%
  filter(brfss.year==nhsii.brfssyear[1]) %>%
  mutate(rate.corrected=rate/(1-hystavg)/100000) %>%
  full_join(data.frame(age.group=seq(40, 85)), by="age.group") %>%
  arrange(age.group) %>%
  fill(rate.corrected) %>%
  dplyr::select(age.group, rate.corrected) %>%
  filter(age.group>=40 & age.group<=85) %>%
  as.matrix()

# -- c. Cleaning data for PLCO -- 
plco.disease.incidence.rates <- endo.rates %>%
  filter(brfss.year %in% plco.brfssyear & seer.year==plco.seeryear) %>%
  group_by(age.group) %>%
  mutate(hystavg=mean(hyst)) %>%
  ungroup() %>%
  filter(brfss.year==plco.brfssyear[1]) %>%
  mutate(rate.corrected=rate/(1-hystavg)/100000) %>%
  full_join(data.frame(age.group=seq(40, 85)), by="age.group") %>%
  arrange(age.group) %>%
  fill(rate.corrected) %>%
  dplyr::select(age.group, rate.corrected) %>%
  filter(age.group>=40 & age.group<=85) %>%
  as.matrix()

# -- d. Cleaning data for current rates -- 
current.disease.incidence.rates <- endo.rates %>%
  filter(brfss.year %in% current.brfssyear & seer.year==current.seeryear) %>%
  group_by(age.group) %>%
  mutate(hystavg=mean(hyst)) %>%
  ungroup() %>%
  filter(brfss.year==current.brfssyear[1]) %>%
  mutate(rate.corrected=rate/(1-hystavg)/100000) %>%
  full_join(data.frame(age.group=seq(40, 85)), by="age.group") %>%
  arrange(age.group) %>%
  fill(rate.corrected) %>%
  dplyr::select(age.group, rate.corrected) %>%
  filter(age.group>=40 & age.group<=85) %>%
  as.matrix()


# ------------------------- (7) Exporting endo rates ---------------------------
write.csv(endo.rates, "../Analysis/Output/E2C2_Development_2_Endo_All.csv", row.names=F)
write.csv(nhsi.disease.incidence.rates, "../Analysis/Output/E2C2_Development_2_Endo_NHSI.csv", row.names=F)
write.csv(nhsii.disease.incidence.rates, "../Analysis/Output/E2C2_Development_2_Endo_NHSII.csv", row.names=F)
write.csv(plco.disease.incidence.rates, "../Analysis/Output/E2C2_Development_2_Endo_PLCO.csv", row.names=F)
write.csv(current.disease.incidence.rates, "../Analysis/Output/E2C2_Development_2_Endo_Current.csv", row.names=F)

