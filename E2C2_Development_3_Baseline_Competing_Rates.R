################################################################################
#                                                                              #
# Endometrial Cancer Risk Prediction Model: Obtaining Baseline Rates of        #
# COmpeting Risks                                                              #
#                                                                              #
# Author: Joy Shi                                                              #
# Last updated: 09/15/2021                                                     #
#                                                                              #
# Purpose of program: combining CDC mortality data, SEER data and BRFSS data   #
#                     to obtain overall rate of competing risks (mortality,    #
#                     other cancers, hysterectomy)                             #
#                                                                              #
################################################################################


# -------------------------------- (1) Set up ----------------------------------

setwd("D:/Dropbox/PhD Dissertation/Paper 3 - E2C2/Datasets")

if (!require('tidyverse')) install.packages('tidyverse'); library('tidyverse')

start.age <- 45

nhsi.seeryear <- "1989-1993"
nhsii.seeryear <- "2003-2007"
plco.seeryear <- "1996-2000"
current.seeryear <- "2013-2017"

nhsi.cdcyear <- 1988
nhsii.cdcyear <- 2004
plco.cdcyear <- 1997
current.cdcyear <- 2017


# --------------------------- (2) CDC mortality data ---------------------------

mortality.rates <- read.csv("CDC Wonder Mortality Rates.csv") %>% 
  mutate(Mortality.Rate=Rate/100000) %>% rename(CDC.Year=Year) %>%
  full_join(expand.grid(CDC.Year=unique(.$CDC.Year), 
    Age=seq(40, 84)), by=c("CDC.Year", "Age")) %>%
  arrange(CDC.Year, Age) %>%
  fill(Mortality.Rate)


# --------------------------- (3) SEER data cancers -----------------------------

allcancers.rates <- read.csv("SEER All Cancers Incidence Rates.csv") %>%
  rename(Age=Age_group) %>%
  mutate(Age=as.numeric(str_sub(Age, 1, 2))) %>%
  filter(Age>=40) %>%
  mutate(Cancer.Rate=as.numeric(gsub(",", "", as.character(Rate)))/100000) %>%
  full_join(data.frame(expand.grid(Year=c("1989-1993", "1996-2000", "2003-2007", "2013-2017"), 
    Age=seq(40,85))), by=c("Year", "Age")) %>%
  arrange(Year, Age) %>%
  fill(Cancer.Rate) %>%
  dplyr::select(Year, Age, Cancer.Rate) %>%
  rename(SEER.Year=Year)

nhsi.disease.incidence.rates <- read.csv("../Analysis/Output/E2C2_Development_2_Endo_NHSI.csv")
nhsii.disease.incidence.rates <- read.csv("../Analysis/Output/E2C2_Development_2_Endo_NHSII.csv")
plco.disease.incidence.rates <- read.csv("../Analysis/Output/E2C2_Development_2_Endo_PLCO.csv")
current.disease.incidence.rates <- read.csv("../Analysis/Output/E2C2_Development_2_Endo_Current.csv")


# ---------------------- (4) BRFSS data on hysterectomy ------------------------

brfss.rollingavg <- read.table("../Analysis/Output/E2C2_Development_2_Clean_BRFSS.txt") %>%
  group_by(year) %>%
  mutate(hyst.avg=(lag(hyst,2)+lag(hyst)+hyst+lead(hyst)+lead(hyst,2))/5) %>%
  ungroup() %>%
  arrange(age, year) %>%
  group_by(age) %>%
  mutate(hyst.avg=(lag(hyst.avg)+hyst.avg+lead(hyst.avg))/3) %>%
  ungroup()

hyst.incidence <- mortality.rates %>%
  arrange(CDC.Year, Age) %>% 
  mutate(age=as.numeric(as.character(Age))) %>%
  rename(year=CDC.Year) %>%
  full_join(brfss.rollingavg, by=c("age", "year")) %>%
  filter(!is.na(hyst) & age>=40 & age<85) %>%
  mutate(yob=year-age) %>%
  arrange(yob, age) %>%
  group_by(yob) %>%  
  mutate(hyst.lead=lead(hyst.avg)) %>%
  ungroup() %>%
  mutate(hyst.incidence=log((1-hyst.lead)*hyst.avg*exp(-Mortality.Rate*2)/
    (1-hyst.avg)+(1-hyst.lead)*exp(-Mortality.Rate*2))/(-2)-Mortality.Rate) %>%
  arrange(year, age) %>%
  mutate(hyst.iavg=(lag(hyst.incidence,2)+lag(hyst.incidence,1)+hyst.incidence+
    lead(hyst.incidence)+lead(hyst.incidence,2))/5) %>%
  mutate(hyst.iavg=ifelse(hyst.iavg<0, 0, hyst.iavg)) %>%
  dplyr::select(year, age, hyst.iavg)


# --------------------------- (5) Combining rates ------------------------------

nhsi.competingrisks <- full_join(mortality.rates, allcancers.rates, by="Age") %>%
  filter(SEER.Year==nhsi.seeryear) %>%
  filter(CDC.Year==nhsi.cdcyear) %>%
  left_join(data.frame(nhsi.disease.incidence.rates), by=c("Age"="age.group")) %>%
  left_join(hyst.incidence, by=c("Age"="age")) %>%
  filter(year==1990) %>%
  mutate(Competing.Risk=Mortality.Rate+Cancer.Rate+hyst.iavg-rate.corrected) %>%
  dplyr::select(Age, Competing.Risk) %>%
  add_row(Age=85) %>%
  fill(Competing.Risk) %>%
  arrange(-Age) %>%
  fill(Competing.Risk)

nhsii.competingrisks <- full_join(mortality.rates, allcancers.rates, by="Age") %>%
  filter(SEER.Year==nhsii.seeryear) %>%
  filter(CDC.Year==nhsii.cdcyear) %>%
  left_join(data.frame(nhsii.disease.incidence.rates), by=c("Age"="age.group")) %>%
  left_join(hyst.incidence, by=c("CDC.Year"="year", "Age"="age")) %>%
  mutate(Competing.Risk=Mortality.Rate+Cancer.Rate+hyst.iavg-rate.corrected) %>%
  # mutate(Competing.Risk=Mortality.Rate+Cancer.Rate-rate.corrected) %>%
  dplyr::select(Age, Competing.Risk) %>%
  add_row(Age=85) %>%
  fill(Competing.Risk) %>%
  arrange(-Age) %>%
  fill(Competing.Risk)

plco.competingrisks <- full_join(mortality.rates, allcancers.rates, by="Age") %>%
  filter(SEER.Year==plco.seeryear) %>%
  filter(CDC.Year==plco.cdcyear) %>%
  left_join(data.frame(plco.disease.incidence.rates), by=c("Age"="age.group")) %>%
  left_join(hyst.incidence, by=c("Age"="age")) %>%
  filter(year %in% c(1996, 1998)) %>%
  group_by(Age) %>%
  mutate(Competing.Risk=Mortality.Rate+Cancer.Rate+mean(hyst.iavg)-rate.corrected)  %>%
  ungroup() %>%
  # mutate(Competing.Risk=Mortality.Rate+Cancer.Rate-rate.corrected) %>%
  dplyr::select(Age, Competing.Risk) %>%
  add_row(Age=85) %>%
  distinct() %>%
  fill(Competing.Risk) %>%
  arrange(-Age) %>%
  fill(Competing.Risk) %>%
  data.frame()

current.competingrisks <- full_join(mortality.rates, allcancers.rates, by="Age") %>%
  filter(SEER.Year==current.seeryear) %>%
  filter(CDC.Year==current.cdcyear) %>%
  left_join(data.frame(current.disease.incidence.rates), by=c("Age"="age.group")) %>%
  left_join(hyst.incidence, by=c("Age"="age")) %>%
  filter(year==2014|Age>=80) %>%
  distinct(Age, .keep_all=T) %>%
  fill(hyst.iavg) %>%
  mutate(Competing.Risk=Mortality.Rate+Cancer.Rate+hyst.iavg-rate.corrected) %>%
  # mutate(Competing.Risk=Mortality.Rate+Cancer.Rate-rate.corrected) %>%
  dplyr::select(Age, Competing.Risk) %>%
  add_row(Age=85) %>%
  fill(Competing.Risk) %>%
  arrange(-Age) %>%
  fill(Competing.Risk)

# ----------------------- (6) Exporting competing risk -------------------------
write.csv(mortality.rates, "../Analysis/Output/E2C2_Development_3_Competing_Mortality.csv", row.names=F)
write.csv(allcancers.rates, "../Analysis/Output/E2C2_Development_3_Competing_Cancer.csv", row.names=F)
write.csv(nhsi.competingrisks, "../Analysis/Output/E2C2_Development_3_Competing_NHSI.csv", row.names=F)
write.csv(nhsii.competingrisks, "../Analysis/Output/E2C2_Development_3_Competing_NHSII.csv", row.names=F)
write.csv(plco.competingrisks, "../Analysis/Output/E2C2_Development_3_Competing_PLCO.csv", row.names=F)
write.csv(current.competingrisks, "../Analysis/Output/E2C2_Development_3_Competing_Current.csv", row.names=F)
