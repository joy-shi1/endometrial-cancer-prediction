################################################################################
#                                                                              #
# Endometrial Cancer Risk Prediction Model: Creating List of Predictors        #
#                                                                              #
# Author: Joy Shi                                                              #
# Last updated: 09/15/2021                                                     #
#                                                                              #
# Purpose of program: creating list of clinical and genetic predictors to be   #
#                     considered in the model                                  #
#                                                                              #
################################################################################

# -------------------------------- (1) Set up ----------------------------------

setwd("D:/Dropbox/PhD Dissertation/Paper 3 - E2C2/Datasets")



# --------------------- (2) List of clinical predictors ------------------------

predictors <- c("site104", "site105", "site109", "site111", "site202", 
                "site204", "site205", "site206", "site207", "site210", 
                "site303", "site304", "site305", "site306", 
                "site307", "site405", "site406",
                "age.group50", "age.group55", "age.group65", 
                "age.group70", "age.group75", "age.group80", "age.group85",
                "education2", "education3",
                "smoking2", "smoking3", "smoking999",
                "bmi.group1", "bmi.group3", "bmi.group4", "bmi.group5", 
                "parity.group1", "parity.group2", "parity.group3", "parity.group4",
                "menarchage2.group10", "menarchage2.group12", "menarchage2.group14", "menarchage2.group16",
                "anyhrtever", 
                "ehrtever1", "ehrtever999", "ehrtdur.group1", "ehrtdur.group2", "ehrtdur.group3", "ehrtdur.group999",
                "ephrtever1", "ephrtever999", "ephrtdur.group1", "ephrtdur.group2", "ephrtdur.group3", "ephrtdur.group999",
                "ocever", "ocdur.group1", "ocdur.group2", "ocdur.group3", "ocdur.group999",
                "firstbirthage.group2", "firstbirthage.group3", "firstbirthage.group4", "firstbirthage.group5", "firstbirthage.group9",
                "firstbirthage.group999",
                "diabetes1", "diabetes999", "hypertension1", "hypertension999",
                "ocever.bmigroup3", "ocever.bmigroup4", "ocever.bmigroup5",
                "ocdur1.bmigroup3", "ocdur1.bmigroup4", "ocdur1.bmigroup5",
                "ocdur2.bmigroup3", "ocdur2.bmigroup4", "ocdur2.bmigroup5",
                "ocdur3.bmigroup3", "ocdur3.bmigroup4", "ocdur3.bmigroup5",
                "anyhrtever.bmigroup3", "anyhrtever.bmigroup4", "anyhrtever.bmigroup5", 
                "ehrtever.bmigroup3", "ehrtever.bmigroup4", "ehrtever.bmigroup5", 
                "ehrtdur1.bmigroup3", "ehrtdur1.bmigroup4", 
                "ehrtdur2.bmigroup3", "ehrtdur2.bmigroup4", 
                "ehrtdur3.bmigroup3", "ehrtdur3.bmigroup4", 
                "ephrtever.bmigroup3", "ephrtever.bmigroup4", "ephrtever.bmigroup5", 
                "ephrtdur1.bmigroup3", "ephrtdur1.bmigroup4", 
                "ephrtdur2.bmigroup3", "ephrtdur2.bmigroup4", 
                "ephrtdur3.bmigroup3", "ephrtdur3.bmigroup4"
)


# ---------------------- (2) List of genetic predictors ------------------------

endo.snp.info <- read.csv("SNP_Reference.csv", colClasses="character") %>% 
  mutate(SNP37=paste("X", Chr, ".", Position37, sep="")) %>%
  mutate(SNP38=paste("X", Chr, ".", Position38, sep="")) %>%
  mutate(OR=as.numeric(substr(OR, 1, 4))) %>%
  mutate(snp.name=rsID) %>%
  mutate(snp.odds.ratio=as.numeric(OR)) %>%
  mutate(snp.freq=as.numeric(MAF)) %>%
  select(snp.name, snp.odds.ratio, snp.freq)


# ----------------------------- (3) Exporting data -----------------------------

write.csv(predictors, "../Analysis/Output/E2C2_Development_4_Predictors.csv", row.names=F)
write.csv(endo.snp.info, "../Analysis/Output/E2C2_Development_4_SNPs.csv", row.names=F)
