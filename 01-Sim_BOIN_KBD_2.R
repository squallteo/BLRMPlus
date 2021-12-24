rm(list=ls())

library(BOIN)
library(Keyboard)
library(tidyverse)
library(readxl)

scenariodt <- rbind(read_xlsx("FixedScenarios.xlsx", sheet = "Sheet1"), 
                    read_xlsx("FixedScenarios.xlsx", sheet = "All"))

ProvDose <- subset(scenariodt, Scenario=="Steep")$Dose

nsim <- 1000

cohortsize <- 3
# ncohortvec <- c(27, 33, 33, 42, 9)/cohortsize; p.target <- 0.25; p.saf <- 0.16; p.tox <- 0.33 #[0.16, 0.33]
ncohortvec <- c(30, 33, 36, 42, 9)/cohortsize; p.target <- 0.25; p.saf <- 0.2; p.tox <- 0.3 #[0.2, 0.3]

for(i in 1:length(unique(scenariodt$Scenario))){
  ncohort <- ncohortvec[i]
  s <- unique(scenariodt$Scenario)[i]
  
  TrueRates <- subset(scenariodt, Scenario==s)$Rate
  
  oc.single <- get.oc(target =p.target, p.saf = p.saf, p.tox = p.tox,
                      p.true = TrueRates,
                      ncohort = ncohort, cohortsize = cohortsize, ntrial = nsim)
  tt1 <- tibble(Dose = "AllToxic", MTDFreq = oc.single$percentstop/100, Npat = 0)
  tt2 <- tibble(Dose = ProvDose, MTDFreq = oc.single$selpercent/100, Npat = oc.single$npatients)
  summarydt1 <- 
    rbind(tt1, tt2) %>% mutate(Noverall = oc.single$totaln, DLTRate = oc.single$totaltox/oc.single$totaln, 
                               Design = "BOIN", Scenario = s)
  
  
  oc.single <- get.oc.kb(target =p.target, marginL = p.target - p.saf, marginR = p.tox - p.target, 
                         p.true = TrueRates, 
                         ncohort = ncohort, cohortsize = cohortsize, ntrial = nsim)
  tt1 <- tibble(Dose = "AllToxic", MTDFreq = oc.single$percentstop/100, Npat = 0)
  tt2 <- tibble(Dose = ProvDose, MTDFreq = oc.single$selpercent/100, Npat = oc.single$npatients)
  summarydt2 <- 
    rbind(tt1, tt2) %>% mutate(Noverall = oc.single$totaln, DLTRate = oc.single$totaltox/oc.single$totaln, 
                               Design = "Keyboard", Scenario = s)
  
  tt <- rbind(summarydt1, summarydt2)
  if(i==1){outdt <- tt}
  else{outdt <- rbind(outdt, tt)}
}

# 
# write.csv(outdt, "BOIN_KBD_16_33.csv", row.names = F)
# write.csv(outdt, "BOIN_KBD_20_30.csv", row.names = F)





