rm(list=ls())

library(BOIN)
library(Keyboard)
library(tidyverse)
library(readxl)

scenariodt <- rbind(read_xlsx("FixedScenarios.xlsx", sheet = "Sheet1"), 
                    read_xlsx("FixedScenarios.xlsx", sheet = "All"))

ProvDose <- subset(scenariodt, Scenario=="Steep")$Dose

nsim <- 1000

ncohort <- 15
cohortsize <- 3
n.earlystop <- 9

# p.target <- 0.25; p.saf <- 0.16; p.tox <- 0.33 #[0.16, 0.33]
p.target <- 0.25; p.saf <- 0.2; p.tox <- 0.3 #[0.2, 0.3]

for(s in unique(scenariodt$Scenario)){
  
  TrueRates <- subset(scenariodt, Scenario==s)$Rate
  #BOIN with early stop
  oc.single <- get.oc(target =p.target, p.saf = p.saf, p.tox = p.tox,
                      p.true = TrueRates, n.earlystop = n.earlystop,
                      ncohort = ncohort, cohortsize = cohortsize, ntrial = nsim)
  
  tt1 <- tibble(Dose = "AllToxic", MTDFreq = oc.single$percentstop/100, Npat = 0)
  tt2 <- tibble(Dose = ProvDose, MTDFreq = oc.single$selpercent/100, Npat = oc.single$npatients)
  summarydt1 <- 
  rbind(tt1, tt2) %>% mutate(Noverall = oc.single$totaln, DLTRate = oc.single$totaltox/oc.single$totaln, 
                             Design = "BOIN_EarlyStop", Scenario = s)
  
  #BOIN without early stop
  oc.single <- get.oc(target =p.target, p.saf = p.saf, p.tox = p.tox,
                      p.true = TrueRates, n.earlystop = 1000,
                      ncohort = ncohort, cohortsize = cohortsize, ntrial = nsim)
  
  tt1 <- tibble(Dose = "AllToxic", MTDFreq = oc.single$percentstop/100, Npat = 0)
  tt2 <- tibble(Dose = ProvDose, MTDFreq = oc.single$selpercent/100, Npat = oc.single$npatients)
  summarydt2 <- 
    rbind(tt1, tt2) %>% mutate(Noverall = oc.single$totaln, DLTRate = oc.single$totaltox/oc.single$totaln, 
                               Design = "BOIN", Scenario = s)
  
  outdt <- rbind(summarydt1, summarydt2)
  
  if(s==scenariodt$Scenario[1]){
    resultdt <- outdt
  }
  else{
    resultdt <- rbind(resultdt, outdt)
  }
  
}

write.csv(resultdt, "BOINSim.csv", row.names = F)

for(s in unique(scenariodt$Scenario)){
  
  TrueRates <- subset(scenariodt, Scenario==s)$Rate
  #Keyboard with early stop
  oc.single <- get.oc.kb(target =p.target, marginL = p.target - p.saf, marginR = p.tox - p.target, 
                         p.true = TrueRates, n.earlystop = n.earlystop, 
                         ncohort = ncohort, cohortsize = cohortsize, ntrial = nsim)
  
  tt1 <- tibble(Dose = "AllToxic", MTDFreq = oc.single$percentstop/100, Npat = 0)
  tt2 <- tibble(Dose = ProvDose, MTDFreq = oc.single$selpercent/100, Npat = oc.single$npatients)
  summarydt1 <- 
    rbind(tt1, tt2) %>% mutate(Noverall = oc.single$totaln, DLTRate = oc.single$totaltox/oc.single$totaln, 
                               Design = "KBD_EarlyStop", Scenario = s)
  
  #Keyboard without early stop
  oc.single <- get.oc.kb(target =p.target, marginL = p.target - p.saf, marginR = p.tox - p.target, 
                         p.true = TrueRates, 
                         ncohort = ncohort, cohortsize = cohortsize, ntrial = nsim)
  
  tt1 <- tibble(Dose = "AllToxic", MTDFreq = oc.single$percentstop/100, Npat = 0)
  tt2 <- tibble(Dose = ProvDose, MTDFreq = oc.single$selpercent/100, Npat = oc.single$npatients)
  summarydt2 <- 
    rbind(tt1, tt2) %>% mutate(Noverall = oc.single$totaln, DLTRate = oc.single$totaltox/oc.single$totaln, 
                               Design = "KBD", Scenario = s)
  
  outdt <- rbind(summarydt1, summarydt2)
  
  if(s==scenariodt$Scenario[1]){
    resultdt <- outdt
  }
  else{
    resultdt <- rbind(resultdt, outdt)
  }
}

write.csv(resultdt, "KBDSim.csv", row.names = F)