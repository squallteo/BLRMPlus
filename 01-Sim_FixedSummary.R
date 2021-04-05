rm(list=ls())
library(tidyverse)
rootpath <- getwd()

interval <- "16_33"
scenariodt <- read_xlsx("FixedScenarios.xlsx", sheet = "Sheet1")
design <- 0:4

for(s in unique(scenariodt$Scenario)){
  for(d in design){
    filename <- paste(s, "_design", d, ".csv", sep="")
    fullpath <- paste(rootpath, "Results", interval, filename, sep = "/")
    
    resultdt <- read_csv(fullpath)
    
    tt1 <- resultdt %>% select(Dose, MTDFreq, Npat) %>% mutate(Design=d)
    tt2 <- resultdt %>% filter(Dose==min(Dose)) %>% select(Noverall, DLTRate) %>% mutate(Design=d)
    if(d==min(design)){
      dosedt <- rbind(tt1)
      overalldt <- rbind(tt2)
    }
    else{
      dosedt <- rbind(dosedt, tt1)
      overalldt <- rbind(overalldt, tt2)
    }
  }
  out1 <- dosedt %>% mutate(Scenario=s)
  out2 <- overalldt %>% mutate(Scenario=s)
  
  write_csv(out1, paste("Fixed_", s, "_dose.csv", sep=""))
  write_csv(out2, paste("Fixed_", s, "_overall.csv", sep=""))
}

