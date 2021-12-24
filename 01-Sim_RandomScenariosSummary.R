rm(list=ls())
library(tidyverse)
library(readxl)
rootpath <- getwd()

# randclass <- "Clertant"
randclass <- "Paoletti"
interval <- "16_33"; Pint_BLRM <- c(0, 0.16, 0.33, 1)
# interval <- "20_30"; Pint_BLRM <- c(0, 0.20, 0.30, 1)
design <- 0:4

toxdt <- read_csv(paste(randclass,"Class.csv",sep="")) %>% group_by(Sim) %>% 
  mutate(AllToxic=all(Rate>=Pint_BLRM[3]), AllUnder=all(Rate<Pint_BLRM[2]), 
         Target=(Rate >= Pint_BLRM[2] & Rate < Pint_BLRM[3]))
DoseProv <- sort(unique(toxdt$Dose))
#part 1: has MTD
part1 <- toxdt %>% filter(Target) %>% rename(TrueMTD=Dose) %>% select(Sim, TrueMTD)
#part 2: all toxic, MTD = -1
part2 <- toxdt %>% filter(AllToxic & Dose==DoseProv[1]) %>% mutate(TrueMTD=-1) %>% select(Sim, TrueMTD)
#part 3: all under, MTD = 999999
part3 <- toxdt %>% filter(AllUnder & Dose==DoseProv[1]) %>% mutate(TrueMTD=999999) %>% select(Sim, TrueMTD)
#in all other scenarios, neither all toxic nor under, just no MTD in the interval: MTD=999999
MTDTrue <- rbind(part1, part2, part3) %>% ungroup() %>% full_join(tibble(Sim=1:length(unique(toxdt$Sim)))) %>% 
  mutate(TrueMTD=ifelse(is.na(TrueMTD), 999999, TrueMTD)) %>% arrange(Sim, TrueMTD)

for(d in design){
  filename <- paste(randclass, "_design", d, ".RData", sep="")
  fullpath <- paste(rootpath, "Results", paste(randclass,interval,sep="_"), filename, sep = "/")
  
  load(fullpath)
  
  Ntotdt <- resultdt %>% group_by(Sim) %>% summarize_at("Npat", sum) %>% rename(Ntotal = Npat)
  tt1 <- resultdt %>% group_by(Sim, Dose) %>% summarize_at("Npat", sum) %>% left_join(Ntotdt, by = 'Sim') 
  
  pctdt <-  tt1 %>% right_join(MTDTrue, by = "Sim") %>% filter(Dose==TrueMTD) %>% mutate(Pct = Npat/Ntotal) %>% group_by(Sim) %>% 
    summarize_at("Pct", sum) %>% right_join(Ntotdt, by = "Sim") %>% arrange(Sim) %>% mutate(Pct=ifelse(is.na(Pct), 0, Pct))
  print(mean(pctdt$Pct))
  
}

