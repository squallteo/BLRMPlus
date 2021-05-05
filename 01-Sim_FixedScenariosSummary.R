rm(list=ls())
library(tidyverse)
library(readxl)
library(ggpubr)
rootpath <- getwd()

# interval <- "16_33";Pint_BLRM <- c(0, 0.16, 0.33, 1)
interval <- "20_30";Pint_BLRM <- c(0, 0.2, 0.3, 1)
scenariodt <- read_xlsx("FixedScenarios.xlsx", sheet = "Sheet1")
design <- 0:4
maxdose <- 800

for(s in unique(scenariodt$Scenario)){
  for(d in design){
    filename <- paste(s, "_design", d, ".csv", sep="")
    fullpath <- paste(rootpath, "Results", interval, filename, sep = "/")
    
    resultdt <- read_csv(fullpath)
    tt <- resultdt %>% filter(MTDFlag=="*")
    trueMTD <- tt$Dose
    mLoss <-
    resultdt %>% mutate(Dose4Loss = ifelse(Dose==-1, 0, Dose),
                        Loss = abs(Dose4Loss - trueMTD)/maxdose,
                        Loss = ifelse(Dose>99999, 1, Loss)) %>%
      mutate(wLoss = MTDFreq*Loss) %>% summarize_at("wLoss", sum)
    
    tt1 <- resultdt %>% select(Dose, MTDFreq, Npat) %>% mutate(Design=d)
    tt2 <- resultdt %>% filter(Dose==min(Dose)) %>% select(Noverall, DLTRate) %>% mutate(Design=d, mLoss = mLoss)
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
  
  if(s==scenariodt$Scenario[1]){
    bydose <- rbind(out1)
    overall <- rbind(out2)
  }
  else{
    bydose <- rbind(bydose, out1)
    overall <- rbind(overall, out2)
  }
  # write_csv(out1, paste("Fixed_", s, "_dose.csv", sep=""))
  # write_csv(out2, paste("Fixed_", s, "_overall.csv", sep=""))
}

#Generate a 3*2 plot, rows are for scenarios, two columns: MTD accuracy and N treated at MTD
MTDdt <-
scenariodt %>% group_by(Scenario)%>% mutate(MTD=(Rate < Pint_BLRM[3] & Rate >= Pint_BLRM[2])) %>% filter(MTD)

plotdt <- bydose %>% full_join(MTDdt, by = c("Scenario", "Dose")) %>% filter(MTD) %>% 
  group_by(Scenario, Design) %>% summarize_at(c("MTDFreq", "Npat"), sum) %>% mutate(Design = factor(Design), Scenario = factor(Scenario)) %>%
  mutate(Design = fct_recode(Design, 
                             "Orig BLRM" = "0",
                             "Design 1" = "1",
                             "Design 2" = "2",
                             "Design 3" = "3",
                             "Design 4" = "4")
         )

#MTD accuracy plots
MTDlist <- list()
for(s in 1:3){
  scenario <- unique(scenariodt$Scenario)[s]
  p <-
  ggplot(plotdt %>% filter(Scenario==scenario), aes(x=Design, y=MTDFreq, fill = Design)) + 
    geom_bar(stat = "identity", position=position_dodge()) + guides(fill=FALSE) +
    geom_text(aes(label=MTDFreq), vjust=-0.3, size=5) +
    scale_y_continuous(limits=c(0,1), breaks = seq(0,1,0.1), name = "MTD Accuracy") + xlab("") +
    scale_fill_grey() + theme_bw() + ggtitle(scenario) +
    theme(axis.text = element_text(size = 10),
          axis.title = element_text(size = 12),
          plot.title = element_text(hjust = 0.5, vjust = -8, size = 15)
    )
  MTDlist[[s]] <- p
}

#N treated at MTD
NMTDlist <- list()
for(s in 1:3){
  scenario <- unique(scenariodt$Scenario)[s]
  p <-
    ggplot(plotdt %>% filter(Scenario==scenario), aes(x=Design, y=Npat, fill = Design)) + 
    geom_bar(stat = "identity", position=position_dodge()) + guides(fill=FALSE) +
    geom_text(aes(label=Npat), vjust=-0.3, size=5) +
    scale_y_continuous(limits=c(0,15), breaks = seq(0,15,1), name = "Avg Patients Treated at MTD") + xlab("") +
    scale_fill_grey() + theme_bw() + ggtitle(scenario) +
    theme(axis.text = element_text(size = 10),
          axis.title = element_text(size = 12),
          plot.title = element_text(hjust = 0.5, vjust = -8, size = 15)
    )
  NMTDlist[[s]] <- p
}

#arrange them into the same figure
png("16_33.png", width = 900, height = 1200, res = 100)
ggarrange(MTDlist[[1]], NMTDlist[[1]], MTDlist[[2]], NMTDlist[[2]], MTDlist[[3]], NMTDlist[[3]],
          nrow = 3, ncol = 2)
dev.off()

