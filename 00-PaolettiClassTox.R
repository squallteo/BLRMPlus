library(ggplot2)
library(tidyverse)
detach(package:plyr) #Need this to get group_by work properly, due to name collision. We want the group_by in dplyr, not plyr. 
set.seed(113)

source("00-DEFunctions.R")

dosevec <- c(10, 25, 50, 100, 200, 400, 800)
target_prob <- 0.25
Pint_BLRM <- c(0, 0.2, 0.3, 1)
n_display <- 3:20

#############

toxdt <- ToxClassGen(dosevec, target_prob, nscenarios = 1000, muvec = c(0.2, 0.2), sigmavec <- c(0.1, 0.3, 0.4))

if(n_display != 0){
  toxdt %>% filter(Sim %in% n_display) %>%
    ggplot(aes(x=Dose, y=Rate, color=factor(Sim))) + geom_line() + geom_point() + 
    scale_y_continuous(breaks = seq(0, 1,by = 0.1), name = "DLT Rate") +
    geom_hline(yintercept = target_prob) + 
    theme(legend.position = "none")
}

MTDdt <- 
  toxdt %>% mutate(Target = (Rate >= Pint_BLRM[2] & Rate < Pint_BLRM[3])) %>% 
  dplyr::group_by(Sim) %>% summarize(NMTD = sum(Target))
table(MTDdt$NMTD)

RateSummdt <- 
  toxdt %>% group_by(Dose) %>% summarize(Lower = quantile(Rate, 0.025), AvgRate = mean(Rate), Upper = quantile(Rate, 0.975)) 

toxdt %>% ggplot(aes(x=Dose, y=Rate, group=Dose)) + geom_boxplot(outlier.shape = NA) + 
  scale_y_continuous(breaks = seq(0, 1,by = 0.1), name = "DLT Rate") + geom_hline(yintercept = Pint_BLRM[2:3])

write_csv(toxdt, "PaolettiClass.csv")
