library(ggplot2)
library(tidyverse)
# set.seed(712)

source("00-DEFunctions.R")

dosevec <- c(10, 25, 50, 100, 200, 400, 800)
target_prob <- 0.25
target_interval <- c(0.2, 0.3)
n_display <- 15

#############

toxdt <- ToxClassGen(dosevec, target_prob, nscenarios = 1000, muvec = c(0.2, 0.2), sigmavec <- c(0.2, 0.3, 0.4))

if(n_display > 0){
  toxdt %>% filter(Sim <= n_display) %>%
    ggplot(aes(x=Dose, y=Rate, color=factor(Sim))) + geom_line() + geom_point() + 
    scale_y_continuous(breaks = seq(0, 1,by = 0.1), name = "DLT Rate") +
    geom_hline(yintercept = target_prob) + 
    theme(legend.position = "none")
}

MTDdt <- 
  toxdt %>% mutate(Target = (Rate >= target_interval[1] & Rate < target_interval[2])) %>% 
  group_by(Sim) %>% summarize(NMTD = sum(Target))
table(MTDdt$NMTD)

RateSummdt <- 
  toxdt %>% group_by(Dose) %>% summarize(Lower = quantile(Rate, 0.025), AvgRate = mean(Rate), Upper = quantile(Rate, 0.975)) 

toxdt %>% ggplot(aes(x=Dose, y=Rate, group=Dose)) + geom_boxplot(outlier.shape = NA) + 
  scale_y_continuous(breaks = seq(0, 1,by = 0.1), name = "DLT Rate") + geom_hline(yintercept = target_prob)
