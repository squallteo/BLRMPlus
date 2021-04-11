ClertantClassGen <- function(dosevec, target_prob, nscenarios){
  J <- length(dosevec)
  for(s in 1:nscenarios){
    j <- sample(1:J, 1)
    M <- rbeta(1, max(J-j, 0.5), 1)
    B <- target_prob + (1 - target_prob)*M
    
    MTD_index <- 0
    while(j!=MTD_index){
      toxrate <- sort(runif(J, 0, B))
      dev <- abs(toxrate - target_prob)
      MTD_index <- which(dev==min(dev))
    }
    
    tt <- tibble(Sim = s, Dose = dosevec, Rate = toxrate)
    
    if(s==1){toxdt <- rbind(tt)}
    else{toxdt <- rbind(toxdt, tt)}
  }
  return(toxdt)
}

#########################################################
#########################################################
#########################################################
library(ggplot2)
library(tidyverse)
set.seed(113)
dosevec <- c(10, 25, 50, 100, 200, 400, 800)
target_prob <- 0.25
Pint_BLRM <- c(0, 0.16, 0.33, 1)
n_display <- 1:20

toxdt <- ClertantClassGen(dosevec, target_prob, nscenario=1000)

if(n_display != 0){
  clertant_plot <- 
  toxdt %>% filter(Sim %in% n_display) %>%
    ggplot(aes(x=Dose, y=Rate, color=factor(Sim))) + geom_line() + geom_point() + 
    scale_y_continuous(breaks = seq(0, 1,by = 0.1), name = "DLT Rate") +
    scale_x_continuous(breaks = dosevec) +
    geom_hline(yintercept = target_prob) + theme_bw() +
    theme(legend.position = "none")
  print(clertant_plot)
}

MTDdt <- 
  toxdt %>% mutate(Target = (Rate >= Pint_BLRM[2] & Rate < Pint_BLRM[3])) %>% 
  dplyr::group_by(Sim) %>% summarize(NMTD = sum(Target))
table(MTDdt$NMTD)

RateSummdt <- 
  toxdt %>% group_by(Dose) %>% summarize(Lower = quantile(Rate, 0.025), AvgRate = mean(Rate), Upper = quantile(Rate, 0.975)) 

toxdt %>% ggplot(aes(x=Dose, y=Rate, group=Dose)) + geom_boxplot(outlier.shape = NA) + 
  scale_y_continuous(breaks = seq(0, 1,by = 0.1), name = "DLT Rate") + geom_hline(yintercept = Pint_BLRM[2:3])

# write_csv(toxdt, "ClertantClass.csv")
