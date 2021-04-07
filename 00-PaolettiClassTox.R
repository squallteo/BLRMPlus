PaolettiClassGen <- function(dosevec, target_prob, nscenarios, muvec, sigmavec){
  #Generate a class of dose-toxicity scenario originally proposed in Paoletti et al (2004), also elaborated in Liu and Yuan (2015) BOIN paper
  #sigmavec: argument 1 - sd for the deviation from target toxicity rate at MTD; argument 2: sd of deviation to generate DLT rates below MTD; argument 3: sd of deviation to generate DLT rates above MTD
  #muvec: argument 1 - mean of deviation to generate DLT rates below MTD; argument 2 - mean of deviation to generate DLT rates above MTD
  library(tidyverse)
  
  for(s in 1:nscenarios){
    ratevec <- rep(NA, length(dosevec))
    #step 1: randomly select a MTD with equal probabilities
    MTD <- sample(dosevec,1)
    MTD_index <- which(dosevec == MTD)
    
    #step 2: generate the DLT rate at MTD
    ratevec[MTD_index] <- pnorm(rnorm(1, mean = qnorm(target_prob), sd = sigmavec[1]))
    
    #step 3: generate the DLT rate(s) at doses adjacent to MTD, subject to the constraint that DLT rate of MTD is closest to the target probability
    if(MTD_index != 1){
      comp1 <- qnorm(ratevec[MTD_index])
      comp2 <- qnorm(ratevec[MTD_index]) - qnorm(2 * target_prob - ratevec[MTD_index])
      comp3 <- (qnorm(ratevec[MTD_index]) > qnorm(target_prob))
      e1 <- rnorm(1, mean = muvec[1], sd = sigmavec[2])
      ratevec[MTD_index-1] <- pnorm(comp1 - comp2 * comp3 - e1^2)
    }
    if(MTD_index != length(dosevec)){
      comp1 <- qnorm(ratevec[MTD_index])
      comp2 <- qnorm(2 * target_prob - ratevec[MTD_index]) - qnorm(ratevec[MTD_index])
      comp3 <- (qnorm(ratevec[MTD_index]) < qnorm(target_prob))
      e2 <- rnorm(1, mean = muvec[2], sd = sigmavec[3])
      ratevec[MTD_index+1] <- pnorm(comp1 + comp2 * comp3 + e2^2)
    }
    
    #step 4: sequentially generate the DLT rates for remaining dose levels
    if(MTD_index > 2){
      for(j in (MTD_index - 2):1){
        e1 <- rnorm(1, mean = muvec[1], sd = sigmavec[2])
        ratevec[j] <- pnorm(qnorm(ratevec[j+1]) - e1^2)
      }
    }
    if(MTD_index < (length(dosevec)-1)){
      for(j in (MTD_index + 2):length(dosevec)){
        e2 <- rnorm(1, mean = muvec[2], sd = sigmavec[3])
        ratevec[j] <- pnorm(qnorm(ratevec[j-1]) + e2^2)
      }
    }
    
    toxdt <- tibble(Dose = dosevec, Rate = ratevec, Sim = s)
    if(s==1) outdt <- toxdt
    else outdt <- rbind(outdt, toxdt)
  }
  
  return(outdt)
}

########################################################################
########################################################################
########################################################################
########################################################################

library(ggplot2)
library(tidyverse)
detach(package:plyr) #Need this to get group_by work properly, due to name collision. We want the group_by in dplyr, not plyr. 
set.seed(113)

dosevec <- c(10, 25, 50, 100, 200, 400, 800)
target_prob <- 0.25
Pint_BLRM <- c(0, 0.2, 0.3, 1)
n_display <- 3:20

#############

toxdt <- PaolettiClassGen(dosevec, target_prob, nscenarios = 1000, muvec = c(0.2, 0.2), sigmavec <- c(0.1, 0.3, 0.4))

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
