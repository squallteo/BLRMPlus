rm(list=ls())

package2load <- c("R2jags", "tidyverse", "doParallel")
lapply(package2load, require, character.only = TRUE)
detach(package:plyr) #Need this to get group_by work properly, due to name collision. We want the group_by in dplyr, not plyr. 

source("00-BLRM.bugs")
source("00-DEFunctions.R")

DoseProv <- c(10, 25, 50, 100, 200, 400, 800)
DoseRef <- 100

# Pint_BLRM <- c(0, 0.16, 0.33, 1); target_prob <- 0.5
Pint_BLRM <- c(0, 0.2, 0.3, 1); target_prob <- 0.3
Nmax <- 45
ewoc <- 0.3
cohort_size <- 3

#Prior distributions
prior_ab <- c(-0.693, 0, 2, 1, 0)

#dose escalation design
#0: original BLRM with EWOC; 
#1: Babb et al like design, based on toxicity interval
#2: relative dose strength design, based on toxicity interval
#3: Babb et al like design, based on UPM
#4: relative dose strength design, based on UPM
design <- 0:4

#Clertant class tox scenarios and house-keeping
toxdt <- read_csv("ClertantClass.csv") %>% group_by(Sim) %>% mutate(AllToxic=all(Rate>=Pint_BLRM[3]), AllUnder=all(Rate<Pint_BLRM[2]),
                                                                    Target=(Rate >= Pint_BLRM[2] & Rate < Pint_BLRM[3]))
#part 1: has MTD
part1 <- toxdt %>% filter(Target) %>% rename(TrueMTD=Dose) %>% select(Sim, TrueMTD)
#part 2: all toxic, MTD = -1
part2 <- toxdt %>% filter(AllToxic & Dose==DoseProv[1]) %>% mutate(TrueMTD=-1) %>% select(Sim, TrueMTD)
#part 3: all under, MTD = 999999
part3 <- toxdt %>% filter(AllUnder & Dose==DoseProv[1]) %>% mutate(TrueMTD=999999) %>% select(Sim, TrueMTD)
#in all other scenarios, neither all toxic nor under, just no MTD in the interval: MTD=999999
MTDTrue <- rbind(part1, part2, part3) %>% ungroup() %>% full_join(tibble(Sim=1:length(unique(toxdt$Sim)))) %>% 
  mutate(TrueMTD=ifelse(is.na(TrueMTD), 999999, TrueMTD)) %>% arrange(Sim, TrueMTD)

nsim <- length(unique(toxdt$Sim))
# nsim <- 10

ncores <- min(parallel::detectCores()-1, 40)
cl <- makeCluster(ncores)
registerDoParallel(cl)

######################################################
######################################################
######################################################
for(d in design){
  resultdt <-
    foreach(s = 1:nsim, .packages = c("R2jags", "tidyverse", "plyr"), .combine = rbind, .errorhandling = "remove") %dopar% {
      set.seed(s+712)
      Dose_curr <- DoseProv[1]
      cumdt <- tibble(Dose = DoseProv, Ntox = 0, Npat = 0, Current = 0)
      stopcode <- 0
      cohortdt_s <- tibble(Dose = 0, Ntox = 0, Npat = 0, cohort = 0)
      cohort_index <- 1
      MTD <- 0
      
      while(stopcode==0){
        toxrate <- toxdt %>% filter(Dose==Dose_curr & Sim==s)
        Ntox_new <- rbinom(n = 1, size = cohort_size, prob = toxrate$Rate)
        
        newdt <-tibble(Dose = Dose_curr, Ntox_new, Npat_new = cohort_size)
        cohortdt_s <- rbind(cohortdt_s, 
                            tibble(Dose = Dose_curr, Ntox = Ntox_new, Npat = cohort_size, cohort = cohort_index))
        cumdt <-
          cumdt %>% mutate(Current = 0) %>% left_join(newdt, by = c("Dose")) %>% 
          mutate(Ntox = Ntox + ifelse(is.na(Ntox_new), 0, Ntox_new), 
                 Npat = Npat + ifelse(is.na(Npat_new), 0, Npat_new), 
                 Current = 1*!is.na(Npat_new)) %>% arrange(Dose) %>% select(-c("Ntox_new", "Npat_new"))
        admdt <- cumdt %>% filter(Npat > 0)
        
        jags_data <- list(Ntox = admdt$Ntox, Npat = admdt$Npat, DoseAdm = admdt$Dose, NdoseAdm = nrow(admdt),
                          DoseProv = cumdt$Dose, NdoseProv = nrow(cumdt), DoseRef = DoseRef,
                          Prior = prior_ab
        )
        jags_obj <- jags(model.file = BLRM_orig, data = jags_data, parameters.to.save = c("Pr.Tox"), 
                         n.chains = 3, n.burnin = 10000, n.iter = 50000, progress.bar = "none")
        PrTox_mcmc <- as_tibble(jags_obj$BUGSoutput$sims.matrix) %>% select(starts_with("Pr.Tox"))
        
        BLRM_prob <- cbind(cumdt, interval_prob(PrTox_mcmc, Pint_BLRM, DoseProv))
        stopcode <- checkstop_TI(BLRM_prob, target.prob = target_prob, max.subj = Nmax, ewoc = ewoc)
        cohort_index <- cohort_index + 1
        
        #continue to next cohort
        if(stopcode == 0){
          if(d==0){
            action <- action_BLRM(BLRM_prob, ewoc)
          }
          if(d==1){
            action <- action_d1(BLRM_prob, ewoc, f_bnd = 0.25)
          }
          if(d==2){
            action <- action_d2(BLRM_prob, ewoc)
          }
          if(d==3){
            action <- action_d3(BLRM_prob, Pint_BLRM, ewoc, f_bnd = 0.25)
          }
          if(d==4){
            action <- action_d4(BLRM_prob, Pint_BLRM, ewoc)
          }
          
          Dose_next <- DoseProv[which(DoseProv == Dose_curr) + action]
          cumdt <- cumdt %>% mutate(Current = (Dose==Dose_next)*1)
          Dose_curr <- as.numeric(cumdt %>% filter(Current==1) %>% select(Dose))
        }
        #stop for reaching maximum sample size
        if(stopcode == 1){
          MTD <- 999999
        }
        #stop for declaring MTD
        if(stopcode == 2){
          MTD <- as.numeric(cumdt %>% filter(Current==1) %>% select(Dose))
        }
        #stop because all doses are toxic
        if(stopcode == 3){
          MTD <- -1
        }
      } #end of while loop
      
      cohortdt_s <- cohortdt_s %>% filter(Dose>0) %>% mutate(Sim=s, MTD=MTD)
      
    } #end of parallel computing
  
  #output R workspace and spreadsheet
  filename <- paste("Clertant_design", d, sep="")
  save.image(file=paste(filename,".RData",sep=""))
  #################################
  #################################
  #################################
  #summarize results
  #MTD accuracy
  MTDResult <- resultdt %>% filter(cohort==1) %>% select(Sim, MTD)
  MTDHit <- MTDTrue %>% filter(Sim <= nsim) %>% left_join(MTDResult, by="Sim") %>% mutate(Hit=(MTD==TrueMTD)) %>% filter(Hit)
  
  MTDaccuracy <- length(unique(MTDHit$Sim))/nrow(MTDResult)
  
  outputdt <- tibble(MTDaccuracy=MTDaccuracy, FinishedSim=nrow(MTDResult))
  
  #output spreadsheet
  write_csv(outputdt, paste(filename,".csv",sep=""))
  
}

stopCluster(cl)