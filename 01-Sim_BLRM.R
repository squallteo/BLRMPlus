package2load <- c("R2jags", "tidyverse", "readxl", "plyr", "foreach")
lapply(package2load, require, character.only = TRUE)

source("00-JAGSModel.R")
source("00-DEFunctions.R")

set.seed(712)

toxdt <- read_excel("DevData.xlsx", sheet = "toxdt") #toxicity scenario

DoseProv <- toxdt$Dose
DoseRef <- 50
Nmax <- 45
Pint_BLRM <- c(0, 0.2, 0.3, 1)
target_prob <- 0.4
ewoc <- 0.3
cohort_size <- 3

#Prior distributions
prior_ab <- c(-0.693, 0, 2, 1, 0)

nsim <- 2
######################################################
######################################################
######################################################
MTDvec <- rep(0, nsim)

for(s in 1:nsim){
  Dose_curr <- DoseProv[1]
  cumdt <- tibble(Dose = DoseProv, Ntox = 0, Npat = 0, Current = 0)
  stopcode <- 0
  cohortdt_s <- tibble(Dose = 0, Ntox = 0, Npat = 0, cohort = 0)
  cohort_index <- 1
  
  while(stopcode==0){
    toxrate <- toxdt %>% filter(Dose==Dose_curr)
    Ntox_new <- rbinom(n = 1, size = cohort_size, prob = toxrate$prob)
    
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
    stopcode <- checkstop_BLRM(BLRM_prob, target.prob = target_prob, max.subj = Nmax, ewoc = ewoc)
    cohort_index <- cohort_index + 1
    
    #continue to next cohort
    if(stopcode == 0){
      action <- action_BLRM(BLRM_prob, ewoc)
      Dose_next <- DoseProv[which(DoseProv == Dose_curr) + action]
      cumdt <- cumdt %>% mutate(Current = (Dose==Dose_next)*1)
      Dose_curr <- as.numeric(cumdt %>% filter(Current==1) %>% select(Dose))
    }
    #stop for reaching maximum sample size
    if(stopcode == 1){
      MTDvec[s] <- NA
    }
    #stop for declaring MTD
    if(stopcode == 2){
      MTDvec[s] <- as.numeric(cumdt %>% filter(Current==1) %>% select(Dose))
    }
    #stop because all doses are toxic
    if(stopcode == 3){
      MTDvec[s] <- -1
    }
  }
  
  if(s==1){
    simdt <- cumdt %>% mutate(sim = s)
    cohortdt <- cohortdt_s %>% mutate(sim = s)
  }
  else{
    simdt <- rbind(simdt, cumdt %>% mutate(sim = s))
    cohortdt <- rbind(cohortdt, cohortdt_s %>% mutate(sim = s))
  } 
  
}
