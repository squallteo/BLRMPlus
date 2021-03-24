ToxClassGen <- function(dosevec, target_prob, nscenarios, muvec, sigmavec){
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

toxcat <- function(row, Pint, DoseProv){
  Ncat <- length(Pint) - 1
  probdiff <- as.matrix(row) %x% rep(1,Ncat) - t(rep(1, length(DoseProv))) %x% Pint[-length(Pint)]
  tt <- (probdiff > 0)
  apply(tt, 2, function(x) max(which(x==1)))
}

interval_prob <- function(jagsdata, Pint, DoseProv){
  interval_levels <- 1:(length(Pint)-1)
  out <- as_tibble(adply(jagsdata, 1, toxcat, Pint = Pint, DoseProv = DoseProv)) %>% select(starts_with("V"))
  
  #crucial step to convert all columns to factors, in case some categories are not present at a certain dose level
  out_fct = lapply(out, factor, levels = interval_levels)
  tt <- as_tibble(do.call(cbind, lapply(out_fct, function(x) {prop.table(table(x, useNA = "no"))})))
  probdt <- as_tibble(t(tt), .name_repair = "minimal")
  names(probdt) <- c("Punder", "Ptarget", "Pover")
  
  return(probdt)
}

#BLRM functions
checkstop_BLRM = function(probdt, target.prob = 0.5, min.subj.MTD = 6, max.subj = 40, ewoc = 0.25){
  probdt <- probdt %>% mutate(Toxic = (Pover > ewoc)*1)
  Ntotal <- probdt %>% summarize_at("Npat",sum) %>% select(Npat)
  currdt <- probdt %>% filter(Current == 1)
  
  #target interval probability achieved & EWOC okay
  cond1 <- (currdt$Ptarget >= target.prob & !currdt$Toxic)
  #minimum number of subjects achieved at MTD
  cond2 <- (currdt$Npat >= min.subj.MTD)
  #maximum number of subjects per strata reached
  cond3 <- c(Ntotal >= max.subj)
  #all doses are toxic
  cond4 <- probdt %>% summarize_at("Toxic", prod) %>% select(Toxic)
  
  stop4mtd <- (cond1 & cond2)
  stop4tox <- cond4
  stop4ss <- !(stop4mtd | stop4tox) & cond3
  nostop <- !(stop4mtd | stop4ss | stop4tox)
  stopcode <- 0 * nostop + 1 * stop4ss + 2 * stop4mtd + 3 * stop4tox
  
  return(stopcode)
}

action_BLRM = function(probdt, ewoc = 0.25){
  next_dose <- probdt %>% filter(Pover < ewoc) %>% filter(Ptarget == max(Ptarget)) %>% select(Dose)
  curr_dose <- probdt %>% filter(Current==1) %>% select(Dose)
  action <- -1*(next_dose < curr_dose) + 0*(next_dose == curr_dose) + 1*(next_dose > curr_dose)
  
  return(action)
}

#New design 1 functions
#f_bnd: feasibility bound in Babb et al. (1998)
action_d1 = function(probdt, ewoc = 0.25, f_bnd = 0.25){
  curr_dose_idx <- which(probdt$Current==1)
  #rule to override EWOC
  ovrd <- 0
  if(curr_dose_idx != nrow(probdt)){
    
    loss_under <- probdt$Punder[curr_dose_idx] * f_bnd
    loss_over <- probdt$Pover[curr_dose_idx] * (1 - f_bnd)
    ovrd <- (loss_under > loss_over)
    action <- 1
  }
  
  if(ovrd==0){
    next_dose <- probdt %>% filter(Pover < ewoc) %>% filter(Ptarget == max(Ptarget)) %>% select(Dose)
    curr_dose <- probdt %>% filter(Current==1) %>% select(Dose)
    action <- -1*(next_dose < curr_dose) + 0*(next_dose == curr_dose) + 1*(next_dose > curr_dose)
  }
  
  return(action)
}

#New design 2 functions
action_d2 = function(probdt, ewoc = 0.25){
  curr_dose_idx <- which(probdt$Current==1)
  #rule to override EWOC
  ovrd <- 0
  if(curr_dose_idx != nrow(probdt)){
    relstr <- probdt$Dose[curr_dose_idx] / probdt$Dose[curr_dose_idx+1]
    cond1 <- (probdt$Punder[curr_dose_idx] * relstr > probdt$Pover[curr_dose_idx+1])
    cond2 <- (probdt$Pover[curr_dose_idx+1] > ewoc)
    ovrd <- cond1*cond2
    action <- 1
  }
  
  if(ovrd==0){
    next_dose <- probdt %>% filter(Pover < ewoc) %>% filter(Ptarget == max(Ptarget)) %>% select(Dose)
    curr_dose <- probdt %>% filter(Current==1) %>% select(Dose)
    action <- -1*(next_dose < curr_dose) + 0*(next_dose == curr_dose) + 1*(next_dose > curr_dose)
  }
  
  return(action)
}

#New design 3 functions
action_d3 = function(probdt, Pint, ewoc = 0.25){
  browser()
  upmdt <- probdt / diff(Pint_BLRM) %x% t(rep(1, length(DoseProv)))
  
  next_dose <- probdt %>% filter(Pover < ewoc) %>% filter(Ptarget == max(Ptarget)) %>% select(Dose)
  curr_dose <- probdt %>% filter(Current==1) %>% select(Dose)
  action <- -1*(next_dose < curr_dose) + 0*(next_dose == curr_dose) + 1*(next_dose > curr_dose)
  
  return(action)
}

#calculate UPM for design 2
