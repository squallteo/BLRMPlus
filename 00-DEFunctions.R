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

#functions for checking whether the trial should stop
checkstop_TI = function(probdt, target.prob = 0.5, min.subj.MTD = 6, max.subj = 40, ewoc = 0.25){
  probdt <- probdt %>% mutate(Toxic = (Pover > ewoc)*1)
  Ntotal <- probdt %>% summarize_at("Npat",sum) %>% select(Npat)
  currdt <- probdt %>% filter(Current == 1)
  
  max_target <- probdt %>% filter(Npat > 0) %>% summarize_at("Ptarget", max)
  #target interval probability achieved, recommendation is to stay & EWOC okay 
  cond1 <- (currdt$Ptarget >= target.prob & !currdt$Toxic & currdt$Ptarget == max_target)
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

#functions for making dose recommendation
#original BLRM
action_BLRM = function(probdt, ewoc = 0.25){
  next_dose <- probdt %>% filter(Pover < ewoc) %>% filter(Ptarget == max(Ptarget)) %>% select(Dose)
  curr_dose <- probdt %>% filter(Current==1) %>% select(Dose)
  action <- -1*(next_dose < curr_dose) + 0*(next_dose == curr_dose) + 1*(next_dose > curr_dose)
  
  return(action)
}

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

action_d3 = function(probdt, Pint, ewoc = 0.25, f_bnd){
  upmdt <- probdt %>% mutate(UPMunder = Punder/diff(Pint)[1],
                             UPMtarget = Ptarget/diff(Pint)[2],
                             UPMover = Pover/diff(Pint)[3])
  curr_dose_idx <- which(probdt$Current==1)
  #rule to override EWOC
  ovrd <- 0
  if(curr_dose_idx != nrow(probdt)){
    
    loss_under <- upmdt$UPMunder[curr_dose_idx] * f_bnd
    loss_over <- upmdt$UPMover[curr_dose_idx] * (1 - f_bnd)
    ovrd <- (loss_under > loss_over)
    action <- 1
  }
  
  if(ovrd==0){
    next_dose <- upmdt %>% filter(Pover < ewoc) %>% filter(Ptarget == max(Ptarget)) %>% select(Dose)
    curr_dose <- upmdt %>% filter(Current==1) %>% select(Dose)
    action <- -1*(next_dose < curr_dose) + 0*(next_dose == curr_dose) + 1*(next_dose > curr_dose)
  }
  
  return(action)
  
}

action_d4 = function(probdt, Pint, ewoc = 0.25){
  upmdt <- probdt %>% mutate(UPMunder = Punder/diff(Pint)[1],
                             UPMtarget = Ptarget/diff(Pint)[2],
                             UPMover = Pover/diff(Pint)[3])
  curr_dose_idx <- which(probdt$Current==1)
  #rule to override EWOC
  ovrd <- 0
  if(curr_dose_idx != nrow(probdt)){
    relstr <- upmdt$Dose[curr_dose_idx] / upmdt$Dose[curr_dose_idx+1]
    cond1 <- (upmdt$UPMunder[curr_dose_idx] * relstr > upmdt$UPMover[curr_dose_idx+1])
    cond2 <- (upmdt$Pover[curr_dose_idx+1] > ewoc)
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
