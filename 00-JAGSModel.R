BLRM_orig <- function(){
  cova[1,1] <- Prior[3]*Prior[3]
  cova[2,2] <- Prior[4]*Prior[4]
  cova[1,2] <- Prior[3]*Prior[4]*Prior[5]
  cova[2,1] <- cova[1,2]
  prec[1:2,1:2] <- inverse(cova[,])
  log.alphabeta[1:2] ~ dmnorm(Prior[1:2],prec[1:2,1:2])
  
  #likelihood
  for (j in 1:NdoseAdm){
    logit(Pr.Tox1[j]) <- log.alphabeta[1] + exp(log.alphabeta[2]) * log(DoseAdm[j]/DoseRef)
    Ntox[j] ~ dbin(Pr.Tox1[j],Npat[j])
  }
  
  # for each dose: probabilities of toxicity, category probabilities, risks
  for (i in 1:NdoseProv) {
    lin[i] <- log.alphabeta[1] + exp(log.alphabeta[2]) * log(DoseProv[i]/DoseRef) 
    logit(Pr.Tox[i]) <- lin[i]
    Pr.Cat[i,1] <- step(Pint[1]-Pr.Tox[i])
    Pr.Cat[i,2] <- step(Pint[2]-Pr.Tox[i])-step(Pint[1]-Pr.Tox[i])
    Pr.Cat[i,3] <- step(1-Pr.Tox[i])-step(Pint[2]-Pr.Tox[i])
  }
}

