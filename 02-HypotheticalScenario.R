rm(list=ls())

package2load <- c("R2jags", "tidyverse", "plyr")
lapply(package2load, require, character.only = TRUE)

source("00-BLRM.bugs")
source("00-DEFunctions.R")

DoseProv <- c(10, 25, 50, 100, 200, 400, 800)
DoseRef <- 100

Pint_BLRM <- c(0, 0.16, 0.33, 1)
ewoc <- 0.25

prior_ab <- c(-0.693, 0, 2, 1, 0)

cumdt <- tibble(Dose=DoseProv,
                Ntox=rep(0,7), Npat=c(3,3,3,3,0,0,0),
                Current=c(0,0,0,1,0,0,0))
admdt <- cumdt %>% filter(Npat>0)

jags_data <- list(Ntox = admdt$Ntox, Npat = admdt$Npat, DoseAdm = admdt$Dose, NdoseAdm = nrow(admdt),
                  DoseProv = cumdt$Dose, NdoseProv = nrow(cumdt), DoseRef = DoseRef,
                  Prior = prior_ab
)

jags_obj <- jags(model.file = BLRM_orig, data = jags_data, parameters.to.save = c("Pr.Tox"), 
                 n.chains = 3, n.burnin = 10000, n.iter = 50000, progress.bar = "none")
PrTox_mcmc <- as_tibble(jags_obj$BUGSoutput$sims.matrix) %>% select(starts_with("Pr.Tox"))

BLRM_prob <- cbind(cumdt, interval_prob(PrTox_mcmc, Pint_BLRM, DoseProv))

# Dose Ntox Npat Current    Punder     Ptarget        Pover
# 1   10    0    3       0 0.9963333 0.003333333 0.0003333333
# 2   25    0    3       0 0.9893333 0.010333333 0.0003333333
# 3   50    0    3       0 0.9660000 0.031666667 0.0023333333
# 4  100    0    3       1 0.7993333 0.153000000 0.0476666667
# 5  200    0    0       0 0.5036667 0.200000000 0.2963333333
# 6  400    0    0       0 0.3470000 0.174333333 0.4786666667
# 7  800    0    0       0 0.2516667 0.150666667 0.5976666667
