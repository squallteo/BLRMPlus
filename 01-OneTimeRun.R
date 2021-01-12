library(R2jags)
library(tidyverse)
library(readxl)

source("00-JAGSModel.R")

DoseProv_mono <- c(5, 10, 25, 50, 100, 200, 400)
DoseRef <- 50
Pint <- c(0.16, 0.33)
ewoc <- 0.3

#Prior distributions
prior_ab <- c(-0.693, 0, 2, 1, 0)

#cumulative and new data sets
cumdt <- read_excel("DevData.xlsx", sheet = "cumdt")
newdt <- read_excel("DevData.xlsx", sheet = "newdt")

updatedt <- 
  cumdt %>% mutate(Current = 0) %>% left_join(newdt, by = c("Dose")) %>% 
  mutate(Ntox = Ntox + ifelse(is.na(Ntox_new), 0, Ntox_new), 
         Npat = Npat + ifelse(is.na(Npat_new), 0, Npat_new), 
         Current = 1*!is.na(Npat_new)) %>% arrange(Dose) %>% select(-c("Ntox_new", "Npat_new"))

admdt <- updatedt %>% filter(Npat > 0)

jags_data <- list(Ntox = admdt$Ntox, Npat = admdt$Npat, DoseAdm = admdt$Dose, NdoseAdm = nrow(admdt),
                  DoseProv = updatedt$Dose, NdoseProv = nrow(updatedt), DoseRef = DoseRef,
                  Prior = prior_ab, Pint = Pint
)

# param <- c("Pr.Cat")
# jags_obj <- jags(model.file = BLRM_orig, data = jags_data, parameters.to.save = param, 
#                  n.chains = 3, n.burnin = 10000, n.iter = 50000)

jags_obj <- jags(model.file = BLRM_orig, data = jags_data, parameters.to.save = "Pr.Tox", 
                 n.chains = 3, n.burnin = 10000, n.iter = 50000, n.thin = 5)

PrTox <- as_tibble(jags_obj$BUGSoutput$sims.matrix) %>% select(-deviance)

Pint <- c(0, 0.16, 0.33, 1)


# jags_auto <- autojags(jags_obj, Rhat = 1.1, n.thin = 4, n.iter = 40000, n.update = 3)
