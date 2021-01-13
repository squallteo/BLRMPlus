package2load <- c("R2jags", "tidyverse", "readxl", "plyr", "foreach")
lapply(package2load, require, character.only = TRUE)

source("00-JAGSModel.R")
source("00-DEFunctions.R")

DoseProv <- c(5, 10, 25, 50, 100, 200, 400)
DoseRef <- 50
Pint_BLRM <- c(0, 0.2, 0.3, 1)
Pint_kbd <- seq(0, 1, 0.1)
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
                  Prior = prior_ab
)


jags_obj <- jags(model.file = BLRM_orig, data = jags_data, parameters.to.save = c("Pr.Tox"), 
                 n.chains = 3, n.burnin = 10000, n.iter = 50000)

PrTox_mcmc <- as_tibble(jags_obj$BUGSoutput$sims.matrix) %>% select(starts_with("Pr.Tox"))

#compute posterior probabilites of each interval outside BUGS for added programming flexibility
#need to revisit if first and last intervals are discarded
#for BLRM and design 2
tt <- interval_prob(Pint_BLRM, DoseProv)
#calculate UPM for design 2
tt / diff(Pint_BLRM) %x% t(rep(1, length(DoseProv)))

#for design 3
interval_prob(Pint_kbd, DoseProv)



