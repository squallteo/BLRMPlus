package2load <- c("R2jags", "tidyverse", "readxl", "plyr", "foreach")
lapply(package2load, require, character.only = TRUE)

source("00-JAGSModel.R")

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


toxcat <- function(row, Pint, DoseProv){
  Ncat <- length(Pint) - 1
  probdiff <- as.matrix(row) %x% rep(1,Ncat) - t(rep(1, length(DoseProv))) %x% Pint[-length(Pint)]
  tt <- (probdiff > 0)
  apply(tt, 2, function(x) max(which(x==1)))
}

row <- PrTox_mcmc[1,]
Pint <- Pint_kbd


out <- as_tibble(adply(PrTox_mcmc, 1, toxcat, Pint = Pint, DoseProv = DoseProv)) %>% select(starts_with("V")) 

out_BLRM <- as_tibble(adply(PrTox_mcmc, 1, toxcat, Pint = Pint_BLRM, DoseProv = DoseProv)) %>% select(starts_with("V")) 
out_kbd <- as_tibble(adply(PrTox_mcmc, 1, toxcat, Pint = Pint_kbd, DoseProv = DoseProv)) %>% select(starts_with("V")) 


do.call(cbind, lapply(out_kbd, function(x) {prop.table(table(x, useNA = "no"))}))

interval_prob <- do.call(cbind, lapply(out_BLRM, function(x) {prop.table(table(x, useNA = "no"))}))
diff(Pint_BLRM)
diff(Pint_kbd)
