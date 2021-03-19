package2load <- c("R2jags", "tidyverse", "readxl", "plyr", "doParallel")
lapply(package2load, require, character.only = TRUE)

source("00-JAGSModel.R")
source("00-DEFunctions.R")

nsim <- 1000

tt <- read_csv("ToxScenarios.csv", col_names = T)
for(i in 1:nsim){
  if(i==1) toxdt <- tibble(tt, Sim = i)
  else toxdt <- rbind(toxdt, tibble(tt, Sim = i))
}

DoseProv <- c(10, 25, 50, 100, 200, 400, 800)
DoseRef <- 100

Pint_BLRM <- c(0, 0.2, 0.3, 1)
Nmax <- 45
target_prob <- 0.4
ewoc <- 0.3
cohort_size <- 3

#Prior distributions
prior_ab <- c(-0.693, 0, 2, 1, 0)

ncores <- min(parallel::detectCores()-1, 40)
cl <- makeCluster(ncores)
registerDoParallel(cl)