library(tidyverse)

DoseProv <- c(10, 25, 50, 100, 200, 400, 800)
DoseRef <- 100

expit <- function(x) 1/(1+exp(-x))
a <- 0.4
b <- 1.2

ratedt <- tibble(Dose = DoseProv, Rate = expit(log(a) + b*log(DoseProv/DoseRef))) %>% View()

write_csv(ratedt, "ToxScenarios.csv")
