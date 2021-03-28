library(tidyverse)
library(ggplot2)

DoseProv <- c(10, 25, 50, 100, 200, 400, 800)
DoseRef <- 100

expit <- function(x) 1/(1+exp(-x))
loga <- -0.916
b <- 1.2

ratedt <- tibble(Dose = DoseProv, Rate = round(expit(loga + b*log(DoseProv/DoseRef)),3))
ggplot(ratedt,aes(Dose, Rate, label=Rate)) + geom_line() + geom_point() + geom_text(nudge_x = 2)

ratedt <- readxl::read_xlsx("FixedScenarios.xlsx", sheet = "Sheet1")
ggplot(ratedt,aes(Dose, Rate, color = Scenario)) + geom_line() + geom_point()

# write_csv(ratedt, "ToxScenarios.csv")
