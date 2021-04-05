library(tidyverse)
library(ggplot2)

DoseProv <- c(10, 25, 50, 100, 200, 400, 800)
DoseRef <- 100

expit <- function(x) 1/(1+exp(-x))

d <- 1:800
dref <- DoseRef

#steep scenario
loga <- -0.916
b <- 1.2
linpred <- loga + b*log(d/dref)
rate <- expit(linpred)

steepdt <- tibble(Scenario="1: Steep", Dose=d, DLTRate=rate)

#S-shaped scenario
s <- 0.02
ed50 <- 225
linpred <- s*(d-ed50)
rate <- 0.6*expit(linpred)

sshapedt <- tibble(Scenario="2: S-Shaped", Dose=d, DLTRate=rate)

#flat scenario
s <- 2
ed50 <- 700
linpred <- -s*log(ed50/d)
rate <- expit(linpred)

flatdt <- tibble(Scenario="3: Flat", Dose=d, DLTRate=rate)

plotdt <- rbind(steepdt, sshapedt, flatdt)
pointdt <- plotdt %>% filter(Dose %in% DoseProv)

ggplot(plotdt, aes(x=Dose, y=DLTRate, linetype=Scenario)) + geom_line(size=1.5) + 
  geom_point(aes(x=Dose, y=DLTRate), data = pointdt, size=4) + ylab("DLT Rate") +
  scale_x_continuous(breaks = DoseProv) + scale_y_continuous(breaks = seq(0,1,0.1)) + theme_bw() +
  theme(axis.text.x = element_text(angle = 90),
        legend.position = c(0.9,0.2))
  




ratedt <- tibble(Dose = DoseProv, Rate = round(expit(),3))
ggplot(ratedt,aes(Dose, Rate, label=Rate)) + geom_line() + geom_point() + geom_text(nudge_x = 2)

ratedt <- readxl::read_xlsx("FixedScenarios.xlsx", sheet = "Sheet1")
ggplot(ratedt,aes(Dose, Rate, color = Scenario)) + geom_line() + geom_point()

# write_csv(ratedt, "ToxScenarios.csv")
