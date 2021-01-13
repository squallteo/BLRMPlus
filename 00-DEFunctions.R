toxcat <- function(row, Pint, DoseProv){
  Ncat <- length(Pint) - 1
  probdiff <- as.matrix(row) %x% rep(1,Ncat) - t(rep(1, length(DoseProv))) %x% Pint[-length(Pint)]
  tt <- (probdiff > 0)
  apply(tt, 2, function(x) max(which(x==1)))
}

interval_prob <- function(Pint, DoseProv){
  interval_levels <- 1:(length(Pint)-1)
  out <- as_tibble(adply(PrTox_mcmc, 1, toxcat, Pint = Pint, DoseProv = DoseProv)) %>% select(starts_with("V"))
  
  #crucial step to convert all columns to factors, in case some categories are not present at a certain dose level
  out_fct = lapply(out, factor, levels = interval_levels)
  do.call(cbind, lapply(out_fct, function(x) {prop.table(table(x, useNA = "no"))}))
}