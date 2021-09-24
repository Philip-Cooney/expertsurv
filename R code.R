
library(devtools)

library(testthat)
testthat::with_mock(
  check_not_nested = function(path, name) return(),
  usethis::create_package('expertsurv'),
  .env = "usethis"
)
install.packages("here")
create_package("C:/Users/phili/OneDrive/PhD/R packages/expertsurv/")
here::dr_here("C:/Users/phili/OneDrive/PhD/R packages/")

library(survHE)

fs1 <- flexsurvreg(Surv(recyrs, censrec) ~ group, data = bc,
                   dist = "weibull")


mods <- c("exp", "weibull", "gamma", "lnorm", "llogis", "gengamma")
formula <- Surv(recyrs, censrec) ~ as.factor(group)
m1 <- fit.models(formula = formula, data = bc, distr = mods)

View(m1)
print(m1, mod = 5)

plot(m1, add.km = TRUE)+
  scale_x_continuous(expand = c(0, 0), limits = c(0,NA)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA))

model.fit.plot(m1)

#Change to WAIC


psa <- make.surv(fit = m1, nsim = 1000, t = seq(.1, 63))

psa.plot(psa, xlab = "Extrapolated time",
         ylab = "Estimation of the survival curves",
         alpha = 0.2, col = c("dark grey", "black", "green"),
         main = "PSA to survival curves", cex.txt = .95)