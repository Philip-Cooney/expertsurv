# install.packages("survHE", lib ='~/R-packages/')
# install.packages("ggmcmc", lib ='~/R-packages/')

#devtools::install_github("Philip-Cooney/expertsurv", lib =  '~/R-packages/')
#remove.packages("expertsurv", lib = '~/R-packages/')
#remove.packages("survHE", lib = '~/R-packages/')
library("expertsurv", lib.loc = '~/R-packages/')
library("dplyr")

n <- 2
times <- rexp(n, 0.2)
status <- rep(1,n)
time_study <- 200
conflict_prefer("mutate", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")

df1 <- data.frame(time = times, status = status)
data2 <- df1  %>% mutate(time2 = ifelse(time > time_study, time_study, time), status2 = ifelse(time> time_study, 0, status))


plot(survfit(Surv(time2,status2)~1,data=data2))

param_expert_example1 <- list()

param_expert_example1[[1]] <- data.frame(dist = c("norm","t"),
                                         wi = c(1,1),
                                         param1 = c(0.1,0.12),
                                         param2 = c(0.005,0.005),
                                         param3 = c(NA,3))
timepoint_expert <- c(14)

plot_opinion1<- plot_expert_opinion(param_expert_example1[[1]], 
                                    weights = param_expert_example1[[1]]$wi)
#ggsave("Vignette_Example 1 - Expert Opinion.png")

cred_int_val <- cred_int(plot_opinion1,val = "linear pool", interval = c(0.025, 0.975))

undebug(survHE::plot.survHE)
undebug(plot_ggplot_survHE)
example1  <- expertsurv:::fit.models.expert(formula=Surv(time2,status2)~1,
                                            data=data2,
                                            distr=c("wei", "lnorm","llogis", "gompertz",
                                                    "gengamma","exp", "rps", "gamma"),
                                            method="hmc",
                                            pool_type = "linear pool",
                                            iter = 5000,
                                            opinion_type = "survival",
                                            times_expert = timepoint_expert, 
                                            param_expert = param_expert_example1,
                                            a0 = rep(0.001,nrow(df1)))

#undebug(expertsurv:::runHMC)
#debug(expertsurv:::)
undebug(expertsurv:::compute_ICs_stan)


posterior <- as.matrix(example1$models[["log-Normal"]])

library("MASS")
DENS <- kde2d(posterior[,c("alpha")],posterior[,c("meanlog")])
contour(DENS)
filled.contour(DENS,plot.axes = {
  axis(1)
  axis(2)
  contour(DENS,add = TRUE)})

plot(example1, add.km = T, t = 0:30)+
  theme_light()+
  scale_x_continuous(expand = c(0, 0), limits = c(0,NA), breaks=seq(0, 30, 2)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA), breaks=seq(0, 1, 0.05))+
  geom_segment(aes(x = 14, y = cred_int_val[1], xend = 14, yend = cred_int_val[2]))

#Issue with surv HE plot is based on the posterior mean. FINE for MLE however, issue when the parameters have different
#interpretation

psa <- survHE::make.surv(fit = example1, mod = 1, nsim = 5000, t = seq(0, 63))
library("ggmcmc", lib.loc = '~/R-packages/' )
library("coda")
ggmcmc(ggs(example1$models$`log-Normal`), file = "Log-Normal.pdf")
ggmcmc(ggs(example1$models$`Royston-Parmar`), file = "RP.pdf")
ggmcmc(ggs(example1$models$Exponential), file = "Exponential.pdf")
ggmcmc(ggs(as.mcmc(example1$models$`Gen. Gamma`)), file = "Gengamma.pdf")
ggmcmc(ggs(as.mcmc(example1$models$`Gamma`)), file = "Gamma.pdf")
ggmcmc(ggs(as.mcmc(example1$models$`Gompertz`)), file = "Gompertz.pdf")
ggmcmc(ggs(example1$models$`Weibull (AFT)`), file = "Weibull.pdf")

example1$misc$data.stan
time_surv_vec <- seq(1,20, by =0.1)

plot(x = time_surv_vec, y = quant_high, ylim = c(0,1), type = "l")
lines(x = time_surv_vec, y = quant_low)
lines(x = psa$times, y= sim.survHE[1,], col = "red", lty = 2)
lines(x = psa$times, y= sim.survHE[2,], col = "red", lty = 2)

quant_low <- quant_high <- rep(NA, length(time_surv_vec))

for(i in 1:length(time_surv_vec)){
  time_surv <- time_surv_vec[i]
  
  surv_temp<- apply(example1$models$Gompertz$BUGSoutput$sims.matrix[,c("alpha","rate")],1,function(x){
    pgompertz(time_surv,shape =x[1],rate= x[2], lower.tail = F)})
  quant_low[i] <- quantile(surv_temp,c(0.025))
  quant_high[i] <- quantile(surv_temp,c(0.975))
}

sim.survHE<- apply(psa$mat[[1]],1, quantile, probs = c(0.025, 0.975))





head(example1$models$Gompertz$BUGSoutput$sims.matrix)
psa <- survHE::make.surv(fit = example1, mod = 5, nsim = 5000, t = seq(0, 63))

#head(psa[["sim"]][[1]])
survHE::psa.plot(psa)+
  theme_light()+
  scale_x_continuous(expand = c(0, 0), limits = c(0,NA), breaks=seq(0, 50, 2)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA), breaks=seq(0, 1, 0.05))+
  geom_segment(aes(x = 14, y = cred_int_val[1], xend = 14, yend = cred_int_val[2]))



