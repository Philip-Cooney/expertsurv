library("expertsurv")

devtools::load_all()
param_expert_example1 <- list()

param_expert_example1[[1]] <- data.frame(dist = c("norm","t"),
                                         wi = c(0.5,0.5), # Ensure Weights sum to 1
                                         param1 = c(0.1,0.12),
                                         param2 = c(0.1,0.1),
                                         param3 = c(NA,3))

plot_opinion1<- plot_expert_opinion(param_expert_example1[[1]], 
                                    weights = param_expert_example1[[1]]$wi, St_indic = 1)

plot_opinion1+
scale_x_continuous( limits = c(0,.6), breaks=seq(0, .6, 0.05)) + 
  scale_y_continuous(limits = c(0, 5), breaks=seq(0, 5, 0.1))
cred_int_val <- cred_int(plot_opinion1,val = "log pool", interval = c(0.025, 0.975))

#View(plot_opinion1[["data"]])
data2 <- expertsurv::data %>% rename(status = censored) %>% mutate(time2 = ifelse(time > 10, 10, time),
                                                               status2 = ifelse(time> 10, 0, status))

#Set the opinion type to "survival"
# undebug(fit.models.expert)
# undebug(fit.models)
# undebug(expertsurv:::runHMC)
# undebug(expertsurv:::make_data_stan)
# undebug(expertsurv:::compute_ICs_stan)
# 
# undebug(roxygen2::roxygenise)
# undebug(devtools::document)
# devtools::document()
# undebug(roxygen2:::roclet_process)
example1  <- fit.models.expert(formula=Surv(time2,status2)~1,data=data2,
                               distr=c("gga", "wph", "wei"),
                               method="hmc",
                               iter = 1000,
                               pool_type = "log pool", 
                               opinion_type = "survival",
                               times_expert = timepoint_expert, 
                               param_expert = param_expert_example1)

# plot(example1, add.km = T, t = 0:30)
# print(example1)
# model.fit.plot(example1, type = "waic")

undebug(compute_ICs_stan)
undebug(expert_like)
undebug(expert_log_dens)
compute_ICs_stan(example1$models$`Gen. Gamma`, distr3 = "gga", example1$misc$data.stan[[1]])
undebug(lik_wei)
undebug(lik_wph)
undebug(expert_like)
undebug(expert_log_dens)

compute_ICs_stan(example1$models$`Weibull (PH)`, distr3 = "wph", example1$misc$data.stan[[2]])


example1$models$`Weibull (PH)`

param_expert_example2 <- list()

param_expert_example2[[1]] <- data.frame(dist = c("beta"),
                                         wi = c(1), # Ensure Weights sum to 1
                                         param1 = c(1),
                                         param2 = c(1),
                                         param3 = c(NA))


example2  <- fit.models.expert(formula=Surv(time2,status2)~1,data=data2,
                               distr=c("gga", "wph", "wei"),
                               method="hmc",
                               iter = 1000,
                               pool_type = "log pool", 
                               opinion_type = "survival",
                               times_expert = timepoint_expert, 
                               param_expert = param_expert_example2)


expertsurv::print.expertsurv(example1)


plot.expertsurv(example1, add.km = T, t = 0:30)+
  theme_light()+
  scale_x_continuous(expand = c(0, 0), limits = c(0,NA), breaks=seq(0, 30, 2)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA), breaks=seq(0, 1, 0.05))+
  geom_segment(aes(x = 14, y = cred_int_val[1], xend = 14, yend = cred_int_val[2]))


plot.expertsurv(example2, add.km = T, t = 0:30)+
  theme_light()+
  scale_x_continuous(expand = c(0, 0), limits = c(0,NA), breaks=seq(0, 30, 2)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA), breaks=seq(0, 1, 0.05))#+
#  geom_segment(aes(x = 14, y = cred_int_val[1], xend = 14, yend = cred_int_val[2]))

model.fit.plot(example1, type = "dic")
model.fit.plot(example2, type = "dic")
example1$model.fitting$dic[2] -example12$model.fitting$dic[2] 
log(1.78)*2 #DIC difference is approximately 2 times the log density for the pooled expert at the posterior median of the parameters


#Sometimes the improvement in Likelihood due to being close to the expert's opinion is offset by the worse fit to the data

#Need to supress those warning messages

example1_mle  <- fit.models.expert(formula=Surv(time2,status2)~1,data=data2,
                                   distr=c("wph", "gomp", "lnorm", "gga"),
                                   method="mle",
                                   iter = 5000,
                                   pool_type = "log pool", 
                                   opinion_type = "survival",
                                   times_expert = timepoint_expert, 
                                   param_expert = param_expert_example1)


plot(example1_mle, add.km = T, t = 0:30)+
  theme_light()+
  scale_x_continuous(expand = c(0, 0), limits = c(0,NA), breaks=seq(0, 30, 2)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA), breaks=seq(0, 1, 0.05))+
  geom_segment(aes(x = 14, y = cred_int_val[1], xend = 14, yend = cred_int_val[2]))

example1_mle$models$`Gen. Gamma`
model.fit.plot(example1_mle, type = "aic")

#Confirm the AIC
library("flexsurv")


n_par <- 2
distr<- "gomp"

if(distr == "wph"){
  coef_mle  <- example1_mle$models$`Weibull (PH)`[["res"]][, "est"]
  LL_mod <- sum(pweibullPH(data2$time2,shape = coef_mle[1], scale =  coef_mle[2], lower.tail = F, log = T),
                hweibullPH(data2$time2,shape = coef_mle[1], scale =  coef_mle[2], log = T)*data2$status2) 
  St_mle <- pweibullPH(timepoint_expert,shape = coef_mle[1], scale =  coef_mle[2], lower.tail = F, log = F)
  
}else{
  coef_mle  <- example1_mle$models$Gompertz[["res"]][, "est"]
  LL_mod <- sum(pgompertz(data2$time2,shape = coef_mle[1], rate =  coef_mle[2], lower.tail = F, log = T),
                hgompertz(data2$time2,shape = coef_mle[1], rate =  coef_mle[2], log = T)*data2$status2) 
  St_mle <- pgompertz(timepoint_expert,shape = coef_mle[1], rate =  coef_mle[2], lower.tail = F, log = F)

}



expert1<- param_expert_example1[[1]][1,]
expert2<- param_expert_example1[[1]][2,]
Surv_prob <- seq(0, 1, by = 0.01)

dens1 <- dnorm(Surv_prob, expert1$param1,expert1$param2)
dens2 <- expertsurv:::dt.scaled(Surv_prob, mean = expert2$param1,sd = expert2$param2, df= expert2$param3)

log_pool <- dens1^expert1$wi*dens2^expert2$wi
k  <- sfsmisc:::integrate.xy(Surv_prob, log_pool)
log_pool2 <- log_pool/k

dens_eval1 <- dnorm(St_mle, expert1$param1,expert1$param2)
dens_eval2 <- expertsurv:::dt.scaled(St_mle, mean = expert2$param1,sd = expert2$param2, df= expert2$param3)

LL_expert <- log((dens_eval1^expert1$wi)*(dens_eval2^expert2$wi)/k)
(LL_mod+LL_expert)*(-2) +n_par*2

AIC(example1_mle$models$`Weibull (PH)`)
AIC(example1_mle$models$Gompertz)


pweibullPH(timepoint_expert,shape = coef_mle[1], scale =  exp(-4.551401 ), lower.tail = F, log = F)





plot(flexsurvreg(formula=Surv(time2,status2)~1,data=data2,
            dist="gengamma.orig"))


param_expert3 <- list()

#Prior belief of 5 "months" difference in expected survival
param_expert3[[1]] <- data.frame(dist = "norm", wi = 1, param1 = 5, param2 = 0.2, param3 = NA)


survHE.data.model  <- fit.models.expert(formula=Surv(time2,status2)~as.factor(arm),data=data2,
                                        distr=c("wei"),
                                        method="hmc",
                                        iter = 5000,
                                        opinion_type = "mean",
                                        id_trt = 1, # Survival difference is Mean_surv[id_trt]- Mean_surv[id_comp] 
                                        param_expert = param_expert3)
model.fit.plot(survHE.data.model, type = "dic")


survHE.data.model2  <- fit.models.expert(formula=Surv(time2,status2)~as.factor(arm),data=data2,
                                        distr=c("wei"),
                                        method="mle",
                                        iter = 5000,
                                        opinion_type = "mean",
                                        id_trt = 1, # Survival difference is Mean_surv[id_trt]- Mean_surv[id_comp] 
                                        param_expert = param_expert3)


model.fit.plot(survHE.data.model2, type = "aic")

plot(survHE.data.model,survHE.data.model2, add.km = T, t = 0:30)

