
# Expert Opinion R script ----

library("survHE")
library("SHELF")
library("xlsx")
library("abind")
library("crayon")
library("flexsurv")
library("expertsurv")
library("kableExtra")
devtools::load_all("C:/Users/phili/OneDrive/PhD/R packages/expertsurv/")
# remove.packages("expertsurv")


## Example Publication: ELIANA Trial Data ---- 

pathway <- "C:/Users/phili/OneDrive/PhD/R packages/expertsurv/"

Survival.df <- read.xlsx(paste0(pathway, "ELIANA OS.xlsx"), 1) %>% data.frame()
times.risk <- seq(0, to = 34, by = 2)
n.risk <- c(79,76,73,68,67,62,55,52,47,42,39,26,21,14,9,5,2,0)
upper <- sapply(times.risk, function(x){sum(x > Survival.df$time)})

times.risk<- times.risk[!duplicated(upper)]
n.risk <- n.risk[!duplicated(upper)]
upper <- upper[!duplicated(upper)]

upper <- upper[-1]
lower <- upper
lower <- lower[-length(lower)]
lower <- lower +1
lower <- c(1, lower)


write.table(data.frame(Time = Survival.df$time,
                       Survival =Survival.df$surv),
            paste0(pathway,"survival.txt"),
            row.names=T,
            sep="\t")

write.table(data.frame(Interval = 1:length(lower),
                       time = times.risk[-which(times.risk>tail(Survival.df$time, n= 1))],
                       lower =lower, upper = upper,
                       nrisk =  n.risk[-which(times.risk>tail(Survival.df$time, n= 1))]),
            paste0(pathway,"nrisk.txt"),
            row.names=FALSE,
            sep="\t")

#digitize function exports it so it needs to be read back in 

digitise(paste0(pathway,"survival.txt"),
         paste0(pathway,"nrisk.txt"),
         km_output = paste0(pathway, "KMdata.txt"),
         ipd_output = paste0(pathway, "IPDdata.txt"))

digitized_IPD <- read.table(paste0(pathway,"IPDdata.txt"))

colnames(digitized_IPD) <- digitized_IPD[1,]

digitized_IPD <- data.frame(digitized_IPD[-1,])

digitized_IPD <- apply(digitized_IPD,2, as.numeric)
km.fit <- survfit(Surv(time, event)~1, data.frame(digitized_IPD))

#Plot the Kaplan Meier curve
#plot(km.fit)

# Read in Expert opinions
expert.prob.df <- read.xlsx(paste0(pathway, "Expert Surv Probabilities.xlsx"), sheetName = "Sheet1") %>% data.frame()

j <- q <- 1

res.mat <- matrix(nrow = nrow(expert.prob.df)/3, ncol = 3)
for(i in 1:(nrow(expert.prob.df)/3)){
        
        res.mat[j,]    <- expert.prob.df[q:(q+2), 2]   
        q <- q+3
        j <- j +1
        
}

expert.prob.df2 <- cbind(res.mat, rep(1:7, times = 4), rep(2:5, each = 7))
colnames(expert.prob.df2) <- c("0.99", "Mode","0.01", "Expert","Time")

times_act <- 4:5

dfs_expert <- list()
plts_pool <- list()
dfs_pool <- list()
fit.eval <- list()
lower_bound <- 0
upper_bound <- 1

for(i in 1:length(times_act)){

expert.prob.eval <- expert.prob.df2[expert.prob.df2[,5] == times_act[i],1:3]
expert.prob.eval <- expert.prob.eval[,c(3,2,1)]


fit.eval[[i]] <- expertsurv:::fitdist_mod(t(expert.prob.eval)[-2,],
        probs = c(0.005,  0.995), upper = upper_bound, lower = lower_bound ,mode = t(expert.prob.eval)[2,], 
        expertnames = c(1:7))

}


#sum ssq error over all timepoints; so that the distribution is the same for each expert
ssq_total <- fit.eval[[1]]$ssq

for(i in 2:length(times_act)){
  ssq_total <- ssq_total +fit.eval[[i]]$ssq

}

dist_considered <- c("normal","t","gamma", "lognormal", "beta") 


ssq_total<- ssq_total[,dist_considered]

best_fit_overall <-colnames(ssq_total)[apply(ssq_total, 1, function(x){which.min(x)})]


for(i in 1:length(times_act)){
  fit.eval[[i]]$best.fitting[,"best.fit"] <- best_fit_overall
  
}



for(i in 1:length(times_act)){
# Appears to be a Slight bug with this function (SHELF version)
#SHELF:::plotfit(fit.eval, lp = T)

plts_pool[[i]] <- expertsurv:::makePoolPlot(fit= fit.eval[[i]],
                                       xl =lower_bound,
                                       xu =upper_bound,
                                       d = "best",
                                       w = 1,
                                       lwd =1,
                                       xlab = "x",
                                       ylab =expression(f[X](x)),
                                       legend_full = TRUE,
                                       ql = NULL,
                                       qu = NULL,
                                       nx = 200,
                                       addquantile = FALSE,
                                       fs = 12,
                                       expertnames = paste0("Expert ",1:nrow(fit.eval[[i]]$Normal)))+
  theme_bw()

plts_pool[[i]] <- plts_pool[[i]]+
  scale_color_brewer(palette = "Paired")+
  ylim(0,30)+
  xlab("S(t)")+
  theme(legend.title=element_blank())
dfs_pool[[i]] <-  plts_pool[[i]][["data"]]

#best_fit_index  <- apply(fit.eval[[i]]$ssq[,dist_considered], 1, which.min)

best_fit <- fit.eval[[i]]$best.fitting[,"best.fit"]

best_fit_loc <- sapply(best_fit, function(x){which(x  == names(fit.eval[[i]]$ssq))})
fit.eval.dist  <- fit.eval[[i]][best_fit_loc]

pool.df_output <- matrix(nrow = length(best_fit_loc),ncol = 3)
colnames(pool.df_output) <- c("param1", "param2", "param3")

for(j in 1:length(best_fit_loc)){
        pool.df_output[j,
                       1:length(fit.eval.dist[[j]][j,])] <-  as.numeric(as.vector(fit.eval.dist[[j]][j,]))
}
dfs_expert[[i]] <- data.frame(dist = best_fit, wi = 1/nrow(pool.df_output), pool.df_output)
}


ggpubr::ggarrange(plts_pool[[1]]+
            theme(legend.title = element_blank()), plts_pool[[2]]+ ylab("")+
            theme(legend.title = element_blank()),  
          labels = c("Year 4", "Year 5"),
          vjust =1,
          ncol = 2, nrow = 1,common.legend = T,legend="bottom")

ggsave(filename = paste0(pathway, "Year 4 & 5 Distributions.png"),
       width =  9, height = 4.5)


dist_fit <- c()
pool.eval.stan <- list()
pool_type_eval <- "linear pool"

#Expert Opinion
param_expert <- dfs_expert
save(param_expert, file = paste0(pathway, "param_expert.RData"))

digitized_IPD <- data.frame(digitized_IPD)

n.chains <- 3
k1 <- 1
# knots1 <- quantile(log((digitized_IPD %>% filter(event == 1))$time),
#                    seq(0, 1, length = k1 + 2))


mle.ests_rps <- flexsurvspline(Surv(time, event) ~ 1, data=digitized_IPD, k=1, scale="hazard")
mle.ests_lno <- flexsurvreg(Surv(time, event) ~ 1, data=digitized_IPD, dist="lno")
mle.ests_exp <- flexsurvreg(Surv(time, event) ~ 1, data=digitized_IPD, dist="exp")
mle.ests_wph <- flexsurvreg(Surv(time, event) ~ 1, data=digitized_IPD, dist="weibullPH")
mle.ests_llo <- flexsurvreg(Surv(time, event) ~ 1, data=digitized_IPD, dist="llogis")
mle.ests_wei <- flexsurvreg(Surv(time, event) ~ 1, data=digitized_IPD, dist="weibull")
mle.ests_wph <- flexsurvreg(Surv(time, event) ~ 1, data=digitized_IPD, dist="weibullPH")
mle.ests_gom <- flexsurvreg(Surv(time, event) ~ 1, data=digitized_IPD, dist="gompertz")
mle.ests_gam <- flexsurvreg(Surv(time, event) ~ 1, data=digitized_IPD, dist="gamma")
mle.ests_gga<- flexsurvreg(Surv(time, event) ~ 1, data=digitized_IPD, dist="gengamma.orig")


init_fun_gga <- function(...){
  
  res <- as.numeric(mvtnorm::rmvnorm(n = 1, mean = mle.ests_gga$res.t[,1], sigma = mle.ests_gga$cov/20))
  beta_jags = c(-res[2],0)
  
  list(beta_exp = pmin(exp(beta_jags),1000),
       r = min(exp(res[1]),100),
       b = min(exp(res[3]), 100))
}

init_fun_gomp <- function(...){
 
  res <- as.numeric(mvtnorm::rmvnorm(n = 1, mean = mle.ests_gom$res.t[,1], sigma = mle.ests_gom$cov))
  beta <- c(res[2],0)
  
  list(beta_exp = exp(beta),
       alpha = res[1])
}


init_fun_gam <- function(...){
  res <- as.numeric(mvtnorm::rmvnorm(n = 1, mean = mle.ests_gam$res.t[,1], sigma = mle.ests_gam$cov))
  beta <- c(res[2],0)
  
  list(beta_exp = exp(beta),
       alpha = exp(res[1]))
}



#We have to initialize the MCMC sampler at parameters which have positive probability
#we use a Multinormal Density centred at the MLE with the covariance matrix (both obtained from flexsurv) 


init_fun_gga <- function(...){
  
  #res <- as.numeric(mvtnorm::rmvnorm(n = 1, mean = mle.ests_gga$res.t[,1], sigma = mle.ests_gga$cov/20))
  res <- mle.ests_gga$res.t[,1]
  beta_jags = c(-res[2],0)
  
  list(beta_exp = pmin(exp(beta_jags),100),
       r = min(exp(res[1]),100),
       b = min(exp(res[3]), 100))
}

init_fun_gomp <- function(...){
  
  res <- as.numeric(mvtnorm::rmvnorm(n = 1, mean = mle.ests_gom$res.t[,1], sigma = mle.ests_gom$cov))
  beta <- c(res[2],0)
  
  list(beta_exp = exp(beta),
       alpha = res[1])
}

init_fun_gam <- function(...){
  res <- as.numeric(mvtnorm::rmvnorm(n = 1, mean = mle.ests_gam$res.t[,1], sigma = mle.ests_gam$cov))
  beta <- c(res[2],0)
  
  list(beta_exp = exp(beta),
       alpha = exp(res[1]))
}

init_fun_rps <- function(...){list(gamma=as.numeric(mvtnorm::rmvnorm(n = 1, mean = mle.ests_rps$res[,1],
                                                                     sigma = mle.ests_rps$cov)))}
init_fun_lnorm <- function(...){
  #Can generate values outside allowable range
  res <- as.numeric(mvtnorm::rmvnorm(n = 1, mean = mle.ests_lno$res.t[,1], sigma = mle.ests_lno$cov/20))
  
  beta <- c(res[1],0)
  
  list(beta = beta,
       alpha = exp(res[2]))
}

init_fun_llogis <- function(...){
  res <- as.numeric(mvtnorm::rmvnorm(n = 1, mean = mle.ests_llo$res.t[,1], sigma = mle.ests_llo$cov))
  beta <- c(res[2],0)
  
  list(beta_exp = exp(beta),
       alpha = exp(res[1]))
}

init_fun_exp <- function(...){
  # res <- as.numeric(rnorm(n = 1, mean = mle.ests_exp$res.t[,1]),sd =  mle.ests_exp$res.t[,4]/50)
  res <-  as.numeric(mle.ests_exp$res.t[,1])
  list(beta_exp =exp(c(res[1],0)))
}

init_fun_wph <- function(...){
  res <- as.numeric(mvtnorm::rmvnorm(n = 1, mean = mle.ests_wph$res.t[,1], sigma = mle.ests_wph$cov))
  beta <- c(res[2],0)
  
  list(beta_exp = exp(beta),
       alpha = exp(res[1]))
}

init_fun_wei <- function(...){
  res <- as.numeric(mvtnorm::rmvnorm(n = 1, mean = mle.ests_wei$res.t[,1], sigma = mle.ests_wei$cov/10))
  beta <- c(res[2],0)
  list(beta_exp = exp(beta),
       alpha = exp(res[1]))
}


m.all_expert <- expertsurv:::fit.models.expert(formula=Surv(time,event)~1,data=digitized_IPD,
                                               distr=c("exp","llo","lno", "wei","rps", "gam", "gomp", "gga"),method="hmc",
                                               opinion_type = "survival",
                                               param_expert = param_expert,
                                               pool_type = "linear pool",
                                               k = k1,
                                               times_expert = c(4,5)*12,
                                               iter = 10000,
                                               #a0 = rep(0.00000001, nrow(digitized_IPD_no_data)),
                                               chains =n.chains,
                                               thin = 1,
                                               init = list(
                                                 rps = lapply(rep(1, n.chains), init_fun_rps),
                                                 lno = lapply(rep(1, n.chains), init_fun_lnorm),
                                                 exp = lapply(rep(1, n.chains), init_fun_exp),
                                                 wph = lapply(rep(1, n.chains), init_fun_wph),
                                                 llo = lapply(rep(1, n.chains), init_fun_llogis),
                                                 exp = lapply(rep(1, n.chains), init_fun_exp),
                                                 wei = lapply(rep(1, n.chains), init_fun_wei),
                                                 gom = lapply(rep(1, n.chains), init_fun_gomp),
                                                 gam = lapply(rep(1, n.chains), init_fun_gam),
                                                 gga = lapply(rep(1, n.chains), init_fun_gga)),
                                               St_lower = 0,
                                               St_upper = 1)




models <- names(m.all_expert$models)
psa_outuput <- list()

for(i in 1:length(m.all_expert$models)){
  psa <- make.surv(fit = m.all_expert,mod = i, nsim = 1000, t = seq(0, 80, by = 1))
  #psa.plot(psa)
  df_temp  <- t(apply(psa$mat[[1]], 1,quantile, probs = c(0.025, 0.5,.975))) %>% data.frame()
  df_temp$time <- seq(0, 80, by = 1)
  mod_name <- names(m.all_expert$models)[i]
  psa_outuput[[mod_name]] <- df_temp %>% mutate(model = mod_name)
  
}

df_final <- do.call(rbind.data.frame, psa_outuput)

km.fit   <- survfit(Surv(time,event)~1,data=digitized_IPD)


datakm <- data.frame(time =km.fit$time, survival=  km.fit$surv, lower = km.fit$lower, upper = km.fit$upper)
times_expert = c(4,5)*12
scale <- 1
pool_type_eval <- "linear pool"
df.linear <- subset(dfs_pool[[1]], ftype == pool_type_eval) %>% rename(y = x) %>% 
  mutate(x = times_expert[1] + fx*scale)
df.linear2 <- subset(dfs_pool[[2]], ftype == pool_type_eval) %>% rename(y = x) %>% 
  mutate(x = times_expert[2] + fx*scale)


df_final$Models <- factor(df_final$model, levels = unique(df_final$model))

plt.expert <- ggplot(data = df_final, aes(y = X50., x = time, group = Models, colour = Models))+
  geom_line()+
   geom_line(aes(y = X97.5., x = time),linetype="dashed")+
   geom_line(aes(y = X2.5., x = time), linetype="dashed")+
  # geom_errorbar(data = df_final%>% filter(time == 80),aes(ymin=X2.5., ymax=X97.5.), width=8,
  #               position=position_dodge(2))+
  theme_light()+
  scale_x_continuous(expand = c(0, 0), limits = c(0,NA), breaks=seq(0, 80, 5)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA), breaks=seq(0, 1, 0.05))+
  xlab("Time (months)")+
  ylab("Survival")+
  geom_path(data = df.linear, aes(x =x, y =y), colour = "grey", lwd=1.1, inherit.aes = F)+
  geom_path(data = data.frame(grp = rep(times_expert[1], 10), 
                              y = seq(0,1, length.out =10)), aes(x = grp, y= y), colour = "grey", lwd=1.1, inherit.aes = F)+
  geom_polygon(data = df.linear, aes(y = y, x = x), fill = "sky blue", alpha = 0.5,inherit.aes = F)+
  geom_path(data = df.linear2, aes(x =x, y =y), colour = "grey", lwd=1.1, inherit.aes = F)+
  geom_polygon(data = df.linear2, aes(y = y, x = x), fill = "sky blue", alpha = 0.5,inherit.aes = F)+
  geom_path(data = data.frame(grp = rep(times_expert[2], 10), 
                              y = seq(0,1, length.out =10)), aes(x = grp, y= y), colour = "grey", lwd=1.1, inherit.aes = F) +
    geom_step(data = datakm, aes(x = time , y = survival), color = "black",inherit.aes = F) + 
    geom_ribbon(data = datakm, aes(x = time, y = survival, ymin = lower, 
                                   ymax = upper), alpha = 0.2, color = "grey",inherit.aes = F)
plt.expert<-  plt.expert+ scale_fill_discrete(name ="Models")
  


ggsave("plt.png")




# gg.obj2_stan <- ggs(m.all_expert$models$`log-Normal`) %>% filter(Parameter %in%  c("beta[1]","St_expert[1]"))
# ggmcmc(gg.obj2_stan,file="model_simple-diag-stan.pdf")



AUC_res_main <-  matrix(NA, nrow = length(models), ncol = 1)
for(j in 1:length(models)){
  AUC_res_main[j,1] <- mean(summary(m.all_expert, mod = j, t = seq(0, 80, by = 1))[["mean.surv"]])
}


#Plot the distributions 

times_expert = c(4,5)*12
scale <- 1
pool_type_eval <- "linear pool"
df.linear <- subset(dfs_pool[[1]], ftype == pool_type_eval) %>% rename(y = x) %>% 
  mutate(x = times_expert[1] + fx*scale)
df.linear2 <- subset(dfs_pool[[2]], ftype == pool_type_eval) %>% rename(y = x) %>% 
  mutate(x = times_expert[2] + fx*scale)

library(RColorBrewer)
display.brewer.all(colorblindFriendly = TRUE)

plot.survHE(m.all_expert, add.km = T, t = 0:80)+
  scale_color_brewer(palette = "Set1")


plt.expert  <- plot.survHE(m.all_expert, add.km = T, t = 0:80)+
  theme_light()+
  scale_x_continuous(expand = c(0, 0), limits = c(0,NA), breaks=seq(0, 80, 5)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA), breaks=seq(0, 1, 0.05))+
  xlab("Time (months)")+
  geom_path(data = df.linear, aes(x =x, y =y), colour = "grey", lwd=1.1, inherit.aes = F)+
  geom_path(data = data.frame(grp = rep(times_expert[1], 10), 
                              y = seq(0,1, length.out =10)), aes(x = grp, y= y), colour = "grey", lwd=1.1, inherit.aes = F)+
  geom_polygon(data = df.linear, aes(y = y, x = x), fill = "sky blue", alpha = 0.5)+
  geom_path(data = df.linear2, aes(x =x, y =y), colour = "grey", lwd=1.1, inherit.aes = F)+
  geom_polygon(data = df.linear2, aes(y = y, x = x), fill = "sky blue", alpha = 0.5)+
  geom_path(data = data.frame(grp = rep(times_expert[2], 10), 
                              y = seq(0,1, length.out =10)), aes(x = grp, y= y), colour = "grey", 
                              lwd=1.1, inherit.aes = F)






if (add.km == TRUE) {
  datakm = lapply(1:length(survHE_objs), function(i) {
    make_data_surv(survHE_objs[[i]], mods = 1, nsim = 1, 
                   t = t, newdata = newdata, add.km = add.km)[[2]] %>% 
      mutate(object_name = as.factor(names(survHE_objs)[i]))
  }) %>% bind_rows() %>% group_by(object_name, model_name) %>% 
    mutate(mods_id = cur_group_id()) %>% ungroup()
}


# psa <- make.surv(fit = m.all, nsim = 1000, t = seq(.1, 80), mod = 1)
# psa.plot(psa, xlab = "Extrapolated time",
#          ylab = "Estimation of the survival curves")

#model.fit.plot(m.all_expert,type = "dic")

## issues ??
# Gompertz  and Generalized Gamma 
# survHE.mods  <- survHE:::fit.models(formula=Surv(time,event)~1,data=digitized_IPD,
#                                                distr=c("exp","llo","lno", "wei","rps", "gam", "gomp", "gga"),
#                                                method="hmc",
#                                                iter = 5000,
#                                     priors=list(gom=list(a_alpha=0.1,b_alpha=0.1)))


param_expert_vague <- list()

param_expert_vague[[1]] <- data.frame(dist = "beta", wi = 1, param1 = 1, param2 = 1, param2 = NA)

m.all_expert_vague <- expertsurv:::fit.models.expert(formula=Surv(time,event)~1,data=digitized_IPD,
                                               distr=c("exp","llo","lno", "wei","rps", "gam", "gomp", "gga"),method="hmc",
                                               opinion_type = "survival",
                                               param_expert = param_expert_vague,
                                               k = k1,
                                               times_expert = c(4)*12,
                                               iter = 10000,
                                               chains =n.chains,
                                               priors = list(lno = list(sigma_beta = c(5,5)),
                                                             gga = list(sigma_beta = c(1,1))),
                                               init = list(rps = lapply(rep(1, n.chains), init_fun_rps),
                                                           lno = lapply(rep(1, n.chains), init_fun_lnorm),
                                                           exp = lapply(rep(1, n.chains), init_fun_exp)))

plt.vague <- plot(m.all_expert_vague, add.km = T, t = 0:80)+
  theme_light()+
  scale_x_continuous(expand = c(0, 0), limits = c(0,NA), breaks=seq(0, 80, 5)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA), breaks=seq(0, 1, 0.05))+
  xlab("Time (months)")+
  ylab("")



models <- names(m.all_expert_vague$models)
psa_outuput <- list()

for(i in 1:length(m.all_expert_vague$models)){
  psa <- make.surv(fit = m.all_expert_vague,mod = i, nsim = 1000, t = seq(0, 80, by = 1))
  #psa.plot(psa)
  df_temp  <- t(apply(psa$mat[[1]], 1,quantile, probs = c(0.025, 0.5,.975))) %>% data.frame()
  df_temp$time <- seq(0, 80, by = 1)
  mod_name <- names(m.all_expert_vague$models)[i]
  psa_outuput[[mod_name]] <- df_temp %>% mutate(model = mod_name)
  
}

df_final <- do.call(rbind.data.frame, psa_outuput)

plt.vague <- ggplot(data = df_final, aes(y = X50., x = time, group = factor(model), colour = factor(model)))+
  geom_line()+
  # geom_line(aes(y = X97.5., x = time),linetype="dotdash")+
  # geom_line(aes(y = X2.5., x = time), linetype="dotdash")+
  geom_errorbar(data = df_final%>% filter(time == 80),aes(ymin=X2.5., ymax=X97.5.), width=8,
                position=position_dodge(2))+
  theme_light()+
  scale_x_continuous(expand = c(0, 0), limits = c(0,NA), breaks=seq(0, 80, 5)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA), breaks=seq(0, 1, 0.05))+
  xlab("Time (months)")+
  ylab("Survival")+
  geom_step(data = datakm, aes(x = time , y = survival), color = "black",inherit.aes = F) + 
  geom_ribbon(data = datakm, aes(x = time, y = survival, ymin = lower, 
                                 ymax = upper), alpha = 0.2, color = "grey",inherit.aes = F)




### Outputs for Publication ----

DIC_comp <- data.frame(Model =names(m.all_expert$models),
                       DIC_expert = m.all_expert$model.fitting$dic ) %>% arrange(DIC_expert) %>% mutate_if(is.numeric, round,digits = 2)


DIC_comp <- data.frame(Model =names(m.all_expert_vague$models),
           DIC_expert = m.all_expert$model.fitting$dic,
           DIC_vague = m.all_expert_vague$model.fitting$dic) %>% mutate_if(is.numeric, round,digits = 2)




kable(DIC_comp[order(DIC_comp$DIC_expert),], format="latex", row.names = F)
model.fit.plot(m.all_expert_vague,type = "dic")
#95% credible interval for the Year 4 (linear pool)

vals <- seq(0,1, by = 0.01)
quant.vec <- rep(NA, length(vals))
for(i in 1:length(vals)){
  quant.vec[i] <- mean(t(apply(param_expert[[1]], 1, function(x){get_cdf_val(
    dist = x["dist"],
    param1 = x["param1"],
    param2 = x["param2"],
    param3 = x["param3"],
    vals = vals[i])})))
  
}

vals[which.min(abs(quant.vec-0.025))]
vals[which.min(abs(quant.vec-0.975))]


### Leave One Out Analysis ----#


AUC_res <- DIC_res <- matrix(NA, nrow = length(models), ncol = nrow(param_expert[[1]]))


for(i in 1:nrow(param_expert[[1]])){
  param_expert_loo <- list()
  for(k in 1:2){
    param_expert_loo[[k]] <- param_expert[[k]][-i,]
  }
  m.all_expert_loo <- expertsurv:::fit.models.expert(formula=Surv(time,event)~1,data=digitized_IPD,
                                                 distr=c("exp","llo","lno", "wei","rps", "gam", "gomp", "gga"),method="hmc",
                                                 opinion_type = "survival",
                                                 param_expert = param_expert_loo,
                                                 pool_type = "linear pool",
                                                 k = k1,
                                                 times_expert = c(4,5)*12,
                                                 iter = 1000,
                                                 chains =n.chains,
                                                 priors = list(lno = list(sigma_beta = c(5,5)),
                                                               gga = list(sigma_beta = c(1,1))),
                                                 init = list(rps = lapply(rep(1, n.chains), init_fun_rps),
                                                             lno = lapply(rep(1, n.chains), init_fun_lnorm)))
  
  DIC_res[,i] <- m.all_expert_loo$model.fitting$dic
   for(j in 1:length(models)){
    AUC_res[j,i] <- mean(summary(m.all_expert_loo, mod = j, t = seq(0, 80, by = 1))[["mean.surv"]])
  }
 
}


expert_num <- 7
kendall_cor <- rep(NA, expert_num )
library("Cairo")
for(i in 1:expert_num){
  
#https://rstudio-pubs-static.s3.amazonaws.com/297778_5fce298898d64c81a4127cf811a9d486.html
  par(mar=c(1.2,1.2,1.2,1.2))
  Cairo::Cairo( file=paste0(pathway,"DIC without Expert ",i,".png"), type="png", bg="white", dpi= 130, width = 600, height = 600)
plot(x =rank(DIC_comp[,2]), y = apply(DIC_res, 2, rank)[,i],
     xlab = "", ylab = "",
     xaxt = "n",
     xaxt = "n", yaxt = "n")
# axis(1, at = rank(DIC_comp[,2]),
#      labels = models,las = 2, cex.axis=0.7)
text(x=3, y=5, "Rank from 1 (at origin) \n upwards",col="gold4", font=2, cex=0.60)
title(xlab = "DIC ranking with all Experts", line = 0)            # Add x-axis text
title(ylab = paste0("DIC ranking without Expert ", i ), line = 0)  
abline(h=seq(1,nrow(DIC_comp),1), v=seq(1,nrow(DIC_comp), 1), lty=3, col="gray")
text(x=rank(DIC_comp[,2]),  par("usr")[3]-0.7, 
     labels = models, srt = 45, pos = 1, xpd = TRUE, cex = 0.60)
text(par("usr")[1]-0.3, apply(DIC_res, 2, rank)[,i],  
     labels = models, srt = 45, pos = 2, xpd = TRUE,cex = 0.60)
dev.off()

kendall_cor[i] <- cor(rank(DIC_comp[,2]),apply(DIC_res, 2, rank)[,i], method = "kendall" )
}
# Add y-axis text

save(kendall_cor, file = paste0(pathway,"kendall.RData"))

AUC_res2 <- AUC_res
for(i in 1:ncol(AUC_res2)){
  
  AUC_res2[,i] <- round((AUC_res[,i] - AUC_res_main[,1])*100/AUC_res_main[,1], digits = 2)
}

save(AUC_res2, file = paste0(pathway,"AUC_res2.RData"))



# 
# undebug(compute_ICs_stan)
# compute_ICs_stan(m.all_expert_vague$models$Gompertz, "gom",m.all_expert_vague$misc$data.stan[[7]])
figure <- ggpubr::ggarrange(plt.expert, plt.vague, 
                    ncol = 2, nrow = 1,common.legend = TRUE, legend="bottom")
ggsave(paste0(pathway,"Survival functions with Expert Opinion and default priors.png"), width =  10, height = 6.66)


## Example Vignette 1: Opinion on the survival of a treatment ----

data2 <- survHE::data %>% rename(status = censored) %>% mutate(time2 = ifelse(time > 10, 10, time),
                                                              status2 = ifelse(time> 10, 0, status))


param_expert_example1 <- list()

param_expert_example1[[1]] <- data.frame(dist = c("norm","t"),
                                         wi = c(1,1),
                                         param1 = c(0.1,0.12),
                                         param2 = c(0.005,0.005),
                                         param3 = c(NA,3))
timepoint_expert <- c(14)

plot_opinion1<- plot_expert_opinion(param_expert_example1[[1]], 
                    weights = param_expert_example1[[1]]$wi)
ggsave("Vignette_Example 1 - Expert Opinion.png")

cred_int_val <- cred_int(plot_opinion1,val = "linear pool", interval = c(0.025, 0.975))

example1  <- expertsurv:::fit.models.expert(formula=Surv(time2,status2)~1,data=data2,
                                            distr=c("wei", "gomp"),
                                            method="hmc",
                                            iter = 5000,
                                            opinion_type = "survival",
                                            times_expert = timepoint_expert, 
                                            param_expert = param_expert_example1)

model.fit.plot(example1, type = "dic")
ggsave("Vignette_Example 1 - DIC.png")


plot(example1, add.km = T, t = 0:30)+
  theme_light()+
  scale_x_continuous(expand = c(0, 0), limits = c(0,NA), breaks=seq(0, 30, 2)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA), breaks=seq(0, 1, 0.05))+
  geom_segment(aes(x = 14, y = cred_int_val[1], xend = 14, yend = cred_int_val[2]))
ggsave("Vignette_Example 1.png")



## Example Vignette 2: Opinion on the survival of the comparator arm ----
param_expert_example2 <- list()
param_expert_example2[[1]] <- data.frame(dist = c("norm"),
                                         wi = c(1),
                                         param1 = c(0.1),
                                         param2 = c(0.005),
                                         param3 = c(NA))


example2  <- expertsurv:::fit.models.expert(formula=Surv(time2,status2)~as.factor(arm),data=data2,
                                        distr=c("wei"),
                                        method="hmc",
                                        iter = 5000,
                                        opinion_type = "survival",
                                        id_St = 0, 
                                        times_expert = 14, 
                                        param_expert = param_expert_example2)


#debug(compute_ICs_stan)
#undebug(fit.models.expert)
plot.exam_2.expert <- plot(example2, add.km = T, t = 0:30)+
  theme_light()+
  scale_x_continuous(expand = c(0, 0), limits = c(0,NA), breaks=seq(0, 80, 2)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA), breaks=seq(0, 1, 0.05))


param_expert_vague <- list()

param_expert_vague[[1]] <- data.frame(dist = "beta", wi = 1, param1 = 1, param2 = 1, param2 = NA)


example2.vague  <- expertsurv:::fit.models.expert(formula=Surv(time2,status2)~as.factor(arm),data=data2,
                                                     distr=c("wei"),
                                                     method="hmc",
                                                     iter = 5000,
                                                     opinion_type = "survival",
                                                     id_St = 0, 
                                                     times_expert = 20, 
                                                     param_expert = param_expert_vague)
plot.exam_2.vague<- plot(example2.vague, add.km = T, t = 0:30)+
  theme_light()+
  scale_x_continuous(expand = c(0, 0), limits = c(0,NA), breaks=seq(0, 80, 2)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA), breaks=seq(0, 1, 0.05))


figure <- ggpubr::ggarrange(plot.exam_2.expert, plot.exam_2.vague, 
                            ncol = 2, nrow = 1,common.legend = TRUE, legend="bottom")
(figure )
ggsave(paste0("Vignette_Example 2.png"), width =  10, height = 6.66)



## Example 3: Opinion about survival difference ----

param_expert3 <- list()

param_expert3[[1]] <- data.frame(dist = "norm", wi = 1, param1 = 5, param2 = 0.2, param3 = NA)


example3  <- expertsurv:::fit.models.expert(formula=Surv(time2,status2)~as.factor(arm),data=data2,
                                                     distr=c("gom"),
                                                     method="hmc",
                                                     iter = 5000,
                                                     opinion_type = "mean",
                                                     id_trt = 1, 
                                                     param_expert = param_expert3)

plot.exam_3.expert <- plot(example3, add.km = T, t = 0:50)

ggsave(paste0("Vignette_Example 3.png"), width =  10, height = 6.66)



mean_arms <- summary(survHE.data.model, t = c(0:50))

mean_diff <- apply(mean_arms[["mean.surv"]], 1, function(x){x[2]-x[1]})
summary(mean_diff)

plot(density(mean_diff))

## Note on Priors for Gompertz ----

#if we have expert opinion on the mean survival we must 


## Check Approximation vs Manual Integration (rmst_generic)----

cbind(head(survHE.data.model$models$Gompertz[["BUGSoutput"]]$sims.matrix),
      "mean_diff" = apply(head(survHE.data.model$models$Gompertz[["BUGSoutput"]]$sims.matrix), 1, function(x){flexsurv::mean_gompertz(x[2], exp(x[3]+x[4]))-
      flexsurv::mean_gompertz(x[2], exp(x[3]))}))

#To be fixed
#Survival plots - base

#Check Log-Like!
#Warning message for Pd



citation_expert <- citHeader("To cite expertsurv in publications use:")
citEntry(entry = "Article",
         title        = "Utilizing Expert Opinion to inform Extrapolation of Survival Models",
         author       = personList(as.person("Philip Cooney"),
                                   as.person("Arthur White")),
         journal      = "arXiv",
         year         = "2021",
         pages        = "1--13",
         url          = "https://arxiv.org/pdf/2112.02288.pdf",
         
         textVersion  =
           paste("Philp Cooney, Arthur White (2021).",
                 "Utilizing Expert Opinion to inform Extrapolation of Survival Models.",
                 "arXiv, 1-13.",
                 "URL https://arxiv.org/pdf/2112.02288.pdf.")
)

citation("expertsurv")


basecit <- system.file("CITATION", package="base")
source(basecit, echo=TRUE)
readCitationFile(basecit)


