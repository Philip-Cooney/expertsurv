
# Expert Opinion R script ----

library("survHE")
library("SHELF")
library("xlsx")
library("abind")
library("crayon")
library("flexsurv")
library("expertsurv")
library("kableExtra")
devtools::load_all()
# remove.packages("expertsurv")


## Example 1: ELIANA Trial Data ---- 

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
        probs = c(0.01,  0.99), upper = upper_bound, lower = lower_bound ,mode = t(expert.prob.eval)[2,], 
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
                                       expertnames = NULL)

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
dfs_expert[[i]] <- data.frame(dist = best_fit, wi = 1, pool.df_output)
}

dist_fit <- c()
pool.eval.stan <- list()
pool_type_eval <- "linear pool"

#Expert Opinion
param_expert <- dfs_expert


digitized_IPD <- data.frame(digitized_IPD)

n.chains <- 2
k1 <- 1
# knots1 <- quantile(log((digitized_IPD %>% filter(event == 1))$time),
#                    seq(0, 1, length = k1 + 2))


mle.ests_rps <- flexsurvspline(Surv(time, event) ~ 1, data=digitized_IPD, k=1, scale="hazard")
mle.ests_lno <- flexsurvreg(Surv(time, event) ~ 1, data=digitized_IPD, dist="lno")

#We have to initialize the MCMC sampler at parameters which have positive probability
#we use a Multinormal Density centred at the MLE with the covariance matrix (both obtained from flexsurv) 

init_fun_rps <- function(...){list(gamma=as.numeric(mvtnorm::rmvnorm(n = 1, mean = mle.ests_rps$res[,1],
                                                                     sigma = mle.ests_rps$cov)))}

init_fun_lnorm <- function(...){
  res <- as.numeric(mvtnorm::rmvnorm(n = 1, mean = mle.ests_lno$res.t[,1], sigma = mle.ests_lno$cov))
  list(beta =c(res[1],0),
       alpha = res[2])}


m.all_expert <- expertsurv:::fit.models.expert(formula=Surv(time,event)~1,data=digitized_IPD,
                                   distr=c("exp","llo","lno", "wei","rps", "gam", "gomp", "gga"),method="hmc",
                                   opinion_type = "survival",
                                   param_expert = param_expert,
                                   pool_type = "linear pool",
                                   k = k1,
                                   times_expert = c(4,5)*12,
                                   iter = 1000,
                                   chains =n.chains,
                                   priors = list(lno = list(sigma_beta = c(5,5)),
                                                 gga = list(sigma_beta = c(1,1))),
                                   init = list(rps = lapply(rep(1, n.chains), init_fun_rps),
                                               lno = lapply(rep(1, n.chains), init_fun_lnorm)))


#Plot the distributions 

times_expert = c(4,5)*12
scale <- 1
pool_type_eval <- "linear pool"
df.linear <- subset(dfs_pool[[1]], ftype == pool_type_eval) %>% rename(y = x) %>% 
  mutate(x = times_expert[1] + fx*scale)
df.linear2 <- subset(dfs_pool[[2]], ftype == pool_type_eval) %>% rename(y = x) %>% 
  mutate(x = times_expert[2] + fx*scale)



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
                              y = seq(0,1, length.out =10)), aes(x = grp, y= y), colour = "grey", lwd=1.1, inherit.aes = F)




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
                                               iter = 1000,
                                               chains =n.chains,
                                               priors = list(lno = list(sigma_beta = c(5,5)),
                                                             gga = list(sigma_beta = c(1,1))),
                                               init = list(rps = lapply(rep(1, n.chains), init_fun_rps),
                                                           lno = lapply(rep(1, n.chains), init_fun_lnorm)))

plt.vague <- plot(m.all_expert_vague, add.km = T, t = 0:80)+
  theme_light()+
  scale_x_continuous(expand = c(0, 0), limits = c(0,NA), breaks=seq(0, 80, 5)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA), breaks=seq(0, 1, 0.05))+
  xlab("Time (months)")+
  ylab("")


### Outputs for Publication ----


DIC_comp <- data.frame(Model =names(m.all_expert_vague$models),
           DIC_expert = m.all_expert$model.fitting$dic,
           DIC_vague = m.all_expert_vague$model.fitting$dic) %>% mutate_if(is.numeric, round,digits = 2)

kable(DIC_comp[order(DIC_comp$DIC_expert),], format="latex")
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
models <- names(m.all_expert_vague$models)

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


# 
# undebug(compute_ICs_stan)
# compute_ICs_stan(m.all_expert_vague$models$Gompertz, "gom",m.all_expert_vague$misc$data.stan[[7]])
figure <- ggpubr::ggarrange(plt.expert, plt.vague, 
                    ncol = 2, nrow = 1,common.legend = TRUE, legend="bottom")
ggsave(paste0(pathway,"Survival functions with Expert Opinion and default priors.png"), width =  10, height = 6.66)


## Example 2: Opinion on the survival of the comparator arm ----

data2 <- survHE::data %>% rename(status = censored) %>% mutate(time2 = ifelse(time > 10, 10, time),
                                                              status2 = ifelse(time> 10, 0, status))


param_expert2 <- list()

param_expert2[[1]] <- data.frame(dist = "norm", wi = 1, param1 = 0.10, param2 = 0.01, df = 15)


survHE.data.model  <- expertsurv:::fit.models.expert(formula=Surv(time2,status2)~as.factor(arm),data=data2,
                                        distr=c("llo"),
                                        method="hmc",
                                        iter = 1000,
                                        opinion_type = "survival",
                                        id_St = 0, 
                                        times_expert = 20, 
                                        param_expert = param_expert2)


#debug(compute_ICs_stan)
#undebug(fit.models.expert)
plot.exam_1.expert <- plot(survHE.data.model, add.km = T, t = 0:30)+
  theme_light()+
  scale_x_continuous(expand = c(0, 0), limits = c(0,NA), breaks=seq(0, 80, 2)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA), breaks=seq(0, 1, 0.05))


survHE.data.model.vague2  <- expertsurv:::fit.models.expert(formula=Surv(time2,status2)~as.factor(arm),data=data2,
                                                     distr=c("llo"),
                                                     method="hmc",
                                                     iter = 1000,
                                                     opinion_type = "survival",
                                                     id_St = 0, 
                                                     times_expert = 20, 
                                                     param_expert = param_expert_vague)
plot.exam_1.mle<- plot(survHE.data.model.vague2, add.km = T, t = 0:30)+
  theme_light()+
  scale_x_continuous(expand = c(0, 0), limits = c(0,NA), breaks=seq(0, 80, 2)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA), breaks=seq(0, 1, 0.05))


figure <- ggpubr::ggarrange(plot.exam_1.expert, plot.exam_1.mle, 
                            ncol = 2, nrow = 1,common.legend = TRUE, legend="bottom")
(figure )
ggsave(paste0(pathway,"Survival functions Example 2.png"), width =  10, height = 6.66)



## Example 3: Opinion about survival difference ----

param_expert3 <- list()

param_expert3[[1]] <- data.frame(dist = "norm", wi = 1, param1 = 5, param2 = 0.2, df = 15)


survHE.data.model  <- expertsurv:::fit.models.expert(formula=Surv(time2,status2)~as.factor(arm),data=data2,
                                                     distr=c("gom"),
                                                     method="hmc",
                                                     iter = 1000,
                                                     opinion_type = "mean",
                                                     id_trt = 1, 
                                                     param_expert = param_expert3)

plot.exam_2.expert <- plot(survHE.data.model, add.km = T, t = 0:50)

ggsave(paste0(pathway,"Survival function Example 3.png"), width =  10, height = 6.66)



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

