library("survHE")
library("SHELF")
library("xlsx")
library("abind")
library("mixR")
library("flexsurv")
library("devtools")
library("rstan")
library("mixtools")
library("mixR")
library("pkgcond")
library("abind")
library("ggpubr")
# library("rstantools")
# library("rstan")
# library("StanHeaders")
# library("survival")
# library("flexsurv")
# library("expertsurv")
devtools::load_all()
#install
# surv.inp <- system.file("extdata", "survival.txt", package = "survHE")
# 
# nrisk.inp <- system.file("extdata", "nrisk.txt", package = "survHE")
# 
# digitise(surv.inp, nrisk.inp, km_output = paste0("C:/Users/phili/Desktop/","KMdata.txt"),
#          ipd_output  = paste0("C:/Users/phili/Desktop/","IPDdata.txt"))

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

cbind(lower, upper)


write.table(data.frame(Time = Survival.df$time,
                       Survival =Survival.df$surv),
            paste0("C:/Users/phili/Desktop/","survival.txt"),
            row.names=T,
            sep="\t")

write.table(data.frame(Interval = 1:length(lower),
                       time = times.risk[-length(n.risk)],
                       lower =lower, upper = upper,
                       nrisk =  n.risk[-length(n.risk)]),
            paste0("C:/Users/phili/Desktop/","nrisk.txt"),
            row.names=FALSE,
            sep="\t")

# write.table(data.frame(Interval = 1:length(lower),
#                        time = times.risk[-which(times.risk>tail(Survival.df$time, n= 1))],
#                        lower =lower, upper = upper,
#                        nrisk =  n.risk[-which(times.risk>tail(Survival.df$time, n= 1))]),
#             paste0("C:/Users/phili/Desktop/","nrisk.txt"),
#             row.names=FALSE,
#             sep="\t")

digitise(paste0("C:/Users/phili/Desktop/","survival.txt"),
         paste0("C:/Users/phili/Desktop/","nrisk.txt"),
         km_output = "C:/Users/phili/Desktop/KMdata.txt",
         ipd_output = "C:/Users/phili/Desktop/IPDdata.txt")

undebug(digitise)

digitized_IPD <- read.table("C:/Users/phili/Desktop/IPDdata.txt")

colnames(digitized_IPD) <- digitized_IPD[1,]

digitized_IPD <- data.frame(digitized_IPD[-1,])

digitized_IPD <- apply(digitized_IPD,2, as.numeric)
km.fit <- survfit(Surv(time, event)~1, data.frame(digitized_IPD))
cbind(km.fit$time,km.fit$surv)
plot(km.fit)
expert.prob.df <- read.xlsx(paste0(pathway, "Expert Surv Probabilities.xlsx"), sheetName = "Sheet1") %>% data.frame()


j <- q <- 1

res.mat <- matrix(nrow = nrow(expert.prob.df)/3, ncol = 3)
for(i in 1:(nrow(expert.prob.df)/3)){
        
        res.mat[j,]    <- expert.prob.df[q:(q+2), 2]   
        
        q <- q+3
        j <- j +1
        
}

expert.prob.df2 <- cbind(res.mat, rep(1:7, times = 4), rep(2:5, each = 7))
colnames(expert.prob.df2) <- c("0.995", "Mode","0.005", "Expert","Time")

#Exclude Expert 2
#expert.prob.df2 <- expert.prob.df2[which(expert.prob.df2[ , which(colnames(expert.prob.df2) == "Expert")] != 2),]

times_act <- 4:5

dfs_expert <- list()
plts_pool <- list()
dfs_pool <- list()
lower_bound <- 0
upper_bound <- 1

surv.vals <- seq(0,1,by= 0.01)

plot(surv.vals,  y = dt.scaled(surv.vals, 0.6504780, 0.01272621, df = 2), type = "l")

qt.scaled(0.995, 0.6504780, 0.01272621, df = 2)
fit.eval.list <- list()

for(i in 1:length(times_act)){

expert.prob.eval <- expert.prob.df2[expert.prob.df2[,5] == times_act[i],1:3]

expert.prob.eval <- expert.prob.eval[,c(3,2,1)]
fit.eval.list[[i]] <- expertsurv:::fitdist_mod(t(expert.prob.eval)[-2,],
        probs = c(0.005,  0.995), upper = upper_bound, lower = lower_bound ,mode = t(expert.prob.eval)[2,], 
        expertnames = c(1:7))
fit.eval <- fit.eval.list[[i]]
#Slight bug with this function (SHELF version)
#SHELF:::plotfit(fit.eval, lp = T)



plts_pool[[i]] <- expertsurv:::makePoolPlot(fit  = fit.eval,
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
                                       expertnames = paste("Expert ", 1:7))

dfs_pool[[i]] <-  plts_pool[[i]][["data"]]



best_fit_index  <- apply(fit.eval$ssq[,c("normal","t","gamma", "lognormal", "beta")], 1, which.min)
best_fit <- names(fit.eval$ssq[,c("normal","t","gamma", "lognormal", "beta")])[best_fit_index]

best_fit_loc <- sapply(best_fit, function(x){which(x  == names(fit.eval$ssq))})
fit.eval.dist  <- fit.eval[best_fit_loc]

pool.df_output <- matrix(nrow = length(best_fit_loc),ncol = 3)
colnames(pool.df_output) <- c("param1", "param2", "param3")

for(j in 1:length(best_fit_loc)){
        pool.df_output[j,1:length(fit.eval.dist[[j]][j,])] <-  as.numeric(as.vector(fit.eval.dist[[j]][j,]))
}
dfs_expert[[i]] <- data.frame(dist = best_fit, wi = 1, pool.df_output)
}


saveRDS(dfs_expert, file = "dfs_expert.rds")
dist_fit <- c()
pool.eval.stan <- list()
pool_type_eval <- "linear pool"

for(i in 1:length(dfs_expert)){
        
        pool.eval.stan[[i]]  <- expertsurv:::makePoolPlot.Data(pool.df = dfs_expert[[i]], 
                                                  pool_type = pool_type_eval,
                                                  add_hist = F, plt_other_dists = F, max_mix_eval = 6)
        
        
        pool.df <- subset(dfs_pool[[i]], ftype == pool_type_eval) %>% rename(Dist = ftype)
        pool.eval.stan[[i]]$plot.fit.mixnorm <- pool.eval.stan[[i]]$plot.fit.mixnorm+
               
                geom_line(data =pool.df,
                          aes(x = x, y = fx, colour = Dist),  lty = 2)+
                xlim(c(0,1))+
                ylim(c(0, max(pool.df$fx)*1.1))
              
        dist_fit <- c(dist_fit, names(pool.eval.stan[[i]])[2])
        if(i == 1){
                param_expert <-  pool.eval.stan[[i]][[2]]
                
        }else{
                param_expert <-  abind(param_expert,pool.eval.stan[[i]][[2]],along = 1)
                
        }
        
        rownames(param_expert) <- dist_fit
        
}



plts_pool[[1]]+ ggtitle("Best fitting Expert Survival Opinions Year 4")
ggsave("Year 4 Expert Opinion.png", width =  10, height = 10)


plts_pool[[2]]+ ggtitle("Best fitting Expert Survival Opinions Year 5")
ggsave("Year 5 Expert Opinion.png", width =  10, height = 10)



ggarrange(plts_pool[[1]]+
            theme(legend.title = element_blank()), plts_pool[[2]]+ ylab("")+
            theme(legend.title = element_blank()),  
          labels = c("Year 4", "Year 5"),
          ncol = 2, nrow = 1,common.legend = T,legend="bottom")

ggsave(filename = "Year 4 & 5 Distributions.png",
       width =  9, height = 4.5)

ggsave(plot =pool.eval.stan[[1]]$plot.fit.mixnorm+
       labs(caption = "Dotted line indicates actual linear pool, solid line is the mixture normal."),
       filename = "Year 4 Linear Pool.png",
       width =  5, height = 5)
       
       


tmpfun <- get("make_sim_hmc", envir = asNamespace("survHE"))
environment(make_sim_hmc) <- environment(tmpfun)
attributes(make_sim_hmc) <- attributes(tmpfun)  
assignInNamespace("make_sim_hmc", make_sim_hmc, ns="survHE")


digitized_IPD<- data.frame(digitized_IPD)
# m.gamma <- survHE::fit.models(formula=Surv(time,event)~1,data=digitized_IPD,
#                               distr="rps",method="hmc", k = 2)
# Fits a parametric model
#issue with print for gomp and rps
# m.gamma <- survHE::fit.models(formula=Surv(time,event)~1,data=digitized_IPD,
#                         distr=c("exp","llo","lno","rps", "wei","lnorm"),
#                         method="hmc", priors=list(gom=list(a_alpha=0.1,b_alpha=0.1)))

k <- rps.knot <- 1


mle.ests_rps <- flexsurvspline(Surv(time, event) ~ 1, data=digitized_IPD, k=rps.knot, scale="hazard")
mle.ests_lno <- flexsurvreg(Surv(time, event) ~ 1, data=digitized_IPD, dist="lno")

#init_fun <- function(...){list(gamma=apply(gamma_params,1,function(x){rnorm(1,mean = x[1], sd = x[2]/2)}))}

init_fun_rps <- function(...){list(gamma=as.numeric(mvtnorm::rmvnorm(n = 1, mean = mle.ests_rps$res[,1],
                                                                     sigma = mle.ests_rps$cov)))}


init_fun_lnorm <- function(...){
        res <- as.numeric(mvtnorm::rmvnorm(n = 1, mean = mle.ests_lno$res.t[,1], sigma = mle.ests_lno$cov))
        list(beta =c(res[1],0),
             alpha = res[2])
        
        }

n.chains <- 3


k1 <- 1

knots1 <- quantile(log((digitized_IPD %>% filter(event == 1))$time),
                   seq(0, 1, length = k1 + 2))

# Will give standard estimates (ie without prior)
# param_expert[,,2] <- 10
# param_expert[1,5,2] <- -999.2

m.all <- expertsurv:::expert_surv2(formula=Surv(time,event)~1,data=digitized_IPD,
                        distr=c("exp","llo","lno", "wei","rps", "gam", "gomp", "gga"),method="hmc",
                        opinion_type = "survival",
                        param_expert = param_expert,
                        k = k1,
                        times_expert = c(4,5)*12,
                        iter = 10000,
                        chains =n.chains,
                        priors = list(lno = list(sigma_beta = c(5,5)),
                                      gga = list(sigma_beta = c(1,1))),
                        init = list(rps = lapply(rep(1, n.chains), init_fun_rps),
                                    lno = lapply(rep(1, n.chains), init_fun_lnorm)))

m.all$models

saveRDS(m.all, file = "m_all.rds")


# library(shinystan)
# 
# shinystan::launch_shinystan(m.all$models$`Royston-Parmar`)



mle.ests_rps <- flexsurvspline(Surv(time, event) ~ 1, data=digitized_IPD, k=2, scale="hazard")

k2 <- 2

knots2 <- quantile(log((digitized_IPD %>% filter(event == 1))$time),
                   seq(0, 1, length = k2 + 2))

m.rps <- expertsurv:::expert_surv2(formula=Surv(time,event)~1,data=digitized_IPD,
                                   distr=c("rps"),method="hmc",
                                   opinion_type = "survival",
                                   param_expert = param_expert,
                                   k = k2,
                                   times_expert = c(4,5)*12,
                                   iter = 2000,
                                   chains =n.chains,
                                   priors = list(lno = list(sigma_beta = c(5,5))),
                                   init = list(rps = lapply(rep(1, n.chains), init_fun_rps)))

m.rps$models
plt <- survHE:::plot.survHE(m.all, add.km = T, t = 0:80, mods = 1:4)


surv.rps <- apply(data.frame(extract(m.all$models$`Royston-Parmar`, "gamma")), 1, FUN =
             function(x){flexsurv::psurvspline(0:80,x,knots = knots1, lower.tail = F )})


gen.gamma.sims <- m.all$models$`Gen. Gamma`$BUGSoutput$sims.matrix[,c("mu", "sigma","Q")]
surv.gengamma <- apply(gen.gamma.sims, 1, function(x){pgengamma(0:80, x[1], x[2], x[3], lower.tail = F)})

gomp.sims <- m.all$models$Gompertz$BUGSoutput$sims.matrix[,c("alpha", "rate")]
surv.gomp <- apply(gomp.sims, 1, function(x){pgompertz(0:80, x[1], x[2], lower.tail = F)})

gamma.sims <- m.all$models$Gamma$BUGSoutput$sims.matrix[,c("alpha", "rate")]
surv.gamma <- apply(gamma.sims, 1, function(x){pgamma(0:80, x[1], x[2], lower.tail = F)})


flex.weibull <- flexsurvreg(Surv(time, event)~1, data = digitized_IPD, dist = "weibull")
 
mle.weibull <- data.frame(summary(flex.weibull, t =1:100))
 mle.weibull$time <- as.numeric(mle.weibull$time)


g <- ggplot_build(plt)
unique_cols <- unique(g$data[[1]] %>% pull(colour))

times_expert = c(4,5)*12
scale <- 1
pool_type_eval <- "linear pool"
df.linear <- subset(dfs_pool[[1]], ftype == pool_type_eval) %>% rename(y = x) %>% 
        mutate(x = times_expert[1] + fx*scale)
df.linear2 <- subset(dfs_pool[[2]], ftype == pool_type_eval) %>% rename(y = x) %>% 
  mutate(x = times_expert[2] + fx*scale)
#install.packages("wesanderson")
#library(wesanderson)
#pal <- wes_palette("Zissou1", 100, type = "continuous")
final.plt_expert <- plt+ #geom_point(data =mle.weibull, aes(x = time, y = est), colour = "blue", inherit.aes = T)+
         geom_line(data = data.frame(time = 0:80, est = rowMeans(surv.rps)),
            aes(x = time, y = est), colour = "black")+
        geom_line(data = data.frame(time = 0:80, est = rowMeans(surv.gamma)),
                   aes(x = time, y = est), colour = "orange")+
        geom_line(data = data.frame(time = 0:80, est = rowMeans(surv.gomp)),
                   aes(x = time, y = est), colour = "brown")+
        geom_line(data = data.frame(time = 0:80, est = rowMeans(surv.gengamma)),
                  aes(x = time, y = est), colour = "royalblue2")+
        scale_colour_manual(name = 'Models', 
                            values =c('Exponential'=unique_cols[1],
                                      'log-Logistic'=unique_cols[2],
                                      "log-Normal"= unique_cols[3],
                                      "Weibull (AFT)"= unique_cols[4],
                                      "Royston-Parmar"= "black",
                                      "Gamma"= "orange",
                                      "Gompertz" = "brown", 
                                      "Gen Gamma"= "royalblue2"))+
        geom_path(data = df.linear, aes(x =x, y =y), colour = "grey", lwd=1.1, inherit.aes = F)+
        geom_path(data = data.frame(grp = rep(times_expert[1], 10), 
                                    y = seq(0,1, length.out =10)), aes(x = grp, y= y), colour = "grey", lwd=1.1, inherit.aes = F)+
        geom_polygon(data = df.linear, aes(y = y, x = x), fill = "sky blue", alpha = 0.5)+
        geom_path(data = df.linear2, aes(x =x, y =y), colour = "grey", lwd=1.1, inherit.aes = F)+
        geom_polygon(data = df.linear2, aes(y = y, x = x), fill = "sky blue", alpha = 0.5)+
        geom_path(data = data.frame(grp = rep(times_expert[2], 10), 
                              y = seq(0,1, length.out =10)), aes(x = grp, y= y), colour = "grey", lwd=1.1, inherit.aes = F)+
  
  theme(legend.position = "bottom")+
   theme_light()+
  scale_x_continuous(expand = c(0, 0), limits = c(0,NA), breaks=seq(0, 80, 5)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA), breaks=seq(0, 1, 0.05))+
  xlab("Time (months)")
ggsave("Expert Elicited Survival Functions.png", width =  10, height = 10)







plt+
geom_line(data = data.frame(time = 0:80, est = rowMeans(p)),
          aes(x = time, y = est), colour = "orange")

final.plt


#Calculate WAIC I can't actually use WAIC

digitized_IPD <- data.frame(digitized_IPD)
t <- digitized_IPD$time
event <- digitized_IPD$event
event2 <- ifelse(digitized_IPD$event ==1, 0,1)
mix_norm_log <- function(surv, x_df){
  log(sum(exp(apply(x_df,1,function(x){dnorm(surv, mean = x[1], sd = x[2], log = T)})+ log(x_df[,3]))))
  
}

df_param <- data.frame(extract(m.all$models$Exponential, "rate"))
LL <- apply(df_param,1, 
            function(x){dexp(t,rate = x[1],log = T)*event+
                        pexp(t,rate = x[1],log = T,lower.tail = F)*event2})

LL.expert <- apply(df_param,1, function(x){
  mix_norm_log(pexp(times_expert[1],rate = x,lower.tail = F),m.all$misc$data.stan[[1]]$param_expert[1,,])+
    mix_norm_log(pexp(times_expert[2],rate = x,lower.tail = F),m.all$misc$data.stan[[1]]$param_expert[2,,])})
  
#mean(colSums(LL))


# LL.2 <- apply(df_param,1, 
#             function(x){hexp(t,rate = x,log = T)*event+
#                 pexp(t,rate = x,log = T,lower.tail = F)})
#mean(colSums(LL.2))
#undebug(compute_ICs_stan)
#compute_ICs_stan(m.all$models$Exponential, "exp",m.all$misc$data.stan[[1]] )
undebug(survHE:::lik_exp)
#survHE:::lik_exp()

LL_hat <- dexp(t,rate = colMeans(df_param),log = T)*event+
  pexp(t,rate = colMeans(df_param),log = T,lower.tail = F)*event2

LL.expert_hat <- 
  mix_norm_log(pexp(times_expert[1],rate = colMeans(df_param),lower.tail = F),m.all$misc$data.stan[[1]]$param_expert[1,,])+
    mix_norm_log(pexp(times_expert[2],rate = colMeans(df_param),lower.tail = F),m.all$misc$data.stan[[1]]$param_expert[2,,])



pd_1 <- 2*(sum(LL_hat)+LL.expert_hat -mean(colSums(LL)+LL.expert)) 

DIC.exp <- -2*(sum(LL_hat)+LL.expert_hat) + 2*pd_1
#waic.exp <- loo::waic(t(LL))

# library(shinystan)
# 
# shinystan::launch_shinystan(m.all$models$`Weibull (AFT)`)
df_param <- data.frame(extract(m.all$models$`Weibull (AFT)`, c("scale", "alpha")))

LL <- apply(df_param,1, 
            function(x){dweibull(t,scale = x[1], shape = x[2], log = T)*event+
                            pweibull(t,scale = x[1], shape = x[2],log = T,lower.tail = F)*event2
            })

LL.expert <- apply(df_param,1, function(x){
  mix_norm_log(pweibull(times_expert[1],scale = x[1],shape = x[2],lower.tail = F),m.all$misc$data.stan[[1]]$param_expert[1,,])+
    mix_norm_log(pweibull(times_expert[2],scale = x[1],shape = x[2],lower.tail = F),m.all$misc$data.stan[[1]]$param_expert[2,,])})


LL.expert_hat <- 
  mix_norm_log(pweibull(times_expert[1],scale = colMeans(df_param)[1],shape = colMeans(df_param)[2],
                        lower.tail = F),m.all$misc$data.stan[[1]]$param_expert[1,,])+
  mix_norm_log(pweibull(times_expert[2],scale = colMeans(df_param)[1],shape = colMeans(df_param)[2],
                        lower.tail = F),m.all$misc$data.stan[[1]]$param_expert[2,,])


# 
# pd_1 <- 2*(sum(LL_hat)+LL.expert_hat -mean(colSums(LL)+LL.expert)) 
# 
# DIC.weib <- -2*(sum(LL_hat)+LL.expert_hat) + 2*pd_1


#waic.weib <- loo::waic(t(LL))

df_hat <- colMeans(df_param)
#df_hat <- apply(df_param, 2, median)

LL_hat <- dweibull(t,scale = df_hat[1], shape = df_hat[2], log = T)*event+
  pweibull(t,scale = df_hat[1], shape = df_hat[2],log = T,lower.tail = F)*event2

pd_1 <- 2*(sum(LL_hat) - mean(colSums(LL))) 

DIC.weib <- -2*sum(LL_hat) + 2*pd_1

# compute_ICs_stan(m.all$models$`Weibull (AFT)`, "wei",m.all$misc$data.stan[[4]] )
# 
# LL_hat <- dweibull(t,rate = colMeans(df_param),log = T)*event+
#   pexp(t,rate = colMeans(df_param),log = T,lower.tail = F)*event2
# 
# pd_1 <- 2*(sum(LL_hat) -mean(colSums(LL))) 
df_param <- data.frame(extract(m.all$models$`log-Normal`, c("meanlog", "alpha")))
LL <- apply(df_param,1, 
            function(x){dlnorm(t,meanlog = x[1], sdlog = x[2], log = T)*event+
                            plnorm(t,meanlog = x[1], sdlog = x[2],log = T,lower.tail = F)*event2    
            })
#waic.lnorm <- loo::waic(t(LL))


LL.expert <- apply(df_param,1, function(x){
  mix_norm_log(plnorm(times_expert[1],meanlog = x[1],sdlog = x[2],lower.tail = F),m.all$misc$data.stan[[1]]$param_expert[1,,])+
    mix_norm_log(plnorm(times_expert[2],meanlog = x[1],sdlog = x[2],lower.tail = F),m.all$misc$data.stan[[1]]$param_expert[2,,])})


LL_hat <- dlnorm(t,meanlog = colMeans(df_param)[1], sdlog = colMeans(df_param)[2], log = T)*event+
                plnorm(t,meanlog = colMeans(df_param)[1], sdlog = colMeans(df_param)[2],log = T,lower.tail = F)*event2    
          


LL.expert_hat <- 
  mix_norm_log(plnorm(times_expert[1],meanlog = colMeans(df_param)[1],sdlog = colMeans(df_param)[2],
                       lower.tail = F),m.all$misc$data.stan[[1]]$param_expert[1,,])+
  mix_norm_log(plnorm(times_expert[2],meanlog = colMeans(df_param)[1],sdlog = colMeans(df_param)[2],
                       lower.tail = F),m.all$misc$data.stan[[1]]$param_expert[2,,])

pd_1 <- 2*(sum(LL_hat)+LL.expert_hat -mean(colSums(LL)+LL.expert)) 
DIC.lnorm <- -2*(sum(LL_hat)+LL.expert_hat) + 2*pd_1




df_param <- data.frame(extract(m.all$models$`log-Logistic`, c("alpha", "rate")))
LL <- apply(df_param,1, 
            function(x){dllogis(t,shape = x[1], scale = x[2], log = T)*event+
                            pllogis(t,shape = x[1], scale = x[2],log = T,lower.tail = F)*event2})
#waic.llogis <- loo::waic(t(LL))
LL_hat <- dllogis(t,shape = colMeans(df_param)[1], scale = colMeans(df_param)[2], log = T)*event+
                pllogis(t,shape = colMeans(df_param)[1], scale = colMeans(df_param)[2],log = T,lower.tail = F)*event2


LL.expert <- apply(df_param,1, function(x){
  mix_norm_log(pllogis(times_expert[1],shape = x[1],scale = x[2],lower.tail = F),m.all$misc$data.stan[[1]]$param_expert[1,,])+
    mix_norm_log(pllogis(times_expert[2],shape = x[1],scale = x[2],lower.tail = F),m.all$misc$data.stan[[1]]$param_expert[2,,])})


LL.expert_hat <- 
  mix_norm_log(pllogis(times_expert[1],shape = colMeans(df_param)[1],scale = colMeans(df_param)[2],
                         lower.tail = F),m.all$misc$data.stan[[1]]$param_expert[1,,])+
  mix_norm_log(pllogis(times_expert[2],shape = colMeans(df_param)[1],scale = colMeans(df_param)[2],
                         lower.tail = F),m.all$misc$data.stan[[1]]$param_expert[2,,])

pd_1 <- 2*(sum(LL_hat)+LL.expert_hat -mean(colSums(LL)+LL.expert)) 
DIC.llogis <- -2*(sum(LL_hat)+LL.expert_hat) + 2*pd_1




df_param <- m.all$models$Gompertz$BUGSoutput$sims.matrix[,c("alpha", "rate")]
LL <- apply(df_param,1, 
            function(x){dgompertz(t,shape = x[1], rate = x[2], log = T)*event+
                            pgompertz(t,shape = x[1], rate = x[2],log = T,lower.tail = F)*event2    
            })

LL_hat <- dgompertz(t,shape = colMeans(df_param)[1], rate = colMeans(df_param)[2], log = T)*event+
                pgompertz(t,shape = colMeans(df_param)[1], rate = colMeans(df_param)[2],log = T,lower.tail = F)*event2    
           
#waic.gomp <- loo::waic(t(LL))

LL.expert <- apply(df_param,1, function(x){
  mix_norm_log(pgompertz(times_expert[1],shape = x[1],rate = x[2],lower.tail = F),m.all$misc$data.stan[[1]]$param_expert[1,,])+
    mix_norm_log(pgompertz(times_expert[2],shape = x[1],rate = x[2],lower.tail = F),m.all$misc$data.stan[[1]]$param_expert[2,,])})


LL.expert_hat <- 
  mix_norm_log(pgompertz(times_expert[1],shape = colMeans(df_param)[1],rate = colMeans(df_param)[2],
                        lower.tail = F),m.all$misc$data.stan[[1]]$param_expert[1,,])+
  mix_norm_log(pgompertz(times_expert[2],shape = colMeans(df_param)[1],rate = colMeans(df_param)[2],
                        lower.tail = F),m.all$misc$data.stan[[1]]$param_expert[2,,])

pd_1 <- 2*(sum(LL_hat)+LL.expert_hat -mean(colSums(LL)+LL.expert)) 
DIC.gomp <- -2*(sum(LL_hat)+LL.expert_hat) + 2*pd_1


df_param <- m.all$models$Gamma$BUGSoutput$sims.matrix[,c("alpha", "rate")]
LL <- apply(df_param,1, 
            function(x){dgamma(t,shape = x[1], rate = x[2], log = T)*event+
                            pgamma(t,shape = x[1], rate = x[2],log = T,lower.tail = F)*event2    
            })
#waic.gam <- loo::waic(t(LL))

LL_hat <- dgamma(t,shape = colMeans(df_param)[1], rate = colMeans(df_param)[2], log = T)*event+
                pgamma(t,shape = colMeans(df_param)[1], rate = colMeans(df_param)[2],log = T,lower.tail = F)*event2    
            
LL.expert <- apply(df_param,1, function(x){
  mix_norm_log(pgamma(times_expert[1],shape = x[1],rate = x[2],lower.tail = F),m.all$misc$data.stan[[1]]$param_expert[1,,])+
    mix_norm_log(pgamma(times_expert[2],shape = x[1],rate = x[2],lower.tail = F),m.all$misc$data.stan[[1]]$param_expert[2,,])})


LL.expert_hat <- 
  mix_norm_log(pgamma(times_expert[1],shape = colMeans(df_param)[1],rate = colMeans(df_param)[2],
                         lower.tail = F),m.all$misc$data.stan[[1]]$param_expert[1,,])+
  mix_norm_log(pgamma(times_expert[2],shape = colMeans(df_param)[1],rate = colMeans(df_param)[2],
                         lower.tail = F),m.all$misc$data.stan[[1]]$param_expert[2,,])

pd_1 <- 2*(sum(LL_hat)+LL.expert_hat -mean(colSums(LL)+LL.expert)) 
DIC.gam <- -2*(sum(LL_hat)+LL.expert_hat) + 2*pd_1



df_param <- m.all$models$`Gen. Gamma`$BUGSoutput$sims.matrix[,c("mu", "sigma", "Q")]

LL <- apply(df_param,1, 
            function(x){dgengamma(t,mu = x[1], sigma = x[2],Q = x[3], log = T)*event+
                            pgengamma(t,mu = x[1], sigma = x[2],Q = x[3],log = T,lower.tail = F)*event2    
            })
#waic.gengamma <- loo::waic(t(LL))

LL_hat <- dgengamma(t,mu = colMeans(df_param)[1], sigma = colMeans(df_param)[2],Q = colMeans(df_param)[3], log = T)*event+
                pgengamma(t,mu = colMeans(df_param)[1], sigma = colMeans(df_param)[2],Q = colMeans(df_param)[3],log = T,lower.tail = F)*event2    
           


LL.expert <- apply(df_param,1, function(x){
  mix_norm_log(pgengamma(times_expert[1],mu = x[1],sigma = x[2],Q = x[3],lower.tail = F),m.all$misc$data.stan[[1]]$param_expert[1,,])+
    mix_norm_log(pgengamma(times_expert[2],mu = x[1],sigma = x[2],Q = x[3],lower.tail = F),m.all$misc$data.stan[[1]]$param_expert[2,,])})


LL.expert_hat <- 
  mix_norm_log(pgengamma(times_expert[1],
                         mu = colMeans(df_param)[1],
                         sigma = colMeans(df_param)[2],
                         Q = colMeans(df_param)[3],
                      lower.tail = F),m.all$misc$data.stan[[1]]$param_expert[1,,])+
  mix_norm_log(pgengamma(times_expert[2],
                         mu = colMeans(df_param)[1],
                         sigma = colMeans(df_param)[2],
                         Q = colMeans(df_param)[3],
                      lower.tail = F),m.all$misc$data.stan[[1]]$param_expert[2,,])

pd_1 <- 2*(sum(LL_hat)+LL.expert_hat -mean(colSums(LL)+LL.expert)) 
DIC.gengamma <- -2*(sum(LL_hat)+LL.expert_hat) + 2*pd_1




df_param <- data.frame(extract(m.all$models$`Royston-Parmar`, "gamma"))
LL <- apply(df_param,1, 
      function(x){dsurvspline(t,gamma = x,knots = knots1,log = T)*event+
                  psurvspline(t,gamma = x,knots = knots1,log = T,lower.tail = F)*event2    
                })

LL_hat <- dsurvspline(t,gamma = colMeans(df_param),knots = knots1,log = T)*event+
                psurvspline(t,gamma = colMeans(df_param),knots = knots1,log = T,lower.tail = F)*event2    
         

LL.expert <- apply(df_param,1, function(x){
  mix_norm_log( psurvspline(times_expert[1],gamma = x,knots = knots1,lower.tail = F),m.all$misc$data.stan[[1]]$param_expert[1,,])+
    mix_norm_log(psurvspline(times_expert[2],gamma = x,knots = knots1,lower.tail = F),m.all$misc$data.stan[[1]]$param_expert[2,,])})


LL.expert_hat <- 
  mix_norm_log( psurvspline(times_expert[1],gamma = colMeans(df_param),knots = knots1,lower.tail = F),m.all$misc$data.stan[[1]]$param_expert[1,,])+
    mix_norm_log(psurvspline(times_expert[2],gamma = colMeans(df_param),knots = knots1,lower.tail = F),m.all$misc$data.stan[[1]]$param_expert[2,,])


pd_1 <- 2*(sum(LL_hat)+LL.expert_hat -mean(colSums(LL)+LL.expert)) 
DIC.rps <- -2*(sum(LL_hat)+LL.expert_hat) + 2*pd_1



#waic.rps <- loo::waic(t(LL))




#waic.rps2 <- loo::waic(t(LL))


waic = c(waic.exp$estimates[3,1],
         waic.llogis$estimates[3,1],
         waic.lnorm$estimates[3,1],
         waic.weib$estimates[3,1],
         waic.rps$estimates[3,1],
         waic.gam$estimates[3,1],
         waic.gomp$estimates[3,1],
         waic.gengamma$estimates[3,1])

DIC = c(DIC.exp,
         DIC.llogis,
         DIC.lnorm,
         DIC.weib,
         DIC.rps,
         DIC.gam,
         DIC.gomp,
         DIC.gengamma)



fit.df <- data.frame(Model = m.all$misc$model_name, m.all$model.fitting, DIC = DIC)

fit.df <- fit.df[order(fit.df[,"DIC"]),]


Model_Names <- rep(NA, nrow(fit.df))

model_full_names <- load_availables()$hmc

names(model_full_names)[which(names(model_full_names)=="RP")] <- "Royston-Parmar (1-knot)"
names(model_full_names)[which(names(model_full_names)=="WeibullAF")] <- "Weibull ATF"
names(model_full_names)[which(names(model_full_names)=="logLogistic")] <- "Log-Logistic"
names(model_full_names)[which(names(model_full_names)=="logNormal")] <- "Log-Normal"


for( i in 1:nrow(fit.df)){
  Model_Names[i] <-   names(model_full_names)[grep(fit.df[i, "Model"], model_full_names)]
}

fit.df.final <- data.frame(Models = Model_Names, DIC = fit.df[, "DIC"])

fit.df.final$DIC <- round(fit.df.final$DIC, digits = 2)

library(kableExtra)

fit.df.final %>% kable(format = 'latex', caption = "Survival models ordered by DIC (lower is better)")%>%
  kable_classic_2(full_width = F)


    
res.array <- survHE:::load_availables()$mle[which(survHE:::load_availables()$mle %in% m.rps[["misc"]][["model_name"]])]
res.array<- res.array[match( m.rps[["misc"]][["model_name"]],res.array)]

names(res.array)

#When k is not added the print is weird!
# Issue with plotting survival....
# DIC for RPS looks dodgy


undebug(make_data_stan)



print.survHE(m.gamma)
undebug(print.survHE)
library(survHE)
m.gamma2 <- survHE::fit.models(formula=Surv(time,event)~1,data=digitized_IPD,
                              distr="gam",method="mle")

plot(m.gamma)
undebug(lik_gga)
undebug(survHE::print.survHE)
undebug(survHE:::format_output_fit.models)

undebug(print.survHE)



  model_code <- '
  functions {
    real stud_t(real x, real nu, real mu , real sigma) {
      return student_t_lpdf(x|nu, mu, sigma);
    }
  }
'
  
  library("rstan")
  expose_stan_functions(stanc(model_code = model_code))
 
  
mean <- 5
sd <- 1
#ncp <- 2 not needed!
df <- 1
x <- 5
stud_t(x, df, mean, sd)

stats::dt((x - mean)/sd, df, ncp = ncp, log = TRUE) -log(sd)

dt.scaled <- function (x, df, mean = 0, sd = 1, ncp, log = FALSE){
  if (!log) 
    stats::dt((x - mean)/sd, df, ncp = ncp, log = FALSE)/sd
  else stats::dt((x - mean)/sd, df, ncp = ncp, log = TRUE) - 
    log(sd)
}


dt.scaled(x, df, mean, sd , log = T)



