
# Expert Opinion R script ----

# library("survHE",lib.loc ="~/R-packages/")
# library("SHELF",lib.loc ="~/R-packages/")
library("xlsx")
library("abind")
library("crayon")
library("flexsurv")
library("survHE")
#library("expertsurv", lib.loc ="~/R-packages/")
library("kableExtra")
# remove.packages("expertsurv")


## Example Publication: ELIANA Trial Data ---- 

pathway <- "C:/Users/phili/OneDrive/PhD/R packages/expertsurv"

Survival.df <- read.xlsx(paste0(pathway, "/data/ELIANA OS.xlsx"), 1) %>% data.frame()
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
            paste0(pathway,"/data/survival.txt"),
            row.names=T,
            sep="\t")

write.table(data.frame(Interval = 1:length(lower),
                       time = times.risk[-which(times.risk>tail(Survival.df$time, n= 1))],
                       lower =lower, upper = upper,
                       nrisk =  n.risk[-which(times.risk>tail(Survival.df$time, n= 1))]),
            paste0(pathway,"/data/nrisk.txt"),
            row.names=FALSE,
            sep="\t")

#digitize function exports it so it needs to be read back in 

digitise(paste0(pathway,"/data/survival.txt"),
         paste0(pathway,"/data/nrisk.txt"),
         km_output = paste0(pathway, "/data/KMdata.txt"),
         ipd_output = paste0(pathway, "/data/IPDdata.txt"))

digitized_IPD <- read.table(paste0(pathway,"/data/IPDdata.txt"))

colnames(digitized_IPD) <- digitized_IPD[1,]

digitized_IPD <- data.frame(digitized_IPD[-1,])

digitized_IPD <- apply(digitized_IPD,2, as.numeric)
km.fit <- survfit(Surv(time, event)~1, data.frame(digitized_IPD))

#Plot the Kaplan Meier curve
#plot(km.fit)

# Read in Expert opinions
expert.prob.df <- read.xlsx(paste0(pathway, "/data/Expert Surv Probabilities.xlsx"), sheetName = "Sheet1") %>% data.frame()

j <- q <- 1

res.mat <- matrix(nrow = nrow(expert.prob.df)/3, ncol = 3)
for(i in 1:(nrow(expert.prob.df)/3)){
  
  res.mat[j,]    <- expert.prob.df[q:(q+2), 2]   
  q <- q+3
  j <- j +1
  
}

expert.prob.df2 <- cbind(res.mat, rep(1:7, times = 4), rep(2:5, each = 7))
colnames(expert.prob.df2) <- c("0.995", "Mode","0.005", "Expert","Time")

times_act <- 4:5

dfs_expert <- list()
plts_pool <- list()
dfs_pool <- list()
fit.eval <- list()
lower_bound <- 0
upper_bound <- 1

time_eval <- c(24, max(km_surv$time))

km_surv<- summary(km.fit, t = time_eval)
km_surv_df <- data.frame(time = km_surv$time, surv = km_surv$surv, lower = km_surv$lower, upper= km_surv$upper)

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
    #scale_color_brewer(palette = "Paired")+
    scale_x_continuous(expand = c(0, 0), limits = c(0,NA)) + 
    scale_y_continuous(expand = c(0, 0), limits = c(0, 30))+
    xlab("S(t)")+
    #geom_errorbar(data = km_surv_df %>% filter(time == time_eval[i]),aes(xmin=lower, xmax=upper,y = 1),inherit.aes = F, col = "purple")+
    theme(legend.title=element_blank(),plot.margin = margin(0.75,0.1,0.1,0.1, "cm"))
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

cred_int_val <- cred_int(plts_pool[[1]],val = "linear pool", interval = c(0.025, 0.975))

ggpubr::ggarrange(plts_pool[[1]]+
                    theme(legend.title = element_blank()), plts_pool[[2]]+ ylab("")+
                    theme(legend.title = element_blank()),  
                  labels = c("Year 4", "Year 5"),
                  vjust =1,
                  ncol = 2, nrow = 1,common.legend = T,legend="bottom")



ggsave(filename = paste0(pathway, "plots/Year 2 & 3 Distributions.png"),
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

n.iter <- 10000

mle.ests_rps <- flexsurvspline(Surv(time, event) ~ 1, data=digitized_IPD, k=k1, scale="hazard")

init_fun_rps <- function(...){list(gamma=as.numeric(mvtnorm::rmvnorm(n = 1, mean = mle.ests_rps$res[,1],
                                                                     sigma = mle.ests_rps$cov)))}


m.all_expert <- expertsurv:::fit.models.expert(formula=Surv(time,event)~1,data=digitized_IPD,
                                               distr=c("exp","llo","lno", "wei","rps", "gam", "gomp", "gga"),method="hmc",
                                               opinion_type = "survival",
                                               param_expert = param_expert,
                                               pool_type = "linear pool",
                                               k = k1,
                                               times_expert = c(4,5)*12,
                                               iter = n.iter,
                                               chains =n.chains,
                                               init = list(rps = lapply(rep(1, n.chains), init_fun_rps)))

expertsurv:::model.fit.plot(m.all_expert, type = "dic", xlim = c(270, NA))
models <- names(m.all_expert$models)
psa_outuput <- list()

for(i in 1:length(m.all_expert$models)){
  psa <- survHE:::make.surv(fit = m.all_expert,mod = i, nsim = n.iter, t = seq(0, 80, by = 1))
  #psa.plot(psa)
  df_temp  <- t(apply(psa$mat[[1]], 1,quantile, probs = c(0.025, 0.5,.975))) %>% data.frame()
  df_temp$time <- seq(0, 80, by = 1)
  mod_name <- names(m.all_expert$models)[i]
  psa_outuput[[mod_name]] <- df_temp %>% mutate(model = mod_name)
  
}

df_final_expert <- do.call(rbind.data.frame, psa_outuput)

km.fit   <- survfit(Surv(time,event)~1,data=digitized_IPD)


datakm <- data.frame(time =km.fit$time, survival=  km.fit$surv, lower = km.fit$lower, upper = km.fit$upper)
times_expert = c(4,5)*12
scale <- 1
pool_type_eval <- "linear pool"
df.linear <- subset(dfs_pool[[1]], ftype == pool_type_eval) %>% rename(y = x) %>% 
  mutate(x = times_expert[1] + fx*scale)
df.linear2 <- subset(dfs_pool[[2]], ftype == pool_type_eval) %>% rename(y = x) %>% 
  mutate(x = times_expert[2] + fx*scale)


df_final_expert$Models <- factor(df_final_expert$model, levels = unique(df_final_expert$model))

AUC_res_main <-  matrix(NA, nrow = length(models), ncol = 1)
for(j in 1:length(models)){
  AUC_res_main[j,1] <- mean(summary(m.all_expert, mod = j, t = seq(0, 80, by = 1))[["mean.surv"]])
}



param_expert_vague <- list()

param_expert_vague[[1]] <- data.frame(dist = "beta", wi = 1, param1 = 1, param2 = 1, param3 = NA)


m.all_expert_vague <- expertsurv:::fit.models.expert(formula=Surv(time,event)~1,data=digitized_IPD,
                                                     distr=c("exp","llo","lno", "wei","rps", "gam", "gomp", "gga"),method="hmc",
                                                     opinion_type = "survival",
                                                     param_expert = param_expert_vague,
                                                     pool_type = "linear pool",
                                                     k = k1,
                                                     times_expert = c(4)*12,
                                                     iter = 10000,
                                                     chains =n.chains,
                                                     thin = 1,
                                                     init = list(rps = lapply(rep(1, n.chains), init_fun_rps)))



models <- names(m.all_expert_vague$models)
psa_outuput <- list()

for(i in 1:length(m.all_expert_vague$models)){
  psa <- make.surv(fit = m.all_expert_vague,mod = i, nsim = n.iter/2, t = seq(0, 80, by = 1))
  #psa.plot(psa)
  df_temp  <- t(apply(psa$mat[[1]], 1,quantile, probs = c(0.025, 0.5,.975))) %>% data.frame()
  df_temp$time <- seq(0, 80, by = 1)
  mod_name <- names(m.all_expert_vague$models)[i]
  psa_outuput[[mod_name]] <- df_temp %>% mutate(model = mod_name)
  
}

df_final_vague <- do.call(rbind.data.frame, psa_outuput)
df_final_vague$Models <- factor(df_final_vague$model, levels = unique(df_final_vague$model))



surv60_expert<- df_final_expert %>% filter(time ==60) %>%  mutate(CredInt = X97.5.- X2.5.) %>% group_by(Models) %>% summarize(CredInt = mean(CredInt))
surv60_vague <- df_final_vague %>% filter(time ==60) %>%  mutate(CredInt_vague = X97.5.- X2.5.) %>% group_by(Models) %>% summarize(CredInt_vague = mean(CredInt_vague))
surv_60_mods <- surv60_expert %>% dplyr::select(Models, CredInt) %>% left_join(surv60_vague, by = "Models")%>% mutate(diff=CredInt_vague -CredInt ) %>% arrange(desc(diff))

mod_select <- surv_60_mods[1:3,"Models"] %>%pull(Models) %>% as.character()

plt.expert <- ggplot(data = df_final_expert, aes(y = X50., x = time, group = Models, colour = Models))+
  geom_line()+
  geom_line(data = df_final_expert %>% filter(Models %in% mod_select),aes(y = X97.5., x = time),linetype="dotdash")+
  geom_line(data = df_final_expert %>% filter(Models %in% mod_select),aes(y = X2.5., x = time), linetype="dotdash")+
  geom_line(data = df_final_expert ,aes(y = X97.5., x = time),linetype="dotdash")+
  geom_line(data = df_final_expert ,aes(y = X2.5., x = time), linetype="dotdash")+
  
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


plt.vague <- ggplot(data = df_final_vague, aes(y = X50., x = time, group = Models, colour = Models))+
  geom_line()+
  geom_line(data = df_final_vague ,aes(y = X97.5., x = time),linetype="dotdash")+
  geom_line(data = df_final_vague,aes(y = X2.5., x = time), linetype="dotdash")+
  # geom_errorbar(data = df_final%>% filter(time == 80),aes(ymin=X2.5., ymax=X97.5.), width=8,
  #               position=position_dodge(2))+
  theme_light()+
  scale_x_continuous(expand = c(0, 0), limits = c(0,NA), breaks=seq(0, 80, 5)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA), breaks=seq(0, 1, 0.05))+
  xlab("Time (months)")+
  ylab("Survival")+
  geom_step(data = datakm, aes(x = time , y = survival), color = "black",inherit.aes = F) + 
  geom_ribbon(data = datakm, aes(x = time, y = survival, ymin = lower, 
                                 ymax = upper), alpha = 0.2, color = "grey",inherit.aes = F)




### Outputs for Publication ----


# figure <- ggpubr::ggarrange(plt.vague, plt.expert,
#                             ncol = 2, nrow = 1,common.legend = TRUE, legend="bottom")
# ggsave(paste0(pathway,"Survival functions with Expert Opinion and default priors.png"), width =  10, height = 6.66)
figure <- ggpubr::ggarrange(plt.vague, plt.expert,
                            ncol = 2, nrow = 1,common.legend = TRUE, legend="bottom")
ggsave(paste0(pathway,"Survival functions with Expert Opinion and default priors - without Cred.png"), width =  10, height = 6.66)
ggsave(paste0(pathway,"Survival functions with Expert Opinion and default priors - all Cred.png"), width =  10, height = 6.66)


DIC_comp <- data.frame(Model =names(m.all_expert$models),
                       DIC_expert = m.all_expert$model.fitting$dic,
                       DIC_vague = m.all_expert_vague$model.fitting$dic) %>% arrange(DIC_expert) %>% mutate_if(is.numeric, round,digits = 2)


DIC_comp <- data.frame(Model =names(m.all_expert_vague$models),
                       DIC_expert = m.all_expert$model.fitting$dic,
                       DIC_vague = m.all_expert_vague$model.fitting$dic) %>% mutate_if(is.numeric, round,digits = 2)



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

