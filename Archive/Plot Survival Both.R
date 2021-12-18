m.all.mle <- survHE::fit.models(formula=Surv(time,event)~1,data=digitized_IPD,
                                   distr=c("exp","llo","lno", "wei","rps", "gam", "gomp", "gga"),method="mle")



plt.data <- survHE:::plot.survHE(m.all.mle, add.km = T, t = 0:80, mods = 1:8)

plt.data<- plt.data+
  scale_colour_manual(name = 'Models', 
                      values =c('Exponential'=unique_cols[1],
                                'log-Logistic'=unique_cols[2],
                                "log-Normal"= unique_cols[3],
                                "Weibull (AFT)"= unique_cols[4],
                                "Royston-Parmar"= "black",
                                "Gamma"= "orange",
                                "Gompertz" = "brown", 
                                "Gen Gamma"= "royalblue2"))+
  theme(legend.position = "bottom")+
  theme_light()+
  scale_x_continuous(expand = c(0, 0), limits = c(0,NA), breaks=seq(0, 80, 5)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA), breaks=seq(0, 1, 0.05))+
  xlab("Time (months)")+
  ylab("")


figure <- ggarrange(final.plt_expert, plt.data, 
                    ncol = 2, nrow = 1,common.legend = TRUE, legend="bottom")

ggsave("Expert Elicited Survival Functions and MLES.png", width =  10, height = 6.66)
m.all.mle$model.fitting


fit.df <- data.frame(Model = m.all$misc$model_name,  DIC = DIC, BIC = m.all.mle$model.fitting$bic)

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

fit.df.final <- data.frame(Models = Model_Names, DIC = fit.df[, "DIC"], BIC = fit.df$BIC)

fit.df.final$DIC <- round(fit.df.final$DIC, digits = 2)
fit.df.final$BIC <- round(fit.df.final$BIC, digits = 2)

library(kableExtra)

fit.df.final %>% kable(format = 'latex', caption = "Survival models ordered by DIC (lower is better)")%>%
  kable_classic_2(full_width = F)




cbind(c("exp","llo","lno", "wei","rps", "gam", "gomp", "gga"),m.all.mle$model.fitting )