library(truncnorm)
library(magrittr)
library(tidyverse)
library(latex2exp)

wi <- 0.5
space <- seq(0,1,by =0.001)

mu_1 <- 0.5
mu_2 <- 0.4
sd_1 <- 0.05
sd_2 <- 0.05

d_mu0<- dtruncnorm(space,a = 0,b = 1, mean = mu_1, sd = sd_1)
#sfsmisc::integrate.xy(space,d_mu0)
d_mu1<- dtruncnorm(space,a = 0,b = 1, mean = mu_2, sd = sd_2)


#Linear Pooling
d_linear <- (d_mu1*wi)+(d_mu0*wi)
#k <- sfsmisc::integrate.xy(space,d_comb) #Integrates to 1



#Logarithm Pooling
d_log <- (d_mu1^wi)*(d_mu0^wi)
k <- sfsmisc::integrate.xy(space,d_log)
d_log <- d_log/k


# "Bayesian Updating" - Likelihood * Prior  (Normal-Normal) : Same as the product between two Normal pdfs
pooled_mu <- (mu_1+mu_2)/2
pooled_sd <- sqrt((1/(sd_1^2)+1/(sd_2^2))^-1)

which(d_mu1==d_mu0)

expert.df <- data.frame(St = space, A = d_mu0, B = d_mu1,d_linear,d_log)

gathered.df<- tidyr::gather(expert.df, key = "Expert", value = "pdf", - St ) %>% 
  mutate(Expert = ifelse(Expert == "d_linear", "Linear pool",ifelse(Expert == "d_log","Log pool",Expert) ))

ggplot(gathered.df,aes(y = pdf, x = St,group = Expert))+
  geom_line(aes(color=Expert,linetype = Expert),size = 1)+
  scale_linetype_manual(values=c("solid","solid","twodash", "longdash"))+
  theme(legend.position="top")+
  theme_bw()
  
ggsave("Expert Opinion.png",width = 10, height = 5)  
  
  

plot(y = d_mu0,x = space, type = "l", col ="red", ylim = c(0,max(d_mu0,d_mu1,d_final)))
lines(y = d_mu1,x = space,  col ="blue")
lines(y = d_linear, x = space, col = "orange")
#lines(y = dnorm(space,pooled_mu,pooled_sd),x = space,col = "green")
lines(y = d_log, x = space, col = "pink")

#Posterior pooling
prod_dnorm <- dnorm(space, mean = mu_2, sd = sd_2)*dnorm(space, mean = mu_1, sd = sd_1)
k2<- sfsmisc::integrate.xy(space,prod_dnorm)
plot(y =prod_dnorm/k2, x = space)





### Gamma example

seq_gamma <- seq(0,5,0.01)

a1 <- 8
b1 <- 10
a2 <- 20
b2 <- 10
wi <- 0.5

# Prior


gamma1_prior <- dgamma(seq_gamma,a1,b1) #Prior 1
gamma2_prior <-  dgamma(seq_gamma,a2,b2)#Prior 2

# Log Pooling

gamma_comb_log_prior <- dgamma(seq_gamma,wi*(a1+a2),wi*(b1+b2)) 

#Log Pooling manually
d_log <- (gamma1_prior^wi)*(gamma2_prior^wi)
k <- sfsmisc::integrate.xy(seq_gamma,d_log)
d_log <- d_log/k

# Linear Pool

gamma_comb_lin_prior <- wi*(gamma1_prior+gamma2_prior)


png(file="Expert Opinion Prior.png",
    width=600, height=350)
#Both types of pooling are the same 
plot(y = gamma_comb_log_prior, seq_gamma, type = "l" , ylim = c(0, 2), col = "purple", 
     xlab = TeX('$\\lambda$'),ylab = "Density", xlim = c(0,4))
#points(y = d_log, seq_gamma,  col = "blue", cex = 0.02)
lines(y = gamma1_prior, seq_gamma,  col = "blue", lty = 2)
lines(y = gamma2_prior, seq_gamma,  col = "red", lty = 2)
lines(y = gamma_comb_lin_prior, seq_gamma, col = "green")
legend("topright", title = "Priors:", legend = c("Logarithmic pool", "Linear pool", 
                                                     paste0("Expert 1 G(8,10)"), "Expert 2 G(20,10)"),
       col = c("purple", "green", "blue", "red"), lty = c(1,1,2,2))
dev.off()

#Data
cum_time <- 7
n <- 10

#Scaled Likelihood 
gamma_grid <- (seq_gamma^n)*exp(-seq_gamma*cum_time)
k_like <- sfsmisc::integrate.xy(seq_gamma,gamma_grid)
gamma_grid_norm <- gamma_grid/k_like


# Posterior Log Pooling
gamma_post <- dgamma(seq_gamma,n+wi*(a1+a2),cum_time+wi*(b1+b2))


#Validation!
# d_log <- (dgamma(seq_gamma,n+a2,cum_time+b2)^wi)*(dgamma(seq_gamma,n+a1,cum_time+b1)^wi)
# k <- sfsmisc::integrate.xy(seq_gamma,d_log)
# d_log <- d_log/k


# Posterior Linear Pooling - Have to evaluate the density manually.
# And find the normalization constant.

# Approach 1: Pool the prior first 

gamma_den_unorm <- gamma_grid*wi*(gamma1_prior+gamma2_prior)
k2 <- sfsmisc::integrate.xy(seq_gamma,gamma_den_unorm)
gamma_den_norm <- gamma_den_unorm/k2


#Posterior
# Approach 2: Evaluate the individual posteriors first 
gamma1_post <- dgamma(seq_gamma,a1+n,b1+ cum_time)
gamma2_post <-  dgamma(seq_gamma,a2+n,b2 + cum_time)
gamma_den_unorm_post <- gamma_grid*wi*(gamma1_post+gamma2_post)
k3 <- sfsmisc::integrate.xy(seq_gamma,gamma_den_unorm_post)
gamma_den_norm_post <- gamma_den_unorm_post/k3




#Posterior Plots
png(file="Expert Opinion Posterior.png",
width=600, height=350)
plot(y = gamma_post, seq_gamma, type = "l", 
     xlab = TeX('$\\lambda$'),ylab = "Density",xlim = c(0,4),ylim = c(0,2)) # Log Pooling
#lines(y = d_log, seq_gamma, col = "green")
#lines(y = gamma_grid_norm, seq_gamma, col = "red", cex= 0.25) #Scaled likelihood
lines(y = gamma_den_norm,seq_gamma, col = "brown" )
lines(y = gamma_den_norm_post,seq_gamma, col = "orange" )
legend("topright", title = "Posterior:", legend = c("Logarithmic pool",
                                                    "Linear pool on prior densities",
                                                    "Linear pool on individual posterior densities"),
        col = c("black", "brown", "orange"), lty = 1)
dev.off()


dens_eval <- gamma_comb_lin_prior
dens_eval <- gamma_comb_log_prior
total_integral <- sfsmisc::integrate.xy(seq_gamma, dens_eval)


partial_integral <- rep(NA, length(dens_eval))
partial_integral[1] <- 0
for(i in 2:length(dens_eval)){
  partial_integral[i] <- sfsmisc::integrate.xy(seq_gamma[1:i], dens_eval[1:i])/total_integral
}
interval <- c(0.025, 0.975)
limits <- c(seq_gamma[which.min(abs(partial_integral - interval[1]))],
            seq_gamma[which.min(abs(partial_integral - interval[2]))])
names(limits) <- c("lower", "upper")
print(limits)
qgamma(c(0.025, 0.975), 14, 10)


