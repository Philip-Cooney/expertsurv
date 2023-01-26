
data.example = data.frame( times = flexsurv::rgompertz(100, 1.2,.9),
                           status = 1)
data.example$status[1:5] <-0
# Fits a parametric model

m.weibull <- survHE::fit.models(formula=Surv(times,status)~1,data=data.example,
                                distr="weibull",method="hmc")
gompertz.jags <- "
data{
for(i in 1:n){
zeros[i] <- 0
}
for(i in 1:n_time_expert){
zero2[i] <- 0
 }
}

model{


for(i in 1:n){
zeros[i] ~ dpois(zero.mean[i])

log_h[i] = log(mu[i]) + (alpha*t[i]);
log_S[i] = -mu[i]/alpha * (exp(alpha*t[i]) - 1);

LL[i] <- log_h[i]*d[i] + log_S[i] + log(a0[i])
zero.mean[i] <- -LL[i] + C

}

alpha ~ dgamma(a_alpha,b_alpha)

for(i in 1:H){
 prec_beta[i] <- 1/(sigma_beta[i])^2
 beta[i] ~ dnorm(mu_beta[i],prec_beta[i])
}

linpred <- X %*%beta
for( i in 1:n){
mu[i] <- exp(linpred[i])
}


for (i in 1:n_time_expert){
  zero2[i] ~ dpois(phi_con[i])
 
 
    St_expert[i,1] = exp(-mu[i]/alpha * (exp(alpha*time_expert[i]) - 1))*St_indic
    mean_surv_trt[i] <- (1/alpha)*exp(mu[id_trt]/alpha)*(exp(loggam(s))*(1-pgamma(mu[id_trt]/alpha,s,1)))
    mean_surv_comp[i] <- (1/alpha)*exp(mu[id_comp]/alpha)*(exp(loggam(s))*(1-pgamma(mu[id_comp]/alpha,s,1)))
    St_expert[i,2] <- (mean_surv_trt[i] - mean_surv_comp[i])*(ifelse(St_indic == 1, 0,1))
    
    St_expert_final[i] <- sum(St_expert[i,])
    
    for(j in 1:num_param[i]){
     expert_dens_mix[i,j] <- dnorm(St_expert_final[i], param_expert[i,j,1],pow(param_expert[i,j,2],-2))*param_expert[i,j,3]
    
    }
    
    expert_dens[i,1] <-  sum(expert_dens_mix[i,1:num_param[i]])
    expert_dens[i,2] <-  dnorm(St_expert_final[i], param_expert[i,1,1],pow(param_expert[i,1,2],-2))
    expert_dens[i,3] <-  dt(St_expert_final[i],    param_expert[i,1,1],param_expert[i,1,2],3) 
    expert_dens[i,4] <- dgamma(St_expert_final[i], max(param_expert[i,1,1],0.001),param_expert[i,1,2])
    expert_dens[i,5] <- dlnorm(St_expert_final[i], param_expert[i,1,1],param_expert[i,1,2])
    expert_dens[i,6] <- dbeta(St_expert_final[i], max(param_expert[i,1,1], 0.01),param_expert[i,1,2])
 
    phi_con[i] <- -log(expert_dens[i,data_dist_ind[i]]) + C

 } 
 
rate = exp(beta[1])
C <- 10000
s <- 0.0001

}"

stan.data <- m.weibull[["misc"]][["data.stan"]][[1]]
num_param <- c()

for(i in 1:dim(param_expert)[1]){
  num_param[i] <-   sum(param_expert[i,,1] != -999.2)
}
lk_up_dist <- c("mixture","norm", "t", "gamma", "lnorm","beta")

data_dist_ind<- as.numeric(sapply(rownames(param_expert), function(x){which(x== lk_up_dist)}))

stan.data$num_param <- num_param;
stan.data$data_dist_ind <- data_dist_ind;
stan.data$max_param <- max(num_param)
stan.data$n_time_expert <- length(times_act)
stan.data$param_expert <-param_expert[,1:max(num_param),];
stan.data$time_expert <- c(4,5)
stan.data$St_indic <- 1
stan.data$id_trt <- 1 
stan.data$id_comp <- 1
stan.data$a0 <- rep(1,nrow(data.example))
rjags::load.module("mix")
rjags::list.modules()

library(rjags)


gompertz.mod <-R2jags::jags(model.file = textConnection(gompertz.jags),
                        data=stan.data,
                        n.chains=2,
                        parameters.to.save = c("alpha","beta", "rate"),
                        n.iter = 2000,
                        n.thin = 5,
                        n.burnin = 10,
                        jags.module = c("glm","dic","mix"))


gomp.mle <- flexsurvreg(Surv(times,status)~1, data = data.example, dist = "gompertz")
plot(gomp.mle)

m.gomp <- survHE::fit.models(formula=Surv(times,status)~1,data=data.example,
                                 distr="gomp",method="hmc", priors = list(gom=list(a_alpha=0.1,b_alpha=0.1)))

View(m.gomp)





#Generalized Gamma function
gen.gamma.jags <- "
data{
for(i in 1:n){
    zeros[i] <- 0
    d2[i] <- ifelse(d[i] == 1, 0,1)
  }
for(i in 1:n_time_expert){
    zero2[i] <- 0
 }
}

model{


for(i in 1:H){
 prec_beta[i] <- 1/(sigma_beta[i])^2
 beta[i] ~ dnorm(mu_beta[i],prec_beta[i])
}

linpred <- X %*%beta
for( i in 1:n){
lambda[i] <- exp(linpred[i])
}

for(i in 1:n){

log_d[i]  <- log(dgen.gamma(t[i],r,lambda[i],b))
log_S[i]  <- log(1-pgen.gamma(t[i],r,lambda[i],b))
LL[i] <- log_d[i]*d[i] + log_S[i]*d2[i] + log(a0[i])
zeros[i] ~ dpois(zero.mean[i])
zero.mean[i] <- -LL[i] + C

}

for (i in 1:n_time_expert){
  zero2[i] ~ dpois(phi_con[i])

    St_expert[i,1] =  (1-pgen.gamma(time_expert[i],r,lambda[i],b))*St_indic
    mean_surv_trt[i] <- exp(loggam((b*r +1)/b) - loggam(r))/lambda[2]
    mean_surv_comp[i] <- exp(loggam((b*r +1)/b) - loggam(r))/lambda[1]
    St_expert[i,2] <- (mean_surv_trt[i] - mean_surv_comp[i])*(ifelse(St_indic == 1, 0,1))
    
    St_expert_final[i] <- sum(St_expert[i,])
    
    for(j in 1:num_param[i]){
     expert_dens_mix[i,j] <- dnorm(St_expert_final[i], param_expert[i,j,1],pow(param_expert[i,j,2],-2))*param_expert[i,j,3]
    
    }
    
    expert_dens[i,1] <-  sum(expert_dens_mix[i,1:num_param[i]])
    expert_dens[i,2] <-  dnorm(St_expert_final[i], param_expert[i,1,1],pow(param_expert[i,1,2],-2))
    expert_dens[i,3] <-  dt(St_expert_final[i],    param_expert[i,1,1],param_expert[i,1,2],3) 
    expert_dens[i,4] <- dgamma(St_expert_final[i], max(param_expert[i,1,1],0.001),param_expert[i,1,2])
    expert_dens[i,5] <- dlnorm(St_expert_final[i], param_expert[i,1,1],param_expert[i,1,2])
    expert_dens[i,6] <- dbeta(St_expert_final[i], max(param_expert[i,1,1], 0.01),param_expert[i,1,2])
 
    phi_con[i] <- -log(expert_dens[i,data_dist_ind[i]]) + C

 } 


r ~ dgamma(a_alpha,b_alpha);
b ~ dgamma(a_alpha,b_alpha);

C <- 10000
sigma <- 1/(b*pow(r,0.5))
Q <- pow(r,-0.5)
mu <- -beta[1] + (log(r)/b)


}"

gengamma.mod <-R2jags::jags(model.file = textConnection(gen.gamma.jags),
                            data=stan.data,
                            n.chains=2,
                            parameters.to.save = c("mu","beta",  "sigma","Q"),
                            n.iter = 2000,
                            n.thin = 5,
                            n.burnin = 10,
                            jags.module = c("glm","dic","mix"))

gengamma.mle <- flexsurvreg(Surv(times,status)~1, data = data.example, dist = "gengamma")




#Gamma JAGS function
gamma.jags <- "
data{
for(i in 1:n){
    zeros[i] <- 0
    d2[i] <- ifelse(d[i] == 1, 0,1)
  }
for(i in 1:n_time_expert){
    zero2[i] <- 0
 }
}

model{

for(i in 1:H){
 prec_beta[i] <- 1/(sigma_beta[i])^2
 beta[i] ~ dnorm(mu_beta[i],prec_beta[i])
}

linpred <- X %*%beta
for( i in 1:n){
lambda[i] <- exp(linpred[i])
}

for(i in 1:n){

log_d[i]  <- log(dgamma(t[i],alpha,lambda[i]))
log_S[i]  <- log(1-pgamma(t[i],alpha, lambda[i]))
LL[i] <- log_d[i]*d[i] + log_S[i]*d2[i] + log(a0[i])
zeros[i] ~ dpois(zero.mean[i])
zero.mean[i] <- -LL[i] + C

}

for (i in 1:n_time_expert){
  zero2[i] ~ dpois(phi_con[i])

    St_expert[i,1] =  (1-pgamma(time_expert[i],alpha,lambda[i]))*St_indic
    mean_surv_trt[i] <- alpha/lambda[id_trt]
    mean_surv_comp[i] <- alpha/lambda[id_comp]
    St_expert[i,2] <- (mean_surv_trt[i] - mean_surv_comp[i])*(ifelse(St_indic == 1, 0,1))
    
    St_expert_final[i] <- sum(St_expert[i,])
    
    for(j in 1:num_param[i]){
     expert_dens_mix[i,j] <- dnorm(St_expert_final[i], param_expert[i,j,1],pow(param_expert[i,j,2],-2))*param_expert[i,j,3]
    
    }
    
    expert_dens[i,1] <-  sum(expert_dens_mix[i,1:num_param[i]])
    expert_dens[i,2] <-  dnorm(St_expert_final[i], param_expert[i,1,1],pow(param_expert[i,1,2],-2))
    expert_dens[i,3] <-  dt(St_expert_final[i],    param_expert[i,1,1],param_expert[i,1,2],3) 
    expert_dens[i,4] <- dgamma(St_expert_final[i], max(param_expert[i,1,1],0.001),param_expert[i,1,2])
    expert_dens[i,5] <- dlnorm(St_expert_final[i], param_expert[i,1,1],param_expert[i,1,2])
    expert_dens[i,6] <- dbeta(St_expert_final[i], max(param_expert[i,1,1], 0.01),param_expert[i,1,2])
 
    phi_con[i] <- -log(expert_dens[i,data_dist_ind[i]]) + C

 } 


alpha ~ dgamma(a_alpha,b_alpha);

rate <- exp(beta[1]);

C <- 10000

}"
gamma.mod <-R2jags::jags(model.file = textConnection(gamma.jags),
                            data=stan.data,
                            n.chains=2,
                            parameters.to.save = c("alpha","beta","rate"),
                            n.iter = 2000,
                            n.thin = 5,
                            n.burnin = 10,
                            jags.module = c("glm","dic","mix"))


m.gamma <- survHE::fit.models(formula=Surv(times,status)~1,data=data.example,
                                 distr="gamma",method="hmc")
