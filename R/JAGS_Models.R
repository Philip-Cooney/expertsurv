Gompertz.jags <- "
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

alpha1 ~ dnorm(mu_beta[1],1/(sigma_beta[1])^2)
alpha2 ~ dgamma(a_alpha, b_alpha)
alpha <-   alpha1*ifelse(St_indic == 1, 1,0) +alpha2*ifelse(St_indic == 0, 1,0)


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
 
 
    St_expert_temp[i,1] = exp(-mu[id_St]/alpha * (exp(alpha*time_expert[i]) - 1))*St_indic
    mean_surv_trt[i] <- (1/alpha)*exp(mu[id_trt]/alpha)*(exp(loggam(s))*(1-pgamma(mu[id_trt]/max(alpha,0.000000001),s,1)))
    mean_surv_comp[i] <- (1/alpha)*exp(mu[id_comp]/alpha)*(exp(loggam(s))*(1-pgamma(mu[id_comp]/max(alpha,0.0000001),s,1)))
    St_expert_temp[i,2] <- (mean_surv_trt[i] - mean_surv_comp[i])*(ifelse(St_indic == 1, 0,1))
    
    St_expert[i] <- sum(St_expert_temp[i,])

    for(j in 1:n_experts[i]){
    expert_dens[j,1,i] <-  dnorm(St_expert[i], param_expert[j,3,i],pow(param_expert[j,4,i],-2))
    expert_dens[j,2,i] <-  dt(St_expert[i],    param_expert[j,3,i],pow(param_expert[j,4,i],-2),max(param_expert[j,5,i],1)) 
    expert_dens[j,3,i] <-  dgamma(St_expert[i], max(param_expert[j,3,i],0.001),param_expert[j,4,i])
    expert_dens[j,4,i] <-  dlnorm(St_expert[i], param_expert[j,3,i],param_expert[j,4,i])
    expert_dens[j,5,i] <-  dbeta(St_expert[i], max(param_expert[j,3,i], 0.01),param_expert[j,4,i])
    phi_temp[j,i] <- equals(pool_type,1)*(expert_dens[j,param_expert[j,1,i],i]*param_expert[j,2,i])+equals(pool_type,0)*(expert_dens[j,param_expert[j,1,i],i]^param_expert[j,2,i])
    }
 
  phi_con[i] <- -log(sum(phi_temp[,i]))*equals(pool_type,1) +  -log(prod(phi_temp[,i]))*equals(pool_type,0) + C #+test_St_expert[i]*C

 } 
 
rate = exp(beta[1])
C <- 10000
s <- 0.0001

}"

Gamma.jags <- "
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

    St_expert_temp[i,1] =  (1-pgamma(time_expert[i],alpha,lambda[id_St]))*St_indic
    mean_surv_trt[i] <- alpha/lambda[id_trt]
    mean_surv_comp[i] <- alpha/lambda[id_comp]
    St_expert_temp[i,2] <- (mean_surv_trt[i] - mean_surv_comp[i])*(ifelse(St_indic == 1, 0,1))
    
    St_expert[i] <- sum(St_expert_temp[i,])
    
    for(j in 1:n_experts[i]){
    expert_dens[j,1,i] <-  dnorm(St_expert[i], param_expert[j,3,i],pow(param_expert[j,4,i],-2))
    expert_dens[j,2,i] <-  dt(St_expert[i],    param_expert[j,3,i],pow(param_expert[j,4,i],-2),max(param_expert[j,5,i],1)) 
    expert_dens[j,3,i] <-  dgamma(St_expert[i], max(param_expert[j,3,i],0.001),param_expert[j,4,i])
    expert_dens[j,4,i] <-  dlnorm(St_expert[i], param_expert[j,3,i],param_expert[j,4,i])
    expert_dens[j,5,i] <-  dbeta(St_expert[i], max(param_expert[j,3,i], 0.01),param_expert[j,4,i])
    phi_temp[j,i] <- equals(pool_type,1)*(expert_dens[j,param_expert[j,1,i],i]*param_expert[j,2,i])+equals(pool_type,0)*(expert_dens[j,param_expert[j,1,i],i]^param_expert[j,2,i])
    }
 
 phi_con[i] <- -log(sum(phi_temp[,i]))*equals(pool_type,1) +  -log(prod(phi_temp[,i]))*equals(pool_type,0) + C #+test_St_expert[i]*C

 } 


alpha ~ dgamma(a_alpha,b_alpha);

rate <- exp(beta[1]);

C <- 10000

}"
GenGamma.jags <- "

data{

    for(i in 1:n){
    zero[i] <- 0
    }
}


model{


for(i in 1:H){
 prec_beta[i] <- 1/(sigma_beta[i])^2
 beta_jags[i] ~ dnorm(mu_beta[i],prec_beta[i])
 
}

linpred <- X %*%beta_jags

for( i in 1:n){
lambda[i] <- exp(linpred[i])
}


for(i in 1:n){
is.censored[i]~dinterval(t_jags[i],t_cen[i])
t_jags[i] ~ dgen.gamma(r,lambda[i],b)

}

for (i in 1:n_time_expert){
  zero[i] ~ dpois(phi_con[i])

    St_expert_temp[i,1] =  (1-pgen.gamma(time_expert[i],r,lambda[id_St],b))*St_indic
    mean_surv_trt[i] <- exp(loggam((b*r +1)/b) - loggam(r))/lambda[id_comp]
    mean_surv_comp[i] <- exp(loggam((b*r +1)/b) - loggam(r))/lambda[id_trt]
    St_expert_temp[i,2] <- (mean_surv_trt[i] - mean_surv_comp[i])*(ifelse(St_indic == 1, 0,1))
    
    St_expert[i] <- sum(St_expert_temp[i,])
    
    
    for(j in 1:n_experts[i]){
    expert_dens[j,1,i] <-  dnorm(St_expert[i], param_expert[j,3,i],pow(param_expert[j,4,i],-2))
    expert_dens[j,2,i] <-  dt(St_expert[i],    param_expert[j,3,i],pow(param_expert[j,4,i],-2),max(param_expert[j,5,i],1)) 
    expert_dens[j,3,i] <-  dgamma(St_expert[i], max(param_expert[j,3,i],0.001),param_expert[j,4,i])
    expert_dens[j,4,i] <-  dlnorm(St_expert[i], param_expert[j,3,i],param_expert[j,4,i])
    expert_dens[j,5,i] <-  dbeta(St_expert[i], max(param_expert[j,3,i], 0.01),param_expert[j,4,i])
    phi_temp[j,i] <- equals(pool_type,1)*(expert_dens[j,param_expert[j,1,i],i]*param_expert[j,2,i])+equals(pool_type,0)*(expert_dens[j,param_expert[j,1,i],i]^param_expert[j,2,i])
    }
 
 
  phi_con[i] <- -log(sum(phi_temp[,i]))*equals(pool_type,1) +  -log(prod(phi_temp[,i]))*equals(pool_type,0) + C #+test_St_expert[i]*C
 } 


r ~ dgamma(a_alpha,b_alpha);
b ~ dgamma(a_alpha,b_alpha);

C <- 10000
sigma <- 1/(b*pow(r,0.5))
Q <- pow(r,-0.5)
mu <- -beta_jags[1] + (log(r)/b)

beta[1] <- mu

for(i in 2:H){
beta[i] <- beta_jags[i]

}


}"
