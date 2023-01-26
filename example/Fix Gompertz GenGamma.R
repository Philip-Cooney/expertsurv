df_comb<- digitized_IPD

colnames(df_comb) <- c("time", "status", "id")
obs_events <- df_comb[which(df_comb$status==1),"time"]
obs_censors <- df_comb[which(df_comb$status==0),"time"]
km.df <- survfit(Surv(time, status)~1,data = df_comb)
plot(km.df)
data_new <- list()
df_jags <- df_comb[,c("time","status", "id")]
df_jags$t <- df_jags$time




is.na(df_jags$t)<-df_jags$status==0
df_jags$is.censored<-1-df_jags$status
df_jags$t.cen<-df_jags$time+df_jags$status


modelinits <- list(list(t = tinits1),
                   list(t = tinits2))
data_jags <- list(N = nrow(df_jags),
                  t.cen = df_jags$t.cen,
                  is.censored = df_jags$is.censored,
                  t = df_jags$t,
                  id = df_jags$id,
                  num_id = length(unique(df_jags$id)),
                  time = df_jags$time,
                  status = df_jags$status)



#Generalized Gamma function
gen.gamma_surv_mean <- "

#http://www.mas.ncl.ac.uk/~nmf16/teaching/mas8391/survival.pdf
data{
    for(i in 1:n_time_expert){
    zero[i] <- 0
    }
    for(i in 1:n){
        is.censored[i] <- ifelse(d[i]==0, 1, 0)
        t_jags[i] <- ifelse(is.censored[i] ==1, NA, t[i]) 
        t_cens[i] <- t[i]+d[i]
    }

}

model{
  for(i in 1:N){
is.censored[i]~dinterval(t_jags[i],t.cen[i])
t[i] ~ dgen.gamma(r,lambda[id[i]],b)
  }

for(i in 1:num_id){
lambda[i] ~ dunif(0,5)
}
r ~ dunif(0,5)
b ~ dunif(0,5)

sigma <- 1/(b*pow(r,0.5))
Q <- pow(r,-0.5)
mu <- -log(lambda) + (log(r)/b)


}"



#Generalized Gamma function
GenGamma.jags <- "

data{
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
is.censored[i]~dinterval(t_jags[i],t_cen[i])
t_jags[i] ~ dgen.gamma(r,lambda[i],b)

}

for (i in 1:n_time_expert){
  zero2[i] ~ dpois(phi_con[i])

    St_expert[i,1] =  (1-pgen.gamma(time_expert[i],r,lambda[id_St],b))*St_indic
    mean_surv_trt[i] <- exp(loggam((b*r +1)/b) - loggam(r))/lambda[id_comp]
    mean_surv_comp[i] <- exp(loggam((b*r +1)/b) - loggam(r))/lambda[id_trt]
    St_expert[i,2] <- (mean_surv_trt[i] - mean_surv_comp[i])*(ifelse(St_indic == 1, 0,1))
    
    St_expert_final[i] <- sum(St_expert[i,])
    
    for(j in 1:num_param[i]){
     expert_dens_mix[i,j] <- dnorm(St_expert_final[i], param_expert[i,j,1],pow(param_expert[i,j,2],-2))*param_expert[i,j,3]
    
    }
    
    expert_dens[i,1] <-  sum(expert_dens_mix[i,1:num_param[i]])
    expert_dens[i,2] <-  dnorm(St_expert_final[i], param_expert[i,1,1],pow(param_expert[i,1,2],-2))
    expert_dens[i,3] <-  dt(St_expert_final[i],    pow(param_expert[i,1,2],-2),param_expert[i,1,2],3) 
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




tinits1 <-data.stan$t + max(data.stan$t)
is.na(tinits1)<-data.stan$d ==1

data.stan$is.censored <- ifelse(data.stan$d==0, 1, 0)
data.stan$t_jags <- ifelse(data.stan$is.censored ==1, NA, data.stan$t) 
data.stan$t_cen <- data.stan$t+data.stan$d

modelinits <- function(){list(t = tinits1)}

gen.gamma.jags <-R2jags::jags(model.file = textConnection(GenGamma.jags),
                              data=data.stan,
                              inits=modelinits,
                              n.chains=2,
                              parameters.to.save = c("lambda","b","r", "mu", "Q",
                                                     "sigma"),
                              n.iter = 100, n.thin = 1)

#Generalized Gamma function
gen.gamma_surv_mean <- "
data{
    zeros <- 0
    
    for(i in 1:n){
    is.censored[i] <- ifelse(d[1] == 0, 1, 0)
    }
    
    for(i in 1:n_time_expert){
    zero[i] <- 0
    }
}

model{
  for(i in 1:N){
is.censored[i]~dinterval(t[i],t.cen[i])
t[i] ~ dgen.gamma(r,lambda[id[i]],b)
  }

for(i in 1:num_id){
lambda[i] ~ dunif(0,5)
}
r ~ dunif(0,5)
b ~ dunif(0,5)

sigma <- 1/(b*pow(r,0.5))
Q <- pow(r,-0.5)
mu <- -log(lambda) + (log(r)/b)


}"

#Convert BUGS parametrization to Flexsurvreg
lambda = 0.2
b = 0.7
r = 2
# Mean of the distribution with the flexsurv paratmerization
gamma((b*r+1)/b)/(lambda * gamma(r))

mu <- -log(lambda) + log(r)/b 
sigma <- 1/(b*sqrt(r))
Q <- sqrt(1/r)


mean_gengamma(mu,sigma,Q)

#r = k
#b = shape 
#a = scale = 1/lambda
modelinits[[1]]$r =1
modelinits[[1]]$b =1
modelinits[[1]]$lambda =1/mean(data.stan$t)

modelinits[[2]]$r =1
modelinits[[2]]$b =1
modelinits[[2]]$lambda =1/mean(data.stan$t)

  flexsurv:::parse.dist("gengamma.orig")
  flexsurvreg(Surv(time, event)~1, data = digitized_IPD, dist = "gengamma")
  flexsurvreg(Surv(time, event)~1, data = digitized_IPD, dist = "gengamma.orig")
gen.gamma.jags <-R2jags::jags(model.file = textConnection(gen.gamma_surv_mean),
                              data=data_jags,
                              inits=modelinits,
                              n.chains=2,
                              parameters.to.save = c("lambda","b","r", "mu", "Q",
                                                     "sigma"),
                              n.iter = 10000, n.thin = 1)

m.gamma2 <- survHE::fit.models(formula=Surv(time,event)~1,data=digitized_IPD,
                               distr="gga",method="mle")


#Generalized Gamma function
GenGamma.jags <- "
data{
for(i in 1:n){
    ones[i] <- 0
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
#zeros[i] ~ dpois(zero.mean[i])
#zero.mean[i] <- -LL[i] + C


p[i] <- exp(LL[i]) / C

ones[i] ~ dbern(p[i])

}


r ~ dgamma(a_alpha,b_alpha);
b ~ dgamma(a_alpha,b_alpha);

C <- 1000000
sigma <- 1/(b*pow(r,0.5))
Q <- pow(r,-0.5)
mu <- -beta[1] + (log(r)/b)


}"


dgengamma.orig(data.stan$t[1], 
               shape =mod$BUGSoutput$sims.matrix[1,"b"],
               scale =1/ mod$BUGSoutput$sims.matrix[1,"lambda[1]"],
               k= mod$BUGSoutput$sims.matrix[1,"r"], log = T)

#r = k
#b = shape 
#a = scale = 1/lambda

mod$BUGSoutput$sims.matrix[1,]

data.stan <- make_data_stan(formula=Surv(time,event)~1,data=digitized_IPD,
                            distr="gga",exArgs)

data.stan$t[18]

data.stan$d[18]
data.stan$mu_beta

data.stan$sigma_beta

1/(data.stan$sigma_beta[1])^2


b <- 1.2921
lambda <- exp(-4.451)
t <- data.stan$t[1]
r <- 0.519

log((b*lambda^(b*r)*t^(b*r-1)*exp(-(lambda*t)^b))/gamma(r))

#View(mod$BUGSoutput$sims.matrix)

mod <- R2jags:::jags( model.file = textConnection(GenGamma.jags),
                      data = data.stan,
                      parameters.to.save = c("beta", "b", "r", "lambda[1]"),
                      n.iter = 10000, n.thin = 1)
flexsurvreg(Surv(time, event)~1, data = digitized_IPD, dist = "gengamma")
flexsurvreg(Surv(time, event)~1, data = digitized_IPD, dist = "gengamma.orig")



data.stan <- make_data_stan(formula=Surv(time,event)~1,data=digitized_IPD,
                            distr="gom",exArgs)


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

alpha ~ dnorm(0,0.001)

for(i in 1:H){
 prec_beta[i] <- 1/(sigma_beta[i])^2
 beta[i] ~ dnorm(mu_beta[i],prec_beta[i])
}

linpred <- X %*%beta
for( i in 1:n){
mu[i] <- exp(linpred[i])
}


rate = exp(beta[1])
C <- 10000
s <- 0.0001

}"
#Need inits
inits <- function(){
  
  beta = c(log(1/mean(data.stan$t)),rep(0,data.stan$H -1))
  list(alpha = 0.001, beta = rep(log(1/mean(data.stan$t)),data.stan$H)) 
}
mod <- R2jags:::jags( model.file = textConnection(Gompertz.jags),
                      data = data.stan,
                      parameters.to.save = c("alpha", "beta"),
                      inits = inits,
                      n.iter = 10000, n.thin = 1)

Gompertz.jags<- "
data{
for(i in 1:n){
zero[i] <- 0}
}
model{
  C <- 10000
  for(i in 1:n){
    logHaz[i] <- (log(b)+ a*t[i])*d[i]
    logSurv[i] <- (-b/a)*(exp(a*t[i])-1)
    LL[i] <- logHaz[i]+ logSurv[i]
    Like[i] <- exp(LL[i])
    invLik[i] <-pow(Like[i],-1)
    zero[i] ~ dpois(zero.mean[i])
    zero.mean[i] <- -logHaz[i]-logSurv[i] + C
  }
  a ~ dnorm(0,0.01)
  b ~ dunif(0,5)
}"


inits <- function(){
  list(a = 0.001, b = 1/mean(data.stan$t))
}

#Need to add inits
mod2 <- R2jags:::jags( model.file = textConnection(Gompertz.jags),
                      data = data.stan,
                      parameters.to.save = c("a", "b"),
                      inits = inits,
                      n.iter = 10000, n.thin = 1)

flexsurv:::parse.dist("gomp")

digitized_IPD

c(1, 1/mean(t), 1)
flexsurv:::parse.dist("gengamma.orig")