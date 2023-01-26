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

alpha ~ dnorm(mu_beta[1],1/(sigma_beta[1])^2)

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
 
 
    St_expert[i,1] = exp(-mu[id_St]/alpha * (exp(alpha*time_expert[i]) - 1))*St_indic
    mean_surv_trt[i] <- (1/alpha)*exp(mu[id_trt]/alpha)*(exp(loggam(s))*(1-pgamma(mu[id_trt]/max(alpha,0.000000001),s,1)))
    mean_surv_comp[i] <- (1/alpha)*exp(mu[id_comp]/alpha)*(exp(loggam(s))*(1-pgamma(mu[id_comp]/max(alpha,0.0000001),s,1)))
    St_expert[i,2] <- (mean_surv_trt[i] - mean_surv_comp[i])*(ifelse(St_indic == 1, 0,1))
    
    St_expert_final[i] <- sum(St_expert[i,])
    
   # for(j in 1:num_param[i]){
    # expert_dens_mix[i,j] <- dnorm(St_expert_final[i], param_expert[i,j,1],pow(param_expert[i,j,2],-2))*param_expert[i,j,3]
   #}
   
    for(j in 1:7){
    expert_dens[j,1,i] <-  0  #sum(expert_dens_mix[i,1:num_param[i]])
    expert_dens[j,2,i] <-  dnorm(St_expert_final[i], param_expert[j,1,i],pow(param_expert[j,2,i],-2))
    expert_dens[j,3,i] <-  dt(St_expert_final[i],    param_expert[j,1,i],pow(param_expert[j,2,i],-2),3) 
    expert_dens[j,4,i] <-  dgamma(St_expert_final[i], max(param_expert[j,1,i],0.001),param_expert[j,2,i])
    expert_dens[j,5,i] <-  dlnorm(St_expert_final[i], param_expert[j,1,i],param_expert[j,2,i])
    expert_dens[j,6,i] <-  dbeta(St_expert_final[i], max(param_expert[j,1,i], 0.01),param_expert[j,2,i])
    phi_temp[j,i] <- expert_dens[j,data_dist_ind[j,1,i],i]/7
    }
 
    phi_con[i] <- -log(sum(phi_temp[,i])) + C

 } 
 
rate = exp(beta[1])
C <- 10000
s <- 0.0001

}"

dfs_expert <- readRDS(file = "dfs_expert.rds")
m.all <- readRDS("m_all.rds")

dfs_model <- dfs_expert
names(dfs_model)
for(i in 1:2){
  
  for(j in 1:7){
    if(  dfs_model[[i]]$dist[j] == "t"){
      dfs_model[[i]]$dist[j] = 3
    }
    if(  dfs_model[[i]]$dist[j] == "beta"){
      dfs_model[[i]]$dist[j] = 6
    }
  }
  
  dfs_model[[i]]$dist <- as.numeric(dfs_model[[i]]$dist )
  
  dfs_model[[i]] <- dfs_model[[i]][,-2]
}



model.data$param_expert <- NULL

array_final <- array( dim = c(dim( dfs_model[[1]]),2))
for(i in 1:2){
  dfs_model[[i]][ is.na(dfs_model[[i]])] <- 0
  array_final[,,i] <- as.matrix(dfs_model[[i]])
  
}

dt.scaled <- function (x, df, mean = 0, sd = 1, ncp, log = FALSE){
  if (!log) 
    stats::dt((x - mean)/sd, df, ncp = ncp, log = FALSE)/sd
  else stats::dt((x - mean)/sd, df, ncp = ncp, log = TRUE) - 
    log(sd)
}
expert_log_dens <- function(St, df){
  
  like <- c() 
  for(i in 1:nrow(df)){
    
    if(df$dist[i] == "t"){
      like <- c(like,dt.scaled(St, df$param3[i], df$param1[i], df$param2[i], log = F))
    }
    
    if(df$dist[i] == "t"){
      like <- c(like,dbeta(St, df$param1[i], df$param2[i],  log = F))
    }
    
    
  }
  return(log(sum(like)/nrow(df)))
  
  
}

digitized_IPD <- data.frame(digitized_IPD)
t <- digitized_IPD$time
event <- digitized_IPD$event
event2 <- ifelse(digitized_IPD$event ==1, 0,1)
times_expert = c(4,5)*12

names(m.all$models)
model.data <- m.all$misc$data.stan[[7]]
model.data$param_expert <-array_final[,-1,]
model.data$data_dist_ind <- array_final[,1,,drop=FALSE]
data.stan <- model.data

d3 <- "gom"
chains <- 3
iter <- 10000
thin <- 1
if (d3 %in% c("gam", "gga", "gom")){
  data.jags <- data.stan
  if(d3 %in% c( "gom")){
    parameters.to.save_jags = c("alpha","beta", "rate")
    
    #Inits as per flexsurvreg (reparameterized)
    modelinits <- function(){
      beta = c(log(1/mean(data.jags$t)*runif(1,0.8,1.5)),rep(0,data.jags$H -1))
      list(alpha = 0.001, beta = beta) 
    }
    
  }else if(d3 == "gga"){ #(d3 == "gga")
    parameters.to.save_jags = c("Q","sigma", "beta", "r", "b","mu")
    tinits1 <-data.jags$t + max(data.jags$t)
    is.na(tinits1)<-data.jags$d ==1
    data.jags$is.censored <- ifelse(data.jags$d==0, 1, 0)
    data.jags$t_jags <- ifelse(data.jags$is.censored ==1, NA, data.jags$t) 
    data.jags$t_cen <- data.jags$t+data.jags$d
    modelinits <- function(){list(t_jags = tinits1)}
    #Stop JAGS Warning messages
    data.jags <- data.jags[names(data.jags) %!in% c("t", "d", "a0")]
    
    
  }else{ #"gam",
    parameters.to.save_jags = c("alpha","beta", "rate")
    modelinits <- NULL
  }
  
  
}
  data.jags <- data.jags[names(data.jags) %!in% "max_param"]
  
  d <- "Gompertz"
model.gomp <-R2jags::jags(model.file = textConnection(get(paste0(d,".jags"))),
                     data=data.jags,
                     n.chains=chains,
                     inits=modelinits,
                     parameters.to.save = parameters.to.save_jags,
                     n.iter = iter*5,
                     n.thin = thin,
                     n.burnin = iter,
                     jags.module = c("glm","dic"))


df_param <- model.gomp$BUGSoutput$sims.matrix[,c("alpha", "rate")]
LL <- apply(df_param,1, 
            function(x){dgompertz(t,shape = x[1], rate = x[2], log = T)*event+
                pgompertz(t,shape = x[1], rate = x[2],log = T,lower.tail = F)*event2    
            })

LL_hat <- dgompertz(t,shape = colMeans(df_param)[1], rate = colMeans(df_param)[2], log = T)*event+
  pgompertz(t,shape = colMeans(df_param)[1], rate = colMeans(df_param)[2],log = T,lower.tail = F)*event2    

#waic.gomp <- loo::waic(t(LL))

LL.expert <- apply(df_param,1, function(x){
  expert_log_dens(pgompertz(times_expert[1],shape = x[1],rate = x[2],lower.tail = F),dfs_expert[[1]])+
    expert_log_dens(pgompertz(times_expert[2],shape = x[1],rate = x[2],lower.tail = F),dfs_expert[[2]])})


LL.expert_hat <- 
  expert_log_dens(pgompertz(times_expert[1],shape = colMeans(df_param)[1],rate = colMeans(df_param)[2],
                         lower.tail = F),dfs_expert[[1]])+
  expert_log_dens(pgompertz(times_expert[2],shape = colMeans(df_param)[1],rate = colMeans(df_param)[2],
                         lower.tail = F),dfs_expert[[2]])

pd_1 <- 2*(sum(LL_hat)+LL.expert_hat -mean(colSums(LL)+LL.expert)) 
DIC.gomp <- -2*(sum(LL_hat)+LL.expert_hat) + 2*pd_1

#Gamma JAGS function
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

    St_expert[i,1] =  (1-pgamma(time_expert[i],alpha,lambda[id_St]))*St_indic
    mean_surv_trt[i] <- alpha/lambda[id_trt]
    mean_surv_comp[i] <- alpha/lambda[id_comp]
    St_expert[i,2] <- (mean_surv_trt[i] - mean_surv_comp[i])*(ifelse(St_indic == 1, 0,1))
    
    St_expert_final[i] <- sum(St_expert[i,])
    
  # for(j in 1:num_param[i]){
    # expert_dens_mix[i,j] <- dnorm(St_expert_final[i], param_expert[i,j,1],pow(param_expert[i,j,2],-2))*param_expert[i,j,3]
   #}
   
    for(j in 1:7){
    expert_dens[j,1,i] <-  0  #sum(expert_dens_mix[i,1:num_param[i]])
    expert_dens[j,2,i] <-  dnorm(St_expert_final[i], param_expert[j,1,i],pow(param_expert[j,2,i],-2))
    expert_dens[j,3,i] <-  dt(St_expert_final[i],    param_expert[j,1,i],pow(param_expert[j,2,i],-2),3) 
    expert_dens[j,4,i] <-  dgamma(St_expert_final[i], max(param_expert[j,1,i],0.001),param_expert[j,2,i])
    expert_dens[j,5,i] <-  dlnorm(St_expert_final[i], param_expert[j,1,i],param_expert[j,2,i])
    expert_dens[j,6,i] <-  dbeta(St_expert_final[i], max(param_expert[j,1,i], 0.01),param_expert[j,2,i])
    phi_temp[j,i] <- expert_dens[j,data_dist_ind[j,1,i],i]/7
    }
 
    phi_con[i] <- -log(sum(phi_temp[,i])) + C

 } 


alpha ~ dgamma(a_alpha,b_alpha);

rate <- exp(beta[1]);

C <- 10000

}"


names(m.all$models)
model.data <- m.all$misc$data.stan[[6]]
model.data$param_expert <-array_final[,-1,]
model.data$data_dist_ind <- array_final[,1,,drop=FALSE]
data.stan <- model.data

d3 <- "gam"
chains <- 3
iter <- 1000
thin <- 1
if (d3 %in% c("gam", "gga", "gom")){
  data.jags <- data.stan
  if(d3 %in% c( "gom")){
    parameters.to.save_jags = c("alpha","beta", "rate")
    
    #Inits as per flexsurvreg (reparameterized)
    modelinits <- function(){
      beta = c(log(1/mean(data.jags$t)*runif(1,0.8,1.5)),rep(0,data.jags$H -1))
      list(alpha = 0.001, beta = beta) 
    }
    
  }else if(d3 == "gga"){ #(d3 == "gga")
    parameters.to.save_jags = c("Q","sigma", "beta", "r", "b","mu")
    tinits1 <-data.jags$t + max(data.jags$t)
    is.na(tinits1)<-data.jags$d ==1
    data.jags$is.censored <- ifelse(data.jags$d==0, 1, 0)
    data.jags$t_jags <- ifelse(data.jags$is.censored ==1, NA, data.jags$t) 
    data.jags$t_cen <- data.jags$t+data.jags$d
    modelinits <- function(){list(t_jags = tinits1)}
    #Stop JAGS Warning messages
    data.jags <- data.jags[names(data.jags) %!in% c("t", "d", "a0")]
    
    
  }else{ #"gam",
    parameters.to.save_jags = c("alpha","beta", "rate")
    modelinits <- NULL
  }
  
  
}
data.jags <- data.jags[names(data.jags) %!in% "max_param"]

d <- "Gamma"
model.gamma <-R2jags::jags(model.file = textConnection(get(paste0(d,".jags"))),
                     data=data.jags,
                     n.chains=chains,
                     inits=modelinits,
                     parameters.to.save = parameters.to.save_jags,
                     n.iter = iter*5,
                     n.thin = thin,
                     n.burnin = iter,
                     jags.module = c("glm","dic"))




df_param <- model.gamma$BUGSoutput$sims.matrix[,c("alpha", "rate")]
LL <- apply(df_param,1, 
            function(x){dgamma(t,shape = x[1], rate = x[2], log = T)*event+
                pgamma(t,shape = x[1], rate = x[2],log = T,lower.tail = F)*event2    
            })
#waic.gam <- loo::waic(t(LL))

LL_hat <- dgamma(t,shape = colMeans(df_param)[1], rate = colMeans(df_param)[2], log = T)*event+
  pgamma(t,shape = colMeans(df_param)[1], rate = colMeans(df_param)[2],log = T,lower.tail = F)*event2    

LL.expert <- apply(df_param,1, function(x){
  expert_log_dens(pgamma(times_expert[1],shape = x[1],rate = x[2],lower.tail = F),dfs_expert[[1]])+
    expert_log_dens(pgamma(times_expert[2],shape = x[1],rate = x[2],lower.tail = F),dfs_expert[[2]])})


LL.expert_hat <- 
  expert_log_dens(pgamma(times_expert[1],shape = colMeans(df_param)[1],rate = colMeans(df_param)[2],
                      lower.tail = F),dfs_expert[[1]])+
  expert_log_dens(pgamma(times_expert[2],shape = colMeans(df_param)[1],rate = colMeans(df_param)[2],
                      lower.tail = F),dfs_expert[[2]])

pd_1 <- 2*(sum(LL_hat)+LL.expert_hat -mean(colSums(LL)+LL.expert)) 
DIC.gam <- -2*(sum(LL_hat)+LL.expert_hat) + 2*pd_1


#Generalized Gamma function
#Generalized Gamma function
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

    St_expert[i,1] =  (1-pgen.gamma(time_expert[i],r,lambda[id_St],b))*St_indic
    mean_surv_trt[i] <- exp(loggam((b*r +1)/b) - loggam(r))/lambda[id_comp]
    mean_surv_comp[i] <- exp(loggam((b*r +1)/b) - loggam(r))/lambda[id_trt]
    St_expert[i,2] <- (mean_surv_trt[i] - mean_surv_comp[i])*(ifelse(St_indic == 1, 0,1))
    
    St_expert_final[i] <- sum(St_expert[i,])
    
    
  # for(j in 1:num_param[i]){
    # expert_dens_mix[i,j] <- dnorm(St_expert_final[i], param_expert[i,j,1],pow(param_expert[i,j,2],-2))*param_expert[i,j,3]
   #}
   
    for(j in 1:7){
    expert_dens[j,1,i] <-  0  #sum(expert_dens_mix[i,1:num_param[i]])
    expert_dens[j,2,i] <-  dnorm(St_expert_final[i], param_expert[j,1,i],pow(param_expert[j,2,i],-2))
    expert_dens[j,3,i] <-  dt(St_expert_final[i],    param_expert[j,1,i],pow(param_expert[j,2,i],-2),3) 
    expert_dens[j,4,i] <-  dgamma(St_expert_final[i], max(param_expert[j,1,i],0.001),param_expert[j,2,i])
    expert_dens[j,5,i] <-  dlnorm(St_expert_final[i], param_expert[j,1,i],param_expert[j,2,i])
    expert_dens[j,6,i] <-  dbeta(St_expert_final[i], max(param_expert[j,1,i], 0.01),param_expert[j,2,i])
    phi_temp[j,i] <- expert_dens[j,data_dist_ind[j,1,i],i]/7
    }
 
    phi_con[i] <- -log(sum(phi_temp[,i])) + C
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
names(m.all$models)
model.data <- m.all$misc$data.stan[[8]]
model.data$param_expert <-array_final[,-1,]
model.data$data_dist_ind <- array_final[,1,,drop=FALSE]
data.stan <- model.data

d3 <- "gga"
chains <- 3
iter <- 1000
thin <- 1
if (d3 %in% c("gam", "gga", "gom")){
  data.jags <- data.stan
  if(d3 %in% c( "gom")){
    parameters.to.save_jags = c("alpha","beta", "rate")
    
    #Inits as per flexsurvreg (reparameterized)
    modelinits <- function(){
      beta = c(log(1/mean(data.jags$t)*runif(1,0.8,1.5)),rep(0,data.jags$H -1))
      list(alpha = 0.001, beta = beta) 
    }
    
  }else if(d3 == "gga"){ #(d3 == "gga")
    parameters.to.save_jags = c("Q","sigma", "beta", "r", "b","mu")
    tinits1 <-data.jags$t + max(data.jags$t)
    is.na(tinits1)<-data.jags$d ==1
    data.jags$is.censored <- ifelse(data.jags$d==0, 1, 0)
    data.jags$t_jags <- ifelse(data.jags$is.censored ==1, NA, data.jags$t) 
    data.jags$t_cen <- data.jags$t+data.jags$d
    modelinits <- function(){list(t_jags = tinits1)}
    #Stop JAGS Warning messages
    data.jags <- data.jags[names(data.jags) %!in% c("t", "d", "a0")]
    
    
  }else{ #"gam",
    parameters.to.save_jags = c("alpha","beta", "rate")
    modelinits <- NULL
  }
  
  
}
#data.jags <- data.jags[names(data.jags) %!in% "max_param"]

d <- "GenGamma"
model.gengamma <-R2jags::jags(model.file = textConnection(get(paste0(d,".jags"))),
                     data=data.jags,
                     n.chains=chains,
                     inits=modelinits,
                     parameters.to.save = parameters.to.save_jags,
                     n.iter = iter*5,
                     n.thin = thin,
                     n.burnin = iter,
                     jags.module = c("glm","dic"))



df_param <- model.gengamma$BUGSoutput$sims.matrix[,c("mu", "sigma", "Q")]

LL <- apply(df_param,1, 
            function(x){dgengamma(t,mu = x[1], sigma = x[2],Q = x[3], log = T)*event+
                pgengamma(t,mu = x[1], sigma = x[2],Q = x[3],log = T,lower.tail = F)*event2    
            })
#waic.gengamma <- loo::waic(t(LL))

LL_hat <- dgengamma(t,mu = colMeans(df_param)[1], sigma = colMeans(df_param)[2],Q = colMeans(df_param)[3], log = T)*event+
  pgengamma(t,mu = colMeans(df_param)[1], sigma = colMeans(df_param)[2],Q = colMeans(df_param)[3],log = T,lower.tail = F)*event2    



LL.expert <- apply(df_param,1, function(x){
  expert_log_dens(pgengamma(times_expert[1],mu = x[1],sigma = x[2],Q = x[3],lower.tail = F),
               dfs_expert[[1]])+
    expert_log_dens(pgengamma(times_expert[2],mu = x[1],sigma = x[2],Q = x[3],lower.tail = F),
                 dfs_expert[[2]])})


LL.expert_hat <- 
  expert_log_dens(pgengamma(times_expert[1],
                         mu = colMeans(df_param)[1],
                         sigma = colMeans(df_param)[2],
                         Q = colMeans(df_param)[3],
                         lower.tail = F),dfs_expert[[1]])+
  expert_log_dens(pgengamma(times_expert[2],
                         mu = colMeans(df_param)[1],
                         sigma = colMeans(df_param)[2],
                         Q = colMeans(df_param)[3],
                         lower.tail = F),dfs_expert[[2]])

pd_1 <- 2*(sum(LL_hat)+LL.expert_hat -mean(colSums(LL)+LL.expert)) 
DIC.gengamma <- -2*(sum(LL_hat)+LL.expert_hat) + 2*pd_1

