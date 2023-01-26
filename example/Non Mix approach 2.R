Weibull.example <- "// Weibull survival model (with PH parameterisation)

functions {
  // Defines the log hazard
  // Defines the log hazard
  vector log_haz (vector t, real shape, vector scale) {
    vector[num_elements(t)] log_haz;
    log_haz = log(shape)+(shape-1)*log(t ./ scale)-log(scale);
    return log_haz;
  }
  
  // Defines the log survival
  vector log_S (vector t, real shape, vector scale) {
    vector[num_elements(t)] log_S;
    for (i in 1:num_elements(t)) {
      log_S[i] = -pow((t[i]/scale[i]),shape);
    }
    return log_S;
  }
  
  
    // Defines the log survival indvidual
  real log_Sind (real t, real shape, real scale) {
	real log_Sind;
      log_Sind = -pow((t/scale),shape);
    return log_Sind;
  }
  
  
      // Defines difference in expected survival
  real Surv_diff ( real shape, real scale_trt, real scale_comp) {
	real Surv_diff;
      Surv_diff = (scale_trt-scale_comp)*tgamma(1 +1/shape);
    return Surv_diff;
  }
  
  // Defines the sampling distribution
  real surv_weibullAF_lpdf (vector t, vector d, real shape, vector scale, vector a0) {
    vector[num_elements(t)] log_lik;
    real prob;
    log_lik = d .* log_haz(t,shape,scale) + log_S(t,shape,scale);
    prob = dot_product(log_lik, a0);
    return prob;
  }
  
  real log_density_dist(real[ , ] params, 
                        real x,int num_expert, int pool_type){
    
    // Evaluates the log density for a range of distributions
    
    real dens[num_expert];
    
    for(i in 1:num_expert){
    if(params[i,1] == 1){
      if(pool_type == 1){
        dens[i] = exp(normal_lpdf(x|params[i,3], params[i,4]))*params[i,2]; /// Only require the log density is correct to a constant of proportionality
      }else{
        dens[i] = exp(normal_lpdf(x|params[i,3], params[i,4]))^params[i,2]; /// Only require the log density is correct to a constant of proportionality
      }
      
    }else if(params[i,1] == 2){
       if(pool_type == 1){
          dens[i] = exp(student_t_lpdf(x|params[i,5],params[i,3], params[i,4]))*params[i,2];
        }else{
          dens[i] = exp(student_t_lpdf(x|params[i,5],params[i,3], params[i,4]))^params[i,2];
       }
      
    }else if(params[i,1] == 3){
        if(pool_type == 1){
            dens[i] = exp(gamma_lpdf(x|params[i,3], params[i,4]))*params[i,2];
        }else{
            dens[i] = exp(gamma_lpdf(x|params[i,3], params[i,4]))^params[i,2];
        }  
            
    }else if(params[i,1] == 4){
    
        if(pool_type == 1){
            dens[i] = exp(lognormal_lpdf(x|params[i,3], params[i,4]))*params[i,2];
        }else{
             dens[i] = exp(lognormal_lpdf(x|params[i,3], params[i,4]))^params[i,2];
        }
        
    }else if(params[i,1] == 5){
     if(pool_type == 1){
            dens[i] = exp(beta_lpdf(x|params[i,3], params[i,4]))*params[i,2];
        }else{
            dens[i] = exp(beta_lpdf(x|params[i,3], params[i,4]))^params[i,2];
        }
      
      
      } 

    }
    
      
    return(log(sum(dens)));
    
  }
  
}

data {
  int n;                  // number of observations
  vector[n] t;            // observed times
  vector[n] d;            // censoring indicator (1=observed, 0=censored)
  int H;                  // number of covariates
  matrix[n,H] X;          // matrix of covariates (with n rows and H columns)
  vector[H] mu_beta;	  // mean of the covariates coefficients
  vector<lower=0> [H] sigma_beta;   // sd of the covariates coefficients
  real<lower=0> a_alpha;
  real<lower=0> b_alpha;
  vector[n] a0;  
  int n_time_expert;
  int<lower = 0, upper = 1> St_indic; // 1 Expert opinion on survival @ timepoint ; 0 Expert opinion on survival difference
  
  int id_St;
  int id_trt;
  int id_comp;
  
  int n_experts[n_time_expert];
  int pool_type;

  real param_expert[max(n_experts),5,n_time_expert];
  vector[St_indic ? n_time_expert : 0] time_expert;


  
  
}

parameters {
  vector[H] beta;         // Coefficients in the linear predictor (including intercept)
  real<lower=0> alpha;    // shape parameter
}

transformed parameters {
  vector[n] linpred;
  vector[n] mu;
  vector[n_time_expert] St_expert;

  linpred = X*beta;
  for (i in 1:n) {
    mu[i] = exp(linpred[i]);
  }
  
   for (i in 1:n_time_expert){

    if(St_indic == 1){
 
    St_expert[i] = exp(log_Sind(time_expert[i],alpha,mu[id_St]));
    }else 
	  St_expert[i] = Surv_diff(alpha,mu[id_trt],mu[id_comp]);
  
  }
  
 
}

model {
  alpha ~ gamma(a_alpha,b_alpha);
  beta ~ normal(mu_beta,sigma_beta);
  t ~ surv_weibullAF(d,alpha,mu, a0);
  
  for (i in 1:n_time_expert){
      
      target += log_density_dist(param_expert[,,i],
                                 St_expert[i],
                                 n_experts[i],
                                 pool_type);
 

  } 
  
}

generated quantities {
  real scale;                // scale parameter
  scale = exp(beta[1]);
}"



dfs_expert <- readRDS(file = "dfs_expert.rds")
m.all <- readRDS("m_all.rds")

library("rstan")

model.data <- m.all$misc$data.stan[[4]]
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

model.data$param_expert <-array_final
library("rstan")

fit.weibull <- stan(model_code = Weibull.example, model_name = "WeibullAF", 
                    data =  stan.data, iter = 2000, chains = 3)

model <- fit.weibull
distr3 <- "wei"
compute_ICs_stan(model,distr3, stan.data)

fit.weibull <- stan(model_code = Weibull.example, model_name = "WeibullAF", 
            data =  model.data, iter = 2000, chains = 3)

df_beta <- data.frame(extract(fit, c("beta[1]")))
colMeans(exp(df_beta))


colMeans(data.frame(extract(fit, c("St_expert"))))
m.all$models$`Weibull (AFT)`




#View(fit)


df_St <- data.frame(extract(fit, c("St_expert")))
colMeans(df_St)

#df_beta <- data.frame(extract(fit, c("alpha", "scale")))
df_beta <- data.frame(extract(fit, c("beta")))
colMeans(df_beta)









eval_dens_pool(x.eval,pool.df,pool_type = "log pool")


df_param <- data.frame(extract(fit.weibull, c("scale", "alpha")))

t <- digitized_IPD$time
event <- digitized_IPD$event
event2 <- ifelse(digitized_IPD$event ==1, 0,1)
times_expert = c(4,5)*12


LL <- apply(df_param,1, 
            function(x){dweibull(t,scale = x[1], shape = x[2], log = T)*event+
                pweibull(t,scale = x[1], shape = x[2],log = T,lower.tail = F)*event2
            })


expert_log_dens(0.2, df)


LL.expert <- apply(df_param,1, function(x){
  expert_log_dens(pweibull(times_expert[1],scale = x[1],shape = x[2],lower.tail = F),
                  dfs_expert[[1]])+
    expert_log_dens(pweibull(times_expert[2],scale = x[1],shape = x[2],lower.tail = F),
                    dfs_expert[[2]])})


LL.expert_hat <- 
  expert_log_dens(pweibull(times_expert[1],scale = colMeans(df_param)[1],shape = colMeans(df_param)[2],
                           lower.tail = F),dfs_expert[[1]])+
  expert_log_dens(pweibull(times_expert[2],scale = colMeans(df_param)[1],shape = colMeans(df_param)[2],
                           lower.tail = F),dfs_expert[[2]])



df_hat <- colMeans(df_param)
#df_hat <- apply(df_param, 2, median)

LL_hat <- dweibull(t,scale = df_hat[1], shape = df_hat[2], log = T)*event+
  pweibull(t,scale = df_hat[1], shape = df_hat[2],log = T,lower.tail = F)*event2

pd_1 <- 2*(sum(LL_hat) - mean(colSums(LL))) 

DIC.weib <- -2*sum(LL_hat) + 2*pd_1



lognorm.example <- "// log-Normal survival model

functions {
  // Defines the log survival
  vector log_S (vector t, vector mean, real sd) {
    vector[num_elements(t)] log_S;
    for (i in 1:num_elements(t)) {
      log_S[i] = log(1-Phi((log(t[i])-mean[i])/sd));
    }
    return log_S;
  }
  
  // Defines the log survival indvidual
  real log_Sind (real t, real mean, real sd) {
    real log_Sind;
    log_Sind = log(1-Phi((log(t)-mean)/sd));
    return log_Sind;
  }
  
  // Defines difference in expected survival
  real Surv_diff (real mean_trt, real mean_comp, real sd ) {
    real Surv_diff;
    real Surv_trt;
    real Surv_comp;
    
    Surv_trt = exp(mean_trt + 0.5 * pow(sd,2));
    Surv_comp = exp(mean_comp + 0.5 * pow(sd,2));
    
    Surv_diff = Surv_trt - Surv_comp;
    return Surv_diff;
  }
  
  
  // Defines the log hazard
  vector log_h (vector t, vector mean, real sd) {
    vector[num_elements(t)] log_h;
    vector[num_elements(t)] ls;
    ls = log_S(t,mean,sd);
    for (i in 1:num_elements(t)) {
      log_h[i] = lognormal_lpdf(t[i]|mean[i],sd) - ls[i];
    }
    return log_h;
  }
  
  // Defines the sampling distribution
  real surv_lognormal_lpdf (vector t, vector d, vector mean, real sd, vector a0) {
    vector[num_elements(t)] log_lik;
    real prob;
    log_lik = d .* log_h(t,mean,sd) + log_S(t,mean,sd);
    prob = dot_product(log_lik, a0);
    return prob;
  }
  
  real log_density_dist(real[ , ] params, 
                        real x,int num_expert, int pool_type){
    
    // Evaluates the log density for a range of distributions
    
    real dens[num_expert];
    
    for(i in 1:num_expert){
    if(params[i,1] == 1){
      if(pool_type == 1){
        dens[i] = exp(normal_lpdf(x|params[i,3], params[i,4]))*params[i,2]; /// Only require the log density is correct to a constant of proportionality
      }else{
        dens[i] = exp(normal_lpdf(x|params[i,3], params[i,4]))^params[i,2]; /// Only require the log density is correct to a constant of proportionality
      }
      
    }else if(params[i,1] == 2){
       if(pool_type == 1){
          dens[i] = exp(student_t_lpdf(x|params[i,5],params[i,3], params[i,4]))*params[i,2];
        }else{
          dens[i] = exp(student_t_lpdf(x|params[i,5],params[i,3], params[i,4]))^params[i,2];
       }
      
    }else if(params[i,1] == 3){
        if(pool_type == 1){
            dens[i] = exp(gamma_lpdf(x|params[i,3], params[i,4]))*params[i,2];
        }else{
            dens[i] = exp(gamma_lpdf(x|params[i,3], params[i,4]))^params[i,2];
        }  
            
    }else if(params[i,1] == 4){
    
        if(pool_type == 1){
            dens[i] = exp(lognormal_lpdf(x|params[i,3], params[i,4]))*params[i,2];
        }else{
             dens[i] = exp(lognormal_lpdf(x|params[i,3], params[i,4]))^params[i,2];
        }
        
    }else if(params[i,1] == 5){
     if(pool_type == 1){
            dens[i] = exp(beta_lpdf(x|params[i,3], params[i,4]))*params[i,2];
        }else{
            dens[i] = exp(beta_lpdf(x|params[i,3], params[i,4]))^params[i,2];
        }
      
      
      } 

    }
    
      
    return(log(sum(dens)));
    
  }

  
  
}

data {
  int n;                  // number of observations
  vector[n] t;            // observed times
  vector[n] d;            // censoring indicator (1=observed, 0=censored)
  int H;                  // number of covariates
  matrix[n,H] X;          // matrix of covariates (with n rows and H columns)
  vector[H] mu_beta;	    // mean of the covariates coefficients
  vector<lower=0> [H] sigma_beta;   // sd of the covariates coefficients
  real a_alpha;			      // lower bound for the sd of the data			  
  real b_alpha;			      // upper bound for the sd of the data
  
  vector[n] a0; //Power prior for the observations
  int<lower = 0, upper = 1> St_indic; // 1 Expert opinion on survival @ timepoint ; 0 Expert opinion on survival difference
  int n_time_expert;

  int id_St;
  int id_trt;
  int id_comp;
  
  int n_experts[n_time_expert];
  int pool_type;

  real param_expert[max(n_experts),5,n_time_expert];
  vector[St_indic ? n_time_expert : 0] time_expert;
  
}

parameters {
  vector[H] beta;         // Coefficients in the linear predictor (including intercept)
  real<lower=0> alpha;    // log-sd parameter
}

transformed parameters {
  vector[n] linpred;
  vector[n] mu;
  vector[n_time_expert] St_expert;
  
  linpred = X*beta;
  for (i in 1:n) {
    mu[i] = linpred[i];
  }
  
  for (i in 1:n_time_expert){
    if(St_indic == 1){
      
      St_expert[i] = exp(log_Sind(time_expert[i],mu[id_St[1]],alpha));
      
    }else{    
      St_expert[i] = Surv_diff(mu[id_trt],mu[id_comp], alpha);  
    }
  }
  
}

model {
  alpha ~ uniform(a_alpha,b_alpha);
  beta ~ normal(mu_beta,sigma_beta);
  t ~ surv_lognormal(d,X*beta,alpha, a0);
  
  for (i in 1:n_time_expert){
     

     target += log_density_dist(param_expert[,,i],
                                 St_expert[i],
                                 n_experts[i],
                                 pool_type);
 
  }

}

generated quantities {
  real meanlog;
  meanlog = beta[1];
}"



model.data <- m.all$misc$data.stan[[3]]
model.data$param_expert <-array_final



fit.lognorm <- stan(model_code = lognorm.example, model_name = "LogNormal", 
            data =  model.data, iter = 2000, chains = 3)
m.all$models$`log-Normal`

# pd_1 <- 2*(sum(LL_hat) -mean(colSums(LL))) 
df_param <- data.frame(extract(fit.lognorm, c("meanlog", "alpha")))
LL <- apply(df_param,1, 
            function(x){dlnorm(t,meanlog = x[1], sdlog = x[2], log = T)*event+
                plnorm(t,meanlog = x[1], sdlog = x[2],log = T,lower.tail = F)*event2    
            })
#waic.lnorm <- loo::waic(t(LL))


LL.expert <- apply(df_param,1, function(x){
  expert_log_dens(plnorm(times_expert[1],meanlog = x[1],sdlog = x[2],lower.tail = F),dfs_expert[[1]])+
    expert_log_dens(plnorm(times_expert[2],meanlog = x[1],sdlog = x[2],lower.tail = F),dfs_expert[[2]])})


LL.expert_hat <- 
  expert_log_dens(plnorm(times_expert[1],meanlog = colMeans(df_param)[1],sdlog = colMeans(df_param)[2],
                      lower.tail = F),dfs_expert[[1]])+
  expert_log_dens(plnorm(times_expert[2],meanlog = colMeans(df_param)[1],sdlog = colMeans(df_param)[2],
                      lower.tail = F),dfs_expert[[2]])

pd_1 <- 2*(sum(LL_hat)+LL.expert_hat -mean(colSums(LL)+LL.expert)) 
DIC.lnorm <- -2*(sum(LL_hat)+LL.expert_hat) + 2*pd_1

names(m.all$models)
#m.all$misc$data.stan[[3]]

exponential.example <- "// Exponential survival model

functions {
  // Defines the log hazard
  vector log_h (vector t, vector rate) {
    vector[num_elements(t)] log_h;
    log_h = log(rate);
    return log_h;
  }
  
  // Defines the log survival
  vector log_S (vector t, vector rate) {
    vector[num_elements(t)] log_S;
    log_S = -rate .* t;
    return log_S;
  }
  
  // Defines the log survival indvidual
  real log_Sind (real t, real rate) {
	real log_Sind;
      log_Sind = -rate .* t;
    return log_Sind;
  }  
  
  // Defines difference in expected survival
  real Surv_diff ( real rate_trt, real rate_comp) {
	real Surv_diff;
      Surv_diff = 1/rate_trt - 1/rate_comp;
    return Surv_diff;
  }
    
  // Defines the sampling distribution
  real surv_exponential_lpdf (vector t, vector d, vector rate, vector a0) {
    vector[num_elements(t)] log_lik;
    real prob;
    log_lik = d .* log_h(t,rate) + log_S(t,rate);
    prob = dot_product(log_lik, a0);
    return prob;
  } 
  

   real log_density_dist(real[ , ] params, 
                        real x,int num_expert, int pool_type){
    
    // Evaluates the log density for a range of distributions
    
    real dens[num_expert];
    
    for(i in 1:num_expert){
    if(params[i,1] == 1){
      if(pool_type == 1){
        dens[i] = exp(normal_lpdf(x|params[i,3], params[i,4]))*params[i,2]; /// Only require the log density is correct to a constant of proportionality
      }else{
        dens[i] = exp(normal_lpdf(x|params[i,3], params[i,4]))^params[i,2]; /// Only require the log density is correct to a constant of proportionality
      }
      
    }else if(params[i,1] == 2){
       if(pool_type == 1){
          dens[i] = exp(student_t_lpdf(x|params[i,5],params[i,3], params[i,4]))*params[i,2];
        }else{
          dens[i] = exp(student_t_lpdf(x|params[i,5],params[i,3], params[i,4]))^params[i,2];
       }
      
    }else if(params[i,1] == 3){
        if(pool_type == 1){
            dens[i] = exp(gamma_lpdf(x|params[i,3], params[i,4]))*params[i,2];
        }else{
            dens[i] = exp(gamma_lpdf(x|params[i,3], params[i,4]))^params[i,2];
        }  
            
    }else if(params[i,1] == 4){
    
        if(pool_type == 1){
            dens[i] = exp(lognormal_lpdf(x|params[i,3], params[i,4]))*params[i,2];
        }else{
             dens[i] = exp(lognormal_lpdf(x|params[i,3], params[i,4]))^params[i,2];
        }
        
    }else if(params[i,1] == 5){
     if(pool_type == 1){
            dens[i] = exp(beta_lpdf(x|params[i,3], params[i,4]))*params[i,2];
        }else{
            dens[i] = exp(beta_lpdf(x|params[i,3], params[i,4]))^params[i,2];
        }
      
      
      } 

    }
    
      
    return(log(sum(dens)));
    
  }

  
}

data {
  int n;                  // number of observations
  vector[n] t;            // observed times
  vector[n] d;            // censoring indicator (1=observed, 0=censored)
  int H;                  // number of covariates
  matrix[n,H] X;          // matrix of covariates (with n rows and H columns)
  vector[H] mu_beta;	    // mean of the covariates coefficients
  vector<lower=0> [H] sigma_beta;   // sd of the covariates coefficients
  
  vector[n] a0; //Power prior for the observations
  
  int n_time_expert;
  int<lower = 0, upper = 1> St_indic; // 1 Expert opinion on survival @ timepoint ; 0 Expert opinion on survival difference
  
  int id_St;
  int id_trt;
  int id_comp;
  
  int n_experts[n_time_expert];
  int pool_type;

  real param_expert[max(n_experts),5,n_time_expert];
  vector[St_indic ? n_time_expert : 0] time_expert;

  
}

parameters {
  vector[H] beta;         // Coefficients in the linear predictor (including intercept)
}

transformed parameters {
  vector[n] linpred;
  vector[n] mu;
  vector[n_time_expert] St_expert;
 
  
  linpred = X*beta;
  for (i in 1:n) {
    mu[i] = exp(linpred[i]);  // Rate parameter
  }
  
  for (i in 1:n_time_expert){
  if(St_indic == 1){
	
		St_expert[i] = exp(log_Sind(time_expert[i],mu[id_St]));

	 }else{    
	   St_expert[i] = Surv_diff(mu[id_trt],mu[id_comp]]);
  	}
  }
}

model {
  beta ~ normal(mu_beta,sigma_beta);
  t ~ surv_exponential(d,mu, a0);
  
  
  for (i in 1:n_time_expert){
     
     target += log_density_dist(param_expert[,,i],
                                 St_expert[i],
                                 n_experts[i],
                                 pool_type);
  }
  
  
}

generated quantities {
  real rate;                // rate parameter
  rate = exp(beta[1]);
}"




model.data <- m.all$misc$data.stan[[1]]
model.data$param_expert <-array_final



fit.Exponential <- stan(model_code = exponential.example, model_name = "Exponential", 
            data =  model.data, iter = 2000, chains = 3)


df_param <- data.frame(extract(fit.Exponential, "rate"))
LL <- apply(df_param,1, 
            function(x){dexp(t,rate = x[1],log = T)*event+
                pexp(t,rate = x[1],log = T,lower.tail = F)*event2})

LL.expert <- apply(df_param,1, function(x){
  expert_log_dens(pexp(times_expert[1],rate = x,lower.tail = F),dfs_expert[[1]])+
    expert_log_dens(pexp(times_expert[2],rate = x,lower.tail = F),dfs_expert[[2]])})


LL_hat <- dexp(t,rate = colMeans(df_param),log = T)*event+
  pexp(t,rate = colMeans(df_param),log = T,lower.tail = F)*event2

LL.expert_hat <- 
  expert_log_dens(pexp(times_expert[1],rate = colMeans(df_param),lower.tail = F),
               dfs_expert[[1]])+
  expert_log_dens(pexp(times_expert[2],rate = colMeans(df_param),lower.tail = F),
               dfs_expert[[2]])



pd_1 <- 2*(sum(LL_hat)+LL.expert_hat -mean(colSums(LL)+LL.expert)) 

DIC.exp <- -2*(sum(LL_hat)+LL.expert_hat) + 2*pd_1
#


loglogistic.example <- "// log-Logistic survival model

functions {
  // Defines the log hazard
  vector log_h (vector t, real shape, vector scale) {
    vector[num_elements(t)] log_h;
    for (i in 1:num_elements(t)) {
      log_h[i] = log(shape)-log(scale[i])+(shape-1)*(log(t[i])-log(scale[i]))-log(1+pow((t[i]/scale[i]),shape));
    }
    return log_h;
  }
  
  // Defines the log survival
  vector log_S (vector t, real shape, vector scale) {
    vector[num_elements(t)] log_S;
    for (i in 1:num_elements(t)) {
      log_S[i] = -log(1+pow((t[i]/scale[i]),shape));
    }
    return log_S;
  }
  
  // Defines the  survival indvidual
  real Sind (real t, real shape, real scale) {
    real Sind;
    Sind = exp(-log(1+pow((t/scale),shape)));
    return Sind;
  }
  
  // Defines difference in expected survival
  real Surv_diff ( real shape, real scale_trt, real scale_comp) {
    real Surv_diff;
    real b;
    b = pi()/shape;
    Surv_diff = (scale_trt-scale_comp)*b/sin(b);
    return Surv_diff;
  }
  
  
  
  // Defines the sampling distribution
  real surv_loglogistic_lpdf (vector t, vector d, real shape, vector scale, vector a0) {
    vector[num_elements(t)] log_lik;
    real prob;
    log_lik = d .* log_h(t,shape,scale) + log_S(t,shape,scale);
    prob = dot_product(log_lik,a0);
    return prob;
  }
  

   real log_density_dist(real[ , ] params, 
                        real x,int num_expert, int pool_type){
    
    // Evaluates the log density for a range of distributions
    
    real dens[num_expert];
    
    for(i in 1:num_expert){
    if(params[i,1] == 1){
      if(pool_type == 1){
        dens[i] = exp(normal_lpdf(x|params[i,3], params[i,4]))*params[i,2]; /// Only require the log density is correct to a constant of proportionality
      }else{
        dens[i] = exp(normal_lpdf(x|params[i,3], params[i,4]))^params[i,2]; /// Only require the log density is correct to a constant of proportionality
      }
      
    }else if(params[i,1] == 2){
       if(pool_type == 1){
          dens[i] = exp(student_t_lpdf(x|params[i,5],params[i,3], params[i,4]))*params[i,2];
        }else{
          dens[i] = exp(student_t_lpdf(x|params[i,5],params[i,3], params[i,4]))^params[i,2];
       }
      
    }else if(params[i,1] == 3){
        if(pool_type == 1){
            dens[i] = exp(gamma_lpdf(x|params[i,3], params[i,4]))*params[i,2];
        }else{
            dens[i] = exp(gamma_lpdf(x|params[i,3], params[i,4]))^params[i,2];
        }  
            
    }else if(params[i,1] == 4){
    
        if(pool_type == 1){
            dens[i] = exp(lognormal_lpdf(x|params[i,3], params[i,4]))*params[i,2];
        }else{
             dens[i] = exp(lognormal_lpdf(x|params[i,3], params[i,4]))^params[i,2];
        }
        
    }else if(params[i,1] == 5){
     if(pool_type == 1){
            dens[i] = exp(beta_lpdf(x|params[i,3], params[i,4]))*params[i,2];
        }else{
            dens[i] = exp(beta_lpdf(x|params[i,3], params[i,4]))^params[i,2];
        }
      
      
      } 

    }
    
      
    return(log(sum(dens)));
    
  }  
}

data {
  int n;                  // number of observations
  vector[n] t;            // observed times
  vector[n] d;            // censoring indicator (1=observed, 0=censored)
  int H;                  // number of covariates
  matrix[n,H] X;          // matrix of covariates (with n rows and H columns)
  vector[H] mu_beta;	  // mean of the covariates coefficients
  vector<lower=0> [H] sigma_beta;   // sd of the covariates coefficients
  real<lower=0> a_alpha;
  real<lower=0> b_alpha;
  vector[n] a0; //Power prior for the observations
  
  int<lower = 0, upper = 1> St_indic; // 1 Expert opinion on survival @ timepoint ; 0 Expert opinion on survival difference
  int n_time_expert;
  int id_St;
  int id_trt;
  int id_comp;
  
  int n_experts[n_time_expert];
  int pool_type;

  real param_expert[max(n_experts),5,n_time_expert];
  vector[St_indic ? n_time_expert : 0] time_expert;
  
  
}

parameters {
  vector[H] beta;         // Coefficients in the linear predictor (including intercept)
  real<lower=1> alpha;    // shape parameter - Constrainted to be greater than 1 so that the expected value exists
}

transformed parameters {
  vector[n] linpred;
  vector[n] mu;
  vector[n_time_expert] St_expert;
  
  linpred = X*beta;
  for (i in 1:n) {
    mu[i] = exp(linpred[i]);
  }
  
  for (i in 1:n_time_expert){
    if(St_indic == 1){
      
      St_expert[i] = Sind(time_expert[i],alpha,mu[id_St]);
      
    }else{    
      St_expert[i] = Surv_diff(alpha,mu[id_trt],mu[id_comp]);  
    }
  }
  
}

model {
  alpha ~ gamma(a_alpha,b_alpha);
  beta ~ normal(mu_beta,sigma_beta);
  t ~ surv_loglogistic(d,alpha,mu, a0);
  
  for (i in 1:n_time_expert){
     
     target += log_density_dist(param_expert[,,i],
                                 St_expert[i],
                                 n_experts[i],
                                 pool_type);
  }
  
}

generated quantities {
  real rate;                // rate parameter
  rate = exp(beta[1]);
}"


names(m.all$models)
model.data <- m.all$misc$data.stan[[2]]
model.data$param_expert <-array_final



fit.loglogistic <- stan(model_code = loglogistic.example, model_name = "Loglogistic Example", 
            data =  model.data, iter = 2000, chains = 3)



df_param <- data.frame(extract(fit.loglogistic, c("alpha", "rate")))
LL <- apply(df_param,1, 
            function(x){dllogis(t,shape = x[1], scale = x[2], log = T)*event+
                pllogis(t,shape = x[1], scale = x[2],log = T,lower.tail = F)*event2})
#waic.llogis <- loo::waic(t(LL))
LL_hat <- dllogis(t,shape = colMeans(df_param)[1], scale = colMeans(df_param)[2], log = T)*event+
  pllogis(t,shape = colMeans(df_param)[1], scale = colMeans(df_param)[2],log = T,lower.tail = F)*event2


LL.expert <- apply(df_param,1, function(x){
  expert_log_dens(pllogis(times_expert[1],shape = x[1],scale = x[2],lower.tail = F),
               dfs_expert[[1]])+
    expert_log_dens(pllogis(times_expert[2],shape = x[1],scale = x[2],lower.tail = F),
                 dfs_expert[[2]])})


LL.expert_hat <- 
  expert_log_dens(pllogis(times_expert[1],shape = colMeans(df_param)[1],scale = colMeans(df_param)[2],
                       lower.tail = F),dfs_expert[[1]])+
  expert_log_dens(pllogis(times_expert[2],shape = colMeans(df_param)[1],scale = colMeans(df_param)[2],
                       lower.tail = F),dfs_expert[[2]])

pd_1 <- 2*(sum(LL_hat)+LL.expert_hat -mean(colSums(LL)+LL.expert)) 
DIC.llogis <- -2*(sum(LL_hat)+LL.expert_hat) + 2*pd_1






RPS.example <- "// Royston-Parmar splines model 

functions {
  real rps_lpdf(vector t, vector d, vector gamma, matrix B, matrix DB, vector linpred, vector a0) {
    // t = vector of observed times
    // d = event indicator (=1 if event happened and 0 if censored)
    // gamma = M+2 vector of coefficients for the flexible part
    // B = matrix of basis
    // DB = matrix of derivatives for the basis
    // linpred = fixed effect part
    vector[num_elements(t)] eta;
    vector[num_elements(t)] eta_prime;
    vector[num_elements(t)] log_lik;
    real lprob;
    
    eta = B*gamma + linpred;
    eta_prime = DB*gamma;
    log_lik = d .* (-log(t) + log(eta_prime) + eta) - exp(eta);
    lprob = dot_product(log_lik, a0);
    return lprob;
  }
  
  
  
  
  real Sind( vector gamma, row_vector B, real linpred) {
    // t = vector of observed times
    // gamma = M+2 vector of coefficients for the flexible part
    // B = row_vector of basis
    // linpred = fixed effect part
    real eta;
    real Sind;
    
    eta = B*gamma + linpred;
    Sind = exp(-exp(eta));
    return Sind;
  }
  
  
  real log_density_dist(real[ , ] params, 
                        real x,int num_expert, int pool_type){
    
    // Evaluates the log density for a range of distributions
    
    real dens[num_expert];
    
    for(i in 1:num_expert){
    if(params[i,1] == 1){
      if(pool_type == 1){
        dens[i] = exp(normal_lpdf(x|params[i,3], params[i,4]))*params[i,2]; /// Only require the log density is correct to a constant of proportionality
      }else{
        dens[i] = exp(normal_lpdf(x|params[i,3], params[i,4]))^params[i,2]; /// Only require the log density is correct to a constant of proportionality
      }
      
    }else if(params[i,1] == 2){
       if(pool_type == 1){
          dens[i] = exp(student_t_lpdf(x|params[i,5],params[i,3], params[i,4]))*params[i,2];
        }else{
          dens[i] = exp(student_t_lpdf(x|params[i,5],params[i,3], params[i,4]))^params[i,2];
       }
      
    }else if(params[i,1] == 3){
        if(pool_type == 1){
            dens[i] = exp(gamma_lpdf(x|params[i,3], params[i,4]))*params[i,2];
        }else{
            dens[i] = exp(gamma_lpdf(x|params[i,3], params[i,4]))^params[i,2];
        }  
            
    }else if(params[i,1] == 4){
    
        if(pool_type == 1){
            dens[i] = exp(lognormal_lpdf(x|params[i,3], params[i,4]))*params[i,2];
        }else{
             dens[i] = exp(lognormal_lpdf(x|params[i,3], params[i,4]))^params[i,2];
        }
        
    }else if(params[i,1] == 5){
     if(pool_type == 1){
            dens[i] = exp(beta_lpdf(x|params[i,3], params[i,4]))*params[i,2];
        }else{
            dens[i] = exp(beta_lpdf(x|params[i,3], params[i,4]))^params[i,2];
        }
      
      
      } 

    }
    
      
    return(log(sum(dens)));
    
  }
  
}

data {
  int<lower=1> n;                   // number of observations
  int<lower=0> M;                   // number of internal knots for the splines model
  int<lower=1> H;                   // number of covariates in the (time-independent) linear predictor
  vector<lower=0>[n] t;             // observed times (including censored values)
  vector<lower=0,upper=1>[n] d;     // censoring indicator: 1 if fully observed, 0 if censored
  matrix[n,H] X;                    // matrix of covariates for the (time-independent) linear predictor
  matrix[n,M+2] B;                  // matrix with basis
  matrix[n,M+2] DB;                 // matrix with derivatives of the basis
  vector[H] mu_beta;                // mean of the covariates coefficients
  vector<lower=0> [H] sigma_beta;   // sd of the covariates coefficients
  vector[M+2] mu_gamma;             // mean of the splines coefficients
  vector<lower=0>[M+2] sigma_gamma; // sd of the splines coefficients
  
  int n_time_expert;
  int<lower = 0, upper = 1> St_indic; // 1 Expert opinion on survival @ timepoint ; 0 Expert opinion on survival difference
  
  int id_St;
  int id_trt;
  int id_comp;
  
  int n_experts[n_time_expert];
  int pool_type;

  real param_expert[max(n_experts),5,n_time_expert];
  vector[St_indic ? n_time_expert : 0] time_expert;
  
  matrix[n_time_expert,M+2] B_expert;                  // matrix with basis for experts times
  vector[n] a0; //Power prior for the observations  
}


parameters {
  vector[M+2] gamma;
  vector[H] beta;
}


transformed parameters{
  vector[n] linpred;
  vector[n] mu;
  vector[n_time_expert] St_expert;
  
  linpred = X*beta;
  for (i in 1:n) {
    mu[i] = linpred[i];
  }
  for (i in 1:n_time_expert){
    St_expert[i] = Sind(gamma, row(B_expert,i),mu[id_St[1]]);
    
  }	
  
}

model {
  // Priors
  gamma ~ normal(mu_gamma,sigma_gamma);
  beta ~ normal(mu_beta,sigma_beta);
  
  // Data model
  t ~ rps(d,gamma,B,DB,X*beta, a0);
  
  for (i in 1:n_time_expert){
     
     target += log_density_dist(param_expert[,,i],
                                 St_expert[i],
                                 n_experts[i],
                                 pool_type);
  }
}"

names(m.all$models)
model.data <- m.all$misc$data.stan[[5]]
model.data$param_expert <-array_final



fit.RPS <- stan(model_code = RPS.example, model_name = "RPS Example", 
                        data =  model.data, iter = 2000, chains = 3)


df_param <- data.frame(extract(fit.RPS, "gamma"))

k1 <- 1

knots1 <- quantile(log((digitized_IPD %>% filter(event == 1))$time),
                   seq(0, 1, length = k1 + 2))
LL <- apply(df_param,1, 
            function(x){dsurvspline(t,gamma = x,knots = knots1,log = T)*event+
                psurvspline(t,gamma = x,knots = knots1,log = T,lower.tail = F)*event2    
            })

LL_hat <- dsurvspline(t,gamma = colMeans(df_param),knots = knots1,log = T)*event+
  psurvspline(t,gamma = colMeans(df_param),knots = knots1,log = T,lower.tail = F)*event2    


LL.expert <- apply(df_param,1, function(x){
  expert_log_dens( psurvspline(times_expert[1],gamma = x,knots = knots1,lower.tail = F),
                dfs_expert[[1]])+
    expert_log_dens(psurvspline(times_expert[2],gamma = x,knots = knots1,lower.tail = F),
                 dfs_expert[[2]])})


LL.expert_hat <- 
  expert_log_dens( psurvspline(times_expert[1],gamma = colMeans(df_param),knots = knots1,lower.tail = F),
                dfs_expert[[1]])+
  expert_log_dens(psurvspline(times_expert[2],gamma = colMeans(df_param),knots = knots1,lower.tail = F),
               dfs_expert[[2]])


pd_1 <- 2*(sum(LL_hat)+LL.expert_hat -mean(colSums(LL)+LL.expert)) 
DIC.rps <- -2*(sum(LL_hat)+LL.expert_hat) + 2*pd_1

