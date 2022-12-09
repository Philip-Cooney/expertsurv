// log-Normal survival model

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
    
      
      
      if(pool_type == 1){
      return(log(sum(dens)));
      }else{
      return(log(prod(dens)));
      }
    
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
      
      St_expert[i] = exp(log_Sind(time_expert[i],mu[id_St],alpha));
      
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
}
