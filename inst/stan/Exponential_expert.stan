// Exponential survival model

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
	   St_expert[i] = Surv_diff(mu[id_trt],mu[id_comp]);
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
}
