// Weibull survival model (with PH parameterisation)

functions {
  // Defines the log hazard
  // Defines the log hazard
  vector log_haz (vector t, real shape, vector scale) {
    vector[num_elements(t)] log_haz;
    //log_haz = log(shape)+(shape-1)*log(t ./ scale)-log(scale);
        for(i in 1:num_elements(t)){
          log_haz[i] = log(shape) +log(scale[i])+  (shape-1)*log(t[i]);
        }
    return log_haz;
  }
  
  // Defines the log survival
  vector log_S (vector t, real shape, vector scale) {
    vector[num_elements(t)] log_S;
    for (i in 1:num_elements(t)) {
      log_S[i] = -scale[i]*pow(t[i],shape);
    }
    return log_S;
  }
  
  
    // Defines the log survival indvidual
  real log_Sind (real t, real shape, real scale) {
	real log_Sind;
      log_Sind = -scale*pow(t,shape);
    return log_Sind;
  }
  
  
      // Defines difference in expected survival
  real Surv_diff ( real shape, real scale_trt, real scale_comp) {
	real Surv_diff;
	real scale_trt_mod;
	real scale_comp_mod;
	
	scale_trt_mod = pow(scale_trt,-1/shape);
	scale_comp_mod = pow(scale_comp,-1/shape);
	
  Surv_diff = (scale_trt_mod-scale_comp_mod)*tgamma(1 +1/shape);
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
  vector[H] beta_lower;	    //lower bound for the covariate coefficients
  vector[H] beta_upper;   // upper bound for the covariate coefficients
  real<lower=0> alpha_lower;
  real<lower=0> alpha_upper;
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

  real St_lower;
  real St_upper;

  int<lower = 0, upper = 1> expert_plus_data;

}

parameters {
  vector[H] beta_exp;         // Coefficients in the linear predictor (including intercept)
  real<lower=0> alpha;    // shape parameter
}

transformed parameters {
  vector[H] beta;
  vector[n] linpred;
  vector[n] mu;
  vector<lower=St_lower,upper=St_upper>[n_time_expert] St_expert;
 
  for(i in 1:H){
    beta[i] = log(beta_exp[i]); // Transformation to make it linear on the outcome scale (i.e. exp)
  }


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
  alpha ~ uniform(alpha_lower,alpha_upper);
   beta_exp ~ uniform(beta_lower, beta_upper);
  
  if(expert_plus_data){
     t ~ surv_weibullAF(d,alpha,mu, a0);
  }

  
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
}
