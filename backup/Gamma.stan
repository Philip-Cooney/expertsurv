// Gamma survival model
functions {
	
 
  real gamma_St(vector t, real alpha,  vector beta, vector a0) {
    vector[num_elements(t)] log_lik;
	vector[num_elements(t)] prob;
    real lprob;
    // Constructs the log-density for each observation
	// don't use gamma_lpdf, use the function derived by yourself
	
    for (i in 1:num_elements(t)) {
      prob[i] = log(1-gamma_cdf(t[i], alpha, beta[i]));
	  //prob[i] = alpha*log(beta[i])-lgamma(alpha)+ (alpha-1)*log(t[i]) - beta*t[i]
	  
    }
    // And the total log-density (as a sum of the individual terms)
    lprob = dot_product(prob, a0);
    return lprob;
  }
  
  
   real gamma_dens(vector t, real alpha,  vector beta, vector a0) {
    vector[num_elements(t)] log_lik;
	vector[num_elements(t)] prob;
    real lprob;
    // Constructs the log-density for each observation
	// don't use gamma_lpdf, use the function derived by yourself
	
    for (i in 1:num_elements(t)) {
      //prob[i] = gamma_lpdf(t[i]| alpha, beta[i]);
	  prob[i] = alpha*log(beta[i])-lgamma(alpha)+ (alpha-1)*log(t[i]) - beta[i]*t[i];
	  
    }
    // And the total log-density (as a sum of the individual terms)
    lprob =  dot_product(prob, a0);
    return lprob;
  }
  
  
  
	
}

data {
  int<lower=1> n_obs;                     // number of observed cases
  int<lower=0> n_cens;                    // number of censored cases
  vector<lower=0>[n_obs] t;               // fully observed times
  vector<lower=0>[n_cens] d;              // observed censoring times
  int<lower=1> H;                         // number of covariates (including intercept)
  matrix[n_obs,H] X_obs;                  // matrix of categorical covariates for the valid cases (0/1 => dummy variables)
  matrix[n_cens,H] X_cens;                // matrix of categorical covariates for the censored cases (0/1 => dummy variables)
  vector[H] mu_beta;                      // vector of means for the covariates
  vector<lower=0>[H] sigma_beta;          // vector of sd for the covariates
  real<lower=0> a_alpha;		  		        // first parameter for the shape distribution
  real<lower=0> b_alpha;		  		        // second parameter for the shape distribution
  
  int n_time_expert;
  vector[n_time_expert] mu_St;
  vector[n_time_expert] sigma_St;
  vector[n_time_expert] time_expert;
  int id_St;
  
  int St_indic; // 1 Expert opinion on survival @ timepoint ; 2 Expert opinion on survival difference
  int id_trt;
  int id_comp;
  real mu_diff;
  real<lower = 0> sigma_diff;
  
  
  vector<lower=0>[n_obs] a0_obs;              
  vector<lower=0>[n_cens]a0_cens ;              

  
}

parameters {
  real<lower=0> alpha;                    // shape of the Gamma distribution
  vector[H] beta;                         // coefficients for the covariates
  //vector<lower=1>[n_cens] cens;           // censoring variable (latent)
}

transformed parameters {
  vector[n_obs] loglambda_obs;            // loglinear predictor for the observed cases
  vector[n_cens] loglambda_cens;          // loglinear predictor for the censored cases
  vector[n_obs] lambda_obs;               // rescaled predictor (rate) for the observed cases
  vector[n_cens] lambda_cens;             // rescaled predictor (rate) for the censored cases
  vector[n_time_expert] St_expert;
  real Surv_diff_val;
 

  loglambda_cens = X_cens*beta;
  for (i in 1:n_cens) {
    lambda_cens[i] = exp(loglambda_cens[i]);
  }
  loglambda_obs = X_obs*beta;
  for (i in 1:n_obs) {
    lambda_obs[i] = exp(loglambda_obs[i]);
  }
      
  if(St_indic == 1){
	for (i in 1:n_time_expert){
		St_expert[i] = 1-gamma_cdf(time_expert[i], alpha, lambda_obs[id_St]); 
	 }
	}else{    
	Surv_diff_val = alpha/lambda_obs[id_trt] -alpha/lambda_obs[id_comp];
  	}  
}

model {
  // Prior distributions
  alpha ~ gamma(a_alpha,b_alpha);
  beta ~ normal(mu_beta,sigma_beta);
  // Data model
  if(n_cens>0){
   target += gamma_St(d, alpha, lambda_cens,a0_cens);
  }
  target += gamma_dens(t, alpha,lambda_obs,a0_obs);
  
  if (St_indic ==1) {
	for (i in 1:n_time_expert){
	target += normal_lpdf(St_expert[i] | mu_St[i], sigma_St[i]);
	}
  }else{
    target += normal_lpdf(Surv_diff_val | mu_diff, sigma_diff);
	}
 }

generated quantities {
  real rate;
  rate = exp(beta[1]);
}
