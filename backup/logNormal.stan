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
  real Surv_diff ( real sd, real mean_trt, real mean_comp) {
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
}

parameters {
  vector[H] beta;         // Coefficients in the linear predictor (including intercept)
  real<lower=0> alpha;    // log-sd parameter
}

transformed parameters {
  vector[n] linpred;
  vector[n] mu;
  vector[n_time_expert] St_expert;
  real Surv_diff_val;
  
  linpred = X*beta;
  for (i in 1:n) {
    mu[i] = linpred[i];
  }
  
  if(St_indic == 1){
  for (i in 1:n_time_expert){
    St_expert[i] = exp(log_Sind(time_expert[i],alpha,mu[id_St]));
	}
  }else{    
	Surv_diff_val = Surv_diff(alpha,mu[id_trt],mu[id_comp]);  
  }
}

model {
  alpha ~ uniform(a_alpha,b_alpha);
  beta ~ normal(mu_beta,sigma_beta);
  t ~ surv_lognormal(d,X*beta,alpha, a0);
  
  if (St_indic ==1) {
	for (i in 1:n_time_expert){
	target += normal_lpdf(St_expert[i] | mu_St[i], sigma_St[i]);
	}
  }else{
    target += normal_lpdf(Surv_diff_val | mu_diff, sigma_diff);
	}
}

generated quantities {
  real meanlog;
  meanlog = beta[1];
}