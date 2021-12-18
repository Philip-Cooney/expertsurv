// Gompertz survival model

functions {
  // Defines the log hazard
  vector log_h (vector t, real shape, vector rate) {
    vector[num_elements(t)] log_h;
    log_h = log(rate) + (shape * t);
    return log_h;
  }
  
  // Defines the log survival
  vector log_S (vector t, real shape, vector rate) {
    vector[num_elements(t)] log_S;
    for (i in 1:num_elements(t)) {
      log_S[i] = -rate[i]/shape * (exp(shape * t[i]) - 1);
    }
    return log_S;
  }
  
    // Defines the log survival indvidual
  real Sind (real t, real shape, real scale) {
	real Sind;
      Sind = exp(-pow((t/scale),shape));
    return Sind;
  }
 
  // Defines the sampling distribution
  real surv_gompertz_lpdf (vector t, vector d, real shape, vector rate, vector a0) {
    vector[num_elements(t)] log_lik;
    real prob;
    log_lik = d .* log_h(t,shape,rate) + log_S(t,shape,rate);
    prob = dot_product(log_lik, a0);
    return prob;
  }
  
  
}

data {
  int n;                          // number of observations
  vector[n] t;                    // observed times
  vector[n] d;                    // censoring indicator (1=observed, 0=censored)
  int H;                          // number of covariates
  matrix[n,H] X;                  // matrix of covariates (with n rows and H columns)
  vector[H] mu_beta;	            // mean of the covariates coefficients
  vector<lower=0> [H] sigma_beta; // sd of the covariates coefficients
  real<lower=0> a_alpha;          // mu_alpha
  real<lower=0> b_alpha;          // sigma_alpha
  
  int n_time_expert;
  vector[n_time_expert] mu_St;
  vector[n_time_expert] sigma_St;
  vector[n_time_expert] time_expert;
  int id_St;
  vector[n] a0; //Power prior for the observations  
}

parameters {
  vector[H] beta;                 // Coefficients in the linear predictor (including intercept)
  real<lower=0> alpha;            // shape parameter
}

transformed parameters {
  vector[n] linpred;
  vector[n] mu;
  vector[n_time_expert] St_expert;
  real Surv_diff_val;
  
  linpred = X*beta;
  for (i in 1:n) {
    mu[i] = exp(linpred[i]);
  }
  
    for (i in 1:n_time_expert){
    St_expert[i] = Sind(time_expert[i],alpha,mu[id_St]);
	
	}

}
  
model {
  alpha ~ gamma(a_alpha,b_alpha);
  beta ~ normal(mu_beta,sigma_beta);
  t ~ surv_gompertz(d,alpha,mu, a0);
  
	for (i in 1:n_time_expert){
	target += normal_lpdf(St_expert[i] | mu_St[i], sigma_St[i]);
	}

  
}

generated quantities {
  real rate;                        // rate parameter
  rate = exp(beta[1]);
}