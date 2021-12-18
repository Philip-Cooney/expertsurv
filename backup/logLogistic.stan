// log-Logistic survival model

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
  
      // Defines the log survival indvidual
  real Sind (real t, real shape, real scale) {
	real Sind;
      Sind = -log(1+pow((t/scale),shape));
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
  real<lower=1> alpha;    // shape parameter - Constrainted to be greater than 1 so that the expected value exists
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
  
    if(St_indic == 1){
  for (i in 1:n_time_expert){
    St_expert[i] = Sind(time_expert[i],alpha,mu[id_St]);
	}
  }else{    
	Surv_diff_val = Surv_diff(alpha,mu[id_trt],mu[id_comp]);  
  }
  
}

model {
  alpha ~ gamma(a_alpha,b_alpha);
  beta ~ normal(mu_beta,sigma_beta);
  t ~ surv_loglogistic(d,alpha,mu, a0);
  
    if (St_indic ==1) {
	for (i in 1:n_time_expert){
	target += normal_lpdf(St_expert[i] | mu_St[i], sigma_St[i]);
	}
  }else{
    target += normal_lpdf(Surv_diff_val | mu_diff, sigma_diff);
	}
  
}

generated quantities {
  real rate;                // rate parameter
  rate = exp(beta[1]);
}