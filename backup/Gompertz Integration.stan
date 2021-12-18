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
 

  //In order to integrate we need this custom integrand function defined exactly as below 
real integrand(real x,real xc,real[] theta,	real[] x_r,	int[] x_i ){

// separate the parameters
	real shape = theta[1];
	real rate = theta[2];

return exp(-rate/shape * (exp(shape * x) - 1));

}

  real expec_surv(real[] theta, data real[] x_r) {
		int x_i[0];
		return integrate_1d(integrand, 0, positive_infinity(),
		theta, x_r, x_i, 1e-8);
}  


  // Defines difference in expected survival
  real Surv_diff ( real shape, real rate_trt, real rate_comp, data real[] x_r) {
	real Surv_diff;
	real Surv_trt;
	real Surv_comp;
	real theta_trt[2];
	real theta_comp[2];
	
	theta_trt[1] = shape;
	theta_trt[2] = rate_trt;
	theta_comp[1] = shape;
	theta_comp[2] = rate_comp;
	
	Surv_trt = expec_surv(theta_trt, x_r);
	Surv_comp = expec_surv(theta_comp, x_r);
	
    Surv_diff = Surv_trt - Surv_comp;
    return Surv_diff;
  }
  
  
  // Defines the sampling distribution
  real surv_gompertz_lpdf (vector t, vector d, real shape, vector rate) {
    vector[num_elements(t)] log_lik;
    real prob;
    log_lik = d .* log_h(t,shape,rate) + log_S(t,shape,rate);
    prob = sum(log_lik);
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
  
  int St_indic; // 1 Expert opinion on survival @ timepoint ; 2 Expert opinion on survival difference
  int id_trt;
  int id_comp;
  real mu_diff;
  real<lower = 0> sigma_diff; 
  
  
}


transformed data {
  real x_r[0];

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
  
    if(St_indic == 1){
  for (i in 1:n_time_expert){
    St_expert[i] = Sind(time_expert[i],alpha,mu[id_St]);
	}
  }else{    
	Surv_diff_val = Surv_diff(alpha, mu[id_trt], mu[id_comp],x_r);
  }
}
  
model {
  alpha ~ gamma(a_alpha,b_alpha);
  beta ~ normal(mu_beta,sigma_beta);
  t ~ surv_gompertz(d,alpha,mu);
  
    if (St_indic ==1) {
	for (i in 1:n_time_expert){
	target += normal_lpdf(St_expert[i] | mu_St[i], sigma_St[i]);
	}
  }else{
    target += normal_lpdf(Surv_diff_val | mu_diff, sigma_diff);
	} 
  
}

generated quantities {
  real rate;                        // rate parameter
  rate = exp(beta[1]);
}