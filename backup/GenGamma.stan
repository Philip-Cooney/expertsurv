// Generalised Gamma model for survival analysis - uses the parameterisation of flexsurv 

// First defines the Generalised Gamma model for fully observed and censored cases

functions {
  real gen_gamma_lpdf(vector x, vector mu, real sigma, real Q, vector a0) {
    // Uses the same parameterisation as flexsurv
    // mu = location
    // sigma = scale
    // Q = shape
    vector[num_elements(x)] prob;
    real lprob;
    vector[num_elements(x)] w;
    // Constructs the log-density for each observation
    w = ((log(x)-mu))/sigma;
    for (i in 1:num_elements(x)) {
      prob[i] = -log(sigma*x[i])+log(fabs(Q))+pow(Q,-2)*log(pow(Q,-2))+pow(Q,-2)*(Q*w[i]-exp(Q*w[i]))-lgamma(pow(Q,-2));
    }
    // And the total log-density (as a sum of the individual terms)
    lprob = dot_product(prob, a0);
    return lprob;
  }
  
  real gen_gamma_cens_lpdf(vector x, vector mu, real sigma, real Q, vector u, vector a0) {
    // Rescales the distribution accounting for right censoring
    vector[num_elements(x)] prob;
    real lprob;
    vector[num_elements(x)] w;
    vector[num_elements(x)] tr;
    // Constructs the log-density for each observation
    tr = x .* u;
    w = ((log(tr)-mu))/sigma;
    for (i in 1:num_elements(x)) {
      prob[i] = log(u[i])-log(sigma*tr[i])+log(fabs(Q))+pow(Q,-2)*log(pow(Q,-2))+pow(Q,-2)*(Q*w[i]-exp(Q*w[i]))-lgamma(pow(Q,-2));
    }
    // And the total log-density (as a sum of the individual terms)
    lprob = dot_product(prob, a0);
    return lprob;
  }
  
  
  real Sind(real time, real mu, real sigma, real Q){
	real Sind;
	real qi;
	real w;
	real expnu;
	
	qi = 1/(Q*Q);
	w = (log(time)-mu)/sigma;
	expnu = exp(Q*w)*qi;
	
	Sind = 1- gamma_cdf(expnu,qi,1);
	return Sind;
  
  }
  
    // Defines difference in expected survival
  real Surv_diff ( real sigma,real Q, real mu_trt, real mu_comp) {
	real Surv_diff;
	real a_trt;
	real a_comp;
	real p;
	real d;
	real mean_trt;
	real mean_comp;
	
	a_trt = exp(mu_trt + (sigma*2*log(Q))/Q);
	a_comp = exp(mu_comp + (sigma*2*log(Q))/Q);
	p = Q/sigma;
    d = 1/(sigma*Q);
	
	mean_trt = (a_trt*tgamma((d+1)/p))/(tgamma(d/p));
	mean_comp = (a_comp*tgamma((d+1)/p))/(tgamma(d/p));
	
	Surv_diff = mean_trt - mean_comp;
    return Surv_diff;
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
  real mu_Q;                              // mean for the parameter Q
  real<lower=0> sigma_Q;                  // sd for the parameter Q
  real<lower=0> a_sigma;                  // first parameter for the scale distribution
  real<lower=0> b_sigma;                  // second parameter for the scale distribution
  
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
  real Q;                                 // shape of the Generalised Gamma distribution
  real<lower=0> sigma;                    // scale of the Generalised Gamma distribution
  vector[H] beta;                         // coefficients for the covariates
  vector<lower=1>[n_cens] cens;           // censoring variable (latent)
 
}

transformed parameters{
  vector[n_obs] linpred_obs;               // rescaled predictor (mu) for the observed cases
  vector[n_cens] linpred_cens;             // rescaled predictor (mu) for the censored cases
  vector[n_time_expert] St_expert;
  real Surv_diff_val;

  linpred_cens = X_cens*beta;
  linpred_obs = X_obs*beta;
      
  if(St_indic == 1){
	for (i in 1:n_time_expert){
		St_expert[i] = Sind(time_expert[i],linpred_obs[id_St],sigma,Q); 
	 }
	}else{    
	Surv_diff_val = Surv_diff(sigma,Q,linpred_obs[id_trt],linpred_obs[id_comp]);
  	}  

}

model {
  // Prior distributions
  Q ~ normal(mu_Q,sigma_Q);
  sigma ~ gamma(a_sigma,b_sigma);
  beta ~ normal(mu_beta,sigma_beta);
  // Data model
    if(n_cens>0){
  cens ~ gen_gamma_cens(X_cens*beta,sigma,Q,d, a0_cens);
   }
  t ~ gen_gamma(X_obs*beta,sigma,Q, a0_obs);
  
  if (St_indic ==1) {
	for (i in 1:n_time_expert){
	target += normal_lpdf(St_expert[i] | mu_St[i], sigma_St[i]);
	}
  }else{
    target += normal_lpdf(Surv_diff_val | mu_diff, sigma_diff);
	}   
}

generated quantities {
  real mu;
  mu = beta[1];
}