// Royston-Parmar splines model 

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
  vector[n_time_expert] mu_St;
  vector[n_time_expert] sigma_St;
  vector[n_time_expert] time_expert;
  int id_St;
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
    St_expert[i] = Sind(gamma, row(B_expert,i),mu[id_St]);
	
 }	
	
}

model {
  // Priors
  gamma ~ normal(mu_gamma,sigma_gamma);
  beta ~ normal(mu_beta,sigma_beta);
  
  // Data model
  t ~ rps(d,gamma,B,DB,X*beta, a0);
  
  	for (i in 1:n_time_expert){
	target += normal_lpdf(St_expert[i] | mu_St[i], sigma_St[i]);
	}
 }