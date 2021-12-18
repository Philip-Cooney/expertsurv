// Generalised F model for survival analysis - uses the parameterisation of flexsurv 

// First defines the Generalised F model for fully observed and censored cases
functions {
    real genf_lpdf(vector x, vector mu, real sigma, real Q, real P,vector a0) {
      // mu = location (\in R)
      // sigma = scale (>0)
      // Q = first shape (\in R)
      // P = second shape (>0)
      vector[num_elements(x)] prob;
      real lprob;
      real tmp;
      real delta;
      real s1; 
      real s2;
      vector[num_elements(x)] expw;
      tmp = pow(Q, 2) + 2*P;
      delta = pow(tmp, 0.5);
      s1 = 2/(tmp + Q*delta);
      s2 = 2/(tmp - Q*delta);
      for (i in 1:num_elements(x)) {
        expw[i] = pow(x[i],delta/sigma)*exp(-mu[i]*delta/sigma);
        prob[i] = log(delta) + s1/sigma * delta * (log(x[i]) - mu[i]) + s1 * (log(s1) - log(s2)) - log(sigma * x[i]) - 
                  (s1 + s2) * log(1 + s1 * expw[i]/s2) - lbeta(s1, s2);
      }
      lprob = dot_product(prob, a0);
      return lprob;
    }
    

	
	real Sind(real time, real mu, real sigma, real Q, real P){
	real tmp;
	real delta;
	real s1;
	real s2;
	real qb;
	real expw;
	real Sind;
	
	tmp = Q * Q + 2 * P;
	delta = sqrt(tmp);
	s1 = 2 / (tmp + Q*delta);
	s2 = 2 / (tmp - Q*delta);
	expw = pow(time,(delta/sigma))*exp(-mu*delta/sigma);
	qb = s2/(s2 + s1*expw);
	Sind = beta_cdf(qb,s2,s1);
	return Sind;
		
    }	
	
	
	real genf_cens_lpdf(vector x, vector mu, real sigma, real Q, real P,  vector a0) {
      vector[num_elements(x)] prob;
      real lprob;
     
	 for (i in 1:num_elements(x)) {
        
        prob[i] = log(Sind(x[i], mu[i], sigma, Q, P));
      }
      lprob=dot_product(prob, a0);
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
  real<lower=0> a_sigma;		 		          // first parameter for the shape parameter
  real<lower=0> b_sigma;		 		          // second parameter for the shape parameter
  real mu_P;						                  // mean for the parameter logP
  real<lower=0> sigma_P;		  		        // sd for the parameter logP
  real mu_Q; 				  			              // mean for the parameter Q
  real<lower=0> sigma_Q;		 		          // sd for the parameter Q
  
  int n_time_expert;
  vector[n_time_expert] mu_St;
  vector[n_time_expert] sigma_St;
  vector[n_time_expert] time_expert;
  int id_St;
    
  
  vector<lower=0>[n_obs] a0_obs;              
  vector<lower=0>[n_cens]a0_cens ;              

}

parameters {
  real<lower=0> sigma;                    // scale of the Generalised F distribution
  real Q;                                 // first degree of freedom
  real logP;                              // second degree of freedom (log-scale)
  vector[H] beta;                         // coefficients for the covariates
  //vector<lower=1>[n_cens] cens;           // censoring variable (latent)
}

transformed parameters {
  vector[n_obs] linpred_obs;               // rescaled predictor (mu) for the observed cases
  vector[n_cens] linpred_cens;             // rescaled predictor (mu) for the censored cases
  vector[n_time_expert] St_expert;
 
  real<lower=0> P;                        // second degree of freedom (natural scale)
  P = exp(logP);
  
  linpred_cens = X_cens*beta;
  linpred_obs = X_obs*beta;
  
  	for (i in 1:n_time_expert){
		St_expert[i] = Sind(time_expert[i],linpred_obs[id_St],sigma,Q,P); 
	 }
  
}

model {
  // Prior distributions
  sigma ~ gamma(a_sigma,b_sigma);
  logP ~ normal(mu_P,sigma_P);
  Q ~ normal(mu_Q,sigma_Q);
  beta ~ normal(mu_beta,sigma_beta);
  // Data model
   if(n_cens>0){
  d ~ genf_cens(X_cens*beta,sigma,Q,P, a0_cens);
    }
  t ~ genf(X_obs*beta,sigma,Q,P, a0_obs);
  
  	for (i in 1:n_time_expert){
	target += normal_lpdf(St_expert[i] | mu_St[i], sigma_St[i]);
	}
  
}

generated quantities {
  real mu;
  mu = beta[1];
}