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
    real Sind_rtn;

    eta = B*gamma + linpred;
    Sind_rtn = exp(-exp(eta));
    return Sind_rtn;
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
  int<lower = 0, upper = 1> St_indic; // 1 Expert opinion on survival @ timepoint ; 0 Expert opinion on survival difference

  int id_St;
  int id_trt;
  int id_comp;

  int n_experts[n_time_expert];
  int pool_type;

  real param_expert[max(n_experts),5,n_time_expert];
  vector[St_indic ? n_time_expert : 0] time_expert;

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

     target += log_density_dist(param_expert[,,i],
                                 St_expert[i],
                                 n_experts[i],
                                 pool_type);
  }
}

