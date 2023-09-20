// log-Logistic survival model

functions {
  // Defines the log hazard
  vector log_h (vector t, real shape, vector scale) {
    vector[num_elements(t)] log_h_rtn;
    for (i in 1:num_elements(t)) {
      log_h_rtn[i] = log(shape)-log(scale[i])+(shape-1)*(log(t[i])-log(scale[i]))-log(1+pow((t[i]/scale[i]),shape));
    }
    return log_h_rtn;
  }

  // Defines the log survival
  vector log_S (vector t, real shape, vector scale) {
    vector[num_elements(t)] log_S_rtn;
    for (i in 1:num_elements(t)) {
      log_S_rtn[i] = -log(1+pow((t[i]/scale[i]),shape));
    }
    return log_S_rtn;
  }

  // Defines the  survival indvidual
  real Sind (real t, real shape, real scale) {
    real Sind_rtn;
    Sind_rtn = exp(-log(1+pow((t/scale),shape)));
    return Sind_rtn;
  }

  // Defines difference in expected survival
  real Surv_diff ( real shape, real scale_trt, real scale_comp) {
    real Surv_diff_rtn;
    real b;
    b = pi()/shape;
    Surv_diff_rtn = (scale_trt-scale_comp)*b/sin(b);
    return Surv_diff_rtn;
  }



  // Defines the sampling distribution
  real surv_loglogistic_lpdf (vector t, vector d, real shape, vector scale, vector a0) {
    vector[num_elements(t)] log_lik;
    real prob;
    log_lik = d .* log_h(t,shape,scale) + log_S(t,shape,scale);
    prob = dot_product(log_lik,a0);
    return prob;
  }


   real log_density_dist(array[ , ] real params,
                        real x,int num_expert, int pool_type){

    // Evaluates the log density for a range of distributions

    array[num_expert] real dens;

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
  vector[H] mu_beta;	  // mean of the covariates coefficients
  vector<lower=0> [H] sigma_beta;   // sd of the covariates coefficients
  real<lower=0> a_alpha;
  real<lower=0> b_alpha;
  vector[n] a0; //Power prior for the observations

  int<lower = 0, upper = 1> St_indic; // 1 Expert opinion on survival @ timepoint ; 0 Expert opinion on survival difference
  int n_time_expert;
  int id_St;
  int id_trt;
  int id_comp;

  array[n_time_expert] int n_experts;
  int pool_type;

  array[max(n_experts),5,n_time_expert] real param_expert;
  vector[St_indic ? n_time_expert : 0] time_expert;


}

parameters {
  vector[H] beta;         // Coefficients in the linear predictor (including intercept)
  real<lower=1> alpha;    // shape parameter - Constrainted to be greater than 1 so that the expected value exists
}

transformed parameters {
  vector[n] linpred;
  vector[n] mu;
  vector[n_time_expert] St_expert;

  linpred = X*beta;
  for (i in 1:n) {
    mu[i] = exp(linpred[i]);
  }

  for (i in 1:n_time_expert){
    if(St_indic == 1){

      St_expert[i] = Sind(time_expert[i],alpha,mu[id_St]);

    }else{
      St_expert[i] = Surv_diff(alpha,mu[id_trt],mu[id_comp]);
    }
  }

}

model {
  alpha ~ gamma(a_alpha,b_alpha);
  beta ~ normal(mu_beta,sigma_beta);
  t ~ surv_loglogistic(d,alpha,mu, a0);

  for (i in 1:n_time_expert){

     target += log_density_dist(param_expert[,,i],
                                 St_expert[i],
                                 n_experts[i],
                                 pool_type);
  }

}

generated quantities {
  real rate;                // rate parameter
  rate = exp(beta[1]);
}
