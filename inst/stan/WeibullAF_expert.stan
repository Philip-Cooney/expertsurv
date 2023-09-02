// Weibull survival model (with AFT parameterisation)

functions {
  // Defines the log hazard
  vector log_haz (vector t, real shape, vector scale) {
    vector[num_elements(t)] log_haz_rtn;
    log_haz_rtn = log(shape)+(shape-1)*log(t ./ scale)-log(scale);
    return log_haz_rtn;
  }

  // Defines the log survival
  vector log_S (vector t, real shape, vector scale) {
    vector[num_elements(t)] log_S_rtn;
    for (i in 1:num_elements(t)) {
      log_S_rtn[i] = -pow((t[i]/scale[i]),shape);
    }
    return log_S_rtn;
  }


    // Defines the log survival indvidual
  real log_Sind (real t, real shape, real scale) {
	real log_Sind_rtn;
      log_Sind_rtn = -pow((t/scale),shape);
    return log_Sind_rtn;
  }


      // Defines difference in expected survival
  real Surv_diff ( real shape, real scale_trt, real scale_comp) {
	real Surv_diff_rtn;
      Surv_diff_rtn = (scale_trt-scale_comp)*tgamma(1 +1/shape);
    return Surv_diff_rtn;
  }

  // Defines the sampling distribution
  real surv_weibullAF_lpdf (vector t, vector d, real shape, vector scale, vector a0) {
    vector[num_elements(t)] log_lik;
    real prob;
    log_lik = d .* log_haz(t,shape,scale) + log_S(t,shape,scale);
    prob = dot_product(log_lik, a0);
    return prob;
  }

  // real log_mixnorm_dens(real[ ] mu, real[ ] sd, real[ ] prob,  real x) {
  //
  //   // Evaluates the log-density based on a mixture normal
  //   // x refers to the Survival or expected difference in survival
  //
  //   //real log_dens;
  //   int n;
  //   vector[num_elements(mu)] log_dens_i;
  //
  //   n = num_elements(mu);
  //   for(i in 1:n){
  //     log_dens_i[i]= normal_lpdf(x| mu[i], sd[i]) + log(prob[i]);
  //   }
  //
  //   return log(sum(exp(log_dens_i)));
  // }
  //


     // Evaluates the log density for a range of distributions
  real log_density_dist(array[ , ] real params,
                        real x,int num_expert, int pool_type){


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
  vector[H] mu_beta;	    // mean of the covariates coefficients
  vector<lower=0> [H] sigma_beta;   // sd of the covariates coefficients
  real<lower=0> a_alpha;
  real<lower=0> b_alpha;

  int n_time_expert;
  int<lower = 0, upper = 1> St_indic; // 1 Expert opinion on survival @ timepoint ; 0 Expert opinion on survival difference

  int id_St;
  int id_trt;
  int id_comp;

  array[n_time_expert] int n_experts;
  int pool_type;

  array[max(n_experts),5,n_time_expert] real param_expert;
  vector[St_indic ? n_time_expert : 0] time_expert;


   vector[n] a0; //Power prior for the observations

}

parameters {
  vector[H] beta;         // Coefficients in the linear predictor (including intercept)
  real<lower=0> alpha;    // shape parameter

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

    St_expert[i] = exp(log_Sind(time_expert[i],alpha,mu[id_St]));
    }else
	  St_expert[i] = Surv_diff(alpha,mu[id_trt],mu[id_comp]);

  }


}

model {
  alpha ~ gamma(a_alpha,b_alpha);
  beta ~ normal(mu_beta,sigma_beta);
  t ~ surv_weibullAF(d,alpha,mu, a0);

for (i in 1:n_time_expert){
     target += log_density_dist(param_expert[,,i],
                                 St_expert[i],
                                 n_experts[i],
                                 pool_type);
    }
}
generated quantities {
  real scale;                // scale parameter
  scale = exp(beta[1]);

 }

