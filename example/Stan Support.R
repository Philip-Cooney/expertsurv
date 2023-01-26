library(rstan)
library(flexsurv)
## gengamm and GEN F don't work with a0 or defining the surival function
GenGamma.model <- "// Generalised Gamma model for survival analysis - uses the parameterisation of flexsurv 
//https://dev.to/martinmodrak/optional-parametersdata-in-stan-4o33

// First defines the Generalised Gamma model for fully observed and censored cases

functions {

 real gamma_cdf_man(real s, real x,int k){
  real part1;
  vector[k+1] res;
  
  part1 = x^s*(tgamma(s))*exp(-x);
  
  for(i in 0:k){
    res[i+1] =  x^i/(tgamma(s+i+1));
  }

  return part1*sum(res)/tgamma(s);
  
  }

  real gen_gamma_lden(real x, real mu, real sigma, real Q) {
    // Uses the same parameterisation as flexsurv
    // mu = location
    // sigma = scale
    // Q = shape
    real lden;
     real w;
    // Constructs the log-density for each observation
    w = ((log(x)-mu))/sigma;
  
    lden = -log(sigma*x)+log(fabs(Q))+pow(Q,-2)*log(pow(Q,-2))+pow(Q,-2)*(Q*w-exp(Q*w))-lgamma(pow(Q,-2));
    
    return lden;
    
  }
  
    real Surv_gengamma(real time, real mu, real sigma, real Q){
    real Sind;
    real qi;
    real w;
    real expnu;
    
    qi = 1/(Q*Q);
    w = (log(time)-mu)/sigma;
    expnu = exp(Q*w)*qi;
    
    Sind =  1- gamma_cdf_man(qi,expnu,10);
    return Sind;
    
  }
  
  
  // Defines the sampling distribution
  real gengamma_lpdf (real t, real d, real d2,  real mu, real sigma, real Q) {
    real log_lik;
    // Issue relates to the inclusion of the survival function;=
    log_lik = d*gen_gamma_lden(t,mu,sigma, Q) + d2*log(Surv_gengamma(t,mu,sigma, Q));
    //log_lik = d*gen_gamma_lden(t,mu,sigma, Q) ;
    return log_lik;
  }
  
  
}

data {
  int<lower=1> n;                      // number of observations
  vector[n] t;               // fully observed times
  vector[n] d;              // censoring indicators
  int<lower=1> H;                         // number of covariates (including intercept)
  matrix[n,H] X;                  // matrix of categorical covariates for the valid cases (0/1 => dummy variables)
  vector[H] mu_beta;                      // vector of means for the covariates
  vector<lower=0>[H] sigma_beta;          // vector of sd for the covariates
  real mu_Q;                              // mean for the parameter Q
  real<lower=0> sigma_Q;                  // sd for the parameter Q
  real<lower=0> a_sigma;                  // first parameter for the scale distribution
  real<lower=0> b_sigma;                  // second parameter for the scale distribution
  
}

transformed data{
vector[n] d2;              // censoring indicators
for (i in 1:n){
  if(d[i] == 1){ 
   d2[i] = 0;
  }else{
   d2[i] =1;
  }
}


}

parameters {
  real Q;                                 // shape of the Generalised Gamma distribution
  real<lower=0> sigma;                    // scale of the Generalised Gamma distribution
  vector[H] beta;                         // coefficients for the covariates
 }

transformed parameters{
  vector[n] linpred;               // rescaled predictor (mu)

  linpred = X*beta;
  
}

model {
  // Prior distributions
  Q ~ normal(mu_Q,sigma_Q);
  sigma ~ gamma(a_sigma,b_sigma);
  beta ~ normal(mu_beta,sigma_beta);
  
  for(i in 1:n){
    t[i] ~ gengamma(d[i],d2[i], linpred[i],sigma,Q);
  }

}

generated quantities {
  real mu;
  mu = beta[1];
}
"


set.seed(123)

data.example = data.frame( times = flexsurv::rgengamma(100, 1.2,.9, .7),
                           status = 1)
#data.example$status[1:5] <-0



stan.data.example <- list()

stan.data.example$t <- data.example$times
stan.data.example$d <- data.example$status
stan.data.example$n <- length(data.example$times)


stan.data.example$H <- 2

stan.data.example$X <- matrix(c(rep(1,length(data.example$times)), 
                                rep(0,length(data.example$times))), ncol = 2)

stan.data.example$mu_beta <- c(0,0)
stan.data.example$sigma_beta <- c(100,100)
stan.data.example$a_sigma <- 0.1
stan.data.example$b_sigma <- 0.1
stan.data.example$mu_Q<- 0
stan.data.example$sigma_Q <- 100


fit <- rstan::stan(model_code = GenGamma.model, model_name = "GenGammaTest", 
                   data = stan.data.example, iter = 2000, chains = 3)

flexsurvreg(Surv(times,status)~1,data =data.example, dist = "gengamma")


stan.funcs <- "
functions {
  real gen_gamma_lden(real x, real mu, real sigma, real Q) {
    // Uses the same parameterisation as flexsurv
    // mu = location
    // sigma = scale
    // Q = shape
    real lden;
     real w;
    // Constructs the log-density for each observation
    w = ((log(x)-mu))/sigma;
  
    lden = -log(sigma*x)+log(fabs(Q))+pow(Q,-2)*log(pow(Q,-2))+pow(Q,-2)*(Q*w-exp(Q*w))-lgamma(pow(Q,-2));
    
    return lden;
    
  }
  
    real Surv_gengamma(real time, real mu, real sigma, real Q){
    real Sind;
    real qi;
    real w;
    real expnu;
    
    qi = 1/(Q*Q);
    w = (log(time)-mu)/sigma;
    expnu = exp(Q*w)*qi;
    
    Sind =  1- gamma_cdf(expnu, qi,1);
    return Sind;
    
  }
  
  
  // Defines the sampling distribution
  real gengamma_lpdf (real t, real d, real d2,  real mu, real sigma, real Q) {
    real log_lik;
    //log_lik = d*gen_gamma_lden(t,mu,sigma, Q) + d2*log(Surv_gengamma(t,mu,sigma, Q));
    log_lik = d*gen_gamma_lden(t,mu,sigma, Q) ;
    return log_lik;
  }
  
  

}"

expose_stan_functions(stanc(model_code = stan.funcs))



time <- 1
mu <- -0.5
sigma <- 0.7
Q <- 0.8

Surv_gengamma(time, mu,sigma,Q)

exp(gen_gamma_lden(time, mu,sigma,Q))
dgengamma(time, mu,sigma,Q)
pgengamma(time, mu,sigma,Q, lower =F)



