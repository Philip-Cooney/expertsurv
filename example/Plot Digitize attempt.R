library("survHE")
library("SHELF")
library("xlsx")
library("abind")
library("mixR")
library("flexsurv")
# surv.inp <- system.file("extdata", "survival.txt", package = "survHE")
# 
# nrisk.inp <- system.file("extdata", "nrisk.txt", package = "survHE")
# 
# digitise(surv.inp, nrisk.inp, km_output = paste0("C:/Users/phili/Desktop/","KMdata.txt"),
#          ipd_output  = paste0("C:/Users/phili/Desktop/","IPDdata.txt"))

pathway <- "C:/Users/phili/OneDrive/PhD/R packages/expertsurv/"


Survival.df <- read.xlsx(paste0(pathway, "ELIANA OS.xlsx"), 1) %>% data.frame()


times.risk <- seq(0, to = 34, by = 2)
n.risk <- c(79,76,73,68,67,62,55,52,47,42,39,26,21,14,9,5,2,0)

upper <- sapply(times.risk, function(x){sum(x > Survival.df$time)})

times.risk<- times.risk[!duplicated(upper)]
n.risk <- n.risk[!duplicated(upper)]
upper <- upper[!duplicated(upper)]


upper <- upper[-1]

lower <- upper
lower <- lower[-length(lower)]
lower <- lower +1

lower <- c(1, lower)

cbind(lower, upper)


write.table(data.frame(Time = Survival.df$time,
                       Survival =Survival.df$surv),
            paste0("C:/Users/phili/Desktop/","survival.txt"),
            row.names=T,
            sep="\t")

write.table(data.frame(Interval = 1:length(lower),
                       time = times.risk[-which(times.risk>tail(Survival.df$time, n= 1))],
                       lower =lower, upper = upper,
                       nrisk =  n.risk[-which(times.risk>tail(Survival.df$time, n= 1))]),
            paste0("C:/Users/phili/Desktop/","nrisk.txt"),
            row.names=FALSE,
            sep="\t")

digitise(paste0("C:/Users/phili/Desktop/","survival.txt"),
         paste0("C:/Users/phili/Desktop/","nrisk.txt"),
         km_output = "C:/Users/phili/Desktop/KMdata.txt",
         ipd_output = "C:/Users/phili/Desktop/IPDdata.txt")

undebug(digitise)

digitized_IPD <- read.table("C:/Users/phili/Desktop/IPDdata.txt")

colnames(digitized_IPD) <- digitized_IPD[1,]

digitized_IPD <- data.frame(digitized_IPD[-1,])

digitized_IPD <- apply(digitized_IPD,2, as.numeric)
km.fit <- survfit(Surv(time, event)~1, data.frame(digitized_IPD))
cbind(km.fit$time,km.fit$surv)
plot(km.fit)
expert.prob.df <- read.xlsx(paste0(pathway, "Expert Surv Probabilities.xlsx"), sheetName = "Sheet1") %>% data.frame()

j <- q <- 1

res.mat <- matrix(nrow = nrow(expert.prob.df)/3, ncol = 3)
for(i in 1:(nrow(expert.prob.df)/3)){
        
        res.mat[j,]    <- expert.prob.df[q:(q+2), 2]   
        
        q <- q+3
        j <- j +1
        
}

expert.prob.df2 <- cbind(res.mat, rep(1:7, times = 4), rep(2:5, each = 7))
colnames(expert.prob.df2) <- c("0.99", "Mode","0.01", "Expert","Time")

times_act <- 4:5

dfs_expert <- list()
plts_pool <- list()
dfs_pool <- list()
lower_bound <- 0
upper_bound <- 1

for(i in 1:length(times_act)){

expert.prob.eval <- expert.prob.df2[expert.prob.df2[,5] == times_act[i],1:3]

expert.prob.eval <- expert.prob.eval[,c(3,2,1)]
fit.eval <- fitdist_mod(t(expert.prob.eval)[-2,],
        probs = c(0.01,  0.99), upper = upper_bound, lower = lower_bound ,mode = t(expert.prob.eval)[2,], 
        expertnames = c(1:7))

#Slight bug with this function (SHELF version)
#SHELF:::plotfit(fit.eval, lp = T)



plts_pool[[i]] <- makePoolPlot(fit  = fit.eval,
                                       xl =lower_bound,
                                       xu =upper_bound,
                                       d = "best",
                                       w = 1,
                                       lwd =1,
                                       xlab = "x",
                                       ylab =expression(f[X](x)),
                                       legend_full = TRUE,
                                       ql = NULL,
                                       qu = NULL,
                                       nx = 200,
                                       addquantile = FALSE,
                                       fs = 12,
                                       expertnames = NULL)

dfs_pool[[i]] <-  plts_pool[[i]][["data"]]



best_fit_index  <- apply(fit.eval$ssq[,c("normal","t","gamma", "lognormal", "beta")], 1, which.min)
best_fit <- names(fit.eval$ssq[,c("normal","t","gamma", "lognormal", "beta")])[best_fit_index]

best_fit_loc <- sapply(best_fit, function(x){which(x  == names(fit.eval$ssq))})
fit.eval.dist  <- fit.eval[best_fit_loc]

pool.df_output <- matrix(nrow = length(best_fit_loc),ncol = 3)
colnames(pool.df_output) <- c("param1", "param2", "param3")

for(j in 1:length(best_fit_loc)){
        pool.df_output[j,1:length(fit.eval.dist[[j]][j,])] <-  as.numeric(as.vector(fit.eval.dist[[j]][j,]))
}
dfs_expert[[i]] <- data.frame(dist = best_fit, wi = 1, pool.df_output)
}

dist_fit <- c()
pool.eval.stan <- list()
pool_type_eval <- "linear pool"

for(i in 1:length(dfs_expert)){
        
        pool.eval.stan[[i]]  <- makePoolPlot.Data(pool.df = dfs_expert[[i]], 
                                                  pool_type = pool_type_eval,
                                                  add_hist = F, plt_other_dists = F, max_mix_eval = 5)
        
        
        pool.df <- subset(dfs_pool[[i]], ftype == pool_type_eval) %>% rename(Dist = ftype)
        pool.eval.stan[[i]]$plot.fit.mixnorm <- pool.eval.stan[[i]]$plot.fit.mixnorm+
               
                geom_line(data =pool.df,
                          aes(x = x, y = fx, colour = Dist),  lty = 2)+
                xlim(c(0,1))+
                ylim(c(0, max(pool.df$fx)*1.1))
              
        dist_fit <- c(dist_fit, names(pool.eval.stan[[i]])[2])
        if(i == 1){
                param_expert <-  pool.eval.stan[[i]][[2]]
                
        }else{
                param_expert <-  abind(param_expert,pool.eval.stan[[i]][[2]],along = 1)
                
        }
        
        rownames(param_expert) <- dist_fit
        
}



# Put This into formula



library("survHE")

data(bc)

data.example = data.frame( times = flexsurv::rgengamma(100, 1.2,.9, .7),
                           status = 1)
data.example$status[1:5] <-0
# Fits a parametric model
m.gengamma <- survHE::fit.models(formula=Surv(times,status)~1,data=data.example,
                        distr="gengamma",method="hmc")

m.gengamma$models$`Gen. Gamma`

m.weibull <- survHE::fit.models(formula=Surv(times,status)~1,data=data.example,
                                 distr="weibull",method="hmc")


num_param <- c()

stan.data.weibull <- m.weibull[["misc"]][["data.stan"]][[1]]


stan.data <- m.gengamma[["misc"]][["data.stan"]][[1]]


set.seed(123)

data.example = data.frame( times = flexsurv::rgengamma(100, 1.2,.9, .7),
                           status = 1)
#data.example$status[1:5] <-0



stan.data.example <- list()

stan.data.example$t <- data.example$times
stan.data.example$d <- data.example$status
stan.data.example$n <- numeric(length(data.example$times))


stan.data.example$H <- 2

stan.data.example$X <- matrix(c(rep(1,length(data.example$times)), 
                                 rep(0,length(data.example$times))), ncol = 2)

stan.data.example$mu_beta <- c(0,0)
stan.data.example$sigma_beta <- c(100,100)
stan.data.example$a_sigma <- numeric(0.1)
stan.data.example$b_sigma <- numeric(0.1)
stan.data.example$mu_Q<- 0
stan.data.example$sigma_Q <- 100


# stan.data$X <- stan.data.weibull$X
# 
# stan.data$t <- stan.data.weibull$t
# stan.data$d <- stan.data.weibull$d
# stan.data$n <- stan.data.weibull$n
# 


for(i in 1:dim(param_expert)[1]){
        num_param[i] <-   sum(param_expert[i,,1] != -999.2)
}

lk_up_dist <- c("mixture","norm", "t", "gamma", "lnorm","beta")

data_dist_ind<- as.numeric(sapply(rownames(param_expert), function(x){which(x== lk_up_dist)}))

stan.data$num_param <- c(2,2)#num_param;
stan.data$data_dist_ind <- c(2,2)#data_dist_ind;
stan.data$max_param <- max(num_param)
stan.data$n_time_expert <- length(times_act)
stan.data$param_expert <-param_expert[,1:max(num_param),];
times_expert <- c(4,5)


stan.data
if(!is.null(times_expert)){
  stan.data$n_time_expert <- length(times_expert)
  stan.data$time_expert <- times_expert
}else{
  stan.data$n_time_expert <- 1
  stan.data$time_expert <- as.numeric(0) # dummy data 
  
}
distr3 = "gga"
id_St = 0
opinion_type = "survival"
if(opinion_type == "survival"){
  stan.data$St_indic<- 1
  stan.data$id_comp <-  numeric(0)
  stan.data$id_trt <-  numeric(0)
  
  if(distr3 %in% c("gam", "gga", "gef")){
    
    stan.data$id_St <-  as.array(min(which(stan.data[["X_obs"]][,2] == id_St)))
    
  }else{
    
    stan.data$id_St <-  as.array(min(which(stan.data[["X"]][,2] == id_St)))
    
  }
  
}else{
  
  stan.data$St_indic <- 0
  
  stan.data$id_St <-  numeric(0)
  
}
stan.data$a0 <- rep(1, nrow(stan.data$X))

stan.data$a0_obs <- rep(1, nrow(stan.data$X_obs))
stan.data$a0_cens <- rep(1, nrow(stan.data$X_cens))


stanmodelcode_valid <- "

functions {
   real log_mixnorm_dens(real[ ] mu, real[ ] sd, real[ ] prob,  real x) {
    
    // Evaluates the log-density based on a mixture normal
    // x refers to the Survival or expected difference in survival
    
    //real log_dens;
    int n;
    vector[num_elements(mu)] log_dens_i;
    
    n = num_elements(mu);
    for(i in 1:n){
        log_dens_i[i]= normal_lpdf(x| mu[i], sd[i]) + log(prob[i]);
    }

     return log(sum(exp(log_dens_i)));
    }
  
  real log_density_dist(int data_dist_ind, real[ ] param1, 
                        real[ ] param2, real[ ] param3, real x){
     
     real log_dens;                    
     if(data_dist_ind == 1){
      log_dens = log_mixnorm_dens(param1, param2, param3,x);
     }else if(data_dist_ind == 2){
      
      log_dens = normal_lpdf(x|param1[1], param2[1]);
      
      }else if(data_dist_ind == 3){
      
      log_dens = student_t_lpdf(x|param3[1],param1[1], param2[1]);
      
      }else if(data_dist_ind == 4){
      
      log_dens = gamma_lpdf(x|param1[1], param2[1]);
      
      }else if(data_dist_ind == 5){
      
      log_dens = lognormal_lpdf(x|param1[1], param2[1]);
      
      }else if(data_dist_ind == 6){
      
      log_dens = beta_lpdf(x|param1[1], param2[1]);
      
      } 
     
     return(log_dens);
  
  }
  

	
}

data {

  int n_time_expert;
  int max_param;
  real param_expert[n_time_expert,max_param,3];
  int num_param[n_time_expert];
  int data_dist_ind[n_time_expert];
  
}

parameters {
 vector<lower=0,upper=1>[n_time_expert] St_expert;
}


model {
  
  St_expert ~uniform(0,1);
	
	
	for (i in 1:n_time_expert){
	
	target += log_density_dist(data_dist_ind[i], 
	                           param_expert[i,1:num_param[i],1],
	                           param_expert[i,1:num_param[i],2],
	                           param_expert[i,1:num_param[i],3],
	                           St_expert[i]);
	 }

}


"

## gengamm and GEN F don't work with a0 or defining the surival function
stanmodelcode_valid2 <- "// Generalised Gamma model for survival analysis - uses the parameterisation of flexsurv 
//https://dev.to/martinmodrak/optional-parametersdata-in-stan-4o33

// First defines the Generalised Gamma model for fully observed and censored cases

functions {
  real gen_gamma_lden(real x, real mu, real sigma, real Q) {
    // Uses the same parameterisation as flexsurv
    // mu = location
    // sigma = scale
    // Q = shape
    real lprob;
   
    real w;
    // Constructs the log-density for each observation
    w = ((log(x)-mu))/sigma;
  
    lprob = -log(sigma*x)+log(fabs(Q))+pow(Q,-2)*log(pow(Q,-2))+pow(Q,-2)*(Q*w-exp(Q*w))-lgamma(pow(Q,-2));
    
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
  
  
  // Defines the sampling distribution
  real surv_gengamma_lpdf (real t, real d, real d2,  real mu, real sigma, real Q, real a0) {
    real log_lik;
    real lprob;
    log_lik = d*gen_gamma_lden(t,mu,sigma, Q) + d2*Sind(t,mu,sigma, Q);
    lprob = log_lik*a0;
    return lprob;
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
  
  real log_mixnorm_dens(real[ ] mu, real[ ] sd, real[ ] prob,  real x) {
    
    // Evaluates the log-density based on a mixture normal
    // x refers to the Survival or expected difference in survival
    
    //real log_dens;
    int n;
    vector[num_elements(mu)] log_dens_i;
    
    n = num_elements(mu);
    for(i in 1:n){
      log_dens_i[i]= normal_lpdf(x| mu[i], sd[i]) + log(prob[i]);
    }
    
    return log(sum(exp(log_dens_i)));
  }
  
  real log_density_dist(int data_dist_ind, real[ ] param1, 
                        real[ ] param2, real[ ] param3, real x){
    
    // Evaluates the log density for a range of distributions
    
    real log_dens;                    
    if(data_dist_ind == 1){
      log_dens = log_mixnorm_dens(param1, param2, param3,x);
    }else if(data_dist_ind == 2){
      
      log_dens = normal_lpdf(x|param1[1], param2[1]);
      
    }else if(data_dist_ind == 3){
      
      log_dens = student_t_lpdf(x|param3[1],param1[1], param2[1]);
      
    }else if(data_dist_ind == 4){
      
      log_dens = gamma_lpdf(x|param1[1], param2[1]);
      
    }else if(data_dist_ind == 5){
      
      log_dens = lognormal_lpdf(x|param1[1], param2[1]);
      
    }else if(data_dist_ind == 6){
      
      log_dens = beta_lpdf(x|param1[1], param2[1]);
      
    } 
    
    return(log_dens);
    
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
  
  int n_time_expert;
  int<lower = 0, upper = 1> St_indic; // 1 Expert opinion on survival @ timepoint ; 0 Expert opinion on survival difference
  
  int id_St[St_indic ? 1 : 0];
  int id_trt[St_indic ? 0 : 1];
  int id_comp[St_indic ? 0 : 1];
  
  vector[n] a0;              
  
  int max_param;
  real param_expert[n_time_expert,max_param,3];
  vector[St_indic ? n_time_expert : 0] time_expert;
  int num_param[n_time_expert];
  int data_dist_ind[n_time_expert];
  
  
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
  //vector<lower=1>[n_cens] cens;           // censoring variable (latent)
  
}

transformed parameters{
  vector[n] linpred;               // rescaled predictor (mu) for the observed cases
   vector[n_time_expert] St_expert;
  vector[n_time_expert] St_expert_diff;
  vector[n_time_expert] St_expert_final;
  
  linpred = X*beta;
  
  
  for (i in 1:n_time_expert){
  if(St_indic == 1){
      St_expert[i] = Sind(time_expert[i],linpred[id_St[1]],sigma,Q); 
  }else{
      St_expert[i] = Surv_diff(sigma,Q,linpred[id_trt[1]],linpred[id_comp[1]]);
  
    }
  }  
  
}

model {
  // Prior distributions
  Q ~ normal(mu_Q,sigma_Q);
  sigma ~ gamma(a_sigma,b_sigma);
  beta ~ normal(mu_beta,sigma_beta);
  
  for(i in 1:n){
    t[i] ~ surv_gengamma(d[i],d2[i], linpred[i],sigma,Q, a0[i]);
  }

 
  for (i in 1:n_time_expert){
      target += log_density_dist(data_dist_ind[i], 
                                 param_expert[i,1:num_param[i],1],
                                 param_expert[i,1:num_param[i],2],
                                 param_expert[i,1:num_param[i],3],
                                 St_expert[i]);
    } 
}

generated quantities {
  real mu;
  mu = beta[1];
}
"



## gengamm and GEN F don't work with a0 or defining the surival function


stanmodelcode_valid2<- "// Generalised Gamma model for survival analysis - uses the parameterisation of flexsurv 

// First defines the Generalised Gamma model for fully observed and censored cases

functions {
  real gen_gamma_lpdf(real x, real mu, real sigma, real Q, real a0) {
    // Uses the same parameterisation as flexsurv
    // mu = location
    // sigma = scale
    // Q = shape
    real prob;
    real lprob;
   real w;
    // Constructs the log-density for each observation
    w = ((log(x)-mu))/sigma;
   
      prob = -log(sigma*x)+log(fabs(Q))+pow(Q,-2)*log(pow(Q,-2))+pow(Q,-2)*(Q*w-exp(Q*w))-lgamma(pow(Q,-2));
   
    // And the total log-density (as a sum of the individual terms)
    lprob = prob*a0;
    return lprob;
  }
  
  real gen_gamma_cens_lpdf(real x, real mu, real sigma, real Q, real u, real a0) {
    // Rescales the distribution accounting for right censoring
    real prob;
    real lprob;
    real w;
    real tr;
    // Constructs the log-density for each observation
    tr = x * u;
    w = ((log(tr)-mu))/sigma;
    prob = log(u)-log(sigma*tr)+log(fabs(Q))+pow(Q,-2)*log(pow(Q,-2))+pow(Q,-2)*(Q*w-exp(Q*w))-lgamma(pow(Q,-2));
   
    // And the total log-density (as a sum of the individual terms)
    lprob = prob*a0;
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
  
    real log_mixnorm_dens(real[ ] mu, real[ ] sd, real[ ] prob,  real x) {
    
    // Evaluates the log-density based on a mixture normal
    // x refers to the Survival or expected difference in survival
    
    //real log_dens;
    int n;
    vector[num_elements(mu)] log_dens_i;
    
    n = num_elements(mu);
    for(i in 1:n){
      log_dens_i[i]= normal_lpdf(x| mu[i], sd[i]) + log(prob[i]);
    }
    
    return log(sum(exp(log_dens_i)));
  }
  
  real log_density_dist(int data_dist_ind, real[ ] param1, 
                        real[ ] param2, real[ ] param3, real x){
    
    // Evaluates the log density for a range of distributions
    
    real log_dens;                    
    if(data_dist_ind == 1){
      log_dens = log_mixnorm_dens(param1, param2, param3,x);
    }else if(data_dist_ind == 2){
      
      log_dens = normal_lpdf(x|param1[1], param2[1]);
      
    }else if(data_dist_ind == 3){
      
      log_dens = student_t_lpdf(x|param3[1],param1[1], param2[1]);
      
    }else if(data_dist_ind == 4){
      
      log_dens = gamma_lpdf(x|param1[1], param2[1]);
      
    }else if(data_dist_ind == 5){
      
      log_dens = lognormal_lpdf(x|param1[1], param2[1]);
      
    }else if(data_dist_ind == 6){
      
      log_dens = beta_lpdf(x|param1[1], param2[1]);
      
    } 
    
    return(log_dens);
    
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
  int<lower = 0, upper = 1> St_indic; // 1 Expert opinion on survival @ timepoint ; 0 Expert opinion on survival difference
  
  int id_St[St_indic ? 1 : 0];
  int id_trt[St_indic ? 0 : 1];
  int id_comp[St_indic ? 0 : 1];
  
 
  int max_param;
  real param_expert[n_time_expert,max_param,3];
  vector[St_indic ? n_time_expert : 0] time_expert;
  int num_param[n_time_expert];
  int data_dist_ind[n_time_expert];
  
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
 
  linpred_cens = X_cens*beta;
  linpred_obs = X_obs*beta;
  
  for (i in 1:n_time_expert){
  if(St_indic == 1){
      St_expert[i] = Sind(time_expert[i],linpred_obs[id_St[1]],sigma,Q); 
  }else{
      St_expert[i] = Surv_diff(sigma,Q,linpred_obs[id_trt[1]],linpred_obs[id_comp[1]]);
  
    }
  }  
}

model {
  // Prior distributions
  Q ~ normal(mu_Q,sigma_Q);
  sigma ~ gamma(a_sigma,b_sigma);
  beta ~ normal(mu_beta,sigma_beta);
  // Data model
  if(n_cens>0){
  for(i in 1:n_cens)
    cens[i] ~ gen_gamma_cens(linpred_cens[i],sigma,Q,d[i], a0_cens[i]);
  }
  
  for(i in 1:n_obs){
    t[i] ~ gen_gamma(linpred_obs[i],sigma,Q,a0_obs[i]);
  }

  
for (i in 1:n_time_expert){
      target += log_density_dist(data_dist_ind[i], 
                                 param_expert[i,1:num_param[i],1],
                                 param_expert[i,1:num_param[i],2],
                                 param_expert[i,1:num_param[i],3],
                                 St_expert[i]);
    } 

}

generated quantities {
  real mu;
  mu = beta[1];
}"



stanmodelcode_valid2<- "// Generalised Gamma model for survival analysis - uses the parameterisation of flexsurv 

// First defines the Generalised Gamma model for fully observed and censored cases

functions {
real gen_gamma_lpdf(vector x, vector mu, real sigma, real Q) {
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
    lprob = sum((prob));
    return lprob;
  }
  
  real gen_gamma_cens_lpdf(vector x, vector mu, real sigma, real Q, vector u) {
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
    lprob = sum((prob));
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
  
    real log_mixnorm_dens(real[ ] mu, real[ ] sd, real[ ] prob,  real x) {
    
    // Evaluates the log-density based on a mixture normal
    // x refers to the Survival or expected difference in survival
    
    //real log_dens;
    int n;
    vector[num_elements(mu)] log_dens_i;
    
    n = num_elements(mu);
    for(i in 1:n){
      log_dens_i[i]= normal_lpdf(x| mu[i], sd[i]) + log(prob[i]);
    }
    
    return log(sum(exp(log_dens_i)));
  }
  
  real log_density_dist(int data_dist_ind, real[ ] param1, 
                        real[ ] param2, real[ ] param3, real x){
    
    // Evaluates the log density for a range of distributions
    
    real log_dens;                    
    if(data_dist_ind == 1){
      log_dens = log_mixnorm_dens(param1, param2, param3,x);
    }else if(data_dist_ind == 2){
      
      log_dens = normal_lpdf(x|param1[1], param2[1]);
      
    }else if(data_dist_ind == 3){
      
      log_dens = student_t_lpdf(x|param3[1],param1[1], param2[1]);
      
    }else if(data_dist_ind == 4){
      
      log_dens = gamma_lpdf(x|param1[1], param2[1]);
      
    }else if(data_dist_ind == 5){
      
      log_dens = lognormal_lpdf(x|param1[1], param2[1]);
      
    }else if(data_dist_ind == 6){
      
      log_dens = beta_lpdf(x|param1[1], param2[1]);
      
    } 
    
    return(log_dens);
    
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
  int<lower = 0, upper = 1> St_indic; // 1 Expert opinion on survival @ timepoint ; 0 Expert opinion on survival difference
  
  int id_St[St_indic ? 1 : 0];
  int id_trt[St_indic ? 0 : 1];
  int id_comp[St_indic ? 0 : 1];
  
 
  int max_param;
  real param_expert[n_time_expert,max_param,3];
  vector[St_indic ? n_time_expert : 0] time_expert;
  int num_param[n_time_expert];
  int data_dist_ind[n_time_expert];
  
  //vector<lower=0>[n_obs] a0_obs;              
  //vector<lower=0>[n_cens]a0_cens ;              
  
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
 
  linpred_cens = X_cens*beta;
  linpred_obs = X_obs*beta;
  
  for (i in 1:n_time_expert){
  if(St_indic == 1){
      //St_expert[i] = Sind(time_expert[i],linpred_obs[id_St[1]],sigma,Q); 
 
      St_expert[i] = 0.5;
  }else{
   // St_expert[i] = Surv_diff(sigma,Q,linpred_obs[id_trt[1]],linpred_obs[id_comp[1]]);
  
  }
  }  
}

model {
  // Prior distributions
  Q ~ normal(mu_Q,sigma_Q);
  sigma ~ gamma(a_sigma,b_sigma);
  beta ~ normal(mu_beta,sigma_beta);
  // Data model
  if(n_cens>0){
  cens ~ gen_gamma_cens(linpred_cens,sigma,Q,d);
  }
  
  t ~ gen_gamma(linpred_obs,sigma,Q);
 

  
 for (i in 1:n_time_expert){
       target += log_density_dist(data_dist_ind[i], 
                                  param_expert[i,1:num_param[i],1],
                                  param_expert[i,1:num_param[i],2],
                                  param_expert[i,1:num_param[i],3],
                                  St_expert[i]);
     } 

}

generated quantities {
  real mu;
  mu = beta[1];
}"
##https://onlinelibrary.wiley.com/doi/epdf/10.1111/j.1467-9574.1991.tb01312.x

# 
# fit <- rstan::stan(model_code = stanmodelcode_valid, model_name = "DensityEval", 
#             data = stan.data, iter = 50000, chains = 3,control = list(adapt_delta = 0.99))

fit <- rstan::stan(model_code = stanmodelcode_valid2, model_name = "GenGammaTest", 
            data = stan.data, iter = 2000, chains = 3)


stan.data$mu_St <- c(0.5, 0.4)
stan.data$sigma_St <- c(0.1, 0.1)

stan.data.test <- stan.data


test <- '// Generalised Gamma model for survival analysis - uses the parameterisation of flexsurv 

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
}'

stan.data.test$id_St <- 50


stan.data.test$id_trt <- 5
stan.data.test$id_comp <- 1
stan.data.test$mu_diff <- 1
stan.data.test$sigma_diff <- 1
fit.test <- rstan::stan(model_code = test, model_name = "GenGammaTest", 
                   data = stan.data.test, iter = 2000, chains = 3)


#plot(fit, show_density = TRUE)

matrix_of_draws <- as.matrix(fit)

density.1 <- density(matrix_of_draws[,1])
density.df1 <- data.frame(x = density.1$x, y = density.1$y) 

density.2 <- density(matrix_of_draws[,2])
density.df2 <- data.frame(x = density.2$x, y = density.2$y) 

dens_eval <- seq(0,1, by = 0.01)
mix.norm.dens1 <- rowSums(apply(param_expert[1,1:num_param[1],],1, function(x){dnorm(dens_eval, x[1], x[2])*x[3]}))

# k_norm <- sfsmisc::integrate.xy(x = dens_eval,fx = mix.norm.dens1)
# k_norm <- sfsmisc::integrate.xy(x = density.df1$x,fx = density.df1$y)

#mix.norm.dens.df1 <- data.frame(x = dens_eval, y = mix.norm.dens1)

pool.eval.stan[[1]]$plot.fit.mixnorm+
        geom_line(data = density.df1, aes(x, y), colour = "green") #+
        #geom_line(data = mix.norm.dens.df1, aes(x, y), colour = "orange")

#Validated Function is correct
#exp(log_density_dist(1,param_expert[1,1:num_param[1],1],param_expert[1,1:num_param[1],2],param_expert[1,1:num_param[1],3], 0.2))

pool.eval.stan[[2]]$plot.fit.mixnorm+
        geom_line(data = density.df2, aes(x, y), colour = "green")


data.stan.example$time_expert <- times_act



v <- matrix(c(30, 40, 50, 20, 25, 35), 3, 2)
p <- c(0.25, 0.5, 0.75)
SHELF:::gamma.error


