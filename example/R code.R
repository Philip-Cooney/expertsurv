library("devtools")
library("rstan")

library("SHELF")
library("mixtools")
library("mixR")
library("pkgcond")
library("abind")
# library("rstantools")
# library("rstan")
# library("StanHeaders")
# library("survival")
# library("flexsurv")
# library("expertsurv")
load_all()
#install

#library("survHE")

#expertsurv:::stanmodels

#https://cran.r-project.org/web/packages/rstantools/vignettes/minimal-rstan-package.html
#https://discourse.mc-stan.org/t/unable-to-build-package-with-rstan-unable-to-load-module-error/14306/3
# example(source) # defines the sourceDir() function
# roxygen2::roxygenize(load_code = sourceDir)
# 
# 
# example(source) # defines the sourceDir() function
# roxygen2::roxygenize()
# remove.packages(expert)
# devtools::install(quick = T)


#use_rstan(pkgdir = "C:/Users/phili/OneDrive/PhD/R packages/expertsurv/")
## Master function



mods <- c("exp", "weibull", "gamma", "lnorm", "llogis", "gengamma")

mods <- c("gengamma")

#formula <- Surv(recyrs, censrec) ~ as.factor(group)

formula <- Surv(recyrs, censrec) ~ 1
mu_St <- 0.30
sigma_St <- 0.02
time_expert <- 5

expert.stan.weibull<- expert_surv2(formula, data = bc,
                                         distr =mods, method = "hmc",
                                         mu_St = mu_St,
                                         sigma_St = sigma_St,
                                         time_expert = time_expert,
                                         id_St = 0, 
                                         opinion_type = "survival")#
plot(expert.stan.weibull, add.km = T)

plot(expert.stan.weibull$models$Exponential, pars = "St_expert")

expert_df <- data.frame(expert = c(1,1,2,2),
              times_expert = c(1,5,1,5),
              dist = c("norm", "norm", "norm", "norm"),
              param1 = c(0.5,0.25, 0.55, .35),
              param2 = c(0.05, 0.05, 0.05, 0.05),
              weights = c(0.3, 1, 0.4, 1))



unique_times <- c(1, 5)
max.quant.vals <- 5
max.expert <- 2
max.timepoints  <- length(unique_times)

expert_quant_array <- array(dim = c(5,2,2))
expert_prob_matrix <- matrix(NA, nrow = max.quant.vals, ncol = max.timepoints)
expert_weig_matrix <- matrix(0.5, nrow = max.expert, ncol = max.timepoints)

v1 <- matrix(c(.30, .40, .50, .20, .25, .35), 3, 2)
p1 <- c(0.25, 0.5, 0.75)
v2 <- matrix(c(.1, .20, .30, .40, 0.5, 0.1, 0.25, 0.29, 0.45, 0.55), 5, 2)
p2 <- c(0.15, 0.27, 0.32, 0.56, 0.79)

for(i in 1:max.timepoints){
  current_v_mat <- get(paste0("v",i))
  expert_quant_array[1:dim(current_v_mat)[1],1:dim(current_v_mat)[2],i] <- current_v_mat 
 
  current_p_vec <- get(paste0("p",i))
  expert_prob_matrix[1:length(current_p_vec),i] <- current_p_vec 
  
}


quantile_list <- list()
quantile_list$v_array <- expert_quant_array
quantile_list$p_mat <- expert_prob_matrix
quantile_list$weights <- NULL
quantile_list$times <- unique_times

expert_final <- expert_pooling(expert_quant_list = quantile_list)

dist_fit <- c()
pool.eval.stan <- list()
for(i in 1:length(expert_final$dfs_expert)){
  
  pool.eval.stan[[i]]  <- makePoolPlot.Data(pool.df = expert_final$dfs_expert[[i]], 
                                            pool_type = "linear pool",
                                            add_hist = F)
  
  pool.eval.stan[[i]]$plot.fit.mixnorm <- pool.eval.stan[[i]]$plot.fit.mixnorm+
    geom_line(data = subset(expert_final$dfs_pool[[i]], ftype == "linear pool"),
              aes(x = x, y = fx, colour = ftype),  lty = 2)
  
  dist_fit <- c(dist_fit, names(pool.eval.stan[[i]])[2])
  if(i == 1){
    param_expert <-  pool.eval.stan[[i]][[2]]
    
  }else{
    param_expert <-  abind(param_expert,pool.eval.stan[[i]][[2]],along = 1)
    
  }
  
}


data.stan.example <- expert.stan.weibull[["misc"]][["data.stan"]][[1]]
data.stan.example$param_expert <-param_expert;
num_param <- c()

for(i in 1:dim(param_expert)[1]){
num_param[i] <-   sum(param_expert[i,,][1,] != -999.2)
}

lk_up_dist <- c("mixture","norm", "t", "gamma", "lnorm","beta")

data_dist_ind<- as.numeric(sapply(dist_fit, function(x){which(x== lk_up_dist)}))

data.stan.example$num_param <- num_param;

data.stan.example$data_dist_ind <- data_dist_ind;
data.stan.example$max_param <- max(num_param)
data.stan.example$n_time_expert <- 2
data.stan.example$time_expert <- c(3,4)

stanmodelcode <- "

functions {
  // Defines the log hazard
  vector log_h (vector t, vector rate) {
    vector[num_elements(t)] log_h;
    log_h = log(rate);
    return log_h;
  }
  
  // Defines the log survival
  vector log_S (vector t, vector rate) {
    vector[num_elements(t)] log_S;
    log_S = -rate .* t;
    return log_S;
  }
  
  // Defines the log survival indvidual
  real log_Sind (real t, real rate) {
	real log_Sind;
      log_Sind = -rate .* t;
    return log_Sind;
  }  
  
  // Defines difference in expected survival
  real Surv_diff ( real rate_trt, real rate_comp) {
	real Surv_diff;
      Surv_diff = 1/rate_trt - 1/rate_comp;
    return Surv_diff;
  }
    
  // Defines the sampling distribution
  real surv_exponential_lpdf(vector t, vector d, vector rate, vector a0) {
    vector[num_elements(t)] log_lik;
    real prob;
    log_lik = d .* log_h(t,rate) + log_S(t,rate);
    prob = dot_product(log_lik, a0);
    return prob;
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
  int n;                  // number of observations
  vector[n] t;            // observed times
  vector[n] d;            // censoring indicator (1=observed, 0=censored)
  int H;                  // number of covariates
  matrix[n,H] X;          // matrix of covariates (with n rows and H columns)
  vector[H] mu_beta;	    // mean of the covariates coefficients
  vector<lower=0> [H] sigma_beta;   // sd of the covariates coefficients
  
  vector[n] a0; //Power prior for the observations
  
  int n_time_expert;
  // vector[n_time_expert] mu_St;
  //vector[n_time_expert] sigma_St;
  vector[n_time_expert] time_expert;
  int id_St;
  
  int St_indic; // 1 Expert opinion on survival @ timepoint ; 2 Expert opinion on survival difference
  int id_trt;
  int id_comp;
  real mu_diff;
  real<lower = 0> sigma_diff;
  
  int max_param;
  real param_expert[n_time_expert,max_param,3];
  int num_param[n_time_expert];
  int data_dist_ind[n_time_expert];
  
}

parameters {
  vector[H] beta;         // Coefficients in the linear predictor (including intercept)
}

transformed parameters {
  vector[n] linpred;
  vector[n] mu;
  vector<lower=0,upper=1>[n_time_expert] St_expert;
  real Surv_diff_val;
  
  linpred = X*beta;
  for (i in 1:n) {
    mu[i] = exp(linpred[i]);  // Rate parameter
  }
  
  if(St_indic == 1){
	for (i in 1:n_time_expert){
		St_expert[i] = exp(log_Sind(time_expert[i],mu[id_St]));
	 }
	}else{    
	Surv_diff_val = Surv_diff(mu[id_trt],mu[id_comp]);
  	}
  
}

model {
  beta ~ normal(mu_beta,sigma_beta);
  t ~ surv_exponential(d,mu, a0);
  
 if (St_indic ==1) {
	for (i in 1:n_time_expert){
	target += log_density_dist(data_dist_ind[i], 
	                           param_expert[i,1:num_param[i],1],
	                           param_expert[i,1:num_param[i],2],
	                           param_expert[i,1:num_param[i],3],
	                           St_expert[i]);
	 }
  }else{
    target += normal_lpdf(Surv_diff_val | mu_diff, sigma_diff);
	}
  
  
}

generated quantities {
  real rate;                // rate parameter
  rate = exp(beta[1]);
} 
"


fit <- stan(model_code = stanmodelcode, model_name = "Exponential2", 
            data = data.stan.example, iter = 2000, chains = 3)

plot(expert.stan.weibull, add.km = T)+
  geom_line(data = df_surv, aes(x = time, y = surv))


post_exp <- exp(summary(fit, pars = "beta", probs = c(0.50))$summary[1,"50%"])

df_surv <- data.frame(time = seq(0, 6, by= 0.1), surv = exp(-seq(0, 6, by= 0.1)*post_exp) )


summary(fit)

plot(fit, pars  ="St_expert")
summary(fit)

# fit <- stan(file = "C:/Users/phili/Desktop/exponential.stan", model_name = "Exponential2", 
#             data = data.stan.example, iter = 2012, chains = 3)
# mod <- stan_model(model_code = "C:/Users/phili/Desktop/exponential.stan", verbose = TRUE)
# mod <- stan_model(model_code = stanmodelcode, verbose = TRUE)





rstan::sampling(mod, data = data.stan.example)

data_dist <- c("mixture","normal", "t", "gamma", "lognormal", "beta")
data_dist_ind<- rep(0, length(data_dist))

data_dist_ind[which(best_fit == data.dist)] <- 1


data_dist_ind <- c(0,1,0,0,0,0)
param1 <- c(1)
prob <- c(NA)
param2 <- c(2)
x <- 0.2

log_density_dist(data_dist_ind, param1, 
                 param2, prob, x)

log_density_dist(t(data_dist_ind), t(param1), 
                 t(param2), t(prob), x)

log(dnorm(x,param1,param2))





polydens <- "functions {
  
  real polydens_lpdf(vector coef, vector power,  real x) {
    
    // Constructs the log-density based on a custom spline
    // x refers to the Survival or expected difference in survival
    
    real dens;
    real dens_bound;
    int n;
    vector[num_elements(coef)] x_pow;
    
    n = num_elements(coef);
    for(i in 1:n){
        x_pow[i]= pow(x, power[i]);
    }

    dens = dot_product(coef,x_pow);
    
    if(dens <0.0){
    dens_bound = 0.0;
    }else{
    dens_bound = dens;
    }
    
    return log(dens_bound);
  }
	
}"


library(rstan)

mixnorm_dens <- "functions {
  
    real log_mixnorm_dens(row_vector mu, row_vector sd, row_vector prob,  real x) {
    
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
  
  real log_density_dist(vector data_dist_ind, row_vector param1, 
                        row_vector param2, row_vector param3, real x){
     
     real log_dens;                    
     if(data_dist_ind[1] == 1){
      log_dens = log_mixnorm_dens(param1, param2, param3,x);
     }else if(data_dist_ind[2] == 1){
      
      log_dens = normal_lpdf(x|param1[1], param2[1]);
      
      }else if(data_dist_ind[3] == 1){
      
      log_dens = student_t_lpdf(x|param3[1],param1[1], param2[1]);
      
      }else if(data_dist_ind[4] == 1){
      
      log_dens = gamma_lpdf(x|param1[1], param2[1]);
      
      }else if(data_dist_ind[5] == 1){
      
      log_dens = lognormal_lpdf(x|param1[1], param2[1]);
      
      }else if(data_dist_ind[6] == 1){
      
      log_dens = beta_lpdf(x|param1[1], param2[1]);
      
      } 
     
     return(log_dens);
  
  }
  
	
}"

expose_stan_functions(stanc(model_code = mixnorm_dens))


library(rstan)
expose_stan_functions(stanc(model_code = polydens))
exp(polydens_lpdf(mod$coefficients,seq(0,length(mod$coefficients)-1), 0.5 ))



data.stan.example$sigma_beta <- c(1,1,1)
data.stan.example$mu_beta <- c(1,1,1)
data.stan.example$a_alpha <- c(1);
data.stan.example$b_alpha <- c(1);



weibull_expert_stan <- rstan::sampling(expertsurv:::stanmodels$WeibullPH_expert, data = data.stan.example
)

devtools::load_all()
View(m1)
print(m1, mod = 1)

plot(m1, add.km = TRUE)+
  scale_x_continuous(expand = c(0, 0), limits = c(0,NA)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA))

model.fit.plot(m1)

#Change to WAIC

psa <- make.surv(fit = expert.stan.weibull, nsim = 1000, t = seq(.1, 63))

psa.plot(psa, xlab = "Extrapolated time",
         ylab = "Estimation of the survival curves",
         alpha = 0.2, col = c("dark grey", "black", "green"),
         main = "PSA to survival curves", cex.txt = .95)






library("survHE")
library("flexsurv")
library("rstan")
library("dplyr")
pathway <-"C:/Users/phili/OneDrive/PhD/R codes/Stan Expert Codes/"
source(paste0(pathway,"SurvHE functions.R"))
rstan_options(auto_write = TRUE)

set.seed(123)
n <- 300
times <- rgengamma(n, mu = 1, sigma = 1, Q = 2)
df <- data.frame(time_event = times)


censoring_time <-1
df <- mutate(df, time = ifelse(time_event > censoring_time, censoring_time,time_event))
df <- mutate(df, status = ifelse(time_event > censoring_time, 0, 1)) %>% mutate(id = 1)

times2 <- rgengamma(n, mu = 1.5, sigma = 1, Q = 2)
df2 <- data.frame(time_event = times2)
df2 <- mutate(df2, time = ifelse(time_event > censoring_time, censoring_time,time_event))
df2 <- mutate(df2, status = ifelse(time_event > censoring_time, 0, 1)) %>% mutate(id = 2)

df_comb <- bind_rows(df, df2)
plot(survfit(Surv(time,status)~id,df_comb))

mods <- c("rps")
#formula <- Surv(time, status) ~as.factor(id)

formula <- Surv(time, status) ~1


mu_St <- c(0.05,0.01)
sigma_St <- c(0.01,0.01)
time_expert <- c(7,9)
mu_St <- 0.1866384#(upper_surv+lower_surv)/2
sigma_St <- 0.010#(upper_surv-lower_surv)/(2*1.96)
time_expert <- 7


init_fun <- function(...){list(alpha = 0.5, beta=c(0,1))}
init_fun()


expert.stan.weibull<- survHE::fit.models(formula, data = df_comb,
                                         distr = "wei", method = "hmc",
                                         mu_St = mu_St,
                                         sigma_St = sigma_St,
                                         time_expert = time_expert,
                                         id_St = 0, 
                                         St_indic = 1)#






fs1 <- flexsurvreg(Surv(recyrs, censrec) ~ group, data = bc,
                   dist = "weibull")


mods <- c("exp", "weibull", "gamma", "lnorm", "llogis", "gengamma")
formula <- Surv(recyrs, censrec) ~ as.factor(group)
m1 <- fit.models(formula = formula, data = bc, distr = mods)

View(m1)
print(m1, mod = 5)

plot(m1, add.km = TRUE)+
  scale_x_continuous(expand = c(0, 0), limits = c(0,NA)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA))

model.fit.plot(m1)

#Change to WAIC


psa <- make.surv(fit = m1, nsim = 1000, t = seq(.1, 63))

psa.plot(psa, xlab = "Extrapolated time",
         ylab = "Estimation of the survival curves",
         alpha = 0.2, col = c("dark grey", "black", "green"),
         main = "PSA to survival curves", cex.txt = .95)