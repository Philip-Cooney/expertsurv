rstan_options(javascript = FALSE)



fit.models <- function (formula = NULL, data, distr = NULL, method = "mle", 
                        ...){
  call = match.call()
  exArgs <- list(...)
  exArgs$formula <- formula
  exArgs$data = data
  exArgs$call = call
  if (is.null(formula)) {
    stop("You need to specify a model 'formula', e.g. 'formula=Surv(time,event)~treat'")
  }
  
  method <- tolower(method)
  if (!method %in% c("hmc", "inla", "mle")) {
    stop("Methods available for use are 'mle', 'hmc' or 'inla'")
  }
  method = check_distributions(method, distr)
  
  print(method)
  
  if (method == "mle") {
    res <- format_output_fit.models(lapply(distr, function(x) runMLE(x, 
                                                                     exArgs)), method, distr, formula, data)
  }
  if (method == "inla") {
    res <- format_output_fit.models(lapply(distr, function(x) runINLA(x, 
                                                                      exArgs)), method, distr, formula, data)
  }
  
  if (method == "hmc") {
    res <- format_output_fit.models(lapply(distr, function(x) runHMC(x, 
                                                                     exArgs)), method, distr, formula, data)
  }
  return(res)
}



runHMC <- function (x, exArgs) 
{
  if (!isTRUE(requireNamespace("rstan", quietly = TRUE))) {
    stop("You need to install the R package 'rstan'. Please run in your R terminal:\n install.packages('rstan')")
  }
  formula <- exArgs$formula
  data = exArgs$data
  availables <- load_availables()
  d3 <- manipulate_distributions(x)$distr3
  method <- "hmc"
  if (exists("chains", where = exArgs)) {
    chains <- exArgs$chains
  }
  else {
    chains <- 2
  }
  if (exists("iter", where = exArgs)) {
    iter <- exArgs$iter
  }
  else {
    iter <- 2000
  }
  if (exists("warmup", where = exArgs)) {
    warmup <- exArgs$warmup
  }
  else {
    warmup <- floor(iter/2)
  }
  if (exists("thin", where = exArgs)) {
    thin <- exArgs$thin
  }
  else {
    thin <- 1
  }
  if (exists("control", where = exArgs)) {
    control <- exArgs$control
  }
  else {
    control <- list(NULL)
  }
  if (exists("seed", where = exArgs)) {
    seed <- exArgs$seed
  }
  else {
    seed <- sample.int(.Machine$integer.max, 1)
  }
  if (exists("pars", where = exArgs)) {
    pars <- exArgs$pars
  }
  else {
    pars <- c("lambda_cens", "lambda_obs", "cens", 
              "d", "lp__", "loglambda_cens", 
              "loglambda_obs", "mu", "logP", 
              "linpred")
  }
  if (exists("include", where = exArgs)) {
    include <- exArgs$include
  }
  else {
    include <- FALSE
  }
  if (exists("k", where = exArgs)) {
    k <- exArgs$k
  }
  else {
    k <- 0
  }
  if (exists("cores", where = exArgs)) {
    cores <- exArgs$cores
  }
  else {
    cores <- 1
  }
  if (exists("init", where = exArgs)) {
    init <- exArgs$init
  }
  else {
    init = "random"
  }
  if (exists("save.stan", where = exArgs)) {
    save.stan <- exArgs$save.stan
  }
  else {
    save.stan = FALSE
  }
  if (exists("refresh", where = exArgs)) {
    refresh = exArgs$refresh
  }
  else {
    refresh = max(iter/10, 1)
  }
  d <- names(availables[[method]][match(d3, availables[[method]])])
  dso <- stanmodels[[d]]
  
  
  #save(d, file = paste0(pathway, "Extra Output/d.RData"))
  #save(dso, file = paste0(pathway, "Extra Output/dso.RData"))
  
  data.stan <- make_data_stan(formula, data, d3, exArgs)
  tic <- proc.time()
  model <- rstan::sampling(dso, data.stan, chains = chains, 
                           iter = iter, warmup = warmup, thin = thin, seed = seed, 
                           control = control, pars = pars, include = include, cores = cores, 
                           init = init, refresh = refresh)
  toc <- proc.time() - tic
  time_survHE <- toc[3]
  time_stan <- sum(rstan::get_elapsed_time(model))
  ics <- compute_ICs_stan(model, d3, data.stan)
  if (save.stan) {
    model_code <- attr(model@stanmodel, "model_code")
    con <- paste0(d, ".stan")
    writeLines(model_code, con = con)
    cat(paste0("Model code saved to the file: ", con, 
               "\n"))
  }
  model_name <- d3
  list(model = model, aic = ics$aic, bic = ics$bic, dic = ics$dic, 
       dic2 = ics$dic2, time2run = pmin(time_survHE, time_stan), 
       data.stan = data.stan, save.stan = save.stan, model_name = model_name)
}


tmpfun <- get("runHMC", envir = asNamespace("survHE"))
environment(runHMC) <- environment(tmpfun)
attributes(runHMC) <- attributes(tmpfun)  # don't know if this is really needed
assignInNamespace("runHMC", runHMC, ns="survHE")




make_data_stan <- function (formula, data, distr3, exArgs = globalenv()){
  availables <- load_availables()
  method <- "hmc"
  formula_temp <- update(formula, paste(all.vars(formula, data)[1], 
                                        "~", all.vars(formula, data)[2], "+."))
  mf <- as_tibble(model.frame(formula_temp, data)) %>% rename(time = 1, 
                                                              event = 2) %>% rename_if(is.factor, .funs = ~gsub("as.factor[( )]", 
                                                                                                                "", .x)) %>% rename_if(is.factor, .funs = ~gsub("[( )]", 
                                                                                                                                                                "", .x)) %>% bind_cols(as_tibble(model.matrix(formula_temp, 
                                                                                                                                                                                                              data)) %>% select(contains("Intercept"))) %>% select(time, 
                                                                                                                                                                                                                                                                   event, contains("Intercept"), everything()) %>% 
    tibble::rownames_to_column("ID")
  
  ####Code Change Here
  
  if (distr3 %in% c("gam", "gga", "gef")) {
    data.stan <- list(t = (mf %>% filter(event == 1))$time, 
                      d = (mf %>% filter(event == 0))$time, n_obs = mf %>% 
                        filter(event == 1) %>% with(nrow(.)), n_cens = mf %>% 
                        filter(event == 0) %>% with(nrow(.)))
    data.stan$X_obs <- matrix(model.matrix(formula, data)[(mf %>% 
                                                             filter(event == 1))$ID, ], nrow = data.stan$n_obs)
    data.stan$H = ncol(data.stan$X_obs)
    
    if(!is.null(exArgs$a0)){
      # print(which(data$event == 1))
      data.stan$a0_obs <-  exArgs$a0[which(mf$event == 1)]
      
    }else{
      data.stan$a0_obs <- rep(1,length(mf$event == 1))
    }
    
    if(data.stan$n_cens >0){
      data.stan$X_cens <- matrix(model.matrix(formula, data)[(mf %>% 
                                                                filter(event == 0))$ID, ], nrow = data.stan$n_cens)
      
      if(!is.null(exArgs$a0)){
        data.stan$a0_cens <-  exArgs$a0[which(mf$event == 0)]
      }else{
        # data.stan$a0_cens <- rep(1,length(data$event == 0))
      }
      
      
    }else{
      data.stan$X_cens <- matrix(0, nrow = 0, ncol = data.stan$H)
      data.stan$a0_cens <- matrix(0, nrow = 0, ncol = 1)
    }
    
    if (data.stan$H == 1) {
      data.stan$X_obs <- cbind(data.stan$X_obs, rep(0, 
                                                    data.stan$n_obs))
      data.stan$X_cens <- cbind(data.stan$X_cens, rep(0, 
                                                      data.stan$n_cens))
      data.stan$H <- ncol(data.stan$X_obs)
    }
  }
  ######
  #print(data.stan$a0_obs)
  
  
  
  
  if (distr3 %in% c("exp", "gom", "wei", 
                    "wph", "llo", "lno", "pow")) {
    data.stan <- list(t = (mf$time), d = mf$event, n = nrow(mf), 
                      X = matrix(model.matrix(formula, data), nrow = nrow(mf)), 
                      H = ncol(model.matrix(formula, data)))
    if (data.stan$H == 1) {
      data.stan$X <- cbind(data.stan$X, rep(0, data.stan$n))
      data.stan$H <- ncol(data.stan$X)
    }
  }
  if (distr3 == "rps") {
    if (exists("k", where = exArgs)) {
      k <- exArgs$k
    }
    else {
      k <- 0
    }
    knots <- quantile(log((mf %>% filter(event == 1))$time), 
                      seq(0, 1, length = k + 2))
    B <- flexsurv::basis(knots, log(mf$time))
    B_expert <- flexsurv::basis(knots, log(exArgs$time_expert))
    DB <- flexsurv::dbasis(knots, log(mf$time))
    mm <- model.matrix(formula, data)[, -1]
    if (length(mm) < 1) {
      mm <- matrix(rep(0, nrow(mf)), nrow = nrow(mf), ncol = 2)
    }
    if (is.null(dim(mm))) {
      mm <- cbind(mm, rep(0, length(mm)))
    }
    data.stan <- list(t = mf$time, d = mf$event, n = nrow(mf), 
                      M = k, X = mm, H = ncol(mm), B = B, DB = DB, mu_gamma = rep(0,k + 2),
                      sigma_gamma = rep(5, k + 2), knots = knots, B_expert = B_expert)
  }
  data.stan$mu_beta = rep(0, data.stan$H)
  if (distr3 %in% c("gef", "gga", "lno")) {
    data.stan$sigma_beta <- rep(100, data.stan$H)
  }
  else {
    data.stan$sigma_beta <- rep(5, data.stan$H)
  }
  if (distr3 == "gef") {
    data.stan$a_sigma = data.stan$b_sigma = 0.1
    data.stan$mu_P = 0
    data.stan$sigma_P = 0.5
    data.stan$mu_Q = 0
    data.stan$sigma_Q = 2.5
  }
  else if (distr3 == "gga") {
    data.stan$a_sigma = data.stan$b_sigma = 0.1
    data.stan$mu_Q = 0
    data.stan$sigma_Q = 100
  }
  else if (distr3 %in% c("gam", "llo", "wei", 
                         "wph", "gom")) {
    data.stan$a_alpha = data.stan$b_alpha = 0.1
  }
  else if (distr3 %in% c("lno")) {
    data.stan$a_alpha = 0
    data.stan$b_alpha = 5
  }
  d <- names(availables[[method]][match(distr3, availables[[method]])])
  priors <- list()
  if (exists("priors", where = exArgs)) {
    abbrs = manipulate_distributions(names(exArgs$priors))$distr3
    pos = grep(distr3, abbrs)
    if (length(pos) > 0) {
      priors = exArgs$priors[[pos]]
    }
  }
  if (!is.null(priors$mu_beta)) {
    data.stan$mu_beta = priors$mu_beta
  }
  if (!is.null(priors$sigma_beta)) {
    data.stan$sigma_beta <- priors$sigma_beta
  }
  if (!is.null(priors$mu_gamma) & distr3 == "rps") {
    data.stan$mu_gamma <- priors$mu_gamma
  }
  if (!is.null(priors$sigma_gamma) & distr3 == "rps") {
    data.stan$sigma_gamma <- priors$sigma_gamma
  }
  if (!is.null(priors$a_sigma)) {
    data.stan$a_sigma = priors$a_sigma
  }
  if (!is.null(priors$b_sigma)) {
    data.stan$b_sigma = priors$b_sigma
  }
  if (!is.null(priors$mu_P)) {
    data.stan$mu_P = priors$mu_P
  }
  if (!is.null(priors$sigma_P)) {
    data.stan$sigma_P = priors$sigma_P
  }
  if (!is.null(priors$mu_Q)) {
    data.stan$mu_Q = priors$mu_Q
  }
  if (!is.null(priors$sigma_Q)) {
    data.stan$sigma_Q = priors$sigma_Q
  }
  if (!is.null(priors$a_alpha)) {
    data.stan$a_alpha = priors$a_alpha
  }
  if (!is.null(priors$b_alpha)) {
    data.stan$b_alpha = priors$b_alpha
  }
  
  data.stan$St_indic <- exArgs$St_indic 
  
  if(!is.null(exArgs$mu_St)){
    
    if (distr3 %in% c("gam", "gga", "gef")){
      data.stan$id_St <-  min(which(data.stan[["X_obs"]][,2] == exArgs$id_St))
      
    }else{
      data.stan$id_St <-  min(which(data.stan[["X"]][,2] == exArgs$id_St))
      
    }
    
    
    if(length(exArgs$mu_St) == 1){
      data.stan$mu_St <- as.array(exArgs$mu_St)
      data.stan$sigma_St <- as.array(exArgs$sigma_St)
      data.stan$time_expert <- as.array(exArgs$time_expert)
      data.stan$n_time_expert <- length(exArgs$time_expert)
      
    }else{
      data.stan$mu_St <- exArgs$mu_St
      data.stan$sigma_St <- exArgs$sigma_St
      data.stan$time_expert <- exArgs$time_expert
      data.stan$n_time_expert <- length(exArgs$time_expert)
    }
    
  }
  
  if(!is.null(exArgs$St_indic)){
    if(exArgs$St_indic == 0){
      
      #Dummy data so stan runs
      data.stan$n_time_expert <-  1
      data.stan$id_St <- 0
      data.stan$mu_St <- as.array(0)
      data.stan$sigma_St <- as.array(0)
      data.stan$time_expert <- as.array(0)
      
      if (distr3 %in% c("gam", "gga", "gef")){
        data.stan$id_trt <-  min(which(data.stan[["X_obs"]][,2] == 1))
        data.stan$id_comp <-  min(which(data.stan[["X_obs"]][,2] == 0))
      }else{
        data.stan$id_trt <-  min(which(data.stan[["X"]][,2] == 1))
        data.stan$id_comp <-  min(which(data.stan[["X"]][,2] == 0)) 
      }
      
      
      data.stan$mu_diff <- exArgs$mu_diff
      data.stan$sigma_diff <- exArgs$sigma_diff
    }else{
      data.stan$id_trt <-  0
      data.stan$id_comp <- 0
      data.stan$mu_diff <- 0
      data.stan$sigma_diff <- 0
    }
    
    
    
  }
  
  #Power prior
  if (distr3 %in% c("exp", "gom", "wei", 
                    "wph", "llo", "lno", "pow", "rps")) { # 
    
    if(!is.null(exArgs$a0)){
      data.stan$a0 <- exArgs$a0
    }else{
      data.stan$a0 <- rep(1, nrow(data))
    }
  }
  
  #data.stan$exArgs <- exArgs
  save(data.stan, file = paste0(pathway, "Extra Output/datastan.RData"))
  #
  data.stan	
}



tmpfun <- get("make_data_stan", envir = asNamespace("survHE"))
environment(make_data_stan) <- environment(tmpfun)
attributes(make_data_stan) <- attributes(tmpfun)  # don't know if this is really needed
assignInNamespace("make_data_stan", make_data_stan, ns="survHE")



tmpfun <- get("fit.models", envir = asNamespace("survHE"))
environment(fit.models) <- environment(tmpfun)
attributes(fit.models) <- attributes(tmpfun)  # don't know if this is really needed
assignInNamespace("fit.models", fit.models, ns="survHE")


stanmodels <- survHE:::stanmodels

stan.code.weibullAF <- stan_model(file = paste0(pathway, "WeibullAF.stan"))
stanmodels[["WeibullAF"]]<- stan.code.weibullAF
stan.code.exponential <- stan_model(file = paste0(pathway, "Exponential.stan"))
stanmodels[["Exponential"]]<- stan.code.exponential

stan.code.gamma <- stan_model(file = paste0(pathway, "Gamma.stan"))
stanmodels[["Gamma"]]<- stan.code.gamma

stan.code.weibullPH <- stan_model(file = paste0(pathway, "WeibullPH.stan"))
stanmodels[["WeibullPH"]]<- stan.code.weibullPH

stan.code.logLogistic <- stan_model(file = paste0(pathway, "logLogistic.stan"))
stanmodels[["logLogistic"]]<- stan.code.logLogistic

stan.code.logNormal <- stan_model(file = paste0(pathway, "logNormal.stan"))
stanmodels[["logNormal"]]<- stan.code.logNormal

stan.code.gompertz <- stan_model(file = paste0(pathway, "Gompertz.stan"))
stanmodels[["Gompertz"]]<- stan.code.gompertz


stan.code.GenGamma <- stan_model(file = paste0(pathway, "GenGamma.stan"))
stanmodels[["GenGamma"]]<- stan.code.GenGamma

stan.code.GenF <- stan_model(file = paste0(pathway, "GenF.stan"))
stanmodels[["GenF"]]<- stan.code.GenF


stan.code.RP <- stan_model(file = paste0(pathway, "RP.stan"))
stanmodels[["RP"]]<- stan.code.RP


tmpfun <- get("stanmodels", envir = asNamespace("survHE"))
environment(stanmodels) <- environment(tmpfun)
attributes(stanmodels) <- attributes(tmpfun)  # don't know if this is really needed
assignInNamespace("stanmodels", stanmodels, ns="survHE")

# load(file = "C:/Users/phili/Desktop/dso.RData")
# load(file = "C:/Users/phili/Desktop/datastan.RData")
# model <- rstan::sampling(stan.code.RP, data.stan, chains = 1)

tmpfun <- get("check_distributions", envir = asNamespace("survHE"))
environment(check_distributions) <- environment(tmpfun)
attributes(check_distributions) <- attributes(tmpfun)  # don't know if this is really needed
assignInNamespace("check_distributions", check_distributions, ns="survHE")



#Function 
stopped = T

#stop("Don't need to run the validation code")

if(stopped == T){
  
}else{
  
  #https://stackoverflow.com/questions/59900118/integrate-1d-syntax-and-interactively-debugging-it
  
  code2 = "functions{

real integrand(real x,  
real xc,
real[] theta,
real[] x_r,
int[] x_i ){

// separate the parameters
real shape = theta[1];
real rate = theta[2];

return(exp(-rate/shape * (exp(shape * x) - 1)));

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
  
}"

rescale_stats_hmc_gam <- function (table, x){
  rate <- matrix(table[grep("rate", rownames(table)), 
  ], ncol = 4)
  rownames(rate) <- "rate"
  shape <- matrix(table[grep("alpha", rownames(table)), 
  ], ncol = 4)
  rownames(shape) <- "shape"
  effects = add_effects_hmc(table,x)
  res <- rbind(shape, rate, effects)
  if (is.null(dim(res))) {
    names(res) <- c("mean", "se", "L95%", 
                    "U95%")
  }
  else {
    colnames(res) <- c("mean", "se", "L95%", 
                       "U95%")
  }
  return(res)
}


tmpfun <- get("rescale_stats_hmc_gam", envir = asNamespace("survHE"))
environment(rescale_stats_hmc_gam) <- environment(tmpfun)
attributes(rescale_stats_hmc_gam) <- attributes(tmpfun)  
assignInNamespace("rescale_stats_hmc_gam", rescale_stats_hmc_gam, ns="survHE")

#Validation

#expose_stan_functions(stanc(model_code = code2))
#expec_surv(theta = c(2,1), x_r = double()) # yields 1

#library(flexsurv)
#flexsurv:::mean_gompertz(shape = 2 , rate = 1)
# Surv_diff(shape = as.array(6.704507e-06),
#           rate_trt = as.array(exp(-0.236889 -0.4293956)),
#           rate_comp = as.array(0.7890789),
#           x_r = as.array(0))



code_surv <- "functions{
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
}"

# expose_stan_functions(stanc(model_code = code_surv))
# 
# Sind(time = as.double(1),
#            mu = as.double(1),
#            sigma = as.double(1),
#            Q = as.double(1))

mu <- 1
sigma <- 0.5
Q <- 0.5

mean_gengamma(mu,sigma,Q)

a <- exp(mu + (sigma*2*log(Q))/Q)
p <- Q/sigma
d <- 1/(sigma*Q)

(a*gamma((d+1)/p))/(gamma(d/p))




Surv.genf <- "functions{
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
}"


expose_stan_functions(stanc(model_code = Surv.genf))

Sind(time = as.double(1),
     mu = as.double(1),
     sigma = as.double(1),
     Q = as.double(1),
     P = as.double(1))


pgenf(1,mu = 1,1,1,1, lower.tail = F)

Surv.rps <-"functions{
real Sind(vector gamma, row_vector B, real linpred) {
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
}"
expose_stan_functions(stanc(model_code = Surv.rps))

#data.stan[["knots"]]
#The first and the last knots are boundary knots, kmin and kmax
knots.exam <- c(9.226319954, -2.384776187, -0.831818545, -0.002183787)

q <- 0.1
gamma.res <- c(1,1,1,1)
psurvspline(q = q,gamma =gamma.res,knots = knots.exam, lower.tail = F)
Sind(gamma.res,flexsurv::basis(knots.exam, log(q)),0)


cuber <- function(x){
  if (x < 0){
    return(0)
  }else{
    return(x^3)
  }
}


q <- 0.1
gamma.res <- c(1,1,1,1)

psurvspline(q = q,gamma =gamma.res,knots = knots.exam, lower.tail = F)

d <- flexsurv:::dbase.survspline(q = q, gamma = gamma.res, knots = knots.exam, 
                                 scale = "hazard")

#ind just checks if q (survival time) is not NA or 0 



flexsurv::basis(knots.exam,log(3.648244e-01))
flexsurv::basis(knots.exam, c(log(3.648244e-01),log(3.648244e-01)))




flexsurv:::dbase.survspline
flexsurv::psurvspline
flexsurv:::Slink
#The Basis functions are as follows;
# 1 is refers to the intercept
# 2 is log(time)
# 3 is the first internal knot
# 4 is the second (last) internal knot

log_t <- log(3.648244e-01)
log_t >data.stan[["knots"]]

#3rd knot
k_max <- tail(data.stan[["knots"]], 1)
k_min <- head(data.stan[["knots"]], 1)
k_m <- data.stan[["knots"]][2]

lambda_m <- (k_max-k_m)/(k_max - k_min)
cuber(log_t - k_m) -lambda_m*cuber(log_t - k_min) - (1-lambda_m)*(log_t - k_max)


# We get eta from the basis functions 
# and multiply by gamma adding covariates as appropriate
# we then link eta to survival by the following expression 
# St = exp(-exp(eta))



flexsurv:::basis_vector(data.stan[["knots"]], log(3.648244e-01))

basis_vector(data.stan[["knots"]], log(3.648244e-01))


sourceCpp(code = '
           #include <Rcpp.h>
            using namespace Rcpp;
            
            namespace {
              inline double cuber(const double x) {
                if (x <= 0) {
                  return 0;
                } else {
                  return x * x * x;
                }
              }
              
              inline double dCuber(const double x) {
                if (x <= 0) {
                  return 0;
                } else {
                  return 3 * x * x;
                }
              }
            }
          
          // [[Rcpp::export]]
          Rcpp::NumericVector
          basis_vector(const Rcpp::NumericVector& knots,
                       const Rcpp::NumericVector& x) {
            if (knots.size() < 2) {
              throw std::runtime_error("Bad knots.");
            }
            Rcpp::NumericMatrix result(x.size(), knots.size());
            result(Rcpp::_, 0) = Rcpp::rep(1, x.size());
            result(Rcpp::_, 1) = x;
            Rcout << "The value of x : " << x ;


            for (R_xlen_t ind=0; ind < knots.size() - 2; ++ind) {
              const double last = *(knots.end() - 1);
              const double first = *(knots.begin());
              const double lam = (last - knots[ind + 1]) / (last - first);
              result(Rcpp::_, ind + 2) =
                Rcpp::sapply(x - knots[ind + 1], cuber)
              - lam * Rcpp::sapply(x - first, cuber)
              - (1 - lam) * Rcpp::sapply(x - last, cuber);
 
            }


            return result;
          }')




rps_lpdf <-"functions{
real rps_lpdf(vector t, vector d, vector gamma, matrix B, matrix DB, vector linpred) {
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
  lprob = sum(log_lik);
  return lprob;
  }
}"
expose_stan_functions(stanc(model_code = rps_lpdf))
beta_inits <- c(0,-0.5)
linpred <- data.stan$X%*%beta_inits

gamma_inits <- c(0.789,0.744,0.027,-0.05)
rps_lpdf(t =data.stan$t,d = data.stan$d, gamma= gamma_inits, 
         B = data.stan$B, DB = data.stan$DB,linpred =linpred)


dim(linpred)



gomp_surv <-"functions{

// Defines difference in expected survival
real Surv_diff ( real shape, real rate_trt, real rate_comp) {
  real Surv_diff;
  real Surv_trt;
  real Surv_comp;
  real s;
  s = 0.0001;
  
  Surv_trt  = (1/shape)*exp(rate_trt/shape)*tgamma(s)*(1-gamma_cdf(rate_trt/shape,s,1));
  Surv_comp  = (1/shape)*exp(rate_comp/shape)*tgamma(s)*(1-gamma_cdf(rate_comp/shape,s,1));
  
  Surv_diff = Surv_trt - Surv_comp;
  return Surv_diff;
}
}"


expose_stan_functions(stanc(model_code = gomp_surv))

#Surv_diff(1,1,3)

gamma_lpdf2 <- "functions {
  
  
  
  real gamma2_lpdf(vector t, real alpha,  vector beta, vector a0) {
    vector[num_elements(t)] log_lik;
    vector[num_elements(t)] prob;
    real lprob;
    // Constructs the log-density for each observation
    // don't use gamma_lpdf, use the function derived by yourself
	
    for (i in 1:num_elements(t)) {
      prob[i] = gamma_lpdf(t[i]| alpha, beta[i]);
	  //prob[i] = alpha*log(beta[i])-lgamma(alpha)+ (alpha-1)*log(t[i]) - beta*t[i]
	  
    }
    // And the total log-density (as a sum of the individual terms)
    lprob = dot_product(prob, a0);
    return lprob;
  }
  
  
	
}"

expose_stan_functions(stanc(model_code = gamma_lpdf2))


t = c(1,2,3)
alpha = as.array(1)
beta = c(1,2,3)
a0 = c(1,1,0)
gamma2_lpdf(t =t ,alpha = alpha, beta = beta, a0)


sum(dgamma(t,alpha, beta, log = T)*a0)

}







