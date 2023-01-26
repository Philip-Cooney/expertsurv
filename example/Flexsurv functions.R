
library(flexsurv)
library(mvtnorm)
library(survival)
library(rstan)


logLikFactory <- function (Y, X = 0, weights, bhazard, rtrunc, dlist, inits, dfns,
                           aux, mx, fixedpars = NULL, expert_opinion){
  pars <- inits
  npars <- length(pars)
  nbpars <- length(dlist$pars)
  insert.locations <- setdiff(seq_len(npars), fixedpars)
  event <- Y[, "status"] == 1
  event.times <- Y[event, "time1"]
  left.censor <- Y[!event, "time2"]
  right.censor <- Y[!event, "time1"]
  event.weights <- weights[event]
  no.event.weights <- weights[!event]
  par.transform <- buildTransformer(inits, nbpars, dlist)
  aux.pars <- buildAuxParms(aux, dlist)
  default.offset <- rep.int(0, length(event.times))
  do.hazard <- any(bhazard > 0)
  loglik <- rep.int(0, nrow(Y))
  function(optpars, ...) {
    pars[insert.locations] <- optpars
    raw.pars <- pars
    pars <- as.list(pars)
    pars.event <- pars.nevent <- pars
    if (npars > nbpars) {
      beta <- raw.pars[(nbpars + 1):npars]
      for (i in dlist$pars) {
        pars[[i]] <- pars[[i]] + X[, mx[[i]], drop = FALSE] %*%
          beta[mx[[i]]]
        pars.event[[i]] <- pars[[i]][event]
        pars.nevent[[i]] <- pars[[i]][!event]
      }
    }
    fnargs <- c(par.transform(pars), aux.pars)
    fnargs.event <- c(par.transform(pars.event), aux.pars)
    fnargs.nevent <- c(par.transform(pars.nevent), aux.pars)
    dargs <- fnargs.event
    dargs$x <- event.times
    dargs$log <- TRUE
    logdens <- do.call(dfns$d, dargs)
    if (any(!event)) {
      pmaxargs <- fnargs.nevent
      pmaxargs$q <- left.censor
      pmax <- do.call(dfns$p, pmaxargs)
      pmax[pmaxargs$q == Inf] <- 1
      pargs <- fnargs.nevent
      pargs$q <- right.censor
      pmin <- do.call(dfns$p, pargs)
    }
    targs <- fnargs
    targs$q <- Y[, "start"]
    plower <- do.call(dfns$p, targs)
    targs$q <- rtrunc
    pupper <- do.call(dfns$p, targs)
    pupper[rtrunc == Inf] <- 1
    pobs <- pupper - plower
    if (do.hazard) {
      pargs <- fnargs.event
      pargs$q <- event.times
      pminb <- do.call(dfns$p, pargs)
      loghaz <- logdens - log(1 - pminb)
      offseti <- log(1 + bhazard[event]/exp(loghaz) * weights[event])
    }
    else {
      offseti <- default.offset
    }
    loglik[event] <- (logdens * event.weights) + offseti
    if (any(!event))
      loglik[!event] <- (log(pmax - pmin) * no.event.weights)
    loglik <- loglik - log(pobs) * weights
    #browser()
    pargs.surv <- fnargs.nevent
    #Expert time
    
    if(!is.null(expert_opinion)){
      expert_opinion
      pargs.surv$q <- expert_opinion$time
      psurv_expert <- 1 - do.call(dfns$p, pargs.surv)
      
      param_expert <- expert_opinion$param_expert
      # mu_surv = c(0.2)#, 0.199999, 0.199998)
      # sd_surv = rep(0.05,length(mu_surv))
      #browser()
      LL_expert <- rep(NA, length(expert_opinion$time))
      for(q in 1:length(expert_opinion$time)){
       param_expert_mat  <- abind::adrop(param_expert[,,q, drop = F],drop = 3)
        LL_expert[q] <-  log_density_dist(param_expert_mat,psurv_expert[q],nrow(param_expert_mat), expert_opinion$pool)
      }
      
      #df_expert <- data.frame(Surv = psurv_expert, mu_surv = expert_opinion$mu_surv , sd_surv = expert_opinion$sd_surv)
      #LL_expert<- apply(df_expert, 1, function(x){dnorm(x["Surv"], mean = x["mu_surv"], sd = x["sd_surv"], log = T)})
      
      ret <- -sum(loglik, LL_expert)
    }else{
      ret <- -sum(loglik)
    }
    
    attr(ret, "indiv") <- loglik
    #Surv_vals <- do.call(dfns$p, list(q =c(10,12)))
    #print(dfns$p)
    #print(pargs)
    #There will always be a darg but not always a parg, however you would never need an opinion if everyone died!
    #print(pargs.surv)
    #print(psurv_expert)
    ret
  }
}



getOption("flexsurv.test.analytic.derivatives")



flexsurvreg <- function (formula, anc = NULL, data, weights, bhazard, rtrunc,
                         subset, na.action, dist, inits, fixedpars = NULL, dfns = NULL,
                         aux = NULL, cl = 0.95, integ.opts = NULL, sr.control = survreg.control(),
                         hessian = TRUE, hess.control = NULL, expert_opinion = NULL, ...)
{
  
  call <- match.call()
  if (missing(dist))
    stop("Distribution \"dist\" not specified")
  dlist <- parse.dist(dist)
  dfns <- form.dp(dlist, dfns, integ.opts)
  parnames <- dlist$pars
  check.formula(formula, dlist)
  anc <- anc_from_formula(formula, anc, dlist)
  ancnames <- setdiff(parnames, dlist$location)
  forms <- c(location = get.locform(formula, ancnames), anc)
  names(forms)[[1]] <- dlist$location
  indx <- match(c("formula", "data", "weights",
                  "bhazard", "rtrunc", "subset", "na.action"),
                names(call), nomatch = 0)
  if (indx[1] == 0)
    stop("A \"formula\" argument is required")
  temp <- call[c(1, indx)]
  temp[[1]] <- as.name("model.frame")
  f2 <- concat.formulae(formula, forms)
  temp[["formula"]] <- f2
  if (missing(data))
    temp[["data"]] <- environment(formula)
  m <- eval(temp, parent.frame())
  m <- droplevels(m)
  attr(m, "covnames") <- attr(f2, "covnames")
  attr(m, "covnames.orig") <- intersect(colnames(m),
                                        attr(f2, "covnames.orig"))
  Y <- check.flexsurv.response(model.extract(m, "response"))
  mml <- mx <- vector(mode = "list", length = length(dlist$pars))
  names(mml) <- names(mx) <- c(dlist$location, setdiff(dlist$pars,
                                                       dlist$location))
  for (i in names(forms)) {
    mml[[i]] <- model.matrix(forms[[i]], m)
    mx[[i]] <- length(unlist(mx)) + seq_len(ncol(mml[[i]][,
                                                          -1, drop = FALSE]))
  }
  X <- compress.model.matrices(mml)
  weights <- model.extract(m, "weights")
  if (is.null(weights))
    weights <- m$"(weights)" <- rep(1, nrow(X))
  bhazard <- model.extract(m, "bhazard")
  if (is.null(bhazard))
    bhazard <- rep(0, nrow(X))
  rtrunc <- model.extract(m, "rtrunc")
  if (is.null(rtrunc))
    rtrunc <- rep(Inf, nrow(X))
  dat <- list(Y = Y, m = m, mml = mml)
  ncovs <- length(attr(m, "covnames.orig"))
  ncoveffs <- ncol(X)
  nbpars <- length(parnames)
  npars <- nbpars + ncoveffs
  if (missing(inits) && is.null(dlist$inits))
    stop("\"inits\" not supplied, and no function to estimate them found in the custom distribution list")
  if (missing(inits) || any(is.na(inits))) {
    yy <- ifelse(Y[, "status"] == 3 & is.finite(Y[,
                                                  "time2"]), (Y[, "time1"] + Y[, "time2"])/2,
                 Y[, "time1"])
    wt <- yy * weights * length(yy)/sum(weights)
    dlist$inits <- expand.inits.args(dlist$inits)
    inits.aux <- c(aux, list(forms = forms, data = if (missing(data)) NULL else data,
                             weights = temp$weights, control = sr.control, counting = (attr(model.extract(m,
                                                                                                          "response"), "type") == "counting")))
    auto.inits <- dlist$inits(t = wt, mf = m, mml = mml,
                              aux = inits.aux)
    if (!missing(inits) && any(is.na(inits)))
      inits[is.na(inits)] <- auto.inits[is.na(inits)]
    else inits <- auto.inits
  }
  if (!is.numeric(inits))
    stop("initial values must be a numeric vector")
  nin <- length(inits)
  if (nin < npars && ncoveffs > 0)
    inits <- c(inits, rep(0, length = npars - nin))
  else if (nin > npars) {
    inits <- inits[1:npars]
    warning("Initial values are a vector length > ",
            npars, ": using only the first ", npars)
  }
  else if (nin < nbpars) {
    stop("Initial values are a vector length ", nin,
         ", but distribution has ", nbpars, " parameters")
  }
  for (i in 1:nbpars) inits[i] <- dlist$transforms[[i]](inits[i])
  outofrange <- which(is.nan(inits) | is.infinite(inits))
  if (any(outofrange)) {
    plural <- if (length(outofrange) > 1)
      "s"
    else ""
    stop("Initial value", plural, " for parameter",
         plural, " ", paste(outofrange, collapse = ","),
         " out of range")
  }
  names(inits) <- c(parnames, colnames(X))
  check.fixedpars(fixedpars, npars)
  #print((is.logical(fixedpars) && fixedpars == TRUE) || (is.numeric(fixedpars) &&
  #identical(as.vector(fixedpars), 1:npars)))
  if ((is.logical(fixedpars) && fixedpars == TRUE) || (is.numeric(fixedpars) &&
                                                       identical(as.vector(fixedpars), 1:npars))) {
    minusloglik <- minusloglik.flexsurv(inits, Y = Y, X = X,
                                        weights = weights, bhazard = bhazard, rtrunc = rtrunc,
                                        dlist = dlist, inits = inits, dfns = dfns, aux = aux,
                                        mx = mx)
    res.t <- matrix(inits, ncol = 1)
    inits.nat <- inits
    for (i in 1:nbpars) inits.nat[i] <- dlist$inv.transforms[[i]](inits[i])
    res <- matrix(inits.nat, ncol = 1)
    dimnames(res) <- dimnames(res.t) <- list(names(inits),
                                             "est")
    ret <- list(res = res, res.t = res.t, npars = 0, loglik = -as.vector(minusloglik),
                logliki = attr(minusloglik, "indiv"))
  }
  else {
    optpars <- inits[setdiff(1:npars, fixedpars)]
    optim.args <- list(...)
    if (is.null(optim.args$method)) {
      optim.args$method <- "BFGS"
    }
    
    #print(dfns$deriv)
    #gr <- if (dfns$deriv)
    #      Dminusloglik.flexsurv
    #else NULL
    gr <- NULL
    optim.args <- c(optim.args, list(par = optpars, fn = logLikFactory(Y = Y,
                                                                       X = X, weights = weights, bhazard = bhazard, rtrunc = rtrunc,
                                                                       inits = inits, dlist = dlist, dfns = dfns, aux = aux,
                                                                       mx = mx, fixedpars = fixedpars, expert_opinion = expert_opinion), gr = gr, Y = Y,
                                     X = X, weights = weights, bhazard = bhazard, rtrunc = rtrunc,
                                     dlist = dlist, inits = inits, dfns = dfns, aux = aux,
                                     mx = mx, fixedpars = fixedpars, hessian = hessian))
    opt <- do.call("optim", optim.args)
    est <- opt$par
    if (hessian && all(!is.na(opt$hessian)) && all(!is.nan(opt$hessian)) &&
        all(is.finite(opt$hessian)) && all(eigen(opt$hessian)$values >
                                           0)) {
      cov <- .hess_to_cov(opt$hessian, hess.control$tol.solve,
                          hess.control$tol.evalues)
      se <- sqrt(diag(cov))
      if (!is.numeric(cl) || length(cl) > 1 || !(cl > 0) ||
          !(cl < 1))
        stop("cl must be a number in (0,1)")
      lcl <- est - qnorm(1 - (1 - cl)/2) * se
      ucl <- est + qnorm(1 - (1 - cl)/2) * se
    }
    else {
      if (hessian)
        warning("Optimisation has probably not converged to the maximum likelihood - Hessian is not positive definite. ")
      cov <- lcl <- ucl <- se <- NA
    }
    res <- cbind(est = inits, lcl = NA, ucl = NA, se = NA)
    res[setdiff(1:npars, fixedpars), ] <- cbind(est, lcl,
                                                ucl, se)
    colnames(res) <- c("est", paste(c("L", "U"),
                                    round(cl * 100), "%", sep = ""), "se")
    res.t <- res
    for (i in 1:nbpars) {
      res[i, 1:3] <- dlist$inv.transforms[[i]](res[i, 1:3])
      if (identical(body(dlist$transforms[[i]]), body(log)))
        res[i, "se"] <- exp(res.t[i, "est"]) *
          res.t[i, "se"]
      else if (identical(body(dlist$transforms[[i]]), body(logh)))
        res[i, "se"] <- dexph(res.t[i, "est"]) *
          res.t[i, "se"]
      else if (!identical(dlist$transforms[[i]], identity))
        res[i, "se"] <- NA
    }
    minusloglik <- minusloglik.flexsurv(res.t[, "est"],
                                        Y = Y, X = X, weights = weights, bhazard = bhazard,
                                        rtrunc = rtrunc, dlist = dlist, inits = inits, dfns = dfns,
                                        aux = aux, mx = mx, expert_opinion = expert_opinion)
    ret <- list(res = res, res.t = res.t, cov = cov, coefficients = res.t[,
                                                                          "est"], npars = length(est), fixedpars = fixedpars,
                optpars = setdiff(1:npars, fixedpars), loglik = -opt$value,
                logliki = attr(minusloglik, "indiv"), cl = cl,
                opt = opt)
  }
  ret <- c(list(call = call, dlist = dlist, aux = aux, ncovs = ncovs,
                ncoveffs = ncoveffs, mx = mx, basepars = 1:nbpars, covpars = if (ncoveffs >
                                                                                 0) (nbpars + 1):npars else NULL, AIC = -2 * ret$loglik +
                  2 * ret$npars, data = dat, datameans = colMeans(X),
                N = nrow(dat$Y), events = sum(dat$Y[, "status"] ==
                                                1), trisk = sum(dat$Y[, "time"]), concat.formula = f2,
                all.formulae = forms, dfns = dfns), ret)
  if (isTRUE(getOption("flexsurv.test.analytic.derivatives")) &&
      (dfns$deriv)) {
    if (is.logical(fixedpars) && fixedpars == TRUE) {
      optpars <- inits
      fixedpars = FALSE
    }
    ret$deriv.test <- deriv.test(optpars = optpars, Y = Y,
                                 X = X, weights = weights, bhazard = bhazard, rtrunc = rtrunc,
                                 dlist = dlist, inits = inits, dfns = dfns, aux = aux,
                                 mx = mx, fixedpars = fixedpars)
  }
  class(ret) <- "flexsurvreg"
  
  
  ret
}

minusloglik.flexsurv <- function (optpars, Y, X = 0, weights, bhazard, rtrunc, dlist,
                                  inits, dfns, aux, mx, fixedpars = NULL, expert_opinion){
  logLikFactory(Y = Y, X = X, weights = weights, bhazard = bhazard,
                rtrunc = rtrunc, dlist = dlist, inits = inits, dfns = dfns,
                aux = aux, mx = mx, fixedpars = fixedpars, expert_opinion = expert_opinion)(optpars)
}





model_code_density <-
  'functions {
  real log_density_dist(matrix params,
                        real x,int num_expert, int pool_type){


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
'
rstan::expose_stan_functions(stanc(model_code = model_code_density))


make_data_expert <- function(param_expert, times_expert = NULL){
  n.experts <- c()
  
  for(i in 1:length(param_expert)){
    n.experts <- c(n.experts, nrow(param_expert[[i]]))
  }
  
  data_dist_ind <- num_param <- matrix(-999.2,nrow = max(n.experts), ncol =  length(times_expert))
  expert.array <- array(-999.2,dim = c(max(n.experts),5,length(times_expert)))
  
  for(i in 1:length(times_expert)){
    lk_up_dist <- c("norm", "t", "gamma", "lnorm","beta")
    dist_fit <- param_expert[[i]][,1]
    if(length(dist_fit) - length(expert.array[,1,i])){
      dist_fit <- c(dist_fit, rep(-999.2,length(dist_fit) - length(expert.array[,1,i])))
    }
    expert.array[,1,i] <- as.numeric(sapply(dist_fit, function(x){which(x==lk_up_dist)}))
    weight_vec <- param_expert[[i]][,2]
    expert.array[1:length(weight_vec),2,i] <- weight_vec
    expert.array[1:nrow(param_expert[[i]][,3:5]),3:5,i] <- as.matrix(param_expert[[i]][,3:5])
  }
  
  #Stan does not allow NA
  expert.array[is.na(expert.array)] <- -999.2
  
  if(!is.null(times_expert)){
    n_time_expert <- length(times_expert)
    time_expert <- as.array(times_expert)
  }else{
    n_time_expert <- 1
    time_expert <- numeric(0) #This produces an array of size 0
    #https://dev.to/martinmodrak/optional-parametersdata-in-stan-4o33
    
    if (distr3 %in% c("gam", "gga", "gom")){
      time_expert <- 1 # Has to be defined for JAGS
      
    }
    
    
  }
  
  expert.array
  
}


tmpfun <- get("logLikFactory", envir = asNamespace("flexsurv"))
environment(logLikFactory) <- environment(tmpfun)
attributes(logLikFactory) <- attributes(tmpfun)
assignInNamespace("logLikFactory", logLikFactory, ns="flexsurv")


tmpfun <- get("flexsurvreg", envir = asNamespace("flexsurv"))
environment(flexsurvreg) <- environment(tmpfun)
attributes(flexsurvreg) <- attributes(tmpfun)
assignInNamespace("flexsurvreg", flexsurvreg, ns="flexsurv")

tmpfun <- get("minusloglik.flexsurv", envir = asNamespace("flexsurv"))
environment(minusloglik.flexsurv) <- environment(tmpfun)
attributes(minusloglik.flexsurv) <- attributes(tmpfun)
assignInNamespace("minusloglik.flexsurv", minusloglik.flexsurv, ns="flexsurv")
