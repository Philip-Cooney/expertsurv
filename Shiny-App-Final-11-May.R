
#Fixes required for expertsurv implemented in this file!!
#Fixes - fix error_mod_normal
#Fix DLSsurvspline was accessed in the wrong order with test.deriv
#Fix this mess in runMLE where id_trt and id_St don't work -- Fixed - pargSurv/expert wasn't reading the correct ids.
#Fix the lik_lno function the likelihood was incorrect for the compute ICS_stan



## Load Packages
list.of.packages <- need<-c("expertsurv", "dplyr", "shiny", "ggplot2", "ggfortify", "shinyWidgets", "survminer","shinycssloaders", "shinyjs", "shinyMatrix", "readxl") #needed libraries

res <- lapply(list.of.packages, require, character.only = TRUE)
not.loaded <-   list.of.packages[which(sapply(res, unlist) ==F)]
not.loaded2 <- lapply(not.loaded, require, character.only = TRUE)
not.installed <-   not.loaded2[which(sapply(not.loaded2, unlist) ==F)]
#load the packages
if(length(not.installed)) install.packages(not.loaded)

m_default_gen <- function(){ #Might add arguments to this function
  
  m_default <- matrix(nrow = 3, ncol = 2)
  colnames(m_default) <- c("Cum Prob", "Expert_1")
  #rownames(m_default) <- rep("Expert_1",3)
  m_default[,1] <- c(0.025, 0.5, 0.975)
  return(m_default)
}



#######################################################################################
return_pooled_info <- function(input_mat, St_indic = 1,dist = "best", mode =NULL){
  #dist_considered <- c("normal","t","gamma", "lognormal", "beta") 
  
  if(St_indic == 1){
    lower_bound = 0
    upper_bound = 1
  }else{
    lower_bound = -Inf
    upper_bound = Inf
  }
  
  fit.eval <- expertsurv:::fitdist_mod(input_mat[,2:ncol(input_mat), drop = F],
                                       probs = input_mat[,1], upper = upper_bound, lower = lower_bound, 
                                       expertnames = paste0("Expert_",1:(ncol(input_mat)-1)),
                                       mode = mode)
  
  plts_pool <- expertsurv:::makePoolPlot(fit= fit.eval,
                                         xl =lower_bound,
                                         xu =upper_bound,
                                         d = dist,
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
                                         expertnames = paste0("Expert_",1:(ncol(input_mat)-1)),
                                         St_indic =St_indic)
  
  dfs_pool <-  plts_pool[["data"]]
  if(dist == "best"){
    selc_fit <- fit.eval$best.fitting[,"best.fit"]
  }else{
    selc_fit <- rep(dist, length(fit.eval$best.fitting[,"best.fit"]))
  }
  selc_fit_loc <- sapply(selc_fit, function(x){which(x  == names(fit.eval$ssq))})
  
  pool.df_output <- matrix(nrow = length(selc_fit),ncol = 3)
  colnames(pool.df_output) <- c("param1", "param2", "param3")
  
  for(j in 1:length(selc_fit_loc)){
    pool.df_output[j,1:length(fit.eval[[selc_fit_loc[j]]][j,])] <-  as.numeric(as.vector(fit.eval[[selc_fit_loc[j]]][j,]))
  }
  dfs_expert <- data.frame(dist = names(selc_fit_loc), wi = 1/nrow(pool.df_output), pool.df_output)
  
  return(list(dfs_expert, plts_pool))
}

normal.error_mod <- function (parameters, values, probabilities, weights, mode){
  res1 <- sum(weights * (pnorm(values, parameters[1], exp(parameters[2])) - 
                           probabilities)^2)
  if (!is.null(mode)) {
    res1 <- res1 + (parameters[1] - mode)^2
  }
  return(res1)
}

t.error_mod <- function (parameters, values, probabilities, weights, degreesfreedom, 
                         mode){ 
  res1 <- sum(weights * (stats::pt((values - parameters[1])/exp(parameters[2]), 
                                   degreesfreedom) - probabilities)^2)
  if (!is.null(mode)) {
    res1 <- res1 + (parameters[1] - mode)^2
  }
  return(res1)
}
tmpfun <- get("t.error_mod", envir = asNamespace("expertsurv"))
environment(t.error_mod) <- environment(tmpfun)
attributes(t.error_mod) <- attributes(tmpfun)
assignInNamespace("t.error_mod", t.error_mod, ns="expertsurv")   


tmpfun <- get("normal.error_mod", envir = asNamespace("expertsurv"))
environment(normal.error_mod) <- environment(tmpfun)
attributes(normal.error_mod) <- attributes(tmpfun)
assignInNamespace("normal.error_mod", normal.error_mod, ns="expertsurv")   


gamma.error_mod <- function (parameters, values, probabilities, weights, mode){
  res1 <- sum(weights * (stats::pgamma(values, exp(parameters[1]), 
                                       exp(parameters[2])) - probabilities)^2)
  if (!is.null(mode)) {
    res1 <- res1 + ((exp(parameters[1]) - 1)/exp(parameters[2]) - 
                      mode)^2
  }
  return(res1)
}

tmpfun <- get("gamma.error_mod", envir = asNamespace("expertsurv"))
environment(gamma.error_mod) <- environment(tmpfun)
attributes(gamma.error_mod) <- attributes(tmpfun)
assignInNamespace("gamma.error_mod", gamma.error_mod, ns="expertsurv")   

lognormal.error_mod <-function (parameters, values, probabilities, weights, mode){
  res1 <- sum(weights * (stats::plnorm(values, parameters[1], 
                                       exp(parameters[2])) - probabilities)^2)
  if (!is.null(mode)) {
    res1 <- res1 + (exp(parameters[1] - exp(parameters[2])^2) - 
                      mode)^2
  }
  return(res1)
}
tmpfun <- get("lognormal.error_mod", envir = asNamespace("expertsurv"))
environment(lognormal.error_mod) <- environment(tmpfun)
attributes(lognormal.error_mod) <- attributes(tmpfun)
assignInNamespace("lognormal.error_mod", lognormal.error_mod, ns="expertsurv")   


beta.error_mod <- function (parameters, values, probabilities, weights, mode){
  res1 <- sum(weights * (stats::pbeta(values, exp(parameters[1]), 
                                      exp(parameters[2])) - probabilities)^2)
  if (!is.null(mode)) {
    res1 <- res1 + ((exp(parameters[1]) - 1)/(exp(parameters[1]) + 
                                                exp(parameters[2]) - 2) - mode)^2
  }
  return(res1)
}

tmpfun <- get("beta.error_mod", envir = asNamespace("expertsurv"))
environment(beta.error_mod) <- environment(tmpfun)
attributes(beta.error_mod) <- attributes(tmpfun)
assignInNamespace("beta.error_mod", beta.error_mod, ns="expertsurv")   


fitdist_mod <- function (vals, probs, lower = -Inf, upper = Inf, weights = 1, 
          tdf = 3, expertnames = NULL, excludelog.mirror = TRUE, mode = NULL){
  logt.error <- utils::getFromNamespace("logt.error", "SHELF")
  gamma.error <- utils::getFromNamespace("gamma.error", "SHELF")
  lognormal.error <- utils::getFromNamespace("lognormal.error", 
                                             "SHELF")
  logt.error <- utils::getFromNamespace("logt.error", "SHELF")
  makeGroupPlot <- utils::getFromNamespace("makeGroupPlot", 
                                           "SHELF")
  makeLinearPoolPlot <- utils::getFromNamespace("makeLinearPoolPlot", 
                                                "SHELF")
  makeSingleExpertPlot <- utils::getFromNamespace("makeSingleExpertPlot", 
                                                  "SHELF")
  expertdensity <- utils::getFromNamespace("expertdensity", 
                                           "SHELF")
  if (is.matrix(vals) == F) {
    vals <- matrix(vals, nrow = length(vals), ncol = 1)
  }
  if (is.matrix(probs) == F) {
    probs <- matrix(probs, nrow = nrow(vals), ncol = ncol(vals))
  }
  if (is.matrix(weights) == F) {
    weights <- matrix(weights, nrow = nrow(vals), ncol = ncol(vals))
  }
  if (length(lower) == 1) {
    lower <- rep(lower, ncol(vals))
  }
  if (length(upper) == 1) {
    upper <- rep(upper, ncol(vals))
  }
  if (length(tdf) == 1) {
    tdf <- rep(tdf, ncol(vals))
  }
  n.experts <- ncol(vals)
  normal.parameters <- matrix(NA, n.experts, 2)
  t.parameters <- matrix(NA, n.experts, 3)
  mirrorgamma.parameters <- gamma.parameters <- matrix(NA, 
                                                       n.experts, 2)
  mirrorlognormal.parameters <- lognormal.parameters <- matrix(NA, 
                                                               n.experts, 2)
  mirrorlogt.parameters <- logt.parameters <- matrix(NA, n.experts, 
                                                     3)
  beta.parameters <- matrix(NA, n.experts, 2)
  ssq <- matrix(NA, n.experts, 9)
  colnames(ssq) <- c("normal", "t", "gamma", "lognormal", "logt", 
                     "beta", "mirrorgamma", "mirrorlognormal", "mirrorlogt")
  if (n.experts > 1 & n.experts < 27 & is.null(expertnames)) {
    expertnames <- paste("expert.", LETTERS[1:n.experts], 
                         sep = "")
  }
  if (n.experts > 27 & is.null(expertnames)) {
    expertnames <- paste("expert.", 1:n.experts, sep = "")
  }
  limits <- data.frame(lower = lower, upper = upper)
  row.names(limits) <- expertnames
  for (i in 1:n.experts) {
   # if (min(probs[, i]) > 0.4) {
   #  stop("smallest elicited probability must be less than 0.4")
   # }
    if (min(probs[, i]) < 0 | max(probs[, i]) > 1) {
      stop("probabilities must be between 0 and 1")
    }
  #  if (max(probs[, i]) < 0.6) {
  #    stop("largest elicited probability must be greater than 0.6")
  #  }
    if (min(vals[, i]) < lower[i]) {
      stop("elicited parameter values cannot be smaller than lower parameter limit")
    }
    if (max(vals[, i]) > upper[i]) {
      stop("elicited parameter values cannot be greater than upper parameter limit")
    }
    if (tdf[i] <= 0) {
      stop("Student-t degrees of freedom must be greater than 0")
    }
    if (min(probs[-1, i] - probs[-nrow(probs), i]) < 0) {
      stop("probabilities must be specified in ascending order")
    }
    if (min(vals[-1, i] - vals[-nrow(vals), i]) <= 0) {
      stop("parameter values must be specified in ascending order")
    }
    inc <- (probs[, i] > 0) & (probs[, i] < 1)
    minprob <- min(probs[inc, i])
    maxprob <- max(probs[inc, i])
    minvals <- min(vals[inc, i])
    maxvals <- max(vals[inc, i])
    q.fit <- stats::approx(x = probs[inc, i], y = vals[inc, 
                                                       i], xout = c(0.4, 0.5, 0.6))$y
    l <- q.fit[1]
    u <- q.fit[3]
    minq <- stats::qnorm(minprob)
    maxq <- stats::qnorm(maxprob)
    m <- (minvals * maxq - maxvals * minq)/(maxq - minq)
    v <- ((maxvals - minvals)/(maxq - minq))^2
    normal.fit <- stats::optim(c(m, 0.5 * log(v)), normal.error_mod, 
                               values = vals[inc, i], probabilities = probs[inc, 
                                                                            i], weights = weights[inc, i], mode = mode[i])
    normal.parameters[i, ] <- c(normal.fit$par[1], exp(normal.fit$par[2]))
    ssq[i, "normal"] <- normal.fit$value
    t.fit <- stats::optim(c(m, 0.5 * log(v)), t.error_mod, 
                          values = vals[inc, i], probabilities = probs[inc, 
                                                                       i], weights = weights[inc, i], degreesfreedom = tdf[i], 
                          mode = mode[i])
    t.parameters[i, 1:2] <- c(t.fit$par[1], exp(t.fit$par[2]))
    t.parameters[i, 3] <- tdf[i]
    ssq[i, "t"] <- t.fit$value
    if (lower[i] > -Inf) {
      vals.scaled1 <- vals[inc, i] - lower[i]
      m.scaled1 <- m - lower[i]
      gamma.fit <- stats::optim(c(log(m.scaled1^2/v), log(m.scaled1/v)), 
                                gamma.error_mod, values = vals.scaled1, probabilities = probs[inc, 
                                                                                              i], weights = weights[inc, i], mode = mode[i])
      gamma.parameters[i, ] <- exp(gamma.fit$par)
      ssq[i, "gamma"] <- gamma.fit$value
      std <- ((log(u - lower[i]) - log(l - lower[i]))/1.35)
      mlog <- (log(minvals - lower[i]) * maxq - log(maxvals - 
                                                      lower[i]) * minq)/(maxq - minq)
      lognormal.fit <- stats::optim(c(mlog, log(std)), 
                                    lognormal.error_mod, values = vals.scaled1, probabilities = probs[inc, 
                                                                                                      i], weights = weights[inc, i], mode = mode[i])
      lognormal.parameters[i, 1:2] <- c(lognormal.fit$par[1], 
                                        exp(lognormal.fit$par[2]))
      ssq[i, "lognormal"] <- lognormal.fit$value
      logt.fit <- stats::optim(c(log(m.scaled1), log(std)), 
                               logt.error, values = vals.scaled1, probabilities = probs[inc, 
                                                                                        i], weights = weights[inc, i], degreesfreedom = tdf[i])
      logt.parameters[i, 1:2] <- c(logt.fit$par[1], exp(logt.fit$par[2]))
      logt.parameters[i, 3] <- tdf[i]
      ssq[i, "logt"] <- logt.fit$value
    }
    if ((lower[i] > -Inf) & (upper[i] < Inf)) {
      vals.scaled2 <- (vals[inc, i] - lower[i])/(upper[i] - 
                                                   lower[i])
      m.scaled2 <- (m - lower[i])/(upper[i] - lower[i])
      v.scaled2 <- v/(upper[i] - lower[i])^2
      alp <- abs(m.scaled2^3/v.scaled2 * (1/m.scaled2 - 
                                            1) - m.scaled2)
      bet <- abs(alp/m.scaled2 - alp)
      if (identical(probs[inc, i], (vals[inc, i] - lower[i])/(upper[i] - 
                                                              lower[i]))) {
        alp <- bet <- 1
      }
      beta.fit <- stats::optim(c(log(alp), log(bet)), beta.error_mod, 
                               values = vals.scaled2, probabilities = probs[inc, 
                                                                            i], weights = weights[inc, i], mode = mode[i])
      beta.parameters[i, ] <- exp(beta.fit$par)
      ssq[i, "beta"] <- beta.fit$value
    }
    if (upper[i] < Inf) {
      valsMirrored <- upper[i] - vals[inc, i]
      probsMirrored <- 1 - probs[inc, i]
      mMirrored <- upper[i] - m
      mirrorgamma.fit <- stats::optim(c(log(mMirrored^2/v), 
                                        log(mMirrored/v)), gamma.error, values = valsMirrored, 
                                      probabilities = probsMirrored, weights = weights[inc, 
                                                                                       i])
      mirrorgamma.parameters[i, ] <- exp(mirrorgamma.fit$par)
      ssq[i, "mirrorgamma"] <- mirrorgamma.fit$value
      mlogMirror <- (log(upper[i] - maxvals) * (1 - minq) - 
                       log(upper[i] - minvals) * (1 - maxq))/(maxq - 
                                                                minq)
      stdMirror <- ((log(upper[i] - l) - log(upper[i] - 
                                               u))/1.35)
      mirrorlognormal.fit <- optim(c(mlogMirror, log(stdMirror)), 
                                   lognormal.error, values = valsMirrored, probabilities = probsMirrored, 
                                   weights = weights[inc, i])
      mirrorlognormal.parameters[i, 1:2] <- c(mirrorlognormal.fit$par[1], 
                                              exp(mirrorlognormal.fit$par[2]))
      ssq[i, "mirrorlognormal"] <- mirrorlognormal.fit$value
      mirrorlogt.fit <- stats::optim(c(log(mMirrored), 
                                       log(stdMirror)), logt.error, values = valsMirrored, 
                                     probabilities = probsMirrored, weights = weights[inc, 
                                                                                      i], degreesfreedom = tdf[i])
      mirrorlogt.parameters[i, 1:2] <- c(mirrorlogt.fit$par[1], 
                                         exp(mirrorlogt.fit$par[2]))
      mirrorlogt.parameters[i, 3] <- tdf[i]
      ssq[i, "mirrorlogt"] <- mirrorlogt.fit$value
    }
  }
  dfn <- data.frame(normal.parameters)
  names(dfn) <- c("mean", "sd")
  row.names(dfn) <- expertnames
  dft <- data.frame(t.parameters)
  names(dft) <- c("location", "scale", "df")
  row.names(dft) <- expertnames
  dfg <- data.frame(gamma.parameters)
  names(dfg) <- c("shape", "rate")
  row.names(dfg) <- expertnames
  dfmirrorg <- data.frame(mirrorgamma.parameters)
  names(dfmirrorg) <- c("shape", "rate")
  row.names(dfmirrorg) <- expertnames
  dfln <- data.frame(lognormal.parameters)
  names(dfln) <- c("mean.log.X", "sd.log.X")
  row.names(dfln) <- expertnames
  dfmirrorln <- data.frame(mirrorlognormal.parameters)
  names(dfmirrorln) <- c("mean.log.X", "sd.log.X")
  row.names(dfmirrorln) <- expertnames
  dflt <- data.frame(logt.parameters)
  names(dflt) <- c("location.log.X", "scale.log.X", "df.log.X")
  row.names(dflt) <- expertnames
  dfmirrorlt <- data.frame(mirrorlogt.parameters)
  names(dfmirrorlt) <- c("location.log.X", "scale.log.X", "df.log.X")
  row.names(dfmirrorlt) <- expertnames
  dfb <- data.frame(beta.parameters)
  names(dfb) <- c("shape1", "shape2")
  row.names(dfb) <- expertnames
  ssq <- data.frame(ssq)
  row.names(ssq) <- expertnames
  if (excludelog.mirror) {
    reducedssq <- ssq[, c("normal", "t", "gamma", "lognormal", 
                          "beta")]
    index <- apply(reducedssq, 1, which.min)
    best.fitting <- data.frame(best.fit = names(reducedssq)[index])
  }
  else {
    index <- apply(ssq, 1, which.min)
    best.fitting <- data.frame(best.fit = names(ssq)[index])
  }
  row.names(best.fitting) <- expertnames
  vals <- data.frame(vals)
  names(vals) <- expertnames
  probs <- data.frame(probs)
  names(probs) <- expertnames
  fit <- list(Normal = dfn, Student.t = dft, Gamma = dfg, Log.normal = dfln, 
              Log.Student.t = dflt, Beta = dfb, mirrorgamma = dfmirrorg, 
              mirrorlognormal = dfmirrorln, mirrorlogt = dfmirrorlt, 
              ssq = ssq, best.fitting = best.fitting, vals = t(vals), 
              probs = t(probs), limits = limits)
  class(fit) <- "elicitation"
  fit
}
tmpfun <- get("fitdist_mod", envir = asNamespace("expertsurv"))
environment(fitdist_mod) <- environment(tmpfun)
attributes(fitdist_mod) <- attributes(tmpfun)
assignInNamespace("fitdist_mod", fitdist_mod, ns="expertsurv") 



logLikFactory <- function (Y, X = 0, weights, bhazard, rtrunc, dlist, inits, 
                           dfns, aux, mx, fixedpars = NULL, expert_opinion){
  
  
  #Fix need DLSsurvspline
  
  
  buildTransformer <- utils::getFromNamespace("buildTransformer", 
                                              "flexsurv")
  buildAuxParms <- utils::getFromNamespace("buildAuxParms", 
                                           "flexsurv")
  check.flexsurv.response <- utils::getFromNamespace("check.flexsurv.response", 
                                                     "flexsurv")
  tsfn <- utils::getFromNamespace("tsfn", "flexsurv")
  DLdsurvspline <- utils::getFromNamespace("DLdsurvspline", 
                                           "flexsurv")
  flexsurv.splineinits.cox <- utils::getFromNamespace("flexsurv.splineinits.cox", 
                                                      "flexsurv")
  parse.dist <- utils::getFromNamespace("parse.dist", "flexsurv")
  form.dp <- utils::getFromNamespace("form.dp", "flexsurv")
  check.formula <- utils::getFromNamespace("check.formula", 
                                           "flexsurv")
  anc_from_formula <- utils::getFromNamespace("anc_from_formula", 
                                              "flexsurv")
  get.locform <- utils::getFromNamespace("get.locform", "flexsurv")
  concat.formulae <- utils::getFromNamespace("concat.formulae", 
                                             "flexsurv")
  compress.model.matrices <- utils::getFromNamespace("compress.model.matrices", 
                                                     "flexsurv")
  expand.inits.args <- utils::getFromNamespace("expand.inits.args", 
                                               "flexsurv")
  .hess_to_cov <- utils::getFromNamespace(".hess_to_cov", 
                                          "flexsurv")
  logh <- utils::getFromNamespace("logh", "flexsurv")
  dexph <- utils::getFromNamespace("dexph", "flexsurv")
  deriv.test <- utils::getFromNamespace("deriv.test", "flexsurv")
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
      offseti <- log(1 + bhazard[event]/exp(loghaz) * 
                       weights[event])
    }
    else {
      offseti <- default.offset
    }
    loglik[event] <- (logdens * event.weights) + offseti
    if (any(!event)) 
      loglik[!event] <- (log(pmax - pmin) * no.event.weights)
    loglik <- loglik - log(pobs) * weights
    pargs.expert <- fnargs
    
    if (!is.null(expert_opinion)) {
      #browser()
      if (expert_opinion$St_indic == 1) {
        
        pargs.expert <- lapply(pargs.expert, function(x) {
          x[expert_opinion$id_St]
        })
        
        if(!is.null(pargs.expert$knots)){
          pargs.expert[["knots"]] <- fnargs$knots
          pargs.expert[["scale"]] <- fnargs$scale
          pargs.expert[["timescale"]] <- fnargs$timescale
        }
          
        pargs.expert$q <- expert_opinion$times
        psurv_expert <- 1 - do.call(dfns$p, pargs.expert)
      }
      else{
        pargs.rmst <- pargs.expert
        pargs.rmst$q <- NULL
        pargs.rmst <- lapply(pargs.rmst, function(x) {
          x[c(expert_opinion$id_comp, expert_opinion$id_trt)]
        })
        
        if(!is.null(pargs.expert$knots)){
          pargs.expert[["knots"]] <- fnargs$knots
          pargs.expert[["scale"]] <- fnargs$scale
          pargs.expert[["timescale"]] <- fnargs$timescale
        }
        psurv_expert <- diff(do.call(dfns$mean, pargs.rmst))
      }
      LL_expert <- rep(NA, length(expert_opinion$times))
      for (q in 1:length(expert_opinion$times)) {
        param_expert_mat <- abind::adrop(expert_opinion$param_expert[, 
                                                                     , q, drop = F], drop = 3)
        
        #Will need to fix for CRAN!!!
        LL_expert[q] <- expertsurv:::expert_log_dens(psurv_expert[q], 
                                        df = param_expert_mat, expert_opinion$pool, 
                                        k_norm = expert_opinion$k_norm[q], St_indic = expert_opinion$St_indic)
      }
      ret <- -sum(loglik, LL_expert)
    }
    else {
      ret <- -sum(loglik)
    }
    attr(ret, "indiv") <- loglik
    ret
  }
}

tmpfun <- get("logLikFactory", envir = asNamespace("expertsurv"))
environment(logLikFactory) <- environment(tmpfun)
attributes(logLikFactory) <- attributes(tmpfun)
assignInNamespace("logLikFactory", logLikFactory, ns="expertsurv")




runMLE <- function (x, exArgs){
  formula <- exArgs$formula
  data = exArgs$data
  availables <- load_availables()
  d3 <- manipulate_distributions(x)$distr3
  x <- manipulate_distributions(x)$distr
  expert_opinion_flex <- list()
  
  
  
  formula_temp <- update(formula, paste(all.vars(formula, 
                                                 data)[1], "~", all.vars(formula, data)[2], "+."))
  mf <- tibble::as_tibble(model.frame(formula_temp, data)) %>% 
    dplyr::rename(time = 1, event = 2) %>% rename_if(is.factor, 
                                                     .funs = ~gsub("as.factor[( )]", "", .x)) %>% dplyr::rename_if(is.factor, 
                                                                                                                   .funs = ~gsub("[( )]", "", .x)) %>% dplyr::bind_cols(tibble::as_tibble(stats::model.matrix(formula_temp, 
                                                                                                                                                                                                              data)) %>% dplyr::select(contains("Intercept"))) %>% 
    dplyr::select(time, event, contains("Intercept"), 
                  everything()) %>% tibble::rownames_to_column("ID")
  
  
  if (ncol(mf) == 4) {
    expert_opinion_flex$id_St <- 1
  }
  else if (ncol(mf) == 5) {
    if (exArgs$opinion_type == "survival") {
      expert_opinion_flex$id_St <- min(which(mf[, 5] == exArgs$id_St))
    }
    else {
      expert_opinion_flex$id_trt <- min(which(mf[, 5] == exArgs$id_trt))
      if (length(unique(mf[, 5] %>% pull())) == 2) {
        expert_opinion_flex$id_comp <- min(which(mf[, 5] != exArgs$id_trt))
      }
      else {
        expert_opinion_flex$id_comp <- min(which(mf[, 5] == exArgs$id_comp))
      }
    }
  }
  else {
    message("We do not allow more than one covariate (i.e. treatment) in the analysis - although it is technically possible")
    stop()
  }
  
  
  if (exArgs$opinion_type != "survival") {
    times <- 99999
    expert_opinion_flex$St_indic <- 0
  }
  
  if (exArgs$opinion_type == "survival") {
    expert_opinion_flex$St_indic <- 1
    times <- exArgs$times_expert
  }
  expert_opinion_flex$param_expert <- make_data_expert(exArgs$param_expert, 
                                                       times)
  expert_opinion_flex$times <- times
  if (exArgs$pool_type == "linear pool") {
    expert_opinion_flex$pool <- 1
  }
  else {
    expert_opinion_flex$pool <- 0
  }
  if (expert_opinion_flex$pool == 0) {
    expert_opinion_flex$k_norm <- get_k_norm(exArgs$param_expert)
  }
  else {
    expert_opinion_flex$k_norm <- NULL
  }
  if (exists("method_mle", where = exArgs)) {
    method_mle <- exArgs$method_mle
  }
  else {
    method_mle <- "BFGS"
  }
  tic <- proc.time()
  if (x == "survspline") {
    if (exists("bhazard", where = exArgs)) {
      bhazard <- exArgs$bhazard
    }
    else {
      bhazard <- NULL
    }
    if (exists("weights", where = exArgs)) {
      weights <- exArgs$weights
    }
    else {
      weights <- NULL
    }
    if (exists("subset", where = exArgs)) {
      subset <- exArgs$subset
    }
    else {
      subset <- NULL
    }
    if (exists("knots", where = exArgs)) {
      knots <- exArgs$knots
    }
    else {
      knots <- NULL
    }
    if (exists("k", where = exArgs)) {
      k <- exArgs$k
    }
    else {
      k <- 0
    }
    if (exists("bknots", where = exArgs)) {
      bknots <- exArgs$bknots
    }
    else {
      bknots <- NULL
    }
    if (exists("scale", where = exArgs)) {
      scale <- exArgs$scale
    }
    else {
      scale <- "hazard"
    }
    if (exists("timescale", where = exArgs)) {
      timescale <- exArgs$scale
    }
    else {
      timescale <- "log"
    }
    suppressWarnings({
      #browser()
      model_mle <- flexsurvspline(formula = formula, data = data, 
                                  k = k, knots = knots, bknots = bknots, scale = scale, 
                                  timescale = timescale, expert_opinion = NULL, 
                                  method = method_mle)
      model <- flexsurvspline(formula = formula, data = data, 
                              k = k, knots = knots, bknots = bknots, scale = scale, 
                              timescale = timescale, expert_opinion = expert_opinion_flex, 
                              method = method_mle, inits = model_mle$res[,1])
    })
  }
  else {
    suppressWarnings({
      model_mle <- flexsurvreg(formula = formula, data = data, 
                               dist = x, expert_opinion = NULL, method = method_mle)
      model <- flexsurvreg(formula = formula, data = data, 
                           dist = x, expert_opinion = expert_opinion_flex, 
                           method = method_mle, inits = model_mle$res[,1])
    })
  }
  toc <- proc.time() - tic
  model_name <- d3
  list(model = model, aic = model$AIC, bic = -2 * model$loglik + 
         model$npars * log(model$N), dic = NULL, time2run = toc[3], 
       model_name = model_name)
}


tmpfun <- get("runMLE", envir = asNamespace("expertsurv"))
environment(runMLE) <- environment(tmpfun)
attributes(runMLE) <- attributes(tmpfun)
assignInNamespace("runMLE", runMLE, ns="expertsurv")



flexsurvspline <-function (formula, data, weights, bhazard, rtrunc, subset, k = 0, 
          knots = NULL, bknots = NULL, scale = "hazard", timescale = "log", 
          expert_opinion = NULL, ...){
  buildTransformer <- utils::getFromNamespace("buildTransformer", 
                                              "flexsurv")
  buildAuxParms <- utils::getFromNamespace("buildAuxParms", 
                                           "flexsurv")
  check.flexsurv.response <- utils::getFromNamespace("check.flexsurv.response", 
                                                     "flexsurv")
  tsfn <- utils::getFromNamespace("tsfn", "flexsurv")
  DLdsurvspline <- utils::getFromNamespace("DLdsurvspline", 
                                           "flexsurv")
  flexsurv.splineinits.cox <- utils::getFromNamespace("flexsurv.splineinits.cox", 
                                                      "flexsurv")
  flexsurv.splineinits <- utils::getFromNamespace("flexsurv.splineinits", 
                                                  "flexsurv")
  parse.dist <- utils::getFromNamespace("parse.dist", "flexsurv")
  form.dp <- utils::getFromNamespace("form.dp", "flexsurv")
  check.formula <- utils::getFromNamespace("check.formula", 
                                           "flexsurv")
  anc_from_formula <- utils::getFromNamespace("anc_from_formula", 
                                              "flexsurv")
  get.locform <- utils::getFromNamespace("get.locform", "flexsurv")
  concat.formulae <- utils::getFromNamespace("concat.formulae", 
                                             "flexsurv")
  compress.model.matrices <- utils::getFromNamespace("compress.model.matrices", 
                                                     "flexsurv")
  expand.inits.args <- utils::getFromNamespace("expand.inits.args", 
                                               "flexsurv")
  .hess_to_cov <- utils::getFromNamespace(".hess_to_cov", "flexsurv")
  logh <- utils::getFromNamespace("logh", "flexsurv")
  dexph <- utils::getFromNamespace("dexph", "flexsurv")
  deriv.test <- utils::getFromNamespace("deriv.test", "flexsurv")
  DLSsurvspline <- utils::getFromNamespace("DLSsurvspline", "flexsurv")
  call <- match.call()
  indx <- match(c("formula", "data", "weights", "bhazard", 
                  "rtrunc", "subset", "na.action"), names(call), nomatch = 0)
  if (indx[1] == 0) 
    stop("A \"formula\" argument is required")
  temp <- call[c(1, indx)]
  temp[[1]] <- as.name("model.frame")
  f2 <- stats::as.formula(gsub("(.*~).*", "\\1 1", Reduce(paste, 
                                                          deparse(formula))))
  environment(f2) <- environment(formula)
  temp[["formula"]] <- f2
  if (missing(data)) 
    temp[["data"]] <- environment(formula)
  if (missing(data)) 
    data <- environment(formula)
  m <- eval(temp, parent.frame())
  Y <- check.flexsurv.response(stats::model.extract(m, "response"))
  dtimes <- Y[, "stop"][Y[, "status"] == 1]
  if (is.null(knots)) {
    is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) abs(x - 
                                                                       round(x)) < tol
    if (is.null(k)) 
      stop("either \"knots\" or \"k\" must be specified")
    if (!is.numeric(k)) 
      stop("k must be numeric")
    if (!is.wholenumber(k) || (k < 0)) 
      stop("number of knots \"k\" must be a non-negative integer")
    knots <- stats::quantile(tsfn(dtimes, timescale), seq(0, 
                                                          1, length = k + 2)[-c(1, k + 2)])
  }
  else {
    if (!is.numeric(knots)) 
      stop("\"knots\" must be a numeric vector")
    minlogtime <- min(tsfn(Y[, "stop"], timescale))
    if (any(knots <= minlogtime)) {
      badknots <- knots[knots < min(tsfn(Y[, "stop"], timescale))]
      plural <- if (length(badknots) > 1) 
        "s"
      else ""
      stop(sprintf("knot%s %s less than or equal to minimum %stime", 
                   plural, paste(badknots, collapse = ", "), (if (timescale == 
                                                                  "log") 
                     "log "
                     else "")))
    }
    maxlogtime <- max(tsfn(Y[, "stop"], timescale))
    if (any(knots >= maxlogtime)) {
      badknots <- knots[knots > maxlogtime]
      plural <- if (length(badknots) > 1) 
        "s"
      else ""
      stop(sprintf("knot%s %s greater than or equal to maximum %stime", 
                   plural, paste(badknots, collapse = ", "), (if (timescale == 
                                                                  "log") 
                     "log "
                     else "")))
    }
  }
  if (is.null(bknots)) {
    if (length(dtimes) > 0) {
      bt <- dtimes
    }
    else {
      bt <- c(Y[, "time1"], Y[, "time2"], Y[, "time"])
      bt <- bt[is.finite(bt)]
    }
    bknots <- c(min(tsfn(bt, timescale)), max(tsfn(bt, timescale)))
    if (bknots[1] == bknots[2]) 
      warning("minimum and maximum log death times are the same: knot and boundary knot locations should be supplied explicitly")
  }
  else if (!is.numeric(bknots) || (length(bknots) != 2)) 
    stop("boundary knots should be a numeric vector of length 2")
  knots <- c(bknots[1], knots, bknots[2])
  nk <- length(knots)
  custom.fss <- list(name = "survspline", pars = c(paste0("gamma", 
                                                          0:(nk - 1))), location = c("gamma0"), transforms = rep(c(identity), 
                                                                                                                 nk), inv.transforms = rep(c(identity), nk), inits = flexsurv.splineinits)
  aux <- list(knots = knots, scale = scale, timescale = timescale)
  dfn <- flexsurv::unroll.function(flexsurv::dsurvspline, gamma = 0:(nk - 
                                                                       1))
  pfn <- flexsurv::unroll.function(flexsurv::psurvspline, gamma = 0:(nk - 
                                                                       1))
  rfn <- flexsurv::unroll.function(flexsurv::rsurvspline, gamma = 0:(nk - 
                                                                       1))
  hfn <- flexsurv::unroll.function(flexsurv::hsurvspline, gamma = 0:(nk - 
                                                                       1))
  Hfn <- flexsurv::unroll.function(flexsurv::Hsurvspline, gamma = 0:(nk - 
                                                                       1))
  qfn <- flexsurv::unroll.function(flexsurv::qsurvspline, gamma = 0:(nk - 
                                                                       1))
  meanfn <- flexsurv::unroll.function(flexsurv::mean_survspline, 
                                      gamma = 0:(nk - 1))
  rmstfn <- flexsurv::unroll.function(flexsurv::rmst_survspline, 
                                      gamma = 0:(nk - 1))
  Ddfn <- if (scale == "normal") 
    NULL
  else flexsurv::unroll.function(DLdsurvspline, gamma = 0:(nk - 
                                                             1))
  DSfn <- if (scale == "normal") 
    NULL
  else flexsurv::unroll.function(DLSsurvspline, gamma = 0:(nk - 
                                                             1))
  args <- c(list(formula = formula, data = data, dist = custom.fss, 
                 dfns = list(d = dfn, p = pfn, r = rfn, h = hfn, H = Hfn, 
                             rmst = rmstfn, mean = meanfn, q = qfn, DLd = Ddfn, 
                             DLS = DSfn, deriv = !(scale == "normal")), aux = aux), 
            list(...))
  fpold <- args$fixedpars
  args$fixedpars <- TRUE
  args$weights <- temp$weights
  args$bhazard <- temp$bhazard
  args$rtrunc <- temp$rtrunc
  args$subset <- temp$subset
  if (is.null(expert_opinion)) {
    expert_opinion <- vector(mode = "list", length = 1)
    names(expert_opinion) <- "expert_opinion"
    args <- append(args, expert_opinion)
  }
  else {
    args[["expert_opinion"]] <- expert_opinion
  }
  if (is.infinite(do.call("flexsurvreg", args)$loglik)) {
    args$dist$inits <- flexsurv.splineinits.cox
  }
  args$fixedpars <- fpold
  ret <- do.call("flexsurvreg", args)
  ret <- c(ret, list(k = length(knots) - 2, knots = knots, 
                     scale = scale, timescale = timescale))
  ret$call <- call
  class(ret) <- "flexsurvreg"
  ret
}

tmpfun <- get("flexsurvspline", envir = asNamespace("expertsurv"))
environment(flexsurvspline) <- environment(tmpfun)
attributes(flexsurvspline) <- attributes(tmpfun)
assignInNamespace("flexsurvspline", flexsurvspline, ns="expertsurv")



lik_lno <- function (x, linpred, linpred.hat, model, data.stan){
  dist = "lno"
  sigma = as.numeric(rstan::extract(model)$alpha)
  sigma.hat = stats::median(sigma)
  logf <- matrix(unlist(lapply(1:nrow(linpred), function(i) {
    data.stan$d * log(flexsurv::hlnorm(data.stan$t, (linpred[i, 
    ]), sigma[i])) + log(1 - stats::plnorm(data.stan$t, 
                                           (linpred[i, ]), sigma[i]))
  })), nrow = nrow(linpred), byrow = T)
  logf.hat <- matrix(data.stan$d * log(flexsurv::hlnorm(data.stan$t, 
                                                        (linpred.hat), sigma.hat)) + log(1 - stats::plnorm(data.stan$t, 
                                                                                                           (linpred.hat), sigma.hat)), nrow = 1)
  logf.expert <- rep(NA, nrow(linpred))
  if (data.stan$St_indic == 1) {
    index_vec <- data.stan$id_St
  }
  else {
    index_vec <- c(data.stan$id_trt, data.stan$id_comp)
  }
  for (i in 1:nrow(linpred)) {
    logf.expert[i] <- expert_like(data.stan, dist_surv = dist, 
                                  param1 = linpred[i, index_vec],
                                  param2 = sigma[i])
  }
  logf.hat.expert <- expert_like(data.stan, dist_surv = dist,
                                 param1 = linpred.hat[1, index_vec],
                                 param2 = sigma.hat)
  npars <- 2 + sum(1 - apply(data.stan$X, 2, function(x) all(x == 
                                                               0)))
  list(logf = logf, logf.hat = logf.hat, npars = npars, f = NULL, 
       f.bar = NULL, s = NULL, s.bar = NULL, logf.expert = logf.expert, 
       logf.hat.expert = logf.hat.expert)
}

tmpfun <- get("lik_lno", envir = asNamespace("expertsurv"))
environment(lik_lno) <- environment(tmpfun)
attributes(lik_lno) <- attributes(tmpfun)
assignInNamespace("lik_lno", lik_lno, ns="expertsurv")





#Slight tweak in this function to take the upper limit from the highest mean value


makePoolPlot <- function (fit, xl, xu, d = "best", w = 1, lwd = 1, xlab = "x", 
          ylab = expression(f[X](x)), legend_full = TRUE, ql = NULL, 
          qu = NULL, nx = 500, addquantile = FALSE, fs = 12, expertnames = NULL, 
          St_indic){
  logt.error <- utils::getFromNamespace("logt.error", "SHELF")
  gamma.error <- utils::getFromNamespace("gamma.error", "SHELF")
  lognormal.error <- utils::getFromNamespace("lognormal.error", 
                                             "SHELF")
  logt.error <- utils::getFromNamespace("logt.error", "SHELF")
  makeGroupPlot <- utils::getFromNamespace("makeGroupPlot", 
                                           "SHELF")
  makeLinearPoolPlot <- utils::getFromNamespace("makeLinearPoolPlot", 
                                                "SHELF")
  makeSingleExpertPlot <- utils::getFromNamespace("makeSingleExpertPlot", 
                                                  "SHELF")
  expertdensity <- utils::getFromNamespace("expertdensity", 
                                           "SHELF")
  lpname <- c("linear pool", "log pool")
  expert <- ftype <- NULL
  n.experts <- nrow(fit$vals)
  if (length(d) == 1) {
    d <- rep(d, n.experts)
  }
  if (is.null(expertnames)) {
    if (n.experts < 27) {
      expertnames <- LETTERS[1:n.experts]
    }
    if (n.experts > 26) {
      expertnames <- 1:n.experts
    }
  }
  nxTotal <- nx + length(c(ql, qu))
  x <- matrix(0, nxTotal, n.experts)
  fx <- x
  if (min(w) < 0 | max(w) <= 0) {
    stop("expert weights must be non-negative, and at least one weight must be greater than 0.")
  }
  if (length(w) == 1) {
    w <- rep(w, n.experts)
  }
  weight <- matrix(w/sum(w), nxTotal, n.experts, byrow = T)
  sd.norm <- rep(NA, n.experts)
  for (i in 1:n.experts) {
  }
  if (is.infinite(xl) || is.infinite(xu)) {
    if (St_indic == 1) {
      xl <- 0
      xu <- 1
    }
    else {
      min.mean.index <- which.min(fit$Normal$mean)
      min.sd.index <- which.min(fit$Normal$sd)
      
      max.mean.index <- which.max(fit$Normal$mean)
      max.sd.index <- which.max(fit$Normal$sd)
      xl <- qnorm(0.001, fit$Normal[min.mean.index, 1], 
                  fit$Normal[min.sd.index, 2])
      xu <- qnorm(0.999, fit$Normal[max.mean.index, 1], 
                  fit$Normal[max.sd.index, 2])
    }
  }
  for (i in 1:n.experts) {
    densitydata <- expertdensity(fit, d[i], ex = i, xl, 
                                 xu, ql, qu, nx)
    x[, i] <- densitydata$x
    if (St_indic == 1) {
      k_trunc <- integrate.xy(x = x[, 1], fx = densitydata$fx)
    }
    else {
      k_trunc <- 1
    }
    fx[, i] <- densitydata$fx/k_trunc
  }
  fx.lp <- apply(fx * weight, 1, sum)
  if (any(is.infinite(fx^weight))) {
    warning("Print Non finite density for log pooling - Results invalid")
  }
  fx.logp <- apply(fx^weight, 1, prod)
  k_norm <- integrate.xy(x = x[, 1], fx = fx.logp)
  fx.logp <- fx.logp/k_norm
  df1 <- data.frame(x = rep(x[, 1], n.experts + 2), fx = c(as.numeric(fx), 
                                                           fx.lp, fx.logp), expert = factor(c(rep(expertnames, 
                                                                                                  each = nxTotal), rep("linear pool", nxTotal), rep("log pool", 
                                                                                                                                                    nxTotal)), levels = c(expertnames, "linear pool", "log pool")), 
                    ftype = factor(c(rep("individual", nxTotal * n.experts), 
                                     rep("linear pool", nxTotal), rep("log pool", nxTotal)), 
                                   levels = c("individual", "linear pool", "log pool")))
  df1$expert <- factor(df1$expert, levels = c(expertnames, 
                                              "linear pool", "log pool"))
  if (legend_full) {
    cols <- (scales::hue_pal())(n.experts + 2)
    linetypes <- c(rep("dashed", n.experts), "solid", "solid")
    sizes <- lwd * c(rep(0.5, n.experts), 1.5, 1.5)
    names(cols) <- names(linetypes) <- names(sizes) <- c(expertnames, 
                                                         lpname)
    p1 <- ggplot(df1, aes(x = x, y = fx, colour = expert, 
                          linetype = expert, size = expert)) + scale_colour_manual(values = cols, 
                                                                                   breaks = c(expertnames, lpname)) + scale_linetype_manual(values = linetypes, 
                                                                                                                                            breaks = c(expertnames, lpname)) + scale_size_manual(values = sizes, 
                                                                                                                                                                                                 breaks = c(expertnames, lpname))
  }
  else {
    p1 <- ggplot(df1, aes(x = x, y = fx, colour = ftype, 
                          linetype = ftype, size = ftype)) + scale_linetype_manual(name = "distribution", 
                                                                                   values = c("dashed", "solid", "solid")) + scale_size_manual(name = "distribution", 
                                                                                                                                               values = lwd * c(0.5, 1.5, 1.5)) + scale_color_manual(name = "distribution", 
                                                                                                                                                                                                     values = c("black", "red", "blue"))
  }
  if (legend_full) {
    for (i in 1:n.experts) {
      if (d[i] == "hist") {
        p1 <- p1 + geom_step(data = subset(df1, expert == 
                                             expertnames[i]), aes(colour = expert))
      }
      else {
        p1 <- p1 + geom_line(data = subset(df1, expert == 
                                             expertnames[i]), aes(colour = expert))
      }
    }
  }
  else {
    for (i in 1:n.experts) {
      if (d[i] == "hist") {
        p1 <- p1 + geom_step(data = subset(df1, expert == 
                                             expertnames[i]), aes(colour = ftype))
      }
      else {
        p1 <- p1 + geom_line(data = subset(df1, expert == 
                                             expertnames[i]), aes(colour = ftype))
      }
    }
  }
  if (length(unique(d)) == 1 & d[1] == "hist") {
    p1 <- p1 + geom_step(data = subset(df1, expert == lpname), 
                         aes(colour = expert))
  }
  else {
    p1 <- p1 + geom_line(data = subset(df1, expert == lpname[1]), 
                         aes(colour = expert))
    p1 <- p1 + geom_line(data = subset(df1, expert == lpname[2]), 
                         aes(colour = expert))
  }
  p1 <- p1 + labs(x = xlab, y = ylab)
  if ((!is.null(ql)) & (!is.null(qu)) & addquantile) {
    if (legend_full) {
      ribbon_col <- (scales::hue_pal())(n.experts + 2)[n.experts + 
                                                         2]
    }
    else {
      ribbon_col <- "red"
    }
    p1 <- p1 + geom_ribbon(data = with(df1, subset(df1, 
                                                   x <= ql & expert == lpname[1])), aes(ymax = fx, 
                                                                                        ymin = 0), alpha = 0.2, show.legend = FALSE, colour = NA, 
                           fill = ribbon_col) + geom_ribbon(data = with(df1, 
                                                                        subset(df1, x >= qu & expert == lpname[2])), aes(ymax = fx, 
                                                                                                                         ymin = 0), alpha = 0.2, show.legend = FALSE, colour = NA, 
                                                            fill = ribbon_col)
    p1 <- p1 + geom_ribbon(data = with(df1, subset(df1, 
                                                   x <= ql & expert == lpname[2])), aes(ymax = fx, 
                                                                                        ymin = 0), alpha = 0.2, show.legend = FALSE, colour = NA, 
                           fill = ribbon_col) + geom_ribbon(data = with(df1, 
                                                                        subset(df1, x >= qu & expert == lpname[2])), aes(ymax = fx, 
                                                                                                                         ymin = 0), alpha = 0.2, show.legend = FALSE, colour = NA, 
                                                            fill = ribbon_col)
  }
  if (lpname[1] == "marginal") {
    p1 <- p1 + theme(legend.title = element_blank())
  }
  p1 + theme(text = element_text(size = fs))
}

tmpfun <- get("makePoolPlot", envir = asNamespace("expertsurv"))
environment(makePoolPlot) <- environment(tmpfun)
attributes(makePoolPlot) <- attributes(tmpfun)
assignInNamespace("makePoolPlot", makePoolPlot, ns="expertsurv")

###############################################################################


# df <- expertsurv::data %>% rename(status = censored) %>% mutate(time2 = ifelse(time > 10, 10, time),
#                                                                 status2 = ifelse(time> 10, 0, status))



options(spinner.color="#0275D8", spinner.color.background="#ffffff", spinner.size=2)

`%!in%` = Negate(`%in%`)


ui <- fluidPage(useShinyjs(),
               # tags$style('.container-fluid {
               #               background-color: #7b8cde;
               #}'),
                titlePanel("ShinyExpertsurv"),
                sidebarPanel(
  wellPanel(
    #fluidRow(column(3, downloadButton("report", "Download report")),
     #column(2, offset = 1, actionButton("exit", "Quit"))),
    fileInput('df_upload', 'Choose Excel data file to upload',
              accept = c(".xlsx")),
    p("Data should have the following columns: time and status. If your data has two treatment arms please include an arm column with numeric values indicating the treatments."),
    numericInput("n_expert", "Number of Experts", value = 1, min = 1),
    numericInput("n_timepoint", "Number of Timepoints", value = 1,min = 1,  max = 2),
     numericInput("scale1", "Scale for density", value = 1),
    numericInput("xlim", "Limit of x-axis on Kaplan-Meier curve", value = round(10,#max(df$time)*2,
                                                          digits = 0)),
    selectInput(inputId ="pool_type_eval", label = "Pooling approach for experts", 
                choices = c("Linear Pool" = "linear pool",
                            "Logarithmic Pool"= "log pool"), 
                selected = "linear pool"),
    selectInput(inputId ="dist_select", label = "Select the best fitting distribution for Expert Pooling", 
                choices = c("Best Fitting" = "best",
                            "Normal"= "normal",
                            "T-distribution" = "t",
                            "Gamma" = "gamma",
                            "Log-Normal" = "lognormal",
                            "Beta" = "beta"), 
                selected = "best"),
    actionButton(paste0('update_expert'), "Plot/Update Expert Opinions")
    
  ),
  
  hr(),
  
  tabsetPanel(id = "Timepoints",
    tabPanel("Timepoints1",
             numericInput(paste0("time1"), label= "Timepoint", value= 1),
              textInput('quant_vec1', 'Enter a Vector of Quantiles', "0.025,0.5,0.975"),
              matrixInput(
                inputId = "matrix1",
                value = m_default_gen(),
                class = "numeric",
                cols = list(names = TRUE,
                            editableNames = FALSE),
                rows = list(names = FALSE,
                            editableNames = FALSE)),
             helpText("Enter the judgements in the table below,
                            one column per expert. Enter quantile values 
                        corresponding to the cumulative probabilities."),
              
              plotOutput(paste0("expert_plot1"))),
             
      tabPanel("Timepoints2",
               numericInput(paste0("time2"), label= "Timepoint", value= 1),
               textInput('quant_vec2', 'Enter a Vector of Quantiles', "0.025,0.5,0.975"),
               matrixInput(
                 inputId = "matrix2",
                 value = m_default_gen(),
                 class = "numeric",
                 cols = list(names = TRUE,
                             editableNames = FALSE),
                 rows = list(names = FALSE,
                             editableNames = FALSE)),
               helpText("Enter the judgements in the table below,
                            one column per expert. Enter quantile values 
                        corresponding to the cumulative probabilities."),
                plotOutput(paste0("expert_plot2")))
  )),
  mainPanel(
    #withSpinner(tableOutput('tb'), type = 2),
    h3("Kaplan-Meier Survival Plot"),
    plotOutput(paste0("plot_km_expert1")),
    selectInput("opinion_type", label = "Choose opinion type", 
                choices = c("Survival at timepoint(s)" = "survival",
                            "Mean difference between survival"= "mean",
                            "No expert opinion" = "no_expert"), 
                  selected = "survival"),
  
    selectInput("stat_type", label = "Choose statistical approach", 
                choices = c("Frequentist" = "mle","Bayesian" = "hmc"), 
                selected = "mle"),
    numericInput("id_trt", label = "Select treatment ID corresponding to expert opinion", value = 1),
    
    
    pickerInput(
      inputId = "param_mod", 
      label = "Choose models:", 
      choices = c("Exponential" = "exp",
                  "Weibull" = "wei",
                  "Gompertz" = "gomp",
                  "Log-Logistic"= "llo",
                  "Log-normal" = "lno",
                  "Generalized-Gamma" = "gga",
                  "Royston-Parmar" = "rps"), 
      options = list(
        `actions-box` = TRUE, 
        size = 10,
        `selected-text-format` = "count > 3"
      ), 
      multiple = TRUE,
      selected  = c("exp", "wei")
    ),
    selectInput("incl_psa", label = "Include Statistical Uncertainty in Plots", 
                choices = c("Yes" = "yes",
                            "No"= "no"), 
                selected = "no"),  
    
   actionButton("run_analysis", "Run Analysis"),
   #br(),
   #h4("Output of the Statistical Analysis"),
   selectInput("gof_type", label = "Choose goodness of fit measure", 
               choices = c("AIC" = "aic","BIC" = "bic"), 
               selected = "AIC"),
    plotOutput("plot_gof"),
   # selectInput("outFormat", label = "Report format",
   #             choices = list('html' = "html_document",
   #                            'pdf' = "pdf_document",
   #                            'Word' = "word_document"))    #textOutput("txt"),
   
   textInput('file_name', 'Enter a File name to save output', "Output-File"),
   actionButton("save_output", "Save current files")
  )
  
)



server <- function(input, output, session) {

  value <- reactiveValues(
                          m_default = m_default_gen(),
                          n_expert_prev = 1,
                          quant_vec2 = NULL,
                          id_trt = NULL) #Up to max timepoints
  
  
  observeEvent(input$stat_type,{
    if(input$stat_type == "mle"){
      updateSelectInput(session,"gof_type",choices = c("AIC" = "aic", "BIC" = "bic"))
    }else{
      updateSelectInput(session,"gof_type",choices = c("WAIC" = "waic", "PML" = "pml"))
    }

  })
  
  observeEvent(input$df_upload,{
    inFile <- input$df_upload
    df_upload <- readxl::read_excel(inFile$datapath)
    
    if(is.null(df_upload[["arm"]])){
      #browser()
      result.km <- survfit(Surv(time, status) ~ 1, data = df_upload, conf.type="log-log")
      km.data <- data.frame(cbind(result.km[[c("time")]],
                                  result.km[[c("surv")]],
                                  result.km[[c("upper")]],
                                  result.km[[c("lower")]],
                                  arm = 1))
      
      updateSelectInput(session,"opinion_type",choices = c("Survival at timepoint(s)" = "survival",
                                                           "No expert opinion" = "no_expert"))
      hide("id_trt") #hide id_trt panel
      value$id_trt <- NULL
      
      df_upload$arm <- 1
    }else{
      #browser()
      show("id_trt") #hide id_trt panel
      updateNumericInput(inputId = "id_trt", min =min(df_upload$arm),max = max(df_upload$arm), value = max(df_upload$arm))
      km.data <- NULL
      for(i in unique(df_upload$arm)){
        df_temp <- df_upload %>% filter(arm == i)
        result.km_temp <- survfit(Surv(time, status) ~ 1, data = df_temp, conf.type="log-log")
        km.data_temp <- data.frame(cbind(result.km_temp[[c("time")]],
                                         result.km_temp[[c("surv")]],
                                         result.km_temp[[c("upper")]],
                                         result.km_temp[[c("lower")]],
                                         arm = i))
        
        km.data <-  rbind(km.data,km.data_temp)
      }
      
     updateSelectInput(session,"opinion_type",
                       choices = c("Survival at timepoint(s)" = "survival",
                                   "Mean difference between survival"= "mean",
                                   "No expert opinion" = "no_expert"),
                        selected = "survival")
      
    }
   
    colnames(km.data) <- c("Time", "Survival", "upper", "lower", "arm")
    
    value$km.data <- km.data
    value$df_upload <- df_upload
    value$id_trt <- input$id_trt
    #Need to adjust for arm
    
    plot_fit <- ggplot(value$km.data, aes(x = Time,y =Survival, col = factor(arm)))+
      geom_step()+
      ylim(0,1)+
      xlim(0, input$xlim)+
      geom_step(aes(x  = Time, y =upper, col = factor(arm)))+
      geom_step(aes(x  = Time, y =lower, col = factor(arm)))+
      theme_light()#+
      #scale_x_continuous(expand = c(0, 0))+#, breaks=seq(0, 30, 2)) + 
      #scale_y_continuous(expand = c(0, 0))#, breaks=seq(0, 1, 0.05))

    

    output$plot_km_expert1<- renderPlot(plot_fit)
    
  })
  
  observeEvent(input$n_timepoint, {
    
    if(input$n_timepoint > 1){
      showTab(inputId = "Timepoints", target = "Timepoints1")
      showTab(inputId = "Timepoints", target = "Timepoints2")
    }
    if(input$n_timepoint == 1){
      showTab(inputId = "Timepoints", target = "Timepoints1")
      hideTab(inputId = "Timepoints", target = "Timepoints2")
    }
  })
  
  observeEvent(input$opinion_type,{
    if(input$opinion_type == "survival"){
      
      if(input$n_timepoint > 1){
        showTab(inputId = "Timepoints", target = "Timepoints1")
        showTab(inputId = "Timepoints", target = "Timepoints2")
      }
      if(input$n_timepoint == 1){
        showTab(inputId = "Timepoints", target = "Timepoints1")
        hideTab(inputId = "Timepoints", target = "Timepoints2")
      }
      
    }
    
    if(input$opinion_type == "mean"){
      hideTab(inputId = "Timepoints", target = "Timepoints2")
    }
    
    if(input$opinion_type == "no_expert"){
      hideTab(inputId = "Timepoints", target = "Timepoints1")
      hideTab(inputId = "Timepoints", target = "Timepoints2")
    }
    
    if(input$opinion_type == "survival" | input$opinion_type == "mean"){
      show("n_expert")
      show("n_timepoint")
      show("scale1")
      show("pool_type_eval")
      show("dist_select")
      if(is.null(value$id_trt)){
        hide("id_trt")
      }else{
        hide("id_trt")
      }
      

    }else{ #No Expert Opinion
      hide("n_expert")
      hide("n_timepoint")
      hide("scale1")
      hide("pool_type_eval")
      hide("dist_select")
      hide("id_trt")
    }

  })
  
  observeEvent(input$id_trt,{
    value$id_trt <- input$id_trt
  })
  
  observeEvent(input$n_expert, {
   # browser()
    if(input$n_expert == 1){
      shinyjs::hideElement(id = "pool_type_eval")
      hide("pool_type_eval")
      hide("dist_select")
    }else{
      shinyjs::showElement(id = "pool_type_eval")
      show("pool_type_eval")
      show("dist_select")
    }
    
    for(i in 1:2){ #Modify this force it to me 2 which is the max number of timepoints
      mat_exist <- input[[paste0("matrix",i)]]
      if(input$n_expert > value$n_expert_prev){
        extra_cols <- input$n_expert - value$n_expert_prev 
        mat_bind <- matrix(nrow = nrow(mat_exist), ncol = extra_cols)
        mat_exist <- cbind(mat_exist,mat_bind)
      }else if(input$n_expert == value$n_expert_prev){
      } else{
        mat_exist <- mat_exist[,1:(input$n_expert+1),drop = F]
      }
      colnames(mat_exist) <- c("Cum Prob", paste0("Expert_", 1:input$n_expert))
      updateMatrixInput(session, paste0("matrix",i), value=mat_exist )
    }
    value$n_expert_prev <- input$n_expert

  })
  
  
  toListen <- reactive({
    list(input$quant_vec1,input$quant_vec2)
  })

  observeEvent(toListen(),{
   
    for(i in 1:input$n_timepoint){#max number of quant_vec
#browser()
      if(!is.null(input[[paste0("quant_vec",i)]])){
        quant_vec_temp <- input[[paste0("quant_vec",i)]]
        quant_num <- as.numeric(unlist(strsplit(quant_vec_temp,",")))
        if(length(quant_num)==0){
          new_mat <-   matrix(ncol = input$n_expert +1, nrow = 1) # Handle case when nothing is entered
          colnames(new_mat) <- c("Cum Prob", paste0("Expert_",1:input$n_expert ))

        }else{
          mat_exist <- input[[paste0("matrix",i)]]
          retain_quant_index <-which(mat_exist[,1] %in% quant_num)
          retain_quant <- mat_exist[retain_quant_index,1]
          change_quant_index <- which(quant_num %!in% retain_quant)
          new_mat <- matrix(ncol = input$n_expert +1, nrow = length(quant_num))
          colnames(new_mat) <- c("Cum Prob", paste0("Expert_",1:input$n_expert ))

          if(length(retain_quant_index)>0){
            new_mat[1:length(retain_quant_index),] <-mat_exist[retain_quant_index,]

          }
          if(length(change_quant_index)>0){
            new_mat[(length(retain_quant_index)+1):nrow(new_mat),1] <- quant_num[change_quant_index]
          }
        }
        updateMatrixInput(session, paste0("matrix",i), value=new_mat)

      }
    }})

  
  observeEvent(input$update_expert,{
    #pool_type_eval <- "linear pool"
    times_expert_vec <- c()
    df.linear_all <- NULL
    param_expert <- list()
   
    
    plot_fit <- ggplot(value$km.data, aes(x = Time,y =Survival, col = factor(arm)))+
      geom_step()+
      ylim(0,1)+
      xlim(0, input$xlim)+
      geom_step(aes(x  = Time, y =upper, col = factor(arm)))+
      geom_step(aes(x  = Time, y =lower, col = factor(arm)))+
      theme_light()#+
      #scale_x_continuous(expand = c(0, 0))+#, breaks=seq(0, 30, 2)) + 
      #scale_y_continuous(expand = c(0, 0))
    
    if(!any(is.na(input[["matrix1"]][,2]))){ # If Expert opinions are not NA values
      
    
     for(i in 1:input$n_timepoint){ #Update
       
       #want to allow for zeros
      output_pool <- return_pooled_info(input[[paste0("matrix",i)]], St_indic = 0,dist = input$dist_select, mode =NULL)
      if(input$opinion_type == "survival"){
        output_pool[[2]] <- output_pool[[2]]+ xlim(c(0,1)) #If survival we want to truncate.
      }
   
      
      output[[paste0("expert_plot",i)]] <- renderPlot(output_pool[[2]])
      
      times_expert = input[[paste0("time",i)]]
      times_expert_vec <- c(times_expert_vec, times_expert)
      
      df.linear <- subset(output_pool[[2]]$data, ftype == input$pool_type_eval) %>% rename(y = x) %>% 
        mutate(x = times_expert + fx*input[[paste0("scale1")]], 
               times_expert = times_expert)
      df.linear_all <- rbind(df.linear_all, df.linear)
      
      output_pool[[1]][,"dist"] <- gsub("normal", "norm", output_pool[[1]][,"dist"])
      
      param_expert[[i]] <- output_pool[[1]]
     }
      
      value$param_expert <- param_expert
      value$timepoint_expert <- times_expert_vec
      value$df.linear_all <- df.linear_all
      
      if(input$opinion_type == "survival"){
        plot_fit <- plot_fit+
          geom_ribbon(data = df.linear_all, aes(x = x, y = y, xmin= x, xmax =times_expert, group=times_expert), 
                      fill = "sky blue", alpha = 0.5, colour = "grey")
      }

      
    }
    
    output$plot_km_expert1<- renderPlot(plot_fit)
    

    
  })
  
  observeEvent(input$run_analysis, {

    #browser()
    if(length(unique(value$df_upload[["arm"]]))==1){
      formula_text <- "Surv(time,status)~1"
    }else{
      formula_text <- "Surv(time,status)~factor(arm)"
    }
    
    if(!is.null(value$param_expert)& input$opinion_type != "no_expert"){
      
      mod_fit  <- fit.models.expert(formula=as.formula(formula_text),data=value$df_upload,
                                    distr=input$param_mod,
                                    method=input$stat_type,
                                    pool_type = input$pool_type_eval,#"log pool", 
                                    opinion_type = input$opinion_type,
                                    times_expert = value$timepoint_expert, 
                                    param_expert = value$param_expert,
                                    id_trt = input$id_trt,
                                    id_St  = input$id_trt,
                                    k = 1)
  
        value$mod_fit <- mod_fit
    }
    
    if(input$opinion_type == "no_expert"){
      
      
      param_expert_vague <- list()
      param_expert_vague[[1]] <- data.frame(dist = "beta", wi = 1, param1 = 1, param2 = 1, param2 = NA)
      
      mod_fit  <- fit.models.expert(formula=as.formula(formula_text),data=value$df_upload,
                                    distr=input$param_mod,
                                    method=input$stat_type,
                                    pool_type = input$pool_type_eval,#"log pool", 
                                    opinion_type = "survival",
                                    times_expert = 2, 
                                    param_expert = param_expert_vague,
                                    k = 1,
                                    id_St  = 1)
      
      value$mod_fit <- mod_fit
    }

  })
  
  value$plot_km_expert1 <- eventReactive(value$mod_fit,{
    
    if(input$incl_psa == "yes"){
      #browser()
      models <- names(value$mod_fit$models)
      psa_outuput <- list()
      
      for(i in 1:length(value$mod_fit$models)){
        psa <- make.surv(fit = value$mod_fit,mod = i, nsim = 1000, t = seq(0,input$xlim,length.out = 1000))
        df_temp  <- t(apply(psa$mat[[1]], 1,quantile, probs = c(0.025, 0.5,.975))) %>% data.frame()
        df_temp$time <- seq(0,input$xlim,length.out = 1000)
        mod_name <- names(value$mod_fit$models)[i]
        psa_outuput[[mod_name]] <- df_temp %>% mutate(model = mod_name)
      }
      
      df_final_plot <- do.call(rbind.data.frame, psa_outuput)
      df_final_plot$models <- factor(df_final_plot$model, levels = unique(df_final_plot$model))
      #browser()
    if(input$opinion_type == "survival"){
      ggplot(data = df_final_plot, aes(y = X50., x = time, group = models, colour = models))+
        geom_line()+
        geom_line(data = df_final_plot ,aes(y = X97.5., x = time),linetype="dotdash")+
        geom_line(data = df_final_plot ,aes(y = X2.5., x = time), linetype="dotdash")+
        geom_step(data =value$km.data, mapping = aes(x = Time,y =Survival, col = factor(arm)),inherit.aes = FALSE)+
        geom_ribbon(data = value$df.linear_all, aes(x = x, y = y, xmin= x, xmax =times_expert, group=times_expert), 
                    fill = "sky blue", alpha = 0.5, colour = "grey")
    }else{
      ggplot(data = df_final_plot, aes(y = X50., x = time, group = models, colour = models))+
        geom_line()+
        geom_line(data = df_final_plot ,aes(y = X97.5., x = time),linetype="dotdash")+
        geom_line(data = df_final_plot ,aes(y = X2.5., x = time), linetype="dotdash")+
        geom_step(value$km.data, aes(x = Time,y =Survival, col = factor(arm)),inherit.aes = FALSE)+
        geom_step(value$km.data, aes(x = Time,y =lower, col = factor(arm)),inherit.aes = FALSE)+
        geom_step(value$km.data, aes(x = Time,y =upper, col = factor(arm)),inherit.aes = FALSE)
      }
      
    }else{
     
    if(input$opinion_type == "survival"){
    plot(value$mod_fit, add.km = TRUE,t = seq(0,input$xlim,length.out = 1000))+
      geom_ribbon(data = value$df.linear_all, aes(x = x, y = y, xmin= x, xmax =times_expert, group=times_expert), 
                  fill = "sky blue", alpha = 0.5, colour = "grey")
      
    }else{
      plot(value$mod_fit, add.km = TRUE,t = seq(0,input$xlim,length.out = 1000))
      
    } 
      
    }  
    })
  
  value$plot_gof <- eventReactive(value$mod_fit,{
    #browser()
    model.fit.plot(value$mod_fit,type = input$gof_type)
  })

  observeEvent(input$update_expert, {
   # browser()
    hide("plot_gof")
  })
  
  
  observeEvent(input$run_analysis, {
    
  output$plot_km_expert1 <- renderPlot(
    value$plot_km_expert1())
    
  output$plot_gof <- renderPlot(
       value$plot_gof())
   
  show("plot_gof")
    
  })
  
  observeEvent(input$save_output,{
    #browser()
    saveRDS(list(model = value$mod_fit, surv_plt = value$plot_km_expert1(), gof_plt = value$plot_gof()),
         file = paste0(input$file_name,".rds"))
    
    #readRDS(file = paste0(input$file_name,".rds"))
  })
  
  
}

shinyApp(ui, server)
