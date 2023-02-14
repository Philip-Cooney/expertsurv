#' logLikFactory
#'
#' @param Y 
#' @param X 
#' @param weights 
#' @param bhazard 
#' @param rtrunc 
#' @param dlist 
#' @param inits 
#' @param dfns 
#' @param aux 
#' @param mx 
#' @param fixedpars 
#' @param expert_opinion 
#'
#' @import flexsurv
#' @importFrom abind adrop
#' @noRd
logLikFactory <- function (Y, X = 0, weights, bhazard, rtrunc, dlist, inits, dfns,
                           aux, mx, fixedpars = NULL, expert_opinion){
  
  
  
  #NAMESPACE HACK FOR CRAN; won't let me use flexsurv::: :-(
  
  buildTransformer <-utils::getFromNamespace("buildTransformer", "flexsurv")
  buildAuxParms <-utils::getFromNamespace("buildAuxParms", "flexsurv")
  check.flexsurv.response <-utils::getFromNamespace("check.flexsurv.response", "flexsurv")
  tsfn <-utils::getFromNamespace("tsfn", "flexsurv")
  DLdsurvspline <-utils::getFromNamespace("DLdsurvspline", "flexsurv")
  flexsurv.splineinits.cox <-utils::getFromNamespace("flexsurv.splineinits.cox", "flexsurv")
  parse.dist <-utils::getFromNamespace("parse.dist", "flexsurv")
  form.dp<-utils::getFromNamespace("form.dp", "flexsurv")
  check.formula<-utils::getFromNamespace("check.formula", "flexsurv")
  anc_from_formula<-utils::getFromNamespace("anc_from_formula", "flexsurv")
  get.locform<-utils::getFromNamespace("get.locform", "flexsurv")
  concat.formulae<-utils::getFromNamespace("concat.formulae", "flexsurv")
  compress.model.matrices<-utils::getFromNamespace("compress.model.matrices", "flexsurv")
  expand.inits.args<-utils::getFromNamespace("expand.inits.args", "flexsurv")
  .hess_to_cov<-utils::getFromNamespace(".hess_to_cov", "flexsurv")
  logh<-utils::getFromNamespace("logh", "flexsurv")
  dexph<-utils::getFromNamespace("dexph", "flexsurv")
  deriv.test<-utils::getFromNamespace("deriv.test", "flexsurv")
  
  
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
      #expert_opinion
      
      if(expert_opinion$St_indic ==1){
        pargs.surv$q <- expert_opinion$times
        psurv_expert <- 1 - do.call(dfns$p, pargs.surv)
      }else{
        pargs.rmst <- pargs.surv
        pargs.rmst$q <- NULL
        pargs.rmst <-  lapply(pargs.rmst, function(x){x[c(expert_opinion$id_comp,expert_opinion$id_trt)]})
        psurv_expert<- diff(do.call(dfns$mean, pargs.rmst))
      }
      #print(psurv_expert)

      LL_expert <- rep(NA, length(expert_opinion$times))
      for(q in 1:length(expert_opinion$times)){
       param_expert_mat  <- abind::adrop(expert_opinion$param_expert[,,q, drop = F],drop = 3)
       #LL_expert[q] <-  log_density_dist(param_expert_mat,psurv_expert[q],nrow(param_expert_mat), expert_opinion$pool)
       LL_expert[q] <- expert_log_dens(psurv_expert[q], df= param_expert_mat, expert_opinion$pool, k_norm = expert_opinion$k_norm[q],
                                                    St_indic =expert_opinion$St_indic)

      }
      ret <- -sum(loglik, LL_expert)
    }else{
      ret <- -sum(loglik)
    }
    
    attr(ret, "indiv") <- loglik
    #There will always be a darg but not always a parg, however you would never need an opinion if everyone died!
    ret
  }
}

flexsurvspline <- function (formula, data, weights, bhazard, rtrunc, subset, k = 0, 
          knots = NULL, bknots = NULL, scale = "hazard", timescale = "log", expert_opinion =NULL,
          ...){
  
  #NAMESPACE HACK FOR CRAN; won't let me use flexsurv::: :-(
  
  buildTransformer <-utils::getFromNamespace("buildTransformer", "flexsurv")
  buildAuxParms <-utils::getFromNamespace("buildAuxParms", "flexsurv")
  check.flexsurv.response <-utils::getFromNamespace("check.flexsurv.response", "flexsurv")
  tsfn <-utils::getFromNamespace("tsfn", "flexsurv")
  DLdsurvspline <-utils::getFromNamespace("DLdsurvspline", "flexsurv")
  flexsurv.splineinits.cox <-utils::getFromNamespace("flexsurv.splineinits.cox", "flexsurv")
  flexsurv.splineinits <-utils::getFromNamespace("flexsurv.splineinits", "flexsurv")
  parse.dist <-utils::getFromNamespace("parse.dist", "flexsurv")
  form.dp<-utils::getFromNamespace("form.dp", "flexsurv")
  check.formula<-utils::getFromNamespace("check.formula", "flexsurv")
  anc_from_formula<-utils::getFromNamespace("anc_from_formula", "flexsurv")
  get.locform<-utils::getFromNamespace("get.locform", "flexsurv")
  concat.formulae<-utils::getFromNamespace("concat.formulae", "flexsurv")
  compress.model.matrices<-utils::getFromNamespace("compress.model.matrices", "flexsurv")
  expand.inits.args<-utils::getFromNamespace("expand.inits.args", "flexsurv")
  .hess_to_cov<-utils::getFromNamespace(".hess_to_cov", "flexsurv")
  logh<-utils::getFromNamespace("logh", "flexsurv")
  dexph<-utils::getFromNamespace("dexph", "flexsurv")
  deriv.test<-utils::getFromNamespace("DLSsurvspline ", "flexsurv")
  DLSsurvspline <-utils::getFromNamespace("deriv.test", "flexsurv")
  
  
  
  call <- match.call()
  indx <- match(c("formula", "data", "weights", "bhazard", 
                  "rtrunc", "subset", "na.action"), names(call), nomatch = 0)
  if (indx[1] == 0) 
    stop("A \"formula\" argument is required")
  temp <- call[c(1, indx)]
  temp[[1]] <- as.name("model.frame")
  f2 <- stats::as.formula(gsub("(.*~).*", "\\1 1", Reduce(paste, deparse(formula))))
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
    knots <- stats::quantile(tsfn(dtimes, timescale), seq(0, 1, 
                                                   length = k + 2)[-c(1, k + 2)])
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
  dfn <- flexsurv::unroll.function(flexsurv::dsurvspline, gamma = 0:(nk - 1))
  pfn <- flexsurv::unroll.function(flexsurv::psurvspline, gamma = 0:(nk - 1))
  rfn <- flexsurv::unroll.function(flexsurv::rsurvspline, gamma = 0:(nk - 1))
  hfn <- flexsurv::unroll.function(flexsurv::hsurvspline, gamma = 0:(nk - 1))
  Hfn <- flexsurv::unroll.function(flexsurv::Hsurvspline, gamma = 0:(nk - 1))
  qfn <- flexsurv::unroll.function(flexsurv::qsurvspline, gamma = 0:(nk - 1))
  meanfn <- flexsurv::unroll.function(flexsurv::mean_survspline, gamma = 0:(nk - 
                                                          1))
  rmstfn <- flexsurv::unroll.function(flexsurv::rmst_survspline, gamma = 0:(nk - 
                                                          1))
  Ddfn <- if (scale == "normal") 
    NULL
  else flexsurv::unroll.function(DLdsurvspline, gamma = 0:(nk - 1))
  DSfn <- if (scale == "normal") 
    NULL
  else flexsurv::unroll.function(DLSsurvspline, gamma = 0:(nk - 1))
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
  
  if(is.null(expert_opinion)){
    expert_opinion <- vector(mode = "list", length = 1)
    names(expert_opinion) <- "expert_opinion"
    args <- append(args,expert_opinion)
  }else{
    args[["expert_opinion"]] <- expert_opinion
  }

  #browser()
  if (is.infinite(do.call("flexsurvreg", args)$loglik)) {
    args$dist$inits <- flexsurv.splineinits.cox
  }
  args$fixedpars <- fpold
  #browser()
  ret <- do.call("flexsurvreg", args)
  ret <- c(ret, list(k = length(knots) - 2, knots = knots, 
                     scale = scale, timescale = timescale))
  ret$call <- call
  class(ret) <- "flexsurvreg"
  ret
}

flexsurvreg <- function (formula, anc = NULL, data, weights, bhazard, rtrunc,
                         subset, na.action, dist, inits, fixedpars = NULL, dfns = NULL,
                         aux = NULL, cl = 0.95, integ.opts = NULL, sr.control = survival::survreg.control(),
                         hessian = TRUE, hess.control = NULL, expert_opinion = NULL, ...){
  
  #NAMESPACE HACK FOR CRAN; won't let me use flexsurv::: :-(
  
  buildTransformer <-utils::getFromNamespace("buildTransformer", "flexsurv")
  buildAuxParms <-utils::getFromNamespace("buildAuxParms", "flexsurv")
  check.flexsurv.response <-utils::getFromNamespace("check.flexsurv.response", "flexsurv")
  tsfn <-utils::getFromNamespace("tsfn", "flexsurv")
  DLdsurvspline <-utils::getFromNamespace("DLdsurvspline", "flexsurv")
  flexsurv.splineinits.cox <-utils::getFromNamespace("flexsurv.splineinits.cox", "flexsurv")
  parse.dist <-utils::getFromNamespace("parse.dist", "flexsurv")
  form.dp<-utils::getFromNamespace("form.dp", "flexsurv")
  check.formula<-utils::getFromNamespace("check.formula", "flexsurv")
  anc_from_formula<-utils::getFromNamespace("anc_from_formula", "flexsurv")
  get.locform<-utils::getFromNamespace("get.locform", "flexsurv")
  concat.formulae<-utils::getFromNamespace("concat.formulae", "flexsurv")
  compress.model.matrices<-utils::getFromNamespace("compress.model.matrices", "flexsurv")
  expand.inits.args<-utils::getFromNamespace("expand.inits.args", "flexsurv")
  .hess_to_cov<-utils::getFromNamespace(".hess_to_cov", "flexsurv")
  logh<-utils::getFromNamespace("logh", "flexsurv")
  dexph<-utils::getFromNamespace("dexph", "flexsurv")
  deriv.test<-utils::getFromNamespace("deriv.test", "flexsurv")
  check.fixedpars<-utils::getFromNamespace("check.fixedpars", "flexsurv")
  
  
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
  Y <- check.flexsurv.response(stats::model.extract(m, "response"))
  mml <- mx <- vector(mode = "list", length = length(dlist$pars))
  names(mml) <- names(mx) <- c(dlist$location, setdiff(dlist$pars,
                                                       dlist$location))
  for (i in names(forms)) {
    mml[[i]] <- stats::model.matrix(forms[[i]], m)
    mx[[i]] <- length(unlist(mx)) + seq_len(ncol(mml[[i]][,
                                                          -1, drop = FALSE]))
  }
  X <- compress.model.matrices(mml)
  weights <- stats::model.extract(m, "weights")
  if (is.null(weights))
    weights <- m$"(weights)" <- rep(1, nrow(X))
  bhazard <- stats::model.extract(m, "bhazard")
  if (is.null(bhazard))
    bhazard <- rep(0, nrow(X))
  rtrunc <- stats::model.extract(m, "rtrunc")
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
                             weights = temp$weights, control = sr.control, counting = (attr(stats::model.extract(m,
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
                                        mx = mx, expert_opinion = expert_opinion)
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
      lcl <- est - stats::qnorm(1 - (1 - cl)/2) * se
      ucl <- est + stats::qnorm(1 - (1 - cl)/2) * se
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

gen_pathway<- function(string){
  folder_bits <- strsplit(string, "\\/")[[1]]
  
  file_extensions <- c("R", "png", "xlsx", "csv", "jpg", "svg", "txt", "pptx","pdf")
  file_extensions_grep <- paste0(paste0("\\.",file_extensions), collapse ="|")
  
  
  exclude.last <- grepl(file_extensions_grep,folder_bits[length(folder_bits)], ignore.case = T)
  
  
  if(exclude.last){
    folder_num <- length(folder_bits) -1
    
  }else{
    folder_num <- length(folder_bits) 
  }
  
  for(i in 1:folder_num){
    folder_path_temp <- paste0(folder_bits[1:i], collapse = "/")
    
    if(!file.exists(folder_path_temp)){
      dir.create(folder_path_temp)
    }
  }
  
  return(string) #Need to return string for actual pasting
  
}
pow <- function(x,y){
  x^y
}
get_k_norm <- function(opinion_list, St_indic = 1){ # Only required if log-pooling
  k_norm <- rep(NA, length(opinion_list))
  
  if(St_indic ==1){
    a <- 0
    b <- 1
  }else{
    a <- -Inf
    b <- +Inf
  }
  
  for(i in 1:length(opinion_list)){
    
    opinion_list[[i]]$dist <- stringr::str_replace_all( opinion_list[[i]]$dist, "normal", "norm") 
    opinion_list[[i]]$dist <- stringr::str_replace_all( opinion_list[[i]]$dist, "lognorm", "lnorm") 
    
    pool.df <- opinion_list[[i]]
    
    if(St_indic ==1){
      min_quant <- 0
      max_quant <- 1
    }else{
      quant.vec <- t(apply(pool.df, 1, function(x){get_quant_val(
        dist = x["dist"],
        param1 = x["param1"],
        param2 = x["param2"],
        param3 = x["param3"],
        probs = c(0.001,0.025,0.5,0.975,0.999))}))
      # central.cauchy <- mean(quant.vec[,3])#mean
      # sd.cauchy <- max(apply(quant.vec,1, function(x){(x[4]-x[2])/4})) #sd
      min_quant <- min(quant.vec)
      max_quant <- max(quant.vec)
    }
    x.eval <- seq(min_quant, max_quant, length.out = 100)
    #Find the range of x.eval to integrate over
    # x.eval <- seq(0,1, by = 0.001)
    dens_eval<- eval_dens_pool(x.eval,opinion_list[[i]],  pool_type = "log pool", St_indic =St_indic)
    k_norm[i] <- integrate.xy(x = x.eval,fx = dens_eval)
   }
  k_norm
}

#' Helper function to run the survival models using MLE and flexsurv
#' 
#' @param x a string containing the name of the model
#' to be fitted
#' #' @param exArgs a list of extra arguments passed from the main 'fit.models' 
#' function
#' @note Something will go here
#' @author Gianluca Baio
#' @seealso fit.models
#' @references Baio (2020). survHE
#' @keywords Parametric survival models Maximum likelihood estimation
#' @import tibble
#' @import dplyr
#' @import stats 
#' @import flexsurv
#' @noRd 
runMLE <- function (x, exArgs){
  #browser()
  formula <- exArgs$formula
  data = exArgs$data
  availables <- load_availables()
  d3 <- manipulate_distributions(x)$distr3
  x <- manipulate_distributions(x)$distr
  
  expert_opinion_flex<-  list()
  
  
  ## -- Create Model Frame
  
  if(exArgs$opinion_type != "survival"){# Survival Difference
   #browser() 
  
  formula_temp <- update(formula, paste(all.vars(formula, data)[1], 
                                        "~", all.vars(formula, data)[2], "+."))
  mf <- tibble::as_tibble(model.frame(formula_temp, data)) %>% 
    dplyr::rename(time = 1,event = 2) %>% rename_if(is.factor, .funs = ~gsub("as.factor[( )]","", .x)) %>% 
    dplyr::rename_if(is.factor, .funs = ~gsub("[( )]","", .x)) %>% 
    dplyr::bind_cols(tibble::as_tibble(stats::model.matrix(formula_temp,data)) %>% dplyr::select(contains("Intercept"))) %>%
    dplyr::select(time,event, contains("Intercept"), everything()) %>% tibble::rownames_to_column("ID")
  
 if(ncol(mf) == 5){
      expert_opinion_flex$id_trt <- min(which(mf[,5] == exArgs$id_trt)) 
      if(length(unique(mf[,5] %>% pull()))==2){
        expert_opinion_flex$id_comp <- min(which(mf[,5] != exArgs$id_trt)) 
      }else{
        expert_opinion_flex$id_comp <- min(which(mf[,5] == exArgs$id_comp))  
      }
 
    #put the number in  could put in a combination of numbers
  }else{
    print("We do not allow more than one covariate (i.e. treatment) in the analysis - although it is technically possible")
    stop()
  }
  
  times <- 99999
  expert_opinion_flex$St_indic <- 0
}
  
  
  if(exArgs$opinion_type == "survival"){
    expert_opinion_flex$St_indic <- 1
    times <- exArgs$times_expert 
  }#else{
  # stop("Expert Opinion Currently only on Survival probabilities;")
  # }
  # 
  expert_opinion_flex$param_expert <- make_data_expert(exArgs$param_expert, times)
  expert_opinion_flex$times <- times
  if(exArgs$pool_type == "linear pool"){
    expert_opinion_flex$pool <- 1 #linear pool is 1; log is 0
  }else{
    expert_opinion_flex$pool <- 0
  }
  
  #Only currently implemented for survival probabilities 
  
  if(expert_opinion_flex$pool == 0){
    expert_opinion_flex$k_norm <-  get_k_norm(exArgs$param_expert)
  }else{
    expert_opinion_flex$k_norm <-  NULL
  }
  
  if(exists("method_mle", where = exArgs)) {
    method_mle <- exArgs$method_mle
  }else{
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
    model_mle <- flexsurvspline(formula = formula, 
                                         data = data, k = k, knots = knots, bknots = bknots, 
                                         scale = scale, timescale = timescale,expert_opinion = NULL , method = method_mle)
    
    model <-flexsurvspline(formula = formula, 
                                      data = data, k = k, knots = knots, bknots = bknots, 
                                      scale = scale, timescale = timescale,expert_opinion = expert_opinion_flex , method = method_mle,inits = model_mle$res[,1])
    })
    }
  else {
    suppressWarnings({
    model_mle <- flexsurvreg(formula = formula, data = data, 
                                      dist = x, expert_opinion = NULL,method = method_mle)
    
    model <- flexsurvreg(formula = formula, data = data, 
                                   dist = x, expert_opinion = expert_opinion_flex,method = method_mle,inits = model_mle$res[,1])
    })
    }
  toc <- proc.time() - tic
  model_name <- d3
  list(model = model, aic = model$AIC, bic = -2 * model$loglik + 
         model$npars * log(model$N), dic = NULL, time2run = toc[3], 
       model_name = model_name)
}
