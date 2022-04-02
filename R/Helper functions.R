# FUNCTION ----
`%!in%` = Negate(`%in%`)
gamma.error_mod <- function (parameters, values, probabilities, weights, mode){
  sum(weights * (pgamma(values, exp(parameters[1]), exp(parameters[2])) - 
                   probabilities)^2) + 
    ((exp(parameters[1])-1)/exp(parameters[2]) -mode)^2 #Mode
}

beta.error_mod <- function (parameters, values, probabilities, weights, mode) 
{
  sum(weights * (pbeta(values, exp(parameters[1]), exp(parameters[2])) - 
                   probabilities)^2) +
    ((exp(parameters[1])-1)/(exp(parameters[1])+exp(parameters[2])-2) -mode)^2 
}


normal.error_mod <- function (parameters, values, probabilities, weights, mode){
  sum(weights * (pnorm(values, parameters[1], exp(parameters[2])) - 
                   probabilities)^2)+
    (parameters[1] -mode)^2
}


t.error_mod <- function (parameters, values, probabilities, weights, degreesfreedom, mode){
  sum(weights * (pt((values - parameters[1])/exp(parameters[2]), 
                    degreesfreedom) - probabilities)^2)+
    (parameters[1] -mode)^2
}


lognormal.error_mod <- function (parameters, values, probabilities, weights, mode){
  sum(weights * (plnorm(values, parameters[1], exp(parameters[2])) - 
                   probabilities)^2)+
    (exp(parameters[1]-exp(parameters[2])^2) -mode)^2
}

dt.scaled <- function (x, df, mean = 0, sd = 1, ncp, log = FALSE){
  if (!log) 
    stats::dt((x - mean)/sd, df, ncp = ncp, log = FALSE)/sd
  else stats::dt((x - mean)/sd, df, ncp = ncp, log = TRUE) - 
    log(sd)
}


expert_log_dens <- function(x, df, pool_type, k_norm = NULL){
  if(is.null(dim(df))){ #corced to vector
    df <-matrix(df, nrow = 1, ncol = length(df))
  }
    
  like <- rep(NA,nrow(df)) 
  for(i in 1:nrow(df)){
    
    if(df[i,1] == 1){ # 1 equal normal
      if(pool_type == 1){
        like[i] <- dnorm(x, df[i,3], df[i,4], log = F)*df[i,2]  
      }else{
        like[i] <- dnorm(x, df[i,3], df[i,4], log = F)^df[i,2] 
      }
      
    }
    
    
    if(df[i,1] == 2){ # 2 equal t
      if(pool_type == 1){
        like[i] <- dt.scaled(x, df[i,5], df[i,3], df[i,4], log = F)*df[i,2]  
      }else{
        like[i] <- dt.scaled(x, df[i,5], df[i,3], df[i,4], log = F)^df[i,2] #
      }
      
    }
    
    if(df[i,1] == 3){ # 3 equal gamma
      if(pool_type == 1){
        like[i] <- dgamma(x, df[i,3], df[i,4],  log = F)*df[i,2]   
      }else{
        like[i] <- dgamma(x,  df[i,3], df[i,4],   log = F)^df[i,2]   
      }
      
    }
    
    if(df[i,1] == 4){ # 4 equal lnorm
      if(pool_type == 1){
        like[i] <- dlnorm(x,  df[i,3], df[i,4],   log = F)*df[i,2]   
      }else{
        like[i] <- dlnorm(x,  df[i,3], df[i,4],  log = F)^df[i,2]   
      }
      
    }
    
    
    if(df[i,1] == 5){# 5 = beta
      if(pool_type == 1){
        like[i] <- dbeta(x,  df[i,3], df[i,4],   log = F)*df[i,2]  
      }else{
        like[i] <- dbeta(x,  df[i,3], df[i,4],   log = F)^df[i,2]  
      }
    }
    
    
    
  }  
  if(pool_type == 1){
    return(log(sum(like)/sum(df[,2])))
  }else{
    return(log(prod(like)/k_norm))
  }
  
  
  
  
}





fitdist_mod <- function (vals, probs, lower = -Inf, upper = Inf, weights = 1, 
                     tdf = 3, expertnames = NULL, excludelog.mirror = TRUE, mode = NULL){
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
  colnames(ssq) <- c("normal", "t", "gamma", 
                     "lognormal", "logt", "beta", "mirrorgamma", 
                     "mirrorlognormal", "mirrorlogt")
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
    if (min(probs[, i]) > 0.4) {
      stop("smallest elicited probability must be less than 0.4")
    }
    if (min(probs[, i]) < 0 | max(probs[, i]) > 1) {
      stop("probabilities must be between 0 and 1")
    }
    if (max(probs[, i]) < 0.6) {
      stop("largest elicited probability must be greater than 0.6")
    }
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
    q.fit <- approx(x = probs[inc, i], y = vals[inc, i], 
                    xout = c(0.4, 0.5, 0.6))$y
    l <- q.fit[1]
    u <- q.fit[3]
    minq <- qnorm(minprob)
    maxq <- qnorm(maxprob)
    m <- (minvals * maxq - maxvals * minq)/(maxq - minq)
    v <- ((maxvals - minvals)/(maxq - minq))^2
    normal.fit <- optim(c(m, 0.5 * log(v)), normal.error_mod, 
                        values = vals[inc, i], probabilities = probs[inc, 
                                                                     i], weights = weights[inc, i], mode = mode[i])
    normal.parameters[i, ] <- c(normal.fit$par[1], exp(normal.fit$par[2]))
    ssq[i, "normal"] <- normal.fit$value
    t.fit <- optim(c(m, 0.5 * log(v)), t.error_mod, values = vals[inc,i], 
                   probabilities = probs[inc, i], weights = weights[inc,i], degreesfreedom = tdf[i], mode = mode[i])
    t.parameters[i, 1:2] <- c(t.fit$par[1], exp(t.fit$par[2]))
    t.parameters[i, 3] <- tdf[i]
    ssq[i, "t"] <- t.fit$value
    if (lower[i] > -Inf) {
      vals.scaled1 <- vals[inc, i] - lower[i]
      
      m.scaled1 <- m - lower[i]
      gamma.fit <- optim(c(log(m.scaled1^2/v), log(m.scaled1/v)), 
                         gamma.error_mod, values = vals.scaled1, probabilities = probs[inc, 
                                                                                   i], weights = weights[inc, i], mode = mode[i])
      gamma.parameters[i, ] <- exp(gamma.fit$par)
      ssq[i, "gamma"] <- gamma.fit$value
      std <- ((log(u - lower[i]) - log(l - lower[i]))/1.35)
      mlog <- (log(minvals - lower[i]) * maxq - log(maxvals - 
                                                      lower[i]) * minq)/(maxq - minq)
      lognormal.fit <- optim(c(mlog, log(std)), lognormal.error_mod, 
                             values = vals.scaled1, probabilities = probs[inc, 
                                                                          i], weights = weights[inc, i], mode = mode[i])
      lognormal.parameters[i, 1:2] <- c(lognormal.fit$par[1], 
                                        exp(lognormal.fit$par[2]))
      ssq[i, "lognormal"] <- lognormal.fit$value
      logt.fit <- optim(c(log(m.scaled1), log(std)), SHELF:::logt.error, 
                        values = vals.scaled1, probabilities = probs[inc, 
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
      beta.fit <- optim(c(log(alp), log(bet)), beta.error_mod, 
                        values = vals.scaled2, probabilities = probs[inc, 
                                                                     i], weights = weights[inc, i], mode = mode[i])
      beta.parameters[i, ] <- exp(beta.fit$par)
      ssq[i, "beta"] <- beta.fit$value
    }
    if (upper[i] < Inf) {
      valsMirrored <- upper[i] - vals[inc, i]
      probsMirrored <- 1 - probs[inc, i]
      mMirrored <- upper[i] - m
      mirrorgamma.fit <- optim(c(log(mMirrored^2/v), log(mMirrored/v)), 
                               SHELF:::gamma.error, values = valsMirrored, probabilities = probsMirrored, 
                               weights = weights[inc, i])
      mirrorgamma.parameters[i, ] <- exp(mirrorgamma.fit$par)
      ssq[i, "mirrorgamma"] <- mirrorgamma.fit$value
      mlogMirror <- (log(upper[i] - maxvals) * (1 - minq) - 
                       log(upper[i] - minvals) * (1 - maxq))/(maxq - 
                                                                minq)
      stdMirror <- ((log(upper[i] - l) - log(upper[i] - 
                                               u))/1.35)
      mirrorlognormal.fit <- optim(c(mlogMirror, log(stdMirror)), 
                                   SHELF:::lognormal.error, values = valsMirrored, probabilities = probsMirrored, 
                                   weights = weights[inc, i])
      mirrorlognormal.parameters[i, 1:2] <- c(mirrorlognormal.fit$par[1], 
                                              exp(mirrorlognormal.fit$par[2]))
      ssq[i, "mirrorlognormal"] <- mirrorlognormal.fit$value
      mirrorlogt.fit <- optim(c(log(mMirrored), log(stdMirror)), 
                              SHELF:::logt.error, values = valsMirrored, probabilities = probsMirrored, 
                              weights = weights[inc, i], degreesfreedom = tdf[i])
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
  names(dflt) <- c("location.log.X", "scale.log.X", 
                   "df.log.X")
  row.names(dflt) <- expertnames
  dfmirrorlt <- data.frame(mirrorlogt.parameters)
  names(dfmirrorlt) <- c("location.log.X", "scale.log.X", 
                         "df.log.X")
  row.names(dfmirrorlt) <- expertnames
  dfb <- data.frame(beta.parameters)
  names(dfb) <- c("shape1", "shape2")
  row.names(dfb) <- expertnames
  ssq <- data.frame(ssq)
  row.names(ssq) <- expertnames
  if (excludelog.mirror) {
    reducedssq <- ssq[, c("normal", "t", "gamma", 
                          "lognormal", "beta")]
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



plotfit <- function (fit, d = "best", xl = -Inf, xu = Inf, ql = NA, qu = NA, 
          lp = FALSE, ex = NA, sf = 3, ind = TRUE, lpw = 1, fs = 12, 
          lwd = 1, xlab = "x", ylab = expression(f[X](x)), legend_full = TRUE, 
          percentages = FALSE, returnPlot = FALSE){
  if (d == "beta" & (min(fit$limits) == -Inf | max(fit$limits) == 
                     Inf)) {
    stop("Parameter limits must be finite to fit a beta distribution")
  }
  if (d == "gamma" & min(fit$limits) == -Inf) {
    stop("Lower parameter limit must be finite to fit a (shifted) gamma distribution")
  }
  if (d == "lognormal" & min(fit$limits) == -Inf) {
    stop("Lower parameter limit must be finite to fit a (shifted) log normal distribution")
  }
  if (d == "logt" & min(fit$limits) == -Inf) {
    stop("Lower parameter limit must be finite to fit a (shifted) log t distribution")
  }
  if (is.na(ql) == F & (ql < 0 | ql > 1)) {
    stop("Lower feedback quantile must be between 0 and 1")
  }
  if (is.na(qu) == F & (qu < 0 | qu > 1)) {
    stop("Upper feedback quantile must be between 0 and 1")
  }
  theme_set(theme_grey(base_size = fs))
  theme_update(plot.title = element_text(hjust = 0.5))
  if (nrow(fit$vals) > 1 & is.na(ex) == T & lp == F) {
    if (xl == -Inf & min(fit$limits[, 1]) > -Inf) {
      xl <- min(fit$limits[, 1])
    }
    if (xu == Inf & max(fit$limits[, 2]) < Inf) {
      xu <- max(fit$limits[, 2])
    }
    p1 <- suppressWarnings(SHELF:::makeGroupPlot(fit, xl, xu, d, 
                                         lwd, xlab, ylab, expertnames = rownames(fit$Normal)))
    #print(p1)
    if (returnPlot) {
      return(p1)
    }
  }
  if (nrow(fit$vals) > 1 & lp == T) {
    if (xl == -Inf & min(fit$limits[, 1]) > -Inf) {
      xl <- min(fit$limits[, 1])
    }
    if (xl == -Inf & min(fit$limits[, 1]) == -Inf) {
      f1 <- feedback(fit, quantiles = 0.01, dist = d)
      xl <- min(f1$expert.quantiles)
    }
    if (xu == Inf & max(fit$limits[, 2]) < Inf) {
      xu <- max(fit$limits[, 2])
    }
    if (xu == Inf & max(fit$limits[, 2]) == Inf) {
      f2 <- feedback(fit, quantiles = 0.99, dist = d)
      xu <- max(f2$expert.quantiles)
    }
    p1 <- SHELF:::makeLinearPoolPlot(fit, xl, xu, d, lpw, lwd, xlab, 
                             ylab, legend_full, expertnames = rownames(fit$Normal))
    #print(p1)
    if (returnPlot) {
      return(p1)
    }
  }
  if (nrow(fit$vals) > 1 & is.na(ex) == F) {
    if (xl == -Inf & fit$limits[ex, 1] > -Inf) {
      xl <- fit$limits[ex, 1]
    }
    if (xu == Inf & fit$limits[ex, 2] < Inf) {
      xu <- fit$limits[ex, 2]
    }
    p1 <- suppressWarnings(SHELF:::makeSingleExpertPlot(fit, d, 
                                                xl, xu, ql, qu, sf, ex = ex, lwd, xlab, ylab, percentages))
    #print(p1)
    if (returnPlot) {
      return(p1)
    }
  }
  if (nrow(fit$vals) == 1) {
    p1 <- suppressWarnings(SHELF:::makeSingleExpertPlot(fit, d, 
                                                xl, xu, ql, qu, sf, ex = 1, lwd, xlab, ylab, percentages))
    #print(p1)
    if (returnPlot) {
      return(p1)
    }
  }
}


`%!in%` <- Negate(`%in%`)

dt.scaled <- function (x, df, mean = 0, sd = 1, ncp, log = FALSE){
  if (!log) 
    stats::dt((x - mean)/sd, df, ncp = ncp, log = FALSE)/sd
  else stats::dt((x - mean)/sd, df, ncp = ncp, log = TRUE) - 
    log(sd)
}

qt.scaled <- function (p, df, mean = 0, sd = 1, ncp, lower.tail = TRUE, log.p = FALSE){
  mean + sd * stats::qt(p, df, ncp = ncp, log.p = log.p)
}


rt.scaled <- function (n, df, mean = 0, sd = 1, ncp) {
  mean + sd * stats::rt(n, df, ncp = ncp)
}
#plot(x= x.eval, y = dt((x.eval-0.6031418)/0.1232016,3))



makePoolPlot.Data <- function(pool.df,  pool_type = "linear pool", add_hist = TRUE, max_mix_eval = 3, plt_other_dists = TRUE){
  
  pool.df$dist <- stringr::str_replace_all(pool.df$dist, "normal", "norm") 
  pool.df$dist <- stringr::str_replace_all(pool.df$dist, "lognorm", "lnorm") 
 quant.vec <- t(apply(pool.df, 1, function(x){get_quant_val(
    dist = x["dist"],
    param1 = x["param1"],
    param2 = x["param2"],
    param3 = x["param3"],
    probs = c(0.001,0.025,0.5,0.975,0.999))}))
 
 central.cauchy <- mean(quant.vec[,3])#mean
 sd.cauchy <- max(apply(quant.vec,1, function(x){(x[4]-x[2])/4})) #sd
 
 min_quant <- min(quant.vec)
 max_quant <- max(quant.vec)
  
  x.eval <- seq(min_quant, max_quant, length.out = 100)
  if(pool_type == "log pool"){
    dens.eval <- eval_dens_pool(x.eval,pool.df,pool_type = "log pool")
    k_norm <- sfsmisc::integrate.xy(x = x.eval,fx = dens.eval)
    dens.eval <- dens.eval/k_norm
  }else{
    dens.eval <- eval_dens_pool(x.eval,pool.df,pool_type = "linear pool")
  }
  

  M <- max(dens.eval)
  M_x <- x.eval[which.max(dens.eval)]
  #plot(y = dens.eval, x = x.eval, type = "l")
  
  
  Nsim=100000
  u <- runif(Nsim)

  y= rt.scaled(Nsim, df =1, mean = central.cauchy, sd = sd.cauchy ) #generation from g
  
  
  if(pool_type == "log pool"){
    x=y[u<(eval_dens_pool(y,pool.df, pool_type)/k_norm)/(M*dt.scaled(y, df =1, mean = central.cauchy, sd = sd.cauchy))]
  }else{
    x=y[u<(eval_dens_pool(y,pool.df, pool_type))/(M*dt.scaled(y, df =1, mean = central.cauchy, sd = sd.cauchy))]
  }
  
  
  # fit2 <- mixfit(x, ncomp = 2) 
  # plot(fit2, ps ="ggplot2")
  
  #hist(x, breaks = seq(min(x),max(x), length.out = 200), freq = F)
  #lines(x = x.eval, y = dens.eval, col = "red")
  
  
  #Approximate log and linear pool with a Mixture of Normals
  #Fit it against the distribution as we can easily get the quantiles from a mixture of normal
  #compare against ssq and pick the best fitting distribution
  
  
  probs_pool <- seq(0.01,0.99, by =0.05)
  #Calcualte SSQ for the Mixture distribution
  quant_pool <- quantile(x, probs =probs_pool)
  
  fit.dist <- SHELF::fitdist(vals = quant_pool, probs = probs_pool)
  #See which fits best in terms of AIC.
  fit.mix.list <- list()
  
  AIC.vec <- rep(Inf,max_mix_eval)
  
    for(k in 2:max_mix_eval){ #Max number of mixture distributions
  
     fit.mix.list[[k]] <-  eval(parse(text = paste0('mixfit(x, ncomp = ',k,', family="normal")')))
     AIC.vec[k] <- fit.mix.list[[k]]$aic
    
    
    # fit.mix.list[[2]] <- mixfit(x, ncomp = 2, family= "normal")
    # AIC.vec[1] <- fit.mix.list[[2]]$aic
    # fit.mix.list[[3]] <- mixfit(x, ncomp = 3, family= "normal")
    # AIC.vec[2] <- fit.mix.list[[3]]$aic
    # 
    # fit.mix.list[[4]] <- mixfit(x, ncomp = 4, family= "normal")
    # AIC.vec[3] <- fit.mix.list[[4]]$aic
    # 
  }
  #print(AIC.vec)
  fit.mixnorm <- fit.mix.list[[which.min(AIC.vec)]]
  
  ssq_df_eval <- fit.dist$ssq[names(fit.dist$ssq) %in% c("normal", "t", "gamma", "lognormal", "beta")]
  
  finite.dists <- names(ssq_df_eval)[sapply(ssq_df_eval, is.finite)]
  
  plots.eval <- list()
  #http://www.sthda.com/english/wiki/ggplot2-colors-how-to-change-colors-automatically-and-manually
  for(i in 1:length(finite.dists)){
    plots.eval[[finite.dists[i]]] <- plotfit(fit = fit.dist, d = finite.dists[i], returnPlot = T)
  }

  plot.fit.mixnorm <- plot(fit.mixnorm, 
                           breaks =100, 
                           add_hist =add_hist)
  
  if(plt_other_dists){
    comb_df <- NULL
    
    for(i in 1:length(finite.dists)){ 
      df_temp <- plots.eval[[i]]$data
      df_temp$Dist <- finite.dists[i]
      comb_df <- rbind(comb_df, df_temp)
      
    }
    
    plot.fit.mixnorm <- plot.fit.mixnorm+
      geom_line(data = comb_df, aes(x = x, y = fx, colour = Dist), inherit.aes = F)+
      ggtitle("All Distributions")
    
  }


  #See which is the best fit
  
  ssq_mix_eval <- ssq_mix(fit.mixnorm, values = quant_pool, probs = probs_pool)
  ratio <- fit.dist$ssq[fit.dist$best.fitting[[1]]]/ssq_mix_eval
  
  if(ratio < 3){
    best_fit <- fit.dist$best.fitting[[1]]
  }else{
    best_fit <- "mixture"
  }
  
  # Extract the parameters of the best fitting distribution
  
  # Stan can't handle missing data even if it is not used.
  param_expert <- array(-999.2,dim = c(1,max_mix_eval,3)) 
  
  if(best_fit == "mixture"){
    num_param <- rep(length(fit.mixnorm$pi), 3)
    
    param_expert[1, 1:num_param[1],1] <- fit.mixnorm$mu
    param_expert[1, 1:num_param[2],2] <- fit.mixnorm$sd
    param_expert[1, 1:num_param[3],3] <- fit.mixnorm$pi
    
  }else if(best_fit == "t"){
    num_param <- rep(1,3)
    param_expert[1, 1,1] <- fit.dist[["Student.t"]][[1]]
    param_expert[1, 1,2] <- fit.dist[["Student.t"]][[2]]
    param_expert[1, 1,3] <- fit.dist[["Student.t"]][[3]]
    
  }else{
     
  best_fit_index <- which(names(fit.dist$ssq) ==  best_fit)
  param_expert[1, 1,1] <- fit.dist[[best_fit_index]][[1]]
  param_expert[1, 1,2] <- fit.dist[[best_fit_index]][[2]]
  
  }
  best_fit <- stringr::str_replace_all(best_fit, "normal", "norm") 
  best_fit <- stringr::str_replace_all(best_fit, "lognorm", "lnorm") 
  
  list_final <- list(plot.fit.mixnorm, param_expert)
  
  names(list_final) <- c("plot.fit.mixnorm",best_fit)
  
 return(list_final)
}



eval_dens_pool <- function(x.eval, pool.df, pool_type){

  #Ensure weights sum to 1
    pool.df$wi <- pool.df$wi/sum(pool.df$wi)

  dens.vec <- apply(pool.df, 1, function(x){get_density(
    dist = x["dist"],
    param1 = x["param1"],
    param2 = x["param2"],
    param3 = x["param3"],
    x = x.eval)})

  if(pool_type == "log pool"){
    
    if(is.matrix(dens.vec)){
      return(apply(dens.vec, 1, function(x){prod(x^pool.df$wi)}))
    }else{
      #scaled by an arbitary constant which we don't need to know for candidate density evaluation
      return(prod(dens.vec^pool.df$wi))
    }
    
  }else{
    if(is.matrix(dens.vec)){
      return(apply(dens.vec, 1, function(x){sum(x*pool.df$wi)}))
      
    }else{
      
      return(sum(dens.vec*pool.df$wi))
    }
    
  }
  
}


ssq_mix <- function(object, values, probs){
  df_ssq <- data.frame(pi = object$pi, mu = object$mu, sd = object$sd)
  
  #Evaluate the pnorm individually
  pnorm_eval  <- apply(df_ssq,1, FUN = function(x){pnorm(values,x["mu"],
                                                         x["sd"])})
  pnorm_eval_weighted <- t(pnorm_eval)*df_ssq$pi
  
  #Sum the pnorm then subtract
  return(sum((colSums(pnorm_eval_weighted) - probs)^2))
  
}


makePoolPlot2 <- function(pool.df, x.eval, pool_type = "linear pool"){
  
  if(pool_type == "log pool"){
    dens.eval <- eval_dens_pool(x.eval,pool.df,pool_type = "log pool")
    k_norm<- sfsmisc::integrate.xy(x = x.eval,fx = dens.eval)
    
    dens.eval <- dens.eval/k_norm
  }else{
    dens.eval <- eval_dens_pool(x.eval,pool.df,pool_type = "linear pool")
  }
  
  M <- max(dens.eval)
  M_x <- x.eval[which.max(dens.eval)]
  plot(y = dens.eval, x = x.eval, type = "l")
  
  
  Nsim=250000
  #u=runif(Nsim,max=M) #uniform over (0,M)
  u <- runif(Nsim)
  y=runif(Nsim) #generation from g
  
  if(pool_type == "log pool"){
    x=y[u<(eval_dens_pool(y,pool.df, pool_type)/k_norm)/(M*1)]
  }else{
    x=y[u<(eval_dens_pool(y,pool.df, pool_type))/(M*1)]
  }
  
  
  # fit2 <- mixfit(x, ncomp = 2) 
  # plot(fit2, ps ="ggplot2")
  
  hist(x, breaks = seq(0,1, length.out = 500), freq = F)
  lines(x = x.eval, y = dens.eval, col = "red")
  
  #saveRDS(data.stan.example, file = "C:/Users/phili/Desktop/stanfile.RDS")
  #saveRDS(data.stan.example, file = "C:/Users/phili/Desktop/stanfile.RDS")
  #library(rstan)
  #data.stan.example<- readRDS(file = "C:/Users/phili/Desktop/stanfile.RDS")
  
  
  #Approximate log and linear pool with a Mixture of Normals
  #Fit it against the distribution as we can easily get the quantiles from a mixture of normal
  #compare against ssq and pick the best fitting distribution
  
  
  probs_pool <- seq(0.01,0.99, by =0.05)
  #Calcualte SSQ for the Mixture distribution
  quant_pool <- quantile(x, probs =probs_pool)
  
  fit.dist <- SHELF::fitdist(vals = quant_pool, probs = probs_pool, lower = 0, upper = 1)
  
  fit.mixnorm <- mixfit(x, ncomp = 2, family= "normal")
  
  
  plot.fit.norm <- plotfit(fit = fit.dist, d = "normal", returnPlot = T)
  plot.fit.student_t <- plotfit(fit = fit.dist, d = "t", returnPlot = T)
  plot.fit.gamma <- plotfit(fit = fit.dist, d = "gamma", returnPlot = T)
  plot.fit.lognormal <- plotfit(fit = fit.dist, d = "lognormal", returnPlot = T)
  plot.fit.beta <- plotfit(fit = fit.dist, d = "beta", returnPlot = T)
  
  plot.fit.mixnorm <- plot(fit.mixnorm)
  
  plot.fit.mixnorm+
    ggtitle("All Distributions")+
    geom_line(data = plot.fit.norm$data, aes(x,fx, colour = "red")) +
    geom_line(data = plot.fit.student_t$data, aes(x,fx, colour = "green"))+
    geom_line(data = plot.fit.gamma$data, aes(x,fx, colour = "orange"))+
    geom_line(data = plot.fit.lognormal$data, aes(x,fx, colour = "blue"))+
    geom_line(data = plot.fit.beta$data, aes(x,fx, colour = "purple"))+
    scale_colour_manual(name = 'Dists', 
                        values =c('red'= 'red','green'='green',
                                  "orange"="orange","blue"="blue", "purple"="purple"),
                        labels = c('normal','t', "gamma", "lnorm", "beta")) 
  
  
}



#' Title
#'
#' @param expert_df 
#'
#' @return
#' @export
#'
#' @examples
expert_dens <- function(expert_df, probs =  seq(0.01, 0.98, by = 0.002)){
  
if(length(unique(expert_df$expert)) !=1){ #Only one expert, Don't need to anything
  
  
  if(is.null(expert_df$weights) && is.null(expert_df$wi)){
    warning("No weights given.. assuming equally weighted expert opinion")
    expert_df$weights <- 1
  }
  
  if(!is.null(expert_df$wi)){
    expert_df$weights <- expert_df$wi
  }
  
  
  
  expert_df_sum <- expert_df %>% group_by(times_expert) %>% arrange(times_expert) %>%
                  summarize(sum_weights = sum(weights))
  
  if(any(expert_df_sum$sum_weights != 1)&& any(expert_df$weights != 1)){
    warning("Some weights don't sum to 1.. reweighting")
    
  }
  expert_df <-expert_df %>% left_join(expert_df_sum,"times_expert") %>%
    mutate(weights = weights/sum_weights)
  
  }
  
  expert_density <- apply(expert_df, 1, function(x){get_quant_val(
    dist = x["dist"],
    param1 = x["param1"],
    param2 = x["param2"],
    param3 = x["param3"],
    probs = probs)})
  
  rownames(expert_density) <- probs
  
  list(expert_df = expert_df  %>% dplyr::select(-sum_weights),
       expert_density = expert_density)
  
}


#expert_density <- expert_dens(expert_df)



expert_pooling <- function(expert_density = NULL, expert_quant_list = NULL,
                           lower_bound = -Inf, upper_bound = Inf){

dfs_expert <- list() 
plts_pool <- list()
dfs_pool <- list()


if(!is.null(expert_quant_list)){ # If density

 if(is.null(expert_quant_list$weights_mat)){
   weights_mat <- NULL
 }
  suppress_messages(attach(expert_quant_list))
  
    
max.timepoints  <- length(times)

for(i in 1:max.timepoints){
  
  timepoint <- paste0("Time ",times[i])
  
  fit.eval <- SHELF::fitdist(vals = na.omit(v_array[,,i]),
                      probs = na.omit(p_mat[,i]), lower = lower_bound, upper = upper_bound)
  
  weights <- na.omit(weights_mat[,i])
  
  if(is.null(weights_mat) && ncol(na.omit(v_array[,,i])) == 1){
    weights <- 1 #Only one expert so weights should be 1
  }else if(is.null(weights_mat)){
    warning("No weights assigned assuming equal weights")
    weights <- 1
  }
  
  best_fit_index  <- apply(fit.eval$ssq[,c("normal","t","gamma", "lognormal", "beta")], 1, which.min)
  best_fit <- names(fit.eval$ssq[,c("normal","t","gamma", "lognormal", "beta")])[best_fit_index]
  
  best_fit_loc <- sapply(best_fit, function(x){which(x  == names(fit.eval$ssq))})
  fit.eval.dist  <- fit.eval[best_fit_loc]
  
  pool.df_output <- matrix(nrow = length(best_fit_loc),ncol = 3)
  colnames(pool.df_output) <- c("param1", "param2", "param3")
  
  for(i in 1:length(best_fit_loc)){
    pool.df_output[i,1:length(fit.eval.dist[[i]][i,])] <-  as.numeric(as.vector(fit.eval.dist[[i]][i,]))
  }
  dfs_expert[[timepoint]] <- data.frame(dist = best_fit, wi = weights, pool.df_output)

  plts_pool[[timepoint]] <- makePoolPlot(fit  = fit.eval,
                                         xl =lower_bound,
                                         xu =upper_bound,
                                         d = "best",
                                         w = weights,
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

  dfs_pool[[timepoint]] <-  plts_pool[[timepoint]][["data"]]

  }

}else{

  times <- unique(expert_density$expert_df[,"times_expert"])
  probs <- as.numeric(rownames(expert_density$expert_density))
  
for(i in 1:length(times)){
  
  timepoint <- paste0("Time ",times[i])
  
  index.loc <- which(expert_density$expert_df$times_expert == times[i])
  temp_df <- expert_density$expert_df[index.loc, ]
  temp_dens <- expert_density$expert_density[,index.loc]
  
  v <- temp_dens
  p <- matrix(rep(probs, ncol(temp_dens)), nrow = length(probs), ncol = ncol(temp_dens))
  
  # Need to consider upper and lower bounds
  fit.eval <- SHELF::fitdist(v, p, lower= lower_bound, upper = upper_bound)
  
  if(!is.null(temp_df$weights)){
    weights <- temp_df$weights
  }else{
    weights <- 1
  } 
  
  
  best_fit_index  <- apply(fit.eval$ssq[,c("normal","t","gamma", "lognormal", "beta")], 1, which.min)
  best_fit <- names(fit.eval$ssq[,c("normal","t","gamma", "lognormal", "beta")])[best_fit_index]
  
  best_fit_loc <- sapply(best_fit, function(x){which(x  == names(fit.eval$ssq))})
  fit.eval.dist  <- fit.eval[best_fit_loc]
  
  pool.df_output <- matrix(nrow = length(best_fit_loc),ncol = 3)
  colnames(pool.df_output) <- c("param1", "param2", "param3")
  
  for(i in 1:length(best_fit_loc)){
    pool.df_output[i,1:length(fit.eval.dist[[i]][i,])] <-  as.numeric(as.vector(fit.eval.dist[[i]][i,]))
  }
  dfs_expert[[timepoint]] <- data.frame(dist = best_fit, wi = weights, pool.df_output)
  
  
  plts_pool[[timepoint]] <- makePoolPlot(fit  = fit.eval,
                            xl =lower_bound,
                            xu =upper_bound,
                            d = "best",
                            w = weights,
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
  
  dfs_pool[[timepoint]] <-  plts_pool[[timepoint]][["data"]]
  
  
  }
  
}
 list(dfs_expert =dfs_expert,
      plts_pool =plts_pool,
      dfs_pool = dfs_pool)

}

#myfit <- fitdist(expert_density, probs)

#plotfit(myfit,lp = T)

# Reweights the weights if they don't sum to 1 anyway



get_quant_val <- function(dist,param1, param2, param3 = NULL, probs = seq(0.01, 0.98, by = 0.01)){
  if(dist == "t"){
    probs_eval <- as.numeric(param1) + as.numeric(param2)*stats::qt(as.numeric(probs),as.numeric(param3))
    return(probs_eval)
     
  }else{
    probs <- paste0(probs, collapse = ",")
    
    probs_eval <-  paste0("q",dist,
                          "(c(",probs,"),", param1,
                          ",",param2,")")
    probs_eval <- eval(parse(text = probs_eval))
    return(probs_eval)
  }
 
}

pt.scaled <-function (q, df, mean = 0, sd = 1, ncp, lower.tail = TRUE, log.p = FALSE){
  stats::pt((q - mean)/sd, df, ncp = ncp, log.p = log.p)
}

get_cdf_val <- function(dist,param1, param2, param3 = NULL, vals = seq(0.01, 0.98, by = 0.01)){
  if(dist == "t"){
    probs_eval <- stats::pt((vals -  as.numeric(param1))/as.numeric(param2), as.numeric(param3), log.p = F)
    return(probs_eval)
    
  }else{
    vals <- paste0(vals, collapse = ",")
    
    probs_eval <-  paste0("p",dist,
                          "(c(",vals,"),", param1,
                          ",",param2,")")
    probs_eval <- eval(parse(text = probs_eval))
    return(probs_eval)
  }
  
}


get_density <- function(dist, param1, param2, param3 = NULL, x = seq(0.01, 0.98, by = 0.01)){
  x <- paste0(x, collapse = ",")
  if(dist == "t"){
    #From SHELF reference Student.t Parameters of the fitted t distributions. 
    #Note that (X - location) / scale has a standard t distribution 
    dens_x <-  paste0("d",dist,
                      "((c(",x,")-",param1,")/",param2,",", param3,")/",param2)
     }else{ #
    dens_x <-  paste0("d",dist,
                      "(c(",x,"),", param1,
                      ",",param2,")")
    
  }

  dens_eval <- eval(parse(text = dens_x))
  
  return(dens_eval)
}


#' Credible interval for pooled distribution
#'
#' Returns the interval based on defined quantiles. 
#' The approach used only provides an approximate (although quite accurate) integral.  
#' @param plt_obj A plot object from `plot_expert_opinion`.
#' @param val The name of the opinion for which the interval will be generated.
#' @param interval A vector of the upper and lower probabilities. Default is the standard 95% interval 
#'
#' @return
#' @export
#'
#' @examples \dontrun{
#' cred_int(plot_opinion1,val = "linear pool", interval = c(0.025, 0.975))
#' }
#' 
cred_int <- function(plt_obj, val = "linear pool",interval = c(0.025, 0.975)){
  
  plt_df <- plt_obj$data %>% filter(expert == val) %>% data.frame()
  
  total_integral <- sfsmisc::integrate.xy(plt_df$x, plt_df$fx)
  partial_integral <- rep(NA, nrow(plt_df))
  partial_integral[1] <- 0
  for(i in 2:nrow(plt_df)){
    partial_integral[i] <- sfsmisc::integrate.xy(plt_df$x[1:i], plt_df$fx[1:i])/total_integral
  }
  
  plt_df$cdf <- partial_integral
  
  limits <- c(plt_df$x[which.min(abs(plt_df$cdf - interval[1]))],plt_df$x[which.min(abs(plt_df$cdf - interval[2]))])
  names(limits) <- c("lower", "upper")
  return(limits)
  
}


makePoolPlot <- function (fit, xl, xu, d = "best", w = 1, lwd =1, xlab="x", 
                          ylab=expression(f[X](x)), legend_full = TRUE, 
                          ql = NULL, qu = NULL, nx = 500, addquantile = FALSE, fs = 12, 
                          expertnames = NULL){
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
  sd.norm <- rep(NA,n.experts )
  for (i in 1:n.experts) {
    
  }
  
  if(is.infinite(xl)||is.infinite(xu)){
    max.sd.index <- which.max(fit$Normal$sd)
    
    xl <- qnorm(0.001, fit$Normal[max.sd.index, 1], fit$Normal[max.sd.index,2])
    
    xu <- qnorm(0.999, fit$Normal[max.sd.index, 1], fit$Normal[max.sd.index,2])
    
  }
  

  
  for (i in 1:n.experts) {
    densitydata <- SHELF:::expertdensity(fit, d[i], ex = i, xl, 
                                         xu, ql, qu, nx)
    x[, i] <- densitydata$x
    fx[, i] <- densitydata$fx
  }
  #plot(x =  x[, 1], y = fx[,1])
  #lines(x =  x[, 2], y = fx[,2])
  
  fx.lp <- apply(fx * weight, 1, sum)
  if(any(is.infinite(fx ^ weight))){
    warning("Print Non finite density for log pooling - Results invalid")
  }
  fx.logp <- apply(fx ^ weight, 1, prod)
  k_norm <- sfsmisc::integrate.xy(x = x[,1], fx = fx.logp)
  fx.logp <- fx.logp/k_norm
  
  #Should integrate to 1 after normalization
  #sfsmisc::integrate.xy(x = x[,1], fx = fx.logp)
  
  
  df1 <- data.frame(x = rep(x[, 1], n.experts + 2),
                    fx = c(as.numeric(fx),fx.lp, fx.logp), 
                    expert = factor(c(rep(expertnames, each = nxTotal),
                                      rep("linear pool", nxTotal),rep("log pool", nxTotal)),
                                    levels = c(expertnames, "linear pool", "log pool")), 
                    ftype = factor(c(rep("individual", nxTotal * n.experts), 
                                     rep("linear pool", nxTotal),rep("log pool", nxTotal)),
                                   levels = c("individual","linear pool", "log pool")))
  
  df1$expert <- factor(df1$expert, levels = c(expertnames, 
                                              "linear pool", "log pool"))
  if (legend_full) {
    cols <- (scales::hue_pal())(n.experts + 2)
    linetypes <- c(rep("dashed", n.experts), "solid","solid")
    sizes <- lwd * c(rep(0.5, n.experts), 1.5,1.5)
    names(cols) <- names(linetypes) <- names(sizes) <- c(expertnames, 
                                                         lpname)
    p1 <- ggplot(df1, aes(x = x, y = fx, colour = expert, 
                          linetype = expert, size = expert)) + scale_colour_manual(values = cols, 
                                                                                   breaks = c(expertnames, lpname)) + scale_linetype_manual(values = linetypes, 
                                                                                                                                            breaks = c(expertnames, lpname)) + scale_size_manual(values = sizes, 
                                                                                                                                                                                                 breaks = c(expertnames, lpname))
  }else {
    p1 <- ggplot(df1, aes(x = x, y = fx, colour = ftype, 
                          linetype = ftype, size = ftype)) + scale_linetype_manual(name = "distribution", 
                                                                                   values = c("dashed", "solid","solid")) + scale_size_manual(name = "distribution", 
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
  } else {
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
  }else{
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
                                                   x <= ql & expert == lpname[1])), aes(ymax = fx, ymin = 0), 
                           alpha = 0.2, show.legend = FALSE, colour = NA, fill = ribbon_col) + 
      geom_ribbon(data = with(df1, subset(df1, x >= qu & 
                                            expert == lpname[2])), aes(ymax = fx, ymin = 0), 
                  alpha = 0.2, show.legend = FALSE, colour = NA, 
                  fill = ribbon_col)
    
    p1 <- p1 + geom_ribbon(data = with(df1, subset(df1, 
                                                   x <= ql & expert == lpname[2])), aes(ymax = fx, ymin = 0), 
                           alpha = 0.2, show.legend = FALSE, colour = NA, fill = ribbon_col) + 
      geom_ribbon(data = with(df1, subset(df1, x >= qu & 
                                            expert == lpname[2])), aes(ymax = fx, ymin = 0), 
                  alpha = 0.2, show.legend = FALSE, colour = NA, 
                  fill = ribbon_col)
    
    
  }
  if (lpname[1] == "marginal") {
    p1 <- p1 + theme(legend.title = element_blank())
  }
  p1 + theme(text = element_text(size = fs))
}




#' Plotting Pooled Expert Opinion
#'
#' Returns a ggplot with the individual expert opinions along with the pooled distributions (both linear and logarithmic).
#'
#' @param object Either a object of class elicitation (from `SHELF`) or a dataframe with parameters of the distribution (see Example below).
#' @param xl_plt Optionally set the lower bound for the plot
#' @param xu_plt Optionally set the upper bound for the plot
#' @param weights A vector with the weight of each expert. If omitted, set to equal weights.
#'
#' @return
#' @export
#'
#' @examples \dontrun{ 
#'  expert_df <- data.frame(dist = c("norm","t"), #Distribution Name
#'                          wi = c(1/3,2/3), #Expert weights
#'                          param1 = c(0.3,0.40), #Parameter 1
#'                          param2 = c(0.05,0.05),# Parameter 2
#'                          param3 = c(NA,3)) #Parameter 3: Only t-distribution
#'  plot_expert_opinion(expert_df , weights = expert_df$wi)}
#'                                                         
plot_expert_opinion <- function(object, xl_plt = NULL, xu_plt = NULL, weights = NULL){
  
  
  if(is.null(weights)){
    weights <- 1
  }
  
  
  
  if(class(object) == "elicitation"){
    
    if(is.null(xl_plt)){
      xl_plt <- min(object$limits["lower"])
    }
    if(is.null(xu_plt)){
      xu_plt <- max(object$limits["upper"])
      
    }
    
    plt <- expertsurv:::makePoolPlot(fit= object,
                                     xl =xl_plt,
                                     xu =xu_plt,
                                     d = "best",
                                     w = weights,
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
    
    
  }else{
    
    object$times_expert <- 2 #Just for compatibility
    
    expert_dens_list <- expertsurv:::expert_dens(object, probs =  seq(0.001, 0.99, by = 0.005))
    
    lower <- as.numeric(head(expert_dens_list$expert_density, n = 1)-0.1)
    upper <- as.numeric(tail(expert_dens_list$expert_density, n = 1)+0.1)
    
    # if(is.null(lower) || is.null(upper)){
    #   stop("Upper and lower bounds required for distributions")
    # }
    
    if(is.null(xl_plt)){
      xl_plt <- min(lower)
    }
    if(is.null(xu_plt)){
      xu_plt <- max(upper)
      
    }
    
    
    
    probs_mat <- matrix(as.numeric(rep(rownames(expert_dens_list$expert_density), 
                                       dim(expert_dens_list$expert_density)[2])),
                        ncol = dim(expert_dens_list$expert_density)[2])
    
    fit_shelf  <- SHELF::fitdist(vals = expert_dens_list$expert_density,
                          probs_mat, lower = lower, upper = upper)
    
    plt <- expertsurv:::makePoolPlot(fit= fit_shelf,
                                     xl = xl_plt,
                                     xu = xu_plt,
                                     d = "best",
                                     w = weights,
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
    
    
  }
  
  return(plt+theme_bw())
}


#' Fitting Parametric Survival models with Expert Opinion
#'
#' Implementation of survival models with expert opinion on the survival probabilities or expected difference in survival.
#' Function is equivalent to the `fit.models` in \texttt{survHE} expect for the inclusion of the "expert_type" and "param_expert" arguments. 
#' Worked examples can be found in the [README](https://github.com/Philip-Cooney/expertsurv/blob/master/README.md) file.
#' Note that the default method is "hmc", however, the user may use "mle" or "inla" for analysis without expert opinion.
#'
#' @param formula As per `fit.models` on \texttt{survHE}
#' @param data As per `fit.models` on \texttt{survHE}
#' @param distr As per `fit.models` on \texttt{survHE}. Note Generalized F model is not available for method = "hmc".
#' @param method As per `fit.models` on \texttt{survHE}
#' @param expert_type Either "survival", which indicates expert opinion on the survival function or "mean" (actually anything that does not contain "survival") which represents a belief on difference in survival.
#' @param param_expert A list containing a dataframe for each timepoint (if applicable). Each dataframe should have columns with the following names and each row representing an expert:
#'  \itemize{
#'   \item \strong{dist}: Names of the distribution assigned to each expert which may be "norm", "t", "lnorm", "gamma", "beta".
#'   \item \strong{wi}: Weight of the expert, if all experts = 1 then equal weights.
#'   \item \strong{param1}: First parameter of the distribution (e.g. mean for norm distribution). Parameters as per \texttt{SHELF} package. 
#'   \item \strong{param2}: Second parameter of the distribution.
#'   \item \strong{param3}: Third parameter of the distribution (NA expect for degrees of freedom for t distribution)
#' }
#' @param ... Other arguments may be required depending on the example. See [README](https://github.com/Philip-Cooney/expertsurv/blob/master/README.md) for details.
#'
#' @return
#' @import survHE
#' @importFrom magrittr %>%
#' @export
#' @md
#' 
#' @examples
fit.models.expert <- function(formula = NULL, data, distr = NULL, method = "hmc", 
                        expert_type = "survival", param_expert = NULL, ...){
  exArgs <- list(...)
  #will need to be modified
  exArgs$formula <- formula
  exArgs$data = data
  exArgs$param_expert <- param_expert
  if(!is.null(expert_type) && method != "hmc"){
    print("Expert Opinion is only implemented with hmc method")
    stop()
  }
  
  if(!is.null(expert_type) && is.null(param_expert)){
    print("You have not specified any expert opinions using the param_expert argument")
    stop()
  }
  
  if(!is.null(expert_type) && expert_type != "survival" && any(distr == "rps")){
    print("Mean Difference is not implemented for RPS models")
    stop()
  }
  
  if(method == "hmc" && any(distr == "genf")){
    print("Generalized F models are not with expert opinion")
    stop()
  }
  
 
   fit.models(formula = formula, data = data, distr = distr, method = method, exArgs = exArgs)
}





#' Title
#'
#' @param formula 
#' @param data 
#' @param distr 
#' @param method 
#' @param ... 
#'
#' @return
#' @import survHE
#'
#' @examples
fit.models <- function (formula = NULL, data, distr = NULL, method = "mle", exArgs, 
          ...){

  if (is.null(formula)) {
    stop("You need to specify a model 'formula', e.g. 'formula=Surv(time,event)~treat'")
  }
  method <- tolower(method)
  if (!method %in% c("hmc", "inla", "mle")) {
    stop("Methods available for use are 'mle', 'hmc' or 'inla'")
  }
  survHE:::check_distributions(method, distr)
  if (method == "mle") {
    cat(crayon::blue("No expert opinion implemented with mle approach \n"))
    
    res <- survHE:::format_output_fit.models(lapply(distr, function(x) survHE:::runMLE(x, 
                                                                     exArgs)), method, distr, formula, data)
  }
  if (method == "inla") {
    
    cat(crayon::blue("No expert opinion implemented with inla approach \n"))
    
    res <- survHE:::format_output_fit.models(lapply(distr, function(x) survHE:::runINLA(x, 
                                                                      exArgs)), method, distr, formula, data)
  }
  if (method == "hmc") {
    
    adjust_survHE_func()#Make sure the survHE formulaes are compatible
    res <- format_output_fit.models(lapply(distr, function(x) runHMC(x, 
                                                                     exArgs)), method, distr, formula, data)
  }
  return(res)
}

runHMC <- function (x, exArgs){
  if (!isTRUE(requireNamespace("rstan", quietly = TRUE))) {
    stop("You need to install the R package 'rstan'. Please run in your R terminal:\n install.packages('rstan')")
  }
  formula <- exArgs$formula
  data = exArgs$data
  availables <- survHE:::load_availables()
  d3 <- survHE:::manipulate_distributions(x)$distr3
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
    if(d3%in% names(exArgs$init) ){
      init <- exArgs$init[[d3]] 
    }else{
      init = "random"
    }
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
  data.stan <- make_data_stan(formula, data, d3, exArgs)
  
  tic <- proc.time()
  
  if (d3 %in% c("gam", "gga", "gom")){
    data.jags <- data.stan
    if(d3 %in% c( "gom")){
      parameters.to.save_jags = c("alpha","beta", "rate")
      
      #Inits as per flexsurvreg (reparameterized)
      modelinits <- function(){
        beta = c(log(1/mean(data.jags$t)*runif(1,0.8,1.5)),rep(0,data.jags$H -1))
        list(alpha1 = runif(1,0.001,0.003),alpha2 = runif(1,0.001,0.003), beta = beta) 
      }
       
    }else if(d3 == "gga"){ #(d3 == "gga")
      parameters.to.save_jags = c("Q","sigma", "beta", "r", "b","mu")
      tinits1 <-data.jags$t + max(data.jags$t)
      is.na(tinits1)<-data.jags$d ==1
      data.jags$is.censored <- ifelse(data.jags$d==0, 1, 0)
      data.jags$t_jags <- ifelse(data.jags$is.censored ==1, NA, data.jags$t) 
      data.jags$t_cen <- data.jags$t+data.jags$d
      modelinits <- function(){list(t_jags = tinits1)}
      #Stop JAGS Warning messages
      data.jags <- data.jags[names(data.jags) %!in% c("t", "d", "a0")]
      
      
    }else{ #"gam",
      parameters.to.save_jags = c("alpha","beta", "rate")
      modelinits <- NULL
    }
    data.jags <- data.jags[names(data.jags) %!in% "max_param"]
    cat(paste0(" \n SAMPLING FOR MODEL '",d,"_expert' NOW.  \n"))
    model <-R2jags::jags(model.file = textConnection(get(paste0(d,".jags"))),
                             data=data.jags,
                             n.chains=chains,
                             inits=modelinits,
                             parameters.to.save = c(parameters.to.save_jags,"St_expert"),
                             n.iter = iter*5,
                             n.thin = thin,
                             n.burnin = iter,
                             jags.module = c("glm","dic"))

    
  }else{
    dso <- stanmodels[[paste0(d, "_expert")]]
    model <- rstan::sampling(dso, data.stan, chains = chains, 
                             iter = iter, warmup = warmup, thin = thin, seed = seed, 
                             control = control, pars = pars, include = include, cores = cores, 
                             init = init, refresh = refresh)
    
     time_stan <- sum(rstan::get_elapsed_time(model))
    
  }
  
  toc <- proc.time() - tic
  time_survHE <- toc[3]
  ics <- compute_ICs_stan(model, d3, data.stan)
  
  if (save.stan) {
    if(d3 %in% c("gam", "gga", "gom")){
      
      model_code <- get(paste0(d,".jags"))
      con <- paste0(d, ".txt")
    }else{
      model_code <- attr(model@stanmodel, "model_code")
      con <- paste0(d, ".stan")

    }
    
    writeLines(model_code, con = con)
    cat(paste0("Model code saved to the file: ", con, 
               "\n"))
    
    ## Add in for Jags
  }
  model_name <- d3
  list(model = model, aic = ics$aic, bic = ics$bic, dic = ics$dic, 
       dic2 = ics$dic2,waic = ics$waic, pml = ics$pml,  time2run = time_survHE, 
       data.stan = data.stan, save.stan = save.stan, model_name = model_name)
}




make_data_stan <- function (formula, data, distr3, exArgs = globalenv()){
  availables <- survHE:::load_availables()
  method <- "hmc"
  formula_temp <- update(formula, paste(all.vars(formula, data)[1], 
                                        "~", all.vars(formula, data)[2], "+."))
  mf <- as_tibble(model.frame(formula_temp, data)) %>% 
          dplyr::rename(time = 1,event = 2) %>% rename_if(is.factor, .funs = ~gsub("as.factor[( )]","", .x)) %>% 
          dplyr::rename_if(is.factor, .funs = ~gsub("[( )]","", .x)) %>% 
          bind_cols(as_tibble(model.matrix(formula_temp,data)) %>% dplyr::select(contains("Intercept"))) %>%
          dplyr::select(time,event, contains("Intercept"), everything()) %>% tibble::rownames_to_column("ID")
  
  ####Code Change Here
  
  ######
  #print(data.stan$a0_obs)

  if (distr3 %!in% c("rps")) {
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
    B_expert <- flexsurv::basis(knots, log(exArgs$times_expert))
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
  #if (distr3 %in% c("gga", "lno", "gam")) {
  if (distr3 %in% c("lno")) {
    
      data.stan$sigma_beta <- rep(100, data.stan$H)
  }
  #else {
    data.stan$sigma_beta <- rep(5, data.stan$H)
  #}
  # if (distr3 == "gef") {
  #   data.stan$a_sigma = data.stan$b_sigma = 0.1
  #   data.stan$mu_P = 0
  #   data.stan$sigma_P = 0.5
  #   data.stan$mu_Q = 0
  #   data.stan$sigma_Q = 2.5
  # }
  # else if (distr3 == "gga") {
  #   data.stan$a_sigma = data.stan$b_sigma = 0.1
  #   data.stan$mu_Q = 0
  #   data.stan$sigma_Q = 100
  # }
  if (distr3 %in% c("gam","gom", "gga", "llo", "wei", 
                         "wph")) {
    data.stan$a_alpha = data.stan$b_alpha = 0.1
  }else if(distr3 %in% c("lno")){
    data.stan$a_alpha = 0
    data.stan$b_alpha = 5
  }
  d <- names(availables[[method]][match(distr3, availables[[method]])])
  priors <- list()
  if (exists("priors", where = exArgs)) {
    abbrs = survHE:::manipulate_distributions(names(exArgs$priors))$distr3
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
  
  
  if(exArgs$opinion_type == "survival"){
    data.stan$St_indic <- 1
    #even if survival need to define these (just put as 1)
    data.stan$id_comp <- 1
    data.stan$id_trt <- 1
  }else{
    data.stan$St_indic <- 0
    #even if survival need to define these (just put as 1)
    data.stan$id_St <- 1
   }

  if(ncol(mf) == 4){
    #No covariates
    # Has to be opinion_type survival 
    data.stan$id_St <- 1
    
  }else if(ncol(mf) == 5){
    
    if(exArgs$opinion_type == "survival"){
      data.stan$id_St <- min(which(mf[,5] == exArgs$id_St))
    }else{# Survival Difference
      data.stan$id_trt <- min(which(mf[,5] == exArgs$id_trt)) 
      if(length(unique(mf[,5] %>% pull()))==2){
        data.stan$id_comp <- min(which(mf[,5] != exArgs$id_trt)) 
      }else{
        data.stan$id_comp <- min(which(mf[,5] == exArgs$id_comp))  
      }
      
    }
    #put the number in  could put in a combination of numbers
  }else{
    print("We do not allow more than one covariate (i.e. treatment) in the analysis")
    stop()
  }
  
  
  
  # 
  # 
  # if(exArgs$opinion_type == "survival"){
  #   if(length(unique(data.stan[["X"]][,2])) == 1){
  #     data.stan$id_St <-1
  #   }else{
  #     #Need to work on this function
  #     data.stan$id_St <-  min(which(data.stan[["X"]][,2] == exArgs$id_St))
  #   }
  # 
  #     
  #   data.stan$id_comp <-  1
  #   data.stan$id_trt <-  1
  #   # if (distr3 %in% c("gam", "gga", "gom")){
  #   #   data.stan$id_comp <-  1  # Has to be defined for JAGS
  #   #   data.stan$id_trt <-  1
  #   #   
  #   # }
  #   
  # }else{
  #   
  #     data.stan$id_trt <-  min(which(data.stan[["X"]][,2] == 1))
  #     data.stan$id_comp <-  min(which(data.stan[["X"]][,2] == 0)) 
  #     data.stan$id_St <-  numeric(0)
  #     if (distr3 %in% c("gam", "gga", "gom")){
  #       data.stan$id_St <-  1  # Has to be defined for JAGS
  #       
  #     }
  # 
  # }
    
   
    param_expert <- exArgs$param_expert
    n.experts <- c()
    
    for(i in 1:length(param_expert)){
      n.experts <- c(n.experts, nrow(param_expert[[i]])) 
    }
    
    data_dist_ind <- num_param <- matrix(-999.2,nrow = max(n.experts), ncol =  length(param_expert))
    expert.array <- array(-999.2,dim = c(max(n.experts),5,length(param_expert))) 
    
    for(i in 1:length(param_expert)){
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

    if(!is.null(exArgs$times_expert)){
      data.stan$n_time_expert <- length(exArgs$times_expert)
      data.stan$time_expert <- as.array(exArgs$times_expert)
    }else{
      data.stan$n_time_expert <- 1
      data.stan$time_expert <- numeric(0) #This produces an array of size 0
      #https://dev.to/martinmodrak/optional-parametersdata-in-stan-4o33
      
      if (distr3 %in% c("gam", "gga", "gom")){
        data.stan$time_expert <- 1 # Has to be defined for JAGS
        
      }
      
      
    }
 
    data.stan$param_expert <-expert.array
    data.stan$n_experts <- as.array(n.experts)  
    
    if(is.null(exArgs$pool_type)){
     
      data.stan$pool_type <- 1
      
      }else{
        data.stan$pool_type <- as.numeric(grepl("line", exArgs$pool_type)) 
        
      }
    
    if(data.stan$pool_type == 0){
      k_norm <- rep(NA,data.stan$n_time_expert )
      for(i in 1:data.stan$n_time_expert){
        
        param_expert[[i]]$dist <- stringr::str_replace_all(param_expert[[i]]$dist, "normal", "norm") 
        param_expert[[i]]$dist <- stringr::str_replace_all(param_expert[[i]]$dist, "lognorm", "lnorm") 
        quant.vec <- t(apply(param_expert[[i]], 1, function(x){get_quant_val(
          dist = x["dist"],
          param1 = x["param1"],
          param2 = x["param2"],
          param3 = x["param3"],
          probs = c(0.001,0.025,0.5,0.975,0.999))}))
        
        central.cauchy <- mean(quant.vec[,3])#mean
        sd.cauchy <- max(apply(quant.vec,1, function(x){(x[4]-x[2])/4})) #sd
        
        min_quant <- min(quant.vec)
        max_quant <- max(quant.vec)
        
        x.eval <- seq(min_quant, max_quant, length.out = 100)
          dens.eval <- eval_dens_pool(x.eval,param_expert[[i]],pool_type = "log pool")
          k_norm[i] <- sfsmisc::integrate.xy(x = x.eval,fx = dens.eval)
          #dens.eval <- dens.eval/k_norm
       
      }
      data.stan$k_norm <- k_norm

    }
    
  

    #Power prior
   
    if(!is.null(exArgs$a0)){
      data.stan$a0 <- exArgs$a0
    }else{
      data.stan$a0 <- rep(1, nrow(data))
    }
 
  
  #data.stan$exArgs <- exArgs
  #save(data.stan, file = paste0(pathway, "Extra Output/datastan.RData"))
  #
  data.stan	
}

# model <- fit.weibull
# distr3 <- "weib"
# data.stan <- stan.data

compute_ICs_stan <-function (model, distr3, data.stan){
  if (distr3 %!in% c("gam", "gga", "gom")) {
    beta <- rstan::extract(model)$beta
  }
  else {
    beta <- model$BUGSoutput$sims.matrix[, grep("beta", 
                                                colnames(model$BUGSoutput$sims.matrix))]
  }
  beta.hat <- apply(beta, 2, mean)
  linpred <- beta %*% t(data.stan$X)
  linpred.hat <- beta.hat %*% t(data.stan$X)
  model.eval <- paste0("lik_", distr3)
  out = do.call(what = eval(parse(text = model.eval)), args = list(distr3, 
                                                                   linpred, linpred.hat, model, data.stan))
  logf = out$logf
  logf.hat = out$logf.hat
  npars = out$npars
  logf_comb <- matrix(nrow = nrow(logf), ncol = ncol(logf))
  for (i in 1:nrow(logf)) {
    logf_comb[i, ] <- logf[i, ] + out$logf.expert[i]/ncol(logf)
  }
  tryCatch(suppressWarnings(WAIC <- loo::loo(logf_comb)[["estimates"]][grep("looic", 
                                                                   rownames(loo::loo(logf_comb)[["estimates"]])), "Estimate"]),
           error = function(e)
             print("Cannot Evaluate WAIC"))
  
  if(!exists("WAIC")){
    WAIC <- Inf
  }
  
  PML <- -2 * sum(log(nrow(logf_comb)/colSums(1/exp(logf_comb))))
  loglik <- apply(logf, 1, sum) + out$logf.expert
  loglik.bar <- apply(logf.hat, 1, sum) + out$logf.hat.expert
  D.theta <- -2 * loglik
  D.bar <- -2 * loglik.bar
  pD <- mean(D.theta) - D.bar
  if (pD < 0) {
    warning(paste0("pD is ", round(pD), " for ", distr3, 
                   "; DIC estimates unreliable, use WAIC or PML."))
  }
  pV <- 0.5 * var(D.theta)
  dic <- mean(D.theta) + pD
  dic2 <- mean(D.theta) + pV
  aic <- D.bar + 2 * npars
  bic <- D.bar + npars * log(data.stan$n)
  list(aic = aic, bic = bic, dic = dic, dic2 = dic2, waic = WAIC, 
       pml = PML)
}
           
# hmc = c(Exponential = "exp",Gamma = "gam", GenGamma = "gga", 
#         Gompertz = "gom", 
#         RP = "rps", WeibullAF = "wei", WeibullPH = "wph", 
#         logLogistic = "llo", logNormal = "lno")


format_output_fit.models <- function (output, method, distr, formula, data){
  labs <- survHE:::manipulate_distributions(distr)$labs
  models <- lapply(output, function(x) x$model)
  model.fitting <- list(aic = unlist(lapply(output, function(x) x$aic)), 
                        bic = unlist(lapply(output, function(x) x$bic)), dic = unlist(lapply(output, 
                                                                                             function(x) x$dic)))
  misc <- list(time2run = unlist(lapply(output, function(x) x$time2run)), 
               formula = formula, data = data, model_name = unlist(lapply(output, 
                                                                          function(x) x$model_name)))
  if (any(distr == "polyweibull")) {
    misc$km = lapply(formula, function(f) survHE:::make_KMmake_KM(f, data))
  }
  else {
    misc$km = survHE:::make_KM(formula, data)
  }
  if (method == "hmc") {
    misc$data.stan <- lapply(output, function(x) x$data.stan)
    model.fitting$dic2 <- unlist(lapply(output, function(x) x$dic2))
    model.fitting$waic <- unlist(lapply(output, function(x) x$waic))
    model.fitting$pml <- unlist(lapply(output, function(x) x$pml))
    
  }
  names(models) <- labs
  res <- list(models = models, model.fitting = model.fitting, 
              method = method, misc = misc)
  class(res) <- "survHE"
  return(res)
}

make_sim_hmc <- function (m, t, X, nsim, newdata, dist, summary_stat, ...){
  iter_stan <- m@stan_args[[1]]$iter
  beta = rstan::extract(m)$beta
  if (ncol(X) == 1) {
    beta = beta[, 1,drop = F]
  }
  # if (dist == "rps" & any(grepl("Intercept", colnames(X)))) {
  #   X <- as.matrix(as_tibble(X) %>% select(-`(Intercept)`))
  #   beta = beta[, -ncol(beta)]
  # }
  linpred <- beta %*% t(X)
  sim <- lapply(1:nrow(X), function(x) {
    do.call(paste0("rescale_hmc_", dist), args = list(m, 
                                                      X, linpred[, x]))
  })
  if (nsim > iter_stan) {
    stop("Please select a value for 'nsim' that is less than or equal to the value set in the call to 'fit.models'")
  }
  if (nsim == 1) {
    sim <- lapply(sim, function(x) as.matrix(as_tibble(x) %>% 
                                               summarise_all(summary_stat), nrow = 1, ncol = ncol(x)))
  }
  if (nsim > 1 & nsim < iter_stan) {
    sim <- lapply(sim, function(x) as.matrix(as_tibble(x) %>% 
                                               sample_n(nsim, replace = FALSE), nrow = nsim, ncol = ncol(x)))
  }
  return(sim)
}






get_Surv <- function(dist, time, param1 = NULL, param2 = NULL, param3 = NULL, log = F, data.stan = NULL){
  
  if(dist == "wei"){
      return(pweibull(time, 
                    shape = param1, scale = param2, lower.tail = F, log = log))
  }
  
  if(dist == "wph"){
      return(pweibullPH(time, 
                      shape = param1, scale = param2, lower.tail = F, log = log))
  }
  
  if(dist == "exp"){
    return(pexp(time, rate = param1, lower.tail = F, log = log))
    
  }
  
  if(dist == "gam"){
    return(pgamma(time, shape = param1, rate = param2, lower.tail =  F, log = log))
  }
  
  if(dist == "gga"){
    return(pgengamma(time, mu = param1, sigma = param2, Q = param3, lower.tail =  F, log = log))
  }
  
  if(dist == "gom"){
    return(pgompertz(time, shape = param1, rate = param2, lower.tail =  F, log = log))
  }
  
  if(dist == "lno"){
    return(plnorm(time, 
                  meanlog  = param1, sdlog  = param2, lower.tail = F, log = log))
  }
  
  if(dist == "llo"){
    return(pllogis(time, 
                  shape  = param1, scale  = param2, lower.tail = F, log = log))
  }
  
  if(dist == "rps"){
    #2 ways to do it -- need to check if it is valid
    #eta <-  param1*data.stan$B_expert[which(time == data.stan$time_expert),] + param2
    #return(exp(-exp(eta)))
     return(psurvspline(time, gamma = param1, knots= data.stan$knots,lower.tail = F, log = log, offset = param2 ))
    
  }
  
  
}


get_mean_diff <- function(dist, time, param1 = NULL, param2 = NULL, param3 = NULL, log = F, data.stan = NULL){
  
  if(dist == "wei"){
    return(mean_weibull(shape = param1, scale = param2[1])-mean_weibull(shape = param1, scale = param2[2]))
  }
  
  if(dist == "wph"){
    return(mean_weibullPH(shape = param1, scale = param2[1])-mean_weibull(shape = param1, scale = param2[2]))
  }
  
  if(dist == "exp"){
    return(mean_exp(rate = param1[1])-mean_exp(rate = param1[2]))
    
  }
  
  if(dist == "gam"){
    return(mean_gamma(shape = param1, rate = param2[1])- mean_gamma(shape = param1, rate = param2[2]))
  }
  
  if(dist == "gga"){
    return(mean_gengamma(mu = param1[1], sigma = param2, Q = param3)-mean_gengamma(mu = param1[2], sigma = param2, Q = param3))
  }
  
  if(dist == "gom"){
    return(mean_gompertz(shape = param1, rate = param2[1])-mean_gompertz(shape = param1, rate = param2[2]))
  }
  
  if(dist == "lno"){
    return(mean_lnorm(meanlog  = param1[1], sdlog  = param2)-mean_lnorm(meanlog  = param1[2], sdlog  = param2))
  }
  
  if(dist == "llo"){
    return(mean_llogis(shape  = param1[1], scale  = param2)-mean_llogis(shape  = param1[2], scale  = param2))
  }
  
  # if(dist == "rps"){
  # 
  # }
  
  
}



expert_like <- function(data.stan, dist_surv, param1, param2 =NULL, param3= NULL){
  
  log_lik <- rep(NA, length(data.stan$time_expert))
  
  for(i in 1:length(data.stan$time_expert)){
     
      if(data.stan$St_indic ==1){ #Survival
        
        output <- get_Surv(dist_surv, data.stan$time_expert[i], 
                           param1 =param1, param2 = param2, param3 = param3, data.stan = data.stan)
      }else{# Add code for mean
        
        output <- get_mean_diff(dist_surv,param1 =param1, param2 = param2, param3 = param3, data.stan = data.stan)
        
        
      }
      
    if(data.stan$pool_type == 0){
      log_lik[i]  <-   expert_log_dens(x = output, df = data.stan$param_expert[,,i], pool_type = data.stan$pool_type, k_norm = data.stan$k_norm[i])
    }else{
      log_lik[i]  <- expert_log_dens(x = output, df = data.stan$param_expert[,,i], pool_type = data.stan$pool_type)
      
    }
    
  }
  
  return(sum(log_lik))
  
}


lik_rps <- function (x, linpred, linpred.hat, model, data.stan){
  
  dist <- "rps"
  gamma <- rstan::extract(model)$gamma
  gamma.hat <- apply(gamma, 2, median)
  linpred.hat <- as.numeric(linpred.hat)
  
  # LL<- apply(gamma_iters, 1, function(x){data.stan$d*log(hsurvspline(data.stan$t, gamma = x, knots = m.all$misc$data.stan[[1]]$knots))+
  #     psurvspline(q = data.stan$t, gamma = x, knots =  m.all$misc$data.stan[[1]]$knots, lower.tail = F, log.p =T)})
  # 
  # 
  # 
  # logf <- data.stan$d * (-log(data.stan$t) + log(gamma %*% 
  #                                                  t(data.stan$DB)) + gamma %*% t(data.stan$B) + linpred) - 
  #   exp(gamma %*% t(data.stan$B) + linpred)
  # logf.hat <- t(data.stan$d * (-log(data.stan$t) + log(data.stan$DB %*% 
  #                                                        gamma.hat) + data.stan$B %*% gamma.hat + linpred.hat) - 
  #                 exp(data.stan$B %*% gamma.hat + linpred.hat))
  
  
  logf.hat <- array(dim =  c(1,dim(linpred)[2]))
  
  
  if(all(data.stan$X ==0)){
    
    logf<- apply(gamma, 1, function(x){data.stan$d*log(hsurvspline(data.stan$t, gamma = x, knots = data.stan$knots))+
        psurvspline(q = data.stan$t, gamma = x, knots =  data.stan$knots, lower.tail = F, log.p =T)})
    logf <- t(logf)
  }else{
    logf <- array(dim = dim(linpred))
    #probably can be optimized
    for(i in 1:nrow(logf)){
      
      for(j in 1:ncol(logf)){
        logf[i,j] <- data.stan$d[j]*log(hsurvspline(data.stan$t[j], gamma = gamma[i,], knots = data.stan$knots, offset = linpred[i,j]))+
          psurvspline(q = data.stan$t[j], gamma = gamma[i,], knots =  data.stan$knots, lower.tail = F, log.p =T, offset = linpred[i,j])
        
      }
      
    }
  }
  

  for(i in 1:ncol(logf.hat)){
    logf.hat[i] <- data.stan$d[i]*log(hsurvspline(data.stan$t[i], gamma = gamma.hat, knots = data.stan$knots, offset = linpred.hat[i]))+
      psurvspline(q = data.stan$t[i], gamma = gamma.hat, knots =  data.stan$knots, lower.tail = F, log.p =T, offset = linpred.hat[i])
    
  }
  
  logf.expert <- rep(NA, nrow(linpred))
  
  if(data.stan$St_indic == 1){
    index_vec  <- data.stan$id_St
    for(i in 1:nrow(linpred)){
      logf.expert[i] <-  expert_like(data.stan, dist_surv = dist,param1 = gamma[i,], param2 = linpred[index_vec])
    }
    logf.hat.expert <- expert_like(data.stan, dist_surv = dist, param1 = gamma.hat, param2 = linpred.hat[index_vec])
  }else{
    index_vec <- c(data.stan$id_trt,data.stan$id_comp)
    #Enter code for Difference in survival
  }
  
  
  npars <- length(gamma.hat) + sum(apply(data.stan$X, 2, function(x) 1 - 
                                           all(x == 0)))
  list(logf = logf, logf.hat = logf.hat, npars = npars, f = NULL, 
       f.bar = NULL, s = NULL, s.bar = NULL, logf.expert = logf.expert, logf.hat.expert = logf.hat.expert)
}

lik_exp <- function (x, linpred, linpred.hat, model, data.stan){
  dist = "exp"
  
  logf <- matrix(unlist(lapply(1:nrow(linpred), function(i) {
    data.stan$d * log(hexp(data.stan$t, exp(linpred[i, ]))) + 
      log(1 - pexp(data.stan$t, exp(linpred[i, ])))
  })), nrow = nrow(linpred), byrow = T)
  logf.hat <- matrix(data.stan$d * log(hexp(data.stan$t, exp(linpred.hat))) + 
                       log(1 - pexp(data.stan$t, exp(linpred.hat))), nrow = 1)
  
  logf.expert <- rep(NA, nrow(linpred))
  
  
  if(data.stan$St_indic == 1){
    index_vec  <- data.stan$id_St
  }else{
    index_vec <- c(data.stan$id_trt,data.stan$id_comp)
    
  }
  
  for(i in 1:nrow(linpred)){
      logf.expert[i] <-  expert_like(data.stan, dist_surv = dist,param1 = exp(linpred[i,index_vec]))
    }
  logf.hat.expert <- expert_like(data.stan, dist_surv = dist, param1 = exp(linpred.hat[1,index_vec]))

  npars <- 1 + sum(1 - apply(data.stan$X, 2, function(x) all(x == 0)))
  list(logf = logf, logf.hat = logf.hat, npars = npars, f = NULL, 
       f.bar = NULL, s = NULL, s.bar = NULL, logf.expert = logf.expert, logf.hat.expert = logf.hat.expert)
}

lik_wei <- function (x, linpred, linpred.hat, model, data.stan ){
  
  dist = "wei"
  shape <- alpha <- as.numeric(rstan::extract(model)$alpha)
  shape.hat <- median(shape)
  logf <- matrix(unlist(lapply(1:nrow(linpred), function(i) {
    data.stan$d * log(hweibull(data.stan$t, shape[i], exp(linpred[i, 
    ]))) + log(1 - pweibull(data.stan$t, shape[i], exp(linpred[i, 
    ])))
  })), nrow = nrow(linpred), byrow = T)
  
  
  logf.hat <- matrix(data.stan$d * log(hweibull(data.stan$t, 
                                                shape.hat, exp(linpred.hat))) + log(1 - pweibull(data.stan$t, 
                                                                                                 shape.hat, exp(linpred.hat))), nrow = 1)
  logf.expert <- rep(NA, nrow(linpred))
  
  if(data.stan$St_indic == 1){
    index_vec  <- data.stan$id_St
  }else{
    index_vec <- c(data.stan$id_trt,data.stan$id_comp)
    
  }
    for(i in 1:nrow(linpred)){
      logf.expert[i] <-  expert_like(data.stan, dist_surv = dist, param1 = shape[i], param2 = exp(linpred[i,index_vec]))
    }

  logf.hat.expert <- expert_like(data.stan, dist_surv = dist, param1 = shape.hat, param2 = exp(linpred.hat[1,index_vec]))

  npars <- 2 + sum(1 - apply(data.stan$X, 2, function(x) all(x == 0)))
  
  
  list(logf = logf, logf.hat = logf.hat, npars = npars, f = NULL, 
       f.bar = NULL, s = NULL, s.bar = NULL, logf.expert = logf.expert, logf.hat.expert = logf.hat.expert)
}


lik_lno <- function (x, linpred, linpred.hat, model, data.stan){
  dist = "lno"
  
  sigma = as.numeric(rstan::extract(model)$alpha)
  sigma.hat = median(sigma)
  logf <- matrix(unlist(lapply(1:nrow(linpred), function(i) {
    data.stan$d * log(hlnorm(data.stan$t, (linpred[i, ]), 
                             sigma[i])) + log(1 - plnorm(data.stan$t, (linpred[i, 
                             ]), sigma[i]))
  })), nrow = nrow(linpred), byrow = T)
  logf.hat <- matrix(data.stan$d * log(hlnorm(data.stan$t, 
                                              (linpred.hat), sigma.hat)) + log(1 - plnorm(data.stan$t, 
                                                                                          (linpred.hat), sigma.hat)), nrow = 1)
  logf.expert <- rep(NA, nrow(linpred))

  if(data.stan$St_indic == 1){
    index_vec  <- data.stan$id_St
  }else{
    index_vec <- c(data.stan$id_trt,data.stan$id_comp)
    
  }
    
  for(i in 1:nrow(linpred)){
      logf.expert[i] <-  expert_like(data.stan, dist_surv = dist, param1 = sigma[i], param2 = exp(linpred[i,index_vec]))
  }
  logf.hat.expert <- expert_like(data.stan, dist_surv = dist, param1 = sigma.hat, param2 = exp(linpred.hat[1,index_vec]))

  npars <- 2 + sum(1 - apply(data.stan$X, 2, function(x) all(x ==0)))
  list(logf = logf, logf.hat = logf.hat, npars = npars, f = NULL, 
       f.bar = NULL, s = NULL, s.bar = NULL, logf.expert = logf.expert, logf.hat.expert = logf.hat.expert)
}

lik_llo <- function (x, linpred, linpred.hat, model, data.stan){
  dist = "llo"
  sigma = as.numeric(rstan::extract(model)$alpha)
  sigma.hat = median(sigma)
  logf <- matrix(unlist(lapply(1:nrow(linpred), function(i) {
    data.stan$d * log(hllogis(data.stan$t, sigma[i], exp(linpred[i, 
    ]))) + log(1 - pllogis(data.stan$t, sigma[i], exp(linpred[i, 
    ])))
  })), nrow = nrow(linpred), byrow = T)
  logf.hat <- matrix(data.stan$d * log(hllogis(data.stan$t, 
                                               sigma.hat, exp(linpred.hat))) + log(1 - pllogis(data.stan$t, 
                                                                                               sigma.hat, exp(linpred.hat))), nrow = 1)
  
  logf.expert <- rep(NA, nrow(linpred))
  
  if(data.stan$St_indic == 1){
    index_vec  <- data.stan$id_St
  }else{
    index_vec <- c(data.stan$id_trt,data.stan$id_comp)
    
  }
   
   for(i in 1:nrow(linpred)){
      logf.expert[i] <-  expert_like(data.stan, dist_surv = dist, param1 = sigma[i], param2 = exp(linpred[i,index_vec]))
    }
    logf.hat.expert <- expert_like(data.stan, dist_surv = dist, param1 = sigma.hat, param2 = exp(linpred.hat[1,index_vec]))
  
  
  npars <- 2 + sum(1 - apply(data.stan$X, 2, function(x) all(x == 
                                                               0)))
  list(logf = logf, logf.hat = logf.hat, npars = npars, f = NULL, 
       f.bar = NULL, s = NULL, s.bar = NULL, logf.expert = logf.expert, logf.hat.expert = logf.hat.expert)
}




lik_wph <- function (x, linpred, linpred.hat, model, data.stan){
  dist = "wph"
  shape <- alpha <- as.numeric(rstan::extract(model)$alpha)
  shape.hat = median(shape)
  logf <- matrix(unlist(lapply(1:nrow(linpred), function(i) {
    data.stan$d * log(hweibullPH(data.stan$t, shape[i], exp(linpred[i, 
    ]))) + log(1 - pweibullPH(data.stan$t, shape[i], 
                              exp(linpred[i, ])))
  })), nrow = nrow(linpred), byrow = T)
  logf.hat <- matrix(data.stan$d * log(hweibullPH(data.stan$t, 
                                                  shape.hat, exp(linpred.hat))) + log(1 - pweibullPH(data.stan$t, 
                                                                                                     shape.hat, exp(linpred.hat))), nrow = 1)
  
  logf.expert <- rep(NA, nrow(linpred))
  
  if(data.stan$St_indic == 1){
    index_vec  <- data.stan$id_St
  }else{
    index_vec <- c(data.stan$id_trt,data.stan$id_comp)
    
  }
  
  for(i in 1:nrow(linpred)){
      logf.expert[i] <-  expert_like(data.stan, dist_surv = dist, param1 = shape[i], param2 = exp(linpred[i,index_vec]))
    }
  logf.hat.expert <- expert_like(data.stan, dist_surv = dist, param1 = shape.hat, param2 = exp(linpred.hat[1,index_vec]))

  
  npars <- 2 + sum(1 - apply(data.stan$X, 2, function(x) all(x == 
                                                               0)))
  list(logf = logf, logf.hat = logf.hat, npars = npars, f = NULL, 
       f.bar = NULL, s = NULL, s.bar = NULL, logf.expert = logf.expert, logf.hat.expert = logf.hat.expert)
}



lik_gam <- function (x, linpred, linpred.hat, model, data.stan){
  dist = "gam"
  
  shape <- alpha <-  as.numeric(model$BUGSoutput$sims.matrix[ , grep("alpha",colnames(model$BUGSoutput$sims.matrix))])
  shape.hat <- median(shape)
  
  logf <- matrix(unlist(lapply(1:nrow(linpred), function(i) {
    data.stan$d * log(hgamma(data.stan$t, shape = shape[i], 
                             rate = exp(linpred[i, ]))) + 
      pgamma(q = data.stan$t,shape[i], rate = exp(linpred[i, ]), lower.tail = F, log = T)
  })), nrow = nrow(linpred), byrow = T)

  logf.hat <- matrix(data.stan$d * log(hgamma(data.stan$t, 
                                              shape.hat, exp(linpred.hat))) + 
                       pgamma(data.stan$t,shape.hat, exp(linpred.hat),lower.tail = F,
                              log = T), nrow = 1)  
  
  logf.expert <- rep(NA, nrow(linpred))

  if(data.stan$St_indic == 1){
    index_vec  <- data.stan$id_St
  }else{
    index_vec <- c(data.stan$id_trt,data.stan$id_comp)
    
  }
    
    for(i in 1:nrow(linpred)){
      logf.expert[i] <-  expert_like(data.stan, dist_surv = dist, param1 = shape[i], param2 = exp(linpred[i,index_vec]))
    }
    logf.hat.expert <- expert_like(data.stan, dist_surv = dist, param1 = shape.hat, param2 = exp(linpred.hat[1,index_vec]))

  npars <- 2 + sum(1 - apply(data.stan$X, 2, function(x) all(x == 0)))
  
  list(logf = logf, logf.hat = logf.hat, npars = npars, f = NULL, 
       f.bar = NULL, s = NULL, s.bar = NULL, logf.expert = logf.expert, logf.hat.expert = logf.hat.expert)
}


lik_gom <- function (x, linpred, linpred.hat, model, data.stan){
  dist = "gom"
  shape <- alpha <- as.numeric(model$BUGSoutput$sims.matrix[ , grep("alpha",colnames(model$BUGSoutput$sims.matrix))])
  shape.hat = median(shape)
  logf <- matrix(unlist(lapply(1:nrow(linpred), function(i) {
    data.stan$d * log(hgompertz(data.stan$t, shape = shape[i], 
                                rate = exp(linpred[i, ]))) + 
      pgompertz(data.stan$t,shape[i], rate = exp(linpred[i, ]), lower.tail = F, log = T)
  })), nrow = nrow(linpred), byrow = T)
  logf.hat <- matrix(data.stan$d * log(hgompertz(data.stan$t,shape.hat, exp(linpred.hat))) + 
                       pgompertz(data.stan$t,shape.hat, exp(linpred.hat), lower.tail = F, log = T), nrow = 1)
  
  
  if(data.stan$St_indic == 1){
    index_vec  <- data.stan$id_St
  }else{
    index_vec <- c(data.stan$id_trt,data.stan$id_comp)
    
  }
  logf.expert <- rep(NA, nrow(linpred))
  
  for(i in 1:nrow(linpred)){
      logf.expert[i] <-  expert_like(data.stan, dist_surv = dist, param1 = shape[i], param2 = exp(linpred[i,index_vec]))
    }
    logf.hat.expert <- expert_like(data.stan, dist_surv = dist, param1 = shape.hat, param2 = exp(linpred.hat[1,index_vec]))
  
  npars <- 2 + sum(1 - apply(data.stan$X, 2, function(x) all(x == 0)))
  list(logf = logf, logf.hat = logf.hat, npars = npars, f = NULL, 
       f.bar = NULL, s = NULL, s.bar = NULL, logf.expert = logf.expert, logf.hat.expert = logf.hat.expert)
}



lik_gga <- function (x, linpred, linpred.hat, model, data.stan){
  dist = "gga"
  q = as.numeric(model$BUGSoutput$sims.matrix[ , grep("Q",colnames(model$BUGSoutput$sims.matrix))])
  q.bar = median(q)
  scale = as.numeric(model$BUGSoutput$sims.matrix[ , grep("sigma",colnames(model$BUGSoutput$sims.matrix))])
  scale.bar = median(scale)
  
  d2 <- sapply(data.stan$d,function(x){ifelse(x == 1, 0,1)})
 
  logf <- matrix(unlist(lapply(1:nrow(linpred), function(i) {
      data.stan$d*dgengamma(data.stan$t, 
                                linpred[i, ], scale[i], q[i], log = T) +
      d2*pgengamma(data.stan$t,linpred[i, ], scale[i], q[i], log = T, lower.tail = F)})),
                              nrow = nrow(linpred), byrow = T)
  

  logf.hat <- matrix(data.stan$d*dgengamma(data.stan$t,linpred.hat,scale.bar, q.bar, log = T) +
                     d2*pgengamma(data.stan$t, linpred.hat, scale.bar, q.bar, log = T, lower.tail = F),
                     nrow = 1)
  
  
  if(data.stan$St_indic == 1){
    index_vec  <- data.stan$id_St
    }else{
    index_vec <- c(data.stan$id_trt,data.stan$id_comp)
   }
  
  logf.expert <- rep(NA, nrow(linpred))
  

    for(i in 1:nrow(linpred)){
      logf.expert[i] <-  expert_like(data.stan, dist_surv = dist, param1 = linpred[i,index_vec], param2 =scale[i], param3 = q[i])
    }
    logf.hat.expert <- expert_like(data.stan, dist_surv = dist, param1 = linpred.hat[1,index_vec], param2 = scale.bar, 
                                   param3 = q.bar[index_vec])

  
  npars <- 3 + sum(1 - apply(data.stan$X, 2, function(x) all(x == 0)))
  list(logf = logf, logf.hat = logf.hat, npars = npars, f = NULL, 
       f.bar = NULL, s = NULL, s.bar = NULL, logf.expert = logf.expert, logf.hat.expert = logf.hat.expert)
  
}


get_stats_hmc <- function (x, mod){
  
  if(names(x[["models"]])[mod] %in% c("Gamma", "Gen. Gamma","Gompertz")){
    table = x$models[[mod]]$BUGSoutput$summary[, c("mean", 
                                           "sd", "2.5%", "97.5%")]
  }else{
    table = rstan::summary(x$models[[mod]])$summary[, c("mean", 
                                                        "sd", "2.5%", "97.5%")]
    table = table[-grep("lp__", rownames(table)), ]
  }

 

    if (any(apply(x$misc$data.stan[[1]]$X, 2, function(x) all(x == 
                                                              0)))) {
      table = table[-grep("beta\\[2\\]", rownames(table)), ]
    }

  # if(names(x[["models"]])[mod] %in% c("Gen. Gamma")){
  #   
  #   
  # }
    res = do.call(eval(parse(text=paste0("survHE:::rescale_stats_hmc_", x$misc$model_name[mod]))), 
                  args = list(table = table, x = x))


  return(res)
}


load_availables <- function(){
  availables = list(mle = c(genf = "gef", genf.orig = "gof", 
                            gengamma = "gga", gengamma.orig = "ggo", 
                            exp = "exp", weibull = "wei", weibullPH = "wph", 
                            lnorm = "lno", gamma = "gam", gompertz = "gom", 
                            llogis = "llo", lognormal = "lno", rps = "rps"), 
                    inla = c(exponential = "exp", weibull = "wei", 
                             weibullPH = "wph", lognormal = "lno", 
                             loglogistic = "llo", rps = "rps"), 
                    hmc = c(Exponential = "exp",Gamma = "gam", GenGamma = "gga", 
                            Gompertz = "gom", 
                            RP = "rps", WeibullAF = "wei", WeibullPH = "wph", 
                            logLogistic = "llo", logNormal = "lno"))
  return(availables)
}

print.survHE <-function (x, mod = 1, ...) 
{
  adjust_survHE_func()
  
  exArgs <- list(...)
  availables <- load_availables()
  if (!exists("digits", where = exArgs)) {
    digits = 6
  }
  else {
    digits = exArgs$digits
  }
  if (!exists("original", where = exArgs)) {
    original = FALSE
  }
  else {
    original = exArgs$original
  }
  if (exists("orig", exArgs)) {
    original = exArgs$orig
  }
  if (original == TRUE) {
    do.call(paste0("original_table_", x$method), args = list(x, 
                                                             mod, digits))
  }
  else {
    method_eval <- paste0("get_stats_", x$method)
    res = do.call(method_eval, args = list(x,mod))
    survHE:::format_table(x, mod, res, digits)
  }
}


#### Adjusts SurvHE functions:

#error in a rps function 
adjust_survHE_func <- function(){


make_sim_hmc <- function (m, t, X, nsim, newdata, dist, summary_stat, ...){
  
  if(class(m)== "rjags"){
    iter_stan <- m[["n.iter"]]
    beta <- m$BUGSoutput$sims.matrix[, grep("beta",colnames(m$BUGSoutput$sims.matrix))]
  }else{
    
    iter_stan <- m@stan_args[[1]]$iter
    beta = rstan::extract(m)$beta
  }
  if (ncol(X) == 1) {
    beta = beta[, 1, drop = F] #PC: add drop to stop it coreceing to a vector
  }
  if (dist == "rps" & any(grepl("Intercept", colnames(X)))) {
    X <- as.matrix(as_tibble(X) %>% select(-`(Intercept)`))
    beta = beta[, -ncol(beta)]
  }
  linpred <- beta %*% t(X)
  sim <- lapply(1:nrow(X), function(x) {
    do.call(paste0("rescale_hmc_", dist), args = list(m, 
                                                      X, linpred[, x]))
  })
  if (nsim > iter_stan) {
    stop("Please select a value for 'nsim' that is less than or equal to the value set in the call to 'fit.models'")
  }
  if (nsim == 1) {
    sim <- lapply(sim, function(x) as.matrix(as_tibble(x) %>% 
                                               summarise_all(summary_stat), nrow = 1, ncol = ncol(x)))
  }
  if (nsim > 1 & nsim < iter_stan) {
    sim <- lapply(sim, function(x) as.matrix(as_tibble(x) %>% 
                                               sample_n(nsim, replace = FALSE), nrow = nsim, ncol = ncol(x)))
  }
  return(sim)
}

tmpfun <- get("make_sim_hmc", envir = asNamespace("survHE"))
environment(make_sim_hmc) <- environment(tmpfun)
attributes(make_sim_hmc) <- attributes(tmpfun)  
assignInNamespace("make_sim_hmc", make_sim_hmc, ns="survHE")




rescale_hmc_gam <- function (m, X, linpred){
  if(class(m)== "rjags"){
    shape <- as.numeric(m$BUGSoutput$sims.matrix[,"alpha"])
  }else{
    shape <- as.numeric(rstan::extract(m)$alpha)
  }
 
  rate <- exp(linpred)
  sim <- cbind(shape, rate)
  colnames(sim) <- c("shape", "rate")
  return(sim)
}

tmpfun <- get("rescale_hmc_gam", envir = asNamespace("survHE"))
environment(rescale_hmc_gam) <- environment(tmpfun)
attributes(rescale_hmc_gam) <- attributes(tmpfun)  
assignInNamespace("rescale_hmc_gam", rescale_hmc_gam, ns="survHE")



rescale_hmc_gom <- function (m, X, linpred){
  if(class(m)== "rjags"){
    shape <- as.numeric(m$BUGSoutput$sims.matrix[,"alpha"])
  }else{
    shape <- as.numeric(rstan::extract(m)$alpha)
  }
  rate <- exp(linpred)
  sim <- cbind(shape, rate)
  colnames(sim) <- c("shape", "rate")
  return(sim)
}


tmpfun <- get("rescale_hmc_gom", envir = asNamespace("survHE"))
environment(rescale_hmc_gom) <- environment(tmpfun)
attributes(rescale_hmc_gom) <- attributes(tmpfun)  
assignInNamespace("rescale_hmc_gom", rescale_hmc_gom, ns="survHE")

rescale_hmc_gga<- function (m, X, linpred){
  
  if(class(m)== "rjags"){
    Q <- as.numeric(m$BUGSoutput$sims.matrix[,"Q"])
    sigma <- as.numeric(m$BUGSoutput$sims.matrix[,"sigma"])
  }else{
    Q <- as.numeric(rstan::extract(m)$Q)
    sigma <- as.numeric(rstan::extract(m)$sigma)
  }
  mu <- linpred
  sim <- cbind(mu, sigma, Q)
  colnames(sim) <- c("mu", "sigma", "Q")
  return(sim)
}
tmpfun <- get("rescale_hmc_gga", envir = asNamespace("survHE"))
environment(rescale_hmc_gga) <- environment(tmpfun)
attributes(rescale_hmc_gga) <- attributes(tmpfun)  
assignInNamespace("rescale_hmc_gga", rescale_hmc_gga, ns="survHE")

get_stats_hmc <- function(x, mod){
  
  if(class(x$models[[mod]])== "rjags"){
   table =  x$models[[mod]]$BUGSoutput$summary[,c("mean", 
                                                  "sd", "2.5%", "97.5%")]
  }else{
    table = rstan::summary(x$models[[mod]])$summary[, c("mean", 
                                                        "sd", "2.5%", "97.5%")]
    table = table[-grep("lp__", rownames(table)), ]
    
  }

  if ("X_obs" %in% names(x$misc$data.stan[[1]])) {
    if (any(apply(x$misc$data.stan[[1]]$X_obs, 2, function(x) all(x == 
                                                                  0)))) {
      table = table[-grep("beta\\[2\\]", rownames(table)), 
      ]
    }
  }
  else {
    if (any(apply(x$misc$data.stan[[1]]$X, 2, function(x) all(x == 
                                                              0)))) {
      table = table[-grep("beta\\[2\\]", rownames(table)), 
      ]
    }
  }
  res = do.call(paste0("rescale_stats_hmc_", x$misc$model_name[mod]), 
                args = list(table = table, x = x))
  return(res)
}


tmpfun <- get("get_stats_hmc", envir = asNamespace("survHE"))
environment(get_stats_hmc) <- environment(tmpfun)
attributes(get_stats_hmc) <- attributes(tmpfun)  
assignInNamespace("get_stats_hmc", get_stats_hmc, ns="survHE")


rescale_stats_hmc_gam <- function (table, x){
  rate <- matrix(table[grep("rate", rownames(table)),], ncol = 4)
  rownames(rate) <- "rate"
  shape <- matrix(table[grep("alpha", rownames(table)), 
  ], ncol = 4)
  rownames(shape) <- "shape"
  effects = add_effects_hmc(table, x) #typo in this function
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


get_stats_hmc <- function (x, mod){
  if (class(x$models[[mod]]) == "rjags") {
    table = x$models[[mod]]$BUGSoutput$summary[, c("mean", 
                                                   "sd", "2.5%", "97.5%")]
  }
  else {
    table = rstan::summary(x$models[[mod]])$summary[, c("mean", 
                                                        "sd", "2.5%", "97.5%")]
    table = table[-grep("lp__", rownames(table)), ]
  }
  if ("X_obs" %in% names(x$misc$data.stan[[1]])) {
    if (any(apply(x$misc$data.stan[[1]]$X_obs, 2, function(x) all(x == 0)))) {
      #Error (mod instead of 1)
      #for RPS X matrix can be 0 for both columns
      beta_drop <- which(apply(x$misc$data.stan[[mod]]$X, 2, function(x) all(x == 0)) == TRUE)
      beta_drop <- paste0("beta\\[", beta_drop,"\\]")
      if(length(beta_drop)>1){
        beta_drop <- paste(beta_drop, collapse = "|")
      }
      table <-  table[-grep(beta_drop, rownames(table)), ]
    }
  }
  else {
    if (any(apply(x$misc$data.stan[[1]]$X, 2, function(x) all(x == 0)))) {
      beta_drop <- which(apply(x$misc$data.stan[[mod]]$X, 2, function(x) all(x == 0)) == TRUE)
      beta_drop <- paste0("beta\\[", beta_drop,"\\]")
      if(length(beta_drop)>1){
        beta_drop <- paste(beta_drop, collapse = "|")
      }
       table <-table[-grep(beta_drop, rownames(table)), ]
    }
  }
  res = do.call(paste0("rescale_stats_hmc_", x$misc$model_name[mod]), 
                args = list(table = table, x = x))
  return(res)
}

tmpfun <- get("get_stats_hmc", envir = asNamespace("survHE"))
environment(get_stats_hmc) <- environment(tmpfun)
attributes(get_stats_hmc) <- attributes(tmpfun)  
assignInNamespace("get_stats_hmc", get_stats_hmc, ns="survHE")

}
