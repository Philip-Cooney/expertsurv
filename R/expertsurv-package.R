


#' Incorporating Expert Opinion with Parametric Survival Models
#' 
#' Contains functions to include expert opinion with the parametric models commonly 
#' used in health economic modelling. Theoretical details are described elsewhere \insertCite{Cooney.2023}{expertsurv}.  Borrows many function from the ``survHE`` package \insertCite{Baio.2020}{expertsurv}. 
#' 
#' \tabular{ll}{ Package: \tab expertsurv \cr Type: \tab Package\cr Version: \tab
#' 1.0.0\cr Date: \tab 2023-01-27\cr License: \tab MIT + file LICENSE \cr LazyLoad: \tab
#' yes\cr } Integrate expert opinions on survival and mean differences in survival with common parametric survival models using either a Bayesian or frequentist framework.
#' 
#' @name expertsurv-package
#' 
#' @aliases expertsurv-package expertsurv
#' @docType package
#' @author 
#' Philip Cooney Package Creator, Maintainer
#' @author 
#' Arthur White  Thesis Supervisor
#' @template refs
#' @keywords Expert Opinion Survival Modelling Health Economic Evaluation
#' @examples
#' #Define expert opinion
#' require("dplyr")
#' param_expert_example1 <- list()
#' #1 timepoint and 2 experts with equal weight,
#' #first a normal distribution, second a non-standard t-distribution with
#' #3 degrees of freedom
#' 
#' param_expert_example1[[1]] <- data.frame(dist = c("norm","t"),
#'                                wi = c(0.5,0.5), # Ensure Weights sum to 1
#'                                param1 = c(0.1,0.12),
#'                                param2 = c(0.05,0.05),
#'                                 param3 = c(NA,3))
#' 
#' 
#' timepoint_expert <- 14
#' 
#' data2 <- data %>% rename(status = censored) %>% mutate(time2 = ifelse(time > 10, 10, time),
#'                   status2 = ifelse(time> 10, 0, status))
#'
#' example1  <- fit.models.expert(formula=Surv(time2,status2)~1,data=data2,
#'                               distr=c("wph", "gomp"),
#'                               method="mle",
#'                               pool_type = "log pool", 
#'                               opinion_type = "survival",
#'                               times_expert = timepoint_expert, 
#'                               param_expert = param_expert_example1)
#' 
#'  #Visualize the goodness of fit
#'  model.fit.plot(example1, type = "aic")
#'  #Visualize the survival curve
#'  plot(example1, add.km = TRUE, t = 0:30)
#' 
#' 
#' @references 
#' \insertRef{Baio.2020}{expertsurv}
#' 
#' \insertRef{Cooney.2023}{expertsurv}
#' 
NULL



