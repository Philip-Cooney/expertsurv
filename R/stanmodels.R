
<!-- README.md is generated from README.Rmd. Please edit that file -->

# expert-surv

<!-- badges: start -->
<!-- badges: end -->

The goal of `expertsurv` is to incorporate expert opinion into an
analysis of time to event data. `expertsurv` uses many of the core
functions of the `survHE` package (Baio 2020). Technical details of the
implementation are detailed in (Anon) and will not be
repeated here.

The key function is `fit.models.expert` and operates almost identically
to the `fit.models` function of `survHE`.

## Installation

You can install the released version of expertsurv from
[GitHub](https://github.com/Anon19820/expertsurv) with:

``` r
devtools::install_github("Anon19820/expertsurv")
```

## Expert Opinion on Survival at timepoints

If we have elicited expert opinion of the survival probability at
certain timepoint(s) and assigned distributions to these beliefs, we
encode that information as follows:

``` r
#A param_expert object; which is a list of 
#length equal to the number of timepoints
param_expert_example1 <- list()

#If we have 1 timepoint and 2 experts
#dist is the names of the distributions
#wi is the weight assigned to each expert (usually 1)
#param1, param2, param3 are the parameters of the distribution
#e.g. for norm, param1 = mean, param2 = sd
#param3 is only used for the t-distribution and is the degress of freedom.
#We allow the following distributions:
#c("normal","t","gamma","lognormal","beta") 


param_expert_example1[[1]] <- data.frame(dist = c("norm","t"),
                                         wi = c(0.5,0.5), # Ensure Weights sum to 1
                                         param1 = c(0.1,0.12),
                                         param2 = c(0.005,0.005),
                                         param3 = c(NA,3))
param_expert_example1
#> [[1]]
#>   dist wi param1 param2 param3
#> 1 norm  1   0.10  0.005     NA
#> 2    t  1   0.12  0.005      3

#Naturally we will specify the timepoint for which these probabilities where elicited

timepoint_expert <- 14


#In case we wanted a second timepoint -- Just for illustration

# param_expert_example1[[2]] <- data.frame(dist = c("norm","norm"),
#                                          wi = c(1,1),
#                                          param1 = c(0.05,0.045),
#                                          param2 = c(0.005,0.005),
#                                          param3 = c(NA,NA))
# 
# timepoint_expert <- c(timepoint_expert,18)
```

If we wanted opinions at multiple timepoints we just include append
another list (i.e. param\_expert\_example1\[\[2\]\] with the relevant
parameters) and specify timepoint\_expert as a vector of length 2 with
the second element being the second timepoint.

For details on assigning distributions to elicited probabilities and
quantiles see the `SHELF` package (Oakley 2021) and for an overview on
methodological approaches to eliciting expert opinion see (O’Hagan
2019). We can see both the individual and pooled distributions using the
following code (note that we could have used the output of the `fitdist`
function from `SHELF` if we actually elicited quantiles from an expert):

    plot_opinion1<- plot_expert_opinion(param_expert_example1[[1]], 
                        weights = param_expert_example1[[1]]$wi)
    ggsave("Vignette_Example 1 - Expert Opinion.png")

For the log pool we have a uni-modal distribution  (in contrast to the bi-modal linear pool) which has a 95%
credible interval between 9.0 − 11.9% calculated with the function
below:

<div class="figure">

<img src="Vignette_Example 1 - Expert Opinion.png" alt="Expert beliefs about survival represented as distributions" width="2721" />
<p class="caption">
Expert beliefs about survival represented as distributions
</p>

</div>

    cred_int_val <- cred_int(plot_opinion1,val = "linear pool", interval = c(0.025, 0.975))

We load and fit the data as follows (in this example considering just
the Weibull and Gompertz models), with pool\_type = “log pool”
specifying that we want to use the  logarithmic pooling (rather than default
“linear pool”). We do this as we wish to compare the results to the penalized maximum likelihood estimates in the next section.


    data2 <- survHE::data %>% rename(status = censored) %>% mutate(time2 = ifelse(time > 10, 10, time),
                                                                  status2 = ifelse(time> 10, 0, status))

    #Set the opinion type to "survival"

    example1  <- fit.models.expert(formula=Surv(time2,status2)~1,data=data2,
                                            distr=c("wei", "gomp"),
                                            method="hmc",
                                            iter = 5000,
                                            pool_type = "log pool", 
                                            opinion_type = "survival",
                                            times_expert = timepoint_expert, 
                                            param_expert = param_expert_example1)

Both visual fit and model fit statistics highlight that the Weibull
model is a poor fit to both the expert opinion and data (black line
referring to the 95% confidence region for the experts belief about survival at the timepoint).

    survHE::model.fit.plot(example1, type = "dic")

#N.B. survHE::plot plots the survival function at the posterior mean parameter values
#     while it is more robust to use the entire posterior sample (make.surv), however, in this case both results are similar. 
     plot(example1, add.km = T, t = 0:30)+
      theme_light()+
      scale_x_continuous(expand = c(0, 0), limits = c(0,NA), breaks=seq(0, 30, 2)) + 
      scale_y_continuous(expand = c(0, 0), limits = c(0, NA), breaks=seq(0, 1, 0.05))+
      geom_segment(aes(x = 14, y = cred_int_val[1], xend = 14, yend = cred_int_val[2]))

<div class="figure">

<img src="Vignette_Example 1 - DIC.png" alt="Model Comparison" width="2721" />
<p class="caption">
Model Comparison
</p>

</div>

<div class="figure">

<img src="Vignette_Example 1.png" alt="Survival function with Expert Information" width="2721" />
<p class="caption">
Survival function with Expert Information
</p>

</div>

## Expert Opinion using Penalized Maximum Likelihood

We can also fit the model by Penalized Maximum Likelihood approaches. Full integration of this function with the package is under development and currently is only implemented for survival probabilities (not mean survival difference).

We require to source the functions in "Flexsurv functions v2.R" available in the Github folder. It should be noted that the results will be very similar to the Bayesian approach when the expert opinion is unimodal (as maximum liklelihood produces a point estimate), therefore we use the logarithmic pool which is unimodal. The functions return flexsurvreg objects so the regular plot, summary functions work without modification. We find that the AIC values also favour the Gompertz model by a large factor (not shown).

```
#source("~/Flexsurv functions v2.R")
expert_opinion<-  list()

#Create - the data to add to the flexsurv function
expert_opinion$param_expert <- make_data_expert(param_expert_example1, timepoint_expert)
expert_opinion$times <- timepoint_expert
expert_opinion$pool <- 0 #linear pool is 1; log is 0

#Only currently implemented for survival probabilities 

if(expert_opinion$pool == 0){
  expert_opinion$k_norm <-  get_k_norm(param_expert_example1)
}else{
  expert_opinion$k_norm <-  NULL
}

#Fit the models - returns a flexsurv object so all flexsurv functions work on it.
fit_expert_mle_weibull <- try({flexsurvreg(formula = Surv(time2, status2) ~ 1, 
                                           data = data2, dist="weibullPH",expert_opinion = expert_opinion)},
                              silent = TRUE)

fit_expert_mle_gompertz <- try({flexsurvreg(formula = Surv(time2, status2) ~ 1, 
                                            data = data2, dist="gompertz",expert_opinion = expert_opinion)},
                               silent = TRUE)

gompertz_summary <- summary(fit_expert_mle_gompertz,t = seq(0,30, by = 0.1))

## Plot the outputs
png("MLE-Weibull-Gomp.png", width = 600, height = 400)
plot(fit_expert_mle_weibull, t = seq(0,30, by = 0.1),xlim= c(0,30),ci = F, xlab = "Time", ylab = "Survival" )
lines(x = gompertz_summary[[1]][,"time"], y = gompertz_summary[[1]][,"est"], col = "blue")
segments(x0 = timepoint_expert, y0 =  cred_int_val[1], y1 =  cred_int_val[2], col = "orange")
legend("topright", legend=c("Weibull", "Gompertz"),
       col=c("red", "blue"), lty=1:2, cex=0.8,
       title="Survival analysis with MLE", text.font=4, bg='lightblue')
dev.off()
```
</div>

<div class="figure">

<img src="MLE-Weibull-Gomp.png" alt="Survival function with Expert Information-Penalized Maximum Likelihood" width="2721" />
<p class="caption">
Survival function with Expert Information-Penalized Maximum Likelihood
</p>

</div>

After using this function you and wish to use the regular flexsurv pacakge you should run the following commands:

    unloadNamespace("flexsurv") #Unload flexsurv and associated name spaces
    require("flexsurv") #reload flexsurv
    
    
    
## Expert Opinion on Survival of a comparator arm

In this situation we place an opinion on the comparator arm.

    param_expert_example2[[1]] <- data.frame(dist = c("norm"),
                                             wi = c(1),
                                             param1 = c(0.1),
                                             param2 = c(0.005),
                                             param3 = c(NA))

``` r
#Check the coding of the arm variable
#Comparator is 0, which is our id_St
unique(survHE::data$arm)
#> [1] 0 1
```

    survHE.data.model  <- fit.models.expert(formula=Surv(time2,status2)~as.factor(arm),data=data2,
                                            distr=c("wei"),
                                            method="hmc",
                                            iter = 5000,
                                            opinion_type = "survival",
                                            id_St = 0, 
                                            times_expert = timepoint_expert, 
                                            param_expert = param_expert_example2)

We can remove the impact of expert opinion by running the same model in
the `survHE` package. Alternatively we note that a ℬ(1, 1) distribution
is uniform on the survival probability and does not change the
likelihood.

    param_expert_vague <- list()
    param_expert_vague[[1]] <- data.frame(dist = "beta", wi = 1, param1 = 1, param2 = 1, param2 = NA)

<div class="figure">

<img src="Vignette_Example 2.png" alt="Survival function with Expert Information (left) and Vague prior - no information (right)" width="3000" />
<p class="caption">
Survival function with Expert Information (left) and Vague prior - no information (right)
</p>

</div>

The survival function for “arm 1” has been shifted downwards slightly,
however the covariate for the accelerated time factor has markedly
increased to counteract the lower survival probability for the reference
(arm 0).

## Expert Opinion on Survival Difference

This example illustrates an opinion on the survival difference. For
illustration we use the Gompterz, noting that a negative shape parameter
will lead to a proportion of subjects living forever. Clearly the mean
is not defined in these cases so the code automatically constrains the
shape to be positive.

    param_expert3 <- list()

    #Expert belief of 5 "months" difference in expected survival
    param_expert3[[1]] <- data.frame(dist = "norm", wi = 1, param1 = 5, param2 = 0.2, param3 = NA)


    survHE.data.model  <- fit.models.expert(formula=Surv(time2,status2)~as.factor(arm),data=data2,
                                                         distr=c("gom"),
                                                         method="hmc",
                                                         iter = 5000,
                                                         opinion_type = "mean",
                                                         id_trt = 1, # Survival difference is Mean_surv[id_trt]- Mean_surv[id_comp] 
                                                         param_expert = param_expert3)

<div class="figure">

<img src="Vignette_Example 3.png" alt="Survival difference" width="3000" />
<p class="caption">
Survival difference
</p>

</div>

## Compatability with survHE

As stated in the introduction this package relies on many of the core
functions of the `survHE` package (Baio 2020) (i.e. note the use of
`survHE::model.fit.plot` in the first example, meaning that the function
is sourced directly from `survHE`). In theory a new version of `survHE`
could result in a lack of compatibility with this package, however, any
required changes should be minor. Because the objective of this package
was to fit the models with expert opinion, plot the survival curves and
compare the goodness of fit, these capabilities (which have been
presented in this README) have been tested for compatibility. Other
functions should (in theory) be compatible (again by adding `survHE::`
to the relevant function), however, I have not tested all these
potential use cases. If you run in issues, bugs or just features which
you feel would be useful, please let me know (<anon.email>) and I
will investigate and update as required.

Additionally I have made modifications to some of the `survHE` functions
to accommodate JAGS models (by changing the namespace of the `survHE`
environment). These should have no impact on the operation of `survHE`
and these changes are only invoked when `expertsurv` is loaded. However,
in the situation where you would like to revert to `survHE` functions
during the session, simply run the following:

    unloadNamespace("survHE") #Unload survHE and associated name spaces
    require("survHE") #reload survHE

One practical difference between the packages is the calculation of DIC
(Deviance Information Criterion). In `survHE` the posterior median is
used as the plug-in estimate for the log-likelihood, while we use the
posterior mean as per the definition of DIC by (Spiegelhalter et al.
2002), noting that both estimates should be very similar.

## Survival curves implied by Expert Opinion alone

In some situations it may be of interest to see the range of predicted
survival functions given the expert opinion. The easiest solution is to
simulate two or more observations. In order to remove the effect of
these data points we supply the following argument to the
`fit.models.expert` function which essentially sets the likelihood
contribution to zero for these points:

    a0 = rep(0.001,nrow(df1))

Using Stan and JAGS to simulate these “posteriors” is inefficient and
because of lack of identifiability (due to having no data), Markov Chain
Monte Carlo diagnostics will suggest there is a problem. A more
efficient approach for the Weibull distribution is sketched out below
and (similar to (Ouwens 2018)) would be to:

-   Simulate times from the Survival distribution
-   Simulate values of the shape from a vague distribution
-   Reexpress the scale in terms of the shape

As we can see the 90% credible intervals are very wide, narrowing only
at the timepoint at which there is expert opinion.


    nsims <- 10000
    Surv_samp <- rbeta(nsims, 10, 100)
    ancs <- runif(nsims, 0, 10) #shape
    time_expert <- 14
    loc <- exp((ancs*log(time_expert)-log(-log(Surv_samp)))/ancs)

    time <- c(0:20)
    res <- cbind(ancs,loc)

    St <- apply(res, 1, function(x){pweibull(time,x[1],x[2], lower.tail = F )})
    St_sum <- apply(St, 1, quantile, probs = c(0.1, 0.9), na.rm = T)

    plot(y = St_sum[1,], x = time, type= "l", xlab = "Time", ylab = "St",
    main = "90% interval for survival with St from Beta(10,100)")
    lines(y = St_sum[2,], x = time )

<div class="figure">

<img src="Survival without Data.png" alt="Predicted survival without data" width="480" />
<p class="caption">
Predicted survival without data
</p>

</div>

## Model Diagnostics

As this is a Bayesian analysis convergence diagnostics should be
performed. Poor convergence can be observed for many reasons, however,
because of our use of expert opinion my be a symptom of conflict between
the observed data and the expert’s opinion.

Default priors should work in most situations, but still need to be
considered. At a minimum the Bayesian results without expert opinion
should be compared against the maximum likelihood estimates. If
considerable differences are present the prior distributions should be
investigated.

If there is concern about the impact of the prior the expert opinion implemented through penalized maximum is an alternative approach.

Because the analysis is done in JAGS and Stan we can leverage the
`ggmcmc` package:

    #For Stan Models # Log-Normal, RP, Exponential, Weibull
    ggmcmc(ggs(as.mcmc(example1$models$`Gen. Gamma`)), file = "Gengamma.pdf")

    #For JAGS Models # Gamma, Gompertz, Generalized Gamma
    ggmcmc(ggs(as.mcmc(example1$models$`Gamma`)), file = "Gamma.pdf")

## References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-Baio.2020" class="csl-entry">

Baio, Gianluca. 2020. “<span class="nocase">survHE</span>: Survival
Analysis for Health Economic Evaluation and Cost-Effectiveness
Modeling.” *Journal of Statistical Software* 95 (14): 1–47.
<https://doi.org/10.18637/jss.v095.i14>.

</div>

<div id="ref-OHagan.2019" class="csl-entry">

O’Hagan, Anthony. 2019. “Expert Knowledge Elicitation: Subjective but
Scientific.” *The American Statistician* 73 (sup1): 69–81.
<https://doi.org/10.1080/00031305.2018.1518265>.

</div>

<div id="ref-SHELF.2021" class="csl-entry">

Oakley, Jeremy. 2021. *SHELF: Tools to Support the Sheffield Elicitation
Framework*. <https://CRAN.R-project.org/package=SHELF>.

</div>

<div id="ref-Ouwens.2018" class="csl-entry">

Ouwens, Mario. 2018. “Use of Clinical Opinion in the Estimation of
Survival Extrapolation Distributions.” *ISPOR EU 2018 - Use of Clinical
Opinion in the Estimation of Survival Extrapolation Distributions*.
ISPOR.
<https://www.ispor.org/docs/default-source/presentations/91714pdf.pdf?sfvrsn=a5b2756f_0>.

</div>

<div id="ref-Spiegelhalter.2003" class="csl-entry">

Spiegelhalter, David J., Nicola G. Best, Bradley P. Carlin, and Angelika
Van Der Linde. 2002. “Bayesian Measures of Model Complexity and Fit.”
*Journal of the Royal Statistical Society: Series B (Statistical
Methodology)* 64 (4): 583–639.
https://doi.org/<https://doi.org/10.1111/1467-9868.00353>.

</div>

</div>
