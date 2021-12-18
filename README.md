
<!-- README.md is generated from README.Rmd. Please edit that file -->

# expert-surv

<!-- badges: start -->
<!-- badges: end -->

The goal of `expertsurv` is to incorporate expert opinion into an
analysis of time to event data. `expertsurv` uses many of the core
functions of the `survHE` package (Baio 2020). Technical details of the
implementation are detailed in (Cooney and White 2021) and will not be
repeated here.

The key function is `fit.models.expert` and operates almost identically
to the `fit.models` function of `survHE`.

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
                                         wi = c(1,1),
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
```

For details on assigning distributions to elicited probabilities and
quantiles see the `SHELF` package (Oakley 2021) and for an overview on
methodological approaches to eliciting expert opinion see (O’Hagan
2019).


    data2 <- survHE::data %>% rename(status = censored) %>% mutate(time2 = ifelse(time > 10, 10, time),
                                                                  status2 = ifelse(time> 10, 0, status))

    #Set the opinion type to "survival"

    example1  <- expertsurv:::fit.models.expert(formula=Surv(time2,status2)~1,data=data2,
                                            distr=c("wei"),
                                            method="hmc",
                                            iter = 5000,
                                            opinion_type = "survival",
                                            times_expert = timepoint_expert, 
                                            param_expert = param_expert_example1)

     plot(example1, add.km = T, t = 0:30)+
      theme_light()+
      scale_x_continuous(expand = c(0, 0), limits = c(0,NA), breaks=seq(0, 30, 2)) + 
      scale_y_continuous(expand = c(0, 0), limits = c(0, NA), breaks=seq(0, 1, 0.05))

<div class="figure">

<img src="Vignette_Example 1.png" alt="Survival function with Expert prior" width="2025" />
<p class="caption">
Survival function with Expert prior
</p>

</div>

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

    survHE.data.model  <- expertsurv:::fit.models.expert(formula=Surv(time2,status2)~as.factor(arm),data=data2,
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

<img src="Vignette_Example 2.png" alt="Survival function with Expert prior (left) and Vague prior (right)" width="3000" />
<p class="caption">
Survival function with Expert prior (left) and Vague prior (right)
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

    #Prior belief of 5 "months" difference in expected survival
    param_expert3[[1]] <- data.frame(dist = "norm", wi = 1, param1 = 5, param2 = 0.2, param3 = NA)


    survHE.data.model  <- expertsurv:::fit.models.expert(formula=Surv(time2,status2)~as.factor(arm),data=data2,
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

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-Baio.2020" class="csl-entry">

Baio, Gianluca. 2020. “<span class="nocase">survHE</span>: Survival
Analysis for Health Economic Evaluation and Cost-Effectiveness
Modeling.” *Journal of Statistical Software* 95 (14): 1–47.
<https://doi.org/10.18637/jss.v095.i14>.

</div>

<div id="ref-Cooney.2021" class="csl-entry">

Cooney, Philip, and Arthur White. 2021. “Utilizing Expert Opinion to
Inform Extrapolation of Survival Models.”
<http://arxiv.org/abs/2112.02288>.

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

</div>
