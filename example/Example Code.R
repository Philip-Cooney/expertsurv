library(dplyr)
library(flexsurv)

#https://strimmerlab.github.io/publications/lecture-notes/MATH20802/likelihood-based-confidence-interval-and-likelihood-ratio.html
x_vals <- seq(0, 0.3, by = 0.01)
mu_expert <- 0.2
sd_expert <- 0.1

#plot( x= x_vals, y = dnorm(x_vals, mean = mu_expert, sd = sd_expert, log= F), type = "l", xlab = "Survival", ylab = "Likelihood")
plot( x= x_vals, y = dnorm(x_vals, mean = mu_expert, sd = sd_expert), type = "l", xlab = "Survival", ylab = "Density",
      xaxt='n')
axis(side = 1, at = seq(0, .3, by = 0.05))
axis(side = 3, at = seq(0, .3, by = 0.05), labels = round(-log(seq(0, .3, by = 0.05))/2, 2))
mtext("hazard", side = 3, line = 3) 

lambda_vals <- seq(0, 3, by = 0.001)

set.seed(123)
times <- rexp(10, 1)
n_events <- length(times)

censor.time <- 0.75
df<- data.frame(time_event = times, status = 1) %>% mutate(time = ifelse(time_event <= censor.time,
                                                                         time_event, censor.time),
                                                           status = ifelse(time_event <= censor.time,
                                                                          1,0) )



plot(flexsurv::flexsurvreg(Surv(time, status)~1,data = df, dist = "exp"),t = seq(0, 2, by = 0.1),
     xlim = c(0,2), xlab = "time", ylab = "S(t)")

LL_surv <- dnorm(exp(-lambda_vals*2), mean = mu_expert, sd = sd_expert, log= T)
LL_data <- n_events*log(lambda_vals)-sum(times)*lambda_vals

MLE_bound = max(LL_surv+LL_data) -1.92

plot(x= lambda_vals, y = LL_surv, 
     type = "l", ylab = "Log-Lik", xlab = "hazard", ylim = c(-10,3), col = "blue", xlim = c(0,3))
lines(x= lambda_vals, y = LL_surv+LL_data)
lines(x = lambda_vals,  LL_data, col = "red")
abline(v = length(times)/sum(times), col = "red", lty = 2)
abline(v = round(-log(mu_expert)/2, digits = 2),col = "blue", lty = 2 )
abline(v = lambda_vals[which.max(LL_surv+LL_data)],lty = 2)
#abline(h =MLE_bound ,lty = 2)
legend("topright", legend=c("Data Log-Likelihood", "Expert Contribution", "Total Log-Likelihood"),
       col=c("red", "blue", "black"), lty=1, cex=0.8)


#MLE CI
MLE_val <- max(LL_surv+LL_data)
MLE_lambda <- lambda_vals[which.max(LL_surv+LL_data)]
UpperCI <- min(lambda_vals[which( MLE_val- (LL_surv+LL_data) < 1.92)])
LowerCI <-max(lambda_vals[which( MLE_val- (LL_surv+LL_data) < 1.92)])


plot(flexsurv::flexsurvreg(Surv(time, status)~1,data = df, dist = "exp"),t = seq(0, 2, by = 0.1),
     xlim = c(0,2), xlab = "time", ylab = "S(t)", ci = F)
lines(y = exp(-MLE_lambda*seq(0, 2, by = 0.1)), x = seq(0, 2, by = 0.1), col = "blue")
lines(y = exp(-UpperCI*seq(0, 2, by = 0.1)), x = seq(0, 2, by = 0.1), col = "blue", lty = 2)
lines(y = exp(-LowerCI*seq(0, 2, by = 0.1)), x = seq(0, 2, by = 0.1), col = "blue", lty = 2)
legend("topright", legend=c("S(t) - Data Only", "S(t) - Including Expert Opinion"),
       col=c("red", "blue"), lty=1, cex=0.8)


con_int <- sfsmisc::integrate.xy(fx = exp(LL_surv+LL_data), x = lambda_vals)
#Calculate Bayesian CI
cum_prob <- rep(NA, length(lambda_vals))
for(i in 2:length(cum_prob)){
    cum_prob[i] <- sfsmisc::integrate.xy(fx = exp(LL_surv+LL_data)[1:i], x = lambda_vals[1:i])/con_int
  
}

plot(y = exp(LL_surv+LL_data)/con_int, x= lambda_vals, type = "l", ylab = "density", xlab = "lambda")
abline(v = c(lambda_vals[which.min(abs(cum_prob-0.025))],lambda_vals[which.min(abs(cum_prob-0.975))]), lty = 2)

#Bayes CI
lambda_vals[which.min(abs(cum_prob-0.025))]
lambda_vals[which.min(abs(cum_prob-0.975))]
# 
# shape_eval <- 2 #alpha
# rate_eval <- 2 #beta
# lambda_vals <- seq(0, 40, by = 0.001)
# 
# #mean(rgamma(1000, shape=shape_eval,rate = rate_eval))
# 
# LL_eval<- shape_eval*log(lambda_vals)-rate_eval*lambda_vals
# MLE_eval <- lambda_vals[which.max(LL_eval)]
# 
# plot(LL_eval , x = lambda_vals, type = "l", xlim = c(0, qgamma(.99, shape=shape_eval ,rate = rate_eval )))
# abline(h = max(LL_eval)-1.92)
# MLE_eval <- lambda_vals[which.max(LL_eval)]
# 
# min(lambda_vals[LL_eval>(max(LL_eval)-1.92)])
# max(lambda_vals[LL_eval>(max(LL_eval)-1.92)])
# qgamma(c(0.025,0.975), shape=shape_eval ,rate = rate_eval )
# 
# plot(lambda_vals, y = dgamma(lambda_vals, shape=shape_eval ,rate = rate_eval ))


