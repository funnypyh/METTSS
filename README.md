# METTSS
Sample size determination for two-stage adaptive designs of single arm clinical trials based on median event time test

# Example
library(survival)

library(dplyr)

set.seed(1234)

out <- METTSS_rule(phi0=10, phi1=17, alpha=0.05, beta=0.2, pi=0, surv_dist="weibull")

out
