METTSS_rule = function(phi0,phi1,alpha=0.05, beta=0.2,
                       alpha1=seq(0.05, 0.5, 0.005),beta1=seq(0.001, 0.199, 0.005),
                       pi=0,surv_dist, tau=2){
  
  # phi0 is median survival time under null hypothesis
  # phi1 is median survival time under alternative hypothesis
  # alpha is overall type I error
  # beta is overall type II error
  # alpha1 is a list of type I error for first stage to search
  # beta1 is a list of type II error for first stage to search
  # pi is probability having phi1 based on prior knowledge
  # surv_dist is the survival distribution specified: exponential, weibull,uniform
  # tau is the shape parameter of the weibull distribution
  
  result <- En_fun_search(phi0,phi1,surv_dist="weibull")
  n1 <- result$n1
  n2 <- result$n2
  t1 <- result$t1
  t2 <- result$t2
  result2 <- rule_stg1(phi0,phi1,alpha=0.05, beta=0.2,
                       pi=0,surv_dist="weibull")
  tstar <- result2$t
  
  res <-list(t1, n1, t2, n2, tstar)
  names(res) <- c("t1", "n1", "t2", "n2", "tstar")
  return(res)
  
}