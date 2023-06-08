En_fun_search = function(phi0,phi1,alpha=0.05, beta=0.2,
                         alpha1=seq(0.05, 0.5, 0.005),beta1=seq(0.001, 0.199, 0.005),
                         pi=0,surv_dist, tau=2){
  
  # phi0 is median survival time under null hypothesis
  # phi1 is median survival time under alternative hypothesis
  # alpha1 is a list of type I error for first stage to search
  # beta1 is a list of type II error for first stage to search
  # pi is probability having phi1 based on prior knowledge
  # surv_dist is the survival distribution specified: exponential, weibull,uniform
  # tau is the shape parameter of the weibull distribution
  
  ## use En_fun to calculate for each alpha1 and beta1
  outcome = matrix(unlist(lapply(alpha1, function(x) lapply(beta1, function(y) En_fun(phi0,phi1,alpha, beta, x,y,pi,surv_dist, tau=2)))),ncol = 10,byrow=TRUE)
  outcome <- data.frame(outcome)
  colnames(outcome) <- c("phi0","phi1","n1", "t1", "n2", "t2","n", "EN", "alpha1", "beta1") 
  result = outcome %>% filter(n2>0,t1>0,t2>0) %>% group_by(phi0,phi1) %>% slice(which.min(EN))
  res <- list(phi0=result$phi0,phi1=result$phi1,
              n1=result$n1, t1=result$t1, n2=result$n2, t2=result$t2, 
              n=result$n,EN0=result$EN,alpha1=result$alpha1,beta1=result$beta1)
  return(res)
}
