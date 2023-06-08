rule_stg1 = function(phi0,phi1,alpha=0.05, beta=0.2,
                     pi=0,surv_dist, tau=2){
  
  # phi0 is median survival time under null hypothesis
  # phi1 is median survival time under alternative hypothesis
  # pi is probability having phi1 based on prior knowledge
  # surv_dist is the survival distribution specified: exponential, weibull,uniform
  # tau is the shape parameter of the weibull distribution
  
  if (surv_dist == "exponential"){
    
    ## density function of exponential distribution
    f_y = function(phi){
      (log(2)/phi)*exp(-phi*log(2)/phi)
    }
    
  } else if (surv_dist == "weibull"){
    ## density function of weibull distribution with shape as 2
    f_y = function(phi){
      lambda = phi/((log(2))^(1/tau))
      (tau/lambda)*(phi/lambda)^(tau-1)*exp(-(phi/lambda)^tau)
      
    }
  } else if (surv_dist == "uniform"){
    
    ## density function of uniform distribution with (0,2*phi)
    f_y = function(phi){
      1/(2*phi)
    }
    
  }
  
  ## calculate sample size according to the formula
  n = ceiling((f_y(phi1)/f_y(phi0)*qnorm(p=alpha,lower.tail=FALSE)+qnorm(p=beta,lower.tail=FALSE))^2*((0.5/(f_y(phi1)*(phi1 - phi0)))^2))
  t = 0.5*qnorm(p=alpha,lower.tail=FALSE)/(sqrt(n)*f_y(phi0))+phi0
  
  res <- list(n=n, t=t)
  
  ## return the decision boundaries and sample size for one set of alpha1 and beta1
  return(res) 
} 
