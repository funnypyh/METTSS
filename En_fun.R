## Calculate decision boundaries and sample size
En_fun = function(phi0,phi1,alpha=0.05, beta=0.2,
                  alpha1=seq(0.05, 0.5, 0.005),beta1=seq(0.001, 0.199, 0.005),
                  pi=0,surv_dist, tau=2){
  
  # phi0 is median survival time under null hypothesis
  # phi1 is median survival time under alternative hypothesis
  # alpha1 is one type I error for first stage
  # beta1 is one type II error for first stage
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
  n1 = ceiling((f_y(phi1)/f_y(phi0)*qnorm(p=alpha1,lower.tail=FALSE)+qnorm(p=beta1,lower.tail=FALSE))^2*((0.5/(f_y(phi1)*(phi1 - phi0)))^2))
  t1 = 0.5*qnorm(p=alpha1,lower.tail=FALSE)/(sqrt(n1)*f_y(phi0))+phi0
  n2 = ceiling((f_y(phi1)/f_y(phi0)*qnorm(p=alpha,lower.tail=FALSE)+qnorm(p=beta - beta1,lower.tail=FALSE))^2*((0.5/(f_y(phi1)*(phi1 - phi0)))^2) - n1)
  t2 = 0.5*qnorm(p=alpha,lower.tail=FALSE)/(sqrt(n1+n2)*f_y(phi0))+phi0
  EN = n1 + n2*(pi*(1-beta1)+(1-pi)*alpha1)
  n = n1+n2
  res <- list(phi0,phi1,n1, t1, n2, t2, n,EN,alpha1,beta1)
  
  ## return the decision boundaries and sample size for one set of alpha1 and beta1
  return(res) 
} 
