# 1. comGB2 regression model -------------------
library(Rsolnp)
comGB2Reg.solnp <- function(y, X, control = list(print_level = 1,
                                     pars.seed = 113,
                                     pars.init.n = 5, 
                                     tol.rel = 1.0e-8,
                                     rho = 0.1)){
  # y: observations
  # X: design matrix
  
  print_level <- control$print_level
  n.init <- control$pars.init.n # Number of groups for initialization parameters
  pars.seed <- control$pars.seed
  tol.rel <- control$tol.rel
  rho <- control$rho

  k <- dim(X)[2] # Numbers of parameters
  
  mod.init <- pars.initialization(y = y, X = X)
  beta.init <- mod.init$beta.init
  
  # p1, p2, sigma1, sigma2, delta1, delta2
  # alpha1, alpha2, alpha3, alpha4, alpha5, alpha6
  # sigma <- p*tau
  # delta <- p*nu 
  # for instance, 
  # sigma1 <- p1*tau1; sigma2 <- p2*tau2
  # delta1 <- p1*nu1; delta2 <- p2*nu2
  LLcomGB2.Reg <-  function(pars) {
    k <- dim(X)[2]
    beta <- pars[1:k]
    
    mu2 <- exp(X%*%beta)
    p1 <- exp(pars[(k+1)]);
    p2 <- exp(pars[(k+2)]);
    sigma1 <- exp(pars[(k+3)])
    sigma2 <- exp(pars[(k+4)])
    delta1 <- exp(pars[(k+5)]) # > 1
    delta2 <- exp(pars[(k+6)]) # > 1
    
    nu1 <- (delta1)/p1; tau1 <- sigma1/p1
    nu2 <- (delta2)/p2; tau2 <- sigma2/p2
    
    
    loglike <- dcomGB2(y = y,
                       mu2 = as.vector(mu2),
                       p1 = p1,
                       p2 = p2, tau1 = tau1, tau2 = tau2, nu1 = nu1, nu2 = nu2, log = T)
    ll <- -sum(loglike)
    if(is.na(ll) | is.infinite(ll)) ll <- 100000000000
    
    return(ll)
  }
  
  
  ieqnPars <- function(pars){
    k <- dim(X)[2]
    p1 <- exp(pars[k+1])
    p2 <- exp(pars[k+2])
    sigma1 <- exp(pars[(k+3)])
    sigma2 <- exp(pars[(k+4)])
    delta1 <- exp(pars[(k+5)]) # > 1
    delta2 <- exp(pars[(k+6)]) # > 1
    
    z2 <- pars[(k+5)] # >= 0 
    z3 <- pars[(k+6)] # >= 0 
    return(c(z2, z3))
  }
  alpha.init.mat <- matrix(0, nrow = 6, ncol = n.init)
  row.names(alpha.init.mat) <- c('p1', 'p2', 
                                 'sigma1', 'sigma2', 
                                 'delta1', 'delta2')
  opt.result <- rep(0, n.init) # Define the vector for storing optimization results
  mReg.list <- list()
  set.seed(pars.seed)
  for (j in 1:n.init) {
    # Generate random numbers
    random_nums <- rnorm(6, 0, 1)
    # random_nums <- runif(6, min = -2, max = 2)
    # Determine if the conditions are met
    while ((random_nums[5] < 0) | (random_nums[6] < 0)) {
      # If it does not meet the requirements, regenerate the random number
      random_nums <- rnorm(6, 0, 1)
      # random_nums <- runif(6, min = -2, max = 2)
    }
    alpha.init <- random_nums
    alpha.init.mat[, j] <- alpha.init
    pars.int <- c(beta.init, alpha.init)
    tryCatch({
      suppressWarnings(
      mReg <- solnp(pars = pars.int,
                    fun =  LLcomGB2.Reg,
                    ineqfun = ieqnPars,
                    ineqLB = c(0, 0),
                    ineqUB = c(10, 10),
                    LB = c(rep(-10, times = dim(X)[2]), rep(-10, times = 6)),
                    UB = c(rep(10, times = dim(X)[2]), rep(10, times = 6)),
                    control = list(rho = rho,
                                   tol = tol.rel,
                                  trace = print_level))
      )
      
      
    }, error = function(e) {
      # The code to handle errors
      # You can use the message() function to output error messages
      message('Perform the next parameter initialization')
    })
    mReg.list[[j]] <- mReg
    opt.result[j] <- mReg$values[length(mReg$values)]
  }
  best.alpha.init <- as.vector(alpha.init.mat[, which.min(opt.result)])
  pars.int <- c(beta.init, best.alpha.init)
  mReg <- mReg.list[[which.min(opt.result)]]
  
  
  est <- mReg$pars
  NLL <- mReg$values[length(mReg$values)]
  npars <- k + 6
  AIC <- 2*NLL + 2*npars
  BIC <- 2*NLL + npars*log(length(y))

  names(est)[1:k] <- paste0('beta', 0:(k-1))
  names(est)[(k+1):(npars)] <- paste0('alpha', 1:6)
  
  out <- list(estimate = est,
              Hessian = mReg$hessian,
              NLL = NLL,
              AIC = AIC, BIC = BIC,
              npars = npars)
  return(out)
  
}



summary.comGB2 <- function(model, class = c('original', 
                                            'transformed')){
  pars <- model$estimate
  npars <- length(pars)
  k <- npars - 6

  beta <- pars[1:k]
  p1 <- exp(pars[(k+1)]);
  p2 <- exp(pars[(k+2)]);
  sigma1 <- exp(pars[(k+3)])
  sigma2 <- exp(pars[(k+4)])
  delta1 <- exp(pars[(k+5)]) # > 1
  delta2 <- exp(pars[(k+6)]) # > 1
  
  nu1 <- (delta1)/p1; tau1 <- sigma1/p1
  nu2 <- (delta2)/p2; tau2 <- sigma2/p2
  
  Hessian <- model$Hessian[3:(npars + 2), 3:(npars + 2)]
  cov <- (solve(Hessian))
  se <- sqrt(diag(solve(Hessian)))                 ## standard error 
  if(class == 'original'){
    est <- c(beta, log(c(p1, p2, sigma1, sigma2, delta1, delta2)))
    zvalue <- est/se                                      ##  Z statistic
    pvalue <- ifelse(zvalue>=0, pnorm(zvalue, lower=F)*2, pnorm(zvalue)*2)           ## p value
    
    summarytable <- data.frame(estimates = est, 
                               se = se,
                               zvalue = zvalue, 
                               pvalue = pvalue)
    row.names(summarytable)[1:k] <- paste0('beta', 0:(k-1))
    row.names(summarytable)[(k+1):length(pars)] <- paste0('alpha', 1:6)
    
  } else {
    var.p1 <- p1*cov[1,1]*p1
    var.p2 <- p2*cov[2,2]*p2
    var.tau1 <- t(c(-tau1, tau1))%*%cov[c(1,3), c(1,3)]%*%(c(-tau1, tau1))
    var.tau2 <- t(c(-tau2, tau2))%*%cov[c(2,4), c(2,4)]%*%(c(-tau2, tau2))
    var.nu1 <- t(c(-nu1, nu1))%*%cov[c(1,5), c(1,5)]%*%(c(-nu1, nu1))
    var.nu2 <- t(c(-nu2, nu2))%*%cov[c(2,6), c(2,6)]%*%(c(-nu2, nu2))
    se <- c(se[1:k], sqrt(c(var.p1, var.p2, 
                            var.tau1, var.tau2, 
                            var.nu1, var.nu2)))
    est <- c(beta, (c(p1, p2, tau1, tau2, nu1, nu2)))
    zvalue <- est/se                                            ##  Z statistic
    pvalue <- ifelse(zvalue>=0, pnorm(zvalue, lower=F)*2, pnorm(zvalue)*2)           ## p value
    
    summarytable <- data.frame(estimates = est, 
                               se = se,
                               zvalue = zvalue, 
                               pvalue = pvalue)
    
    row.names(summarytable)[1:k] <- paste0('beta', 0:(k-1))
    row.names(summarytable)[(k+1):length(pars)] <- c('p1', 'p2', 'tau1', 'tau2', 'nu1', 'nu2')
  }
  

  
  out <- list(summarytable = summarytable)
  return(out)
}



