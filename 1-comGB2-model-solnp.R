# 1. comGB2 regression model -------------------
library(Rsolnp) # Load the Rsolnp package for optimization

# This is R code for estimating a regression model using the comGB2 distribution. 
# The comGB2 distribution is a seven-parameter distribution and is often used to model data with heavy tails and skewness.
# The code starts by defining a log-likelihood function (LLcomGB2.Reg) that takes a set of parameters as input and returns the negative log-likelihood of the data given those parameters. 
# The function uses the dcomGB2 function to calculate the log-likelihood of the comGB2 distribution given the data, the predicted values from the regression model (mu2), and the parameter values.
# The code then defines a function (ieqnPars) that imposes inequality constraints on the parameters. These constraints ensure that certain parameters (e.g., delta1 and delta2) are positive. The constraints are implemented using the solnp function, which is a nonlinear optimization solver that allows for inequality constraints.
# The code sets the lower and upper bounds for the parameters using the LB and UB vectors. It then initializes the parameter vector (pars) to the values in pars.int and uses the solnp function to find the parameter values that minimize the negative log-likelihood. The resulting parameter estimates are stored in the mReg$pars vector.
# Finally, the code extracts the estimated beta coefficients (which are the first k elements of the pars vector) and the estimated comGB2 distribution parameters (which are the last six elements of the pars vector). These estimates are returned as a single vector (est).


comGB2Reg.solnp <- function(y, X, control = list(print_level = 1,
                                     pars.seed = 113,
                                     pars.init.n = 5, 
                                     tol.rel = 1.0e-8,
                                     rho = 0.1)){
  # y: observations
  # X: design matrix
  
  print_level <- control$print_level # Get print_level from the control list
  n.init <- control$pars.init.n # Get the number of groups for initialization parameters
  pars.seed <- control$pars.seed # Get the pars.seed from the control list
  tol.rel <- control$tol.rel # Get the relative tolerance from the control list
  rho <- control$rho # Get the rho from the control list
  
  k <- dim(X)[2] # Get the number of parameters
  
  mod.init <- pars.initialization(y = y, X = X) # Get the initial values of parameters using the pars.initialization function
  beta.init <- mod.init$beta.init # Get the initial value of beta
  
  # p1, p2, sigma1, sigma2, delta1, delta2
  # alpha1, alpha2, alpha3, alpha4, alpha5, alpha6
  # sigma <- p*tau
  # delta <- p*nu 
  # for instance, 
  # sigma1 <- p1*tau1; sigma2 <- p2*tau2
  # delta1 <- p1*nu1; delta2 <- p2*nu2
  
  LLcomGB2.Reg <-  function(pars) {
    # Define a function for computing the log-likelihood of the comGB2 model with parameters beta, alpha1, alpha2, alpha3, alpha4, alpha5, and alpha6
    k <- dim(X)[2] # Get the number of parameters
    beta <- pars[1:k] # Get beta from the parameter vector
    
    mu2 <- exp(X%*%beta) # Calculate mu2 based on the formula
    p1 <- exp(pars[(k+1)]); # Get p1 from the parameter vector
    p2 <- exp(pars[(k+2)]); # Get p2 from the parameter vector
    sigma1 <- exp(pars[(k+3)]) # Get sigma1 from the parameter vector
    sigma2 <- exp(pars[(k+4)]) # Get sigma2 from the parameter vector
    delta1 <- exp(pars[(k+5)]) # > 1; Get delta1 from the parameter vector
    delta2 <- exp(pars[(k+6)]) # > 1; Get delta2 from the parameter vector
    
    nu1 <- (delta1)/p1; tau1 <- sigma1/p1 # Calculate nu1 and tau1 based on the formulas
    nu2 <- (delta2)/p2; tau2 <- sigma2/p2 # Calculate nu2 and tau2 based on the formulas
    
    loglike <- dcomGB2(y = y, # Compute the log-likelihood of the comGB2 model based on the given parameters
                       mu2 = as.vector(mu2),
                       p1 = p1,
                       p2 = p2, tau1 = tau1, tau2 = tau2, nu1 = nu1, nu2 = nu2, 
                       log = T)
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


# This function calculates and returns a summary table of the estimated parameters for the comGB2 model. The input parameter model is an object returned by the comGB2 function. The class parameter can be set to either 'original' or 'transformed' to indicate whether to display the estimates and standard errors in the original or transformed scale.
# The function first extracts the parameter estimates and standard errors from the model object, and then calculates the z-values and p-values for each parameter. The z-value is calculated as the ratio of the estimate to the standard error, and the p-value is calculated using a two-sided normal distribution test.
# If class is set to 'original', the function formats the output table to display the estimates, standard errors, z-values, and p-values for the beta and alpha parameters. If class is set to 'transformed', the function additionally calculates the variance-covariance matrix of the transformed parameters and displays the estimates, standard errors, z-values, and p-values for the transformed parameters.
# The output of the function is a list with one element named summarytable, which contains the summary table as a data frame. The row names of the table correspond to the parameter names.

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



