

BGReg.solnp <- function(y, X, control = list(print_level = 1,
                                       pars.seed = 113,
                                       pars.init.n = 5, 
                                       tol.rel = 1.0e-8,
                                       print_level = 1,
                                       rho = 0.1)){
  # y: observations
  # X: design matrix
  
  print_level <- control$print_level
  n.init <- control$pars.init.n # 初始化参数的组数
  pars.seed <- control$pars.seed
  tol.rel <- control$tol.rel
  rho <- control$rho
  
  k <- dim(X)[2] # 参数个数
  
  mod.init <- pars.initialization(y = y, X = X)
  beta.init <- mod.init$beta.init
  
  # p1, p2, sigma1, sigma2, delta1, delta2
  # alpha1, alpha2, alphaa3,lpha4, alpha5, alpha6
  # sigma <- p*tau
  # delta <- p*nu 
  # for instance, 
  # sigma1 <- p1*tau1; sigma2 <- p2*tau2
  # delta1 <- p1*nu1; delta2 <- p2*nu2
  # nu1 <- 1; nu2 <- 0.5
  
  LLcomGB2.Reg <-  function(pars) {
    k <- dim(X)[2]
    beta <- pars[1:k]
    
    mu2 <- exp(X%*%beta)
    p1 <- exp(pars[(k+1)]);
    p2 <- exp(pars[(k+2)]);
    sigma1 <- exp(pars[(k+3)])
    sigma2 <- exp(pars[(k+4)])
    delta1 <- p1 # > 1
    delta2 <- p2/2 # > 1
    
    nu1 <- (delta1)/p1; tau1 <- sigma1/p1
    nu2 <- (delta2)/p2; tau2 <- sigma2/p2
    
    
    loglike <- dcomGB2(y = y,
                          mu2 = as.vector(mu2),
                          p1 = p1,
                          p2 = p2, tau1 = tau1, 
                          tau2 = tau2, nu1 = nu1, nu2 = nu2, log = T)
    
    ll <- -sum(loglike)
    if(is.nan(ll) |is.na(ll) | is.infinite(ll)) ll <- 100000000000
    
    return(ll)
  }
  
  
  ieqnPars <- function(pars){
    k <- dim(X)[2]
    
    p1 <- exp(pars[(k+1)]);
    p2 <- exp(pars[(k+2)]);
    sigma1 <- exp(pars[(k+3)])
    sigma2 <- exp(pars[(k+4)])
    delta1 <- p1 # > 1
    delta2 <- p2/2 # > 1
    
    z2 <- pars[(k+1)] - log(1)# >= 0 
    z3 <- pars[(k+2)] - log(2)
    return(c(z2, z3))
  }
  alpha.init.mat <- matrix(0, nrow = 4, ncol = n.init)
  row.names(alpha.init.mat) <- c('p1', 'p2', 
                                 'sigma1', 'sigma2')
  opt.result <- rep(0, n.init) # 定义存储优化结果的向量
  mReg.list <- list()
  for (j in 1:n.init) {
    # 生成随机数
    random_nums <- rnorm(4, 0, 1)
    # 判断条件是否符合要求
    while ((random_nums[1] < log(1))|(random_nums[2] < log(2))) {
      # 如果不符合要求，重新生成随机数
      random_nums <-  rnorm(4, 0, 1)
    }
    alpha.init <- random_nums
    alpha.init.mat[, j] <- alpha.init
    pars.int <- c(beta.init, alpha.init)
    ieqnPars(pars.int)
    LLcomGB2.Reg(pars.int)
    tryCatch({
      mReg <- solnp(pars = pars.int,
                    fun =  LLcomGB2.Reg,
                    ineqfun = ieqnPars,
                    ineqLB = c(0, 0),
                    ineqUB = c(10, 10),
                    LB = c(rep(-10, times = dim(X)[2]), rep(-10, times = 4)),
                    UB = c(rep(10, times = dim(X)[2]), rep(10, times = 4)),
                    control = list(rho = rho,
                                   tol = tol.rel,
                                   trace = print_level))
    }, error = function(e) {
      # 处理错误的代码
      # 可以使用message()函数来输出错误信息
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
  npars <- k + 4
  AIC <- 2*NLL + 2*npars
  BIC <- 2*NLL + npars*log(length(y))
  

  names(est)[1:k] <- paste0('beta', 0:(k-1))
  names(est)[(k+1):(npars)] <- paste0('alpha', 1:4)
  
  out <- list(estimate = est,
              NLL = NLL,
              AIC = AIC, BIC = BIC,
              npars = npars,
              Hessian = mReg$hessian)
  
  
  return(out)
}


# model <- BGReg.solnp(y = y, 
#                        X = X, 
#                        control = list(pars.init.n = 1))
# model$X
# model$estimate
# model$Hessian
# class <- 'original'
summary.BG <- function(model, class = c('original', 
                                          'transformed')){
  pars <- model$estimate
  npars <- length(pars)
  k <- npars - 4
  
  beta <- pars[1:k]
  p1 <- exp(pars[(k+1)]);
  p2 <- exp(pars[(k+2)]);
  sigma1 <- exp(pars[(k+3)])
  sigma2 <- exp(pars[(k+4)])
  delta1 <- p1 # > 1
  delta2 <- p2/2 # > 1
  
  nu1 <- (delta1)/p1; tau1 <- sigma1/p1
  nu2 <- (delta2)/p2; tau2 <- sigma2/p2
  
  Hessian <- model$Hessian[3:(npars + 2), 3:(npars + 2)]
  cov <- (solve(Hessian))
  sd <- sqrt(diag(solve(Hessian)))                 ## standard error 
  
  if(class == 'original'){
    est <- c(beta, log(c(p1, p2, sigma1, sigma2, delta1, delta2)))
    se <- c(sd[1:k], sd[(k+1):(k+4)], sd[k+1], sd[k+2])
    
    zvalue <- est/se                                      ##  Z statistic
    pvalue <- ifelse(zvalue>=0, pnorm(zvalue, lower=F)*2, pnorm(zvalue)*2)           ## p value
    
    summarytable <- data.frame(estimates = est, 
                               se = se,
                               zvalue = zvalue, 
                               pvalue = pvalue)
    row.names(summarytable)[1:k] <- paste0('beta', 0:(k-1))
    row.names(summarytable)[(k+1):length(est)] <- paste0('alpha', 1:6)
    
  } else {
    var.p1 <- p1*cov[1,1]*p1
    var.p2 <- p2*cov[2,2]*p2
    var.tau1 <- t(c(-tau1, tau1))%*%cov[c(1,3), c(1,3)]%*%(c(-tau1, tau1))
    var.tau2 <- t(c(-tau2, tau2))%*%cov[c(2,4), c(2,4)]%*%(c(-tau2, tau2))
    var.nu1 <- NA
    var.nu2 <- NA
    se <- c(sd[1:k], sqrt(c(var.p1, var.p2, 
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
    row.names(summarytable)[(k+1):length(est)] <- c('p1', 'p2', 'tau1', 'tau2', 'nu1', 'nu2')
  }
  
  
  
  out <- list(summarytable = summarytable)
  return(out)
}

# summary.BGReg(model = model, class = 'original')
# summary.BGReg(model = model, class = 'transformed')




