library(gamlss)
library(GB2)

# Parameter initialization strategy
pars.initialization <- function(y, X){
  mod <- gamlss(y ~ X - 1, family = GA(), method = mixed(2, 100000),
                control = gamlss.control(trace = FALSE)) # fit the Gamma regression
  beta.init <- coefficients(mod, what = 'mu') # extract the estimates
  out <- list(beta.init = beta.init)
  return(out)
}

# define of the composite GB2 distribution --------------------------------------
# 1. density function (using GB2 package)
dcomGB2 <- function(y, mu2, p1, p2 , nu1, nu2, tau1, tau2, log = FALSE){
  g1 <- (p1*nu1 - 1)/(p1*tau1 + 1)
  g2 <- (p2*nu2 - 1)/(p2*tau2 + 1)
  x0 <- mu2*g2^(1/p2)
  mu1 <- mu2*g2^(1/p2)/g1^(1/p1)
  
  f_1 <- dgb2(x = x0, shape1 = p1, scale = mu1, shape2 = nu1, shape3 = tau1)
  F_1 <- pgb2(x = x0, shape1 = p1, scale = mu1, shape2 = nu1, shape3 = tau1)
  f_2 <- dgb2(x = x0, shape1 = p2, scale = mu2, shape2 = nu2, shape3 = tau2)
  F_2 <- pgb2(x = x0, shape1 = p2, scale = mu2, shape2 = nu2, shape3 = tau2)
 
  r <- f_2*F_1*(f_2*F_1+f_1*(1-F_2))^(-1)
  
  index1 <- (y <= x0)
  index2 <- y>x0
  fpart1 <- r*dgb2(x = y, shape1 = p1, scale = mu1, shape2 = nu1, shape3 = tau1)/F_1
  fpart2 <- (1-r)*dgb2(x = y, shape1 = p2, scale = mu2, shape2 = nu2, shape3 = tau2)/(1-F_2)
  
  out <- fpart1 # y <= x0
  out[index2] <- fpart2[index2]  # y > x0
  
  if(log == TRUE){ll = log(out)} else {ll = out}
  return(ll)
}

# 2. distribution function (using gamlss package)
pcomGB2 <- function(y, mu2, p1, p2, nu1, nu2, tau1, tau2){
  g1 <- (p1*nu1 - 1)/(p1*tau1 + 1)
  g2 <- (p2*nu2 - 1)/(p2*tau2 + 1)
  x0 <- mu2*g2^(1/p2)
  mu1 <- mu2*g2^(1/p2)/g1^(1/p1)
  f_1 <- dGB2(x = x0, mu = mu1, sigma = p1, nu = nu1, tau = tau1)
  F_1 <- pGB2(q = x0, mu = mu1, sigma = p1, nu = nu1, tau = tau1)
  f_2 <- dGB2(x = x0, mu = mu2, sigma = p2, nu = nu2, tau = tau2)
  F_2 <- pGB2(q = x0, mu = mu2, sigma = p2, nu = nu2, tau = tau2)
  
  r <- f_2*F_1*(f_2*F_1+f_1*(1-F_2))^(-1)
  
  if(y <= x0){
    Fpart <- r*pGB2(q = y, mu = mu1, sigma = p1, nu = nu1, tau = tau1)/F_1
  } else {
    Fpart <- r + (1-r)*(pGB2(q = y, mu = mu2, sigma = p2, nu = nu2, tau = tau2) - pGB2(q = x0, mu = mu2, sigma = p2, nu = nu2, tau = tau2))/(1-F_2)
  }
  out <- Fpart
  return(out)
}
pcomGB2 <- Vectorize(pcomGB2)

# 3. qunatile function of composite GB2 distribution
qcomGB2 <- function(p, mu2, p1, p2 , nu1, nu2, tau1, tau2){
  g1 <- (p1*nu1 - 1)/(p1*tau1 + 1)
  g2 <- (p2*nu2 - 1)/(p2*tau2 + 1)
  x0 <- mu2*g2^(1/p2)
  mu1 <- mu2*g2^(1/p2)/g1^(1/p1)
  f_1 <- dGB2(x = x0, mu = mu1, sigma = p1, nu = nu1, tau = tau1)
  F_1 <- pGB2(q = x0, mu = mu1, sigma = p1, nu = nu1, tau = tau1)
  f_2 <- dGB2(x = x0, mu = mu2, sigma = p2, nu = nu2, tau = tau2)
  F_2 <- pGB2(q = x0, mu = mu2, sigma = p2, nu = nu2, tau = tau2)
  
  r <- f_2*F_1*(f_2*F_1+f_1*(1-F_2))^(-1)
  
  if(p <= r){
    u1sim <- p/r*F_1
    Isim <- qbeta(p = u1sim, shape1 = nu1, shape2 = tau1)
    zsim <- Isim/(1 - Isim)
    out <- mu1*zsim^(1/p1)
    # out <- qGB2(p = u1sim, mu = mu1, sigma = p1, nu = nu1, tau = tau1)
  } else if (p > r) {
    u2sim <- F_2 + (p - r)/(1 - r)*(1 - F_2)
    out <- qGB2(p = u2sim, mu = mu2, sigma = p2, nu = nu2, tau = tau2)
  }
  
  return(out)
}
qcomGB2 <- Vectorize(qcomGB2)

# 3. random generation for the ComGB2 distribution
rcomGB2 <- function(n, mu2, p1, p2 , nu1, nu2, tau1, tau2){
  g1 <- (p1*nu1 - 1)/(p1*tau1 + 1)
  g2 <- (p2*nu2 - 1)/(p2*tau2 + 1)
  x0 <- mu2*g2^(1/p2)
  mu1 <- mu2*g2^(1/p2)/g1^(1/p1)
  f_1 <- dGB2(x = x0, mu = mu1, sigma = p1, nu = nu1, tau = tau1)
  F_1 <- pGB2(q = x0, mu = mu1, sigma = p1, nu = nu1, tau = tau1)
  f_2 <- dGB2(x = x0, mu = mu2, sigma = p2, nu = nu2, tau = tau2)
  F_2 <- pGB2(q = x0, mu = mu2, sigma = p2, nu = nu2, tau = tau2)
  
  r <- f_2*F_1*(f_2*F_1+f_1*(1-F_2))^(-1)
  
  qsim <- runif(n = n, min = 0, max = 1)
  u1sim <- qsim/r*F_1
  u2sim <- F_2 + (qsim - r)/(1 - r)*(1 - F_2)
  r.vector <- rep(r, times = n)
  cond1 <- qsim <= r.vector
  cond2 <- qsim > r.vector
  y1sim <- qGB2(p = u1sim, mu = mu1, sigma = p1, nu = nu1, tau = tau1)
  y2sim <- qGB2(p = u2sim, mu = mu2, sigma = p2, nu = nu2, tau = tau2)
  
  ysim <- y1sim
  ysim[cond2] <- y2sim[cond2]
  
  return(y)
  
}


# 4.Value-at-risk and TVaR measure for comGB2 distribution
VaR.comGB2 <- qcomGB2 # define the Value-at-risk measure

TVaR.comGB2 <- function(p, mu2, p1, p2 , nu1, nu2, tau1, tau2){
  g1 <- (p1*nu1 - 1)/(p1*tau1 + 1)
  g2 <- (p2*nu2 - 1)/(p2*tau2 + 1)
  x0 <- mu2*g2^(1/p2)
  mu1 <- mu2*g2^(1/p2)/g1^(1/p1)
  f_1 <- dGB2(x = x0, mu = mu1, sigma = p1, nu = nu1, tau = tau1)
  F_1 <- pGB2(q = x0, mu = mu1, sigma = p1, nu = nu1, tau = tau1)
  f_2 <- dGB2(x = x0, mu = mu2, sigma = p2, nu = nu2, tau = tau2)
  F_2 <- pGB2(q = x0, mu = mu2, sigma = p2, nu = nu2, tau = tau2)

  r <- f_2*F_1*(f_2*F_1+f_1*(1-F_2))^(-1)
 
  if(p <= r){
    u1sim <- p/r*F_1
    z <- qbeta(u1sim, shape1 = nu1, shape2 = tau1)
    out <- mu1*beta(nu1 + 1/p1, tau1 - 1/p1)/beta(nu1, tau1)*pbeta(q = z, shape1 = nu1 + 1/p1, shape2 = tau1 - 1/p1)/(1 - u1sim)

  } else if (p > r) {
    u2sim <- F_2 + (p - r)/(1 - r)*(1 - F_2)
    z <- qbeta(u2sim, shape1 = nu2, shape2 = tau2)
    out <- mu2*beta(nu2 + 1/p2, tau2 - 1/p2)/beta(nu2, tau2)*(1 - pbeta(q = z, shape1 = nu2 + 1/p2, shape2 = tau2 - 1/p2))/(1 - u2sim)

  }
  return(out)


}
TVaR.comGB2 <- Vectorize(TVaR.comGB2)



# 5. first moment for comGB2 distribution (expectation)
mean.comGB2 <- function(p, mu2, p1, p2, nu1, nu2, tau1, tau2){
  g1 <- (p1*nu1 - 1)/(p1*tau1 + 1)
  g2 <- (p2*nu2 - 1)/(p2*tau2 + 1)
  x0 <- mu2*g2^(1/p2)
  mu1 <- mu2*g2^(1/p2)/g1^(1/p1)
  f_1 <- dGB2(x = x0, mu = mu1, sigma = p1, nu = nu1, tau = tau1)
  F_1 <- pGB2(q = x0, mu = mu1, sigma = p1, nu = nu1, tau = tau1)
  f_2 <- dGB2(x = x0, mu = mu2, sigma = p2, nu = nu2, tau = tau2)
  F_2 <- pGB2(q = x0, mu = mu2, sigma = p2, nu = nu2, tau = tau2)
  
  r <- f_2*F_1*(f_2*F_1+f_1*(1-F_2))^(-1)
  
  m1 <- (p1*nu1 - 1)/(p1*nu1 + p1*tau1)
  m2 <- (p2*nu2 - 1)/(p2*nu2 + p2*tau2)
  
  part1 <- mu1*beta(nu1 + 1/p1, tau1 - 1/p1)/beta(nu1, tau1)*pbeta(m1, shape1 = nu1 + 1/p1, shape2 = tau1 - 1/p1)/pbeta(m1, shape1 = nu1 , shape2 = tau1)
  part2 <- mu2*beta(nu2 + 1/p2, tau2 - 1/p2)/beta(nu2, tau2)*(1 - pbeta(m2, shape1 = nu2 + 1/p2, shape2 = tau2 - 1/p2))/(1 - pbeta(m2, shape1 = nu2 , shape2 = tau2))
  
  out <- r*part1 + (1 - r)*part2
  
  return(out)
}
mean.comGB2 <- Vectorize(mean.comGB2)

mean.comGB2(mu2 = 1, p1 = c(2), p2 = 2, tau1 = c(2), tau2 = 1, nu1 = c(2), nu2 = 3)



