library(gamlss)
library(GB2)

pars.initialization <- function(y, X){
  # mod <- gamlss(y ~ X - 1, family = GB2(), method = mixed(2, 100000), 
  #               control = gamlss.control(trace = FALSE))
  # beta.init <- coefficients(mod, what = 'mu')
  # alpha.init <- c(coefficients(mod, what = 'sigma'),
  #                 coefficients(mod, what = 'tau'),
  #                 coefficients(mod, what = 'nu'))
  # names(alpha.init) <- c('p', 'tau', 'nu')
  # out <- list(beta.init = beta.init, shape.init = alpha.init)
  # 尝试运行这个程序段，捕捉错误并执行相应的操作
  # tryCatch({
  #   mod <- glm(y ~ X-1, family = Gamma(link = 'log'))
  #   beta.init <- coef(mod)
  # }, error = function(e) {
  #   print("An error occurred, running alternative code.")
  #   mod <- gamlss(y ~ X - 1, family = GA(), method = mixed(2, 100000),
  #                 control = gamlss.control(trace = FALSE))
  #   beta.init <- coefficients(mod, what = 'mu')
  # })
  mod <- gamlss(y ~ X - 1, family = GA(), method = mixed(2, 100000),
                control = gamlss.control(trace = FALSE))
  beta.init <- coefficients(mod, what = 'mu')
  # mod <- glm(y ~ X-1, family = Gamma(link = 'log'))
  # beta.init <- coef(mod)
  out <- list(beta.init = beta.init)
  return(out)
}

# composite GB2 distribution --------------------------------------
# density function
dcomGB2 <- function(y, mu2, p1, p2 , nu1, nu2, tau1, tau2, log = FALSE){
  g1 <- (p1*nu1 - 1)/(p1*tau1 + 1)
  g2 <- (p2*nu2 - 1)/(p2*tau2 + 1)
  x0 <- mu2*g2^(1/p2)
  mu1 <- mu2*g2^(1/p2)/g1^(1/p1)
  f_1 <- dGB2(x = x0, mu = mu1, sigma = p1, nu = nu1, tau = tau1)
  F_1 <- pGB2(q = x0, mu = mu1, sigma = p1, nu = nu1, tau = tau1)
  f_2 <- dGB2(x = x0, mu = mu2, sigma = p2, nu = nu2, tau = tau2)
  F_2 <- pGB2(q = x0, mu = mu2, sigma = p2, nu = nu2, tau = tau2)
  
  r <- f_2*F_1*(f_2*F_1+f_1*(1-F_2))^(-1)
  
  index1 <- (y <= x0)
  index2 <- y>x0
  fpart1 <- r*dGB2(x = y, mu = mu1, sigma = p1, nu = nu1, tau = tau1)/F_1
  fpart2 <- (1-r)*dGB2(x = y, mu = mu2, sigma = p2, nu = nu2, tau = tau2)/(1-F_2)
  
  # y <= x0
  out <- fpart1
  # y > x0
  out[index2] <- fpart2[index2]
  
  if(log == TRUE){ll = log(out)} else {ll = out}
  return(ll)
}

# using gb2 package
dcomGB2_v3 <- function(y, mu2, p1, p2 , nu1, nu2, tau1, tau2, log = FALSE){
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
  
  # y <= x0
  out <- fpart1
  # y > x0
  out[index2] <- fpart2[index2]
  
  if(log == TRUE){ll = log(out)} else {ll = out}
  return(ll)
}


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



VaR.comGB2 <- qcomGB2

TVaR2.comGB2 <- function(p, mu2, p1, p2 , nu1, nu2, tau1, tau2){
  g1 <- (p1*nu1 - 1)/(p1*tau1 + 1)
  g2 <- (p2*nu2 - 1)/(p2*tau2 + 1)
  x0 <- mu2*g2^(1/p2)
  mu1 <- mu2*g2^(1/p2)/g1^(1/p1)
  f_1 <- dGB2(x = x0, mu = mu1, sigma = p1, nu = nu1, tau = tau1)
  F_1 <- pGB2(q = x0, mu = mu1, sigma = p1, nu = nu1, tau = tau1)
  f_2 <- dGB2(x = x0, mu = mu2, sigma = p2, nu = nu2, tau = tau2)
  F_2 <- pGB2(q = x0, mu = mu2, sigma = p2, nu = nu2, tau = tau2)

  r <- f_2*F_1*(f_2*F_1+f_1*(1-F_2))^(-1)
  # VaR <- qcomGB2(p = p, mu2 = mu2, p1 = p1, p2 = p2 , nu1 = nu1, nu2 = nu2, tau1 = tau1, tau2 = tau2)


  # if(p <= r){
  #   u1sim <- p/r*F_1
  #   Isim <- qbeta(p = u1sim, shape1 = nu1, shape2 = tau1)
  #   zsim <- Isim/(1 - Isim)
  #   VaR <- mu1*zsim^(1/p1)
  #
  #   # out <- qGB2(p = u1sim, mu = mu1, sigma = p1, nu = nu1, tau = tau1)
  #
  # } else if (p > r) {
  #   u2sim <- F_2 + (p - r)/(1 - r)*(1 - F_2)
  #   VaR <- qGB2(p = u2sim, mu = mu2, sigma = p2, nu = nu2, tau = tau2)
  # }

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

TVaR2.comGB2 <- Vectorize(TVaR2.comGB2)




TVaR.comGB2 <- function(p, mu2, p1, p2, nu1, nu2, tau1, tau2){
  myfun <- function(p){
    VaR.comGB2(p = p,
               mu2 = mu2, p1 = p1, p2 = p2, tau1 = tau1, tau2 = tau2, nu1 = nu1, nu2 = nu2)

  }
  out <- integrate(myfun, lower = p, upper = 1)$value/(1 - p)

  return(out)


}
TVaR.comGB2 <- Vectorize(TVaR.comGB2)
# TVaR.comGB2.new <- Vectorize(TVaR.comGB2.new)


dcomGB2(c(1),
        mu2 = 1, p1 = c(2), p2 = 2, tau1 = c(2), tau2 = 1, nu1 = c(2), nu2 = 3)


TVaR.comGB2(p = 0.9,
        mu2 = 1, p1 = c(2), p2 = 2, tau1 = c(2), tau2 = 1, nu1 = c(2), nu2 = 3)
# VaR.comGB2(p = 0.9,
#            mu2 = 1, p1 = c(2), p2 = 2, tau1 = c(2), tau2 = 1, nu1 = c(2), nu2 = 3)

TVaR2.comGB2(p = 0.9,
           mu2 = 1, p1 = c(2), p2 = 2, tau1 = c(2), tau2 = 1, nu1 = c(2), nu2 = 3)





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

mean2.comGB2 <- function(mu2, p1, p2, nu1, nu2, tau1, tau2){
  myfun <- function(y){
    y*dcomGB2(y = y,
               mu2 = mu2, p1 = p1, p2 = p2, tau1 = tau1, tau2 = tau2, nu1 = nu1, nu2 = nu2)
    
  }
  out <- integrate(myfun, lower = 0, upper = Inf)$value
  
  return(out)
  
}

mean.comGB2(mu2 = 1, p1 = c(2), p2 = 2, tau1 = c(2), tau2 = 1, nu1 = c(2), nu2 = 3)
mean2.comGB2(mu2 = 1, p1 = c(2), p2 = 2, tau1 = c(2), tau2 = 1, nu1 = c(2), nu2 = 3)





mean.GB2 <- function(p, mu, nu, tau){
  
  part <- mu*beta(nu + 1/p, tau - 1/p)/beta(nu, tau)
  
  out <- part
  
  return(out)
  
}

# TVaR.comGB2.new(p = c(0.9, 0.8),
            # mu2 = 1, p1 = c(2), p2 = 2, tau1 = c(2), tau2 = 1, nu1 = c(2), nu2 = 3)


# pars.int <- c(-1, -1, 0.5, -0.2, -0.2, 2, 5)
# pars <- pars.int
# mu <- pars[1]
# p <- pars[2:3]
# tau <- pars[4:5]
# nu <- pars[6:7]
# mu2 <- exp(mu)
# p1 <- exp(p[1]); p2 <- exp(p[2])
# tau1 <- exp(tau[1]); tau2 <- exp(tau[2])
# nu1 <- exp(nu[1]); nu2 <- exp(nu[2])
# 
# g1 <- (p1*nu1 - 1)/(p1*tau1 + 1)
# g2 <- (p2*nu2 - 1)/(p2*tau2 + 1)
# x0 <- mu2*g2^(1/p2)
# mu1 <- mu2*(g2^((1/p2)))/g1^(1/p1)
# 
# f_1 <- dGB2(x = x0, mu = mu1, sigma = p1, nu = nu1, tau = tau1)
# F_1 <- pGB2(q = x0, mu = mu1, sigma = p1, nu = nu1, tau = tau1)
# f_2 <- dGB2(x = x0, mu = mu2, sigma = p2, nu = nu2, tau = tau2)
# F_2 <- pGB2(q = x0, mu = mu2, sigma = p2, nu = nu2, tau = tau2)
# 
# r <- f_2*F_1/(f_2*F_1+f_1*(1-F_2))
# 
# LLcomGB2(y = y, pars = pars.int)
# dcomGB2(y = y,
#         mu2 = mu2,
#         p1 = p1,
#         p2 = p2, tau1 = tau1, tau2 = tau2, nu1 = nu1, nu2 = nu2,
#         log = F)
# pcomGB2(y = y,
#         mu2 = mu2,
#         p1 = p1,
#         p2 = p2, tau1 = tau1, tau2 = tau2, nu1 = nu1, nu2 = nu2)
# qcomGB2(p = seq(from = 0.01, to = 0.9, length.out = 20),
#         mu2 = mu2,
#         p1 = p1,
#         p2 = p2, tau1 = tau1, tau2 = tau2, nu1 = nu1, nu2 = nu2)

