## case I -----------------------------------------------------------
library(gamlss)
library(Rsolnp)
set.seed(111)
nsample <- 2000
Nsim <- 10
beta.mat <- matrix(0, nrow = Nsim, ncol = 3)
alpha.mat <- matrix(0, nrow = Nsim, ncol = 6)
p1 <- 1; p2 <- 1.5; 
tau1 <- 1.5; tau2 <- 2; 
nu1 <- 2; nu2 <- 2;
sigma1 <- p1*tau1; sigma2 <- p2*tau2
delta1 <- p1*nu1; delta2 <- p2*nu2

pars1.true <- c(c(2, 0.5, 1), 
               log(c(p1, p2, tau1, tau2, nu1, nu2))) # true value
pars2.true <- c(c(2, 0.5, 1), 
                log(c(p1, p2, sigma1, sigma2, delta1, delta2))) # true value

for (j in 1:Nsim){
  x1 <- rnorm(n = nsample, mean = 0, sd = 1)
  x2 <- rnorm(n = nsample, mean = 0, sd = 1)
  X1 <- model.matrix( ~ x1 + x2)
  k <- dim(X1)[2]
  pars <- pars1.true
  beta <- pars[1:k]
  mu2 <- exp(X1%*%beta)
  p1 <- exp(pars[(k+1)]); p2 <- exp(pars[(k+2)])
  tau1 <- exp(pars[(k+3)]); tau2 <- exp(pars[(k+4)])
  nu1 <- exp(pars[(k+5)]); nu2 <- exp(pars[(k+6)])

  p1*nu1 - 1
  p2*nu2 - 1
  ysim <- 0
  for(i in 1:nsample){
    ysim[i] <-  qcomGB2(p = runif(1, 0, 1), mu2 = as.numeric(mu2)[i], 
                         p1 = p1, p2 = p2, nu1 = nu1, nu2 = nu2, tau1 = tau1, tau2 = tau2)
  }
  y <- ysim
  X <- X1


  tryCatch({
    mReg <- comGB2Reg.solnp(y = y, X = X, control = list(pars.init.n = 5, 
                                                         rho = 2, tol.rel = 1.0e-4))
  }, error = function(e) {
    message('Return the next iteration')
  })
  
  estimates <- mReg$estimate
  beta.mat[j,] <- estimates[1:k]
  alpha.mat[j,] <- estimates[(k+1):(k+6)]
}

pars.true <- pars1.true
beta.true <- pars.true[1:3]
alpha.true <- pars.true[4:9]


# ggplot
library(data.table)
dtbeta <- data.frame(beta.mat)
colnames(dtbeta) <- c('beta0', 'beta1', 'beta2')
dtplot <- data.table()
dtplot$y <- c(beta.mat[,1], beta.mat[,2], beta.mat[,3],
              alpha.mat[,1], alpha.mat[,2],alpha.mat[,3], alpha.mat[,4],alpha.mat[,5], alpha.mat[,6])
dtplot$Parameters <- c(rep('beta0', Nsim),
                       rep('beta1', Nsim),
                       rep('beta2', Nsim),
                       rep('alpha0', Nsim),
                       rep('alpha1', Nsim),
                       rep('alpha2', Nsim),
                       rep('alpha3', Nsim),
                       rep('alpha4', Nsim),
                       rep('alpha5', Nsim))
dtplot$group <- as.factor(dtplot$Parameters)
labels = c(expression(paste(beta[0])),
           expression(paste(beta[1])),
           expression(paste(beta[2])),
           expression(paste(alpha[0])),
           expression(paste(alpha[1])),
           expression(paste(alpha[2])),
           expression(paste(alpha[3])),
           expression(paste(alpha[4])),
           expression(paste(alpha[5]))
           )
dtbeta <- dtplot[Parameters %in% c('beta0','beta1', 'beta2')]
dtalpha <- dtplot[!Parameters %in% c('beta0','beta1', 'beta2')]

library(ggplot2)
dtbeta <- data.frame(dtbeta)
fig1basic <- ggplot(data = dtbeta, mapping = aes(y = y, x = Parameters, fill = Parameters)) + 
  theme_bw() + xlab('') + ylab('Estimates') 
fig1 <- fig1basic +  geom_boxplot(outlier.colour = "red", outlier.shape = 1,
                               outlier.size = 1, varwidth = TRUE) + 
        scale_y_continuous(limits = c(0, 3),
                           breaks = seq(0, 3, by = 0.5),
                           expand = c(0, 0)) +
        scale_x_discrete(labels = c(expression(paste(beta[0])),
                                    expression(paste(beta[1])),
                                    expression(paste(beta[2]))
                                    )) +
  scale_fill_discrete(name = "Parameters",
                      labels = c(expression(paste(beta[0])),
                                 expression(paste(beta[1])),
                                 expression(paste(beta[2]))
                      ))



dtalpha <- data.frame(dtalpha)
fig2basic <- ggplot(data = dtalpha, mapping = aes(y = y, x = Parameters, fill = Parameters)) + 
  theme_bw() + xlab('') + ylab('Estimates') 
fig2 <- fig2basic +  geom_boxplot(outlier.colour = "red", outlier.shape = 1,
                                  outlier.size = 1, varwidth = TRUE) + 
  scale_y_continuous(limits = c(-8, 8),
                     breaks = seq(-8, 8, by = 2),
                     expand = c(0, 0)) +
  scale_x_discrete(labels = c(expression(paste(alpha[1])),
                              expression(paste(alpha[2])),
                              expression(paste(alpha[3])),
                              expression(paste(alpha[4])),
                              expression(paste(alpha[5])),
                              expression(paste(alpha[6]))
  )) +
  scale_fill_discrete(name = "Parameters",
                      labels = c(expression(paste(alpha[1])),
                                 expression(paste(alpha[2])),
                                 expression(paste(alpha[3])),
                                 expression(paste(alpha[4])),
                                 expression(paste(alpha[5])),
                                 expression(paste(alpha[6]))
                      ))



library(patchwork)
fig0 <- fig1 + fig2 + plot_layout(ncol = 2)
fig0






