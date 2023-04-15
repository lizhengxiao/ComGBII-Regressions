library(gamlss) # load the gamlss package
library(Rsolnp) # load the Rsolnp package
set.seed(111) # set the seed for reproducibility
nsample <- 2000 # set the sample size
Nsim <- 20 # set the number of simulations
beta.mat <- matrix(0, nrow = Nsim, ncol = 3) # initialize the matrix for storing beta estimates
alpha.mat <- matrix(0, nrow = Nsim, ncol = 6) # initialize the matrix for storing alpha estimates
p1 <- 1; p2 <- 1.5; # set values for p1 and p2
tau1 <- 1.5; tau2 <- 2; # set values for tau1 and tau2
nu1 <- 2; nu2 <- 2; # set values for nu1 and nu2
sigma1 <- p1*tau1; sigma2 <- p2*tau2 # calculate sigma1 and sigma2
delta1 <- p1*nu1; delta2 <- p2*nu2 # calculate delta1 and delta2
pars1.true <- c(c(2, 0.5, 1), 
                log(c(p1, p2, tau1, tau2, nu1, nu2))) # true values of parameters


# Loop over simulations
for (j in 1:Nsim){
  # Generate normally distributed random variables
  x1 <- rnorm(n = nsample, mean = 0, sd = 1) 
  x2 <- rnorm(n = nsample, mean = 0, sd = 1) 
  
  # Create a model matrix with x1 and x2
  X1 <- model.matrix( ~ x1 + x2) 
  
  # Get the number of columns in X1
  k <- dim(X1)[2] 
  
  # Initialize the parameters with true values
  pars <- pars1.true 
  
  # Get beta parameters from pars
  beta <- pars[1:k] 
  
  # Calculate the expected value of mu2
  mu2 <- exp(X1%*%beta) 
  
  # Calculate p1 and p2
  p1 <- exp(pars[(k+1)]); p2 <- exp(pars[(k+2)]) 
  
  # Calculate tau1 and tau2
  tau1 <- exp(pars[(k+3)]); tau2 <- exp(pars[(k+4)]) 
  
  # Calculate nu1 and nu2
  nu1 <- exp(pars[(k+5)]); nu2 <- exp(pars[(k+6)]) 
  
  # Check if p1 and nu1 are greater than 1
  p1*nu1 - 1 
  
  # Check if p2 and nu2 are greater than 1
  p2*nu2 - 1 
  
  # Initialize the simulated response variable
  ysim <- 0 
  
  # Loop over the sample size
  for(i in 1:nsample){ 
    # Simulate response variable using the qcomGB2 function
    ysim[i] <-  qcomGB2(p = runif(1, 0, 1), mu2 = as.numeric(mu2)[i], 
                        p1 = p1, p2 = p2, nu1 = nu1, nu2 = nu2, tau1 = tau1, tau2 = tau2) 
  }
  
  # Set y equal to the simulated response variable
  y <- ysim 
  
  # Set X equal to the model matrix
  X <- X1 
  
  # Use comGB2Reg.solnp to estimate the parameters
  tryCatch({
    mReg <- comGB2Reg.solnp(y = y, X = X, control = list(pars.init.n = 5, 
                                                         rho = 2, tol.rel = 1.0e-4))
  }, error = function(e) {
    message('Return the next loop ')
  })
  
  # Get the estimates from mReg
  estimates <- mReg$estimate
  
  # Store beta estimates in the beta.mat matrix
  beta.mat[j,] <- estimates[1:k]
  
  # Store alpha estimates in the alpha.mat matrix
  alpha.mat[j,] <- estimates[(k+1):(k+6)]
}

pars.true <- pars1.true
beta.true <- pars.true[1:3]
alpha.true <- pars.true[4:9]


# ggplot -------------------------------------------------------------------------------------------------

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






