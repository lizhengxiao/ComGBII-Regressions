library(gamlss)
library(data.table)
library(SMPracticals)

data("danish") # read the data
data <- data.table(danish)
y <- as.numeric(data$danish)
summary(y); mean(y); sd(y)
X = model.matrix( ~ 1, data = data.frame(y = y)) # design matrix

source('1-comGB2-model-solnp.r', local = TRUE)
source('2-GBIIG-model-solnp.r', local = TRUE)
source('3-BIIG-model-solnp.r', local = TRUE)
source('4-BG-model-solnp.r', local = TRUE)
source('5-IBG-model-solnp.r', local = TRUE)
source('6-PG-model-solnp.r', local = TRUE)
source('7-IPG-model-solnp.r', local = TRUE)
source('functions-fit-comGBIIReg.r', local = TRUE)
source('functions-comGB2-distribution.r', local = TRUE)


# =========================================================
# model fitting
# =========================================================
mycontrol <- list(pars.seed = 3223, # Random number seeds
                  pars.init.n = 5, 
                  rho = 0.1, 
                  tol.rel = 1.0e-8)

m.comGB2 <- comGB2Reg.solnp(y = y, X = X, control = mycontrol)
summary.comGB2(model = m.comGB2, class = 'transformed') # shows the estimates for beta_0, p_1, p_2, tau_1, tau_2, nu_1 and nu_2

m.GBIIG <- GBIIGReg.solnp(y = y, X = X, control = mycontrol)
summary.GBIIG(model = m.GBIIG, class = 'transformed')

m.BIIG <- BIIGReg.solnp(y = y, X = X, control = mycontrol)
summary.BIIG(model = m.BIIG, class = 'original')
summary.BIIG(model = m.BIIG, class = 'transformed')

m.BG <- BGReg.solnp(y = y, X = X,  control = mycontrol)
summary.BG(model = m.BG, class = 'transformed')

m.IBG <- IBGReg.solnp(y = y, X = X, control = mycontrol)
summary.IBG(model = m.IBG, class = 'transformed')

m.PG <- PGReg.solnp(y = y, X = X,  control = mycontrol)
summary.PG(model = m.PG, class = 'transformed')

m.IPG <- IPGReg.solnp(y = y, X = X, control = mycontrol)
summary.IPG(model = m.IPG, class = 'transformed')




# model comparison (log likelihood, AIC and BIC sorting) -----------------------


family.name <- c('comGB2', 'GBIIG',
                 'BIIG', 'BG', 'IBG', 'PG', 'IPG')
m.name <- paste0('m.', family.name)
method.name <- paste0('summary.', family.name)
results.mat  <- se.mat <- matrix(0, nrow = 7, ncol = 4 + 7)
rownames(results.mat) <- rownames(se.mat) <- m.name
colnames(results.mat) <- colnames(se.mat) <- c('log(mu2)',
                           'p1', 'p2', 'tau1', 'tau2', 
                           'nu1', 'nu2',
                           'npars','NLL','AIC','BIC')
for(i in 1:7){
  m.temp <- get(m.name[i])
  est.table <- do.call(what = get(method.name[i]), 
                   args = list(model = m.temp, class = 'transformed'))
  m.est <- est.table$summarytable[,1]
  m.se <- est.table$summarytable[,2]
  results.mat[i,] <- c(m.est,
                       m.temp$npars, 
                       m.temp$NLL,  m.temp$AIC, m.temp$BIC)
  se.mat[i,] <- c(m.se,
                  m.temp$npars, 
                  m.temp$NLL,  m.temp$AIC, m.temp$BIC)
  
}
round(results.mat[order(results.mat[,2+7]), ], 3) # NULL sorting

round(results.mat[order(results.mat[,3+7]), ], 3) # AIC sorting
results.mat[order(results.mat[,4+7]), ] # BIC sorting

round(results.mat, 2)


