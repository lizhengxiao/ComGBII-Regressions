comGB2RegM <- function(y, X, family = c('comGB2'),
                       control = list(print_level = 1,
                                      pars.seed = 113,
                                      pars.init.n = 5, 
                                      tol.rel = 1.0e-8,
                                      rho = 0.1)){

  if (family == 'comGB2'){
    mod <- comGB2Reg.solnp(y = y, X = X, 
                     control = control)
  } else if(family == 'GBIIG'){
    mod <- GBIIGReg.solnp(y = y, X = X, 
                     control = control)
  }  else if(family == 'BIIG'){
    mod <- BIIGReg.solnp(y = y, X = X, 
                    control = control)
  }  else if(family == 'BG'){
    mod <- BGReg.solnp(y = y, X = X, 
                    control = control)
  }  else if(family == 'IBG'){
    mod <- IBGReg.solnp(y = y, X = X, 
                    control = control)
  }  else if(family == 'PG'){
    mod <- PGReg.solnp(y = y, X = X, 
                    control = control)
  } else if(family == 'IPG'){
    mod <- IPGReg.solnp(y = y, X = X, 
                 control = control)
  } 
  
  return(mod)
  
  
}