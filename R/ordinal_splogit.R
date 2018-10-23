#'
#' Data setup for fitting
#'

setup.ordinal.nonspatial <- function(formula, data, na.action){
  mf <- match.call(expand.dots = FALSE)
  m <- match(c('formula', 'data', 'na.action'), names(mf), 0)
  mf <- mf[c(1, m)]
  mf <- eval(mf[[3]], parent.frame())
  f <- Formula(formula)
  
  data <- list(
    n = nrow(data),
    y = as.numeric(model.frame(formula(terms(f, lhs = 1, rhs = 0)), data = mf)[,1])
  )
  
  # covariates for direct effects on outcome
  data$X <- model.matrix(f, data = mf, rhs = 1) # fixed effects
  data$N_X <- dim(data$X)[2]
  data$N_obs <- data$n
  data$K <- length(unique(data$y))
  list(formula = formula, data = data)
}

setup.ordinal.timeseries <- function(formula, data, na.action){
  mf <- match.call(expand.dots = FALSE)
  m <- match(c('formula', 'data', 'na.action'), names(mf), 0)
  mf <- mf[c(1, m)]
  mf <- eval(mf[[3]], parent.frame())
  f <- Formula(formula)
  
  # location for spatial effects
  t <- as.numeric(model.frame(formula(terms(f, lhs = 0, rhs = 2)), data = mf)[,1])
  
  data <- list(
    n = length(t),
    t = t,
    y = as.numeric(model.frame(formula(terms(f, lhs = 1, rhs = 0)), data = mf)[,1])
  )
  
  data$distmat <- as.matrix(dist(t)^2)  # distance matrix
  
  # covariates for direct effects on outcome
  data$X <- model.matrix(f, data = mf, rhs = 1) # fixed effects
  data$N_X <- dim(data$X)[2]
  data$N_obs <- data$n
  data$K <- length(unique(data$y))
  list(formula = formula, data = data)
}

#'
#' Simulate data
#' 

simulate.ordinal.timeseries <- function(ord.object, 
                                        beta = rep(1, N_X),
                                        eta = c(0, 1),
                                        a = 1,
                                        b = 1,
                                        r = 1,  # power parameter of the splogit link function
                                        jitter = 1e-6, 
                                        seed = NULL){
  
  if(!is.null(seed)) set.seed(seed)  
  for(nm in names(ord.object$data)){
    assign(nm, ord.object$data[[nm]])
  }
  
  if (b == 0){
    sigma <- a * diag(n) + jitter * diag(n)
  }else{
    sigma <- a * exp(- distmat / b) + jitter * diag(n)
  }
  
  # Gaussian Process, spatial effects
  f <- as.vector(rmvnorm(1, mean = rep(0, n), sigma = sigma))
  # latent variables
  z <- as.vector(X %*% beta) + f + rplogit(n, r)    
  y <- as.vector(do.call(cbind, lapply(1:n, function(i){
    match(z[i], sort(c(z[i], eta)))
  })))
  proportion <- as.data.frame(table(y)/n)
  list(y = y, proportion = proportion, z.latent = z, f = f)
}

###################################################

# random generation from power logit distribution #
rplogit <- function(n, r){
  u <- runif(n, 0, 1)
  if(r <= 1){
    t <- u^(1/r)
    return (r * log(t/(1-t)))
  }else if(r > 1){
    t <- (1-u)^r
    return ((-1/r) * log(t/(1-t)))
  }
}

# density function of power logit distribution #
dplogit <- function(x, location = 0, r){
  if(r <= 1){
    return(plogis((x-location)/r, 0, 1)^(r-1) * dlogis((x-location)/r, 0, 1))
  }else if(r > 1){
    return(plogis(-r*(x-location), 0, 1)^(1/r-1) * dlogis(-r*(x-location), 0, 1))
  }
}

# distribution function of power logit distribution #
pplogit <- function(x, location = 0, r){
  if(r <= 1){
    return(plogis((x-location)/r, 0, 1)^r)
  }else if (r > 1){
    return(1 - plogis(-r*(x-location), 0, 1)^(1/r))
  }
}

#### Compute Model Comparison Criterion: WAIC and LOOIC, LPML, DIC ####

comparison.splogit <- function(fit, dd, dependent = c(FALSE, TRUE)[1]){
  K <- length(unique(dd$Y))
  ## WAIC and LOOIC ##
  ll <- extract_log_lik(fit)
  WAIC <- waic(ll)$waic
  LOOIC <- loo(ll)$looic
  
  ## LPML ##
  ilik <- exp(ll)
  inv_ilik <- 1./ilik
  CPOi <- 1./colMeans(inv_ilik)
  LPML <- sum(log(CPOi))
  
  ## compute DIC: DIC = 2*Dbar - Dhat = Dbar + pD
  # pD = Dbar - Dhat
  # Dbar is average deviance
  # Dhat is deviance at average parameter values
  # pD is effective number of parameters
  
  dev <- (-2) * rowSums(ll)
  Dbar <- mean(dev)
  # extract samples
  e <- extract(fit, inc_warmup = FALSE)
  beta.mean <- colMeans(e$beta)
  cut.mean <- colMeans(e$cuts)[-1]
  r.mean <- mean(e$r)
  parvalue.mean <- c(beta.mean = beta.mean, cut.mean = cut.mean, r.mean = r.mean)
 
  # probability of power logit
  probfn <- function(obj, cutpoints, r){
    if (obj[1] == 1){
      prob <- pplogit(x = -obj[2], r = r)
    }else{
      if (obj[1] == K){
        prob <- 1 - pplogit(x = cutpoints[K-1] - obj[2], r = r)
      }else{
        prob <- pplogit(x = cutpoints[obj[1]] - obj[2], r = r) - pplogit(x = cutpoints[obj[1]-1] - obj[2], r = r)
      }
    }
    return(prob)
  }
  
  if (dependent == FALSE){
    loglik <- function(y, x, parvalue){
      X <- cbind(1, x)
      beta <- parvalue.mean[grepl("beta", names(parvalue.mean))]
      cutpoints <- c(0, parvalue.mean[grepl("cut", names(parvalue.mean))])
      r <- parvalue.mean[grepl("r", names(parvalue.mean))]
      temp <- cbind(y, X %*% beta)
      log_lik <- sum(log(apply(temp, 1, probfn, cutpoints = cutpoints, r = r))) * (-2)
      return(log_lik)
    }
    Dhat <- loglik(dd$Y, dd$x1, parvalue.mean)
    pD <- Dbar - Dhat
    DIC <- Dbar + pD
  }
  
  if (dependent == TRUE){
    f.mean <- colMeans(e$f)
    loglik <- function(y, x, parvalue, f){
      X <- cbind(1, x)
      beta <- parvalue.mean[grepl("beta", names(parvalue.mean))]
      cutpoints <- c(0, parvalue.mean[grepl("cut", names(parvalue.mean))])
      r <- parvalue.mean[grepl("r", names(parvalue.mean))]
      temp <- cbind(y, X %*% beta + f)
      log_lik <- sum(log(apply(temp, 1, probfn, cutpoints = cutpoints, r = r))) * (-2)
      return(log_lik)
    }
    Dhat <- loglik(dd$Y, dd$x1, parvalue.mean, f.mean)
    pD <- Dbar - Dhat
    DIC <- Dbar + pD
  }
  
  return(list(WAIC = WAIC, LOOIC = LOOIC, LPML = LPML, DIC = DIC))
}

comparison.logit <- function(fit, dd, dependent = c(FALSE, TRUE)[1]){
  K <- length(unique(dd$Y))
  ## WAIC and LOOIC ##
  ll <- extract_log_lik(fit)
  WAIC <- waic(ll)$waic
  LOOIC <- loo(ll)$looic
  
  ## LPML ##
  ilik <- exp(ll)
  inv_ilik <- 1./ilik
  CPOi <- 1./colMeans(inv_ilik)
  LPML <- sum(log(CPOi))
  
  ## compute DIC: DIC = 2*Dbar - Dhat = Dbar + pD
  # pD = Dbar - Dhat
  # Dbar is average deviance
  # Dhat is deviance at average parameter values
  # pD is effective number of parameters
  
  dev <- (-2) * rowSums(ll)
  Dbar <- mean(dev)
  # extract samples
  e <- extract(fit, inc_warmup = FALSE)
  beta.mean <- colMeans(e$beta)
  cut.mean <- colMeans(e$cuts)[-1]
  parvalue.mean <- c(beta.mean = beta.mean, cut.mean = cut.mean)
  
  # probability of power logit
  probfn <- function(obj, cutpoints){
    if (obj[1] == 1){
      prob <- plogis(-obj[2], 0, 1)
    }else{
      if (obj[1] == K){
        prob <- 1 - plogis(cutpoints[K-1] - obj[2], 0, 1)
      }else{
        prob <- plogis(cutpoints[obj[1]] - obj[2], 0, 1) - plogis(cutpoints[obj[1]-1] - obj[2], 0, 1)
      }
    }
    return(prob)
  }
  
  if (dependent == FALSE){
    loglik <- function(y, x, parvalue){
      X <- cbind(1, x)
      beta <- parvalue.mean[grepl("beta", names(parvalue.mean))]
      cutpoints <- c(0, parvalue.mean[grepl("cut", names(parvalue.mean))])
      temp <- cbind(y, X %*% beta)
      log_lik <- sum(log(apply(temp, 1, probfn, cutpoints = cutpoints))) * (-2)
      return(log_lik)
    }
    Dhat <- loglik(dd$Y, dd$x1, parvalue.mean)
    pD <- Dbar - Dhat
    DIC <- Dbar + pD
  }
  
  if (dependent == TRUE){
    f.mean <- colMeans(e$f)
    loglik <- function(y, x, parvalue, f){
      X <- cbind(1, x)
      beta <- parvalue.mean[grepl("beta", names(parvalue.mean))]
      cutpoints <- c(0, parvalue.mean[grepl("cut", names(parvalue.mean))])
      r <- parvalue.mean[grepl("r", names(parvalue.mean))]
      temp <- cbind(y, X %*% beta + f)
      log_lik <- sum(log(apply(temp, 1, probfn, cutpoints = cutpoints))) * (-2)
      return(log_lik)
    }
    Dhat <- loglik(dd$Y, dd$x1, parvalue.mean, f.mean)
    pD <- Dbar - Dhat
    DIC <- Dbar + pD
  }
  return(list(WAIC = WAIC, LOOIC = LOOIC, LPML = LPML, DIC = DIC))
}
