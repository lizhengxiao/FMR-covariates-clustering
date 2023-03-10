
# This function provides the ADMM algorithm for Gamma regression, given the auxiliary variable used in ADMM.
innerLayer_ADMM <- function(Y, X, 
                        beta0, betaX, param,
                        z.mat, r.mat,
                        lambda, rho,
                        weights, w = 1, 
                        control = list(tol = 1e-4,
                                       reltol = 1e-4,
                                       gradtol = 1e-4,
                                       iterlim = 10000,
                                       printLevel = 0)){
  # Y: the response variable.
  # X: the design matrix with dims n*(p+1), where p denotes the number of covariates, n is the sample size.
  # beta0, betaX, param: the initial parameters:
    # beta: a matrix with dims 1*H where H is the number of component for mixture distribution
    # betaX: a matrix with dims (p+1)*H
    # param：a list object containing the dispersion parameter in Gamma distribution with dim of 2: (mu, sigma).
  # lambda: the penalty parameter for the L2-norm of betaX.
  # rho: the tuning parameter for ADMM.
  # control: a list of control parameters for ADMM.
  # weights: an vector of ‘prior weights’ to be used in the fitting process.
  # w: the component weights in mixture models.
  # z.mat/r.mat: auxiliary variable used in ADMM, a matrix with dims p*p.
  
  ## parameter initialization
  pX <- ncol(X) - 1 # number of covarites
  n <- length(Y)    # sample size
  y <- Y # response
  int.pars <- c(beta0, betaX, log(param)) # initialization for all parameters
  int.beta <- c(beta0, betaX) # initialization for regression coefficients
  int.param <- log(param) # the log of dispersion parameter in Gamma distribution
  
  # Definition penalty term: L2 normal function
  Regularfun1 <- function(lambda = 1, pars){
    if (lambda != 0){
      L2 <- pars^2
      L2dev <- 2*pars
      Nl <- sum(lambda*L2)     # the penalty function
      Dl <- lambda*L2dev  # the first dev
      myfun <- function(pars){
        out <- lambda*pars^2
        out
      }
      D2l <- numDeriv::jacobian(myfun, x = pars)

    } else {
      Nl = 0
      Dl = rep(0, times = length(pars))
      D2l = matrix(0, nrow = length(pars), ncol = length(pars))
    }

    out <- list(value = Nl,
                deriv = Dl, # the first derivative
                deriv2 = D2l # the second derivative
                )
    out
  }

  # Definition penalty term: covarite clustering function
  Regularfun2 <- function(pars, z.mat, r.mat, rho){
    beta.mat <- matrix(rep(pars, times = pX), byrow = FALSE, nrow = pX)
    Nl <- sum(rho/2*(z.mat + r.mat - beta.mat)^2)
    Dl <- -rho*  (apply(z.mat + r.mat, 1, sum) - pX*pars)
    D2l <- matrix(0, nrow = pX, ncol = pX)
    diag(D2l) <- rho*pX

    out <- list(value = Nl,
                deriv = Dl,
                deriv2 = D2l)
    out
  }
 
  # Definition of the log-likelihood function for Gamma regression model
  fn.beta <- function(modpars, sigma, weights, z.mat, r.mat, lambda, rho, pX, y, X, w) {
    betaall <- modpars[1:(pX+1)]
    beta0 <- betaall[1]
    betaX <- betaall[2:(pX+1)]
    mu <- as.vector(exp(X%*%betaall)) # mu parameter
    # log-likelihood for the Gamma given that the value of sigma
    loglike <- -1*(1/sigma^2)*(2*log(sigma) + log(mu)) + (1/sigma^2 - 1)*log(y) - y/(sigma^2*mu) - lgamma(1/sigma^2)
    
    Reugular.p1 <- Regularfun1(lambda = lambda, pars = betaX)$value
    Reugular.p2 <-  Regularfun2(pars = betaX,
                                 z.mat = z.mat, r.mat = r.mat, rho = rho)$value

    ll <- sum(weights*loglike) - Reugular.p1*w - Reugular.p2

    return(ll)
  }
  # Definition of the first derivative of the log-likelihood function for Gamma regression model
  gn.beta <- function(modpars, weights, sigma, z.mat, r.mat, lambda, rho, pX, y, X, w){
    betaall <- modpars[1:(pX+1)]
    beta0 <- betaall[1]
    betaX <- betaall[2:(pX+1)]
    mu <- as.vector(exp(X%*%betaall)) # mu

    dldmu <- (y - mu)/((sigma^2) * (mu^2)) 
    A <- dldmu*mu*weights
    llderiv <- as.vector(t(A)%*%X)
    
    regderiv.p1 <- Regularfun1(lambda = lambda, pars = betaX)$deriv
    regderiv.p2 <- Regularfun2(pars = betaX,
                               z.mat = z.mat, r.mat = r.mat, rho = rho)$deriv
    regderivmat <- c(0, (regderiv.p1*w + regderiv.p2))

    return(llderiv - regderivmat)
  }
  # Definition of the Hessian matrix of the log-likelihood function for Gamma regression model
  hessn.beta <- function(modpars, sigma, weights, z.mat, r.mat, lambda, rho, pX, y, X, w){
    betaall <- modpars[1:(pX+1)]
    beta0 <- betaall[1]
    betaX <- betaall[2:(pX+1)]
    mu <- as.vector(exp(X%*%betaall)) # mu

    Ipars.mat <- matrix(0, nrow = pX + 1, ncol = pX + 1) # information matrix
    for(i in 1:length(y)){
      d2ldm2 = function(mu, sigma) -1/((sigma^2) * (mu^2)) # hessian for mu (mu)
      Ibb <- - d2ldm2(mu = mu[i], sigma = sigma)
      Ibapq.mat <- matrix(data = Ibb, nrow = 1, ncol = 1, byrow = T)
      xi <- as.vector(X[i,])
      D <- matrix(0, nrow = 1, ncol = pX + 1)
      D[1,1:(pX + 1)]  <- mu[i]*xi
      tempmat <- t(D)%*%Ibapq.mat%*%D*weights[i]
      Ipars.mat <- Ipars.mat + tempmat
    }

    Amat <- matrix(0, nrow = pX+1, ncol = pX+1)
    regderiv2.p1 <- Regularfun1(lambda = lambda, pars = betaX)$deriv2
    regderiv2.p2 <- Regularfun2(pars = betaX,
                                z.mat = z.mat, r.mat = r.mat, rho = rho)$deriv2
    Amat[2:(pX+1), 2:(pX+1)] <- regderiv2.p1*w + regderiv2.p2

    hessmat <- - Ipars.mat # hessian = - information matrix
    return(hessmat - Amat)
  }


  opt1 <- maxLik(start = int.beta,
                              logLik = fn.beta,
                              grad  = gn.beta,
                              hess = hessn.beta,
                              method = "BFGS",
                              sigma = param,
                              weights = weights, z.mat = z.mat, r.mat = r.mat,
                              lambda = lambda, rho = rho, y = Y, pX = pX,
                              X = X, w = w,
                              control = control)
  # update the regression coefficients
  out <- opt1$estimate
  beta0 <- out[1]
  betaX <- out[2:(1+pX)]
  
  # update the dispersion parameter
  # the log-likelihood given that the value of mu
  fn.param <- function(param.value, beta0, betaX, weights, y, X){
    sigma <- exp(param.value)
    mu <- exp(X%*%c(beta0, betaX))
 
    loglike <- -1*(1/sigma^2)*(2*log(sigma) + log(mu)) + (1/sigma^2 - 1)*log(y) - y/(sigma^2*mu) - lgamma(1/sigma^2)
    
    ll <- sum(weights*loglike)
    return(ll)
  }
  # the first dev of the log-likelihood given that the value of mu
  gn.param <- function(param.value, beta0, betaX, weights, y, X){
    sigma <- exp(param.value)
    mu <- exp(X%*%c(beta0, betaX))
    dldd = (2/sigma^3) * ((y/mu) - log(y) + log(mu) + log(sigma^2) - 1 + digamma(1/(sigma^2)))

    out <- sum(dldd*weights)
    return(out)
  }
  
  opt2 <-  maxLik(start = int.param,
                  logLik = fn.param,
                  grad  = gn.param,
                  beta0 = beta0, betaX = betaX,
                  weights = weights,  y = y, X = X, 
                  control = control)
  param <- exp(opt2$estimate)
  
  result <- list(beta0 = beta0, betaX = betaX, param = param)
  return(result)
}


# define the log-likelihood of GA distribution: weights*loglike - w*lambda*L1 - w*v*L2
GA.item <- function(Y, X, beta0, betaX, param, Sims, lambda, v, weights, w){
  pX <- ncol(X) - 1 # length of beta
  n <- length(Y)    # numbers of samples
  y <- Y
  modpars <- c(beta0, betaX, log(param)) # model parameters
  
  betaall <- modpars[1:(pX+1)]
  beta0 <- betaall[1]
  betaX <- betaall[2:(pX+1)]
  mu <- as.vector(exp(X%*%betaall)) # mu
  sigma <- exp(modpars[pX+2]) # sigma, dispersion parameter
  
  loglike <- dgamma(x = y, shape = rep(1/sigma^2, time = length(y)),
                    scale = sigma^2*mu,
                    log = TRUE)
  
  
  L1 <- sum(betaX^2)
  
  q.mat <-  matrix(0, nrow = pX, ncol = pX, byrow = T)
  for(j in 1:pX){
    for(k in 1:pX){
      q.mat[j,k] <- betaX[j] - betaX[k]
    }
  }
  L2 <- 0.5*sum(as.vector(abs(Sims*q.mat)))
  
  ll <- sum(weights*loglike) - w*lambda*L1 - w*v*L2
  
  results <- list(ll = ll,
                  L1 = L1, L2 = L2)
}

