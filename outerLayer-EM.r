# define the FRM with covariate clustering using EM-ADMM Algorithm 

library(maxLik)
library(caret)
library(gamlss)

MixGamma_EM <- function(y, init.w, init.beta0, init.betaX, init.param,
                       X, lambda = 1, v = 1, rho = 2, Sims,
                       control.admm = list(tol.primal = 1, num.print = 5, tol.dual = 1, maxit = 10),
                       control.EM = list(maxit = 10, tol = 1),
                       control.NR = list(tol = 1e-4, reltol = 1e-4, gradtol = 1e-4)){
  # y: the response variable.
  # X: the design matrix with dims n*(p+1), where p denotes the number of covariates, n is the sample size.
  # init.w: the initial parameters for the mixing proportions. a vector of length H.
  # init.beta0: a matrix with dims 1*H where H is the number of component for mixture distribution.
  # init.betaX: a matrix with dims (p+1)*H.
  # init.param: a list object containing the regression coefficients and dispersion parameter in Gamma: (sigma).
  # lambda: the penalty parameter for the L2-norm of betaX.
  # v: the penalty parameter for the similarity.
  # rho: the tuning parameter for ADMM.
  # control.admm: list of control parameters for ADMM.
  # control.EM: list of control parameters for EM.
  # control.NR. list of control parameters for Newton-Raphson optimization
  
  ## model parameters initialization
  T0.EM <- Sys.time()
  n <- length(y)
  pX <- nrow(init.betaX)
  H <- length(init.w)
  w <<- init.w
  beta0.mat <<- init.beta0
  betaX.mat <<- init.betaX
  param.mat <<- init.param
  ## auxiliary variables initialization
  z.arr <- r.arr <- array(0, dim = c(pX, pX, H)) # three-dimensional array
  for (k in 1:H){
    z.arr[,,k] <- matrix(rep(betaX.mat[,k], times = pX), nrow = pX, byrow = FALSE) # j is row, k is col
    r.arr[,,k] <- 0 # initialization
  }
  ## control parameters initialization
  tol.primal = control.admm$tol.primal
  tol.dual = control.admm$tol.dual
  ADMMmaxiters = control.admm$maxit
  EMmaxiters = control.EM$maxit
  EMtol = control.EM$tol
  mycontrol = control.NR
  if (is.null(control.admm$num.print)){
    num.print = 1} else {num.print = control.admm$num.print}
  
  # EM algorithm
  EM.loglike.trace <- c()
  for(iters.EM in 1:EMmaxiters){
    cat(paste('===========================',' The ',iters.EM, 'th iterations EM algorithm', '===========================', '\n', sep = ''))
    cat('Hyperparameters: ',paste(c('gamma:', 'v:', 'rho:'), c(lambda, v, rho), collapse ='; '), '\n')

    # ====================================================================
    # E step
    # beta0: intercept for each GA
    # betaX: coefficients for each GA
    # param: (sigma) for each GA
    # ====================================================================
  
    # 1. the log-likelihood of the mixture models.
    loglike.mGA <- function(beta0.mat, betaX.mat, param.mat, y, n, H, X, w){
      fy.ih <- dGA(x = matrix(y, nrow = n, ncol = H),
                    mu = exp(X %*% rbind(beta0.mat, betaX.mat)),
                    sigma = matrix(param.mat$sigma, nrow = n, ncol = H, byrow = T))
      wfy.ih <- matrix(w, nrow = n, ncol = H, byrow = T) * fy.ih
      fy.i <- apply(wfy.ih, 1, sum)
      loglike <- sum(log(fy.i))
      return(loglike)
    }
    # the current value of log-likelihood
    curr.loglike <- loglike.mGA(beta0.mat = beta0.mat, betaX.mat = betaX.mat,
                           param.mat = param.mat, y = y, n = n, H = H, X = X, w = w)
    # the density function of Gamma distribution for each component.
    fy.ih <- dGA(x = matrix(y, nrow = n, ncol = H),
                  mu = exp(X %*% rbind(beta0.mat, betaX.mat)),
                  sigma = matrix(param.mat$sigma, nrow = n, ncol = H, byrow = T))
    wfy.ih <- matrix(w, nrow = n, ncol = H, byrow = T) * fy.ih
    e.ih <<- wfy.ih / apply(wfy.ih, 1, sum)
    rm(fy.ih, wfy.ih)
    
    # ====================================================================
    # M step
    # ====================================================================
    # update w using the following or
    w <- apply(e.ih, 2, mean)
    # # 2. define the minimization objective function f~(w), for updating the mixing proportions.
    # fn.w <- function(pars, e.ih, betaX.mat, lambda, v, Sims){
    #   eta <- c(pars, 0)
    #   w.h <- exp(eta)/(sum(exp(eta)))
    #   out1 <- sum(log(w.h)* apply(e.ih, 2, sum))
    #   out2 <- sum(w.h*apply(betaX.mat, 2, function(z){sum(z^2)}))
    #   temp <- 0
    #   for(h in 1:H){
    #     betaX.h <- betaX.mat[,h]
    #     w.ih <- w.h[h]
    #     for(j in 1:pX){
    #       for(k in 1:pX){
    #         temp <- temp + w.ih*Sims[j,k]*abs(betaX.h[j] - betaX.h[k])
    #       }
    #     }
    #   }
    #   out3 <- temp
    #   loss <- - out1 + lambda*out2 + v*0.5*out3
    #   return(loss)
    # }
    # 
    # if(H == 2){
    #   eta.h <- optimize(f = fn.w, e.ih = e.ih, betaX.mat = betaX.mat, lambda = lambda, v = v, Sims = Sims, interval = c(-5, 5))$minimum;
    #   eta.h <- c(eta.h, 0)
    # } else if(H == 1){
    #   eta.h <- 0
    # } else {
    #   eta.h <- optim(par = log(w[-1]/(1-w[-1])), fn = fn.w, e.ih = e.ih,betaX.mat = betaX.mat, lambda = lambda, v = v, Sims = Sims)$par;
    #   eta.h <- c(eta.h, 0)
    # }
    # w <<- as.numeric(exp(eta.h)/sum(exp(eta.h))) # update the weight (mixing proportions)
    cat(paste('E-step: mixing weights are: ', paste(round(w,3), collapse = ', '), '\n'))
    
    
    ## Step 3: update the beta0, betaX and param for each GA component
  
    T1.ADMM <- Sys.time()
    current.beta.mat <- rbind(beta0.mat, betaX.mat) # current estimates of beta
    z.jk.mat <- z.kj.mat <- matrix(0, nrow = pX*pX, ncol = H)
    for(h in 1:H){
      cat(paste('M-step:', 'component', h, ':','\n'))
      for(iters.ADMM in 1:ADMMmaxiters){
        # current log-likelihhod
        current.item.h <- GA.item(Y = y, X = X,
                                   beta0 = beta0.mat[h],
                                   betaX = betaX.mat[,h],
                                   param = c(param.mat$sigma[h],
                                             param.mat$nu[h],
                                             param.mat$tau[h]),
                                   w = w[h],
                                   Sims = Sims, lambda = lambda, v = v, weights = e.ih[,h])
        current.z.mat <- z.arr[,,h] # current z mat
        
        # Step 1: update regression parameters in each component using GA.ADMM.seq function
        mod.h <- innerLayer_ADMM(Y = y, X = X,
                          beta0 = beta0.mat[h],
                          betaX = betaX.mat[,h],
                          param = c(param.mat$sigma[h]),
                          z.mat = z.arr[,,h], r.mat = r.arr[,,h],
                          lambda = lambda, rho = rho,
                          w = w[h],
                          weights = e.ih[,h], control = mycontrol)
        beta0.mat[h] <- mod.h$beta0
        betaX.mat[,h] <- mod.h$betaX
        param.mat$sigma[h] <- mod.h$param[1]

        ## Step 2: update z.jk and r.jk
        # convert matrix (z.mat, r.mat) into vector
        z.mat <- z.arr[,,h]
        r.mat <- r.arr[,,h]
        z.jk <- as.vector(t(z.mat)) # z[j,k] j=1, k = 1,2,...
        z.kj <- as.vector((z.mat)) # z[k,j]
        r.jk <- as.vector(t(r.mat))
        r.kj <- as.vector((r.mat))
        s.jk <- as.vector(t(Sims))
        beta.jk <- rep(betaX.mat[,h], each = pX)
        beta.kj <- rep(betaX.mat[,h], times = pX)
        
        hh <- abs(((beta.jk - r.jk - beta.kj + r.kj)))
        theta <- pmax(1 - (v * w[h] / rho * as.numeric(s.jk) / hh), 0.5)
        theta[hh == 0] <- 0.5
        
        # update z.jk and r.jk
        z.jk.new <- theta*(beta.jk - r.jk) + (1 - theta)*(beta.kj - r.kj)
        z.kj.new <- (1-theta)*(beta.jk - r.jk) + theta*(beta.kj - r.kj)
        r.jk.new <- r.jk + (z.jk.new - beta.jk)
        
        
        # update z.arr and r.arr
        z.arr[,,h] <- matrix(z.jk.new, nrow = pX, ncol = pX, byrow = TRUE)
        r.arr[,,h] <- matrix(r.jk.new, nrow = pX, ncol = pX, byrow = TRUE)
        
        rm(theta, hh, z.mat, r.mat)
        
        # calculate the new log-likelihood of GA regression in h component
        new.item.h <- GA.item( Y = y, X = X,
                                beta0 = beta0.mat[h],
                                betaX = betaX.mat[,h],
                                param = c(param.mat$sigma[h],
                                          param.mat$nu[h],
                                          param.mat$tau[h]), 
                                w = w[h],
                                Sims = Sims, lambda = lambda, v = v, weights = e.ih[,h])
        
        new.z.mat <- z.arr[,,h]
        new.beta.j.mat <- matrix(betaX.mat[,h], nrow = pX, ncol = pX, byrow = FALSE)
        
        #  computes primal and dual residuals for checking the convergence
        res.primal <- norm(new.z.mat - new.beta.j.mat, type = 'F') # difference between auxiliary variables and the beta's
        res.dual <- norm(rho*(new.z.mat - current.z.mat), type = 'F') # the changing of auxiliary variables
        
        if (iters.ADMM %% num.print == 0) {
          cat(' Iter',  iters.ADMM, ':',
              paste(c('Loglike:', 'L2-norm:', 'Simil-norm:', 'res.primal:', 'res.dual: '), 
                    round(c(new.item.h$ll, new.item.h$L1, new.item.h$L2, 
                            res.primal, res.dual), 2), collapse ='; '), '\n')
        }
        # the convergence conditions
        if((res.primal < tol.primal)&(res.dual < tol.dual)) break # end when little improvement
        
      }
      z.jk.mat[,h] <- z.jk.new
      z.kj.mat[,h] <- z.kj.new
    }
    
    rm(h)
    T3.ADMM <- Sys.time()
  
    # the updated loglikelihood
    new.loglike <- loglike.mGA(beta0.mat = beta0.mat, betaX.mat = betaX.mat,
                                param.mat = param.mat, y = y, n = n, H = H, X = X, w = w)
    # the updated parameters
    new.beta.mat <- rbind(beta0.mat, betaX.mat)
    # the changes of estimated parameters
    delta.EM <- norm((new.beta.mat - current.beta.mat), type = 'F')
    
    cat(paste(c('EM: Loglike: ', 'delta: '), 
              round(c(new.loglike, delta.EM), 3),  collapse ='; '), '\n')
    
    EM.loglike.trace[iters.EM] <- new.loglike
    
    if(delta.EM >0 & delta.EM < EMtol) break
  }
  
  T1.EM <- Sys.time()
  cat(paste('Running time of EM algorithm: ', 
            round(T1.EM - T0.EM, 3),'\n', sep = ''))
  

  
  return(list(w = w, beta0 = beta0.mat, betaX = betaX.mat, 
              param = param.mat,
              EM.loglike.trace = EM.loglike.trace,
              EM.loglike = new.loglike,
              z.jk.mat = z.jk.mat,
              z.kj.mat = z.kj.mat
              ))

}



