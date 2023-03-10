library(maxLik)
library(caret)
library(data.table)
library(Metrics)

source("0-innerLayer-AMDD.r", encoding = 'UTF-8')
source("0-outerLayer-EM.r", encoding = 'UTF-8')

set.seed(1231)
rho.coef <- 0.5 # correlation between covariates
nsample <- 1000 # sample size
pcluster <- 5 # the number of clusters
beta01 <- 1; beta02 <- 2;  # intercepts
betaX1 <- c(rep(-0.1, each = pcluster), rep(-0.2, each = pcluster))
betaX2 <- c(rep(0.1, each = pcluster), rep(0.2, each = pcluster))
w <- c(0.7, 0.3) # weight
p <- length(betaX1) # the number of regression coefficients
cluster.label <- c(rep(2:3, each = pcluster)) # the label of cluster
S.sim <- matrix(0, nrow = p, ncol = p) # Similarity matrix
for(j in 1:p){
  for(k in 1:p){
    cluster1 <- cluster.label[j]
    cluster2 <- cluster.label[k]
    S.sim[j,k] <- ifelse(cluster1 == cluster2, rho.coef, 0)
  }
}
diag(S.sim) <- 1
# true values of parameters
true.w <- w
true.beta0 <- c(beta01, beta02)
true.betaX <- cbind(betaX1, betaX2)
true.param <- list(sigma = c(0.1, 0.2))
cov.sim <- S.sim*0.2*0.2 # covariance matrix,
# simulated covarties from the multivariate normal distribution
X <- mvtnorm::rmvnorm(nsample,
                      mean = rep(0, length.out = p),
                      sigma = cov.sim)
# regression coefficients for all components
beta0.mat <- c(beta01, beta02)
betaX.mat <- cbind(betaX1, betaX2)
mu1.mat <- exp(beta0.mat[1] + X%*%betaX.mat[,1])
mu2.mat <- exp(beta0.mat[2] + X%*%betaX.mat[,2])
# simulate the samples (y) from the mixture Gamma distribution
ysim <- 0
for(i in 1:nsample){
  usim <- runif(1, 0, 1)
  if(usim < w[1]){
    ysim[i] <- rGA(1, mu = mu1.mat[i], sigma = true.param$sigma[1]) # w1
  } else {
    ysim[i] <- rGA(1, mu = mu2.mat[i], sigma = true.param$sigma[2]) # w2
  }
}
X <- cbind(1, X) # design matrix 

# similar matrix for estimation 
Sims.mat <<- as.matrix(proxy::simil(X[, -1], by_rows = F), nrow = p, ncol = p)
Sims.mat[Sims.mat < 0] <- 0
Sims.mat[is.na(Sims.mat)] <- 0
Sims.mat[Sims.mat < 0.1] <- 0
Sims <- Sims.mat
diag(Sims) <- 1 # similar matrix
mydata <- data.table(ysim, X[,-1])
colnames(mydata) <- c('y', paste0('X', 1:p))


# ameter initialisation strategy
model.init <- function(y, X, H){
  library(glmnet)
  
  init.cluster <- kmeans(y, centers = H)
  label1 <- names(sort(table(init.cluster$cluster), decreasing = T))[1]
  label2 <- names(sort(table(init.cluster$cluster), decreasing = T))[2]
  label1 <- as.numeric(label1)
  label2 <- as.numeric(label2)
  init.w <- as.numeric(sort(table(init.cluster$cluster), decreasing = T) / length(y))
  init.beta0 <- numeric(H); init.betaX <- matrix(0, nrow = ncol(X)-1, ncol = H)
  init.param <- list(sigma = numeric(H))
  mod <- glmnet(x = X[, -1],
                y= y,
                family = Gamma(link = "log"), alpha = 0, lambda = 0.01)
  labels <- c(label1, label2)
  for(h in 1:H){
    mod <- glmnet(x = X[which(init.cluster$cluster == labels[h]), -1],
                  y= y[which(init.cluster$cluster == labels[h])],
                  family = Gamma(link = "log"), alpha = 0, lambda = 0.01)
    init.beta0[h] <- coef(mod)[1]
    init.betaX[, h] <- ifelse(is.na(coef(mod)[-1]), 0, coef(mod)[-1])
    init.param$sigma[h] <- optimize(f = function(s.value){
      -sum(dGA(y[which(init.cluster$cluster == labels[h])],
               mu = exp(X[which(init.cluster$cluster == labels[h]), ] %*% as.numeric(coef(mod))),
               sigma = s.value, log = T))
    }, interval = c(0, 5))$minimum
  }
  out <- list(init.w = init.w, init.beta0 = init.beta0,
              init.betaX = init.betaX,
              init.param = init.param)
  out
}  
 
# fit the model
ytrain <- mydata$y
Xtrain <- X
H <- 2
init.w <- model.init(ytrain, Xtrain, H = H)$init.w
init.beta0 <- model.init(ytrain, Xtrain, H = H)$init.beta0
init.betaX <- model.init(ytrain, Xtrain, H = H)$init.betaX
init.param <- model.init(ytrain, Xtrain, H = H)$init.param
  
mix_cov.cluster_mod <- MixGamma_EM(y = ytrain,
                                       X = Xtrain, 
                                       init.w = init.w,
                                       init.beta0 = init.beta0,
                                       init.betaX = init.betaX,
                                       init.param = init.param,
                                       Sims = Sims,
                                       lambda = 0.001,
                                       v =   1,
                                       rho =  1,
                                       control.admm = list(tol.primal = 0.05, tol.dual = 0.05,
                                                           maxit = 100, num.print = 50),
                                       control.EM = list(maxit = 10, tol = 0.01),
                                       control.NR = list(tol = 1e-4, reltol = 1e-4, gradtol = 1e-4))
betaX1.est <- mix_cov.cluster_mod$betaX[,1]
betaX2.est <- mix_cov.cluster_mod$betaX[,2]
w.est <- mix_cov.cluster_mod$w
  
cbind(betaX1.est, betaX2.est)
   


