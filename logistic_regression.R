#install.packages('bayesplot')
#install.packages('cowplot')

############################################
###### Data Simulation
set.seed(1)

n <- 50 # number of samples
q <- 3 # number of covariates
d <- q+1 # dimensionality of parameters
mu <- 2    # mean for data samples
sigma <- sqrt(2)    #sd for data samples

X <- matrix(1,n,d) # matrix to hold data (inlcudin)

for(i in 1:n){
  X[i,-1] <-rnorm(q,mu,sigma)
}

############################################
#### Model Specifications
n <- 50   # number of samples
m <- 20     # parameter for the binomial response
priormu <- 0 # mean of gaussian prior for betas
priorsigma <- 3 # standard deviation of gaussian prior for betas
sets <- 50 # number of chains of posterior and prior samples
iterations <- 11000
warmup <- 1000 # number of iterations from each chain (5000 will be used as burn in)
N <- iterations - warmup # size of our posterior and prior samples

betas <- sample(c(-1.5,-1,-0.5,0.5,1,1.5), 4, replace = TRUE) # randomly chosen beta values
############################################
#### Generating Response from True Model

p <- 1/(1+exp(-(X%*%betas))) # p values of beta values in TRUE model

y <- numeric(n)
for( i in 1:n){
  y[i] <- rbinom(1,m,p[i])
}

###############################################################
### FUNCTIONS
### We create the functions here that we need for our methods

likelihood <- function(y,n,p){ prod(dbinom(y,m,p))}

loglikelihood <- function(y,n,p){sum(dbinom(y,m,p, log = TRUE))}

prior <- function(x,priormu,priorsigma){prod(dnorm(x,priormu,priorsigma))}

logprior <- function(x,priormu,priorsigma){sum(dnorm(x,priormu,priorsigma, log = TRUE))}

###############################################################
### Data for each model

X1 <- matrix(X[,2], ncol = 1)
X1.1 <- cbind(rep(1,n),X1)    ## .1 represents adding a column of 1s to matrix
X2 <- matrix(X[,3], ncol = 1)
X2.1 <- cbind(rep(1,n),X2)
X3 <- matrix(X[,4], ncol = 1)
X3.1 <- cbind(rep(1,n),X3)
X4 <- X[,c(2,3)]
X4.1 <- cbind(rep(1,n),X4)
X5 <- X[,c(2,4)]
X5.1 <- cbind(rep(1,n),X5)
X6 <- X[,c(3,4)]
X6.1 <- cbind(rep(1,n),X6)
X7 <- X[,c(2,3,4)]
X7.1 <- cbind(rep(1,n),X7)

##############################################################################
##############################################################################
##### MODEL 1 (BETA = (beta0, beta1, 0, 0))
##############################################################################
############################################
#### Sampling from the posterior
set.seed(1)

q <- 1    
d <- q+1
library(rstan)

our_model <- stan_model('logistic_regression.stan')
our_data <- list(n = n, q = q, X= X1, y = y, priormu = priormu, priorsigma = priorsigma, m = m)

our_sample <- sampling(our_model, data = our_data, chains = sets, iter = iterations, warmup = warmup, seed = 1)
samples <- as.array(our_sample)

chains <- list() # chains is a list which hold the posterior simulations
for (i in 1:sets){
  chains[[i]] <- samples[,i,-(d+1)]
}

chainpvals <- list() # chainpvals is a list which holds the p values of our posterior simulations
for(i in 1:sets){
  chainpvals[[i]] <- t(apply(chains[[i]],1,function(x) {1/(1+exp(-(X1.1%*%x)))}))
}

#### Sampling from the prior

Nprior <-  N # number of samples from prior

priorsims <- list() # priorsims is a list which holds our prior simulations
for( i in 1:sets){
  priorsims[[i]] <- matrix(rnorm(Nprior*d,priormu,priorsigma), nrow = Nprior, ncol = d)
}

priorpvals <- list()  # priorpvals is a list which holds the p values of our prior simulations
for(i in 1:sets){
  priorpvals[[i]] <- t(apply(priorsims[[i]],1,function(x) {1/(1+exp(-(X1.1%*%x)))}))
}
############################################
#### METHODS
############################################
### MONTE CARLO ESTIMATOR
set.seed(1)

MClikelihoods <- matrix(0,sets,Nprior)

for (i in 1:sets){
  MClikelihoods[i,] <- apply(priorpvals[[i]],1,likelihood,y = y, n = m)
}

MCloglikelihoods <- matrix(0,sets,Nprior)

for (i in 1:sets){
  MCloglikelihoods[i,] <- apply(priorpvals[[i]],1,loglikelihood,y = y, n = m)
}

MCestimates1 <- numeric(sets)
for (i in 1:sets) { MCestimates1[i] <- log(1/Nprior*sum(MClikelihoods[i,]))}

MCestimate1 <- mean(MCestimates1)
MCsd1 <- sqrt(var(MCestimates1))

MCestimate1
MCsd1

############################################
### HARMONIC MEAN ESTIMATOR
set.seed(1)

HMlikelihoods <- matrix(0,sets,N)
for (i in 1:sets){
  HMlikelihoods[i,] <- apply(chainpvals[[i]],1,likelihood,y = y, n = m)
}

HMestimates1 <- apply(HMlikelihoods,1,function(x) {-log(1/N*sum(1/x))})

HMestimate1 <- mean(HMestimates1)
HMsd1 <- sqrt(var(HMestimates1))

HMestimate1
HMsd1

############################################
### GENERALISED HARMONIC MEAN ESTIMATOR
set.seed(1)
library(mvtnorm)

GHMpriors <- matrix(0,sets,N)
for(i in 1:sets){
  GHMpriors[i,] <- apply(chains[[i]],1,prior,priormu = priormu, priorsigma = priorsigma)
}

GHMproposaldensities <- matrix(0,sets,N)

for(i in 1:sets){
  mean <- apply(chains[[i]],2,mean)
  var <- var(chains[[i]])
  GHMproposaldensities[i,] <- apply(chains[[i]],1,dmvnorm,mean = mean, sigma = var)
}

matrix <- GHMproposaldensities/(HMlikelihoods*GHMpriors)

GHMestimates1 <- apply(matrix,1,function(x) {-log(1/N*sum(x))})
GHMestimate1 <- mean(GHMestimates1)
GHMsd1 <- sqrt(var(GHMestimates1))

GHMestimate1
GHMsd1
############################################################################################################################
## NEWTON-RAFTERY ESTIMATOR
set.seed(1)

delta = 0.2
NRiterations = 100

indices <- matrix(0,sets,delta*N)
for (i in 1:sets) {
  indices[i,] <- sample(N,delta*N)
}

NRsims <- chains

for (i in 1:sets){
  for(j in indices[i,]){
    NRsims[[i]][j,] <- priorsims[[i]][j,]
  }}

NRpvals <- list()
for(i in 1:sets){
  NRpvals[[i]] <- t(apply(NRsims[[i]],1,function(x) {1/(1+exp(-(X1.1%*%x)))}))
}

NRlikelihoods <- matrix(0,sets,N)
for (i in 1:sets){
  NRlikelihoods[i,] <- apply(NRpvals[[i]],1,likelihood,y = y, n = m)
}

NRestimates1 <- exp(HMestimates1)

for (i in 1:sets){
  for (j in 1:NRiterations) {
    NRestimates1[i] <- (sum(NRlikelihoods[i,]*(delta*NRestimates1[i]+(1-delta)*NRlikelihoods[i,])^-1))/(sum((delta*NRestimates1[i]+(1-delta)*NRlikelihoods[i,])^-1))
  }
}

NRestimates1 <- log(NRestimates1)
NRestimate1 <- mean(NRestimates1)
NRsd1 <- sqrt(var(NRestimates1))

NRestimate1
NRsd1

############################################
### LAPLACE-METROPOLIS ESTIMATOR
set.seed(1)

logposteriors <- matrix(0,sets,N)
for ( i in 1:sets) {
  for( j in 1:N) {
    logposteriors[i,j] <- loglikelihood(y,m,chainpvals[[i]][j,])+logprior(chains[[i]][j,],priormu,priorsigma)
  }}

maxlambda <- numeric(sets)
for ( i in 1:sets ) {maxlambda[i] <-  which.max(logposteriors[i,])}

S <- list()
for (i in 1:sets) {S[[i]] <- var(chains[[i]])}

LMestimates1 <- numeric(sets)
for (i in 1:sets) { LMestimates1[i] <- (d/2)*log(2*pi)+0.5*log(det(S[[i]]))+logposteriors[i,][maxlambda[i]]}

LMestimate1 <- mean(LMestimates1)
LMsd1 <- sqrt(var(LMestimates1))

LMestimate1
LMsd1

############################################
### BRIDGE SAMPLING ESTIMATOR
set.seed(1)

BSiterations=1000

N1=N
N2=N

s1=N1/(N1+N2)
s2=N2/(N1+N2)

proposalsims <- list()

for(i in 1:sets){
  mean <- apply(chains[[i]],2,mean)
  var <- var(chains[[i]])
  proposalsims[[i]] <- rmvnorm(N2,mean = mean, sigma = var)
}

proposalpvals <- list()

for(i in 1:sets){
  proposalpvals[[i]] <- t(apply(proposalsims[[i]],1,function(x) {1/(1+exp(-(X1.1%*%x)))}))
}

likelihoods2 <- matrix(0,sets,N2)

for ( i in 1:sets) {
  likelihoods2[i,] <- apply(proposalpvals[[i]],1,likelihood, y = y, n = m)
}

priors2 <- matrix(0,sets,N2)
for ( i in 1:sets) {
  priors2[i,] <- apply(proposalsims[[i]],1,prior, priormu = priormu, priorsigma = priorsigma)
}

proposaldensities2 <- matrix(0,sets,N2)

for(i in 1:sets){
  mean <- apply(chains[[i]],2,mean)
  var <- var(chains[[i]])
  proposaldensities2[i,] <- apply(proposalsims[[i]],1,dmvnorm, mean = mean, sigma = var)
}

likelihoods1 <- HMlikelihoods

priors1 <- GHMpriors

proposaldensities1 <- GHMproposaldensities

x2 <- matrix(0,sets,N2)
for (i in 1:sets){
  x2[i,] <- rep(exp(LMestimate1),N2)
}

x1 <- matrix(0,sets,N1)
for (i in 1:sets){
  x1[i,] <- rep(exp(LMestimate1),N1)
}


for(i in 1:BSiterations){
  matrix2 <- (likelihoods2*priors2)/(s1*likelihoods2*priors2+s2*x2*proposaldensities2)
  matrix1 <- proposaldensities1/(s1*likelihoods1*priors1+s2*x1*proposaldensities1)
  numerator <- apply(matrix2,1,function(x) {1/N2*sum(x)})
  denominator <- apply(matrix1,1,function(x) {1/N1*sum(x)})
  BSestimates1 <- numerator/denominator
  
  x2 <- matrix(0,sets,N2)
  for (i in 1:sets){
    x2[i,] <- rep(exp(BSestimates1[i]),N2)
  }
  
  x1 <- matrix(0,sets,N1)
  for (i in 1:sets){
    x1[i,] <- rep(exp(BSestimates1[i]),N1)
  }
  
}

BSestimates1 <- log(BSestimates1)
BSestimate1 <- mean(BSestimates1)
BSsd1 <- sqrt(var(BSestimate1))

BSestimate1
BSsd1

############################################
### FOURIER INTEGRAL ESTIMATOR
#set.seed(1)
#R <- 20

#centres <- matrix(0,sets,d) ### the centres of our Fourier estimation

#for( i in 1:sets){
#  centres[i,] <- apply(chains[[i]],2,mean) 
#}

#matrix <- matrix(0,sets,N)
#for( x in 1:sets){
#  matrix[x,] <- apply(chains[[x]],1,function(y){prod(sin(R*(centres[x,]-y))/(centres[x,]-y))})
#}

#densityestimates <- apply(matrix,1,function(x){1/(N*pi^d)*sum(x)})

#centrepvals <- t(apply(centres,1,function(x) {1/(1+exp(-(X1.1%*%x)))}))

#likelihoods <- apply(centrepvals,1,likelihood,y = y,n = m)

#priors <- numeric(sets)

#for(i in 1:sets){
#  priors[i] <- prior(centres[i,],priormu = priormu, priorsigma = priorsigma)
#}

#FIestimates <- log(likelihoods*priors)+log(densityestimates)

#FIestimate <- mean(FIestimates)
#FIsd <- sqrt(var(FIestimates))

#print(FIestimate)
#print(FIsd)

############################################
#### MODEL 1 - RESULTS
MCestimate1
MCsd1
HMestimate1
HMsd1
GHMestimate1
GHMsd1
NRestimate1
NRsd1
LMestimate1
LMsd1
BSestimate1
BSsd1
##############################################################################
##############################################################################
##### MODEL 2 (BETA = (beta0, 0 ,beta2, 0))
##############################################################################
##############################################################################
############################################
##### Sampling from the posterior
set.seed(1)

q <- 1    
d <- q+1
library(rstan)

our_model <- stan_model('logistic_regression.stan')
our_data <- list(n = n, q = q, X= X2, y = y, priormu = priormu, priorsigma = priorsigma, m = m)

our_sample <- sampling(our_model, data = our_data, chains = sets, iter = iterations, warmup = warmup, seed = 1)
samples <- as.array(our_sample)

chains <- list() # chains is a list which hold the posterior simulations
for (i in 1:sets){
  chains[[i]] <- samples[,i,-(d+1)]
}

chainpvals <- list() # chainpvals is a list which holds the p values of our posterior simulations
for(i in 1:sets){
  chainpvals[[i]] <- t(apply(chains[[i]],1,function(x) {1/(1+exp(-(X2.1%*%x)))}))
}

#### Sampling from the prior

priorsims <- list() # priorsims is a list which holds our prior simulations
for( i in 1:sets){
  priorsims[[i]] <- matrix(rnorm(Nprior*d,priormu,priorsigma), nrow = Nprior, ncol = d)
}

priorpvals <- list() # priorpvals is a list which holds thh p values of our prior simulations
for(i in 1:sets){
  priorpvals[[i]] <- t(apply(priorsims[[i]],1,function(x) {1/(1+exp(-(X2.1%*%x)))}))
}

############################################
#### METHODS
############################################
### MONTE CARLO ESTIMATOR
set.seed(1)

Nprior = N # number of samples from prior

MClikelihoods <- matrix(0,sets,Nprior)

for (i in 1:sets){
  MClikelihoods[i,] <- apply(priorpvals[[i]],1,likelihood,y = y, n = m)
}

MCestimates2 <- numeric(sets)
for (i in 1:sets) { MCestimates2[i] <- log(1/Nprior*sum(MClikelihoods[i,]))}

MCestimate2 <- mean(MCestimates2)
MCsd2 <- sqrt(var(MCestimates2))

MCestimate2
MCsd2

############################################
### HARMONIC MEAN ESTIMATOR
set.seed(1)

HMlikelihoods <- matrix(0,sets,N)
for (i in 1:sets){
  HMlikelihoods[i,] <- apply(chainpvals[[i]],1,likelihood,y = y, n = m)
}

HMestimates2 <- apply(HMlikelihoods,1,function(x) {-log(1/N*sum(1/x))})

HMestimate2 <- mean(HMestimates2)
HMsd2 <- sqrt(var(HMestimates2))

HMestimate2
HMsd2
############################################
### GENERALISED HARMONIC MEAN ESTIMATOR
set.seed(1)
library(mvtnorm)

GHMpriors <- matrix(0,sets,N)
for(i in 1:sets){
  GHMpriors[i,] <- apply(chains[[i]],1,prior,priormu = priormu, priorsigma = priorsigma)
}

GHMproposaldensities <- matrix(0,sets,N)

for(i in 1:sets){
  mean <- apply(chains[[i]],2,mean)
  var <- var(chains[[i]])
  GHMproposaldensities[i,] <- apply(chains[[i]],1,dmvnorm,mean = mean, sigma = var)
}

matrix <- GHMproposaldensities/(HMlikelihoods*GHMpriors)

GHMestimates2 <- apply(matrix,1,function(x) {-log(1/N*sum(x))})
GHMestimate2 <- mean(GHMestimates2)
GHMsd2 <- sqrt(var(GHMestimates2))

GHMestimate2
GHMsd2
############################################################################################################################
## NEWTON-RAFTERY ESTIMATOR
set.seed(1)

delta = 0.2
NRiterations = 100

indices <- matrix(0,sets,delta*N)
for (i in 1:sets) {
  indices[i,] <- sample(N,delta*N)
}

NRsims <- chains

for (i in 1:sets){
  for(j in indices[i,]){
    NRsims[[i]][j,] <- priorsims[[i]][j,]
  }}

NRpvals <- list()
for(i in 1:sets){
  NRpvals[[i]] <- t(apply(NRsims[[i]],1,function(x) {1/(1+exp(-(X2.1%*%x)))}))
}

NRlikelihoods <- matrix(0,sets,N)
for (i in 1:sets){
  NRlikelihoods[i,] <- apply(NRpvals[[i]],1,likelihood,y = y, n = m)
}


NRestimates2 <- exp(HMestimates2)

for (i in 1:sets){
  for (j in 1:NRiterations) {
    NRestimates2[i] <- (sum(NRlikelihoods[i,]*(delta*NRestimates2[i]+(1-delta)*NRlikelihoods[i,])^-1))/(sum((delta*NRestimates2[i]+(1-delta)*NRlikelihoods[i,])^-1))
  }
}

NRestimates2 <- log(NRestimates2)
NRestimate2 <- mean(NRestimates2)
NRsd2 <- sqrt(var(NRestimates2))

NRestimate2
NRsd2

############################################
### LAPLACE-METROPOLIS ESTIMATOR
set.seed(1)

logposteriors <- matrix(0,sets,N)
for ( i in 1:sets) {
  for( j in 1:N) {
    logposteriors[i,j] <- loglikelihood(y,m,chainpvals[[i]][j,])+logprior(chains[[i]][j,],priormu,priorsigma)
  }}

maxlambda <- numeric(sets)
for ( i in 1:sets ) {maxlambda[i] <-  which.max(logposteriors[i,])}

S <- list()
for (i in 1:sets) {S[[i]] <- var(chains[[i]])}

LMestimates2 <- numeric(sets)
for (i in 1:sets) { LMestimates2[i] <- (d/2)*log(2*pi)+0.5*log(det(S[[i]]))+logposteriors[i,][maxlambda[i]]}

LMestimate2 <- mean(LMestimates2)
LMsd2 <- sqrt(var(LMestimates2))

LMestimate2
LMsd2

############################################
### BRIDGE SAMPLING ESTIMATOR
set.seed(1)

BSiterations=1000

N1=N
N2=N

s1=N1/(N1+N2)
s2=N2/(N1+N2)

proposalsims <- list()

for(i in 1:sets){
  mean <- apply(chains[[i]],2,mean)
  var <- var(chains[[i]])
  proposalsims[[i]] <- rmvnorm(N2,mean = mean, sigma = var)
}

proposalpvals <- list()

for(i in 1:sets){
  proposalpvals[[i]] <- t(apply(proposalsims[[i]],1,function(x) {1/(1+exp(-(X2.1%*%x)))}))
}

likelihoods2 <- matrix(0,sets,N2)

for ( i in 1:sets) {
  likelihoods2[i,] <- apply(proposalpvals[[i]],1,likelihood, y = y, n = m)
}

priors2 <- matrix(0,sets,N2)
for ( i in 1:sets) {
  priors2[i,] <- apply(proposalsims[[i]],1,prior, priormu = priormu, priorsigma = priorsigma)
}

proposaldensities2 <- matrix(0,sets,N2)

for(i in 1:sets){
  mean <- apply(chains[[i]],2,mean)
  var <- var(chains[[i]])
  proposaldensities2[i,] <- apply(proposalsims[[i]],1,dmvnorm, mean = mean, sigma = var)
}

likelihoods1 <- HMlikelihoods

priors1 <- GHMpriors

proposaldensities1 <- GHMproposaldensities

x2 <- matrix(0,sets,N2)
for (i in 1:sets){
  x2[i,] <- rep(exp(LMestimate2),N2)
}

x1 <- matrix(0,sets,N1)
for (i in 1:sets){
  x1[i,] <- rep(exp(LMestimate2),N1)
}


for(i in 1:BSiterations){
  matrix2 <- (likelihoods2*priors2)/(s1*likelihoods2*priors2+s2*x2*proposaldensities2)
  matrix1 <- proposaldensities1/(s1*likelihoods1*priors1+s2*x1*proposaldensities1)
  numerator <- apply(matrix2,1,function(x) {1/N2*sum(x)})
  denominator <- apply(matrix1,1,function(x) {1/N1*sum(x)})
  BSestimates2 <- numerator/denominator
  
  x2 <- matrix(0,sets,N2)
  for (i in 1:sets){
    x2[i,] <- rep(exp(BSestimates2[i]),N2)
  }
  
  x1 <- matrix(0,sets,N1)
  for (i in 1:sets){
    x1[i,] <- rep(exp(BSestimates2[i]),N1)
  }
  
}

BSestimates2 <- log(BSestimates2)
BSestimate2 <- mean(BSestimates2)
BSsd2 <- sqrt(var(BSestimates2))

BSestimate2
BSsd2

############################################
#### MODEL 2 - RESULTS
MCestimate2
MCsd2
HMestimate2
HMsd2
GHMestimate2
GHMsd2
NRestimate2
NRsd2
LMestimate2
LMsd2
BSestimate2
BSsd2

##############################################################################
##############################################################################
##### MODEL 3 (BETA = (beta0, 0, 0, beta3))
##############################################################################
############################################
#### Sampling from the posterior

set.seed(1)

q <- 1    
d <- q+1
library(rstan)

our_model <- stan_model('logistic_regression.stan')
our_data <- list(n = n, q = q, X= X3, y = y, priormu = priormu, priorsigma = priorsigma, m = m)

our_sample <- sampling(our_model, data = our_data, chains = sets, iter = iterations, warmup = warmup, seed = 1)
samples <- as.array(our_sample)

chains <- list()  # chains is a list which hold the posterior simulations
for (i in 1:sets){
  chains[[i]] <- samples[,i,-(d+1)]
}

chainpvals <- list() # chainpvals is a list which holds the p values of our posterior simulations
for(i in 1:sets){
  chainpvals[[i]] <- t(apply(chains[[i]],1,function(x) {1/(1+exp(-(X3.1%*%x)))}))
}

### Sampling from the prior

priorsims <- list() # priorsims is a list which holds our prior simulations
for( i in 1:sets){
  priorsims[[i]] <- matrix(rnorm(Nprior*d,priormu,priorsigma), nrow = Nprior, ncol = d)
}

priorpvals <- list() # priorpvals is a list which holds the p values of our prior simulations
for(i in 1:sets){
  priorpvals[[i]] <- t(apply(priorsims[[i]],1,function(x) {1/(1+exp(-(X3.1%*%x)))}))
}


############################################
#### METHODS
############################################
### MONTE CARLO ESTIMATOR
set.seed(1)

Nprior = N # number of samples from prior

MClikelihoods <- matrix(0,sets,Nprior)

for (i in 1:sets){
  MClikelihoods[i,] <- apply(priorpvals[[i]],1,likelihood,y = y, n = m)
}

MCestimates3 <- numeric(sets)
for (i in 1:sets) { MCestimates3[i] <- log(1/Nprior*sum(MClikelihoods[i,]))}

MCestimate3 <- mean(MCestimates3)
MCsd3<- sqrt(var(MCestimates3))

MCestimate3
MCsd3

############################################
### HARMONIC MEAN ESTIMATOR
set.seed(1)

HMlikelihoods <- matrix(0,sets,N)
for (i in 1:sets){
  HMlikelihoods[i,] <- apply(chainpvals[[i]],1,likelihood,y = y, n = m)
}

HMestimates3 <- apply(HMlikelihoods,1,function(x) {-log(1/N*sum(1/x))})

HMestimate3 <- mean(HMestimates3)
HMsd3 <- sqrt(var(HMestimates3))

HMestimate3
HMsd3

############################################
### GENERALISED HARMONIC MEAN ESTIMATOR
set.seed(1)
library(mvtnorm)

GHMpriors <- matrix(0,sets,N)
for(i in 1:sets){
  GHMpriors[i,] <- apply(chains[[i]],1,prior,priormu = priormu, priorsigma = priorsigma)
}

GHMproposaldensities <- matrix(0,sets,N)

for(i in 1:sets){
  mean <- apply(chains[[i]],2,mean)
  var <- var(chains[[i]])
  GHMproposaldensities[i,] <- apply(chains[[i]],1,dmvnorm,mean = mean, sigma = var)
}

matrix <- GHMproposaldensities/(HMlikelihoods*GHMpriors)

GHMestimates3 <- apply(matrix,1,function(x) {-log(1/N*sum(x))})
GHMestimate3 <- mean(GHMestimates3)
GHMsd3 <- sqrt(var(GHMestimates3))

GHMestimate3
GHMsd3
############################################################################################################################
## NEWTON-RAFTERY ESTIMATOR
set.seed(1)

delta = 0.2
NRiterations = 100

indices <- matrix(0,sets,delta*N)
for (i in 1:sets) {
  indices[i,] <- sample(N,delta*N)
}

NRsims <- chains

for (i in 1:sets){
  for(j in indices[i,]){
    NRsims[[i]][j,] <- priorsims[[i]][j,]
  }}

NRpvals <- list()
for(i in 1:sets){
  NRpvals[[i]] <- t(apply(NRsims[[i]],1,function(x) {1/(1+exp(-(X3.1%*%x)))}))
}

NRlikelihoods <- matrix(0,sets,N)
for (i in 1:sets){
  NRlikelihoods[i,] <- apply(NRpvals[[i]],1,likelihood,y = y, n = m)
}

NRestimates3 <- exp(HMestimates3)

for (i in 1:sets){
  for (j in 1:NRiterations) {
    NRestimates3[i] <- (sum(NRlikelihoods[i,]*(delta*NRestimates3[i]+(1-delta)*NRlikelihoods[i,])^-1))/(sum((delta*NRestimates3[i]+(1-delta)*NRlikelihoods[i,])^-1))
  }
}

NRestimates3 <- log(NRestimates3)
NRestimate3 <- mean(NRestimates3)
NRsd3 <- sqrt(var(NRestimates3))

NRestimate3
NRsd3

############################################
### LAPLACE-METROPOLIS ESTIMATOR
set.seed(1)

logposteriors <- matrix(0,sets,N)
for ( i in 1:sets) {
  for( j in 1:N) {
    logposteriors[i,j] <- loglikelihood(y,m,chainpvals[[i]][j,])+logprior(chains[[i]][j,],priormu,priorsigma)
  }}

maxlambda <- numeric(sets)
for ( i in 1:sets ) {maxlambda[i] <-  which.max(logposteriors[i,])}

S <- list()
for (i in 1:sets) {S[[i]] <- var(chains[[i]])}

LMestimates3 <- numeric(sets)
for (i in 1:sets) { LMestimates3[i] <- (d/2)*log(2*pi)+0.5*log(det(S[[i]]))+logposteriors[i,][maxlambda[i]]}

LMestimate3 <- mean(LMestimates3)
LMsd3 <- sqrt(var(LMestimates3))

LMestimate3
LMsd3

############################################
### BRIDGE SAMPLING ESTIMATOR
set.seed(1)

BSiterations=1000

N1=N
N2=N

s1=N1/(N1+N2)
s2=N2/(N1+N2)

proposalsims <- list()

for(i in 1:sets){
  mean <- apply(chains[[i]],2,mean)
  var <- var(chains[[i]])
  proposalsims[[i]] <- rmvnorm(N2,mean = mean, sigma = var)
}

proposalpvals <- list()

for(i in 1:sets){
  proposalpvals[[i]] <- t(apply(proposalsims[[i]],1,function(x) {1/(1+exp(-(X3.1%*%x)))}))
}

likelihoods2 <- matrix(0,sets,N2)

for ( i in 1:sets) {
  likelihoods2[i,] <- apply(proposalpvals[[i]],1,likelihood, y = y, n = m)
}

priors2 <- matrix(0,sets,N2)
for ( i in 1:sets) {
  priors2[i,] <- apply(proposalsims[[i]],1,prior, priormu = priormu, priorsigma = priorsigma)
}

proposaldensities2 <- matrix(0,sets,N2)

for(i in 1:sets){
  mean <- apply(chains[[i]],2,mean)
  var <- var(chains[[i]])
  proposaldensities2[i,] <- apply(proposalsims[[i]],1,dmvnorm, mean = mean, sigma = var)
}

likelihoods1 <- HMlikelihoods

priors1 <- GHMpriors

proposaldensities1 <- GHMproposaldensities

x2 <- matrix(0,sets,N2)
for (i in 1:sets){
  x2[i,] <- rep(exp(LMestimate3),N2)
}

x1 <- matrix(0,sets,N1)
for (i in 1:sets){
  x1[i,] <- rep(exp(LMestimate3),N1)
}


for(i in 1:BSiterations){
  matrix2 <- (likelihoods2*priors2)/(s1*likelihoods2*priors2+s2*x2*proposaldensities2)
  matrix1 <- proposaldensities1/(s1*likelihoods1*priors1+s2*x1*proposaldensities1)
  numerator <- apply(matrix2,1,function(x) {1/N2*sum(x)})
  denominator <- apply(matrix1,1,function(x) {1/N1*sum(x)})
  BSestimates3 <- numerator/denominator
  
  x2 <- matrix(0,sets,N2)
  for (i in 1:sets){
    x2[i,] <- rep(exp(BSestimates3[i]),N2)
  }
  
  x1 <- matrix(0,sets,N1)
  for (i in 1:sets){
    x1[i,] <- rep(exp(BSestimates3[i]),N1)
  }
  
}

BSestimates3 <- log(BSestimates3)
BSestimate3 <- mean(BSestimates3)
BSsd3 <- sqrt(var(BSestimates3))

BSestimate3
BSsd3


############################################
#### MODEL 3 - RESULTS
MCestimate3
MCsd3
HMestimate3
HMsd3
GHMestimate3
GHMsd3
NRestimate3
NRsd3
LMestimate3
LMsd3
BSestimate3
BSsd3

##############################################################################
##############################################################################
##### MODEL 4 (BETA = (beta0, beta1, beta2, 0))
##############################################################################
############################################
#### Sampling from the posterior
set.seed(1)

q <- 2    
d <- q+1
library(rstan)

our_model <- stan_model('logistic_regression.stan')
our_data <- list(n = n, q = q, X= X4, y = y, priormu = priormu, priorsigma = priorsigma, m = m)

our_sample <- sampling(our_model, data = our_data, chains = sets, iter = iterations, warmup = warmup, seed = 1)
samples <- as.array(our_sample)

chains <- list() # chains is a list which hold the posterior simulations
for (i in 1:sets){
  chains[[i]] <- samples[,i,-(d+1)]
}

chainpvals <- list() # chainpvals is a list which holds the p values of our posterior simulations
for(i in 1:sets){
  chainpvals[[i]] <- t(apply(chains[[i]],1,function(x) {1/(1+exp(-(X4.1%*%x)))}))
}

#### Sampling from the prior

priorsims <- list() # priorsims is a list which holds our prior simulations
for( i in 1:sets){
  priorsims[[i]] <- matrix(rnorm(Nprior*d,priormu,priorsigma), nrow = Nprior, ncol = d)
}

priorpvals <- list() # priorpvals is a list which holds the p values of our prior simulations
for(i in 1:sets){
  priorpvals[[i]] <- t(apply(priorsims[[i]],1,function(x) {1/(1+exp(-(X4.1%*%x)))}))
}

############################################
#### METHODS
############################################
### MONTE CARLO ESTIMATOR
set.seed(1)

Nprior = N # number of samples from prior

MClikelihoods <- matrix(0,sets,Nprior)

for (i in 1:sets){
  MClikelihoods[i,] <- apply(priorpvals[[i]],1,likelihood,y = y, n = m)
}

MCestimates4 <- numeric(sets)
for (i in 1:sets) { MCestimates4[i] <- log(1/Nprior*sum(MClikelihoods[i,]))}

MCestimate4 <- mean(MCestimates4)
MCsd4 <- sqrt(var(MCestimates4))

MCestimate4
MCsd4

############################################
### HARMONIC MEAN ESTIMATOR
set.seed(1)

HMlikelihoods <- matrix(0,sets,N)
for (i in 1:sets){
  HMlikelihoods[i,] <- apply(chainpvals[[i]],1,likelihood,y = y, n = m)
}

HMestimates4 <- apply(HMlikelihoods,1,function(x) {-log(1/N*sum(1/x))})

HMestimate4 <- mean(HMestimates4)
HMsd4 <- sqrt(var(HMestimates4))

HMestimate4
HMsd4

############################################
### GENERALISED HARMONIC MEAN ESTIMATOR
set.seed(1)
library(mvtnorm)

GHMpriors <- matrix(0,sets,N)
for(i in 1:sets){
  GHMpriors[i,] <- apply(chains[[i]],1,prior,priormu = priormu, priorsigma = priorsigma)
}

GHMproposaldensities <- matrix(0,sets,N)

for(i in 1:sets){
  mean <- apply(chains[[i]],2,mean)
  var <- var(chains[[i]])
  GHMproposaldensities[i,] <- apply(chains[[i]],1,dmvnorm,mean = mean, sigma = var)
}

matrix <- GHMproposaldensities/(HMlikelihoods*GHMpriors)

GHMestimates4 <- apply(matrix,1,function(x) {-log(1/N*sum(x))})
GHMestimate4 <- mean(GHMestimates4)
GHMsd4 <- sqrt(var(GHMestimates4))

GHMestimate4
GHMsd4
############################################################################################################################
## NEWTON-RAFTERY ESTIMATOR
set.seed(1)

delta = 0.2
NRiterations = 100

indices <- matrix(0,sets,delta*N)
for (i in 1:sets) {
  indices[i,] <- sample(N,delta*N)
}

NRsims <- chains

for (i in 1:sets){
  for(j in indices[i,]){
    NRsims[[i]][j,] <- priorsims[[i]][j,]
  }}

NRpvals <- list()
for(i in 1:sets){
  NRpvals[[i]] <- t(apply(NRsims[[i]],1,function(x) {1/(1+exp(-(X4.1%*%x)))}))
}

NRlikelihoods <- matrix(0,sets,N)
for (i in 1:sets){
  NRlikelihoods[i,] <- apply(NRpvals[[i]],1,likelihood,y = y, n = m)
}

NRestimates4 <- exp(HMestimates4)

for (i in 1:sets){
  for (j in 1:NRiterations) {
    NRestimates4[i] <- (sum(NRlikelihoods[i,]*(delta*NRestimates4[i]+(1-delta)*NRlikelihoods[i,])^-1))/(sum((delta*NRestimates4[i]+(1-delta)*NRlikelihoods[i,])^-1))
  }
}

NRestimates4<- log(NRestimates4)
NRestimate4 <- mean(NRestimates4)
NRsd4 <- sqrt(var(NRestimates4))

NRestimate4
NRsd4

############################################
### LAPLACE-METROPOLIS ESTIMATOR
set.seed(1)

logposteriors <- matrix(0,sets,N)
for ( i in 1:sets) {
  for( j in 1:N) {
    logposteriors[i,j] <- loglikelihood(y,m,chainpvals[[i]][j,])+logprior(chains[[i]][j,],priormu,priorsigma)
  }}

maxlambda <- numeric(sets)
for ( i in 1:sets ) {maxlambda[i] <-  which.max(logposteriors[i,])}

S <- list()
for (i in 1:sets) {S[[i]] <- var(chains[[i]])}

LMestimates4 <- numeric(sets)
for (i in 1:sets) { LMestimates4[i] <- (d/2)*log(2*pi)+0.5*log(det(S[[i]]))+logposteriors[i,][maxlambda[i]]}

LMestimate4 <- mean(LMestimates4)
LMsd4 <- sqrt(var(LMestimates4))

LMestimate4
LMsd4

############################################
### BRIDGE SAMPLING ESTIMATOR
set.seed(1)

BSiterations=1000

N1=N
N2=N

s1=N1/(N1+N2)
s2=N2/(N1+N2)

proposalsims <- list()

for(i in 1:sets){
  mean <- apply(chains[[i]],2,mean)
  var <- var(chains[[i]])
  proposalsims[[i]] <- rmvnorm(N2,mean = mean, sigma = var)
}

proposalpvals <- list()

for(i in 1:sets){
  proposalpvals[[i]] <- t(apply(proposalsims[[i]],1,function(x) {1/(1+exp(-(X4.1%*%x)))}))
}

likelihoods2 <- matrix(0,sets,N2)

for ( i in 1:sets) {
  likelihoods2[i,] <- apply(proposalpvals[[i]],1,likelihood, y = y, n = m)
}

priors2 <- matrix(0,sets,N2)
for ( i in 1:sets) {
  priors2[i,] <- apply(proposalsims[[i]],1,prior, priormu = priormu, priorsigma = priorsigma)
}

proposaldensities2 <- matrix(0,sets,N2)

for(i in 1:sets){
  mean <- apply(chains[[i]],2,mean)
  var <- var(chains[[i]])
  proposaldensities2[i,] <- apply(proposalsims[[i]],1,dmvnorm, mean = mean, sigma = var)
}

likelihoods1 <- HMlikelihoods

priors1 <- GHMpriors

proposaldensities1 <- GHMproposaldensities

x2 <- matrix(0,sets,N2)
for (i in 1:sets){
  x2[i,] <- rep(exp(LMestimate4),N2)
}

x1 <- matrix(0,sets,N1)
for (i in 1:sets){
  x1[i,] <- rep(exp(LMestimate4),N1)
}


for(i in 1:BSiterations){
  matrix2 <- (likelihoods2*priors2)/(s1*likelihoods2*priors2+s2*x2*proposaldensities2)
  matrix1 <- proposaldensities1/(s1*likelihoods1*priors1+s2*x1*proposaldensities1)
  numerator <- apply(matrix2,1,function(x) {1/N2*sum(x)})
  denominator <- apply(matrix1,1,function(x) {1/N1*sum(x)})
  BSestimates4 <- numerator/denominator
  
  x2 <- matrix(0,sets,N2)
  for (i in 1:sets){
    x2[i,] <- rep(exp(BSestimates4[i]),N2)
  }
  
  x1 <- matrix(0,sets,N1)
  for (i in 1:sets){
    x1[i,] <- rep(exp(BSestimates4[i]),N1)
  }
  
}

BSestimates4 <- log(BSestimates4)
BSestimate4 <- mean(BSestimates4)
BSsd4 <- sqrt(var(BSestimates4))

BSestimate4
BSsd4


############################################
#### MODEL 4 - RESULTS
MCestimate4
MCsd4
HMestimate4
HMsd4
GHMestimate4
GHMsd4
NRestimate4
NRsd4
LMestimate4
LMsd4
BSestimate4
BSsd4

##############################################################################
##############################################################################
##### MODEL 5 (BETA = (beta0, beta1, 0, beta3))
##############################################################################
############################################
#### Sampling from the posterior
set.seed(1)

q <- 2    
d <- q+1
library(rstan)

our_model <- stan_model('logistic_regression.stan')
our_data <- list(n = n, q = q, X= X5, y = y, priormu = priormu, priorsigma = priorsigma, m = m)

our_sample <- sampling(our_model, data = our_data, chains = sets, iter = iterations, warmup = warmup, seed = 1)
samples <- as.array(our_sample)

chains <- list() # chains is a list which hold the posterior simulations
for (i in 1:sets){
  chains[[i]] <- samples[,i,-(d+1)]
}

chainpvals <- list() # chainpvals is a list which holds the p values of our posterior simulations
for(i in 1:sets){
  chainpvals[[i]] <- t(apply(chains[[i]],1,function(x) {1/(1+exp(-(X5.1%*%x)))}))
}

#### Sampling from the prior

priorsims <- list() # priorsims is a list which holds our prior simulations
for( i in 1:sets){
  priorsims[[i]] <- matrix(rnorm(Nprior*d,priormu,priorsigma), nrow = Nprior, ncol = d)
}

priorpvals <- list() # priorpvals is a list which holds the p values of our prior simulations
for(i in 1:sets){
  priorpvals[[i]] <- t(apply(priorsims[[i]],1,function(x) {1/(1+exp(-(X5.1%*%x)))}))
}

############################################
#### METHODS
############################################
### MONTE CARLO ESTIMATOR
set.seed(1)

Nprior = N # number of samples from prior

MClikelihoods <- matrix(0,sets,Nprior)

for (i in 1:sets){
  MClikelihoods[i,] <- apply(priorpvals[[i]],1,likelihood,y = y, n = m)
}

MCestimates5 <- numeric(sets)
for (i in 1:sets) { MCestimates5[i] <- log(1/Nprior*sum(MClikelihoods[i,]))}

MCestimate5 <- mean(MCestimates5)
MCsd5 <- sqrt(var(MCestimates5))

MCestimate5
MCsd5

############################################
### HARMONIC MEAN ESTIMATOR
set.seed(1)

HMlikelihoods <- matrix(0,sets,N)
for (i in 1:sets){
  HMlikelihoods[i,] <- apply(chainpvals[[i]],1,likelihood,y = y, n = m)
}

HMestimates5 <- apply(HMlikelihoods,1,function(x) {-log(1/N*sum(1/x))})

HMestimate5 <- mean(HMestimates5)
HMsd5 <- sqrt(var(HMestimates5))

HMestimate5
HMsd5
############################################
### GENERALISED HARMONIC MEAN ESTIMATOR
set.seed(1)
library(mvtnorm)

GHMpriors <- matrix(0,sets,N)
for(i in 1:sets){
  GHMpriors[i,] <- apply(chains[[i]],1,prior,priormu = priormu, priorsigma = priorsigma)
}

GHMproposaldensities <- matrix(0,sets,N)

for(i in 1:sets){
  mean <- apply(chains[[i]],2,mean)
  var <- var(chains[[i]])
  GHMproposaldensities[i,] <- apply(chains[[i]],1,dmvnorm,mean = mean, sigma = var)
}

matrix <- GHMproposaldensities/(HMlikelihoods*GHMpriors)

GHMestimates5 <- apply(matrix,1,function(x) {-log(1/N*sum(x))})
GHMestimate5 <- mean(GHMestimates5)
GHMsd5 <- sqrt(var(GHMestimates5))

GHMestimate5
GHMsd5
############################################################################################################################
## NEWTON-RAFTERY ESTIMATOR
set.seed(1)

delta = 0.2
NRiterations = 100

indices <- matrix(0,sets,delta*N)
for (i in 1:sets) {
  indices[i,] <- sample(N,delta*N)
}

NRsims <- chains

for (i in 1:sets){
  for(j in indices[i,]){
    NRsims[[i]][j,] <- priorsims[[i]][j,]
  }}

NRpvals <- list()
for(i in 1:sets){
  NRpvals[[i]] <- t(apply(NRsims[[i]],1,function(x) {1/(1+exp(-(X5.1%*%x)))}))
}

NRlikelihoods <- matrix(0,sets,N)
for (i in 1:sets){
  NRlikelihoods[i,] <- apply(NRpvals[[i]],1,likelihood,y = y, n = m)
}


NRestimates5 <- exp(HMestimates5)

for (i in 1:sets){
  for (j in 1:NRiterations) {
    NRestimates5[i] <- (sum(NRlikelihoods[i,]*(delta*NRestimates5[i]+(1-delta)*NRlikelihoods[i,])^-1))/(sum((delta*NRestimates5[i]+(1-delta)*NRlikelihoods[i,])^-1))
  }
}

NRestimates5 <- log(NRestimates5)
NRestimate5 <- mean(NRestimates5)
NRsd5 <- sqrt(var(NRestimates5))

NRestimate5
NRsd5

############################################
### LAPLACE-METROPOLIS ESTIMATOR
set.seed(1)

logposteriors <- matrix(0,sets,N)
for ( i in 1:sets) {
  for( j in 1:N) {
    logposteriors[i,j] <- loglikelihood(y,m,chainpvals[[i]][j,])+logprior(chains[[i]][j,],priormu,priorsigma)
  }}

maxlambda <- numeric(sets)
for ( i in 1:sets ) {maxlambda[i] <-  which.max(logposteriors[i,])}

S <- list()
for (i in 1:sets) {S[[i]] <- var(chains[[i]])}

LMestimates5 <- numeric(sets)
for (i in 1:sets) { LMestimates5[i] <- (d/2)*log(2*pi)+0.5*log(det(S[[i]]))+logposteriors[i,][maxlambda[i]]}

LMestimate5 <- mean(LMestimates5)
LMsd5 <- sqrt(var(LMestimates5))

LMestimate5
LMsd5

############################################
### BRIDGE SAMPLING ESTIMATOR
set.seed(1)

BSiterations=1000

N1=N
N2=N

s1=N1/(N1+N2)
s2=N2/(N1+N2)

proposalsims <- list()

for(i in 1:sets){
  mean <- apply(chains[[i]],2,mean)
  var <- var(chains[[i]])
  proposalsims[[i]] <- rmvnorm(N2,mean = mean, sigma = var)
}

proposalpvals <- list()

for(i in 1:sets){
  proposalpvals[[i]] <- t(apply(proposalsims[[i]],1,function(x) {1/(1+exp(-(X5.1%*%x)))}))
}

likelihoods2 <- matrix(0,sets,N2)

for ( i in 1:sets) {
  likelihoods2[i,] <- apply(proposalpvals[[i]],1,likelihood, y = y, n = m)
}

priors2 <- matrix(0,sets,N2)
for ( i in 1:sets) {
  priors2[i,] <- apply(proposalsims[[i]],1,prior, priormu = priormu, priorsigma = priorsigma)
}

proposaldensities2 <- matrix(0,sets,N2)

for(i in 1:sets){
  mean <- apply(chains[[i]],2,mean)
  var <- var(chains[[i]])
  proposaldensities2[i,] <- apply(proposalsims[[i]],1,dmvnorm, mean = mean, sigma = var)
}

likelihoods1 <- HMlikelihoods

priors1 <- GHMpriors

proposaldensities1 <- GHMproposaldensities

x2 <- matrix(0,sets,N2)
for (i in 1:sets){
  x2[i,] <- rep(exp(LMestimate5),N2)
}

x1 <- matrix(0,sets,N1)
for (i in 1:sets){
  x1[i,] <- rep(exp(LMestimate5),N1)
}


for(i in 1:BSiterations){
  matrix2 <- (likelihoods2*priors2)/(s1*likelihoods2*priors2+s2*x2*proposaldensities2)
  matrix1 <- proposaldensities1/(s1*likelihoods1*priors1+s2*x1*proposaldensities1)
  numerator <- apply(matrix2,1,function(x) {1/N2*sum(x)})
  denominator <- apply(matrix1,1,function(x) {1/N1*sum(x)})
  BSestimates5 <- numerator/denominator
  
  x2 <- matrix(0,sets,N2)
  for (i in 1:sets){
    x2[i,] <- rep(exp(BSestimates5[i]),N2)
  }
  
  x1 <- matrix(0,sets,N1)
  for (i in 1:sets){
    x1[i,] <- rep(exp(BSestimates5[i]),N1)
  }
  
}

BSestimates5<- log(BSestimates5)
BSestimate5 <- mean(BSestimates5)
BSsd5 <- sqrt(var(BSestimates5))

BSestimate5
BSsd5


############################################
#### MODEL 5 - RESULTS
MCestimate5
MCsd5
HMestimate5
HMsd5
GHMestimate5
GHMsd5
NRestimate5
NRsd5
LMestimate5
LMsd5
BSestimate5
BSsd5

##############################################################################
##############################################################################
##### MODEL 6 (BETA = (beta0, 0, beta2, beta3))
##############################################################################
############################################
#### Sampling from the posterior
set.seed(1)

q <- 2    
d <- q+1
library(rstan)

our_model <- stan_model('logistic_regression.stan')
our_data <- list(n = n, q = q, X= X6, y = y, priormu = priormu, priorsigma = priorsigma, m = m)

our_sample <- sampling(our_model, data = our_data, chains = sets, iter = iterations, warmup = warmup, seed = 1)
samples <- as.array(our_sample)

chains <- list() # chains is a list which hold the posterior simulations
for (i in 1:sets){
  chains[[i]] <- samples[,i,-(d+1)]
}

chainpvals <- list() # chainpvals is a list which holds the p values of our posterior simulations
for(i in 1:sets){
  chainpvals[[i]] <- t(apply(chains[[i]],1,function(x) {1/(1+exp(-(X6.1%*%x)))}))
}

#### Sampling from the prior

priorsims <- list() # priorsims is a list which holds our prior simulations
for( i in 1:sets){
  priorsims[[i]] <- matrix(rnorm(Nprior*d,priormu,priorsigma), nrow = Nprior, ncol = d)
}

priorpvals <- list() # priorpvals is a list which holds the p values of our prior simulations
for(i in 1:sets){
  priorpvals[[i]] <- t(apply(priorsims[[i]],1,function(x) {1/(1+exp(-(X6.1%*%x)))}))
}

############################################
#### METHODS
############################################
### MONTE CARLO ESTIMATOR
set.seed(1)

Nprior = N # number of samples from prior

MClikelihoods <- matrix(0,sets,Nprior)

for (i in 1:sets){
  MClikelihoods[i,] <- apply(priorpvals[[i]],1,likelihood,y = y, n = m)
}

MCestimates6 <- numeric(sets)
for (i in 1:sets) { MCestimates6[i] <- log(1/Nprior*sum(MClikelihoods[i,]))}

MCestimate6 <- mean(MCestimates6)
MCsd6 <- sqrt(var(MCestimates6))

MCestimate6
MCsd6

############################################
### HARMONIC MEAN ESTIMATOR
set.seed(1)

HMlikelihoods <- matrix(0,sets,N)
for (i in 1:sets){
  HMlikelihoods[i,] <- apply(chainpvals[[i]],1,likelihood,y = y, n = m)
}

HMestimates6 <- apply(HMlikelihoods,1,function(x) {-log(1/N*sum(1/x))})

HMestimate6 <- mean(HMestimates6)
HMsd6 <- sqrt(var(HMestimates6))

HMestimate6
HMsd6

############################################
### GENERALISED HARMONIC MEAN ESTIMATOR
set.seed(1)
library(mvtnorm)

GHMpriors <- matrix(0,sets,N)
for(i in 1:sets){
  GHMpriors[i,] <- apply(chains[[i]],1,prior,priormu = priormu, priorsigma = priorsigma)
}

GHMproposaldensities <- matrix(0,sets,N)

for(i in 1:sets){
  mean <- apply(chains[[i]],2,mean)
  var <- var(chains[[i]])
  GHMproposaldensities[i,] <- apply(chains[[i]],1,dmvnorm,mean = mean, sigma = var)
}

matrix <- GHMproposaldensities/(HMlikelihoods*GHMpriors)

GHMestimates6 <- apply(matrix,1,function(x) {-log(1/N*sum(x))})
GHMestimate6 <- mean(GHMestimates6)
GHMsd6 <- sqrt(var(GHMestimates6))

GHMestimate6
GHMsd6
############################################################################################################################
## NEWTON-RAFTERY ESTIMATOR
set.seed(1)

delta = 0.2
NRiterations = 100

indices <- matrix(0,sets,delta*N)
for (i in 1:sets) {
  indices[i,] <- sample(N,delta*N)
}

NRsims <- chains

for (i in 1:sets){
  for(j in indices[i,]){
    NRsims[[i]][j,] <- priorsims[[i]][j,]
  }}

NRpvals <- list()
for(i in 1:sets){
  NRpvals[[i]] <- t(apply(NRsims[[i]],1,function(x) {1/(1+exp(-(X6.1%*%x)))}))
}

NRlikelihoods <- matrix(0,sets,N)
for (i in 1:sets){
  NRlikelihoods[i,] <- apply(NRpvals[[i]],1,likelihood,y = y, n = m)
}

NRestimates6 <- exp(HMestimates)

for (i in 1:sets){
  for (j in 1:NRiterations) {
    NRestimates6[i] <- (sum(NRlikelihoods[i,]*(delta*NRestimates6[i]+(1-delta)*NRlikelihoods[i,])^-1))/(sum((delta*NRestimates6[i]+(1-delta)*NRlikelihoods[i,])^-1))
  }
}

NRestimates6 <- log(NRestimates6)
NRestimate6 <- mean(NRestimates6)
NRsd6 <- sqrt(var(NRestimates6))

NRestimate6
NRsd6

############################################
### LAPLACE-METROPOLIS ESTIMATOR
set.seed(1)

logposteriors <- matrix(0,sets,N)
for ( i in 1:sets) {
  for( j in 1:N) {
    logposteriors[i,j] <- loglikelihood(y,m,chainpvals[[i]][j,])+logprior(chains[[i]][j,],priormu,priorsigma)
  }}

maxlambda <- numeric(sets)
for ( i in 1:sets ) {maxlambda[i] <-  which.max(logposteriors[i,])}

S <- list()
for (i in 1:sets) {S[[i]] <- var(chains[[i]])}

LMestimates6 <- numeric(sets)
for (i in 1:sets) { LMestimates6[i] <- (d/2)*log(2*pi)+0.5*log(det(S[[i]]))+logposteriors[i,][maxlambda[i]]}

LMestimate6 <- mean(LMestimates6)
LMsd6 <- sqrt(var(LMestimates6))

LMestimate6
LMsd6

############################################
### BRIDGE SAMPLING ESTIMATOR
set.seed(1)

BSiterations=1000

N1=N
N2=N

s1=N1/(N1+N2)
s2=N2/(N1+N2)

proposalsims <- list()

for(i in 1:sets){
  mean <- apply(chains[[i]],2,mean)
  var <- var(chains[[i]])
  proposalsims[[i]] <- rmvnorm(N2,mean = mean, sigma = var)
}

proposalpvals <- list()

for(i in 1:sets){
  proposalpvals[[i]] <- t(apply(proposalsims[[i]],1,function(x) {1/(1+exp(-(X6.1%*%x)))}))
}

likelihoods2 <- matrix(0,sets,N2)

for ( i in 1:sets) {
  likelihoods2[i,] <- apply(proposalpvals[[i]],1,likelihood, y = y, n = m)
}

priors2 <- matrix(0,sets,N2)
for ( i in 1:sets) {
  priors2[i,] <- apply(proposalsims[[i]],1,prior, priormu = priormu, priorsigma = priorsigma)
}

proposaldensities2 <- matrix(0,sets,N2)

for(i in 1:sets){
  mean <- apply(chains[[i]],2,mean)
  var <- var(chains[[i]])
  proposaldensities2[i,] <- apply(proposalsims[[i]],1,dmvnorm, mean = mean, sigma = var)
}

likelihoods1 <- HMlikelihoods

priors1 <- GHMpriors

proposaldensities1 <- GHMproposaldensities

x2 <- matrix(0,sets,N2)
for (i in 1:sets){
  x2[i,] <- rep(exp(LMestimate6),N2)
}

x1 <- matrix(0,sets,N1)
for (i in 1:sets){
  x1[i,] <- rep(exp(LMestimate6),N1)
}


for(i in 1:BSiterations){
  matrix2 <- (likelihoods2*priors2)/(s1*likelihoods2*priors2+s2*x2*proposaldensities2)
  matrix1 <- proposaldensities1/(s1*likelihoods1*priors1+s2*x1*proposaldensities1)
  numerator <- apply(matrix2,1,function(x) {1/N2*sum(x)})
  denominator <- apply(matrix1,1,function(x) {1/N1*sum(x)})
  BSestimates6 <- numerator/denominator
  
  x2 <- matrix(0,sets,N2)
  for (i in 1:sets){
    x2[i,] <- rep(exp(BSestimates6[i]),N2)
  }
  
  x1 <- matrix(0,sets,N1)
  for (i in 1:sets){
    x1[i,] <- rep(exp(BSestimates6[i]),N1)
  }
  
}

BSestimates6 <- log(BSestimates6)
BSestimate6 <- mean(BSestimates6)
BSsd6 <- sqrt(var(BSestimates6))

BSestimate6
BSsd6


############################################
#### MODEL 6 - RESULTS
MCestimate6
MCsd6
HMestimate6
HMsd6
GHMestimate6
GHMsd6
NRestimate6
NRsd6
LMestimate6
LMsd6
BSestimate6
BSsd6

##############################################################################
##############################################################################
##### MODEL 7 (BETA = (beta0, beta1, beta2, beta3))
##############################################################################
############################################
#### Sampling from the posterior
set.seed(1)

q <- 3    
d <- q+1
library(rstan)

our_model <- stan_model('logistic_regression.stan')
our_data <- list(n = n, q = q, X= X7, y = y, priormu = priormu, priorsigma = priorsigma, m = m)

our_sample <- sampling(our_model, data = our_data, chains = sets, iter = iterations, warmup = warmup, seed = 1)
samples <- as.array(our_sample)

chains <- list() # chains is a list which hold the posterior simulations
for (i in 1:sets){
  chains[[i]] <- samples[,i,-(d+1)]
}

chainpvals <- list() # chainpvals is a list which holds the p values of our posterior simulations
for(i in 1:sets){
  chainpvals[[i]] <- t(apply(chains[[i]],1,function(x) {1/(1+exp(-(X7.1%*%x)))}))
}

#### Sampling from the prior

priorsims <- list() # priorsims is a list which holds our prior simulations
for( i in 1:sets){
  priorsims[[i]] <- matrix(rnorm(Nprior*d,priormu,priorsigma), nrow = Nprior, ncol = d)
}

priorpvals <- list() # priorpvals is a list which holds the p values of our prior simulations
for(i in 1:sets){
  priorpvals[[i]] <- t(apply(priorsims[[i]],1,function(x) {1/(1+exp(-(X7.1%*%x)))}))
}

############################################
#### METHODS
############################################
### MONTE CARLO ESTIMATOR
set.seed(1)

Nprior = N # number of samples from prior

MClikelihoods <- matrix(0,sets,Nprior)

for (i in 1:sets){
  MClikelihoods[i,] <- apply(priorpvals[[i]],1,likelihood,y = y, n = m)
}

MCestimates7 <- numeric(sets)
for (i in 1:sets) { MCestimates7[i] <- log(1/Nprior*sum(MClikelihoods[i,]))}

MCestimate7 <- mean(MCestimates7)
MCsd7 <- sqrt(var(MCestimates7))

MCestimate7
MCsd7

############################################
### HARMONIC MEAN ESTIMATOR
set.seed(1)

HMlikelihoods <- matrix(0,sets,N)
for (i in 1:sets){
  HMlikelihoods[i,] <- apply(chainpvals[[i]],1,likelihood,y = y, n = m)
}

HMestimates7 <- apply(HMlikelihoods,1,function(x) {-log(1/N*sum(1/x))})

HMestimate7 <- mean(HMestimates7)
HMsd7 <- sqrt(var(HMestimates7))

HMestimate7
HMsd7

############################################
### GENERALISED HARMONIC MEAN ESTIMATOR
set.seed(1)
library(mvtnorm)

GHMpriors <- matrix(0,sets,N)
for(i in 1:sets){
  GHMpriors[i,] <- apply(chains[[i]],1,prior,priormu = priormu, priorsigma = priorsigma)
}

GHMproposaldensities <- matrix(0,sets,N)

for(i in 1:sets){
  mean <- apply(chains[[i]],2,mean)
  var <- var(chains[[i]])
  GHMproposaldensities[i,] <- apply(chains[[i]],1,dmvnorm,mean = mean, sigma = var)
}

matrix <- GHMproposaldensities/(HMlikelihoods*GHMpriors)

GHMestimates7 <- apply(matrix,1,function(x) {-log(1/N*sum(x))})
GHMestimate7 <- mean(GHMestimates7)
GHMsd7 <- sqrt(var(GHMestimates7))

GHMestimate7
GHMsd7
############################################################################################################################
## NEWTON-RAFTERY ESTIMATOR
set.seed(1)

delta = 0.2
NRiterations = 100

indices <- matrix(0,sets,delta*N)
for (i in 1:sets) {
  indices[i,] <- sample(N,delta*N)
}

NRsims <- chains

for (i in 1:sets){
  for(j in indices[i,]){
    NRsims[[i]][j,] <- priorsims[[i]][j,]
  }}

NRpvals <- list()
for(i in 1:sets){
  NRpvals[[i]] <- t(apply(NRsims[[i]],1,function(x) {1/(1+exp(-(X7.1%*%x)))}))
}

NRlikelihoods <- matrix(0,sets,N)
for (i in 1:sets){
  NRlikelihoods[i,] <- apply(NRpvals[[i]],1,likelihood,y = y, n = m)
}


NRestimates7 <- exp(HMestimates7)

for (i in 1:sets){
  for (j in 1:NRiterations) {
    NRestimates7[i] <- (sum(NRlikelihoods[i,]*(delta*NRestimates7[i]+(1-delta)*NRlikelihoods[i,])^-1))/(sum((delta*NRestimates7[i]+(1-delta)*NRlikelihoods[i,])^-1))
  }
}

NRestimates7 <- log(NRestimates7)
NRestimate7 <- mean(NRestimates7)
NRsd7 <- sqrt(var(NRestimates7))

NRestimate7
NRsd7

############################################
### LAPLACE-METROPOLIS ESTIMATOR
set.seed(1)

logposteriors <- matrix(0,sets,N)
for ( i in 1:sets) {
  for( j in 1:N) {
    logposteriors[i,j] <- loglikelihood(y,m,chainpvals[[i]][j,])+logprior(chains[[i]][j,],priormu,priorsigma)
  }}

maxlambda <- numeric(sets)
for ( i in 1:sets ) {maxlambda[i] <-  which.max(logposteriors[i,])}

S <- list()
for (i in 1:sets) {S[[i]] <- var(chains[[i]])}

LMestimates7 <- numeric(sets)
for (i in 1:sets) { LMestimates7[i] <- (d/2)*log(2*pi)+0.5*log(det(S[[i]]))+logposteriors[i,][maxlambda[i]]}

LMestimate7 <- mean(LMestimates7)
LMsd7 <- sqrt(var(LMestimates7))

LMestimate7
LMsd7

############################################
### BRIDGE SAMPLING ESTIMATOR
set.seed(1)

BSiterations=1000

N1=N
N2=N

s1=N1/(N1+N2)
s2=N2/(N1+N2)

proposalsims <- list()

for(i in 1:sets){
  mean <- apply(chains[[i]],2,mean)
  var <- var(chains[[i]])
  proposalsims[[i]] <- rmvnorm(N2,mean = mean, sigma = var)
}

proposalpvals <- list()

for(i in 1:sets){
  proposalpvals[[i]] <- t(apply(proposalsims[[i]],1,function(x) {1/(1+exp(-(X7.1%*%x)))}))
}

likelihoods2 <- matrix(0,sets,N2)

for ( i in 1:sets) {
  likelihoods2[i,] <- apply(proposalpvals[[i]],1,likelihood, y = y, n = m)
}

priors2 <- matrix(0,sets,N2)
for ( i in 1:sets) {
  priors2[i,] <- apply(proposalsims[[i]],1,prior, priormu = priormu, priorsigma = priorsigma)
}

proposaldensities2 <- matrix(0,sets,N2)

for(i in 1:sets){
  mean <- apply(chains[[i]],2,mean)
  var <- var(chains[[i]])
  proposaldensities2[i,] <- apply(proposalsims[[i]],1,dmvnorm, mean = mean, sigma = var)
}

likelihoods1 <- HMlikelihoods

priors1 <- GHMpriors

proposaldensities1 <- GHMproposaldensities

x2 <- matrix(0,sets,N2)
for (i in 1:sets){
  x2[i,] <- rep(exp(LMestimate7),N2)
}

x1 <- matrix(0,sets,N1)
for (i in 1:sets){
  x1[i,] <- rep(exp(LMestimate7),N1)
}


for(i in 1:BSiterations){
  matrix2 <- (likelihoods2*priors2)/(s1*likelihoods2*priors2+s2*x2*proposaldensities2)
  matrix1 <- proposaldensities1/(s1*likelihoods1*priors1+s2*x1*proposaldensities1)
  numerator <- apply(matrix2,1,function(x) {1/N2*sum(x)})
  denominator <- apply(matrix1,1,function(x) {1/N1*sum(x)})
  BSestimates7 <- numerator/denominator
  
  x2 <- matrix(0,sets,N2)
  for (i in 1:sets){
    x2[i,] <- rep(exp(BSestimates7[i]),N2)
  }
  
  x1 <- matrix(0,sets,N1)
  for (i in 1:sets){
    x1[i,] <- rep(exp(BSestimates7[i]),N1)
  }
  
}

BSestimates7 <- log(BSestimates7)
BSestimate7 <- mean(BSestimates7)
BSsd7 <- sqrt(var(BSestimates7))

BSestimate7
BSsd7


############################################
#### MODEL 7 - RESULTS
MCestimate7
MCsd7
HMestimate7
HMsd7
GHMestimate7
GHMsd7
NRestimate7
NRsd7
LMestimate7
LMsd7
BSestimate7
BSsd7

############################################
#### COMPLETE RESULTS
############################################
df <- data.frame(
  'Model' = c("M1", "-", "M2", "-","M3", "-","M4", "-","M5", "-","M6", "-","M7", "-"),
  'MC' = c(MCestimate1,MCsd1,MCestimate2,MCsd2,MCestimate3,MCsd3,MCestimate4,MCsd4,MCestimate5,MCsd5,MCestimate6,MCsd6,MCestimate7,MCsd7),
  'HM' = c(HMestimate1,HMsd1,HMestimate2,HMsd2,HMestimate3,HMsd3,HMestimate4,HMsd4,HMestimate5,HMsd5,HMestimate6,HMsd6,HMestimate7,HMsd7),
  'GHM' = c(GHMestimate1,GHMsd1,GHMestimate2,GHMsd2,GHMestimate3,GHMsd3,GHMestimate4,GHMsd4,GHMestimate5,GHMsd5,GHMestimate6,GHMsd6,GHMestimate7,GHMsd7),
  'NR' = c(NRestimate1,NRsd1,NRestimate2,NRsd2,NRestimate3,NRsd3,NRestimate4,NRsd4,NRestimate5,NRsd5,NRestimate6,NRsd6, NRestimate7, NRsd7),
  'LM' = c(LMestimate1,LMsd1,LMestimate2,LMsd2,LMestimate3,LMsd3,LMestimate4,LMsd4,LMestimate5,LMsd5,LMestimate6,LMsd6,LMestimate7,LMsd7),
  'BS' = c(BSestimate1,BSsd1,BSestimate2,BSsd2,BSestimate3,BSsd3,BSestimate4,BSsd4,BSestimate5,BSsd5,BSestimate6,BSsd6,BSestimate7,BSsd7)
)

knitr::kable(df)

####################################################################################
#### Posterior Sample plots
####################################################################################

par(mfrow=c(1,2))
## posterior density of betas
for (i in 1:d){
  name <- paste(expression(beta), i, sep = "")
  plot(density(samples[,1:3,i]),main = name)
}

library(bayesplot)
library(cowplot)
library(ggplot2)

####################################################################################################

color_scheme_set('viridis')

mcmc_trace(samples[1:500,1:3,],  regex_pars = "beta",
           facet_args = list(labeller = ggplot2::label_parsed)) +
  legend_move('top')


trace_plot1 <- mcmc_trace(samples[1:500,1:3,],  regex_pars = "beta\\[[1,2]\\]",
                          facet_args = list(labeller = ggplot2::label_parsed)) 

trace_plot2 <- mcmc_trace(samples[1:500,1:3,],  regex_pars = "beta\\[[3,4]\\]",
                          facet_args = list(labeller = ggplot2::label_parsed))

plot_grid(trace_plot1, trace_plot2, nrow = 2)

########################################################################################################################
color_scheme_set("brightblue")

density_plot1 <- mcmc_areas(samples[,1:3,],  regex_pars = "beta\\[[1]\\]",
)
density_plot2 <- mcmc_areas(samples[,1:3,],  regex_pars = "beta\\[[2]\\]", ylim = c(0,1)
)
density_plot3 <- mcmc_areas(samples[,1:3,],  regex_pars = "beta\\[[3]\\]",
)
density_plot4 <- mcmc_areas(samples[,1:3,],  regex_pars = "beta\\[[4]\\]",
)

plot_grid(density_plot1,density_plot2,density_plot3,density_plot4, nrow = 1)

################################################################################################

color_scheme_set("brightblue")

plot <- mcmc_areas_ridges(samples[,1:3,], regex_pars = "beta", border_size = 0.75)

plot + xlim(c(-0.7,0.4)) + theme(plot.margin = margin(t = -100, r = 0, b = 5, l = 5))
