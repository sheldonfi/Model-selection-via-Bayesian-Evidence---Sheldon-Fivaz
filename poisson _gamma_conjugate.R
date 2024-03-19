## working with a Gamma - Poisson conjugate model
## We set our parameters lambda, a and b for Poisson(lambda) and Gamma(a,b)

#######################################################
### Configurations
sets=50   # the number of sets of data
n = 50    # the number of samples in each set
N=10000
lambda = 2
a = 0.001
b = 0.001
########################################################
### Prior vs Posterior density plots
library(ggplot2)

x <- seq(0,5,0.001)
y <- dgamma(x,a,b)

data <- data.frame(x = x, y = y)
ggplot(data = data, aes(x = x, y = y)) +
  geom_line(color="#69b3a2", size = 1) +
  labs(x = "Values", y = "Density") +
  scale_y_continuous(limits = c(0, 0.02))

plot(x,y, type = 'l', col="#69b3a2", xlab = 'Lambda values', ylab = 'Density')
########################################################
## Useful Functions

logmarglikelihood <- function(data,ybar,n,a,b) {a*log(b)+log(gamma(n*ybar+a))-sum(log(factorial(data)))-log(gamma(a))-(n*ybar+a)*log(n+b)}

likelihood <- function(x,data) {prod(dpois(data,x, log = FALSE))}
loglikelihood <- function(x,data) {sum(dpois(data,x, log = TRUE))}

prior <- function(x) {dgamma(x,a,b, log = FALSE)}
logprior <- function(x) {dgamma(x,a,b, log = TRUE)}

########################################################
### Data simulations
set.seed(1)

sims <- matrix(0,sets,n)

for (i in 1:sets) { sims[i,] <- rpois(n,lambda)}

ybar <- apply(sims,1,mean)

results <- numeric(sets) # holds our true model evidence values
for (i in 1:sets) { results[i] <- logmarglikelihood(sims[i,],ybar[i],n,a,b) }

postsims <- t(apply(sims,1,function(x) {rgamma(N, a + sum(x), b+n)})) # 50 posterior samples

priorsims <- matrix(0,sets,N)
for (i in 1:sets)
priorsims[i,] <- rgamma(N,a,b)
priorsims[priorsims==0] <- .Machine$double.xmin
## Now we use  computational methods to estimate the evidence and compare this to the likelihood
############################################################################################################################
## Monte Carlo estimator
set.seed(1)

likelihoods <- matrix(0,sets,N)

for(i in 1:sets){
likelihoods[i,] <- sapply(priorsims[i,],likelihood,data = sims[i,])
}

MCestimates <- apply(likelihoods,1, function(x) {log(1/N*sum(x))})

MCerrors <- abs(MCestimates-results)
boxplot(MCerrors)
############################################################################################################################   
## Harmonic Mean estimate
set.seed(1)

likelihoods <- matrix(0,sets,N)
for (i in 1:sets){
  likelihoods[i,] <- sapply(postsims[i,],likelihood,data = sims[i,])
}

HMestimates <- apply(likelihoods,1,function(x) {-log(1/N*sum(1/x))})

HMerrors <- abs(HMestimates-results)

boxplot(HMerrors)
############################################################################################################################
## Generalised Harmonic Mean Estimate
set.seed(1)

proposaldensities <- matrix(0,sets,N)
for (i in 1:sets){
  #mean = (a+sum(sims[i,]))/(b+n)
  mean = mean(postsims[i,])
  var = var(postsims[i,])
  proposaldensities[i,] <- sapply(postsims[i,],dnorm,mean = mean, sd = sqrt(var))
}
priors <- t(apply(postsims,1,dgamma, shape = a, rate = b))

matrix <- proposaldensities/(likelihoods*priors)
HMestimates1 <- apply(matrix,1,function(x) {-log(1/N*sum(x))})

HMerrors1 <- abs(HMestimates1-results)

boxplot(HMerrors1)
############################################################################################################################
## Newton- Raftery estimator
set.seed(1)

delta = 0.2
iterations = 1000

indices <- matrix(0,sets,delta*N)
for (i in 1:sets) {
  indices[i,] <- sample(N,delta*N)
}

NRsims <- postsims

for (i in 1:sets){
  for(j in indices[i,]){
    NRsims[i,j] <- priorsims[i,j]
  }}


NRlikelihoods <- matrix(0,sets,N)
for (i in 1:sets){
  NRlikelihoods[i,] <- sapply(NRsims[i,],likelihood,data = sims[i,])
}

NRestimates <- exp(HMestimates1)
sequence <- matrix(0,sets,iterations)
for(i in 1:sets){
  sequence[i,1] <- NRestimates[i]
}

for (i in 1:sets){
  for (j in 1:iterations) {
    NRestimates[i] <- (sum(NRlikelihoods[i,]*(delta*NRestimates[i]+(1-delta)*NRlikelihoods[i,])^-1))/(sum((delta*NRestimates[i]+(1-delta)*NRlikelihoods[i,])^-1))
  sequence[i,j] <- NRestimates[i]
    }
}

NRestimates <- log(NRestimates)

NRerrors <- abs(NRestimates-results)
plot(1:iterations,sequence[1,])
boxplot(NRerrors,MCerrors)
############################################################################################################################
## Laplace-Metropolis estimator
set.seed(1)

logposterior <- matrix(0,sets,N)
for ( i in 1:sets) {
  for( j in 1:N) {
  logposterior[i,j] <- loglikelihood(postsims[i,j],sims[i,])+logprior(postsims[i,j])
}}

maxlambda <- numeric(sets)
for ( i in 1:sets ) {maxlambda[i] <-  which.max(logposterior[i,])}

S <- numeric(sets)
for (i in 1:sets) {S[i] <- var(postsims[i,])}

d <- 1

LMestimates <- numeric(sets)
for (i in 1:sets) { LMestimates[i] <- (d/2)*log(2*pi)+0.5*log(S[i])+logposterior[i,][maxlambda[i]]}

LMestimates

LMerrors <- abs(LMestimates-results)

boxplot(LMerrors)

############################################################################################################################
## Bridge sampling estimator
set.seed(1)

iterations=5

N1=N
N2=N

s1=N1/(N1+N2)
s2=N2/(N1+N2)

postsims <- t(apply(sims,1,function(x) {rgamma(N1, a + sum(x), b+n)}))

proposalsims <- matrix(0,sets,N2)
for (i in 1:sets){
  mean = (a+sum(sims[i,]))/(b+n)
  sd=sqrt((a+sum(sims[i,]))/((b+n)^2))
  proposalsims[i,] <- (rnorm(N,mean = mean, sd = sd))
}

likelihoods2 <- matrix(0,sets,N2)

#for (i in 1:sets){
#  likelihoods2[i,] <- sapply(proposalsims[i,],likelihood,data = sims[i,])
#}

for ( i in 1:sets) {
  for( j in 1:N2) {
    likelihoods2[i,j] <- likelihood(proposalsims[i,j],sims[i,])
  }}


priors2 <- t(apply(proposalsims,1,dgamma, shape = a, rate = b))

proposaldensities2 <- matrix(0,sets,N2)
for (i in 1:sets){
  mean = (a+sum(sims[i,]))/(b+n)
  sd=sqrt((a+sum(sims[i,]))/((b+n)^2))
  proposaldensities2[i,] <- sapply(proposalsims[i,],dnorm,mean = mean, sd = sd)
}

likelihoods1 <- matrix(0,sets,N1)

#for (i in 1:sets){
#  likelihoods1[i,] <- sapply(postsims[i,],likelihood,data = sims[i,])
#}

for ( i in 1:sets) {
  for( j in 1:N1) {
    likelihoods1[i,j] <- likelihood(postsims[i,j],sims[i,])
  }}


priors1 <- t(apply(postsims,1,dgamma, shape = a, rate = b))

proposaldensities1 <- matrix(0,sets,N1)
for (i in 1:sets){
  mean = (a+sum(sims[i,]))/(b+n)
  sd=sqrt((a+sum(sims[i,]))/((b+n)^2))
  proposaldensities1[i,] <- sapply(postsims[i,],dnorm,mean = mean, sd = sd)
}

x2 <- matrix(0,sets,N2)
for (i in 1:sets){
  x2[i,] <- rep(exp(LMestimates[i]),N2)
}

x1 <- matrix(0,sets,N1)
for (i in 1:sets){
  x1[i,] <- rep(exp(LMestimates[i]),N1)
}

  
for(i in 1:iterations){
matrix2 <- (likelihoods2*priors2)/(s1*likelihoods2*priors2+s2*x2*proposaldensities2)
matrix1 <- proposaldensities1/(s1*likelihoods1*priors1+s2*x1*proposaldensities1)
numerator <- apply(matrix2,1,function(x) {1/N2*sum(x)})
denominator <- apply(matrix1,1,function(x) {1/N1*sum(x)})
BSestimates <- numerator/denominator

x2 <- matrix(0,sets,N2)
for (i in 1:sets){
  x2[i,] <- rep(exp(BSestimates[i]),N2)
}

x1 <- matrix(0,sets,N1)
for (i in 1:sets){
  x1[i,] <- rep(exp(BSestimates[i]),N1)
}

}

BSestimates <- log(BSestimates)
BSerrors <- abs(BSestimates-results)
boxplot(BSerrors)
############################################################################################################################
### Fourier Integral Method
set.seed(1)

N <- N
R <- 20


postsims <- t(apply(sims,1,function(x) {rgamma(N, a + sum(x), b+n)}))

postmeans <- apply(postsims,1,mean)

FIestimates <- numeric(sets)
postestimates <- numeric(sets)

for(i in 1:sets){
  postestimates[i] <- 1/(N*pi)*sum(sin(R*(postmeans[i]-postsims[i,]))/(postmeans[i]-postsims[i,]))
}

##for (i in seq(1,100)){
##print(1/(N*pi)*sum(sin(i*(postmeans[1]-postsims[1,]))/(postmeans[1]-postsims[1,])))
##}

dgamma(postmeans[1],a+sum(sims[1,]),b+n)

for( i in 1:sets){
FIestimates[i] <- loglikelihood(postmeans[i],sims[i,])+logprior(postmeans[i])-log(postestimates[i])
}

FIerrors <- abs(FIestimates-results)
boxplot(abs(FIestimates-results))

###################################################################################################  
###### MC, GHM AND LM summaries
par(mar = c(2.2,4.3,2,0.4))
boxplot(MCerrors,HMerrors1, LMerrors, xlab = '', ylab = 'Absolute Error', names = c('Monte Carlo','Gen. Harmonic Mean','Laplace-Metropolis'), outine.col = 'red', pch = 20 , 
        boxfill = "lightblue", whiskcol = 'blue', staplecol = 'blue', border = rgb(0.2,0.4,1), cex.lab = 1.4, cex.main = 1.3)
points(jitter(rep(1, length(MCerrors)), factor = 0.2), MCerrors, col = "yellow", pch = 4, bg = 'transparent')
points(jitter(rep(2, length(HMerrors1)), factor = 0.2), HMerrors1, col = "yellow", pch = 4, bg = 'transparent')
points(jitter(rep(3, length(LMerrors)), factor = 0.2), LMerrors, col = "yellow", pch = 4, bg = 'transparent')

############################################################################################################################
## box plot of errors
boxplot(LMerrors,HMerrors,HMerrors1,NRerrors,MCerrors,BSerrors,FIerrors)
