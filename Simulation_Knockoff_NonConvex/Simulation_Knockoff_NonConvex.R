library(foreach)
library(doParallel)
cores <- detectCores(logical=T)
cl <- makeCluster(cores)
registerDoParallel(cl, cores=cores)
library(knockoff) #install.packages("knockoff")
library(parallel)#install.packages("parallel")
library(ncvreg)#install.packages("ncvreg")
library(doMC)#install.packages("doMC", repos="http://R-Forge.R-project.org")
n = 300          # number of observations
p = 100         # number of variables
k = 60            # number of variables with nonzero coefficients
trials=30
setwd('D:/OneDrive - Macau University of Science and Technology/KnockOff/Code1/MCP/Knockoff MCP')
getwd()
source('stat.nonconvex.R')
source('cv_coeffs_ncvreg.R') # SCAD OR MCP penalty
# Parallel to amplitude, OR correlation rho
Res <- foreach(amplitude=c(2,2.5,3,3.5,4,4.5,5), .combine='rbind',.packages = c("knockoff","doParallel")) %dopar%
  {
  set.seed(123)
  mu = rep(0,p)
  rho = 0
  Sigma = toeplitz(rho^(0:(p-1))) 
  X = matrix(rnorm(n*p),n) %*% chol(Sigma)#Cholesky 
  nonzero = sample(p, k)
  beta = amplitude * (1:p %in% nonzero) / sqrt(n)#%in% 
  y.sample = function(X) X %*% beta + rnorm(n) #  y default gaussion linear model 
  # pr = 1/(1+exp(-X %*% beta)) # y ~ Binomial linear model
  # y = rbinom(n,1,pr)
  y = y.sample(X)
  allPower<-c()
  allFDR<-c()
  fdp = function(selected) sum(beta[selected] == 0) / max(1, length(selected))
  power <- function(selected, beta) {sum(beta[selected] != 0) / sum(beta != 0)}
  for (i in 1:trials) {
    result = knockoff.filter(X, y, offset = 0, statistic = stat.nonconvex)
    allFDR<-c(allFDR,fdp(result$selected)) 
    allPower<-c(allPower,power(result$selected,beta))
  }
  output<-rbind(allFDR,allPower)
  return(output)
  }
write.csv(Res,'result_KN_MCP.csv')
stopImplicitCluster()
stopCluster(cl)

