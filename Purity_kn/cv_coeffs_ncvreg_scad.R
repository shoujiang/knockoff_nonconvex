cv_coeffs_ncvreg_scad<-function (X, y, nlambda = 500, intercept = T, parallel = T, 
          ...) 
{
  X = scale(X)
  n = nrow(X)
  p = ncol(X)
  if (!methods::hasArg(family)) 
    family = "gaussian"
  else family = list(...)$family
  if (!methods::hasArg(lambda)) {
    if (identical(family, "gaussian")) {
      if (!is.numeric(y)) {
        stop("Input y must be numeric.")
      }
      lambda_max = max(abs(t(X) %*% y))/n
      lambda_min = lambda_max/2000
      k = (0:(nlambda - 1))/nlambda
      #BRCA 3: 0.02 0.03 0.019 0.045 0.049 0.07 0.085
      #BRCA 4: 0.05 0.08
      #SKCM 4: 1 1.2 1.125 gamma=15/3/50
      #SKCM 5: 20/5/143
      lambda = lambda_max * (lambda_min/lambda_max)^k
    }
    else {
      lambda = NULL
    }
  }
  no_cores <- detectCores(logical = T)
  cl <- makeCluster(no_cores)
  registerDoParallel(cl, cores=no_cores)
  # cv.ncvreg.fit <- ncvreg::cv.ncvreg(X, y,cluster = cl, penalty="SCAD",family="binomial",lambda = lambda)
  cv.ncvreg.fit <- ncvreg::cv.ncvreg(X, y,cluster = cl, penalty="SCAD",gamma=30,lambda = lambda,max.iter=100000)
  stopImplicitCluster()
  stopCluster(cl)
  as.vector(coef(cv.ncvreg.fit, s = "lambda.min")[2:(p + 1)])
}
