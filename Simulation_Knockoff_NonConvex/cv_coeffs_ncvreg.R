cv_coeffs_ncvreg<-function (X, y, nlambda = 500, intercept = T, parallel = T, 
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
      lambda = lambda_max * (lambda_min/lambda_max)^k
    }
    else {
      lambda = NULL
    }
  }
  no_cores <- detectCores(logical = T)
  cl <- makeCluster(no_cores)
  registerDoParallel(cl, cores=no_cores)
  # penalty: SCAD OR MCP (default LASSO) family: gaussion(default), or binomail
  # cv.ncvreg.fit <- ncvreg::cv.ncvreg(X, y,cluster = cl, penalty="SCAD",family="binomial",lambda = lambda)
  cv.ncvreg.fit <- ncvreg::cv.ncvreg(X, y,cluster = cl, penalty="SCAD",lambda = lambda)
  stopImplicitCluster()
  stopCluster(cl)
  as.vector(coef(cv.ncvreg.fit, s = "lambda.min")[2:(p + 1)])
}
