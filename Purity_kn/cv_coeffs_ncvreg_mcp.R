cv_coeffs_ncvreg_mcp<-function (X, y, nlambda = 500, intercept = T, parallel = T, 
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
      # lambda_max = 10^-2*max(abs(t(X) %*% y))/n
      lambda_max = max(abs(t(X) %*% y))/n
      lambda_min = lambda_max/2000
      # lambda_min = 10^-18*lambda_max/2000
      k = (0:(nlambda - 1))/nlambda
      #skcm 3: 1 gamma=3.5 
      #skcm 4 gamma=10 gamma=15/7/59
      #skcm 5 gamma=20 40/6/160
      #LIHC 5 gamma=20 5/76
      lambda = lambda_max * (lambda_min/lambda_max)^k
    }
    else {
      lambda = NULL
    }
  }
  no_cores <- detectCores(logical = T)
  cl <- makeCluster(no_cores)
  registerDoParallel(cl, cores=no_cores)
  # cv.ncvreg.fit <- ncvreg::cv.ncvreg(X, y,cluster = cl,penalty="MCP",gamma=4,family="binomial",lambda = lambda)
  # cv.ncvreg.fit <- ncvreg::cv.ncvreg(X, y,cluster = cl,penalty="MCP",gamma=20,lambda = lambda,max.iter=100000)
  cv.ncvreg.fit <- ncvreg::cv.ncvreg(X, y,cluster = cl,penalty="MCP",gamma=30,lambda = lambda,max.iter=100000)
  stopImplicitCluster()
  stopCluster(cl)
  as.vector(coef(cv.ncvreg.fit, s = "lambda.min")[2:(p + 1)])
}
