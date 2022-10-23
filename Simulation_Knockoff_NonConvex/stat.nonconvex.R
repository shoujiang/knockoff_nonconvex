stat.nonconvex<-function (X, X_k,y,...) 
{
  if (!requireNamespace("ncvreg", quietly = T)) 
    stop("ncvreg is not installed", call. = F)
  swap = rbinom(ncol(X), 1, 0.5)
  swap.M = matrix(swap, nrow = nrow(X), ncol = length(swap), 
                  byrow = TRUE)
  X.swap = X * (1 - swap.M) + X_k * swap.M
  Xk.swap = X * swap.M + X_k * (1 - swap.M)
  Z = cv_coeffs_ncvreg(cbind(X.swap, Xk.swap), y)
  p = ncol(X)
  orig = 1:p
  W = abs(Z[orig]) - abs(Z[orig + p])
  W = W * (1 - 2 * swap)
}
