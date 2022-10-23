rm(list = ls(all.names = TRUE))
gc()
setwd("D:\\")
library(foreach)
library(doParallel)
cores <- detectCores(logical=T)
cl <- makeCluster(cores)
registerDoParallel(cl, cores=cores)
library(knockoff)
library(parallel)
library(ncvreg)
library("KOBT")

df_SKCM = read.table('df_SKCM.txt',header = T)
data_SKCM = na.omit(df_SKCM)
dataX = data_SKCM[,-c(1:4)]
# dataX = dataX[,colSums(dataX)>0]
dataX = dataX[,apply(dataX, 2, var)!=0]
datay = data_SKCM$purity

X=data.matrix(dataX)
y=datay
cv_nlambda<-function (X, y, nlambda = 500)
{
  X = scale(X)
  n = nrow(X)
  p = ncol(X)
  lambda_max = max(abs(t(X) %*% y))/n
  # lambda_min = 10^-15*lambda_max/2000
  lambda_min = lambda_max/2000
  k = (0:(nlambda - 1))/nlambda
  lambda = lambda_max * (lambda_min/lambda_max)^k
}
nlambda = cv_nlambda(X, y, nlambda = 100)#  
cvob1 <- glmnet::cv.glmnet(X, y,lambda = nlambda)
tmpX <- reduce.dim(fit = cvob1, X = X, bound = 1000)# bound = 500 is well worked
dim(tmpX$sub.X)
X = tmpX$sub.X
dim(X)

###################################################################
t1=proc.time()
# knockoff mcp
source('stat.mcp.R')
source('cv_coeffs_ncvreg_mcp.R')
res_mcp<-c()
X <- scale(X)
SEED = 2022
set.seed(SEED)
trials = 100
for (i in 1:trials) {
  result_mcp = knockoff.filter(X, y, statistic = stat.mcp,fdr=0.1)
  res_mcp_sed = result_mcp$selected
  res_mcp = c(res_mcp,res_mcp_sed)
}
t2=proc.time()
t=t2-t1
print(paste0("running time:",round(t[3][[1]]/60,2),' minuts'))
stopImplicitCluster()
stopCluster(cl)
table(names(res_mcp))
selected_gene_mcp<-table(names(res_mcp))
selected_gene<-selected_gene_mcp[selected_gene_mcp>50]
sort(selected_gene_mcp)
max(selected_gene_mcp)
length(selected_gene)
selected_gene
# save(res_mcp,file = 'ResMCPSKCM.RData')


#######################################################################

#knockoff lasso
# res_lasso<-c()
# for (i in 1:trials) {
#   X = scale(X)
#   my_stat = function(...) stat.glmnet_coefdiff(...,cores = cores)
#   result_LASSO = knockoff.filter(X, y, fdr = 0.1,statistic = my_stat)
#   res_lasso_sed = result_LASSO$selected
#   res_lasso = c(res_lasso,res_lasso_sed)
# }

# rlasso <- result_LASSO$selected
# rlasso1 <- result_LASSO$selected
# rlasso = c(rlasso,rlasso1)
# table(names(res_lasso))
