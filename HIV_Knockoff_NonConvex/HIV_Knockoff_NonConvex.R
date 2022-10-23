drug_class = 'NNRTI' # Possible drug types are 'PI', 'NRTI', and 'NNRTI'. 
base_url = 'http://hivdb.stanford.edu/pages/published_analysis/genophenoPNAS2006'
gene_url = paste(base_url, 'DATA', paste0(drug_class, '_DATA.txt'), sep='/')
tsm_url = paste(base_url, 'MUTATIONLISTS', 'NP_TSM', drug_class, sep='/')

gene_df = read.delim(gene_url, na.string = c('NA', ''), stringsAsFactors = FALSE)
tsm_df = read.delim(tsm_url, header = FALSE, stringsAsFactors = FALSE)
names(tsm_df) = c('Position', 'Mutations')
# Returns rows for which every column matches the given regular expression.
grepl_rows <- function(pattern, df) {
  cell_matches = apply(df, c(1,2), function(x) grepl(pattern, x))
  apply(cell_matches, 1, all)
}

pos_start = which(names(gene_df) == 'P1')
pos_cols = seq.int(pos_start, ncol(gene_df))
valid_rows = grepl_rows('^(\\.|-|[A-Zid]+)$', gene_df[,pos_cols])
gene_df = gene_df[valid_rows,]


# Flatten a matrix to a vector with names from concatenating row/column names.
flatten_matrix <- function(M, sep='.') {
  x <- c(M)
  names(x) <- c(outer(rownames(M), colnames(M),
                      function(...) paste(..., sep=sep)))
  x
}

# Construct preliminary design matrix.
muts = c(LETTERS, 'i', 'd')
X = outer(muts, as.matrix(gene_df[,pos_cols]), Vectorize(grepl))
X = aperm(X, c(2,3,1))
dimnames(X)[[3]] <- muts
X = t(apply(X, 1, flatten_matrix))
mode(X) <- 'numeric'

# Remove any mutation/position pairs that never appear in the data.
X = X[,colSums(X) != 0]

# Extract response matrix.
Y = gene_df[,4:(pos_start-1)]

library(foreach)
library(doParallel)
cores <- detectCores(logical=T)
cl <- makeCluster(cores)
registerDoParallel(cl, cores=cores)
library(knockoff)
library(parallel)
library(ncvreg)
library(doMC)
setwd('E:/OneDrive - Macau University of Science and Technology/Basic/02KnockOff/Code1/HIV')
getwd()
source('stat.SCAD.R')
source('cv_coeffs_ncvreg_SCAD.R')
source('stat.MCP.R')
source('cv_coeffs_ncvreg_MCP.R')
myknockoff <- function (X, y, q) {
  # Log-transform the drug resistance measurements.
  y = log(y)
  
  # Remove patients with missing measurements.
  missing = is.na(y)
  y = y[!missing]
  X = X[!missing,]
  
  # Remove predictors that appear less than 3 times.
  X = X[,colSums(X) >= 3]
  
  # Remove duplicate predictors.
  X = X[,colSums(abs(cor(X)-1) < 1e-4) == 1]
  
  # LASSO knockoff filter
  LASSO_result = knockoff.filter(X, y, fdr=fdr,offset = 0,statistic = stat.lasso_coefdiff)
  LASSO_knockoff_selected = names(LASSO_result$selected)
  
  # SCAD knockoff filter
  SCAD_result = knockoff.filter(X, y,fdr=fdr,offset = 0,statistic = stat.SCAD)
  SCAD_knockoff_selected = names(SCAD_result$selected)
  
  # MCP knockoff filter
  MCP_result = knockoff.filter(X, y, fdr=fdr,offset = 0, statistic = stat.MCP)
  MCP_knockoff_selected = names(MCP_result$selected)
  
  list(LASSO_Knockoff = LASSO_knockoff_selected, SCAD_Knockoff = SCAD_knockoff_selected,MCP_Knockoff = MCP_knockoff_selected)
}

trials=100
fdr = 0.20
Res <- foreach(i=1:trials, .combine='rbind',.packages = c("knockoff","doParallel")) %dopar%
{
results = lapply(Y, function(y) myknockoff(X, y, fdr))
# print(results[1])

get_position <- function(x)
  sapply(regmatches(x, regexpr('[[:digit:]]+', x)), as.numeric)

comparisons <- lapply(results, function(drug) {
  lapply(drug, function(selected) {
    positions = unique(get_position(selected)) # remove possible duplicates
    discoveries = length(positions)
    false_discoveries = length(setdiff(positions, tsm_df$Position))
    list(true_discoveries = discoveries - false_discoveries,
         false_discoveries = false_discoveries,
         fdp = false_discoveries / max(1, discoveries))
  })
})
return(list(comparisons=comparisons))
}
stopImplicitCluster()
stopCluster(cl)
save(Res,file = 'Res.Rdata')

