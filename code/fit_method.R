## get a betahat and a sebetahat for each gene, as well as the associated df, from a count matrix with genes in rows and samples in columns, and a condition vector


library(edgeR)
library(limma)

# Voom transformation + limma fit



# TMM normalization + Voom transformation + limma fit
voom_transform = function(counts, condition){
  dgecounts = edgeR::calcNormFactors(edgeR::DGEList(counts = counts, group = condition))

  design = model.matrix(~ condition)

  v = limma::voom(dgecounts, design, plot=FALSE)

  lim = limma::lmFit(v)
  
  betahat.voom = lim$coefficients[, 2]
  sebetahat.voom = lim$stdev.unscaled[, 2] * lim$sigma
  df.voom = length(condition) - 2
  
  return(list(betahat = betahat.voom, sebetahat = sebetahat.voom, df = df.voom, v = v))
}


##' Simple ordinary least squares.
##'
##' @param log_counts A matrix of numerics. The responses.
##' @param condition A matrix of numerics. The predictors.
get_ols <- function(log_counts, condition) {
  limma_out <- limma::lmFit(log_counts, model.matrix(~condition))
  betahat <- limma_out$coefficients[, 2]
  sebetahat <- limma_out$stdev.unscaled[, 2] * limma_out$sigma
  df <- limma_out$df.residual
  return(list(betahat = betahat, sebetahat = sebetahat, df = df))
}