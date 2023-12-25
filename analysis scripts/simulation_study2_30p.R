

# -------------------------------------------------------------------------
#
# Fitting Propensity in SEM
# Study 2 (30 manifest variables)
#
# -------------------------------------------------------------------------


library(doSNOW)
library(lavaan)
library(clusterGeneration)
library(matrixcalc)
library(foreach)
library(doParallel)
library(data.table)



## load sigmas
sigmas <- as.data.frame(fread("sigmas_30p.csv", sep=";"))


# -------------------------------------------------------------------------

# settings
p <- 30
N <- 500
onlypos <- TRUE
nreps <- nrow(sigmas)

# for parallelization
ncluster <- 20
cluster <- makeCluster(ncluster, outfile = "")
registerDoParallel(cluster)
registerDoSNOW(cluster)

# progress bar
pb <- txtProgressBar(max = nreps, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)


clusterCall(cluster, function(){
  library(lavaan)
  library(clusterGeneration)
  library(matrixcalc)
  library(psych)
  library(EFAtools)
  library(data.table)
})







# functions ---------------------------------------------------------------

getFitMeasures <- function(object){
  
  # for CFA models estimated with lavaan
  if(grepl('lavaan', class(object))){ 
    
    # free parameters
    npar <- object@Fit@npar
    
    # degrees of freedom
    df <- p*(p+1)*0.5 - npar
    
    # sample cov
    S <- object@SampleStats@cov[[1]]
    
    # model-implied cov
    SigmaHat <- object@Fit@Sigma.hat[[1]]
    
    
    
  } else { # for EFA models  
    
    ## free parameters
    npar <- p*nfactors + nfactors*(nfactors+1)*0.5 + p - nfactors^2
    
    ## df
    df <- 0.5 * ((p - nfactors)^2 - p - nfactors)
    
    ## sample cov
    S <- csigma
    
    ## model-implied cov
    lambdaHat <- unclass(object$unrot_loadings)
    thetaHat <- diag(1-object$h2)
    SigmaHat <- lambdaHat %*% t(lambdaHat) + thetaHat
    
  }
  
  
  # minimum of fit function
  Fmin <- sum(diag(S %*% solve(SigmaHat))) + log(det(SigmaHat)) - log(det(S)) - ncol(S)
  
  
  # log-Likelihood
  LL <- -(N/2) * (p * log(2*pi) + log(det(SigmaHat)) + sum(diag(S %*% solve(SigmaHat))))
  
  
  # RMSEA
  rmsea <- sqrt(max(0, (Fmin/df - 1/(N-1))))
  
  # SRMR
  sqrt.d <- 1/sqrt(diag(S))
  D <- diag(sqrt.d, ncol=length(sqrt.d))
  R <- D %*% (S - SigmaHat) %*% D
  e <- p*(p+1)/2
  srmr <-  sqrt( sum(R[lower.tri(R, diag=TRUE)]^2) / e )
  
  # null model
  SigmaHat.null <- diag(1, p, p)
  Fmin.null <- sum(diag(S %*% solve(SigmaHat.null))) + log(det(SigmaHat.null)) - log(det(S)) - ncol(S)
  df.null <- p*(p+1)*0.5 - p
  
  # CFI
  cfi <- 1 - (max(0, (N-1)*Fmin-df)/max(0, (N-1)*Fmin.null-df.null, (N-1)*Fmin-df))
  
  # TLI / NNFI
  tli <- 1 - ((((N-1)*Fmin - df)*df.null)/((((N-1)*Fmin.null - df.null)*df)))
  
  # AIC
  aic <- (-2 * LL) + (2 * npar)
  
  # BIC
  bic <- (-2 * LL) + (log(N) * npar)
  
  
  fit <- c(npar, df, Fmin, LL, rmsea, srmr, cfi, tli, aic, bic)
  names(fit) <- c('npar', 'df', 'Fmin', 'LL', 'rmsea', 'srmr', 'cfi', 'tli', 'aic', 'bic')  
  
  return(fit)
  
}


# check whether model converged

check_convergence <- function(object){
  object@optim[["converged"]]
  
}


# measurement models ------------------------------------------------------

# 1-factor cfa
mcfa1 <- '
g =~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10 + x11 + x12 + x13 + x14 + x15 + x16 + x17 + x18 + x19 + x20 + x21 + x22 + x23 + x24 + x25 + x26 + x27 + x28 + x29 + x30
'

# 5-factor cfa
mcfa5 <- '
f1 =~ x1 + x2 + x3 + x4 + x5 + x6 
f2 =~ x7 + x8 + x9 + x10 + x11 + x12
f3 =~ x13 + x14 + x15 + x16 + x17 + x18
f4 =~ x19 + x20 + x21 + x22 + x23 + x24
f5 =~ x25 + x26 + x27 + x28 + x29 + x30
'

# 6-factor cfa
mcfa6 <- '
f1 =~ x1 + x2 + x3 + x4 + x5 
f2 =~ x6 + x7 + x8 + x9 + x10
f3 =~ x11 + x12 + x13 + x14 + x15
f4 =~ x16 + x17 + x18 + x19 + x20
f5 =~ x21 + x22 + x23 + x24 + x25
f6 =~ x26 + x27 + x28 + x29 + x30
'

# Higher-Order (5 first-order factors)
mho5 <- '
f1 =~ x1 + x2 + x3 + x4 + x5 + x6 
f2 =~ x7 + x8 + x9 + x10 + x11 + x12
f3 =~ x13 + x14 + x15 + x16 + x17 + x18
f4 =~ x19 + x20 + x21 + x22 + x23 + x24
f5 =~ x25 + x26 + x27 + x28 + x29 + x30
g =~ f1 + f2 + f3 + f4 + f5
'

# Higher-order (6 first-order factors)
mho6 <- '
f1 =~ x1 + x2 + x3 + x4 + x5 
f2 =~ x6 + x7 + x8 + x9 + x10
f3 =~ x11 + x12 + x13 + x14 + x15
f4 =~ x16 + x17 + x18 + x19 + x20
f5 =~ x21 + x22 + x23 + x24 + x25
f6 =~ x26 + x27 + x28 + x29 + x30
g =~ f1 + f2 + f3 + f4 + f5 + f6
'

# Bifactor (5 specific factors)
mbif5 <- '
f1 =~ x1 + x2 + x3 + x4 + x5 + x6 
f2 =~ x7 + x8 + x9 + x10 + x11 + x12
f3 =~ x13 + x14 + x15 + x16 + x17 + x18
f4 =~ x19 + x20 + x21 + x22 + x23 + x24
f5 =~ x25 + x26 + x27 + x28 + x29 + x30
g =~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10 + x11 + x12 + x13 + x14 + x15 + x16 + x17 + x18 + x19 + x20 + x21 + x22 + x23 + x24 + x25 + x26 + x27 + x28 + x29 + x30
'

# Bifactor (6 specific factors)
mbif6 <- '
f1 =~ x1 + x2 + x3 + x4 + x5 
f2 =~ x6 + x7 + x8 + x9 + x10
f3 =~ x11 + x12 + x13 + x14 + x15
f4 =~ x16 + x17 + x18 + x19 + x20
f5 =~ x21 + x22 + x23 + x24 + x25
f6 =~ x26 + x27 + x28 + x29 + x30
g =~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10 + x11 + x12 + x13 + x14 + x15 + x16 + x17 + x18 + x19 + x20 + x21 + x22 + x23 + x24 + x25 + x26 + x27 + x28 + x29 + x30
'





# save function -----------------------------------------------------------

save.fun <- function(a, b){
  
  intermediate.res <- rbind(a,b)
  colnames(intermediate.res) <- c('run', 'N', 'nvar',
                                  paste0('cfa1.', c('npar', 'df', 'Fmin', 'LL', 'rmsea', 'srmr', 'cfi', 'tli', 'aic', 'bic', 'check.conv')),
                                  paste0('cfa5.', c('npar', 'df', 'Fmin', 'LL', 'rmsea', 'srmr', 'cfi', 'tli', 'aic', 'bic', 'check.conv')),
                                  paste0('cfa6.', c('npar', 'df', 'Fmin', 'LL', 'rmsea', 'srmr', 'cfi', 'tli', 'aic', 'bic', 'check.conv')),
                                  paste0('ho5.', c('npar', 'df', 'Fmin', 'LL', 'rmsea', 'srmr', 'cfi', 'tli', 'aic', 'bic', 'check.conv')),
                                  paste0('ho6.', c('npar', 'df', 'Fmin', 'LL', 'rmsea', 'srmr', 'cfi', 'tli', 'aic', 'bic', 'check.conv')),
                                  paste0('bif5.', c('npar', 'df', 'Fmin', 'LL', 'rmsea', 'srmr', 'cfi', 'tli', 'aic', 'bic', 'check.conv')),
                                  paste0('bif6.', c('npar', 'df', 'Fmin', 'LL', 'rmsea', 'srmr', 'cfi', 'tli', 'aic', 'bic', 'check.conv')),
                                  'detS', 'logDetS')
  
  fwrite(intermediate.res, file = 'fitmeasures_study2_30p.csv', row.names = F, sep = ';')

  intermediate.res
}





# simulation --------------------------------------------------------------



results <- foreach(crep = 1:nreps, .combine = save.fun, .options.snow = opts) %dopar% {
  

  ## re-create correlation matrix
  v <- as.numeric(sigmas[sigmas$run == crep, grepl('sigma', names(sigmas))])
  csigma <- matrix(0, nrow = p, ncol = p)
  csigma[lower.tri(csigma, diag = T)] <- v
  csigma <- t(csigma)
  csigma[lower.tri(csigma, diag = T)] <- v
  
  colnames(csigma) <- rownames(csigma) <- paste0('x', 1:p)
  
  

  tryCatch({
    
    ## 1-factor cfa
    res.cfa1 <- sem(mcfa1, sample.cov = csigma, sample.nobs = N, sample.cov.rescale = F, std.lv = T, se = 'none', check.lv.names = F, check.start = F, check.gradient = F, check.post = F, check.vcov = F, warn = F)
    
    fit.cfa1 <- getFitMeasures(res.cfa1)
    checks.cfa1 <- check_convergence(res.cfa1)
    
    results.cfa1 <- c(fit.cfa1, checks.cfa1)
    
    
    
    ## 5-factor cfa 
    res.cfa5 <- sem(mcfa5, sample.cov = csigma, sample.nobs = N, sample.cov.rescale = F, std.lv = T, se = 'none', check.lv.names = F, check.start = F, check.gradient = F, check.post = F, check.vcov = F, warn = F)

    fit.cfa5 <- getFitMeasures(res.cfa5)
    checks.cfa5 <- check_convergence(res.cfa5)

    results.cfa5 <- c(fit.cfa5, checks.cfa5)
    
    
    ## 6-factor cfa 
    res.cfa6 <- sem(mcfa6, sample.cov = csigma, sample.nobs = N, sample.cov.rescale = F, std.lv = T, se = 'none', check.lv.names = F, check.start = F, check.gradient = F, check.post = F, check.vcov = F, warn = F)
    
    fit.cfa6 <- getFitMeasures(res.cfa6)
    checks.cfa6 <- check_convergence(res.cfa6)
    
    results.cfa6 <- c(fit.cfa6, checks.cfa6)
    
    
    ## Higher-Order (5 first-order factors) 
    res.ho5 <- sem(mho5, sample.cov = csigma, sample.nobs = N, sample.cov.rescale = F, std.lv = T, se = 'none', check.lv.names = F, check.start = F, check.gradient = F, check.post = F, check.vcov = F, warn = F)
    
    fit.ho5 <- getFitMeasures(res.ho5)
    checks.ho5 <- check_convergence(res.ho5)
    
    results.ho5 <- c(fit.ho5, checks.ho5)
    
    
    ## Higher-Order (6 first-order factors)
    res.ho6 <- sem(mho6, sample.cov = csigma, sample.nobs = N, sample.cov.rescale = F, std.lv = T, se = 'none', check.lv.names = F, check.start = F, check.gradient = F, check.post = F, check.vcov = F, warn = F)
    
    fit.ho6 <- getFitMeasures(res.ho6)
    checks.ho6 <- check_convergence(res.ho6)
    
    results.ho6 <- c(fit.ho6, checks.ho6)
    
    
    ## Bifactor (5 specific factors) 
    res.bif5 <- sem(mbif5, sample.cov = csigma, sample.nobs = N, sample.cov.rescale = F, std.lv = T, se = 'none', check.lv.names = F, check.start = F, check.gradient = F, check.post = F, check.vcov = F, warn = F)
    
    fit.bif5 <- getFitMeasures(res.bif5)
    checks.bif5 <- check_convergence(res.bif5)
    
    results.bif5 <- c(fit.bif5, checks.bif5)
    
    
    ## Bifactor (6 specific factors)
    res.bif6 <- sem(mbif6, sample.cov = csigma, sample.nobs = N, sample.cov.rescale = F, std.lv = T, se = 'none', check.lv.names = F, check.start = F, check.gradient = F, check.post = F, check.vcov = F, warn = F)
    
    fit.bif6 <- getFitMeasures(res.bif6)
    checks.bif6 <- check_convergence(res.bif6)
    
    results.bif6 <- c(fit.bif6, checks.bif6)
    
    # det S
    detS <- det(csigma)
    
    # log det S
    logDetS <- log(det(csigma))
    
  
    ## save results
    
    res <- matrix(c(crep, N, p, results.cfa1, results.cfa5, results.cfa6, results.ho5, results.ho6, results.bif5, results.bif6, detS, logDetS), nrow = 1)
    res


  }, warning = function(w) {
    #print(paste('WARNING: ',w))
  }, error = function(e) {
    # print(paste('ERROR: ',e))
  })
  
  
} # end foreach



stopCluster(cluster)




