

# -------------------------------------------------------------------------
#
# Fitting Propensity in SEM
# Study 1 (16 manifest variables)
#
# -------------------------------------------------------------------------


library(doSNOW)
library(lavaan)
library(clusterGeneration)
library(matrixcalc)
library(foreach)
library(doParallel)
library(data.table)



# settings
p <- 16
N <- 500
onlypos <- TRUE
nreps <- 1e6


# for parallelization
ncluster <- 47
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
m1cfa <- '
g =~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10 + x11 + x12 + x13 + x14 + x15 + x16
'

# 2-factor cfa
m2cfa <- '
f1 =~ x1
f2 =~ x2
f1 =~ x3
f2 =~ x4
f1 =~ x5
f2 =~ x6
f1 =~ x7
f2 =~ x8
f1 =~ x9
f2 =~ x10
f1 =~ x11
f2 =~ x12
f1 =~ x13
f2 =~ x14
f1 =~ x15
f2 =~ x16
'

# 3-factor cfa
m3cfa <- '
f1 =~ x1 + x2 + x3 + x4 + x5
f2 =~ x6 + x7 + x8 + x9 + x10
f3 =~ x11 + x12 + x13 + x14 + x15 + x16
'

# 4-factor cfa
m4cfa <- '
f1 =~ x1 + x2 + x3 + x4
f2 =~ x5 + x6 + x7 + x8
f3 =~ x9 + x10 + x11 + x12
f4 =~ x13 + x14 + x15 + x16
'

# 5-factor cfa
m5cfa <- '
f1 =~ x1 + x2 + x3 
f2 =~ x4 + x5 + x6 
f3 =~ x7 + x8 + x9 
f4 =~ x10 + x11 + x12
f5 =~ x13 + x14 + x15 + x16
'

## in addition
# 2-factor efa
# 3-factor efa
# 4-factor efa
# 5-factor efa



# save function -----------------------------------------------------------

save.fun <- function(a, b){
  
  intermediate.res <- rbind(a,b)
  colnames(intermediate.res) <- c('run', 'N', 'nvar',
                                  paste0('sigma.', 1:(p*(p+1)*0.5)),
                                  paste0('cfa1.', c('npar', 'df', 'Fmin', 'LL', 'rmsea', 'srmr', 'cfi', 'tli', 'aic', 'bic', 'check.conv')),
                                  paste0('cfa2.', c('npar', 'df', 'Fmin', 'LL', 'rmsea', 'srmr', 'cfi', 'tli', 'aic', 'bic', 'check.conv')),
                                  paste0('cfa3.', c('npar', 'df', 'Fmin', 'LL', 'rmsea', 'srmr', 'cfi', 'tli', 'aic', 'bic', 'check.conv')),
                                  paste0('cfa4.', c('npar', 'df', 'Fmin', 'LL', 'rmsea', 'srmr', 'cfi', 'tli', 'aic', 'bic', 'check.conv')),
                                  paste0('cfa5.', c('npar', 'df', 'Fmin', 'LL', 'rmsea', 'srmr', 'cfi', 'tli', 'aic', 'bic', 'check.conv')),
                                  paste0('efa2.', c('npar', 'df', 'Fmin', 'LL', 'rmsea', 'srmr', 'cfi', 'tli', 'aic', 'bic', 'check.conv', 'errorMessage')),
                                  paste0('efa3.', c('npar', 'df', 'Fmin', 'LL', 'rmsea', 'srmr', 'cfi', 'tli', 'aic', 'bic', 'check.conv', 'errorMessage')),
                                  paste0('efa4.', c('npar', 'df', 'Fmin', 'LL', 'rmsea', 'srmr', 'cfi', 'tli', 'aic', 'bic', 'check.conv', 'errorMessage')),
                                  paste0('efa5.', c('npar', 'df', 'Fmin', 'LL', 'rmsea', 'srmr', 'cfi', 'tli', 'aic', 'bic', 'check.conv', 'errorMessage')),
                                  'detS', 'logDetS')  
  

    fwrite(intermediate.res, file = 'fitmeasures_study1_16p.csv', row.names = F, sep = ';')

  intermediate.res
}





# simulation --------------------------------------------------------------


results <- foreach(crep = 1:nreps, .combine = save.fun, .options.snow = opts) %dopar% {
  

  ## create correlation matrix
  csigma <- genPositiveDefMat(dim = p, # dimension of covariance matrix
                              covMethod = "onion", 
                              eta = 1, # uniform
                              rangeVar = c(1,1))$Sigma # unit variance per variable
  

  
  ## function adapted from ockhamSEM (10.1037/met0000422; https://github.com/falkcarl/ockhamSEM)
  if (onlypos == TRUE) {
    csigma = (csigma+1)/2 # ad-hoc correction to ensure positive manifold
  }
  if(!is.symmetric.matrix(csigma)){
    csigma <-round(csigma,5) # ad-hoc fix
  }
  
  while(!is.positive.definite(csigma)){
    csigma <- genPositiveDefMat(dim = p, # dimension of covariance matrix
                                covMethod = "onion", 
                                eta = 1, # uniform
                                rangeVar = c(1,1))$Sigma
    
    if (onlypos) {
      csigma = (csigma+1)/2 # ad-hoc correction to ensure positive manifold
    }
    if(!is.symmetric.matrix(csigma)){
      csigma <-round(csigma,5) # ad-hoc fix
    }
    
  } # end while looop 
  
  colnames(csigma) <- rownames(csigma) <- paste0('x', 1:p)
  
  csigma.vec <- as.numeric(vech(csigma))
  
  
  

  tryCatch({
    
    ## 1-factor cfa
    res.1cfa <- sem(m1cfa, sample.cov = csigma, sample.nobs = N, sample.cov.rescale = F, std.lv = T, se = 'none', check.lv.names = F, check.start = F, check.gradient = F, check.post = F, check.vcov = F, warn = F)

    fit.1cfa <- getFitMeasures(res.1cfa)
    checks.1cfa <- check_convergence(res.1cfa)

    results.1cfa <- c(fit.1cfa, checks.1cfa)
    
    
    ## 2-factor cfa
    res.2cfa <- sem(m2cfa, sample.cov = csigma, sample.nobs = N, sample.cov.rescale = F, std.lv = T, se = 'none', check.lv.names = F, check.start = F, check.gradient = F, check.post = F, check.vcov = F, warn = F)

    fit.2cfa <- getFitMeasures(res.2cfa)
    checks.2cfa <- check_convergence(res.2cfa)
    
    results.2cfa <- c(fit.2cfa, checks.2cfa)
    
    
    ## 3-factor cfa
    res.3cfa <- sem(m3cfa, sample.cov = csigma, sample.nobs = N, sample.cov.rescale = F, std.lv = T, se = 'none', check.lv.names = F, check.start = F, check.gradient = F, check.post = F, check.vcov = F, warn = F)
    
    fit.3cfa <- getFitMeasures(res.3cfa)
    checks.3cfa <- check_convergence(res.3cfa)
    
    results.3cfa <- c(fit.3cfa, checks.3cfa)
    
    
    ## 4-factor cfa
    res.4cfa <- sem(m4cfa, sample.cov = csigma, sample.nobs = N, sample.cov.rescale = F, std.lv = T, se = 'none', check.lv.names = F, check.start = F, check.gradient = F, check.post = F, check.vcov = F, warn = F)
    
    fit.4cfa <- getFitMeasures(res.4cfa)
    checks.4cfa <- check_convergence(res.4cfa)
    
    results.4cfa <- c(fit.4cfa, checks.4cfa)
    
    
    ## 5-factor cfa
    res.5cfa <- sem(m5cfa, sample.cov = csigma, sample.nobs = N, sample.cov.rescale = F, std.lv = T, se = 'none', check.lv.names = F, check.start = F, check.gradient = F, check.post = F, check.vcov = F, warn = F)
    
    fit.5cfa <- getFitMeasures(res.5cfa)
    checks.5cfa <- check_convergence(res.5cfa)
    
    results.5cfa <- c(fit.5cfa, checks.5cfa)
    
    
  
    ## 2-factor EFA
    nfactors <- 2
    res.efa <- NULL
    check <- try(res.efa <- EFA(csigma, n_factors = nfactors, N = N, type = 'psych', method = 'ML', rotation = 'none'), silent = T)
    
    if(class(check) == 'try-error' & is.null(res.efa)){
      errorMessage <- as.character(check)
      results.2efa <- c(rep(NA, 8), errorMessage)
      
    } else{
      
      fit.efa <- getFitMeasures(res.efa)
      
      # check whether model converged
      check.conv.efa <- !res.efa$convergence

      if(class(check) == 'try-error'){errorMessage <- as.character(check)
      } else{errorMessage <- NA}
      
      results.2efa <- c(fit.efa, check.conv.efa, errorMessage)
    }
    
    
    ## 3-factor EFA
    nfactors <- 3
    res.efa <- NULL
    check <- try(res.efa <- EFA(csigma, n_factors = nfactors, N = N, type = 'psych', method = 'ML', rotation = 'none'), silent = T)
    
    if(class(check) == 'try-error' & is.null(res.efa)){
      errorMessage <- as.character(check)
      results.3efa <- c(rep(NA, 8), errorMessage)
      
    } else{
      
      fit.efa <- getFitMeasures(res.efa)
      
      # check whether model converged
      check.conv.efa <- !res.efa$convergence

      if(class(check) == 'try-error'){errorMessage <- as.character(check)
      } else{errorMessage <- NA}
      
      results.3efa <- c(fit.efa, check.conv.efa, errorMessage)
    }
    
    
    ## 4-factor EFA
    nfactors <- 4
    res.efa <- NULL
    check <- try(res.efa <- EFA(csigma, n_factors = nfactors, N = N, type = 'psych', method = 'ML', rotation = 'none'), silent = T)
    
    if(class(check) == 'try-error' & is.null(res.efa)){
      errorMessage <- as.character(check)
      results.4efa <- c(rep(NA, 8), errorMessage)
      
    } else{
      
      fit.efa <- getFitMeasures(res.efa)
      
      # check whether model converged
      check.conv.efa <- !res.efa$convergence

      if(class(check) == 'try-error'){errorMessage <- as.character(check)
      } else{errorMessage <- NA}
      
      results.4efa <- c(fit.efa, check.conv.efa, errorMessage)
    }
    
    
    ## 5-factor EFA
    nfactors <- 5
    res.efa <- NULL
    check <- try(res.efa <- EFA(csigma, n_factors = nfactors, N = N, type = 'psych', method = 'ML', rotation = 'none'), silent = T)
    
    if(class(check) == 'try-error' & is.null(res.efa)){
      errorMessage <- as.character(check)
      results.5efa <- c(rep(NA, 8), errorMessage)
      
    } else{
      
      fit.efa <- getFitMeasures(res.efa)
      
      # check whether model converged
      check.conv.efa <- !res.efa$convergence

      if(class(check) == 'try-error'){errorMessage <- as.character(check)
      } else{errorMessage <- NA}
      
      results.5efa <- c(fit.efa, check.conv.efa, errorMessage)
    }
    
    
    # det S
    detS <- det(csigma)
    
    # log det S
    logDetS <- log(det(csigma))
  
    
    ## save results
    
    res <- matrix(c(crep, N, p, csigma.vec, results.1cfa, results.2cfa, results.3cfa, results.4cfa, results.5cfa, results.2efa, results.3efa, results.4efa, results.5efa, detS, logDetS), nrow = 1)
    res
    


  }, warning = function(w) {
    #print(paste('WARNING: ',w))
  }, error = function(e) {
    # print(paste('ERROR: ',e))
  })
  
  
} # end foreach



stopCluster(cluster)

