

# -------------------------------------------------------------------------
#
# Fitting Propensity in SEM
# Study 1 (30 manifest variables)
#
# -------------------------------------------------------------------------


library(doSNOW)
library(lavaan)
library(clusterGeneration)
library(matrixcalc)
library(foreach)
library(doParallel)
library(data.table)
library(EFAtools)



# settings
p <- 30
N <- 500
onlypos <- TRUE
nreps <- 1e6

# for parallelization
ncluster <- 46
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

# 1 CFA
mcfa1 <- '
g =~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10 + x11 + x12 + x13 + x14 + x15 + x16 + x17 + x18 + x19 + x20 + x21 + x22 + x23 + x24 + x25 + x26 + x27 + x28 + x29 + x30
'

# 2 CFA
mcfa2 <- '
f1 =~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10 + x11 + x12 + x13 + x14 + x15
f2 =~ x16 + x17 + x18 + x19 + x20 + x21 + x22 + x23 + x24 + x25 + x26 + x27 + x28 + x29 + x30
'

# 3 CFA
mcfa3 <- '
f1 =~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10
f2 =~ x11 + x12 + x13 + x14 + x15 + x16 + x17 + x18 + x19 + x20
f3 =~ x21 + x22 + x23 + x24 + x25 + x26 + x27 + x28 + x29 + x30
'

# 4 CFA
mcfa4 <- '
f1 =~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8
f2 =~ x9 + x10 + x11 + x12 + x13 + x14 + x15 + x16
f3 =~ x17 + x18 + x19 + x20 + x21 + x22 + x23 
f4 =~ x24 + x25 + x26 + x27 + x28 + x29 + x30
'

# 5 CFA
mcfa5 <- '
f1 =~ x1 + x2 + x3 + x4 + x5 + x6 
f2 =~ x7 + x8 + x9 + x10 + x11 + x12
f3 =~ x13 + x14 + x15 + x16 + x17 + x18
f4 =~ x19 + x20 + x21 + x22 + x23 + x24
f5 =~ x25 + x26 + x27 + x28 + x29 + x30
'

# 6 CFA
mcfa6 <- '
f1 =~ x1 + x2 + x3 + x4
f6 =~ x5
f2 =~ x6 + x7 + x8 + x9 + x10
f3 =~ x11 + x12 + x13 + x14 + x15
f4 =~ x16 + x17 + x18 + x19 + x20
f5 =~ x21 + x22 + x23 + x24 + x25
f6 =~ x26 + x27 + x28 + x29
f1 =~ x30
'

## in addition
# 2-factor efa
# 3-factor efa
# 4-factor efa
# 5-factor efa
# 6-factor efa





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
                                  paste0('cfa6.', c('npar', 'df', 'Fmin', 'LL', 'rmsea', 'srmr', 'cfi', 'tli', 'aic', 'bic', 'check.conv')),
                                  paste0('efa2.', c('npar', 'df', 'Fmin', 'LL', 'rmsea', 'srmr', 'cfi', 'tli', 'aic', 'bic', 'check.conv', 'errorMessage')),
                                  paste0('efa3.', c('npar', 'df', 'Fmin', 'LL', 'rmsea', 'srmr', 'cfi', 'tli', 'aic', 'bic', 'check.conv', 'errorMessage')),
                                  paste0('efa4.', c('npar', 'df', 'Fmin', 'LL', 'rmsea', 'srmr', 'cfi', 'tli', 'aic', 'bic', 'check.conv', 'errorMessage')),
                                  paste0('efa5.', c('npar', 'df', 'Fmin', 'LL', 'rmsea', 'srmr', 'cfi', 'tli', 'aic', 'bic', 'check.conv', 'errorMessage')),
                                  paste0('efa6.', c('npar', 'df', 'Fmin', 'LL', 'rmsea', 'srmr', 'cfi', 'tli', 'aic', 'bic', 'check.conv', 'errorMessage')),
                                  'detS', 'logDetS')  
  
    fwrite(intermediate.res, file = 'fitmeasures_study1_30p.csv', row.names = F, sep = ';')
    
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
    res.cfa1 <- sem(mcfa1, sample.cov = csigma, sample.nobs = N, sample.cov.rescale = F, std.lv = T, se = 'none', check.lv.names = F, check.start = F, check.gradient = F, check.post = F, check.vcov = F, warn = F)

    fit.cfa1 <- getFitMeasures(res.cfa1)
    checks.cfa1 <- check_convergence(res.cfa1)
    results.cfa1 <- c(fit.cfa1, checks.cfa1)
    
    
    ## 2-factor CFA
    res.cfa2 <- sem(mcfa2, sample.cov = csigma, sample.nobs = N, sample.cov.rescale = F, std.lv = T, se = 'none', check.lv.names = F, check.start = F, check.gradient = F, check.post = F, check.vcov = F, warn = F)
    
    fit.cfa2 <- getFitMeasures(res.cfa2)
    checks.cfa2 <- check_convergence(res.cfa2)
    results.cfa2 <- c(fit.cfa2, checks.cfa2)
    
    
    ## 3-factor CFA
    res.cfa3 <- sem(mcfa3, sample.cov = csigma, sample.nobs = N, sample.cov.rescale = F, std.lv = T, se = 'none', check.lv.names = F, check.start = F, check.gradient = F, check.post = F, check.vcov = F, warn = F)
    
    fit.cfa3 <- getFitMeasures(res.cfa3)
    checks.cfa3 <- check_convergence(res.cfa3)
    results.cfa3 <- c(fit.cfa3, checks.cfa3)
    
    
    ## 4-factor CFA
    res.cfa4 <- sem(mcfa4, sample.cov = csigma, sample.nobs = N, sample.cov.rescale = F, std.lv = T, se = 'none', check.lv.names = F, check.start = F, check.gradient = F, check.post = F, check.vcov = F, warn = F)
    
    fit.cfa4 <- getFitMeasures(res.cfa4)
    checks.cfa4 <- check_convergence(res.cfa4)
    results.cfa4 <- c(fit.cfa4, checks.cfa4)
 
    
    ## 5-factor CFA
    res.cfa5 <- sem(mcfa5, sample.cov = csigma, sample.nobs = N, sample.cov.rescale = F, std.lv = T, se = 'none', check.lv.names = F, check.start = F, check.gradient = F, check.post = F, check.vcov = F, warn = F)
    
    fit.cfa5 <- getFitMeasures(res.cfa5)
    checks.cfa5 <- check_convergence(res.cfa5)
    results.cfa5 <- c(fit.cfa5, checks.cfa5)
    
    
    ## 6-factor CFA
    res.cfa6 <- sem(mcfa6, sample.cov = csigma, sample.nobs = N, sample.cov.rescale = F, std.lv = T, se = 'none', check.lv.names = F, check.start = F, check.gradient = F, check.post = F, check.vcov = F, warn = F)
    
    fit.cfa6 <- getFitMeasures(res.cfa6)
    checks.cfa6 <- check_convergence(res.cfa6)
    results.cfa6 <- c(fit.cfa6, checks.cfa6)
    
    
    
    
    ## 2-factor EFA
    nfactors <- 2
    res.efa2 <- NULL
    check <- try(res.efa2 <- EFA(csigma, n_factors = nfactors, N = N, type = 'psych', method = 'ML', rotation = 'none'), silent = T)
    
    if(class(check) == 'try-error' & is.null(res.efa2)){
      errorMessage <- as.character(check)
      results.efa2 <- c(rep(NA, 8), errorMessage)
      
    } else{
      
      fit.efa2 <- getFitMeasures(res.efa2)
      # check whether model converged
      check.conv.efa2 <- !res.efa2$convergence 

      if(class(check) == 'try-error'){errorMessage <- as.character(check)
      } else{errorMessage <- NA}
      
      results.efa2 <- c(fit.efa2, check.conv.efa2, errorMessage)
    }    
    

    ## 3-factor EFA
    nfactors <- 3
    res.efa3 <- NULL
    check <- try(res.efa3 <- EFA(csigma, n_factors = nfactors, N = N, type = 'psych', method = 'ML', rotation = 'none'), silent = T)
    
    if(class(check) == 'try-error' & is.null(res.efa3)){
      errorMessage <- as.character(check)
      results.efa3 <- c(rep(NA, 8), errorMessage)
      
    } else{
      
       fit.efa3 <- getFitMeasures(res.efa3)
       # check whether model converged
       check.conv.efa3 <- !res.efa3$convergence 

      if(class(check) == 'try-error'){errorMessage <- as.character(check)
      } else{errorMessage <- NA}
      
      results.efa3 <- c(fit.efa3, check.conv.efa3, errorMessage)
    }
    
    
    ## 4-factor EFA
    nfactors <- 4
    res.efa4 <- NULL
    check <- try(res.efa4 <- EFA(csigma, n_factors = nfactors, N = N, type = 'psych', method = 'ML', rotation = 'none'), silent = T)
    
    if(class(check) == 'try-error' & is.null(res.efa4)){
      errorMessage <- as.character(check)
      results.efa4 <- c(rep(NA, 8), errorMessage)
      
    } else{
      
      fit.efa4 <- getFitMeasures(res.efa4)
      # check whether model converged
      check.conv.efa4 <- !res.efa4$convergence 
 
      if(class(check) == 'try-error'){errorMessage <- as.character(check)
      } else{errorMessage <- NA}
      
      results.efa4 <- c(fit.efa4, check.conv.efa4, errorMessage)
    }
    
    
    ## 5-factor EFA
    nfactors <- 5
    res.efa5 <- NULL
    check <- try(res.efa5 <- EFA(csigma, n_factors = nfactors, N = N, type = 'psych', method = 'ML', rotation = 'none'), silent = T)
    
    if(class(check) == 'try-error' & is.null(res.efa5)){
      errorMessage <- as.character(check)
      results.efa5 <- c(rep(NA, 8), errorMessage)
      
    } else{
      
      fit.efa5 <- getFitMeasures(res.efa5)
      # check whether model converged
      check.conv.efa5 <- !res.efa5$convergence 

      if(class(check) == 'try-error'){errorMessage <- as.character(check)
      } else{errorMessage <- NA}
      
      results.efa5 <- c(fit.efa5, check.conv.efa5, errorMessage)
    }
    
    
    ## 6-factor EFA
    nfactors <- 6
    res.efa6 <- NULL
    check <- try(res.efa6 <- EFA(csigma, n_factors = nfactors, N = N, type = 'psych', method = 'ML', rotation = 'none'), silent = T)
    
    if(class(check) == 'try-error' & is.null(res.efa6)){
      errorMessage <- as.character(check)
      results.efa6 <- c(rep(NA, 8), errorMessage)
      
    } else{
      
      fit.efa6 <- getFitMeasures(res.efa6)
      # check whether model converged
      check.conv.efa6 <- !res.efa6$convergence 

      if(class(check) == 'try-error'){errorMessage <- as.character(check)
      } else{errorMessage <- NA}
      
      results.efa6 <- c(fit.efa6, check.conv.efa6, errorMessage)
    }
    
    # det S
    detS <- det(csigma)
    
    # log det S
    logDetS <- log(det(csigma))
    
  
    ## save results
    
    res <- matrix(c(crep, N, p, csigma.vec, results.cfa1, results.cfa2, results.cfa3, results.cfa4, results.cfa5, results.cfa6, results.efa2, results.efa3, results.efa4, results.efa5, results.efa6, detS, logDetS), nrow = 1)
    res
    


  }, warning = function(w) {
    #print(paste('WARNING: ',w))
  }, error = function(e) {
    # print(paste('ERROR: ',e))
  })
  
  
} # end foreach



stopCluster(cluster)

