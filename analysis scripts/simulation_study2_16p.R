

# -------------------------------------------------------------------------
#
# Fitting Propensity in SEM
# Study 2 (16 manifest variables)
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
sigmas <- as.data.frame(fread("sigmas_16p.csv", sep=";"))


# -------------------------------------------------------------------------

# settings
p <- 16
N <- 500
onlypos <- TRUE
nreps <- nrow(sigmas)

# for parallelization
ncluster <- 45
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
g =~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10 + x11 + x12 + x13 + x14 + x15 + x16
'

# 4-factor cfa
mcfa4 <- '
f1 =~ x1 + x2 + x3 + x4
f2 =~ x5 + x6 + x7 + x8
f3 =~ x9 + x10 + x11 + x12
f4 =~ x13 + x14 + x15 + x16
'

# higher-order
mho <- '
f1 =~ x1 + x2 + x3 + x4
f2 =~ x5 + x6 + x7 + x8
f3 =~ x9 + x10 + x11 + x12
f4 =~ x13 + x14 + x15 + x16
g =~ f1 + f2 + f3 + f4
'

# bifactor 
mbif <- '
f1 =~ x1 + x2 + x3 + x4
f2 =~ x5 + x6 + x7 + x8
f3 =~ x9 + x10 + x11 + x12
f4 =~ x13 + x14 + x15 + x16
bf =~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10 + x11 + x12 + x13 + x14 + x15 + x16

f1+ f2+ f3 + f4 ~~ 0*bf
f1 + f2 + f3 ~~ 0*f4
f1 + f2 ~~ 0*f3
f1 ~~ 0*f2
'





# save function -----------------------------------------------------------

save.fun <- function(a, b){
  
  intermediate.res <- rbind(a,b)
  colnames(intermediate.res) <- c('run', 'N', 'nvar',
                                  paste0('cfa1.', c('npar', 'df', 'Fmin', 'LL', 'rmsea', 'srmr', 'cfi', 'tli', 'aic', 'bic', 'check.conv')),
                                  paste0('cfa4.', c('npar', 'df', 'Fmin', 'LL', 'rmsea', 'srmr', 'cfi', 'tli', 'aic', 'bic', 'check.conv')),
                                  paste0('ho.', c('npar', 'df', 'Fmin', 'LL', 'rmsea', 'srmr', 'cfi', 'tli', 'aic', 'bic', 'check.conv')),
                                  paste0('bif.', c('npar', 'df', 'Fmin', 'LL', 'rmsea', 'srmr', 'cfi', 'tli', 'aic', 'bic', 'check.conv')),
                                  'detS', 'logDetS')
  
  fwrite(intermediate.res, file = 'fitmeasures_study2_16p.csv', row.names = F, sep = ';')

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
    
    
    
    ## 4-factor cfa 
    res.cfa4 <- sem(mcfa4, sample.cov = csigma, sample.nobs = N, sample.cov.rescale = F, std.lv = T, se = 'none', check.lv.names = F, check.start = F, check.gradient = F, check.post = F, check.vcov = F, warn = F)

    fit.cfa4 <- getFitMeasures(res.cfa4)
    checks.cfa4 <- check_convergence(res.cfa4)

    results.cfa4 <- c(fit.cfa4, checks.cfa4)
    
    
    ## higher-order
    res.ho <- sem(mho, sample.cov = csigma, sample.nobs = N, sample.cov.rescale = F, std.lv = T, se = 'none', check.lv.names = F, check.start = F, check.gradient = F, check.post = F, check.vcov = F, warn = F)
    
    fit.ho <- getFitMeasures(res.ho)
    checks.ho <- check_convergence(res.ho)
    
    results.ho <- c(fit.ho, checks.ho)
    
    
    ## bifactor
    res.bif <- sem(mbif, sample.cov = csigma, sample.nobs = N, sample.cov.rescale = F, std.lv = T, se = 'none', check.lv.names = F, check.start = F, check.gradient = F, check.post = F, check.vcov = F, warn = F)
    
    fit.bif <- getFitMeasures(res.bif)
    checks.bif <- check_convergence(res.bif)
    
    results.bif <- c(fit.bif, checks.bif)
    
    # det S
    detS <- det(csigma)
    
    # log det S
    logDetS <- log(det(csigma))
 
    
  
    ## save results
    
    res <- matrix(c(crep, N, p, results.cfa1, results.cfa4, results.ho, results.bif, detS, logDetS), nrow = 1)
    res


  }, warning = function(w) {
    #print(paste('WARNING: ',w))
  }, error = function(e) {
    # print(paste('ERROR: ',e))
  })
  
  
} # end foreach



stopCluster(cluster)




