

# -------------------------------------------------------------------------
#
# Fitting Propensity in SEM
# Analysis Study 2
#
# -------------------------------------------------------------------------

library(data.table)
library(effsize)
library(matrixStats)

options(scipen = 9999)


data.9p.raw <- as.data.frame(fread('fitmeasures_study2_9p.csv'))
data.16p.raw <- as.data.frame(fread('fitmeasures_study2_16p.csv'))
data.30p.raw <- as.data.frame(fread('fitmeasures_study2_30p.csv'))


# only converged solutions ------------------------------------------------

# 9 manifest variables
all.converged <- apply(data.9p.raw[, grep('check.conv', names(data.9p.raw), value = TRUE)], 1, sum)
data.9p <- data.9p.raw[all.converged == length(grep('check.conv', names(data.9p.raw))), ]
nrow(data.9p)/nrow(data.9p.raw) # convergence rate

# 16 manifest variables
all.converged <- apply(data.16p.raw[, grep('check.conv', names(data.16p.raw), value = TRUE)], 1, sum)
data.16p <- data.16p.raw[all.converged == length(grep('check.conv', names(data.16p.raw))), ]
nrow(data.16p)/nrow(data.16p.raw) # convergence rate

# 30 manifest variables
all.converged <- apply(data.30p.raw[, grep('check.conv', names(data.30p.raw), value = TRUE)], 1, sum)
data.30p <- data.30p.raw[all.converged == length(grep('check.conv', names(data.30p.raw))), ]
nrow(data.30p)/nrow(data.30p.raw) # convergence rate



# 9 manifest variables ----------------------------------------------------

# models
models.9p <- sub(pattern = '.npar', grep('npar', names(data.9p), value = TRUE), replacement = '')

# free parameters
npar.9p <- data.9p[1, grep('npar', names(data.9p))]; npar.9p

# degrees of freedom
df.9p <- data.9p[1, grep('df', names(data.9p))]; df.9p



# descriptive statistics
measures <- c('Fmin', 'cfi', 'srmr', 'rmsea', 'tli', 'aic', 'bic')
colnames <- c('k', 'df', 'Fmin mean', 'Fmin SD', 'CFI mean', 'CFI SD', 'SRMR mean', 'SRMR SD', 'RMSEA mean', 'RMSEA SD', 'TLI mean', 'TLI SD', 'AIC mean', 'AIC SD', 'BIC mean', 'BIC SD')
res.9p <- as.data.frame(matrix(ncol = length(colnames), nrow = length(npar.9p), dimnames = list(models.9p, colnames)))

res.9p[, 1] <- as.numeric(npar.9p)
res.9p[, 2] <- as.numeric(df.9p)

for(i in 1:length(measures)){
  cdat <- data.9p[, grep(measures[i], names(data.9p))]
  res.9p[, 2*i+1] <- apply(cdat, 2, mean)
  res.9p[, 2*i+2] <- apply(cdat, 2, sd)
}
round(res.9p, 2)

# C*NML
Fmin.9p <- data.9p[, grep('Fmin', names(data.9p))]
logdetS <- data.9p$logDetS
apply(Fmin.9p, 2, function(x) logSumExp(-.5*(x + logdetS)))
round(apply(Fmin.9p, 2, function(x) logSumExp(-.5*(x + logdetS))), 2)

# cfi
cfi.9p <- data.9p[, grep('cfi', names(data.9p))]
apply(cfi.9p, 2, FUN = function(x) sum(x >= .90)/length(x))
apply(cfi.9p, 2, FUN = function(x) sum(x >= .90))


# srmr
srmr.9p <- data.9p[, grep('srmr', names(data.9p))]
apply(srmr.9p, 2, FUN = function(x) sum(x <= .10)/length(x))
apply(srmr.9p, 2, FUN = function(x) sum(x <= .08)/length(x))
apply(srmr.9p, 2, FUN = function(x) sum(x <= .06)/length(x))


# rmsea
rmsea.9p <- data.9p[, grep('rmsea', names(data.9p))]
apply(rmsea.9p, 2, FUN = function(x) sum(x <= .08)/length(x))


# similarity of rmsea distributions
rmsea.models.9p <- names(rmsea.9p)
pairs <- combn(rmsea.models.9p, 2)
compare.rmsea.9p <- matrix(nrow=ncol(pairs), ncol = 3)
colnames(compare.rmsea.9p) <- c('pair', 'cliffs delta', 'kolmogorov-smirnov')

for(i in 1:ncol(pairs)){
  mod.a <- pairs[1,i]
  mod.b <- pairs[2,i]
  # Cliff's delta
  cliff <- cliff.delta(data.9p[,mod.a], data.9p[,mod.b])$estimate
  # Kolmogorov Smirnov
  kolm <- ks.test(data.9p[,mod.a], data.9p[,mod.b])$statistic
  compare.rmsea.9p[i,1] <- gsub(x=paste(mod.a, mod.b, sep = ' vs '), pattern='.rmsea', replacement='')
  compare.rmsea.9p[i, 2:3] <- c(cliff, kolm)
}

compare.rmsea.9p <- as.data.frame(compare.rmsea.9p)
compare.rmsea.9p

round(median(abs(as.numeric(compare.rmsea.9p$`cliffs delta`))), 2)
round(median(abs(as.numeric(compare.rmsea.9p$`kolmogorov-smirnov`))), 2)


# tli
tli.9p <- data.9p[, grep('tli', names(data.9p))]
apply(tli.9p, 2, FUN = function(x) sum(x >= .95))
apply(tli.9p, 2, FUN = function(x) sum(x >= .95)/length(x))

# similarity of TLI distributions
tli.models.9p <- names(tli.9p)
pairs <- combn(tli.models.9p, 2)
compare.tli.9p <- matrix(nrow=ncol(pairs), ncol = 3)
colnames(compare.tli.9p) <- c('pair', 'cliffs delta', 'kolmogorov-smirnov')

for(i in 1:ncol(pairs)){
  mod.a <- pairs[1,i]
  mod.b <- pairs[2,i]
  # Cliff's delta
  cliff <- cliff.delta(data.9p[,mod.a], data.9p[,mod.b])$estimate
  # Kolmogorov Smirnov
  kolm <- ks.test(data.9p[,mod.a], data.9p[,mod.b])$statistic
  compare.tli.9p[i,1] <- gsub(x=paste(mod.a, mod.b, sep = ' vs '), pattern='.tli', replacement='')
  compare.tli.9p[i, 2:3] <- c(cliff, kolm)
}

compare.tli.9p <- as.data.frame(compare.tli.9p)
compare.tli.9p

round(median(abs(as.numeric(compare.tli.9p$`cliffs delta`))), 2)
round(median(abs(as.numeric(compare.tli.9p$`kolmogorov-smirnov`))), 2)



# supplement: CNML
LL.9p <- data.9p[, grep('LL', names(data.9p))]
round(apply(LL.9p, 2, logSumExp))



# 16 manifest variables ---------------------------------------------------

# models
models.16p <- sub(pattern = '.npar', grep('npar', names(data.16p), value = TRUE), replacement = '')

# free parameters
npar.16p <- data.16p[1, grep('npar', names(data.16p))]; npar.16p

# degrees of freedom
df.16p <- data.16p[1, grep('df', names(data.16p))]; df.16p



# descriptive statistics
measures <- c('Fmin', 'cfi', 'srmr', 'rmsea', 'tli', 'aic', 'bic')
colnames <- c('k', 'df', 'Fmin mean', 'Fmin SD', 'CFI mean', 'CFI SD', 'SRMR mean', 'SRMR SD', 'RMSEA mean', 'RMSEA SD', 'TLI mean', 'TLI SD', 'AIC mean', 'AIC SD', 'BIC mean', 'BIC SD')
res.16p <- as.data.frame(matrix(ncol = length(colnames), nrow = length(npar.16p), dimnames = list(models.16p, colnames)))

res.16p[, 1] <- as.numeric(npar.16p)
res.16p[, 2] <- as.numeric(df.16p)

for(i in 1:length(measures)){
  cdat <- data.16p[, grep(measures[i], names(data.16p))]
  res.16p[, 2*i+1] <- apply(cdat, 2, mean)
  res.16p[, 2*i+2] <- apply(cdat, 2, sd)
}
round(res.16p, 2)

# C*NML
Fmin.16p <- data.16p[, grep('Fmin', names(data.16p))]
logdetS <- data.16p$logDetS
apply(Fmin.16p, 2, function(x) logSumExp(-.5*(x + logdetS)))
round(apply(Fmin.16p, 2, function(x) logSumExp(-.5*(x + logdetS))), 2)


# cfi
cfi.16p <- data.16p[, grep('cfi', names(data.16p))]
apply(cfi.16p, 2, FUN = function(x) sum(x >= .90)/length(x))


# srmr
srmr.16p <- data.16p[, grep('srmr', names(data.16p))]
apply(srmr.16p, 2, FUN = function(x) sum(x <= .10)/length(x))
apply(srmr.16p, 2, FUN = function(x) sum(x <= .08)/length(x))
apply(srmr.16p, 2, FUN = function(x) sum(x <= .06)/length(x))


# rmsea
rmsea.16p <- data.16p[, grep('rmsea', names(data.16p))]
apply(rmsea.16p, 2, FUN = function(x) sum(x <= .08)/length(x))

# similarity of rmsea distributions
rmsea.16p <- data.16p[, grep('rmsea', names(data.16p))]
rmsea.models.16p <- names(rmsea.16p)
pairs <- combn(rmsea.models.16p, 2)
compare.rmsea.16p <- matrix(nrow=ncol(pairs), ncol = 3)
colnames(compare.rmsea.16p) <- c('pair', 'cliffs delta', 'kolmogorov-smirnov')

for(i in 1:ncol(pairs)){
  mod.a <- pairs[1,i]
  mod.b <- pairs[2,i]
  # Cliff's delta
  cliff <- cliff.delta(data.16p[,mod.a], data.16p[,mod.b])$estimate
  # Kolmogorov Smirnov
  kolm <- ks.test(data.16p[,mod.a], data.16p[,mod.b])$statistic
  compare.rmsea.16p[i,1] <- gsub(x=paste(mod.a, mod.b, sep = ' vs '), pattern='.rmsea', replacement='')
  compare.rmsea.16p[i, 2:3] <- c(cliff, kolm)
}

compare.rmsea.16p <- as.data.frame(compare.rmsea.16p)
compare.rmsea.16p

round(median(abs(as.numeric(compare.rmsea.16p$`cliffs delta`))), 2)
round(median(abs(as.numeric(compare.rmsea.16p$`kolmogorov-smirnov`))), 2)

# tli
tli.16p <- data.16p[, grep('tli', names(data.16p))]
apply(tli.16p, 2, FUN = function(x) sum(x >= .95))
apply(tli.16p, 2, FUN = function(x) sum(x >= .95)/length(x))

# similarity of TLI distributions
tli.16p <- data.16p[, grep('tli', names(data.16p))]
tli.models.16p <- names(tli.16p)
pairs <- combn(tli.models.16p, 2)
compare.tli.16p <- matrix(nrow=ncol(pairs), ncol = 3)
colnames(compare.tli.16p) <- c('pair', 'cliffs delta', 'kolmogorov-smirnov')

for(i in 1:ncol(pairs)){
  mod.a <- pairs[1,i]
  mod.b <- pairs[2,i]
  # Cliff's delta
  cliff <- cliff.delta(data.16p[,mod.a], data.16p[,mod.b])$estimate
  # Kolmogorov Smirnov
  kolm <- ks.test(data.16p[,mod.a], data.16p[,mod.b])$statistic
  compare.tli.16p[i,1] <- gsub(x=paste(mod.a, mod.b, sep = ' vs '), pattern='.tli', replacement='')
  compare.tli.16p[i, 2:3] <- c(cliff, kolm)
}

compare.tli.16p <- as.data.frame(compare.tli.16p)
compare.tli.16p

round(median(abs(as.numeric(compare.tli.16p$`cliffs delta`))), 2)
round(median(abs(as.numeric(compare.tli.16p$`kolmogorov-smirnov`))), 2)


# supplement: CNML
LL.16p <- data.16p[, grep('LL', names(data.16p))]
round(apply(LL.16p, 2, logSumExp))




# 30 manifest variables ---------------------------------------------------

# models
models.30p <- sub(pattern = '.npar', grep('npar', names(data.30p), value = TRUE), replacement = '')

# free parameters
npar.30p <- data.30p[1, grep('npar', names(data.30p))]; npar.30p

# degrees of freedom
df.30p <- data.30p[1, grep('df', names(data.30p))]; df.30p



# descriptive statistics
measures <- c('Fmin', 'cfi', 'srmr', 'rmsea', 'tli', 'aic', 'bic')
colnames <- c('k', 'df', 'Fmin mean', 'Fmin SD', 'CFI mean', 'CFI SD', 'SRMR mean', 'SRMR SD', 'RMSEA mean', 'RMSEA SD', 'TLI mean', 'TLI SD', 'AIC mean', 'AIC SD', 'BIC mean', 'BIC SD')
res.30p <- as.data.frame(matrix(ncol = length(colnames), nrow = length(npar.30p), dimnames = list(models.30p, colnames)))

res.30p[, 1] <- as.numeric(npar.30p)
res.30p[, 2] <- as.numeric(df.30p)

for(i in 1:length(measures)){
  cdat <- data.30p[, grep(measures[i], names(data.30p))]
  res.30p[, 2*i+1] <- apply(cdat, 2, mean)
  res.30p[, 2*i+2] <- apply(cdat, 2, sd)
}
round(res.30p, 2)


# 5-specifics bifactor vs. 6-specifics bifactor
cliff.delta(data.30p$bif5.Fmin, data.30p$bif6.Fmin)$estimate
ks.test(data.30p$bif5.Fmin, data.30p$bif6.Fmin)$statistic

# C*NML
Fmin.30p <- data.30p[, grep('Fmin', names(data.30p))]
logdetS <- data.30p$logDetS
apply(Fmin.30p, 2, function(x) logSumExp(-.5*(x + logdetS)))
round(apply(Fmin.30p, 2, function(x) logSumExp(-.5*(x + logdetS))), 2)



# cfi
cfi.30p <- data.30p[, grep('cfi', names(data.30p))]
apply(cfi.30p, 2, FUN = function(x) sum(x >= .90)/length(x))


# srmr
srmr.30p <- data.30p[, grep('srmr', names(data.30p))]
apply(srmr.30p, 2, FUN = function(x) sum(x <= .10)/length(x))
apply(srmr.30p, 2, FUN = function(x) sum(x <= .08)/length(x))
apply(srmr.30p, 2, FUN = function(x) sum(x <= .06)/length(x))


# rmsea
rmsea.30p <- data.30p[, grep('rmsea', names(data.30p))]
apply(rmsea.30p, 2, FUN = function(x) sum(x <= .08)/length(x))

# similarity of rmsea distributions
rmsea.30p <- data.30p[, grep('rmsea', names(data.30p))]
rmsea.models.30p <- names(rmsea.30p)
pairs <- combn(rmsea.models.30p, 2)
compare.rmsea.30p <- matrix(nrow=ncol(pairs), ncol = 3)
colnames(compare.rmsea.30p) <- c('pair', 'cliffs delta', 'kolmogorov-smirnov')

for(i in 1:ncol(pairs)){
  mod.a <- pairs[1,i]
  mod.b <- pairs[2,i]
  # Cliff's delta
  cliff <- cliff.delta(data.30p[,mod.a], data.30p[,mod.b])$estimate
  # Kolmogorov Smirnov
  kolm <- ks.test(data.30p[,mod.a], data.30p[,mod.b])$statistic
  compare.rmsea.30p[i,1] <- gsub(x=paste(mod.a, mod.b, sep = ' vs '), pattern='.rmsea', replacement='')
  compare.rmsea.30p[i, 2:3] <- c(cliff, kolm)
}

compare.rmsea.30p <- as.data.frame(compare.rmsea.30p)
compare.rmsea.30p

round(median(abs(as.numeric(compare.rmsea.30p$`cliffs delta`))), 2)
round(median(abs(as.numeric(compare.rmsea.30p$`kolmogorov-smirnov`))), 2)

# tli
tli.30p <- data.30p[, grep('tli', names(data.30p))]
apply(tli.30p, 2, FUN = function(x) sum(x >= .95))


# similarity of TLI distributions
tli.models.30p <- names(tli.30p)
pairs <- combn(tli.models.30p, 2)
compare.tli.30p <- matrix(nrow=ncol(pairs), ncol = 3)
colnames(compare.tli.30p) <- c('pair', 'cliffs delta', 'kolmogorov-smirnov')

for(i in 1:ncol(pairs)){
  mod.a <- pairs[1,i]
  mod.b <- pairs[2,i]
  # Cliff's delta
  cliff <- cliff.delta(data.30p[,mod.a], data.30p[,mod.b])$estimate
  # Kolmogorov Smirnov
  kolm <- ks.test(data.30p[,mod.a], data.30p[,mod.b])$statistic
  compare.tli.30p[i,1] <- gsub(x=paste(mod.a, mod.b, sep = ' vs '), pattern='.tli', replacement='')
  compare.tli.30p[i, 2:3] <- c(cliff, kolm)
}

compare.tli.30p <- as.data.frame(compare.tli.30p)
compare.tli.30p

round(median(abs(as.numeric(compare.tli.30p$`cliffs delta`))), 2)
round(median(abs(as.numeric(compare.tli.30p$`kolmogorov-smirnov`))), 2)


# supplement: CNML
LL.30p <- data.30p[, grep('LL', names(data.30p))]
round(apply(LL.30p, 2, logSumExp))
