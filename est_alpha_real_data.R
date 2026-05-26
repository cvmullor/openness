# Estimate power law's alpha (openness) using linear and exponential fit (10 replicates)
#
# Written by Carlos Valiente-Mullor
#
# Usage:
#       Rscript est_alpha_real_data.R [sp_name]
#
# Input file name 'gene_presence_absence_roary.{species_name}.csv' is expected by default
#
# Example:
#       Rscript est_alpha_real_data.R Enterococcus_faecium


## packages
library(pagoo)
library(matrixStats)
library(minpack.lm)
library(nlstools) # nlsMicrobio

## Functions

# Compute Delta(n) and perform permutations from pan matrix
#   * Adapted from Micropan (doi: 10.1186/s12859-015-0517-0)
#   * Default: n.perm = 100
powerlaw_perm <- function(pan.matrix, n.perm = 100){
  pan.matrix[which(pan.matrix > 0, arr.ind = T)] <- 1
  ng <- dim(pan.matrix)[1]
  nmat <- matrix(0, nrow = nrow(pan.matrix) - 1, ncol = n.perm)
  for(i in 1:n.perm){
    cm <- apply(pan.matrix[sample(nrow(pan.matrix)),], 2, cumsum)
    nmat[,i] <- rowSums((cm == 1)[2:ng,] & (cm == 0)[1:(ng-1),])
    cat(i, "/", n.perm, "\r")
  }
  return(nmat)
}

# Functions to prepare data from permutation matrix (including values=0)
prepare_medians_zero <- function(raref){
  size = dim(raref)[1]+1
  nmat_med = rowMedians(raref)
  x.med <- rep(2:size, times = 1)
  y.med <- as.numeric(nmat_med)
  
  median_curve = as.data.frame(cbind(x.med, y.med))
  colnames(median_curve) = c("N", "new_genes")
  
  return(median_curve)
}

prepare_means_zero <- function(raref){
  size = dim(raref)[1]+1
  nmat_mean = rowMeans(raref)
  x.mean <- rep(2:size, times = 1)
  y.mean <- as.numeric(nmat_mean)
  y.mean[y.mean == 0] <- NA
  
  mean_curve = as.data.frame(cbind(x.mean, y.mean))
  colnames(mean_curve) = c("N", "new_genes")
  
  return(mean_curve)
}

prepare_all_zero <- function(raref, n.perm=100){
  size = dim(raref)[1]+1
  x.all <- rep((2:size), times = n.perm)
  y.all <- as.numeric(raref)
  y.all[y.all == 0] <- NA

  all_curve = as.data.frame(cbind(x.all, y.all))
  colnames(all_curve) = c("N", "new_genes")

  return(all_curve)
}


# Functions to prepare data from permutation matrix (values=0 converted to NA)
prepare_medians <- function(raref){
  size = dim(raref)[1]+1
  nmat_med = rowMedians(raref)
  x.med <- rep(2:size, times = 1)
  y.med <- as.numeric(nmat_med)
  y.med[y.med == 0] <- NA
  
  median_curve = as.data.frame(cbind(x.med, y.med))
  colnames(median_curve) = c("N", "new_genes")
  
  return(median_curve)
}

prepare_means <- function(raref){
  size = dim(raref)[1]+1
  nmat_mean = rowMeans(raref)
  x.mean <- rep(2:size, times = 1)
  y.mean <- as.numeric(nmat_mean)
  y.mean[y.mean == 0] <- NA
  
  mean_curve = as.data.frame(cbind(x.mean, y.mean))
  colnames(mean_curve) = c("N", "new_genes")
  
  return(mean_curve)
}

prepare_all <- function(raref, n.perm=100){
  size = dim(raref)[1]+1
  x.all <- rep((2:size), times = n.perm)
  y.all <- as.numeric(raref)
  y.all[y.all == 0] <- NA

  all_curve = as.data.frame(cbind(x.all, y.all))
  colnames(all_curve) = c("N", "new_genes")

  return(all_curve)
}


# Function to estimate alpha from prepared data using linear model fit
alpha_lm <- function(df){
  # Extract data from data.frame with N + increment(n) values
  xval = df[,1]
  yval = df[,2]
  
  # Fit linear model
  fitHeaps <- try( lm(log(yval) ~ log(xval), na.action = na.exclude) )
  
  # Extract info from model
  if ( class(fitHeaps) != "try-error" ) {
    fit_res = summary(fitHeaps)
    fit_res.alpha = -1*(fit_res$coefficients[2])
    fit_res.se = fit_res$coefficients[4]
    
    # Organize results
    size = round(xval[length(xval)], 0)
    res = c(size, fit_res.alpha, fit_res.se)
  } else {
    res = rep(NA, 3)
  }
  
  # Return results
  names(res) = c("N", "alpha", "SE")
  return(res)
}


# Function to estimate alpha from prepared data using non-linear model fit (nlsLM)
alpha_nls <- function(df){
  # Extract data from data.frame with N + increment(n) values
  xval = df[,1]
  yval = df[,2]
  
  # Control parameters
  control1 <- nls.control(maxiter= 1000, minFactor= 1e-30, warnOnly= FALSE,tol=1e-05)
  
  # Fit power-law model with 'nlsLM' and extract data if no error
  fit.med = try( nlsLM(yval ~ k*xval^(-a), 
                       start = list(a = 1, k = 1),
                       control = control1,
                       na.action = na.exclude) )
  
  if ( class(fit.med) != "try-error" ) {
    # Model parameters
    fit_res       = summary(fit.med)
    fit_res.alpha = fit_res$parameters[1]
    fit_res.se    = fit_res$parameters[3]
  
    # Organize results
    size = round(xval[length(xval)], 0)
    res = c(size, fit_res.alpha, fit_res.se)
  } else {
    res = rep(NA, 3)
  }
  
  # Return results
  names(res) = c("N", "alpha", "SE")
  return(res)
}



########## START ##########

## Arguments
args = commandArgs(trailingOnly=TRUE)
spp = args[1] # species name ("Genus_epithet")

## Paths
inpath = "path/to/input/dir"   # set input dir
outpath= "path/to/output/dir"  # set output dir

## Read input / pan-matrix
myfile_r = paste0(inpath, "gene_presence_absence_roary.", spp, ".csv") # Roary-formatted PAM from Panaroo
pg = roary_2_pagoo(gene_presence_absence_csv = myfile_r) # Panaroo table to pagoo object

## pan.matrix
pan.matrix = pg$pan_matrix # extract pan-matrix

# replicates (=10)
alpha.nls  = rep(0, 10)
alpha.lm   = rep(0, 10)
alpha.nls2 = rep(0, 10)


#############################################
#   Estimate alpha using ALL permutations   #
#############################################

for (i in seq(1, 10)) {
  # Perform rarefaction + permutations
  nmat = powerlaw_perm(pan.matrix)
  
  ## Method A: nls fit (power-law), including values=0
  data.nls = prepare_all_zero(nmat)
  res.nls = alpha_nls(data.nls)
  
  alpha.nls[i] = res.nls[2]
  
  ## Method B: lm fit (linearized power-law), excluding values=0 (NAs)
  data.lm = prepare_all(nmat)
  res.lm = alpha_lm(data.lm)
  
  alpha.lm[i] = res.lm[2]
  
  ## Method C: nls fit (power-law), excluding values=0 (NAs)
  res.nls2 = alpha_nls(data.lm)
  
  alpha.nls2[i] = res.nls2[2]
}

# Summarize data
sum.nls  = round(c(summary(alpha.nls), var(alpha.nls), sd(alpha.nls)), 4)
sum.nls2 = round(c(summary(alpha.nls2), var(alpha.nls2), sd(alpha.nls2)), 4)
sum.lm   = round(c(summary(alpha.lm), var(alpha.lm), sd(alpha.lm)), 4)

method.v = c("nlsLM", "nlsLM", "lm")
zeros.v  = c("include", "exclude", "exclude")

df.temp = rbind(sum.nls, sum.nls2, sum.lm)
df.res = cbind(method.v, zeros.v, df.temp)

header = c("Method", "0_values", "Min", "Q1", "Median", "Mean", "Q3", "Max", "Var", "SD")
colnames(df.res) = header

# Write summaries
mytable = paste0(outpath, spp, ".all.rep_x10_sum.csv")
write.csv(df.res, mytable, row.names = F)

# Write raw data.
df.alpha = rbind(alpha.nls, alpha.nls2, alpha.lm)
rownames(df.alpha) = c("nlsLM_in", "nlsLM_ex", "lm_ex")
colnames(df.alpha) = c("r1", "r2", "r3", "r4", "r5", "r6", "r7", "r8", "r9", "r10")

mytable = paste0(outpath, spp, ".all.rep_x10_data.csv")
write.csv(df.alpha, mytable, row.names = T)


#################################################
#   Estimate alpha using MEAN of permutations   #
#################################################

for (i in seq(1, 10)) {
  # Rarefaction + permutation
  nmat = powerlaw_perm(pan.matrix)
  
  ## Method A: nls fit (power-law), including values=0
  data.nls = prepare_means_zero(nmat)
  res.nls = alpha_nls(data.nls)
  
  alpha.nls[i] = res.nls[2]
  
  ## Method B: lm fit (linearized power-law), excluding values=0 (NAs)
  data.lm = prepare_means(nmat)
  res.lm = alpha_lm(data.lm)
  
  alpha.lm[i] = res.lm[2]
  
  ## Method C: nls fit (power-law), excluding values=0 (NAs)
  res.nls2 = alpha_nls(data.lm)
  
  alpha.nls2[i] = res.nls2[2]
}

# Summarize data
sum.nls  = round(c(summary(alpha.nls), var(alpha.nls), sd(alpha.nls)), 4)
sum.nls2 = round(c(summary(alpha.nls2), var(alpha.nls2), sd(alpha.nls2)), 4)
sum.lm   = round(c(summary(alpha.lm), var(alpha.lm), sd(alpha.lm)), 4)

method.v = c("nlsLM", "nlsLM", "lm")
zeros.v  = c("include", "exclude", "exclude")

df.temp = rbind(sum.nls, sum.nls2, sum.lm)
df.res = cbind(method.v, zeros.v, df.temp)

header = c("Method", "0_values", "Min", "Q1", "Median", "Mean", "Q3", "Max", "Var", "SD")
colnames(df.res) = header

# Write summaries
mytable = paste0(outpath, spp, ".mean.rep_x10_sum.csv")
write.csv(df.res, mytable, row.names = F)

# Write raw data
df.alpha = rbind(alpha.nls, alpha.nls2, alpha.lm)
rownames(df.alpha) = c("nlsLM_in", "nlsLM_ex", "lm_ex")
colnames(df.alpha) = c("r1", "r2", "r3", "r4", "r5", "r6", "r7", "r8", "r9", "r10")

mytable = paste0(outpath, spp, ".mean.rep_x10_data.csv")
write.csv(df.alpha, mytable, row.names = T)


###################################################
#   Estimate alpha using MEDIAN of permutations   #
###################################################

for (i in seq(1, 10)) {
  # Rarefaction + permutation
  nmat = powerlaw_perm(pan.matrix)
  
  ## Method A: nls fit (power-law), including values=0
  data.nls = prepare_medians_zero(nmat)
  res.nls = alpha_nls(data.nls)
  
  alpha.nls[i] = res.nls[2]
  
  ## Method B: lm fit (linearized power-law), excluding values=0 (NAs)
  data.lm = prepare_medians(nmat)
  res.lm = alpha_lm(data.lm)
  
  alpha.lm[i] = res.lm[2]
  
  ## Method C: nls fit (power-law), excluding values=0 (NAs)
  res.nls2 = alpha_nls(data.lm)
  
  alpha.nls2[i] = res.nls2[2]
}

# Summarize data
sum.nls  = round(c(summary(alpha.nls), var(alpha.nls), sd(alpha.nls)), 4)
sum.nls2 = round(c(summary(alpha.nls2), var(alpha.nls2), sd(alpha.nls2)), 4)
sum.lm   = round(c(summary(alpha.lm), var(alpha.lm), sd(alpha.lm)), 4)

method.v = c("nlsLM", "nlsLM", "lm")
zeros.v  = c("include", "exclude", "exclude")

df.temp = rbind(sum.nls, sum.nls2, sum.lm)
df.res = cbind(method.v, zeros.v, df.temp)

header = c("Method", "0_values", "Min", "Q1", "Median", "Mean", "Q3", "Max", "Var", "SD")
colnames(df.res) = header

# Write summaries
mytable = paste0(outpath, spp, ".med.rep_x10_sum.csv")
write.csv(df.res, mytable, row.names = F)

# Write raw data.
df.alpha = rbind(alpha.nls, alpha.nls2, alpha.lm)
rownames(df.alpha) = c("nlsLM_in", "nlsLM_ex", "lm_ex")
colnames(df.alpha) = c("r1", "r2", "r3", "r4", "r5", "r6", "r7", "r8", "r9", "r10")

mytable = paste0(outpath, spp, ".med.rep_x10_data.csv")
write.csv(df.alpha, mytable, row.names = T)

