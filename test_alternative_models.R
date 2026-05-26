# Evaluation of alternative models for pangenome curves
#
# Written by Carlos Valiente-Mullor
#
# Usage:
#       Rscript test_alternative_models.R [sp_name]
#
# Input file name 'gene_presence_absence_roary.{species_name}.csv' is expected by default
#
# Example:
#       Rscript test_alternative_models.R Enterococcus_faecium


## packages
library(pagoo)
library(matrixStats)
library(minpack.lm)
library(nlstools) # nlsMicrobio
library(AICcmodavg)

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

########################################
#   Fit alternative models (MEDIANS)   #
########################################
# * 5 replicates

for (i in seq(1, 5)) {
  ## Perform rarefaction + permutation
  nmat = powerlaw_perm(pan.matrix)
  
  data.med = prepare_medians_zero(nmat)
  
  ## nls fit
  x.med = data.med[,1]
  y.med = data.med[,2]
  
  ## Create objects to save model comparison
  fit.models = list()
  model.names = list()
  idx = 0
  
  # nlsLM control parameters
  control1 <- nls.control(maxiter= 1000, minFactor= 1e-30, warnOnly= FALSE,tol=1e-05)
  
  ## (1) 2-phase exponential decay
  # y = c + k1*exp(-t1*x) + k2*exp(-t2*x)
  
  fit.med = try( nlsLM(y.med ~ 0 + k1*exp(-t1*x.med) + k2*exp(-t2*x.med), 
                       start = list(k1=1, k2=1, t1=0.1, t2=0.1 ),
                       control = control1) )
  
  if ( class(fit.med) != "try-error" ) {
    # Model parameters
    nls_med   = summary(fit.med)
    constant1 = nls_med$parameters[1]
    constant2 = nls_med$parameters[2]
    tau1      = nls_med$parameters[3]
    tau2      = nls_med$parameters[4]
    
    # Add model to list
    idx = idx + 1
    fit.models[[idx]] = fit.med
    model.names[[idx]] = "M1"
    
    # Plot data + fit
    fitted_vals = fitted(fit.med)
    data.fit = as.data.frame(cbind(x.med, fitted_vals))
    
    pr1 = as.data.frame(cbind(y.med, round(fitted_vals, 3)))
    
    pr = pr1
    pr[, 3] = round(pr[,1]-pr[,2], 3)
    pr[, 4] = seq(2, (length(y.med)+1), by=1)
    colnames(pr) = c("real", "fitted", "diff", "index")
    
    
    # plot residuals
    resid = nlsResiduals(fit.med)
    
    # Save plots
    mdln = "model1"
    
    pl1 = ggplot(data.med, aes(x=N, y=new_genes)) + geom_point() + geom_line(data = data.fit, aes(x=x.med, y=fitted_vals), color="red") + ylim(0, max(max(data.fit$fitted_vals+0.5), max(y.med)))
    myplot = paste0(outpath, spp, ".fit_line_", mdln, ".r", i, ".png")
    ggsave(myplot, plot=pl1, width=7.22, height=7.22)
    
    pl2 = ggplot(data.med, aes(x=N, y=new_genes)) + geom_point() + geom_point(data = data.fit, aes(x=x.med, y=fitted_vals), color="red", shape=1) + ylim(0, max(max(data.fit$fitted_vals+0.5), max(y.med)))
    myplot = paste0(outpath, spp, ".fit_point_", mdln, ".r", i, ".png")
    ggsave(myplot, plot=pl2, width=7.22, height=7.22)
    
    pl3 = ggplot(pr, aes(x=real , y=fitted)) + geom_point(shape=1)
    myplot = paste0(outpath, spp, ".fit_diff1_", mdln, ".r", i, ".png")
    ggsave(myplot, plot=pl3, width=7.22, height=7.22)
    
    pl4 = ggplot(pr, aes(x=index, y=diff)) + geom_point(shape=1)
    myplot = paste0(outpath, spp, ".fit_diff2_", mdln, ".r", i, ".png")
    ggsave(myplot, plot=pl4, width=7.22, height=7.22)
    
    myplot = paste0(outpath, spp, ".fit_resid_", mdln, ".r", i, ".png")
    png(myplot, units="px", width = 693.12, height = 693.12) # 7.22 inches = 693.12 px
    plot(resid)
    dev.off()
    
    # Save data
    mytable = paste0(outpath, spp, ".fit_data_", mdln, ".r", i, ".csv")
    write.csv(pr, mytable, row.names = F)
    
    mytable = paste0(outpath, spp, ".fit_param_", mdln, ".r", i, ".csv")
    write.csv(nls_med$parameters, mytable, row.names = T)
  }
  
  
  ## (2) 2-phase power-law 
  # y = c + k1*x^-a1 + k2*x^-a2
  
  # fit
  fit.med = try( nlsLM(y.med ~ 0 + k1*(x.med^-a1) + k2*(x.med^-a2), 
                       start = list(k1=1, k2=1, a1=0.1, a2=0.1 ),
                       control = control1) )

  if ( class(fit.med) != "try-error" ) {
    # Model parameters
    nls_med   = summary(fit.med)
    constant1 = nls_med$parameters[1]
    constant2 = nls_med$parameters[2]
    alpha1    = nls_med$parameters[3]
    alpha2    = nls_med$parameters[4]
    
    # Add model to list
    idx = idx + 1
    fit.models[[idx]] = fit.med
    model.names[[idx]] = "M2"
    
    # Plot data + fit
    fitted_vals = fitted(fit.med)
    data.fit = as.data.frame(cbind(x.med, fitted_vals))
    
    pr2 = as.data.frame(cbind(y.med, round(fitted_vals, 3)))
    
    pr = pr2
    pr[, 3] = round(pr[,1]-pr[,2], 3)
    pr[, 4] = seq(2, (length(y.med)+1), by=1)
    colnames(pr) = c("real", "fitted", "diff", "index")
    
    # plot residuals
    resid = nlsResiduals(fit.med)
    
    # Save plots
    mdln = "model2"
    
    pl1 = ggplot(data.med, aes(x=N, y=new_genes)) + geom_point() + geom_line(data = data.fit, aes(x=x.med, y=fitted_vals), color="red") + ylim(0, max(max(data.fit$fitted_vals+0.5), max(y.med)))
    myplot = paste0(outpath, spp, ".fit_line_", mdln, ".r", i, ".png")
    ggsave(myplot, plot=pl1, width=7.22, height=7.22)
    
    pl2 = ggplot(data.med, aes(x=N, y=new_genes)) + geom_point() + geom_point(data = data.fit, aes(x=x.med, y=fitted_vals), color="red", shape=1) + ylim(0, max(max(data.fit$fitted_vals+0.5), max(y.med)))
    myplot = paste0(outpath, spp, ".fit_point_", mdln, ".r", i, ".png")
    ggsave(myplot, plot=pl2, width=7.22, height=7.22)
    
    pl3 = ggplot(pr, aes(x=real , y=fitted)) + geom_point(shape=1)
    myplot = paste0(outpath, spp, ".fit_diff1_", mdln, ".r", i, ".png")
    ggsave(myplot, plot=pl3, width=7.22, height=7.22)
    
    pl4 = ggplot(pr, aes(x=index, y=diff)) + geom_point(shape=1)
    myplot = paste0(outpath, spp, ".fit_diff2_", mdln, ".r", i, ".png")
    ggsave(myplot, plot=pl4, width=7.22, height=7.22)
    
    myplot = paste0(outpath, spp, ".fit_resid_", mdln, ".r", i, ".png")
    png(myplot, units="px", width = 693.12, height = 693.12) # 7.22 inches = 693.12 px
    plot(resid)
    dev.off()
    
    # Save data
    mytable = paste0(outpath, spp, ".fit_data_", mdln, ".r", i, ".csv")
    write.csv(pr, mytable, row.names = F)
    
    mytable = paste0(outpath, spp, ".fit_param_", mdln, ".r", i, ".csv")
    write.csv(nls_med$parameters, mytable, row.names = T)
  }
  
  
  ## (3) 2-phase pl + exp.decay
  # y = c + r1*x^-a + r2*exp(-k*x)
  
  # fit
  fit.med = try( nlsLM(y.med ~ 0 + k1*(x.med^-a) + k2*exp(-t*x.med), 
                       start = list(k1=1, k2=1, a=1, t=0.1 ),
                       control = control1) )
  
  if ( class(fit.med) != "try-error" ) {
    # Model parameters
    nls_med   = summary(fit.med)
    constant1 = nls_med$parameters[1]
    constant2 = nls_med$parameters[2]
    alpha     = nls_med$parameters[3]
    tau       = nls_med$parameters[4]
    
    # Add model to list
    idx = idx + 1
    fit.models[[idx]] = fit.med
    model.names[[idx]] = "M3"
    
    # Plot data + fit
    fitted_vals = fitted(fit.med)
    data.fit = as.data.frame(cbind(x.med, fitted_vals))
    
    pr3 = as.data.frame(cbind(y.med, round(fitted_vals, 3)))
    
    pr = pr3
    pr[, 3] = round(pr[,1]-pr[,2], 3)
    pr[, 4] = seq(2, (length(y.med)+1), by=1)
    colnames(pr) = c("real", "fitted", "diff", "index")
    
    # plot residuals
    resid = nlsResiduals(fit.med)
    
    # Save plots
    mdln = "model3"
    
    pl1 = ggplot(data.med, aes(x=N, y=new_genes)) + geom_point() + geom_line(data = data.fit, aes(x=x.med, y=fitted_vals), color="red") + ylim(0, max(max(data.fit$fitted_vals+0.5), max(y.med)))
    myplot = paste0(outpath, spp, ".fit_line_", mdln, ".r", i, ".png")
    ggsave(myplot, plot=pl1, width=7.22, height=7.22)
    
    pl2 = ggplot(data.med, aes(x=N, y=new_genes)) + geom_point() + geom_point(data = data.fit, aes(x=x.med, y=fitted_vals), color="red", shape=1) + ylim(0, max(max(data.fit$fitted_vals+0.5), max(y.med)))
    myplot = paste0(outpath, spp, ".fit_point_", mdln, ".r", i, ".png")
    ggsave(myplot, plot=pl2, width=7.22, height=7.22)
    
    pl3 = ggplot(pr, aes(x=real , y=fitted)) + geom_point(shape=1)
    myplot = paste0(outpath, spp, ".fit_diff1_", mdln, ".r", i, ".png")
    ggsave(myplot, plot=pl3, width=7.22, height=7.22)
    
    pl4 = ggplot(pr, aes(x=index, y=diff)) + geom_point(shape=1)
    myplot = paste0(outpath, spp, ".fit_diff2_", mdln, ".r", i, ".png")
    ggsave(myplot, plot=pl4, width=7.22, height=7.22)
    
    myplot = paste0(outpath, spp, ".fit_resid_", mdln, ".r", i, ".png")
    png(myplot, units="px", width = 693.12, height = 693.12) # 7.22 inches = 693.12 px
    plot(resid)
    dev.off()
    
    # Save data
    mytable = paste0(outpath, spp, ".fit_data_", mdln, ".r", i, ".csv")
    write.csv(pr, mytable, row.names = F)
    
    mytable = paste0(outpath, spp, ".fit_param_", mdln, ".r", i, ".csv")
    write.csv(nls_med$parameters, mytable, row.names = T)
  }
  
  
  ## (4) "double" power-law
  # y = c + k*x^(-c*x^b)
  
  fit.med = try( nlsLM(y.med ~  k * x.med^(-c * x.med^b),
                       start = list(k=1, c=1, b = -0.1),
                       control = control1) )
  
  if ( class(fit.med) != "try-error" ) {
    # Model parameters
    nls_med     = summary(fit.med)
    constant_k  = nls_med$parameters[1]
    constant_c  = nls_med$parameters[2]
    beta        = nls_med$parameters[3]
    
    # Add model to list
    idx = idx + 1
    fit.models[[idx]] = fit.med
    model.names[[idx]] = "M4"
    
    # Plot data + fit
    fitted_vals = fitted(fit.med)
    data.fit = as.data.frame(cbind(x.med, fitted_vals))
    
    pr4 = as.data.frame(cbind(y.med, round(fitted_vals, 3)))
    
    pr = pr4
    pr[, 3] = round(pr[,1]-pr[,2], 3)
    pr[, 4] = seq(2, (length(y.med)+1), by=1)
    colnames(pr) = c("real", "fitted", "diff", "index")
    
    # plot residuals
    resid = nlsResiduals(fit.med)
    
    # Save plots
    mdln = "model4"
    
    pl1 = ggplot(data.med, aes(x=N, y=new_genes)) + geom_point() + geom_line(data = data.fit, aes(x=x.med, y=fitted_vals), color="red") + ylim(0, max(max(data.fit$fitted_vals+0.5), max(y.med)))
    myplot = paste0(outpath, spp, ".fit_line_", mdln, ".r", i, ".png")
    ggsave(myplot, plot=pl1, width=7.22, height=7.22)
    
    pl2 = ggplot(data.med, aes(x=N, y=new_genes)) + geom_point() + geom_point(data = data.fit, aes(x=x.med, y=fitted_vals), color="red", shape=1) + ylim(0, max(max(data.fit$fitted_vals+0.5), max(y.med)))
    myplot = paste0(outpath, spp, ".fit_point_", mdln, ".r", i, ".png")
    ggsave(myplot, plot=pl2, width=7.22, height=7.22)
    
    pl3 = ggplot(pr, aes(x=real , y=fitted)) + geom_point(shape=1)
    myplot = paste0(outpath, spp, ".fit_diff1_", mdln, ".r", i, ".png")
    ggsave(myplot, plot=pl3, width=7.22, height=7.22)
    
    pl4 = ggplot(pr, aes(x=index, y=diff)) + geom_point(shape=1)
    myplot = paste0(outpath, spp, ".fit_diff2_", mdln, ".r", i, ".png")
    ggsave(myplot, plot=pl4, width=7.22, height=7.22)
    
    myplot = paste0(outpath, spp, ".fit_resid_", mdln, ".r", i, ".png")
    png(myplot, units="px", width = 693.12, height = 693.12) # 7.22 inches = 693.12 px
    plot(resid)
    dev.off()
    
    # Save data
    mytable = paste0(outpath, spp, ".fit_data_", mdln, ".r", i, ".csv")
    write.csv(pr, mytable, row.names = F)
    
    mytable = paste0(outpath, spp, ".fit_param_", mdln, ".r", i, ".csv")
    write.csv(nls_med$parameters, mytable, row.names = T)
  }
  
  
  ## (5) power-law
  fit.med = try( nlsLM(y.med ~ k*x.med^(-a), 
                       start = list(a = 1, k = 1),
                       control = control1) )
  
  if ( class(fit.med) != "try-error" ) {
    # Model parameters
    nls_med   = summary(fit.med)
    alpha     = nls_med$parameters[1]
    constant  = nls_med$parameters[2]
    
    # Add model to list
    idx = idx + 1
    fit.models[[idx]] = fit.med
    model.names[[idx]] = "M5"
    
    # Plot data + fit
    fitted_vals = fitted(fit.med)
    data.fit = as.data.frame(cbind(x.med, fitted_vals))
    
    pr5 = as.data.frame(cbind(y.med, round(fitted_vals, 3)))
    
    pr = pr5
    pr[, 3] = round(pr[,1]-pr[,2], 3)
    pr[, 4] = seq(2, (length(y.med)+1), by=1)
    colnames(pr) = c("real", "fitted", "diff", "index")
    
    # plot residuals
    resid = nlsResiduals(fit.med)
    
    # Save plots
    mdln = "model5"
    
    pl1 = ggplot(data.med, aes(x=N, y=new_genes)) + geom_point() + geom_line(data = data.fit, aes(x=x.med, y=fitted_vals), color="red") + ylim(0, max(max(data.fit$fitted_vals+0.5), max(y.med)))
    myplot = paste0(outpath, spp, ".fit_line_", mdln, ".r", i, ".png")
    ggsave(myplot, plot=pl1, width=7.22, height=7.22)
    
    pl2 = ggplot(data.med, aes(x=N, y=new_genes)) + geom_point() + geom_point(data = data.fit, aes(x=x.med, y=fitted_vals), color="red", shape=1) + ylim(0, max(max(data.fit$fitted_vals+0.5), max(y.med)))
    myplot = paste0(outpath, spp, ".fit_point_", mdln, ".r", i, ".png")
    ggsave(myplot, plot=pl2, width=7.22, height=7.22)
    
    pl3 = ggplot(pr, aes(x=real , y=fitted)) + geom_point(shape=1)
    myplot = paste0(outpath, spp, ".fit_diff1_", mdln, ".r", i, ".png")
    ggsave(myplot, plot=pl3, width=7.22, height=7.22)
    
    pl4 = ggplot(pr, aes(x=index, y=diff)) + geom_point(shape=1)
    myplot = paste0(outpath, spp, ".fit_diff2_", mdln, ".r", i, ".png")
    ggsave(myplot, plot=pl4, width=7.22, height=7.22)
    
    myplot = paste0(outpath, spp, ".fit_resid_", mdln, ".r", i, ".png")
    png(myplot, units="px", width = 693.12, height = 693.12) # 7.22 inches = 693.12 px
    plot(resid)
    dev.off()
    
    # Save data
    mytable = paste0(outpath, spp, ".fit_data_", mdln, ".r", i, ".csv")
    write.csv(pr, mytable, row.names = F)
    
    mytable = paste0(outpath, spp, ".fit_param_", mdln, ".r", i, ".csv")
    write.csv(nls_med$parameters, mytable, row.names = T)
  }
  
  
  ## (6) Exponential decay
  fit.med = try( nlsLM(y.med ~ 0 + k*exp(-x.med*tau), 
                       start   = list(k=1, tau=1),
                       control = control1) ) 
  
  if ( class(fit.med) != "try-error" ) {
    # Model parameters
    nls_med = summary(fit.med)
    tan.omega = 0 # assuming tg(omega)=0
    constant  = nls_med$parameters[1]
    tau.model = nls_med$parameters[2]
    tau.orig  = 1/nls_med$parameters[2]
    
    # Add model to list
    idx = idx + 1
    fit.models[[idx]] = fit.med
    model.names[[idx]] = "M6"
    
    # Plot data + fit
    fitted_vals = fitted(fit.med)
    data.fit = as.data.frame(cbind(x.med, fitted_vals))
    
    pr6 = as.data.frame(cbind(y.med, round(fitted_vals, 1)))
    
    pr = pr6
    pr[, 3] = round(pr[,1]-pr[,2], 3)
    pr[, 4] = seq(2, (length(y.med)+1), by=1)
    colnames(pr) = c("real", "fitted", "diff", "index")
    ggplot(pr, aes(x=real , y=fitted)) + geom_point(shape=1)
    ggplot(pr, aes(x=index, y=diff)) + geom_point(shape=1)
    
    
    # plot residuals
    resid = nlsResiduals(fit.med)
    
    # Save plots
    mdln = "model6"
    
    pl1 = ggplot(data.med, aes(x=N, y=new_genes)) + geom_point() + geom_line(data = data.fit, aes(x=x.med, y=fitted_vals), color="red") + ylim(0, max(max(data.fit$fitted_vals+0.5), max(y.med)))
    myplot = paste0(outpath, spp, ".fit_line_", mdln, ".r", i, ".png")
    ggsave(myplot, plot=pl1, width=7.22, height=7.22)
    
    pl2 = ggplot(data.med, aes(x=N, y=new_genes)) + geom_point() + geom_point(data = data.fit, aes(x=x.med, y=fitted_vals), color="red", shape=1) + ylim(0, max(max(data.fit$fitted_vals+0.5), max(y.med)))
    myplot = paste0(outpath, spp, ".fit_point_", mdln, ".r", i, ".png")
    ggsave(myplot, plot=pl2, width=7.22, height=7.22)
    
    pl3 = ggplot(pr, aes(x=real , y=fitted)) + geom_point(shape=1)
    myplot = paste0(outpath, spp, ".fit_diff1_", mdln, ".r", i, ".png")
    ggsave(myplot, plot=pl3, width=7.22, height=7.22)
    
    pl4 = ggplot(pr, aes(x=index, y=diff)) + geom_point(shape=1)
    myplot = paste0(outpath, spp, ".fit_diff2_", mdln, ".r", i, ".png")
    ggsave(myplot, plot=pl4, width=7.22, height=7.22)
    
    myplot = paste0(outpath, spp, ".fit_resid_", mdln, ".r", i, ".png")
    png(myplot, units="px", width = 693.12, height = 693.12) # 7.22 inches = 693.12 px
    plot(resid)
    dev.off()
    
    # Save data
    mytable = paste0(outpath, spp, ".fit_data_", mdln, ".r", i, ".csv")
    write.csv(pr, mytable, row.names = F)
    
    mytable = paste0(outpath, spp, ".fit_param_", mdln, ".r", i, ".csv")
    write.csv(nls_med$parameters, mytable, row.names = T)
  }
  
  
  ## Model selection (AIC)
  
  # Test
  aic.res = aictab(cand.set = fit.models, modnames = unlist(model.names))
  aic.res
  
  # Save test results
  mytable = paste0(outpath, spp, ".AIC_exp", ".r", i, ".csv")
  write.csv(aic.res, mytable, row.names = F)

}


########################################
#   Fit alternative models (MEANS)   #
########################################
# * 1 replicates

outpath2 = "path/to/output_2/dir"   # different dir to save output using means

data.med = prepare_means_zero(nmat) # prepare data using means

## nls fit
x.med = data.med[,1]
y.med = data.med[,2]


## Model comparison
fit.models = list()
model.names = list()
idx = 0

control1 <- nls.control(maxiter= 1000, minFactor= 1e-30, warnOnly= FALSE,tol=1e-05)

## (1) 2-phase exponential decay
# y = c + k1*exp(-t1*x) + k2*exp(-t2*x)

fit.med = try( nlsLM(y.med ~ 0 + k1*exp(-t1*x.med) + k2*exp(-t2*x.med), 
                     start = list(k1=1, k2=1, t1=0.1, t2=0.1 ),
                     control = control1) )

if ( class(fit.med) != "try-error" ) {
  # Model parameters
  nls_med   = summary(fit.med)
  constant1 = nls_med$parameters[1]
  constant2 = nls_med$parameters[2]
  tau1      = nls_med$parameters[3]
  tau2      = nls_med$parameters[4]
  
  # Add model to list
  idx = idx + 1
  fit.models[[idx]] = fit.med
  model.names[[idx]] = "M1"
  
  # Plot data + fit
  fitted_vals = fitted(fit.med)
  data.fit = as.data.frame(cbind(x.med, fitted_vals))
  
  pr1 = as.data.frame(cbind(y.med, round(fitted_vals, 3)))
  
  pr = pr1
  pr[, 3] = round(pr[,1]-pr[,2], 3)
  pr[, 4] = seq(2, (length(y.med)+1), by=1)
  colnames(pr) = c("real", "fitted", "diff", "index")
  
  
  # plot residuals
  resid = nlsResiduals(fit.med)
  
  # Save plots
  mdln = "model1"
  
  pl1 = ggplot(data.med, aes(x=N, y=new_genes)) + geom_point() + geom_line(data = data.fit, aes(x=x.med, y=fitted_vals), color="red") + ylim(0, max(max(data.fit$fitted_vals+0.5), max(y.med)))
  myplot = paste0(outpath2, spp, ".fit_line_", mdln, ".png")
  ggsave(myplot, plot=pl1, width=7.22, height=7.22)
  
  pl2 = ggplot(data.med, aes(x=N, y=new_genes)) + geom_point() + geom_point(data = data.fit, aes(x=x.med, y=fitted_vals), color="red", shape=1) + ylim(0, max(max(data.fit$fitted_vals+0.5), max(y.med)))
  myplot = paste0(outpath2, spp, ".fit_point_", mdln, ".png")
  ggsave(myplot, plot=pl2, width=7.22, height=7.22)
  
  pl3 = ggplot(pr, aes(x=real , y=fitted)) + geom_point(shape=1)
  myplot = paste0(outpath2, spp, ".fit_diff1_", mdln, ".png")
  ggsave(myplot, plot=pl3, width=7.22, height=7.22)
  
  pl4 = ggplot(pr, aes(x=index, y=diff)) + geom_point(shape=1)
  myplot = paste0(outpath2, spp, ".fit_diff2_", mdln, ".png")
  ggsave(myplot, plot=pl4, width=7.22, height=7.22)
  
  myplot = paste0(outpath2, spp, ".fit_resid_", mdln, ".png")
  png(myplot, units="px", width = 693.12, height = 693.12) # 7.22 inches = 693.12 px
  plot(resid)
  dev.off()
  
  # Save data
  mytable = paste0(outpath2, spp, ".fit_data_", mdln, ".csv")
  write.csv(pr, mytable, row.names = F)
  
  mytable = paste0(outpath2, spp, ".fit_param_", mdln, ".csv")
  write.csv(nls_med$parameters, mytable, row.names = T)
}


## (2) 2-phase power-law 
# y = c + k1*x^-a1 + k2*x^-a2

# fit
fit.med = try( nlsLM(y.med ~ 0 + k1*(x.med^-a1) + k2*(x.med^-a2), 
                     start = list(k1=1, k2=1, a1=0.1, a2=0.1 ),
                     control = control1) )

if ( class(fit.med) != "try-error" ) {
  # Model parameters
  nls_med   = summary(fit.med)
  constant1 = nls_med$parameters[1]
  constant2 = nls_med$parameters[2]
  alpha1    = nls_med$parameters[3]
  alpha2    = nls_med$parameters[4]
  
  # Add model to list
  idx = idx + 1
  fit.models[[idx]] = fit.med
  model.names[[idx]] = "M2"
  
  # Plot data + fit
  fitted_vals = fitted(fit.med)
  data.fit = as.data.frame(cbind(x.med, fitted_vals))
  
  pr2 = as.data.frame(cbind(y.med, round(fitted_vals, 3)))
  
  pr = pr2
  pr[, 3] = round(pr[,1]-pr[,2], 3)
  pr[, 4] = seq(2, (length(y.med)+1), by=1)
  colnames(pr) = c("real", "fitted", "diff", "index")
  
  # plot residuals
  resid = nlsResiduals(fit.med)
  
  # Save plots
  mdln = "model2"
  
  pl1 = ggplot(data.med, aes(x=N, y=new_genes)) + geom_point() + geom_line(data = data.fit, aes(x=x.med, y=fitted_vals), color="red") + ylim(0, max(max(data.fit$fitted_vals+0.5), max(y.med)))
  myplot = paste0(outpath2, spp, ".fit_line_", mdln, ".png")
  ggsave(myplot, plot=pl1, width=7.22, height=7.22)
  
  pl2 = ggplot(data.med, aes(x=N, y=new_genes)) + geom_point() + geom_point(data = data.fit, aes(x=x.med, y=fitted_vals), color="red", shape=1) + ylim(0, max(max(data.fit$fitted_vals+0.5), max(y.med)))
  myplot = paste0(outpath2, spp, ".fit_point_", mdln, ".png")
  ggsave(myplot, plot=pl2, width=7.22, height=7.22)
  
  pl3 = ggplot(pr, aes(x=real , y=fitted)) + geom_point(shape=1)
  myplot = paste0(outpath2, spp, ".fit_diff1_", mdln, ".png")
  ggsave(myplot, plot=pl3, width=7.22, height=7.22)
  
  pl4 = ggplot(pr, aes(x=index, y=diff)) + geom_point(shape=1)
  myplot = paste0(outpath2, spp, ".fit_diff2_", mdln, ".png")
  ggsave(myplot, plot=pl4, width=7.22, height=7.22)
  
  myplot = paste0(outpath2, spp, ".fit_resid_", mdln, ".png")
  png(myplot, units="px", width = 693.12, height = 693.12) # 7.22 inches = 693.12 px
  plot(resid)
  dev.off()
  
  # Save data
  mytable = paste0(outpath2, spp, ".fit_data_", mdln, ".csv")
  write.csv(pr, mytable, row.names = F)
  
  mytable = paste0(outpath2, spp, ".fit_param_", mdln, ".csv")
  write.csv(nls_med$parameters, mytable, row.names = T)
}


## (3) 2-phase pl + exp.decay
# y = c + r1*x^-a + r2*exp(-k*x)

# fit
fit.med = try( nlsLM(y.med ~ 0 + k1*(x.med^-a) + k2*exp(-t*x.med), 
                     start = list(k1=1, k2=1, a=1, t=0.1 ),
                     control = control1) )

if ( class(fit.med) != "try-error" ) {
  # Model parameters
  nls_med   = summary(fit.med)
  constant1 = nls_med$parameters[1]
  constant2 = nls_med$parameters[2]
  alpha     = nls_med$parameters[3]
  tau       = nls_med$parameters[4]
  
  # Add model to list
  idx = idx + 1
  fit.models[[idx]] = fit.med
  model.names[[idx]] = "M3"
  
  # Plot data + fit
  fitted_vals = fitted(fit.med)
  data.fit = as.data.frame(cbind(x.med, fitted_vals))
  
  pr3 = as.data.frame(cbind(y.med, round(fitted_vals, 3)))
  
  pr = pr3
  pr[, 3] = round(pr[,1]-pr[,2], 3)
  pr[, 4] = seq(2, (length(y.med)+1), by=1)
  colnames(pr) = c("real", "fitted", "diff", "index")
  
  # plot residuals
  resid = nlsResiduals(fit.med)
  
  # Save plots
  mdln = "model3"
  
  pl1 = ggplot(data.med, aes(x=N, y=new_genes)) + geom_point() + geom_line(data = data.fit, aes(x=x.med, y=fitted_vals), color="red") + ylim(0, max(max(data.fit$fitted_vals+0.5), max(y.med)))
  myplot = paste0(outpath2, spp, ".fit_line_", mdln, ".png")
  ggsave(myplot, plot=pl1, width=7.22, height=7.22)
  
  pl2 = ggplot(data.med, aes(x=N, y=new_genes)) + geom_point() + geom_point(data = data.fit, aes(x=x.med, y=fitted_vals), color="red", shape=1) + ylim(0, max(max(data.fit$fitted_vals+0.5), max(y.med)))
  myplot = paste0(outpath2, spp, ".fit_point_", mdln, ".png")
  ggsave(myplot, plot=pl2, width=7.22, height=7.22)
  
  pl3 = ggplot(pr, aes(x=real , y=fitted)) + geom_point(shape=1)
  myplot = paste0(outpath2, spp, ".fit_diff1_", mdln, ".png")
  ggsave(myplot, plot=pl3, width=7.22, height=7.22)
  
  pl4 = ggplot(pr, aes(x=index, y=diff)) + geom_point(shape=1)
  myplot = paste0(outpath2, spp, ".fit_diff2_", mdln, ".png")
  ggsave(myplot, plot=pl4, width=7.22, height=7.22)
  
  myplot = paste0(outpath2, spp, ".fit_resid_", mdln, ".png")
  png(myplot, units="px", width = 693.12, height = 693.12) # 7.22 inches = 693.12 px
  plot(resid)
  dev.off()
  
  # Save data
  mytable = paste0(outpath2, spp, ".fit_data_", mdln, ".csv")
  write.csv(pr, mytable, row.names = F)
  
  mytable = paste0(outpath2, spp, ".fit_param_", mdln, ".csv")
  write.csv(nls_med$parameters, mytable, row.names = T)
}


## (4) "double" power-law
# y = c + k*x^(-c*x^b)

# fit
fit.med = try( nlsLM(y.med ~  k * x.med^(-c * x.med^b), 
                     start = list(k=1, c=1, b = 1),
                     control = control1) )

fit.med = try( nlsLM(y.med ~  k * x.med^(-c * x.med^b), 
                     start = list(k=1, c=1, b = 0.01),
                     control = control1) )

fit.med = try( nlsLM(y.med ~  k * x.med^(-c * x.med^b), 
                     start = list(k=1, c=1, b = -0.1),
                     control = control1) )


if ( class(fit.med) != "try-error" ) {
  # Model parameters
  nls_med     = summary(fit.med)
  constant_k  = nls_med$parameters[1]
  constant_c  = nls_med$parameters[2]
  beta        = nls_med$parameters[3]
  
  # Add model to list
  idx = idx + 1
  fit.models[[idx]] = fit.med
  model.names[[idx]] = "M4"
  
  # Plot data + fit
  fitted_vals = fitted(fit.med)
  data.fit = as.data.frame(cbind(x.med, fitted_vals))
  
  pr4 = as.data.frame(cbind(y.med, round(fitted_vals, 3)))
  
  pr = pr4
  pr[, 3] = round(pr[,1]-pr[,2], 3)
  pr[, 4] = seq(2, (length(y.med)+1), by=1)
  colnames(pr) = c("real", "fitted", "diff", "index")
  
  # plot residuals
  resid = nlsResiduals(fit.med)
  
  # Save plots
  mdln = "model4"
  
  pl1 = ggplot(data.med, aes(x=N, y=new_genes)) + geom_point() + geom_line(data = data.fit, aes(x=x.med, y=fitted_vals), color="red") + ylim(0, max(max(data.fit$fitted_vals+0.5), max(y.med)))
  myplot = paste0(outpath2, spp, ".fit_line_", mdln, ".png")
  ggsave(myplot, plot=pl1, width=7.22, height=7.22)
  
  pl2 = ggplot(data.med, aes(x=N, y=new_genes)) + geom_point() + geom_point(data = data.fit, aes(x=x.med, y=fitted_vals), color="red", shape=1) + ylim(0, max(max(data.fit$fitted_vals+0.5), max(y.med)))
  myplot = paste0(outpath2, spp, ".fit_point_", mdln, ".png")
  ggsave(myplot, plot=pl2, width=7.22, height=7.22)
  
  pl3 = ggplot(pr, aes(x=real , y=fitted)) + geom_point(shape=1)
  myplot = paste0(outpath2, spp, ".fit_diff1_", mdln, ".png")
  ggsave(myplot, plot=pl3, width=7.22, height=7.22)
  
  pl4 = ggplot(pr, aes(x=index, y=diff)) + geom_point(shape=1)
  myplot = paste0(outpath2, spp, ".fit_diff2_", mdln, ".png")
  ggsave(myplot, plot=pl4, width=7.22, height=7.22)
  
  myplot = paste0(outpath2, spp, ".fit_resid_", mdln, ".png")
  png(myplot, units="px", width = 693.12, height = 693.12) # 7.22 inches = 693.12 px
  plot(resid)
  dev.off()
  
  # Save data
  mytable = paste0(outpath2, spp, ".fit_data_", mdln, ".csv")
  write.csv(pr, mytable, row.names = F)
  
  mytable = paste0(outpath2, spp, ".fit_param_", mdln, ".csv")
  write.csv(nls_med$parameters, mytable, row.names = T)
}


## (5) "simple" power-law
fit.med = try( nlsLM(y.med ~ k*x.med^(-a), 
                     start = list(a = 1, k = 1),
                     control = control1) )

if ( class(fit.med) != "try-error" ) {
  # Model parameters
  nls_med   = summary(fit.med)
  alpha     = nls_med$parameters[1]
  constant  = nls_med$parameters[2]
  
  # Add model to list
  idx = idx + 1
  fit.models[[idx]] = fit.med
  model.names[[idx]] = "M5"
  
  # Plot data + fit
  fitted_vals = fitted(fit.med)
  data.fit = as.data.frame(cbind(x.med, fitted_vals))
  
  pr5 = as.data.frame(cbind(y.med, round(fitted_vals, 3)))
  
  pr = pr5
  pr[, 3] = round(pr[,1]-pr[,2], 3)
  pr[, 4] = seq(2, (length(y.med)+1), by=1)
  colnames(pr) = c("real", "fitted", "diff", "index")
  
  # plot residuals
  resid = nlsResiduals(fit.med)
  
  # Save plots
  mdln = "model5"
  
  pl1 = ggplot(data.med, aes(x=N, y=new_genes)) + geom_point() + geom_line(data = data.fit, aes(x=x.med, y=fitted_vals), color="red") + ylim(0, max(max(data.fit$fitted_vals+0.5), max(y.med)))
  myplot = paste0(outpath2, spp, ".fit_line_", mdln, ".png")
  ggsave(myplot, plot=pl1, width=7.22, height=7.22)
  
  pl2 = ggplot(data.med, aes(x=N, y=new_genes)) + geom_point() + geom_point(data = data.fit, aes(x=x.med, y=fitted_vals), color="red", shape=1) + ylim(0, max(max(data.fit$fitted_vals+0.5), max(y.med)))
  myplot = paste0(outpath2, spp, ".fit_point_", mdln, ".png")
  ggsave(myplot, plot=pl2, width=7.22, height=7.22)
  
  pl3 = ggplot(pr, aes(x=real , y=fitted)) + geom_point(shape=1)
  myplot = paste0(outpath2, spp, ".fit_diff1_", mdln, ".png")
  ggsave(myplot, plot=pl3, width=7.22, height=7.22)
  
  pl4 = ggplot(pr, aes(x=index, y=diff)) + geom_point(shape=1)
  myplot = paste0(outpath2, spp, ".fit_diff2_", mdln, ".png")
  ggsave(myplot, plot=pl4, width=7.22, height=7.22)
  
  myplot = paste0(outpath2, spp, ".fit_resid_", mdln, ".png")
  png(myplot, units="px", width = 693.12, height = 693.12) # 7.22 inches = 693.12 px
  plot(resid)
  dev.off()
  
  # Save data
  mytable = paste0(outpath2, spp, ".fit_data_", mdln, ".csv")
  write.csv(pr, mytable, row.names = F)
  
  mytable = paste0(outpath2, spp, ".fit_param_", mdln, ".csv")
  write.csv(nls_med$parameters, mytable, row.names = T)
}


## (6) Exponential decay
fit.med = try( nlsLM(y.med ~ 0 + k*exp(-x.med*tau), 
                     start   = list(k=1, tau=1),
                     control = control1) ) 

if ( class(fit.med) != "try-error" ) {
  # Model parameters
  nls_med = summary(fit.med)
  tan.omega = 0 # assuming tg(omega)=0
  constant  = nls_med$parameters[1]
  tau.model = nls_med$parameters[2]
  tau.orig  = 1/nls_med$parameters[2]
  
  # Add model to list
  idx = idx + 1
  fit.models[[idx]] = fit.med
  model.names[[idx]] = "M6"
  
  # Plot data + fit
  fitted_vals = fitted(fit.med)
  data.fit = as.data.frame(cbind(x.med, fitted_vals))
  
  pr6 = as.data.frame(cbind(y.med, round(fitted_vals, 1)))
  
  pr = pr6
  pr[, 3] = round(pr[,1]-pr[,2], 3)
  pr[, 4] = seq(2, (length(y.med)+1), by=1)
  colnames(pr) = c("real", "fitted", "diff", "index")
  ggplot(pr, aes(x=real , y=fitted)) + geom_point(shape=1)
  ggplot(pr, aes(x=index, y=diff)) + geom_point(shape=1)
  
  # plot residuals
  resid = nlsResiduals(fit.med)
  
  # Save plots
  mdln = "model6"
  
  pl1 = ggplot(data.med, aes(x=N, y=new_genes)) + geom_point() + geom_line(data = data.fit, aes(x=x.med, y=fitted_vals), color="red") + ylim(0, max(max(data.fit$fitted_vals+0.5), max(y.med)))
  myplot = paste0(outpath2, spp, ".fit_line_", mdln, ".png")
  ggsave(myplot, plot=pl1, width=7.22, height=7.22)
  
  pl2 = ggplot(data.med, aes(x=N, y=new_genes)) + geom_point() + geom_point(data = data.fit, aes(x=x.med, y=fitted_vals), color="red", shape=1) + ylim(0, max(max(data.fit$fitted_vals+0.5), max(y.med)))
  myplot = paste0(outpath2, spp, ".fit_point_", mdln, ".png")
  ggsave(myplot, plot=pl2, width=7.22, height=7.22)
  
  pl3 = ggplot(pr, aes(x=real , y=fitted)) + geom_point(shape=1)
  myplot = paste0(outpath2, spp, ".fit_diff1_", mdln, ".png")
  ggsave(myplot, plot=pl3, width=7.22, height=7.22)
  
  pl4 = ggplot(pr, aes(x=index, y=diff)) + geom_point(shape=1)
  myplot = paste0(outpath2, spp, ".fit_diff2_", mdln, ".png")
  ggsave(myplot, plot=pl4, width=7.22, height=7.22)
  
  myplot = paste0(outpath2, spp, ".fit_resid_", mdln, ".png")
  png(myplot, units="px", width = 693.12, height = 693.12) # 7.22 inches = 693.12 px
  plot(resid)
  dev.off()
  
  # Save data
  mytable = paste0(outpath2, spp, ".fit_data_", mdln, ".csv")
  write.csv(pr, mytable, row.names = F)
  
  mytable = paste0(outpath2, spp, ".fit_param_", mdln, ".csv")
  write.csv(nls_med$parameters, mytable, row.names = T)
}

## Model selection (AIC)

# Test
aic.res = aictab(cand.set = fit.models, modnames = unlist(model.names))
aic.res

# Save test results
mytable = paste0(outpath2, spp, ".AIC_exp.csv")
write.csv(aic.res, mytable, row.names = F)
