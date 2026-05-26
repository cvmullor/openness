# Estimate power law's alpha (openness) as a function of N using exponential model fit
#
# Written by Carlos Valiente-Mullor
#
# Usage:
#       Rscript calculate_alphaN_expm.R [sp_name]
#
# Input file name 'gene_presence_absence_roary.{species_name}.csv' is expected by default
#
# Example:
#       Rscript calculate_alphaN_expm.R Enterococcus_faecium


## packages
library(pagoo)
library(micropan)
library(matrixStats)
library(stringr)
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

## Arguments
args = commandArgs(trailingOnly=TRUE)
spp = args[1] # species name ("Genus_epithet")

## Paths
inpath = "path/to/input/dir"   # set input dir
outpath= "path/to/output/dir"  # set output dir

## Read input / pan-matrix
myfile_r = paste0(inpath, "gene_presence_absence_roary.", spp, ".csv") # Roary-formatted PAM from Panaroo
pg = roary_2_pagoo(gene_presence_absence_csv = myfile_r) # Panaroo table to pagoo object


#############################################
#   Estimate alpha(N) - Exponential model   #
#############################################

## pan.matrix
panmat1 = pg$pan_matrix # extract pan-matrix
pan.matrix = panmat1

n.perm = 100 # number of permutations to obtain pangneome curve in each iteration
step_N = 1   # num. of genomes added in each iteration

N = dim(panmat1)[1]
N_values = c(seq(from=3, to=N, by=step_N), N) # Default initial N=3

res.def.all = c()
res.def.mean = c()
res.def.med = c()

res.def.all.in = c()  # *.in: results including zeros
res.def.mean.in = c()
res.def.med.in = c()

data_curve.all = list()
data_curve.mean = list()
data_curve.med = list()

data_curve.all.in = list()
data_curve.mean.in = list()
data_curve.med.in = list()

nmat.or <- powerlaw_perm(pan.matrix, n.perm)

# Control parametrers
control1 <- nls.control(maxiter= 1000, minFactor= 1e-30, warnOnly= FALSE,tol=1e-05)

# Add the next genome (rarefaction) each iteration
for (i in 1:length(N_values)) {
  
  # Random sub-sample of permutations matrix
  if (i < length(N_values)) {
    nmat = nmat.or[seq(N_values[i]-1),] # N-1 because we ignore 1st value
    
  } else if (N_values[i]-1 == N-1) {
    nmat = nmat.or
  }
  
  ####### Estimates excluding values==0 #######
  
  #### All permutations ####
  # Prepare data
  x.all <- rep((2:N_values[i]), times = n.perm)
  y.all <- as.numeric(nmat)
  y.all[y.all == 0] <- NA
  
  # add to object/matrix to plot An curves of every run (N value)
  data_curve.all[[i]] = cbind(x.all, y.all)
  
  # Fit power-law model with 'nlsLM' and extract data if no error
  fit.med = try( nlsLM(y.all ~ k*x.all^(-a), 
                       start = list(a = 1, k = 1),
                       control = control1,
                       na.action = na.exclude) )
  
  if ( class(fit.med) != "try-error" ) {
    # Model parameters
    fit_res       = summary(fit.med)
    fit_res.alpha = fit_res$parameters[1]
    fit_res.se    = fit_res$parameters[3]
    
    res = c("[all,ex] nlsLM", i, N_values[i], fit_res.alpha, fit_res.se)
    names(res) = c("method", "index", "N", "alpha", "SE")
    
    res.def.all = rbind(res.def.all, res)
  }
  
  #### Mean of permutations ####
  
  # Prepare data
  nmat_mean = rowMeans(nmat)
  x.mean <- rep((2:N_values[i]), times = 1)
  y.mean <- as.numeric(nmat_mean)
  y.mean[y.mean == 0] <- NA
  
  # add to object/matrix to plot An curves of every run (N value)
  data_curve.mean[[i]] = cbind(x.mean, y.mean)
  
  # Fit power-law model with 'nlsLM' and extract data if no error
  fit.med = try( nlsLM(y.mean ~ k*x.mean^(-a), 
                       start = list(a = 1, k = 1),
                       control = control1,
                       na.action = na.exclude) )
  
  if ( class(fit.med) != "try-error" ) {
    # Model parameters
    fit_res       = summary(fit.med)
    fit_res.alpha = fit_res$parameters[1]
    fit_res.se    = fit_res$parameters[3]
    
    res = c("[mean,ex] nlsLM", i, N_values[i], fit_res.alpha, fit_res.se)
    names(res) = c("method", "index", "N", "alpha", "SE")
    
    res.def.mean = rbind(res.def.mean, res)
  }
  
  #### Median of permutations ####
  
  # Prepare data
  nmat_med = rowMedians(nmat)
  x.med <- rep(2:N_values[i], times = 1)
  y.med <- as.numeric(nmat_med)
  y.med[y.med == 0] <- NA
  
  # add to object/matrix to plot An curves of every run (N value)
  data_curve.med[[i]] = cbind(x.med, y.med)
  
  # Fit power-law model with 'nlsLM' and extract data if no error
  fit.med = try( nlsLM(y.med ~ k*x.med^(-a), 
                       start = list(a = 1, k = 1),
                       control = control1,
                       na.action = na.exclude) )
  
  if ( class(fit.med) != "try-error" ) {
    # Model parameters
    fit_res       = summary(fit.med)
    fit_res.alpha = fit_res$parameters[1]
    fit_res.se    = fit_res$parameters[3]
    
    res = c("[med,ex] nlsLM", i, N_values[i], fit_res.alpha, fit_res.se)
    names(res) = c("method", "index", "N", "alpha", "SE")
    
    res.def.med = rbind(res.def.med, res)
  }
  
  
  ####### Estimates excluding values==0 #######
  
  #### All permutations ####
  
  # Prepare data
  x.all <- rep((2:N_values[i]), times = n.perm)
  y.all <- as.numeric(nmat)
  
  # Add to object/matrix to plot An curves of every run (N value) 
  data_curve.all.in[[i]] = cbind(x.all, y.all)
  
  # Fit power-law model with 'nlsLM' and extract data if no error
  fit.med = try( nlsLM(y.all ~ k*x.all^(-a), 
                       start = list(a = 1, k = 1),
                       control = control1,
                       na.action = na.exclude) )
  
  if ( class(fit.med) != "try-error" ) {
    # Model parameters
    fit_res       = summary(fit.med)
    fit_res.alpha = fit_res$parameters[1]
    fit_res.se    = fit_res$parameters[3]
    
    res = c("[all,in] nlsLM", i, N_values[i], fit_res.alpha, fit_res.se)
    names(res) = c("method", "index", "N", "alpha", "SE")
    
    res.def.all.in = rbind(res.def.all.in, res)
  }
  
  #### Mean of permutations ####
  
  # Prepare data
  nmat_mean = rowMeans(nmat)
  x.mean <- rep((2:N_values[i]), times = 1)
  y.mean <- as.numeric(nmat_mean)
  
  # add to object/matrix to plot An curves of every run (N value)
  data_curve.mean.in[[i]] = cbind(x.mean, y.mean)
  
  # Fit power-law model with 'nlsLM' and extract data if no error
  fit.med = try( nlsLM(y.mean ~ k*x.mean^(-a), 
                       start = list(a = 1, k = 1),
                       control = control1,
                       na.action = na.exclude) )
  
  if ( class(fit.med) != "try-error" ) {
    # Model parameters
    fit_res       = summary(fit.med)
    fit_res.alpha = fit_res$parameters[1]
    fit_res.se    = fit_res$parameters[3]
    
    res = c("[mean,in] nlsLM", i, N_values[i], fit_res.alpha, fit_res.se)
    names(res) = c("method", "index", "N", "alpha", "SE")
    
    res.def.mean.in = rbind(res.def.mean.in, res)
  }
  
  #### Median of permutations ####
  
  # Prepare data
  nmat_med = rowMedians(nmat)
  x.med <- rep(2:N_values[i], times = 1)
  y.med <- as.numeric(nmat_med)
  
  # add to object/matrix to plot An curves of every run (N value)
  data_curve.med.in[[i]] = cbind(x.med, y.med)
  
  # Fit power-law model with 'nlsLM' and extract data if no error
  fit.med = try( nlsLM(y.med ~ k*x.med^(-a), 
                       start = list(a = 1, k = 1),
                       control = control1,
                       na.action = na.exclude) )
  
  if ( class(fit.med) != "try-error" ) {
    # Model parameters
    fit_res       = summary(fit.med)
    fit_res.alpha = fit_res$parameters[1]
    fit_res.se    = fit_res$parameters[3]
    
    res = c("[med,in] nlsLM", i, N_values[i], fit_res.alpha, fit_res.se)
    names(res) = c("method", "index", "N", "alpha", "SE")
    
    res.def.med.in = rbind(res.def.med.in, res)
  }
  
}

# Datasets
# Exclude zeros
a_perN.all = as.data.frame(res.def.all)
a_perN.all[,2] = as.numeric(a_perN.all[,2])
a_perN.all[,3] = as.numeric(a_perN.all[,3])
a_perN.all[,4] = as.numeric(a_perN.all[,4])
a_perN.all[,5] = as.numeric(a_perN.all[,5])

a_perN.mean = as.data.frame(res.def.mean)
a_perN.mean[,2] = as.numeric(a_perN.mean[,2])
a_perN.mean[,3] = as.numeric(a_perN.mean[,3])
a_perN.mean[,4] = as.numeric(a_perN.mean[,4])
a_perN.mean[,5] = as.numeric(a_perN.mean[,5])

a_perN.med = as.data.frame(res.def.med)
a_perN.med[,2] = as.numeric(a_perN.med[,2])
a_perN.med[,3] = as.numeric(a_perN.med[,3])
a_perN.med[,4] = as.numeric(a_perN.med[,4])
a_perN.med[,5] = as.numeric(a_perN.med[,5])

# Include zeros
a_perN.all.in = as.data.frame(res.def.all.in)
a_perN.all.in[,2] = as.numeric(a_perN.all.in[,2])
a_perN.all.in[,3] = as.numeric(a_perN.all.in[,3])
a_perN.all.in[,4] = as.numeric(a_perN.all.in[,4])
a_perN.all.in[,5] = as.numeric(a_perN.all.in[,5])

a_perN.mean.in = as.data.frame(res.def.mean.in)
a_perN.mean.in[,2] = as.numeric(a_perN.mean.in[,2])
a_perN.mean.in[,3] = as.numeric(a_perN.mean.in[,3])
a_perN.mean.in[,4] = as.numeric(a_perN.mean.in[,4])
a_perN.mean.in[,5] = as.numeric(a_perN.mean.in[,5])

a_perN.med.in = as.data.frame(res.def.med.in)
a_perN.med.in[,2] = as.numeric(a_perN.med.in[,2])
a_perN.med.in[,3] = as.numeric(a_perN.med.in[,3])
a_perN.med.in[,4] = as.numeric(a_perN.med.in[,4])
a_perN.med.in[,5] = as.numeric(a_perN.med.in[,5])


# Write estimate results:
# Exclude zeros
rownames(a_perN.all) = NULL
rownames(a_perN.mean) = NULL
rownames(a_perN.med) = NULL

mytable = paste0(outpath, spp, ".perm_all.nlsLM_ex.alpha_perN.csv")
write.csv(a_perN.all, mytable, row.names = F)

mytable = paste0(outpath, spp, ".perm_mean.nlsLM_ex.alpha_perN.csv")
write.csv(a_perN.mean, mytable, row.names = F)

mytable = paste0(outpath, spp, ".perm_med.nlsLM_ex.alpha_perN.csv")
write.csv(a_perN.med, mytable, row.names = F)

# Include zeros
rownames(a_perN.all.in) = NULL
rownames(a_perN.mean.in) = NULL
rownames(a_perN.med.in) = NULL

mytable = paste0(outpath, spp, ".perm_all.nlsLM_in.alpha_perN.csv")
write.csv(a_perN.all.in, mytable, row.names = F)

mytable = paste0(outpath, spp, ".perm_mean.nlsLM_in.alpha_perN.csv")
write.csv(a_perN.mean.in, mytable, row.names = F)

mytable = paste0(outpath, spp, ".perm_med.nlsLM_in.alpha_perN.csv")
write.csv(a_perN.med.in, mytable, row.names = F)


# Data frames
# Exclude zeros
df.3methods = cbind(a_perN.all$N, a_perN.all$alpha, a_perN.mean$alpha, a_perN.med$alpha)
colnames(df.3methods) = c("N", "alpha.all", "alpha.mean", "alpha.median")

mytable = paste0(outpath, spp, ".perm_summ.nlsLM_ex.alpha_perN.csv")
write.csv(df.3methods, mytable, row.names = F)

# Include zeros
df.3methods.in = cbind(a_perN.all.in$N, a_perN.all.in$alpha, a_perN.mean.in$alpha, a_perN.med.in$alpha)
colnames(df.3methods.in) = c("N", "alpha.all.in", "alpha.mean.in", "alpha.median.in")

mytable = paste0(outpath, spp, ".perm_summ.nlsLM_in.alpha_perN.csv")
write.csv(df.3methods.in, mytable, row.names = F)


#############
##  PLOTS  ##
#############

# Save alpha curve plots
# Exclude zeros
curve.all = ggplot(df.3methods, aes(x=N, y=alpha.all)) + geom_point() + geom_smooth(se=F)
curve.mean = ggplot(df.3methods, aes(x=N, y=alpha.mean)) + geom_point() + geom_smooth(se=F)
curve.med = ggplot(df.3methods, aes(x=N, y=alpha.median)) + geom_point() + geom_smooth(se=F)

myplot.all = paste0(outpath, spp, ".curve_all.nlsLM_ex.alpha_perN.png")
myplot.mean = paste0(outpath, spp, ".curve_mean.nlsLM_ex.alpha_perN.png")
myplot.med = paste0(outpath, spp, ".curve_med.nlsLM_ex.alpha_perN.png")

ggsave(myplot.all, plot=curve.all)
ggsave(myplot.mean, plot=curve.mean)
ggsave(myplot.med, plot=curve.med)

# Include zeros
curve.all.in  = ggplot(df.3methods.in, aes(x=N, y=alpha.all.in)) + geom_point() + geom_smooth(se=F)
curve.mean.in = ggplot(df.3methods.in, aes(x=N, y=alpha.mean.in)) + geom_point() + geom_smooth(se=F)
curve.med.in  = ggplot(df.3methods.in, aes(x=N, y=alpha.median.in)) + geom_point() + geom_smooth(se=F)

myplot.all.in  = paste0(outpath, spp, ".curve_all.nlsLM_in.alpha_perN.png")
myplot.mean.in = paste0(outpath, spp, ".curve_mean.nlsLM_in.alpha_perN.png")
myplot.med.in  = paste0(outpath, spp, ".curve_med.nlsLm_in.alpha_perN.png")

ggsave(myplot.all.in, plot=curve.all.in)
ggsave(myplot.mean.in, plot=curve.mean.in)
ggsave(myplot.med.in, plot=curve.med.in)


# Data frame data curve
# Exclude zeros
data_curve.all = lapply(data_curve.all, as.data.frame)
data_curve.mean = lapply(data_curve.mean, as.data.frame)
data_curve.med = lapply(data_curve.med, as.data.frame)

# Include zeros
data_curve.all.in  = lapply(data_curve.all.in, as.data.frame)
data_curve.mean.in = lapply(data_curve.mean.in, as.data.frame)
data_curve.med.in  = lapply(data_curve.med.in, as.data.frame)

# Log-log curve data
# Exclude zeros
log_curve.all = lapply(data_curve.all, log)
log_curve.mean = lapply(data_curve.mean, log)
log_curve.med = lapply(data_curve.med, log)

# Include zeros
log_curve.all.in = lapply(data_curve.all.in, log)
log_curve.mean.in = lapply(data_curve.mean.in, log)
log_curve.med.in = lapply(data_curve.med.in, log)

# Plot rarefaction curve for each total N value:
# Exclude zeros
rar.all.outpath = paste0(outpath, spp, ".alpha.raref_all.nlsLM_ex.alpha_perN.pdf")
rar.mean.outpath = paste0(outpath, spp, ".alpha.raref_mean.nlsLM_ex.alpha_perN.pdf")
rar.med.outpath = paste0(outpath, spp, ".alpha.raref_med.nlsLM_ex.alpha_perN.pdf")

# Include zeros
rar.all.outpath.in  = paste0(outpath, spp, ".alpha.raref_all.nlsLM_in.alpha_perN.pdf")
rar.mean.outpath.in = paste0(outpath, spp, ".alpha.raref_mean.nlsLM_in.alpha_perN.pdf")
rar.med.outpath.in  = paste0(outpath, spp, ".alpha.raref_med.nlsLM_in.alpha_perN.pdf")


## ALL perms

# Get x, y lims
# ... max and min N, An +/- 0.5
last.v = length(data_curve.all)
x_lim.max = max(data_curve.all[[last.v]]$x.all, na.rm = T)+0.5
x_lim.min = min(data_curve.all[[last.v]]$x.all, na.rm = T)-0.5

y_lim.max = max(data_curve.all[[last.v]]$y.all, na.rm = T)+0.5
y_lim.min = min(data_curve.all[[last.v]]$y.all, na.rm = T)-0.5

# plots
x.text_limit = x_lim.max-(x_lim.max*0.13)
y.text_limit = y_lim.max-1

list4 = list()
for (i in 1:length(data_curve.all)){
  lbl = paste0("\u03b1 = ", round(a_perN.all$alpha[i],3))
  list4[[i]] <- ggplot(data_curve.all[[i]], aes(x=x.all, y=y.all)) + 
    geom_point() + 
    #geom_smooth(method="loess", se=F) +
    xlim(x_lim.min,x_lim.max) + ylim(y_lim.min,y_lim.max) +
    annotate("text", x=x.text_limit, y=y.text_limit, 
             label = lbl, 
             size=7, color="black")
}
# (save w/ cairo_pdf in order to print "alpha" character)
grDevices::cairo_pdf(filename = rar.all.outpath, onefile = T)
list4
dev.off()


## MEAN of perms

# Get x, y lims
# ... max and min N, An +/- 0.5
last.v = length(data_curve.mean)
x_lim.max = max(data_curve.mean[[last.v]]$x.mean, na.rm = T)+0.5
x_lim.min = min(data_curve.mean[[last.v]]$x.mean, na.rm = T)-0.5

y_lim.max = max(data_curve.mean[[last.v]]$y.mean, na.rm = T)+0.5
y_lim.min = min(data_curve.mean[[last.v]]$y.mean, na.rm = T)-0.5

# plots
x.text_limit = x_lim.max-(x_lim.max*0.13)
y.text_limit = y_lim.max-1

list4 = list()
for (i in 1:length(data_curve.mean)){
  lbl = paste0("\u03b1 = ", round(a_perN.mean$alpha[i],3))
  list4[[i]] <- ggplot(data_curve.mean[[i]], aes(x=x.mean, y=y.mean)) + 
    geom_point() + 
    #geom_smooth(se=F, formula = y ~ x^(a_perN.mean$alpha[i])) +
    xlim(x_lim.min,x_lim.max) + ylim(y_lim.min,y_lim.max) +
    annotate("text", x=x.text_limit, y=y.text_limit, 
             label = lbl, 
             size=7, color="black")
}
# (save w/ cairo_pdf in order to print "alpha" character)
grDevices::cairo_pdf(filename = rar.mean.outpath, onefile = T)
list4
dev.off()


## MEDIANS of perms

# Get x, y lims
# ... max and min log(N), log(An) +/- 0.5
last.v = length(data_curve.med)
x_lim.max = max(data_curve.med[[last.v]]$x.med, na.rm = T)+0.5
x_lim.min = min(data_curve.med[[last.v]]$x.med, na.rm = T)-0.5

y_lim.max = max(data_curve.med[[last.v]]$y.med, na.rm = T)+0.5
y_lim.min = min(data_curve.med[[last.v]]$y.med, na.rm = T)-0.5

# plots
x.text_limit = x_lim.max-(x_lim.max*0.13)
y.text_limit = y_lim.max-1

list4 = list()
for (i in 1:length(data_curve.med)){
  lbl = paste0("\u03b1 = ", round(a_perN.med$alpha[i],3))
  list4[[i]] <- ggplot(data_curve.med[[i]], aes(x=x.med, y=y.med)) + 
    geom_point() + 
    #geom_smooth(method="loess", se=F) +
    xlim(x_lim.min,x_lim.max) + ylim(y_lim.min,y_lim.max) +
    annotate("text", x=x.text_limit, y=y.text_limit, 
             label = lbl, 
             size=7, color="black")
}
# (save w/ cairo_pdf in order to print "alpha" character)
grDevices::cairo_pdf(filename = rar.med.outpath, onefile = T)
list4
dev.off()


## ALL perms (in)

# Get x, y lims
# ... max and min N, An +/- 0.5
last.v = length(data_curve.all.in)
x_lim.max = max(data_curve.all.in[[last.v]]$x.all, na.rm = T)+0.5
x_lim.min = min(data_curve.all.in[[last.v]]$x.all, na.rm = T)-0.5

y_lim.max = max(data_curve.all.in[[last.v]]$y.all, na.rm = T)+0.5
y_lim.min = min(data_curve.all.in[[last.v]]$y.all, na.rm = T)-0.5

# plots
x.text_limit = x_lim.max-(x_lim.max*0.13)
y.text_limit = y_lim.max-1

list4 = list()
for (i in 1:length(data_curve.all.in)){
  lbl = paste0("\u03b1 = ", round(a_perN.all.in$alpha[i],3))
  list4[[i]] <- ggplot(data_curve.all.in[[i]], aes(x=x.all, y=y.all)) + 
    geom_point() + 
    #geom_smooth(method="loess", se=F) +
    xlim(x_lim.min,x_lim.max) + ylim(y_lim.min,y_lim.max) +
    annotate("text", x=x.text_limit, y=y.text_limit, 
             label = lbl, 
             size=7, color="black")
}
# (save w/ cairo_pdf in order to print "alpha" character)
grDevices::cairo_pdf(filename = rar.all.outpath.in, onefile = T)
list4
dev.off()


## MEAN of perms (in)

# Get x, y lims
# ... max and min N, An +/- 0.5
last.v = length(data_curve.mean.in)
x_lim.max = max(data_curve.mean.in[[last.v]]$x.mean, na.rm = T)+0.5
x_lim.min = min(data_curve.mean.in[[last.v]]$x.mean, na.rm = T)-0.5

y_lim.max = max(data_curve.mean.in[[last.v]]$y.mean, na.rm = T)+0.5
y_lim.min = min(data_curve.mean.in[[last.v]]$y.mean, na.rm = T)-0.5

# plots
x.text_limit = x_lim.max-(x_lim.max*0.13)
y.text_limit = y_lim.max-1

list4 = list()
for (i in 1:length(data_curve.mean.in)){
  lbl = paste0("\u03b1 = ", round(a_perN.mean.in$alpha[i],3))
  list4[[i]] <- ggplot(data_curve.mean.in[[i]], aes(x=x.mean, y=y.mean)) + 
    geom_point() + 
    #geom_smooth(se=F, formula = y ~ x^(a_perN.mean$alpha[i])) +
    xlim(x_lim.min,x_lim.max) + ylim(y_lim.min,y_lim.max) +
    annotate("text", x=x.text_limit, y=y.text_limit, 
             label = lbl, 
             size=7, color="black")
}
# (save w/ cairo_pdf in order to print "alpha" character)
grDevices::cairo_pdf(filename = rar.mean.outpath.in, onefile = T)
list4
dev.off()


## MEDIANS of perms (in)

# Get x, y lims
# ... max and min log(N), log(An) +/- 0.5
last.v = length(data_curve.med.in)
x_lim.max = max(data_curve.med.in[[last.v]]$x.med, na.rm = T)+0.5
x_lim.min = min(data_curve.med.in[[last.v]]$x.med, na.rm = T)-0.5

y_lim.max = max(data_curve.med.in[[last.v]]$y.med, na.rm = T)+0.5
y_lim.min = min(data_curve.med.in[[last.v]]$y.med, na.rm = T)-0.5

# plots
x.text_limit = x_lim.max-(x_lim.max*0.13)
y.text_limit = y_lim.max-1

list4 = list()
for (i in 1:length(data_curve.med.in)){
  lbl = paste0("\u03b1 = ", round(a_perN.med.in$alpha[i],3))
  list4[[i]] <- ggplot(data_curve.med.in[[i]], aes(x=x.med, y=y.med)) + 
    geom_point() + 
    #geom_smooth(method="loess", se=F) +
    xlim(x_lim.min,x_lim.max) + ylim(y_lim.min,y_lim.max) +
    annotate("text", x=x.text_limit, y=y.text_limit, 
             label = lbl, 
             size=7, color="black")
}
# (save w/ cairo_pdf in order to print "alpha" character)
grDevices::cairo_pdf(filename = rar.med.outpath.in, onefile = T)
list4
dev.off()


# Plot log-log rarefaction curve for each total N value:
# outpaths:
log.all.outpath = paste0(outpath, spp, ".alpha.log-log_all.nlsLM_ex.alpha_perN.pdf")
log.mean.outpath = paste0(outpath, spp, ".alpha.log-log_mean.nlsLM_ex.alpha_perN.pdf")
log.med.outpath = paste0(outpath, spp, ".alpha.log-log_med.nlsLM_ex.alpha_perN.pdf")

log.all.outpath.in = paste0(outpath, spp, ".alpha.log-log_all.nlsLM_in.alpha_perN.pdf")
log.mean.outpath.in = paste0(outpath, spp, ".alpha.log-log_mean.nlsLM_in.alpha_perN.pdf")
log.med.outpath.in = paste0(outpath, spp, ".alpha.log-log_med.nlsLM_in.alpha_perN.pdf")


### (Ex)

## ALL perms

# Get x, y lims
# ... max and min log(N), log(An) +/- 0.5
last.v = length(log_curve.all)
x_lim.max = max(log_curve.all[[last.v]]$x.all, na.rm = T)+0.5
x_lim.min = min(log_curve.all[[last.v]]$x.all, na.rm = T)-0.5

y_lim.max = max(log_curve.all[[last.v]]$y.all, na.rm = T)+0.5
y_lim.min = min(log_curve.all[[last.v]]$y.all, na.rm = T)-0.5

# plots
x.text_limit = x_lim.max-1
y.text_limit = y_lim.max-1

list4 = list()
for (i in 1:length(log_curve.all)){
  lbl = paste0("\u03b1 = ", round(a_perN.all$alpha[i],3))
  list4[[i]] <- ggplot(log_curve.all[[i]], aes(x=x.all, y=y.all)) + 
    geom_point() + 
    geom_smooth(method="lm", se=F) +
    xlim(x_lim.min,x_lim.max) + ylim(y_lim.min,y_lim.max) +
    annotate("text", x=x.text_limit, y=y.text_limit, 
             label = lbl, 
             size=7, color="black")
}
# (save w/ cairo_pdf in order to print "alpha" character)
grDevices::cairo_pdf(filename = log.all.outpath, onefile = T)
list4
dev.off()


## MEAN of perms

# Get x, y lims
# ... max and min log(N), log(An) +/- 0.5
last.v = length(log_curve.mean)
x_lim.max = max(log_curve.mean[[last.v]]$x.mean, na.rm = T)+0.5
x_lim.min = min(log_curve.mean[[last.v]]$x.mean, na.rm = T)-0.5

y_lim.max = max(log_curve.mean[[last.v]]$y.mean, na.rm = T)+0.5
y_lim.min = min(log_curve.mean[[last.v]]$y.mean, na.rm = T)-0.5

# plots
x.text_limit = x_lim.max-1
y.text_limit = y_lim.max-1

list4 = list()
for (i in 1:length(log_curve.mean)){
  lbl = paste0("\u03b1 = ", round(a_perN.mean$alpha[i],3))
  list4[[i]] <- ggplot(log_curve.mean[[i]], aes(x=x.mean, y=y.mean)) + 
    geom_point() + 
    geom_smooth(method="lm", se=F) +
    xlim(x_lim.min,x_lim.max) + ylim(y_lim.min,y_lim.max) +
    annotate("text", x=x.text_limit, y=y.text_limit, 
             label = lbl, 
             size=7, color="black")
}
# (save w/ cairo_pdf in order to print "alpha" character)
grDevices::cairo_pdf(filename = log.mean.outpath, onefile = T)
list4
dev.off()


## MEDIANS of perms

# Get x, y lims
# ... max and min log(N), log(An) +/- 0.5
last.v = length(log_curve.med)
x_lim.max = max(log_curve.med[[last.v]]$x.med, na.rm = T)+0.5
x_lim.min = min(log_curve.med[[last.v]]$x.med, na.rm = T)-0.5

y_lim.max = max(log_curve.med[[last.v]]$y.med, na.rm = T)+0.5
y_lim.min = min(log_curve.med[[last.v]]$y.med, na.rm = T)-0.5

# plots
x.text_limit = x_lim.max-1
y.text_limit = y_lim.max-1

list4 = list()
for (i in 1:length(log_curve.med)){
  lbl = paste0("\u03b1 = ", round(a_perN.med$alpha[i],3))
  list4[[i]] <- ggplot(log_curve.med[[i]], aes(x=x.med, y=y.med)) + 
    geom_point() + 
    geom_smooth(method="lm", se=F) +
    xlim(x_lim.min,x_lim.max) + ylim(y_lim.min,y_lim.max) +
    annotate("text", x=x.text_limit, y=y.text_limit, 
             label = lbl, 
             size=7, color="black")
}
# (save w/ cairo_pdf in order to print "alpha" character)
grDevices::cairo_pdf(filename = log.med.outpath, onefile = T)
list4
dev.off()


### (In) ; solo hacemos las medias, ya que ln(0)=-Inf (solo las medias, a veces, no tienen 0s)

## MEAN of perms

# Get x, y lims
# ... max and min log(N), log(An) +/- 0.5
last.v = length(log_curve.mean.in)
x_lim.max = max(log_curve.mean.in[[last.v]]$x.mean, na.rm = T)+0.5
x_lim.min = min(log_curve.mean.in[[last.v]]$x.mean, na.rm = T)-0.5

y_lim.max = max(log_curve.mean.in[[last.v]]$y.mean, na.rm = T)+0.5
y_lim.min = min(log_curve.mean.in[[last.v]]$y.mean, na.rm = T)-0.5

# plots
x.text_limit = x_lim.max-1
y.text_limit = y_lim.max-1

list4 = list()
for (i in 1:length(log_curve.mean.in)){
  lbl = paste0("\u03b1 = ", round(a_perN.mean.in$alpha[i],3))
  list4[[i]] <- ggplot(log_curve.mean.in[[i]], aes(x=x.mean, y=y.mean)) + 
    geom_point() + 
    geom_smooth(method="lm", se=F) +
    xlim(x_lim.min,x_lim.max) + ylim(y_lim.min,y_lim.max) +
    annotate("text", x=x.text_limit, y=y.text_limit, 
             label = lbl, 
             size=7, color="black")
}
# (save w/ cairo_pdf in order to print "alpha" character)
grDevices::cairo_pdf(filename = log.mean.outpath.in, onefile = T)
list4
dev.off()


## Save info/summary matrix of permutations (total N)
nmat = nmat.or

N = dim(panmat1)[1]
N_values = seq(2, N, by=1)

rownames(nmat) = N_values # N (nº genomas incluidos)
colnames(nmat) = paste0("r", seq(1, 100))

#rms = rowMeans(nmat)
#rmd = rowMedians(nmat)

summ.N = t(apply(nmat, MARGIN=1, summary)) # summary de cada N a lo largo de la curva

# write tables (rarefaction matrix + summary)
summ_file = paste0(outpath, spp, ".raref_summary.csv")
raref_file = paste0(outpath, spp, ".raref_table.csv")

write.csv(summ.N, summ_file)
write.csv(nmat, raref_file)


# Plot boxplots for each N
x.all = rep(N_values, times = n.perm)
y.all = as.numeric(nmat)

df.all = as.data.frame(cbind(x.all, y.all))

pl = ggplot(df.all, aes(x=x.all, y=y.all, group=x.all)) + 
  geom_boxplot(outlier.size = 0.1) +
  xlab("Genomes") + 
  ylab("New genes")

myplot = paste0(outpath, spp, ".boxplot_curve.png")
ggsave(myplot, plot=pl)


