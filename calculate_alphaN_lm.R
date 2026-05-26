# Estimate power law's alpha (openness) as a function of N using linear model fit
#
# Written by Carlos Valiente-Mullor
#
# Usage:
#       Rscript calculate_alphaN_lm.R [sp_name]
#
# Input file name 'gene_presence_absence_roary.{species_name}.csv' is expected by default
#
# Example:
#       Rscript calculate_alphaN_lm.R Enterococcus_faecium


## packages
library(pagoo)
library(micropan)
library(matrixStats)
library(stringr)

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


########################################
#   Estimate alpha(N) - Linear model   #
########################################

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

data_curve.all = list()
data_curve.mean = list()
data_curve.med = list()

nmat.or <- powerlaw_perm(pan.matrix, n.perm)

# Add the next genome (rarefaction) each iteration
for (i in 1:length(N_values)) {
  
  # Random sub-sample of permutations matrix
  if (i < length(N_values)) {
    nmat = nmat.or[seq(N_values[i]-1),] # N-1 because we ignore 1st value
    
  } else if (N_values[i]-1 == N-1) {
    nmat = nmat.or
  }
  
  
  #### All permutations ####
  
  # Prepare data (values==0 --> NA)
  x.all <- rep((2:N_values[i]), times = n.perm)
  y.all <- as.numeric(nmat)
  y.all[y.all == 0] <- NA
  
  # add to object/matrix to plot An curves of every run (N value)
  data_curve.all[[i]] = cbind(x.all, y.all)
  
  # Fit model
  fitHeaps <- try( lm(log(y.all) ~ log(x.all), na.action = na.exclude) )
  
  # Extract alpha, prepare output df
  if ( class(fitHeaps) != "try-error" ) {
    fit_res = summary(fitHeaps)
    fit_res.alpha = -1*(fit_res$coefficients[2])
    fit_res.se = fit_res$coefficients[4]
    
    res = c("[all,lm,exclude]", i, N_values[i], fit_res.alpha, fit_res.se)
    names(res) = c("method", "index", "N", "alpha", "SE")
    
    res.def.all = rbind(res.def.all, res)
  }
  
  
  #### Mean of permutations ####
  
  # Prepare data (values==0 --> NA)
  nmat_mean = rowMeans(nmat)
  x.mean <- rep((2:N_values[i]), times = 1)
  y.mean <- as.numeric(nmat_mean)
  y.mean[y.mean == 0] <- NA
  
  # add to object/matrix to plot An curves of every run (N value)
  data_curve.mean[[i]] = cbind(x.mean, y.mean)
  
  # Fit model
  fitHeaps <- try( lm(log(y.mean) ~ log(x.mean), na.action = na.exclude) )
  
  # Extract alpha, prepare output df
  if ( class(fitHeaps) != "try-error" ) {
    fit_res = summary(fitHeaps)
    fit_res.alpha = -1*(fit_res$coefficients[2])
    fit_res.se = fit_res$coefficients[4]
    
    res = c("[mean,lm,exclude]", i, N_values[i], fit_res.alpha, fit_res.se)
    names(res) = c("method", "index", "N", "alpha", "SE")
    
    res.def.mean = rbind(res.def.mean, res)
  }
  
  #### Median of permutations ####
  
  # Prepare data (values==0 --> NA)
  nmat_med = rowMedians(nmat)
  x.med <- rep(2:N_values[i], times = 1)
  y.med <- as.numeric(nmat_med)
  y.med[y.med == 0] <- NA
  
  # add to object/matrix to plot An curves of every run (N value)
  data_curve.med[[i]] = cbind(x.med, y.med)
  
  # Fit model
  fitHeaps <- try( lm(log(y.med) ~ log(x.med), na.action = na.exclude) )
  
  if ( class(fitHeaps) != "try-error" ) {
    fit_res = summary(fitHeaps)
    fit_res.alpha = -1*(fit_res$coefficients[2])
    fit_res.se = fit_res$coefficients[4]
    
    # results df - alpha values
    res = c("[med,lm,exclude]", i, N_values[i], fit_res.alpha, fit_res.se)
    names(res) = c("method", "index", "N", "alpha", "SE")
    
    res.def.med = rbind(res.def.med, res)
  }
  
}

# Datasets
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


### Write estimate results:
rownames(a_perN.all) = NULL
rownames(a_perN.mean) = NULL
rownames(a_perN.med) = NULL

mytable = paste0(outpath, spp, ".perm_all.Na-ex.alpha_perN.csv")
write.csv(a_perN.all, mytable, row.names = F)

mytable = paste0(outpath, spp, ".perm_mean.Na-ex.alpha_perN.csv")
write.csv(a_perN.mean, mytable, row.names = F)

mytable = paste0(outpath, spp, ".perm_med.Na-ex.alpha_perN.csv")
write.csv(a_perN.med, mytable, row.names = F)


# Data frames
df.3methods = cbind(a_perN.all$N, a_perN.all$alpha, a_perN.mean$alpha, a_perN.med$alpha)
colnames(df.3methods) = c("N", "alpha.all", "alpha.mean", "alpha.median")

mytable = paste0(outpath, spp, ".perm_summ.Na-ex.alpha_perN.csv")
write.csv(df.3methods, mytable, row.names = F)


#############
##  PLOTS  ##
#############

# Save alpha curve plots
curve.all = ggplot(df.3methods, aes(x=N, y=alpha.all)) + geom_point() + geom_smooth(se=F)
curve.mean = ggplot(df.3methods, aes(x=N, y=alpha.mean)) + geom_point() + geom_smooth(se=F)
curve.med = ggplot(df.3methods, aes(x=N, y=alpha.median)) + geom_point() + geom_smooth(se=F)

myplot.all = paste0(outpath, spp, ".curve_all.Na-ex.alpha_perN.png")
myplot.mean = paste0(outpath, spp, ".curve_mean.Na-ex.alpha_perN.png")
myplot.med = paste0(outpath, spp, ".curve_med.Na-ex.alpha_perN.png")

ggsave(myplot.all, plot=curve.all)
ggsave(myplot.mean, plot=curve.mean)
ggsave(myplot.med, plot=curve.med)

# Data frame data curve
data_curve.all = lapply(data_curve.all, as.data.frame)
data_curve.mean = lapply(data_curve.mean, as.data.frame)
data_curve.med = lapply(data_curve.med, as.data.frame)

# Log-log curve data
log_curve.all = lapply(data_curve.all, log)
log_curve.mean = lapply(data_curve.mean, log)
log_curve.med = lapply(data_curve.med, log)

# Plot rarefaction curve for each total N value:
rar.all.outpath = paste0(outpath, spp, ".raref_all.Na-ex.alpha_perN.pdf")
rar.mean.outpath = paste0(outpath, spp, ".raref_mean.Na-ex.alpha_perN.pdf")
rar.med.outpath = paste0(outpath, spp, ".raref_med.Na-ex.alpha_perN.pdf")

pdf(file = rar.all.outpath)
for (i in 1:length(data_curve.all)){
  plot(data_curve.all[[i]])
}
dev.off()

pdf(file = rar.mean.outpath)
for (i in 1:length(data_curve.mean)){
  plot(data_curve.mean[[i]])
}
dev.off()

pdf(file = rar.med.outpath)
for (i in 1:length(data_curve.med)){
  plot(data_curve.med[[i]])
}
dev.off()


# Plot log-log rarefaction curve for each total N value:

## ALL perms

# Get x, y lims
# ... max and min log(N), log(An) +/- 0.5
last.v = length(log_curve.all)
x_lim.max = max(log_curve.all[[last.v]]$x.all, na.rm = T)+0.5
x_lim.min = min(log_curve.all[[last.v]]$x.all, na.rm = T)-0.5

y_lim.max = max(log_curve.all[[last.v]]$y.all, na.rm = T)+0.5
y_lim.min = min(log_curve.all[[last.v]]$y.all, na.rm = T)-0.5


# Plot log-log rarefaction curve for each total N value:
log.all.outpath = paste0(outpath, spp, ".alpha.log-log_all.Na-ex.alpha_perN.pdf")

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


# Plot log-log rarefaction curve for each total N value:
log.mean.outpath = paste0(outpath, spp, ".alpha.log-log_mean.Na-ex.alpha_perN.pdf")

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


# Get x, y lims
# ... max and min log(N), log(An) +/- 0.5
last.v = length(log_curve.med)
x_lim.max = max(log_curve.med[[last.v]]$x.med, na.rm = T)+0.5
x_lim.min = min(log_curve.med[[last.v]]$x.med, na.rm = T)-0.5

y_lim.max = max(log_curve.med[[last.v]]$y.med, na.rm = T)+0.5
y_lim.min = min(log_curve.med[[last.v]]$y.med, na.rm = T)-0.5


# Plot log-log rarefaction curve for each total N value:
log.med.outpath = paste0(outpath, spp, ".alpha.log-log_med.Na-ex.alpha_perN.pdf")

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

