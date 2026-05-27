## Estimate openness from simulated curves with deviations from the power law (I)
#     * Add constant distribution to the tail of a power law curve
#
# Written by Carlos Valiente-Mullor
#
# Usage:
#       Rscript power_law_deviations_1.R.R
#

# Libs
library(reshape2)
library(ggplot2)
library(minpack.lm)

# Paths
outpath="path/to/output/dir"

### 1. Generate power-law distribution ###

# Parameters of the model (power-law)
# C=1000, a=1.5, N=100 -> last Delta(n) value=1
C = 1000   # constant
a = 1.5
sample_size = 100 # N
fixed.index = seq(sample_size) # vector of N values

# ==> Values to simulate open pangenome <==
# a = 0.8
# sample_size = 4000
# fixed.index = seq(sample_size) # vector of N values

# Compute An values with power-law formula
values = fixed.index
An = C * values^(-a)

# Plots to explore data
plot( fixed.index, An ) # plot rarefaction curve
plot( log(fixed.index), log(An) ) # plot log-log curve

# Remove 1st value
An = An[2:sample_size]

## 1.2 Estimate alpha from original curve ##
x.all = values[2:sample_size]
y.all = An

# Linearized method
fitHeaps <- lm(log(y.all) ~ log(x.all), na.action = na.exclude)
linear_hl.all = summary(fitHeaps)
linear_hl.all.alpha = -1*(linear_hl.all$coefficients[2])

linear_hl.all.alpha

# Minimum squares method
control1 = nls.control(maxiter=1000, tol=1e-05, warnOnly = T)
fitM = try( nlsLM(y.all ~ k*x.all^power, start = list(power = 1, k = 1),
                control = control1, na.action=na.exclude))
model_result = summary(fitM)
model_result.alpha = -1*(model_result$parameters[1])
model_result.se = model_result$parameters[3]

model_result.alpha

### 2. Add tail of Delta(n)=1 ###

tail_size = 100
An_lt = c(An, rep(1, tail_size))

values_lt = seq(2, sample_size+tail_size) # index/N including An=1 observations

# Plots to explore data
plot( values_lt, An_lt )           # plot rarefaction curve
plot( log(values_lt), log(An_lt) ) # plot log-log curve


### 3. Compute alpha with total N ###
x.all = values_lt
y.all = An_lt

# Linearized method, lm
fitHeaps <- lm(log(y.all) ~ log(x.all), na.action = na.exclude)
linear_hl.all = summary(fitHeaps)
linear_hl.all.alpha = -1*(linear_hl.all$coefficients[2])

linear_hl.all.alpha

# Minimum squares method, nlsLM

control1 = nls.control(maxiter=1000, tol=1e-05, warnOnly = T)
fitM = try( nlsLM(y.all ~ k*x.all^power, start = list(power = 1, k = 1),
                control = control1, na.action=na.exclude))
model_result = summary(fitM)
model_result.alpha = -1*(model_result$parameters[1])
model_result.se = model_result$parameters[3]

model_result.alpha

### 4. Compute alpha with incremental N values ###
res.def.med = c()
data_curve.med = list()

res.def.med.nls = c()

for (i in seq(2, length(values_lt))) {
  # Data for incremental N values
  x.med = values_lt[1:i]
  y.med = An_lt[1:i]
  
  y.med[y.med == 0] <- NA
  
  data_curve.med[[(i-1)]] = cbind(x.med, y.med)
  
  ## Heaps' law (power law) - linear; ln
  fitHeaps <- try( lm(log(y.med) ~ log(x.med), na.action = na.exclude) )

  if ( class(fitHeaps) != "try-error" ) {
    fit_res = summary(fitHeaps)
    fit_res.alpha = -1*(fit_res$coefficients[2])
    fit_res.se = fit_res$coefficients[4]

    # results df - alpha values
    res = c("Simulated An values (lm)", i, values_lt[i], fit_res.alpha, fit_res.se)
    names(res) = c("method", "index", "N", "alpha", "SE")

    res.def.med = rbind(res.def.med, res)
  }
  
  ## Heaps' law (power law) - nlsLM
  control1 = nls.control(maxiter=1000, tol=1e-05, warnOnly = T)
  fitM = try( nlsLM(y.med ~ k*x.med^power, start = list(power = 1, k = 1),
                  control = control1, na.action=na.exclude))
    
  if ( class(fitHeaps) != "try-error" ) {
    model_result = summary(fitM)
    model_result.alpha = -1*(model_result$parameters[1])
    model_result.se = model_result$parameters[3]
    
    # results df - alpha values
    res = c("Simulated An values (nlsLM)", i, values_lt[i], model_result.alpha, model_result.se)
    names(res) = c("method", "index", "N", "alpha", "SE")
    
    res.def.med.nls = rbind(res.def.med.nls, res)
  }
  
}


# Datasets
a_perN.med = as.data.frame(res.def.med)
a_perN.med[,2] = as.numeric(a_perN.med[,2])
a_perN.med[,3] = as.numeric(a_perN.med[,3])
a_perN.med[,4] = as.numeric(a_perN.med[,4])
a_perN.med[,5] = as.numeric(a_perN.med[,5])

a_perN.med.nls = as.data.frame(res.def.med.nls)
a_perN.med.nls[,2] = as.numeric(a_perN.med.nls[,2])
a_perN.med.nls[,3] = as.numeric(a_perN.med.nls[,3])
a_perN.med.nls[,4] = as.numeric(a_perN.med.nls[,4])
a_perN.med.nls[,5] = as.numeric(a_perN.med.nls[,5])

# Write alpha results
# lm
rownames(a_perN.med) = NULL

mytable = paste0(outpath, "simulated_", a, ".results.alpha_perN_lm.csv")
write.csv(a_perN.med, mytable, row.names = F)

# nlsLM
rownames(a_perN.med.nls) = NULL

mytable = paste0(outpath, "simulated_", a, ".results.alpha_perN_nlsLM.csv")
write.csv(a_perN.med.nls, mytable, row.names = F)


## Plots
# plot alpha-vs-N

# lm
curve.med = ggplot(a_perN.med, aes(x=N, y=alpha)) + geom_point()

myplot.med = paste0(outpath, "simulated_", a, ".plot.alpha_perN_lm.png")
ggsave(myplot.med, plot=curve.med, width = 7.22, height = 7.22)

# nlsLM
curve.med.nls = ggplot(a_perN.med.nls, aes(x=N, y=alpha)) + geom_point()

myplot.med.nls = paste0(outpath, "simulated_", a, ".plot.alpha_perN_nlsLM.png")
ggsave(myplot.med.nls, plot=curve.med.nls, width = 7.22, height = 7.22)


# Log-log curve data + Alpha
log_curve.med = lapply(data_curve.med, log)

# lm
log.med.outpath = paste0(outpath, "simulated_", a, ".alpha_log-log.alpha_perN_lm.pdf")
list4 = list()
for (i in 1:length(log_curve.med)){
  lbl = paste0("\u03b1 = ", round(a_perN.med$alpha[i],3)) +
  list4[[i]] <- ggplot(as.data.frame(log_curve.med[[i]]), aes(x=x.med, y=y.med)) + 
    geom_point() + 
    geom_smooth(method="lm", se=F) +
    xlim(0, 9) + ylim(-0.5, 6.5) +
    annotate("text", x=5, y=5, 
             label = lbl, 
             size=7, color="black")
}

# save w/ cairo_pdf in order to print "alpha" character
grDevices::cairo_pdf(filename = log.med.outpath, onefile = T)
list4
dev.off()

# nlsLM
log.med.outpath = paste0(outpath, "simulated_", a, ".alpha_log-log.alpha_perN_nlsLM.pdf")
list4 = list()
for (i in 1:length(log_curve.med)){
  lbl = paste0("\u03b1 = ", round(a_perN.med.nls$alpha[i],3)) +
  list4[[i]] <- ggplot(as.data.frame(log_curve.med[[i]]), aes(x=x.med, y=y.med)) + 
    geom_point() + 
    geom_smooth(method="lm", se=F) +
    xlim(0, 9) + ylim(-0.5, 6.5) +
    annotate("text", x=5, y=5, 
             label = lbl, 
             size=7, color="black")
}

# (save w/ cairo_pdf in order to print "alpha" character)
grDevices::cairo_pdf(filename = log.med.outpath, onefile = T)
list4
dev.off()



