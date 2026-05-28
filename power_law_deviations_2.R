## Estimate openness from simulated curves with deviations from the power law (II)
#     * Simulated open -> closed pangenome curve.
#
# Written by Carlos Valiente-Mullor
#
# Usage:
#       Rscript power_law_deviations_2.R
#

# Libs
library(reshape2)
library(ggplot2)
library(minpack.lm)

# Paths
outpath="path/to/output/dir/"

### Generate pangenome simulating open -> closed pangenome as N increases ###

# Model of open pangenome "closing"
C = 1000                                              # constant
a = c(0.8, 0.85, 0.9, 0.95, 1, 1.05, 1.1, 1.15, 1.2)  # range of alpha
sample_size = 100                                     # N
fixed.index = seq(sample_size)                        # vector of N values

# build Delta(n) by combining values form different curves
values = fixed.index
An = rep(NA, length(fixed.index))

An[1:11]    = C * values[1:11]^(-a[1])
An[12:22]   = C * values[12:22]^(-a[2])
An[23:33]   = C * values[23:33]^(-a[3])
An[34:44]   = C * values[34:44]^(-a[4])
An[45:55]   = C * values[45:55]^(-a[5])
An[56:66]   = C * values[56:66]^(-a[6])
An[67:77]   = C * values[67:77]^(-a[7])
An[78:88]   = C * values[78:88]^(-a[8])
An[89:100]  = C * values[89:100]^(-a[9])

# Plots to explore data
plot( fixed.index, An )            # plot rarefaction curve
plot( log(fixed.index), log(An) )  # plot log-log curve

# Remove 1st value
An = An[2:sample_size]


### Estimate alpha of the final curve ###
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


### Estimate alpha with incremental N values ###
values_lt = values[2:sample_size]
An_lt = An

res.def.med = c()
data_curve.med = list()

res.def.med.nls = c()

for (i in seq(2, length(values_lt))) {
  # Data for incremental N values
  x.med = values_lt[1:i]
  y.med = An_lt[1:i]
  
  y.med[y.med == 0] <- NA
  
  data_curve.med[[(i-1)]] = cbind(x.med, y.med)
  
  ## Heaps' law (power law) - linear model (lm)
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
  
  ## Heaps' law (power law) - exponential model (nlsLM)
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

# Write results dataset
# lm
rownames(a_perN.med) = NULL

mytable = paste0(outpath, "simulated_open2closed", ".results.alpha_perN_lm.csv")
write.csv(a_perN.med, mytable, row.names = F)

# nlsLM
rownames(a_perN.med.nls) = NULL

mytable = paste0(outpath, "simulated_open2closed", ".results.alpha_perN_nlsLM.csv")
write.csv(a_perN.med.nls, mytable, row.names = F)


## Plots

# lm
curve.med = ggplot(a_perN.med, aes(x=N, y=alpha)) + geom_point()

myplot.med = paste0(outpath, "simulated_open2closed", ".plot.alpha_perN_lm.png")
ggsave(myplot.med, plot=curve.med, width = 7.22, height = 7.22)

# nlsLM
curve.med.nls = ggplot(a_perN.med.nls, aes(x=N, y=alpha)) + geom_point()

myplot.med.nls = paste0(outpath, "simulated_open2closed", ".plot.alpha_perN_nlsLM.png")
ggsave(myplot.med.nls, plot=curve.med.nls, width = 7.22, height = 7.22)

# Log-log curve data + Alpha
log_curve.med = lapply(data_curve.med, log)

# lm
log.med.outpath = paste0(outpath, "simulated_open2closed", ".alpha_log-log.alpha_perN_lm.pdf")
list4 = list()
for (i in 1:length(log_curve.med)){
  lbl = paste0("\u03b1 = ", round(a_perN.med$alpha[i],3))
  #lbl = paste0("a = ", round(a_perN.med$alpha[i],3))
  list4[[i]] <- ggplot(as.data.frame(log_curve.med[[i]]), aes(x=x.med, y=y.med)) + 
    geom_point() + 
    geom_smooth(method="lm", se=F) +
    xlim(0, 5.5) + ylim(-0.5, 6.5) +
    annotate("text", x=5, y=5, 
             label = lbl, 
             size=7, color="black")
}

# save w/ cairo_pdf in order to print "alpha" character
grDevices::cairo_pdf(filename = log.med.outpath, onefile = T)
list4
dev.off()

# nlsLM
log.med.outpath = paste0(outpath, "simulated_open2closed", ".alpha_log-log.alpha_perN_nlsLM.pdf")
list4 = list()
for (i in 1:length(log_curve.med)){
  lbl = paste0("\u03b1 = ", round(a_perN.med.nls$alpha[i],3))
  #lbl = paste0("a = ", round(a_perN.med$alpha[i],3))
  list4[[i]] <- ggplot(as.data.frame(log_curve.med[[i]]), aes(x=x.med, y=y.med)) + 
    geom_point() + 
    geom_smooth(method="lm", se=F) +
    xlim(0, 5.5) + ylim(-0.5, 6.5) +
    annotate("text", x=5, y=5, 
             label = lbl, 
             size=7, color="black")
}

# save w/ cairo_pdf in order to print "alpha" character
grDevices::cairo_pdf(filename = log.med.outpath, onefile = T)
list4
dev.off()



