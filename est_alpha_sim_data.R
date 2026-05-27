## Estimate openness from power-laws with predefined alpha (6 methods)
#
# Written by Carlos Valiente-Mullor
#
# Usage:
#       Rscript est_alpha_sim_data.R
#


# Packages
library(reshape2)
library(minpack.lm)

## Simulate power laws using a rage of alpha values (0.05-2, step 0.05)

# power law parameters
C = 1000   # constant
a_ran = seq(from = 0.05, to = 2, by = 0.05) # range of alpha values
g_ran = 1-a_ran
sample_size = 100 # N

# empty dataset
a_ran.fixed = data.frame(matrix(ncol = length(a_ran), nrow = sample_size))
fixed.index = seq(sample_size)

for (i in 1:length(a_ran)) {
  values = fixed.index
  #n = C * values^(g) # Heap's law to compute number of genes (n)
  An = C * values^(-a_ran[i]) # power-law to compute n increment [Delta(n)]
  a_ran.fixed[,i] = An # add delta(n) to next column.
}

colnames(a_ran.fixed) = as.character(a_ran) # set alpha values as column names

a_ran.fixed_2melt = melt(a_ran.fixed) # convert dataframe
fixed.index_2melt = rep(fixed.index, length(a_ran))
a_ran.fixed_2melt[,1] <- fixed.index_2melt

# exploratory plots
plot(a_ran.fixed_2melt)     
plot(log(a_ran.fixed_2melt))

plot( fixed.index, a_ran.fixed[, "0.8"] ) # plot specific curve
plot( log(fixed.index), log(a_ran.fixed[, "0.8"]) ) # plot log-log specific curve


######################
### Estimate alpha ###
######################

## METHOD 1 pagoo (adapted) ##

g_ran.fixed = data.frame(matrix(ncol = length(g_ran), nrow = sample_size))
fixed.index = seq(sample_size)

for (i in 1:length(g_ran)) {
  values = fixed.index
  n = C * values^(g_ran[i]) # Heap's law to compute n
  g_ran.fixed[,i] = n # add values to next column
}

colnames(g_ran.fixed) = as.character(g_ran) # set alpha values as column names

g_ran.fixed_2melt = melt(g_ran.fixed) # convert dataframe
fixed.index_2melt = rep(fixed.index, length(g_ran))
g_ran.fixed_2melt[,1] <- fixed.index_2melt

# exploratory plots
plot( fixed.index, g_ran.fixed[, "0.8"] ) # plot specific curve
plot( log(fixed.index), log(g_ran.fixed[, "0.8"]) ) # plot log-log specific curve
plot(g_ran.fixed_2melt) # plot each curve
plot(log(g_ran.fixed_2melt))
plot(log10(g_ran.fixed_2melt))


# Adapt pagoo::rarefact (skip permutations)
# Empty df, cols = expected/real alpha, estimated alpha, difference, idx ...
# ... rows: 1 per each alpha value tested
g_est = data.frame(matrix(ncol = 7, nrow = length(g_ran)))

x.all = fixed.index
for (i in 1:length(g_ran)) {
  idx=i
  
  y.all = g_ran.fixed[,idx] # delta(n) from corresponding alpha value
  
  fitHeaps <- lm(log(y.all) ~ log(x.all), na.action = na.exclude)
  linear_hl.all = summary(fitHeaps)
  linear_hl.all.gamma = linear_hl.all$coefficients[2]
  linear_hl.all.se = linear_hl.all$coefficients[4]
  
  linear_hl.all.alpha = 1-(linear_hl.all$coefficients[2])
  
  # Set columns gamma
  g_est[idx, 1] = g_ran[idx]                      # expected("real") alpha
  g_est[idx, 2] = linear_hl.all.gamma             # estimated alpha
  g_est[idx, 3] = g_ran[idx]-linear_hl.all.gamma  # difference
  g_est[idx, 4] = idx                             # index
  
  # Set columns alpha
  g_est[idx, 5] = a_ran[idx]                      # expected("real") alpha
  g_est[idx, 6] = linear_hl.all.alpha             # estimated alpha
  g_est[idx, 7] = a_ran[idx]-linear_hl.all.alpha  # difference
}

relative_error = g_est[,3]/g_est[,1] # relative error
sum.relative_error = c(mean(relative_error), 
                       median(relative_error),
                       min(relative_error),
                       max(relative_error),
                       sd(relative_error))
names(sum.relative_error) = c("mean", "median", "min", "max", "SD")

# plots alpha real -vs- estimated
plot(g_est[,1], g_est[,1]) # expected (expected alpha -vs- expect alpha)
plot(g_est[,1], g_est[,2]) # estimated (expected alpha -vs- observed alpha)

y.comp1 = c(g_est[,1], g_est[,2]) # combine values to compare the 2 curves
x.comp1 = c(g_est[,1], g_est[,1])
plot(x.comp1, y.comp1)


## METHOD 2 micropan::heaps ##

# "objectFun" function (micropan)
objectFun <- function(p, x, y){
  y.hat <- p[1] * x^(-p[2])
  J <- sqrt(sum((y - y.hat)^2))/length(x)
  return(J)
}

# Empty df, cols = expected/real alpha, estimated alpha, difference, idx ...
# ... rows: 1 per each alpha value tested
a_est = data.frame(matrix(ncol = 4, nrow = length(a_ran)))

x.all = fixed.index
for (i in 1:length(a_ran)) {
  idx=i
  
  y.all = a_ran.fixed[,idx] # delta(n) from corresponding alpha value
  p0 <- c(median(y.all[which(x.all == 2)] ), 1)
  
  # fit
  fit <- optim(p0, objectFun, gr = NULL, x.all, y.all, method = "L-BFGS-B", lower = c(0, 0), upper = c(10000, 2))
  a.hat <- fit$par[2]
  
  # Set columns
  a_est[idx, 1] = a_ran[idx]               # expected("real") alpha
  a_est[idx, 2] = a.hat                    # estimated alpha
  a_est[idx, 3] = a_ran[idx]-a.hat         # difference
  a_est[idx, 4] = idx                      # index
}

relative_error = a_est[,3]/a_est[,1] # relative error
sum.relative_error = c(mean(relative_error), 
                       median(relative_error),
                       min(relative_error),
                       max(relative_error),
                       sd(relative_error))
names(sum.relative_error) = c("mean", "median", "min", "max", "SD")

# plots alpha real -vs- estimated
plot(a_est[,1], a_est[,1]) #  expected (expected alpha -vs- expect alpha)
plot(a_est[,1], a_est[,2]) # estimated (expected alpha -vs- observed alpha)

y.comp1 = c(a_est[,1], a_est[,2]) # combine values to compare the 2 curves
x.comp1 = c(a_est[,1], a_est[,1])
plot(x.comp1, y.comp1)


## METHOD 3 (newman et al. 2005) ##

# Empty df, cols = expected/real alpha, estimated alpha, difference, idx ...
# ... rows: 1 per each alpha value tested
a_est = data.frame(matrix(ncol = 4, nrow = length(a_ran)))

n = sample_size
for (i in 1:length(a_ran)) {
  idx=i
  
  y = a_ran.fixed[,idx]
  
  # Formula
  ml.alpha = try( 1 + n*((sum(log(y/y[1]))))^(-1) )
  
  # Set columns
  a_est[idx, 1] = a_ran[idx]               # expected("real") alpha
  a_est[idx, 2] = ml.alpha                 # estimated alpha
  a_est[idx, 3] = a_ran[idx]-ml.alpha      # difference
  a_est[idx, 4] = idx                      # index
}

relative_error = a_est[,3]/a_est[,1] # relative error
sum.relative_error = c(mean(relative_error), 
                       median(relative_error),
                       min(relative_error),
                       max(relative_error),
                       sd(relative_error))
names(sum.relative_error) = c("mean", "median", "min", "max", "SD")

# plots alpha real -vs- estimated
plot(a_est[,1], a_est[,1]) # expected (expected alpha -vs- expect alpha)
plot(a_est[,1], a_est[,2]) # estimated (expected alpha -vs- observed alpha)

y.comp1 = c(a_est[,1], a_est[,2]) # combine values to compare the 2 curves
x.comp1 = c(a_est[,1], a_est[,1])
plot(x.comp1, y.comp1)


## METHOD 4: Heaps' law (power law) - exponential (nlsLM) ##

# Empty df, cols = expected/real alpha, estimated alpha, difference, idx ...
# ... rows: 1 per each alpha value tested
a_est = data.frame(matrix(ncol = 4, nrow = length(a_ran)))
x.all = fixed.index
control1 = nls.control(maxiter=1000, tol=1e-05, warnOnly = T)

for (i in 1:length(a_ran)) {
  idx=i
  
  y.all = a_ran.fixed[,idx] # delta(n) from corresponding alpha value
  
  mod.nls.all = try( nlsLM(y.all ~ k*x.all^power, start = list(power = 1, k = 1),
                           control = control1, na.action=na.exclude))
  
  nls_all = summary(mod.nls.all)
  nls_all.alpha = -1*(nls_all$parameters[1])
  nls_all.se = nls_all$parameters[3]
  
  # Set columns
  a_est[idx, 1] = a_ran[idx]               # expected("real") alpha
  a_est[idx, 2] = nls_all.alpha            # estimated alpha
  a_est[idx, 3] = a_ran[idx]-nls_all.alpha # difference
  a_est[idx, 4] = idx                      # index
}

a_est

relative_error = a_est[,3]/a_est[,1] # relative error
sum.relative_error = c(mean(relative_error), 
                       median(relative_error),
                       min(relative_error),
                       max(relative_error),
                       sd(relative_error))
names(sum.relative_error) = c("mean", "median", "min", "max", "SD")

# plots alpha real -vs- estimated
plot(a_est[,1], a_est[,1]) # expected (expected alpha -vs- expect alpha)
plot(a_est[,1], a_est[,2]) # estimated (expected alpha -vs- observed alpha)

y.comp1 = c(a_est[,1], a_est[,2]) # combine values to compare the 2 curves
x.comp1 = c(a_est[,1], a_est[,1])
plot(x.comp1, y.comp1)


## METHOD 5 Heaps' law (power law) - linear (ln) ##

# Empty df, cols = expected/real alpha, estimated alpha, difference, idx ...
# ... rows: 1 per each alpha value tested
a_est = data.frame(matrix(ncol = 4, nrow = length(a_ran)))

x.all = fixed.index
for (i in 1:length(a_ran)) {
  idx=i
  
  y.all = a_ran.fixed[,idx] # delta(n) from corresponding alpha value
  
  fitHeaps <- lm(log(y.all) ~ log(x.all))
  linear_hl.all = summary(fitHeaps)
  linear_hl.all.alpha = -1*(linear_hl.all$coefficients[2])
  linear_hl.all.se = linear_hl.all$coefficients[4]
  
  # Set columns
  a_est[idx, 1] = a_ran[idx]                      # expected("real") alpha
  a_est[idx, 2] = linear_hl.all.alpha             # estimated alpha
  a_est[idx, 3] = a_ran[idx]-linear_hl.all.alpha  # difference
  a_est[idx, 4] = idx                             # index
}

relative_error = a_est[,3]/a_est[,1] # relative error
sum.relative_error = c(mean(relative_error), 
                       median(relative_error),
                       min(relative_error),
                       max(relative_error),
                       sd(relative_error))
names(sum.relative_error) = c("mean", "median", "min", "max", "SD")

# plots alpha real -vs- estimated
plot(a_est[,1], a_est[,1]) # expected (expected alpha -vs- expect alpha)
plot(a_est[,1], a_est[,2]) # estimated (expected alpha -vs- observed alpha)

y.comp1 = c(a_est[,1], a_est[,2]) # combine values to compare the 2 curves
x.comp1 = c(a_est[,1], a_est[,1])
plot(x.comp1, y.comp1)


## METHOD 6 Heaps' linear [ln(x+1)] ##

# Empty df, cols = expected/real alpha, estimated alpha, difference, idx ...
# ... rows: 1 per each alpha value tested
a_est = data.frame(matrix(ncol = 4, nrow = length(a_ran)))

x.all = fixed.index
for (i in 1:length(a_ran)) {
  idx=i
  
  y.all = a_ran.fixed[,idx] # delta(n) from corresponding alpha value
  
  fitHeaps <- lm(log1p(y.all) ~ log1p(x.all))
  linear_hl.all = summary(fitHeaps)
  linear_hl.all.alpha = -1*(linear_hl.all$coefficients[2])
  linear_hl.all.se = linear_hl.all$coefficients[4]
  
  # Set columns
  a_est[idx, 1] = a_ran[idx]                      # expected("real") alpha
  a_est[idx, 2] = linear_hl.all.alpha             # estimated alpha
  a_est[idx, 3] = a_ran[idx]-linear_hl.all.alpha  # difference
  a_est[idx, 4] = idx                             # index
}


relative_error = a_est[,3]/a_est[,1] # relative error
sum.relative_error = c(mean(relative_error), 
                       median(relative_error),
                       min(relative_error),
                       max(relative_error),
                       sd(relative_error))
names(sum.relative_error) = c("mean", "median", "min", "max", "SD")

# plots alpha real -vs- estimated
plot(a_est[,1], a_est[,1]) #  expected (expected alpha -vs- expect alpha)
plot(a_est[,1], a_est[,2]) # estimated (expected alpha -vs- observed alpha)

y.comp1 = c(a_est[,1], a_est[,2]) # combine values to compare the 2 curves
x.comp1 = c(a_est[,1], a_est[,1])
plot(x.comp1, y.comp1)


