# Estimate genomic fluidity with micropan::fluidity
#
# Usage:
#       Rscript estimate_genomic_fluidity.R [sp_name]
#
# Input file name 'gene_presence_absence_roary.{species_name}.csv' is expected by default
#
# Example:
#       Rscript estimate_genomic_fluidity.R Enterococcus_faecium


# libraries
library(pagoo)
library(micropan)

## Arguments
args = commandArgs(trailingOnly=TRUE)
spp = args[1] # species name ("Genus_epithet")

## Paths
inpath = "path/to/input/dir"   # set input dir
outpath= "path/to/output/dir"  # set output dir

## Read input / pan-matrix
myfile_r = paste0(inpath, "gene_presence_absence_roary.", spp, ".csv") # Roary-formatted PAM from Panaroo
pg = roary_2_pagoo(gene_presence_absence_csv = myfile_r) # Panaroo table to pagoo object

## Binary matrix
mybin = pg$pan_matrix


########################################
#   Calculate genomic fluidity (phi)   #
########################################

# Fluidity function
compute_phi = function(x) {
  fluidity(mybin, n.sim = x)
}

# Estimate for genomes samples N_initial -> N_total
ns_values2 = seq(2, nrow(mybin), by=1)      # increment of 1 genome
phi_list2 = sapply(ns_values2, compute_phi)

phi_col2 = t(rbind(phi_list2, ns_values2))
colnames(phi_col2) = c("Mean", "Std", "Genomes.Sampled")

filename = paste0(outpath, "fluidity_est.", spp, ".tsv")
write.table(phi_col2, file = filename, quote = F, sep = "\t", row.names = F)

# Plot
pl <- ggplot(phi_col2, aes(x = Genomes.Sampled, y = Mean))

plotfile2 = paste0(outpath, "fluidity_plot.", spp, ".png")
png(plotfile2, width = 800, height = 800, units = "px")
pl + geom_ribbon(aes(x = Genomes.Sampled, y = Mean, ymin = Mean-Std, ymax = Mean+Std), fill = "salmon") +
  geom_line(size = 0.5)
dev.off()


