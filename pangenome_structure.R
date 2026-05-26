# Calcualte pangenome structure from presence-absence matrix
#
# Written by Carlos Valiente-Mullor
#
# Usage:
#       Rscript pangenome_structure.R [sp_name]
#
# Input file name 'gene_presence_absence_roary.{species_name}.csv' is expected by default
#
# Example:
#       Rscript pangenome_structure.R Enterococcus_faecium


## packages
library(pagoo)
library(micropan)
library(matrixStats)

## Arguments
args = commandArgs(trailingOnly=TRUE)
spp = args[1] # species name ("Genus_epithet")

## Paths
inpath = "path/to/input/dir"   # set input dir
outpath= "path/to/output/dir"  # set output dir

## Read input / pan-matrix
myfile_r = paste0(inpath, "gene_presence_absence_roary.", spp, ".csv") # Roary-formatted PAM from Panaroo
pg = roary_2_pagoo(gene_presence_absence_csv = myfile_r) # Panaroo table to pagoo object

N_spp = dim(pg$pan_matrix)[1]

## Pangenome info/structure (alt.)
bin = pg$pan_matrix
bin_gfreq = colSums(bin)

dist_gfreq = table(bin_gfreq)

categ1 = as.numeric(names(dist_gfreq)) # N genomes categ. with n genes > 0
categ2 = rep(0, N_spp) # 0s vector for all possible N genomes categ.
idx = 0
for (i in categ1) {
  idx = idx + 1
  categ2[i] = dist_gfreq[idx]
}

# dataframe with all N genomes categ. + nº of genes (including 0s)
df.gfd = as.data.frame(cbind(categ2, seq(1, N_spp, by=1)))
colnames(df.gfd) = c("n_genes", "N_genomes")

# write table
utable1 = paste0(outpath, spp, ".utable.csv")
write.csv(df.gfd, utable1, row.names = F)


# Thresholds
core100 = N_spp
core95 = round(N_spp*0.95, 0)
core90 = round(N_spp*0.90, 0)
cloud5 = round(N_spp*0.05, 0)
unique = 1

# Row selection for each fraction. Calculate pangenome structure
n.total   = sum(df.gfd[,1])
n.core100 = sum(df.gfd[core100, 1])                 # core 100%
n.core95  = sum(df.gfd[core95:core100, 1])          # core 95%
n.core90  = sum(df.gfd[core90:core100, 1])          # core 90%

if (cloud5 <= 1) {
  n.acc100  = sum(df.gfd[2:(core100-1), 1])  # acc 100%
  n.acc95   = sum(df.gfd[2:(core95-1), 1])   # acc 95%
  n.acc90   = sum(df.gfd[2:(core90-1), 1])   # acc 90%

  n.cloud5  = 0                              # cloud 5% (excluding unique) 
  
} else {
  n.acc100  = sum(df.gfd[(cloud5+1):(core100-1), 1])  # acc 100%
  n.acc95   = sum(df.gfd[(cloud5+1):(core95-1), 1])   # acc 95%
  n.acc90   = sum(df.gfd[(cloud5+1):(core90-1), 1])   # acc 90%
  n.cloud5  = sum(df.gfd[2:cloud5, 1])                # cloud 5% (excluding unique) 
}

n.unique  = sum(df.gfd[1, 1])                         # unique strict
n.a100_cloud = n.acc100 + n.cloud5
n.a95_cloud  = n.acc95  + n.cloud5
n.a90_cloud  = n.acc90  + n.cloud5
n.cloud_uni  = n.cloud5 + n.unique


# Pangenome structure results tables
fractions = c("Total", "Core", "Accesory", "Cloud", "Unique", "Acc.+Cloud", "Cloud+Unique")
t_percents  = c("100%", "95%", "90%")
v100 = c(n.total, n.core100, n.acc100, n.cloud5, n.unique, n.a100_cloud, n.cloud_uni)
v95 = c(n.total, n.core95, n.acc95, n.cloud5, n.unique, n.a95_cloud, n.cloud_uni)
v90 = c(n.total, n.core90, n.acc90, n.cloud5, n.unique, n.a90_cloud, n.cloud_uni) 

n.dfp = as.data.frame(rbind(v100, v95, v90)) # table n absolute
rownames(n.dfp) = t_percents
colnames(n.dfp) = fractions

p.dfp = round((n.dfp/n.total)*100, 2) # table %

# write tables
outtable1 = paste0(outpath, spp, ".structure_ngenes.csv")
outtable2 = paste0(outpath, spp, ".structure_percen.csv")
write.csv(n.dfp, outtable1, row.names = T)
write.csv(p.dfp, outtable2, row.names = T)

# U-plot
gg1 = ggplot(df.gfd, aes(x=N_genomes, y=n_genes)) + geom_point() + ylab("Genes") + xlab("Genomes")
gg2 = ggplot(df.gfd, aes(x=N_genomes, y=n_genes)) + geom_bar(stat='identity') + ylab("Genes") + xlab("Genomes")

u1 = paste0(outpath, spp, ".uplot1.png")
ggsave(u1, plot=gg1, width = 7.22, height = 7.22, units = "in")

u2 = paste0(outpath, spp, ".uplot2.png")
ggsave(u2, plot=gg2, width = 7.22, height = 7.22, units = "in")

