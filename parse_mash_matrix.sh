#! /bin/bash

# Adapt Mash triangle's matrix output for analysis with script "cluster_redundant_divergent.py"
# Input file name '{species_name}.mash.csv' is expected by default

spp=$1 # species name ("Genus_epithet")

# Transpose row names
awk 'BEGIN{print ","}{ if (NR!=1) {print $1} }' $spp.mash.tsv | paste -s -d, | sed 's/,,/,/' > $spp.mash.parsed.csv
# Add NAs to diagonal
awk -F "," '{ if (NR!=1) { print $0",NA"} }' $spp.mash.tsv | sed 's/\t/,/g' >> $spp.mash.parsed.csv
