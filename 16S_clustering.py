# Filter by 16S divergence (<98.5%) and aln length
#
# Written by Carlos Valiente-Mullor
#
# Usage:
#       python3 16S_cluster.py [sp_name]
#
# Example:
#       python3 16S_cluster.py Enterococcus_faecium
#
# Files required:
#	· % identity CSV matrix between 16S gene of each genome (BLAST): {sp_name}.16S_percid.csv
#	· Alignment length CSV matrix between 16S gene of each genome (BLAST): {sp_name}.16S_alnlen.csv
#	· List of assembly IDs with no 16S gene sequence found (barrnap): {sp_name}.nohit.txt"


## Libraries
import sys
import os
import pandas as pd
import numpy as np
import networkx as nx


## Arguments
spp = str(sys.argv[1]) # species name ("Genus_epithet")


## Read tables

# Read table % id
myfile = spp + ".16S_percid.csv"  # input file
mper_temp = pd.read_csv(myfile, delimiter=",", index_col=0) # open/read matrix file with module "pandas"

# Read table aln length
myfile = spp + ".16S_alnlen.csv"  # input file
maln_temp = pd.read_csv(myfile, delimiter=",", index_col=0) # open/read matrix file with module "pandas"


## Remove genomes in "no hit" list (no 16S detected)
myfile = spp + ".nohit.txt"

if os.stat(myfile).st_size != 0: # check if nohit.txt is empty (no genomes to remove)
    with open(myfile, 'r') as f:
        noh = [line.strip() for line in f.readlines()] # lines of the file to list w/o newlines
    
    kept_ids = [id for id in mper_temp.index if id not in noh] # seqs to keep (not in "no hit" list)

    # Remove seqs w/o hits in percid table
    mper_c1 = mper_temp.loc[kept_ids, kept_ids]

    # Remove seqs w/o hits in alnlen table
    maln_c1 = maln_temp.loc[kept_ids, kept_ids]

else:
    mper_c1 = mper_temp
    maln_c1 = maln_temp


## Remove genomes with short 16S (<=1300bp)
diag_len = np.diag(maln_c1) # diagonal of aln table
long_idx = np.where(diag_len>1300)[0] # indices of seqs with >1300nt

# Remove seqs <=1300nt from percid matrix
mper_c2 = mper_c1.iloc[long_idx, long_idx]

# List of rejected by aln length (compare original index list to filtered index list, keep the absent)
rejected_len = [id for id in mper_c1.index if id not in mper_c2.index]

newfile = spp + ".16S_len_reject.txt"
with open(newfile, mode='w') as nf:
    nf.write('\n'.join(rejected_len))


# write new % id matrix after removin rejected sequences
newfile = spp + ".16S_percid.parsed.csv"
mper_c2.to_csv(newfile)


## Genome clustering by 16S (threshold=98.65% identity)

# Binary matrix
mper_c2_temp = mper_c2.to_numpy()  # no rownames, numpy object

thres = 98.65

mper_values = np.where(mper_c2_temp >= thres, 1, mper_c2)    # 1 if greater/equal than threshold
mper_values = np.where(mper_c2_temp < thres, 0, mper_values) # 0 if lower than threshold

n = len(mper_values)
for i in range(n):
    for j in range(n):
        if i==j: # diagonal to 0
            mper_values[i,j] = 0

# binary numpy array to pandas dataframe (add row and col names)
rowcol = mper_c2.index  # row and column names from original pandas df

mper_c2_bin = pd.DataFrame(data = mper_values, index = rowcol, columns = rowcol)


## NetworkX graph

# make graph
gr_mper = nx.from_pandas_adjacency(mper_c2_bin)

# write gml graph file
gname = spp + ".16S=" + str(thres) + ".gml"
nx.write_gml(gr_mper, gname)

# Save file with genome IDs in each cluster (component)
gr_cc = nx.connected_components(gr_mper)

cn = 0 # cluster counter
list_sizes = []  # cluster sizes
list_numbrs = [] # cluster numbers
for id_list in gr_cc:
    cn += 1
    list_numbrs.append(cn)

    csize = len(id_list)
    list_sizes.append(csize)

    newfile = spp + ".16S=" + str(thres) + ".cluster" + str(cn) + ".csize=" + str(csize) + ".txt"
    with open(newfile, mode='w') as myfile:
        myfile.write('\n'.join(id_list))

# Write summary of cluster sizes
cdict = {"cluster": list_numbrs, "size": list_sizes}
cdf = pd.DataFrame(cdict, columns=["cluster", "size"])

newfile = spp + ".16S=" + str(thres) + ".graph_csizes.csv"
cdf.to_csv(newfile, index = False)

# Write summary spp + nº of clusters
newfile = spp + ".16S=" + str(thres) + ".summ_graph_nclust.tsv"
with open(newfile, 'w') as myfile:
    spp_comp = spp + "\t" + str(cn)
    myfile.write(spp_comp)
