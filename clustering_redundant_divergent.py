# Mash-distance based filtering of redundant genomes and clustering at the species level (ANI~96%).
#
# Written by Carlos Valiente-Mullor
#
# Usage:
#       python3 cluster_redundant_divergent.py [sp_name]
#
# Example:
#       python3 cluster_redundant_divergent.py Enterococcus_faecium
#
# Files required:
#   * Mash distance matrix (output of 'mash triangle') obtained with script 'parse_mash_matrix.sh'
#   * Assembly quality summary in CSV format. Filename expected by default: '{species_name}.quality_summary.csv
#  
# Assembly quality CSV structure and fields:
#   - Header ("assembly_id,contigs,L90,Ns,level,reference)
#   - Each row after header includes quality data of a single unique assembly
#   - Assembly ID: unique assembly GenBank identifier (GCA_*)
#   - Contigs: # of contigs.
#   - L90: as calculated by 'assembly-stats' (https://github.com/sanger-pathogens/assembly-stats)
#   - Ns: N count as calculated by 'assembly-stats'
#   - Level: assembly level as specified in GenBank (Complete, Chromosome, Scaffold or Contig)
#   - Reference: binary field (0/1) indicating if an assembly is the species' reference genome.


## Libraries
import sys
import pandas as pd
import numpy as np
import networkx as nx


## Functions
def get_id(string):
    id = 'GCA_' + string.split("_")[1]
    return id

# Arguments
spp = str(sys.argv[1]) # species name ("Genus_epithet")


## Read Mash matrix

# input
myfile = spp + ".mash.parsed.csv"

# open/read matrix file with module "pandas"
m = pd.read_csv(myfile, delimiter=",", index_col=0) # type: dataframe

# lower to symmetric matrix (necessary for graph)
m2 = m.to_numpy()  # no rownames, numpy object

thres = 0.0001 # ANI~99.99%

m_values = np.where(m2 <= thres, 1, m)        # binary
m_values = np.where(m2 > thres, 0, m_values)  # binary

n = len(m_values)
for i in range(n):
    for j in range(n):
        if i==j: # diagonal to 0
            m_values[i,j] = 0
        else: # set upper matrix (NaN) with values in lower matrix
            m_values[i,j] = m_values[j,i]

# binary numpy array to pandas dataframe (add row and col names)
rowcol = list(m.columns.values)        # row and column names from original pandas df
rowcol_id = list(map(get_id, rowcol))  # only assembly IDs

m_sym = pd.DataFrame(data = m_values, index = rowcol_id, columns = rowcol_id)


## Read metadata and remove L90>100 from matrix
#.................. --> !!!!!! [indicar estructura de "tablafinal.csv" (ver script)]

# input
myfile = spp + ".quality_summary.csv"

# open/read matrix file with module "pandas"
meta = pd.read_csv(myfile, delimiter=",", header=0) # type: dataframe

# Replace values (.loc function)
meta.loc[meta.level == "Complete", "contigs"] = 1
meta.loc[meta.level == " Chromosome", "contigs"] = 1
meta.loc[meta.level == "Complete", "L90"] = 1
meta.loc[meta.level == " Chromosome", "L90"] = 1
meta.loc[meta.level == "Complete", "contigs"] = 1

meta.loc[meta.reference == 1, "reference"] = True
meta.loc[meta.reference == 0, "reference"] = False

meta.loc[meta.level == "Complete", "level"] = 1
meta.loc[meta.level == " Chromosome", "level"] = 2
meta.loc[meta.level == "Scaffold", "level"] = 3
meta.loc[meta.level == "Contig", "level"] = 4

# Remove sequences with L90 over 100
ids_L90 = meta.loc[~(meta["L90"] > 100), "assembly_id"].tolist()  # get IDs of sequences to keep (L90<=100)

m_sym_L90 = m_sym.loc[ids_L90, ids_L90] # clean matrix


## NetworkX graph

# make graph
gr = nx.from_pandas_adjacency(m_sym_L90)

# write gml graph file
gname = spp + ".graph_redundant.gml"
nx.write_gml(gr, gname)

# iterable object with connected components
gr_cc = nx.connected_components(gr)

# Loop for subsetting metadata + sorting + selecting representatives
sortby = ["contigs", "L90", "Ns", "level"] # columns to sort by

representative = []
clustn = []
i = 0

for id_list in gr_cc:
    i += 1

    sub_meta = meta[meta["assembly_id"].isin(id_list)] # subset metadata table by genomes in a component
    sub_meta_sorted = sub_meta.sort_values(sortby) # sort subset table by columns
    sub_meta_sorted.index = list(range(len(sub_meta_sorted))) # remove original indexes (subset keeps them), replace by new indexes starting at 0

    # add representative + cluster nº to lists
    clustn.append(i)
    representative.append(sub_meta_sorted.iloc[0, 0]) # access 1st row, 1st col (ID of the best ranked genome))

    # include spp reference
    if any(sub_meta_sorted.reference):  # checks if there is any "True" value in iterable
        if len(sub_meta_sorted)>1:      # checks if there is >1 element in component [no tiene sentido seguir si solo hay 1, que será la referencia]
            ref_idx = sub_meta_sorted.index[sub_meta_sorted["reference"]].tolist()[0]

            if ref_idx != 0:
                first_qual = sub_meta_sorted.iloc[0, 1:5].tolist()
                ref_qual = sub_meta_sorted.iloc[ref_idx, 1:5].tolist()

                if first_qual == ref_qual:
                    representative.pop()
                    representative.append(sub_meta_sorted.iloc[ref_idx, 0])

                else:
                    clustn.append(i)
                    representative.append(sub_meta_sorted.iloc[ref_idx, 0])


## write list of kept genomes after redundance filter

newfile = spp + ".kept.txt" 

with open(newfile, 'w') as f:
    for line in representative:
        f.write("%s\n" % line)


## Keep only "representative" (non-redundant) genomes in binary matrix

thres = 0.04 # ANI~96%

m_values2 = np.where(m2 <= thres, 1, m)         # binary
m_values2 = np.where(m2 > thres, 0, m_values2)  # binary

n = len(m_values2)
for i in range(n):
    for j in range(n):
        if i==j: # diagonal to 0
            m_values2[i,j] = 0
        else: # set upper matrix (NaN) with values in lower matrix
            m_values2[i,j] = m_values2[j,i]

# binary numpy array to pandas dataframe (add row and col names)
#rowcol = list(m.columns.values)        # row and column names from original pandas df
#rowcol_id = list(map(get_id, rowcol))  # only assembly IDs

m_sym2 = pd.DataFrame(data = m_values2, index = rowcol_id, columns = rowcol_id)

# Partimos de la nueva matriz binaria de mash con el umbral de divergencia (M>0.04; ANI<96%)
m_nored = m_sym2.loc[representative, representative] # clean matrix


## NetworkX graph w/o redundant genomes
# make graph
gr_nored = nx.from_pandas_adjacency(m_nored)

# write gml graph file
gname = spp + ".graph_divergent.gml"
nx.write_gml(gr_nored, gname)

# iterable object with connected components
gr_nored_cc = nx.connected_components(gr_nored)

cn = 0
for id_list in gr_nored_cc:
    cn += 1
    newfile = spp + ".cluster" + str(cn) + ".mash.txt"
    
    with open(newfile, mode='w') as myfile:
        myfile.write('\n'.join(id_list))
