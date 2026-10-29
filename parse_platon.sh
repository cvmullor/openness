#! /bin/bash

spp=$1

# Check if Platon has detected/removed the described plasmids (NCBI) of a species; if not, create a new (corrected) .chromosome.fasta
# [spp].csv (NCBI metadata) should be in the same folder.

# Assembly and plasmid (replicons) info
awk -F"," '{ print $6","$10 }' $spp.csv \
 | grep "plasmid" \
 | sed 's/"//g' \
 | tr ";," "\n" \
 | grep -v -E 'chromosome|complete|contig' > only_plasmids.$spp.temp

# One line for each plasmid
echo "spp,assembly,pls_name,accession,contig" > described_plasmids.$spp.csv

awk '/^G/ { assembly = $0 } !/^G/ { print assembly","$0}' only_plasmids.$spp.temp \
 | sed 's/ plasmid/plasmid/g' \
 | sed 's/:/,/g' \
 | sed 's|\/|,|g' \
 | awk -v sp=$spp '{ print sp","$0 }' >> described_plasmids.$spp.csv

rm only_plasmids.$spp.temp

# Check if platon .tsv summary exists
echo "spp,assembly,pls_name,accession,contig,removed" > removed_plasmids.$spp.csv

awk -F"," '{ if (NR!=1) {print $2" "$5} }' described_plasmids.$spp.csv \
 | while read assid contig
    do
        dir="plasmid_contigs.$spp/platon.$spp.$assid/*.tsv"
        if [ -e $dir ]  # -e checks if anything with that name exists
            then
                grep_res=$(grep $contig $dir | awk '{ print $1 }')
                if [ -z $grep_res ] # -z checks if variables is empty (no contig match in tsv file - not detected as a plasmid)
                    then
                        dir_chr="plasmid_contigs.$spp/platon.$spp.$assid/*.chromosome.fasta"
                        dir_pls="plasmid_contigs.$spp/platon.$spp.$assid/*.plasmid.fasta"

                        # Plasmids corrected file
                        grep -A1 $contig $dir_chr >> plasmid_contigs.$spp/platon.$spp.$assid/$assid.plasmid.corrected.fasta
                        
                        # Chromosome corrected file
                        # [https://www.biostars.org/p/157811/]
                        if [ ! -e "plasmid_contigs.$spp/platon.$spp.$assid/$assid.chromosome.corrected.fasta" ]
                            then
                                cp $dir_chr plasmid_contigs.$spp/platon.$spp.$assid/$assid.chromosome.corrected.fasta

                                awk -v id=$contig '{if (substr($0,1) == ">"id) censor=1; else if (substr($0,1,1) == ">") censor=0; if (censor==0) print $0}' plasmid_contigs.$spp/platon.$spp.$assid/$assid.chromosome.corrected.fasta > plasmid_contigs.$spp/platon.$spp.$assid/$assid.chromosome.corrected.temp
                                rm plasmid_contigs.$spp/platon.$spp.$assid/$assid.chromosome.corrected.fasta
                                mv plasmid_contigs.$spp/platon.$spp.$assid/$assid.chromosome.corrected.temp plasmid_contigs.$spp/platon.$spp.$assid/$assid.chromosome.corrected.fasta

                            else
                                awk -v id=$contig '{if (substr($0,1) == ">"id) censor=1; else if (substr($0,1,1) == ">") censor=0; if (censor==0) print $0}' plasmid_contigs.$spp/platon.$spp.$assid/$assid.chromosome.corrected.fasta > plasmid_contigs.$spp/platon.$spp.$assid/$assid.chromosome.corrected.temp
                                rm plasmid_contigs.$spp/platon.$spp.$assid/$assid.chromosome.corrected.fasta
                                mv plasmid_contigs.$spp/platon.$spp.$assid/$assid.chromosome.corrected.temp plasmid_contigs.$spp/platon.$spp.$assid/$assid.chromosome.corrected.fasta                        
                        fi
                    # Create a new table with the plasmids not detected by Platon
                    grep $contig described_plasmids.$spp.csv | awk '{ print $0",yes" }' >> removed_plasmids.$spp.csv
                fi
        fi
    done
