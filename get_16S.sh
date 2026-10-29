#! /bin/bash

spp=$1  # spp_name

# Get 16S sequences
cat $spp.GCAlist.sorted.txt | while read id;
 do
	# get 16S coords
	barrnap -q --kingdom bac $spp.fna_f2/$id*.fna \
	 | grep "16S" \
	 | awk -v gca=$id '{ print gca"\t"$0"\t"$5-$4 }' \
	 | sort -k11 >> $spp.16S_coord.tsv  # sort by length rRNAs from the same genome (11th field: end-start)
	
	# check if barrnap detects 16S
	check=$(grep $id $spp.16S_coord.tsv)
	if [ -z "$check" ]  # quotes --> variable = only one string (necessary to avoid error "too much arguments")
	 then
		echo "barrnap can't predict 16S in "$id
	else
		# Check strand (if s.start > or < s.end)
		pos1=$(grep $id $spp.16S_coord.tsv | awk '{print $5}' | head -n1)
		pos2=$(grep $id $spp.16S_coord.tsv | awk '{print $6}' | head -n1)

		if [ $pos1 -gt $pos2 ]
		then
			pos=$pos2"-"$pos1
		elif [ $pos2 -gt $pos1 ]
		then
			pos=$pos1"-"$pos2
		fi

		# Check contig and extract 16S seq
		contig=$(grep $id $spp.16S_coord.tsv | awk '{print $2}' | head -n1)
		first_c=$(head -n1 $spp.fna_f2/$id*.fna | sed 's/>//g')

		if [ $contig == $first_c ]
		then
			extractseq --sequence $spp.fna_f2/$id*.fna -reg $pos -separate -outseq $id.$spp.16S.fa
		else
			awk -v seq=$contig -v RS='>' '$1 == seq {print RS $0}' $spp.fna_f2/$id*.fna > $id.$contig.temp  # extract specific contig from fna
			extractseq --sequence $id.$contig.temp -reg $pos -separate -outseq $id.$spp.16S.fa
		fi
		
		# one line
		awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);} END {printf("\n");}' < $id.$spp.16S.fa > $id.$spp.16S.temp
		echo ">"$id > newhead.$id # header
		# change headers
		sed '1{/^$/d}' $id.$spp.16S.temp \
		| awk 'NR%2==0' \
		| paste -d'\n' newhead.$id - > $id.$spp.16S.fa
	fi
 done

# Make database
cat *.16S.fa > $spp.16S.db.fa
makeblastdb -in $spp.16S.db.fa -dbtype nucl

mkdir $spp.db
mv *.16S.db.fa* $spp.db

# Nº of sequences in database
ndb=$(grep -c ">" $spp.db/$spp.16S.db.fa)
target_thres=$(( $ndb+100 ))

# Blast queries
cat $spp.GCAlist.sorted.txt | while read id;
 do
	# Check if there are genomes whitout predicted 16S
	fafile=$id.$spp.16S.fa
	if [ ! -f "$fafile" ]  # check if file does not exist
	 then
	 	echo $id > $id.percid.col.temp
		echo $id > $id.alnlen.col.temp

		cat $spp.GCAlist.sorted.txt | sed -e 's|GCA_[0-9]*.[0-9]|NA|g' >> $id.percid.col.temp
		cat $spp.GCAlist.sorted.txt | sed -e 's|GCA_[0-9]*.[0-9]|NA|g' >> $id.alnlen.col.temp

	else
		blastn -outfmt 7 -max_hsps 1 -max_target_seqs $target_thres -db $spp.db/$spp.16S.db.fa -query $id.$spp.16S.fa -out $id.$spp.16S.blast.tsv
		
		grep -v "#" $id.$spp.16S.blast.tsv | sort -k2 | awk '{ print $2"\t"$3 }' > $id.percid.temp # % identity
		grep -v "#" $id.$spp.16S.blast.tsv | sort -k2 | awk '{ print $2"\t"$4 }' > $id.alnlen.temp # aln length

		echo $id > $id.percid.col.temp
		echo $id > $id.alnlen.col.temp

		# Check if a hit between two sequences is missing
		total_seq=$(wc -l $spp.GCAlist.sorted.txt | awk '{ print $1 }')
		total_hit=$(wc -l $id.percid.temp | awk '{ print $1 }')

		if [ $total_hit -lt $total_seq ]
		then
			cat $spp.GCAlist.sorted.txt | while read gca
			do
				percid=$(grep $gca $id.percid.temp)
				if [ -z "$percid" ]  # check if variable is empty
				then
					echo "NA" >> $id.percid_NA.temp
				else
					echo $percid | awk '{ print $2 }' >> $id.percid_NA.temp
				fi

				alnlen=$(grep $gca $id.alnlen.temp)
				if [ -z "$alnlen" ]  # check if variable is empty
				then
					echo "NA" >> $id.alnlen_NA.temp
				else
					echo $alnlen | awk '{ print $2 }' >> $id.alnlen_NA.temp
				fi

			done
			cat $id.percid_NA.temp >> $id.percid.col.temp
			cat $id.alnlen_NA.temp >> $id.alnlen.col.temp

		else
			awk '{ print $2 }' $id.percid.temp >> $id.percid.col.temp
			awk '{ print $2 }' $id.alnlen.temp >> $id.alnlen.col.temp
		fi
	fi
 done

# Matrices
echo "" > row_names.percid.temp
cat $spp.GCAlist.sorted.txt >> row_names.percid.temp
cp row_names.percid.temp row_names.alnlen.temp

for file in *.percid.col.temp  # % identity
do
	cat row_names.percid.temp \
	 | paste -d"," - $file > temp
	cp temp row_names.percid.temp
done
rm temp
mv row_names.percid.temp $spp.16S_percid.csv

for file in *.alnlen.col.temp  # aln length
do
	cat row_names.alnlen.temp \
	 | paste -d"," - $file > temp
	cp temp row_names.alnlen.temp
done
rm temp
mv row_names.alnlen.temp $spp.16S_alnlen.csv

# Organize
rm *.temp newhead.* # remove

mkdir $spp.16S.fastas
mv *.16S.fa $spp.16S.fastas

mkdir $spp.blast_outs
mv *.16S.blast.tsv $spp.blast_outs

mkdir $spp.16S_outputs
mv $spp.db $spp.16S_outputs
mv $spp.16S.fastas $spp.16S_outputs
mv $spp.blast_outs $spp.16S_outputs
mv $spp.GCAlist.sorted.txt $spp.16S_outputs
mv $spp.16S_coord.tsv $spp.16S_outputs
mv $spp.16S_percid.csv $spp.16S_outputs
mv $spp.16S_alnlen.csv $spp.16S_outputs

