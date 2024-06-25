###########################################################################
## Author: Youtao Lu <luyoutao@sas.upenn.edu>
##  
## Copyright (c) 2021-2023, Youtao Lu and Junhyong Kim, Department of Biology, University of Pennsylvania
## Copyright (c) 2021-2023, James Eberwine, Perelman School of Medicine, University of Pennsylvania
## Copyright (c) 2021-2023, David Issadore, Penn Engineering, University of Pennsylvania
## 
## This Source Code is distributed under Creative Commons Attribution License 4.0 (CC BY).
###########################################################################
## mouse strains from Mouse Genomes Project
cat Data/SangerMouse/assembly/*.fna.gz > Data/strains/17strains.chrM.fa.gz
gzip -cd Data/strains/17strains.chrM.fa.gz > Data/strains/17strains.chrM.fa
perl -i.bak -pe 's/>([^ ]+) (.*) strain (.+) mitochondrion/>\3 \2 strain \3 mitochondrion, \1/' Data/strains/17strains.chrM.fa
## manually add mm10.mito
grep '>' Data/strains/17strains.chrM.fa | awk '{ print $1 }' | tr -d '>' | parallel -j 1 -k "perl -se 'BEGIN { print \$v, qq(\t) } while (<>) { if (/^\$v/) { @F=split /\s+/; \$s .= \$F[1] } } END { print \$s, qq(\n) }' -- -v={1} Report/strains/clustalo-17strains.clustal_num" > Report/strains/clustalo-17strains_str.tsv
echo -n -e "alignment\t" >> Report/strains/clustalo-17strains_str.tsv
perl -e 'while (<>) { if (/\s+\*/) { chomp; s/^\s{17}//; $s .= $_ } } END { print $s, "\n" }' Report/strains/clustalo-17strains.clustal_num >> Report/strains/clustalo-17strains_str.tsv

## mouse species from NCBI organelle genomes
cat Data/species/Mus_*.chrM.fa > Data/species/8species.chrM.fa
## manually add mm10.mito
grep '>' Data/species/8species.chrM.fa | awk '{ print $1 }' | tr -d '>' | parallel -j 1 -k "perl -se 'BEGIN { print \$v, qq(\t) } while (<>) { if (/^\$v/) { @F=split /\s+/; \$s .= \$F[1] } } END { print \$s, qq(\n) }' -- -v={1} Report/species/clustalo-8species.clustal_num" > Report/species/clustalo-8species_str.tsv
echo -n -e "alignment\t" >> Report/species/clustalo-8species_str.tsv
perl -e 'while (<>) { if (/\s+\*/) { chomp; s/^\s{22}//; $s .= $_ } } END { print $s, "\n" }' Report/species/clustalo-8species.clustal_num >> Report/species/clustalo-8species_str.tsv
## add 12 whitespaces to the last row where every position is divergent

## castaneus vs domesticus
grep '>' Data/strains/castaneus_domesticus.chrM.fa | awk '{ print $1 }' | tr -d '>' | parallel -j 1 -k "perl -se 'BEGIN { print \$v, qq(\t) } while (<>) { if (/^\$v/) { @F=split /\s+/; \$s .= \$F[1] } } END { print \$s, qq(\n) }' -- -v={1} Report/strains/clustalo-castaneus_domesticus.clustal_num" > Report/strains/clustalo-castaneus_domesticus_str.tsv
echo -n -e "alignment\t" >> Report/strains/clustalo-castaneus_domesticus_str.tsv
perl -e 'while (<>) { if (/\s+\*/) { chomp; s/^\s{29}//; $s .= $_ } } END { print $s, "\n" }' Report/strains/clustalo-castaneus_domesticus.clustal_num >> Report/strains/clustalo-castaneus_domesticus_str.tsv
