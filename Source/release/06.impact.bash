###########################################################################
## Author: Youtao Lu <luyoutao@sas.upenn.edu>
##  
## Copyright (c) 2021-2023, Youtao Lu and Junhyong Kim, Department of Biology, University of Pennsylvania
## Copyright (c) 2021-2023, James Eberwine, Perelman School of Medicine, University of Pennsylvania
## Copyright (c) 2021-2023, David Issadore, Penn Engineering, University of Pennsylvania
## 
## This Source Code is distributed under Creative Commons Attribution License 4.0 (CC BY).
###########################################################################
## all SNVs
vep --offline --cache --fasta Data/mm10.mito/chrM.fa --species mus_musculus --force_overwrite --everything -i Report/release/SNVs/impact/noctrl_vars.txt -o Report/release/SNVs/impact/noctrl_vep.txt
