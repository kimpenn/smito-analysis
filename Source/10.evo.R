###########################################################################
## Author: Youtao Lu <luyoutao@sas.upenn.edu>
##  
## Copyright (c) 2021-2023, Youtao Lu and Junhyong Kim, Department of Biology, University of Pennsylvania
## Copyright (c) 2021-2023, James Eberwine, Perelman School of Medicine, University of Pennsylvania
## Copyright (c) 2021-2023, David Issadore, Penn Engineering, University of Pennsylvania
## 
## This Source Code is distributed under Creative Commons Attribution License 4.0 (CC BY).
###########################################################################
source("Source/functions.R")
library("data.table")
library("ggplot2")

chrmbases_properties <- fread("Report/artifact/chrmbases_properties.csv.gz")
chrmgenes <- fread("Report/artifact/chrmgenes.csv")
chrmgenes <- rbind(chrmgenes,
    chrmbases_properties[is_in_range == "Y" & is_in_primer == "N" & SNVID == "SNV8", list(symbol = "D-loop", start = min(pos), end = max(pos), strand = "*", nbases_covered = .N)], 
    chrmbases_properties[pos %in% 3843:3844, list(symbol = "intergenic", start = min(pos), end = max(pos), strand = "*", nbases_covered = .N)]
)
chrmgenes_cds <- chrmgenes[c(6, 10, 16, 19, 21, 22, 23, 25, 27, 28, 32, 33, 35), ]
chrmbases_properties[, is_cds := ifelse(sapply(pos, function(p) { any(sapply(1:chrmgenes_cds[,.N], function(x) { b <- chrmgenes_cds[x]; p >= b[, start] & p <= b[, end] })) }), "Y", "N")]
chrmbases_properties[, proteins := sapply(pos, function(p) { x <- sapply(1:chrmgenes_cds[,.N], function(x) { b <- chrmgenes_cds[x]; ifelse(p >= b[, start] & p <= b[, end], b[, symbol], "") }); x <- x[x!=""]; paste0(x, collapse = ",") })]

chrmbases_assayed <- chrmbases_properties[is_in_range == "Y" & is_in_primer == "N"]


genus_mclust <- fread(file = "Report/genus/clustalo-mm10_rn6_base.csv")
names(genus_mclust)[1] <- "pos"
names(genus_mclust)[2] <- "ref"
setkey(genus_mclust, pos)
genus_mclust_varonly <- subset(genus_mclust, alignment != "*" & ref != "-")
genus_mclust_varonly <- genus_mclust_varonly[, -c("alignment")]
fwrite(genus_mclust_varonly, file = "Report/evo/genus_mclust_varonly.csv")
genus_mclust_varonly <- fread(file = "Report/evo/genus_mclust_varonly.csv")

species_mclust <- fread(file = "Report/species/clustalo-8species_base.csv")
names(species_mclust)[1] <- "pos"
names(species_mclust)[which(names(species_mclust) == "Mus_musculus")] <- "ref"
setkey(species_mclust, pos)
setcolorder(species_mclust, c(1, 9, 2:8, 10))
species_mclust_varonly <- subset(species_mclust, alignment != "*" & ref != "-")
species_mclust_varonly <- species_mclust_varonly[, -c("alignment")]
dim(species_mclust_varonly)
## [1] 4240    9
fwrite(species_mclust_varonly, file = "Report/evo/species_mclust_varonly.csv")
species_mclust_varonly <- fread(file = "Report/evo/species_mclust_varonly.csv")


strains_mclust <- fread(file = "Report/strains/clustalo-17strains_base.csv")
names(strains_mclust)[1] <- "pos"
names(strains_mclust)[which(names(strains_mclust) == "C57BL/6J")] <- "ref"
setkey(strains_mclust, pos)
setcolorder(strains_mclust, c(1, 18, 2:17, 19))
strains_mclust_varonly <- subset(strains_mclust, alignment != "*" & ref != "-")
strains_mclust_varonly <- strains_mclust_varonly[, -c("alignment")]
dim(strains_mclust_varonly)
## [1] 1494   18
fwrite(strains_mclust_varonly, file = "Report/evo/strains_mclust_varonly.csv")


castdom_mclust <- fread(file = "Report/population/clustalo-castaneus_domesticus_base.csv")
names(castdom_mclust)[1] <- "pos"
names(castdom_mclust)[which(names(castdom_mclust) == "C57BL/6J")] <- "ref"
setkey(castdom_mclust, pos)
setcolorder(castdom_mclust, c(1, 31, 2:30, 32))
castdom_mclust_varonly <- subset(castdom_mclust, alignment != "*" & ref != "-")
castdom_mclust_varonly <- castdom_mclust_varonly[, -c("alignment")]
dim(castdom_mclust_varonly)
## [1] 722  31

X_cast <- data.table(t(castdom_mclust_varonly[, castaneus_TW3:castaneus_KM][, apply(.SD, 1, function(x) table(factor(x, levels = c("A", "C", "G", "T", "-"))))]))
X_dom <- data.table(t(castdom_mclust_varonly[, domesticus_Ker_Gui05:domesticus_HB_4242][, apply(.SD, 1, function(x) table(factor(x, levels = c("A", "C", "G", "T", "-"))))]))
X_castdom <- data.table(X_cast, X_dom)
kcast <- X_castdom[, 1:5][, apply(.SD, 1, function(x) any(x == 9))]
kdom <- X_castdom[, 6:10][, apply(.SD, 1, function(x) any(x == 20))]
castdom_div_pos <- X_castdom[, castdom_mclust_varonly[kcast & kdom, pos]]
castdom_poly_pos <- X_castdom[, castdom_mclust_varonly[xor(kcast, kdom), pos]]

castdom_mclust_varonly[, type := ""]
castdom_mclust_varonly[pos %in% castdom_div_pos, type := "div"]
castdom_mclust_varonly[pos %in% castdom_poly_pos, type := "poly"]
fwrite(castdom_mclust_varonly, file = "Report/evo/population_mclust_varonly.csv")
castdom_mclust_varonly <- fread(file = "Report/evo/population_mclust_varonly.csv")

###########################################################################
## Ti/Tv in genus, species, strains and population data
###########################################################################
everymuts_vep_unique <- fread("Report/SNVs/impact/everymuts_vep_unique.csv.gz")
classes <- c("nonsynonymous", "synonymous", "tRNA", "rRNA", "D-loop", "intergenic")

genus_basediff <- copy(genus_mclust_varonly)
genus_basediff[, alt := toupper(Rattus_norvegicus)]
genus_basediff[, alt := sub('-', "del", alt)]
genus_basediff[, ref := toupper(ref)]
genus_basediff[, mut := paste0(ref, '>', alt)]
genus_basediff[, posmut := paste0(pos, ':', mut)]
genus_basediff <- genus_basediff[!duplicated(posmut)]
genus_basediff <- everymuts_vep_unique[, c("posmut", "class")][genus_basediff, on = "posmut"]
genus_transmuts_byclass <- sapply(classes, function(cl) dcast.data.table(genus_basediff[class == cl, .N, by = c("ref", "alt")], ref ~ alt, value.var = "N"), simplify = FALSE)
genus_titv_byclass <- t(sapply(genus_transmuts_byclass, function(X) Genetics$titv(Tools$dt2df(X))))
genus_transmuts <- dcast.data.table(genus_basediff[, .N, by = c("ref", "alt")], ref ~ alt, value.var = "N")
genus_titv_byclass <- rbind(genus_titv_byclass, Genetics$titv(Tools$dt2df(genus_transmuts)))
rownames(genus_titv_byclass) <- c(classes, "average")
write.csv(genus_titv_byclass, file = "Report/evo/genus_titv_byclass.csv")
genus_titv_byclass <- fread(file = "Report/evo/genus_titv_byclass.csv")
genus_basediff_assayed <- everymuts_vep_unique[, c("posmut", "class")][genus_basediff[pos %in% chrmbases_assayed[, pos]], on = "posmut"]
genus_transmuts_byclass_assayed <- sapply(classes[classes %in% unique(genus_basediff_assayed[, class])], function(cl) dcast.data.table(genus_basediff_assayed[class == cl, .N, by = c("ref", "alt")], ref ~ alt, value.var = "N"), simplify = FALSE)
genus_titv_byclass_assayed <- t(sapply(genus_transmuts_byclass_assayed, function(X) Genetics$titv(Tools$dt2df(X))))
genus_titv_byclass_assayed <- as.data.frame(genus_titv_byclass_assayed)[classes, ]
rownames(genus_titv_byclass_assayed) <- classes
genus_titv_byclass_assayed <- as.matrix(genus_titv_byclass_assayed)
genus_transmuts_assayed <- dcast.data.table(genus_basediff_assayed[, .N, by = c("ref", "alt")], ref ~ alt, value.var = "N")
genus_titv_byclass_assayed <- rbind(genus_titv_byclass_assayed, Genetics$titv(Tools$dt2df(genus_transmuts_assayed)))
rownames(genus_titv_byclass_assayed) <- c(classes, "average")
write.csv(genus_titv_byclass_assayed, file = "Report/evo/genus_titv_byclass_assayed.csv")
genus_titv_byclass_assayed <- fread(file = "Report/evo/genus_titv_byclass_assayed.csv")

species_basediff <- copy(species_mclust_varonly)
species_basediff <- species_basediff[, apply(.SD, 1, function(x) paste(x[1], ">", x[-1], sep = "")), by = "pos"]
names(species_basediff)[2] <- "mut"
species_basediff[, mut := toupper(mut)]
species_basediff[, mut := sub('-', "del", mut)]
species_basediff[, ref := sapply(strsplit(mut, '>'), '[', 1)]
species_basediff[, alt := sapply(strsplit(mut, '>'), '[', 2)]
species_basediff <- species_basediff[ref != alt]
species_basediff[, posmut := paste0(pos, ':', mut)]
species_basediff <- species_basediff[!duplicated(posmut)]
species_basediff <- everymuts_vep_unique[, c("posmut", "class")][species_basediff, on = "posmut"]
species_transmuts_byclass <- sapply(classes, function(cl) dcast.data.table(species_basediff[class == cl, .N, by = c("ref", "alt")], ref ~ alt, value.var = "N"), simplify = FALSE)
species_titv_byclass <- t(sapply(species_transmuts_byclass, function(X) Genetics$titv(Tools$dt2df(X))))
species_transmuts <- dcast.data.table(species_basediff[, .N, by = c("ref", "alt")], ref ~ alt, value.var = "N")
species_titv_byclass <- rbind(species_titv_byclass, Genetics$titv(Tools$dt2df(species_transmuts)))
rownames(species_titv_byclass) <- c(classes, "average")
write.csv(species_titv_byclass, file = "Report/evo/species_titv_byclass.csv")
species_titv_byclass <- fread(file = "Report/evo/species_titv_byclass.csv")
species_basediff_assayed <- everymuts_vep_unique[, c("posmut", "class")][species_basediff[pos %in% chrmbases_assayed[, pos]], on = "posmut"]
species_transmuts_byclass_assayed <- sapply(classes[classes %in% unique(species_basediff_assayed[, class])], function(cl) dcast.data.table(species_basediff_assayed[class == cl, .N, by = c("ref", "alt")], ref ~ alt, value.var = "N"), simplify = FALSE)
species_titv_byclass_assayed <- t(sapply(species_transmuts_byclass_assayed, function(X) Genetics$titv(Tools$dt2df(X))))
species_titv_byclass_assayed <- as.data.frame(species_titv_byclass_assayed)[classes, ]
rownames(species_titv_byclass_assayed) <- classes
species_titv_byclass_assayed <- as.matrix(species_titv_byclass_assayed)
species_transmuts_assayed <- dcast.data.table(species_basediff_assayed[, .N, by = c("ref", "alt")], ref ~ alt, value.var = "N")
species_titv_byclass_assayed <- rbind(species_titv_byclass_assayed, Genetics$titv(Tools$dt2df(species_transmuts_assayed)))
rownames(species_titv_byclass_assayed) <- c(classes, "average")
write.csv(species_titv_byclass_assayed, file = "Report/evo/species_titv_byclass_assayed.csv")
species_titv_byclass_assayed <- fread(file = "Report/evo/species_titv_byclass_assayed.csv")


strains_basediff <- copy(strains_mclust_varonly)
strains_basediff <- strains_basediff[, apply(.SD, 1, function(x) paste(x[1], ">", x[-1], sep = "")), by = "pos"]
names(strains_basediff)[2] <- "mut"
strains_basediff[, mut := toupper(mut)]
strains_basediff[, mut := sub('-', "del", mut)]
strains_basediff[, ref := sapply(strsplit(mut, '>'), '[', 1)]
strains_basediff[, alt := sapply(strsplit(mut, '>'), '[', 2)]
strains_basediff <- strains_basediff[ref != alt]
strains_basediff[, posmut := paste0(pos, ':', mut)]
strains_basediff <- strains_basediff[!duplicated(posmut)]
strains_basediff <- everymuts_vep_unique[, c("posmut", "class")][strains_basediff, on = "posmut"]
strains_transmuts_byclass <- sapply(classes[classes %in% unique(strains_basediff[, class])], function(cl) dcast.data.table(strains_basediff[class == cl, .N, by = c("ref", "alt")], ref ~ alt, value.var = "N"), simplify = FALSE)
strains_titv_byclass <- t(sapply(strains_transmuts_byclass, function(X) Genetics$titv(Tools$dt2df(X))))
strains_titv_byclass <- as.data.frame(strains_titv_byclass)[classes, ]
rownames(strains_titv_byclass) <- classes
strains_titv_byclass <- as.matrix(strains_titv_byclass)
strains_transmuts <- dcast.data.table(strains_basediff[, .N, by = c("ref", "alt")], ref ~ alt, value.var = "N")
strains_titv_byclass <- rbind(strains_titv_byclass, Genetics$titv(Tools$dt2df(strains_transmuts)))
rownames(strains_titv_byclass) <- c(classes, "average")
write.csv(strains_titv_byclass, file = "Report/evo/strains_titv_byclass.csv")
strains_titv_byclass <- fread(file = "Report/evo/strains_titv_byclass.csv")
strains_basediff_assayed <- everymuts_vep_unique[, c("posmut", "class")][strains_basediff[pos %in% chrmbases_assayed[, pos]], on = "posmut"]
strains_transmuts_byclass_assayed <- sapply(classes[classes %in% unique(strains_basediff_assayed[, class])], function(cl) dcast.data.table(strains_basediff_assayed[class == cl, .N, by = c("ref", "alt")], ref ~ alt, value.var = "N"), simplify = FALSE)
strains_titv_byclass_assayed <- t(sapply(strains_transmuts_byclass_assayed, function(X) Genetics$titv(Tools$dt2df(X))))
strains_titv_byclass_assayed <- as.data.frame(strains_titv_byclass_assayed)[classes, ]
rownames(strains_titv_byclass_assayed) <- classes
strains_titv_byclass_assayed <- as.matrix(strains_titv_byclass_assayed)
strains_transmuts_assayed <- dcast.data.table(strains_basediff_assayed[, .N, by = c("ref", "alt")], ref ~ alt, value.var = "N")
strains_titv_byclass_assayed <- rbind(strains_titv_byclass_assayed, Genetics$titv(Tools$dt2df(strains_transmuts_assayed)))
rownames(strains_titv_byclass_assayed) <- c(classes, "average")
write.csv(strains_titv_byclass_assayed, file = "Report/evo/strains_titv_byclass_assayed.csv")
strains_titv_byclass_assayed <- fread(file = "Report/evo/strains_titv_byclass_assayed.csv")


castdom_basediff <- copy(castdom_mclust_varonly)
castdom_basediff <- castdom_basediff[, apply(.SD, 1, function(x) paste(x[1], ">", x[-1], sep = "")), by = "pos"]
names(castdom_basediff)[2] <- "mut"
castdom_basediff[, mut := toupper(mut)]
castdom_basediff[, mut := sub('-', "del", mut)]
castdom_basediff[, ref := sapply(strsplit(mut, '>'), '[', 1)]
castdom_basediff[, alt := sapply(strsplit(mut, '>'), '[', 2)]
castdom_basediff <- castdom_basediff[ref != alt]
castdom_basediff[, posmut := paste0(pos, ':', mut)]
castdom_basediff <- castdom_basediff[!duplicated(posmut)]
castdom_basediff <- everymuts_vep_unique[, c("posmut", "class")][castdom_basediff, on = "posmut"]
castdom_transmuts_byclass <- sapply(classes, function(cl) dcast.data.table(castdom_basediff[class == cl, .N, by = c("ref", "alt")], ref ~ alt, value.var = "N"), simplify = FALSE)
castdom_titv_byclass <- t(sapply(castdom_transmuts_byclass, function(X) Genetics$titv(Tools$dt2df(X))))
castdom_transmuts <- dcast.data.table(castdom_basediff[, .N, by = c("ref", "alt")], ref ~ alt, value.var = "N")
castdom_titv_byclass <- rbind(castdom_titv_byclass, Genetics$titv(Tools$dt2df(castdom_transmuts)))
rownames(castdom_titv_byclass) <- c(classes, "average")
write.csv(castdom_titv_byclass, file = "Report/evo/population_titv_byclass.csv")
castdom_titv_byclass <- fread(file = "Report/evo/population_titv_byclass.csv")
castdom_basediff_assayed <- everymuts_vep_unique[, c("posmut", "class")][castdom_basediff[pos %in% chrmbases_assayed[, pos]], on = "posmut"]
castdom_transmuts_byclass_assayed <- sapply(classes[classes %in% unique(castdom_basediff_assayed[, class])], function(cl) dcast.data.table(castdom_basediff_assayed[class == cl, .N, by = c("ref", "alt")], ref ~ alt, value.var = "N"), simplify = FALSE)
castdom_titv_byclass_assayed <- t(sapply(castdom_transmuts_byclass_assayed, function(X) Genetics$titv(Tools$dt2df(X))))
castdom_titv_byclass_assayed <- as.data.frame(castdom_titv_byclass_assayed)[classes, ]
rownames(castdom_titv_byclass_assayed) <- classes
castdom_titv_byclass_assayed <- as.matrix(castdom_titv_byclass_assayed)
castdom_transmuts_assayed <- dcast.data.table(castdom_basediff_assayed[, .N, by = c("ref", "alt")], ref ~ alt, value.var = "N")
castdom_titv_byclass_assayed <- rbind(castdom_titv_byclass_assayed, Genetics$titv(Tools$dt2df(castdom_transmuts_assayed)))
rownames(castdom_titv_byclass_assayed) <- c(classes, "average")
write.csv(castdom_titv_byclass_assayed, file = "Report/evo/population_titv_byclass_assayed.csv")
castdom_titv_byclass_assayed <- fread(file = "Report/evo/population_titv_byclass_assayed.csv")

## whole genome
genus_titv_byclass_dt <- data.table(data = "genera", class = rownames(genus_titv_byclass), genus_titv_byclass)
species_titv_byclass_dt <- data.table(data = "species", class = rownames(species_titv_byclass), species_titv_byclass)
strains_titv_byclass_dt <- data.table(data = "strains", class = rownames(strains_titv_byclass), strains_titv_byclass)
populations_titv_byclass_dt <- data.table(data = "population", class = rownames(castdom_titv_byclass), castdom_titv_byclass)
evo_titv_byclass <- rbindlist(list(genus_titv_byclass_dt, species_titv_byclass_dt, strains_titv_byclass_dt, populations_titv_byclass_dt))
evo_titv_byclass[, data := factor(data, levels = c("genera", "species", "strains", "population"))]

pdf("Report/evo/evo_titv_byclass.pdf", width = 6, height = 6)
ggplot(evo_titv_byclass[class != "intergenic"], aes(x = class, fill = data, y = r)) + geom_bar(stat = "identity", position = position_dodge2(preserve = "single")) + scale_x_discrete(limits = c("rRNA", "tRNA", "synonymous", "nonsynonymous", "average", "D-loop"), labels = c(rRNA = "rRNA", tRNA = "tRNA", synonymous = "syn", nonsynonymous = "nonsyn", average = "average", "D-loop" = "D-loop")) + theme_classic(base_size = 16) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + xlab("") + ylab("Ti/Tv") + scale_fill_brewer(palette = "Set1") + scale_y_continuous()
dev.off()

## SMITO covered
genus_titv_byclass_assayed_dt <- data.table(data = "genera", class = rownames(genus_titv_byclass_assayed), genus_titv_byclass_assayed)
species_titv_byclass_assayed_dt <- data.table(data = "species", class = rownames(species_titv_byclass_assayed), species_titv_byclass_assayed)
strains_titv_byclass_assayed_dt <- data.table(data = "strains", class = rownames(strains_titv_byclass_assayed), strains_titv_byclass_assayed)
populations_titv_byclass_assayed_dt <- data.table(data = "population", class = rownames(castdom_titv_byclass_assayed), castdom_titv_byclass_assayed)
evo_titv_byclass_assayed <- rbindlist(list(genus_titv_byclass_assayed_dt, species_titv_byclass_assayed_dt, strains_titv_byclass_assayed_dt, populations_titv_byclass_assayed_dt))
evo_titv_byclass_assayed[, data := factor(data, levels = c("genera", "species", "strains", "population", "SMITO"))]

pdf("Report/evo/evo_titv_byclass_assayed.pdf", width = 6.5, height = 6)
ggplot(evo_titv_byclass_assayed[class != "intergenic"], aes(x = class, fill = data, y = r)) + geom_bar(stat = "identity", position = position_dodge2(preserve = "single")) + scale_x_discrete(limits = c("rRNA", "tRNA", "synonymous", "nonsynonymous", "average", "D-loop"), labels = c(rRNA = "rRNA", tRNA = "tRNA", synonymous = "syn", nonsynonymous = "nonsyn", average = "average", "D-loop" = "D-loop")) + theme_classic(base_size = 16) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + xlab("") + ylab("Ti/Tv") + scale_fill_brewer(palette = "Set1") + scale_y_continuous()
dev.off()

###########################################################################
## Check the codon table
###########################################################################
everymuts_vep <- fread(file = "Report/SNVs/impact/everymuts_vep.csv.gz")
everymuts_vep_synonymous <- everymuts_vep[Feature_type == "Transcript" & Consequence == "synonymous_variant"]
everymuts_vep_nonsynonymous <- everymuts_vep[Feature_type == "Transcript" & Consequence %in% 
    c("frameshift_variant", 
      "frameshift_variant,start_lost", 
      "frameshift_variant,stop_lost", 
      "incomplete_terminal_codon_variant,coding_sequence_variant", 
      "missense_variant", 
      "protein_altering_variant,incomplete_terminal_codon_variant", 
      "start_lost", 
      "stop_gained", 
      "stop_gained,start_lost", 
      "stop_lost", 
      "stop_retained_variant")]
everymuts_vep_synonymous[, class := "synonymous"]
everymuts_vep_nonsynonymous[, class := "nonsynonymous"]
everymuts_vep_cds <- rbind(everymuts_vep_synonymous, everymuts_vep_nonsynonymous)
everymuts_vep_codons <- everymuts_vep_cds[, list(symbol = SYMBOL, codon = sapply(strsplit(Codons, "/"), function(x) toupper(x[1]))), by = c("pos", "STRAND", "CDS_position")]
everymuts_vep_codons[, table(codon)]
## codon
##  AAA  AAC  AAG  AAT  ACA  ACC  ACG  ACT  AGC  AGT  ATA  ATC  ATG  ATT  CAA  CAC
## 1212 1296   24  720 1884 1032   72  708  420  168 2604 1668  360 2808  936  756
##  CAG  CAT  CCA  CCC  CCG  CCT  CGA  CGC  CGG  CGT  CTA  CTC  CTG  CTT  GAA  GAC
##   36  408 1560  408   24  360  420  216   36   96 3192  768  324 1044  960  516
##  GAG  GAT  GCA  GCC  GCG  GCT  GGA  GGC  GGG  GGT  GTA  GTC  GTG  GTT    T  TAA
##  132  372 1152  996   84  564 1308  468  348  420  888  396  132  600   12   96
##  TAC  TAG  TAT  TCA  TCC  TCG  TCT  TGA  TGC  TGG  TGT  TTA  TTC  TTG  TTT
##  720   24  780 1764  564   48  504 1152  216   84  120 1548 1596  180 1308
everymuts_vep_codons[codon == "T"]
##       pos STRAND CDS_position  symbol codon
##  1:  9390      1          784  mt-Co3     T
##  2:  9390      1          784  mt-Co3     T
##  3:  9390      1          784  mt-Co3     T
##  4:  9390      1          784  mt-Co3     T
##  5: 11544      1         1378  mt-Nd4     T
##  6: 11544      1         1378  mt-Nd4     T
##  7: 11544      1         1378  mt-Nd4     T
##  8: 11544      1         1378  mt-Nd4     T
##  9: 15288      1         1144 mt-Cytb     T
## 10: 15288      1         1144 mt-Cytb     T
## 11: 15288      1         1144 mt-Cytb     T
## 12: 15288      1         1144 mt-Cytb     T

## mt-Co3, mt-Nd4, mt-Cytb have the stop codon TAA by polyadenylation. 
everymuts_vep_codons[codon == "T", codon:="TAA"]
everymuts_vep_codons <- unique(everymuts_vep_codons)
everymuts_vep_codons[, strand := ifelse(STRAND == 1, "+", ifelse(STRAND == -1, "-", ""))]
everymuts_vep_codons[, STRAND := NULL]

###########################################################################
## MK test on the population data
###########################################################################
vertmtcodon <- fread(file = "Report/evo/vertmtcodon.csv")

castdom_mclust_varonly <- vertmtcodon[castdom_mclust_varonly, on = "pos"]
setcolorder(castdom_mclust_varonly, c(8:42, 1:7))

castdom_mclust_varonly[, CDS_position:= as.integer(CDS_position)]
castdom_mclust_varonly[, codon_index := round((CDS_position+1)/3)]

setcolorder(castdom_mclust_varonly, c(1:2, 43, 3:42))
fwrite(castdom_mclust_varonly, file = "Report/evo/population_mclust_varonly.csv")
castdom_mclust_varonly <- fread(file = "Report/evo/population_mclust_varonly.csv")

castdom_mclust_varonly_ns_percodon <- castdom_mclust_varonly[!is.na(symbol) & type != ""][, list(sites = paste0(unique(pos), collapse = ","), types = paste0(unique(type), collapse = ","), AA = unique(AA), syn = mean(S), nonsyn = mean(N)), by = c("symbol", "strand", "codon_index")]
fwrite(castdom_mclust_varonly_ns_percodon, file = "Report/evo/population_mclust_varonly_ns_percodon.csv")
castdom_mclust_varonly_ns_percodon <- fread(file = "Report/evo/population_mclust_varonly_ns_percodon.csv")

castdom_mclust_varonly_ns_pergene <- castdom_mclust_varonly_ns_percodon[types %in% c("div", "poly")][, list(S = sum(syn), N = sum(nonsyn)), keyby = c("symbol", "types")]
fwrite(castdom_mclust_varonly_ns_pergene, file = "Report/evo/population_mclust_varonly_ns_pergene.csv")
castdom_mclust_varonly_ns_pergene <- fread(file = "Report/evo/population_mclust_varonly_ns_pergene.csv")

K <- data.table(expand.grid(symbol = castdom_mclust_varonly_ns_pergene[, unique(symbol)], types = c("div", "poly"), stringsAsFactors = FALSE))[order(symbol, types)]
castdom_mclust_varonly_ns_pergene <- castdom_mclust_varonly_ns_pergene[K, on = c("symbol", "types")]
castdom_mclust_varonly_ns_pergene[is.na(S), S:=0]
castdom_mclust_varonly_ns_pergene[is.na(N), N:=0]
castdom_mclust_varonly_mkpval_pergene <- castdom_mclust_varonly_ns_pergene[, list(p = chisq.test(matrix(unlist(.SD), ncol = 2))$p.value), by = c("symbol"), .SDcols = c("S", "N")]
fwrite(castdom_mclust_varonly_mkpval_pergene, file = "Report/evo/population_mclust_varonly_mkpval_pergene.csv")
castdom_mclust_varonly_mkpval_pergene <- fread(file = "Report/evo/population_mclust_varonly_mkpval_pergene.csv")
