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
library("ggsignif")

snv_info <- fread("Data/snv_loci_v2.csv")
snvIDs <- snv_info[, SNVID]
mito_barcodes <- fread("Data/mito_barcodes.csv")
mitoIDs <- mito_barcodes[, ID] 

chrmbases_properties <- fread("Report/artifact/chrmbases_properties.csv.gz")
chrmbases <- chrmbases_properties[, ref]
nchrmbases <- length(chrmbases)
nchrmbases
## [1] 16299

MitoInfo <- fread("Report/metadata/MitoInfo.csv")
dim(MitoInfo)
## [1] 1717   19
MitoInfo[, ExptID := factor(ExptID)]
MitoInfo[, MitoID := factor(MitoID, levels = mitoIDs)]
MitoInfo[, CellID := factor(CellID)]

CellInfo <- fread("Report/metadata/CellInfo.csv")
dim(CellInfo)
## [1] 102  12
MouseInfo <- fread("Report/metadata/MouseInfo.csv")
dim(MouseInfo)
## [1] 13  2

###########################################################################
## Prepare the input for VEP
###########################################################################
support_byposmut <- fread("Report/SNVs/filter/basediffperc_cutdemux_sub500k_q30_unstranded_highdepth_highaf_qcfltd_support_byposmut.csv")
posmuts <- support_byposmut[, posmut]
noctrl_vars_start <- support_byposmut[, pos]
noctrl_vars_end <- noctrl_vars_start

noctrl_vars_vep <- data.table(
    chr = rep("MT", nrow(support_byposmut)),
    start = noctrl_vars_start, 
    end = noctrl_vars_end, 
    allele = support_byposmut[, paste0(ref, '/', ifelse(alt == "del", '-', alt))], 
    strand = rep("+", nrow(support_byposmut)), 
    id = support_byposmut[, posmut]
)
fwrite(noctrl_vars_vep, file = "Report/SNVs/impact/noctrl_vars.txt", col.names = FALSE, sep = "\t")

###########################################################################
## convert the VEP output into data table
###########################################################################
vep <- fread("Report/SNVs/impact/noctrl_vep.txt", skip = 85, header = TRUE, sep = "\t")
setnames(vep, "#Uploaded_variation", "Uploaded_variation")
vep <- vep[Uploaded_variation %in% posmuts]
notes <- vep[, list(lapply(strsplit(Extra, ";"), function(x) { y <- strsplit(x, "="); n <- sapply(y, '[', 1); m <- sapply(y, '[', 2); structure(m, names = n) }))][, V1]
cnames <- Reduce(union, lapply(notes, names))
values <- t(sapply(notes, function(kv) kv[cnames]))
colnames(values) <- cnames
vep[, Extra := NULL]
vep <- data.table(vep, values)
vep <- support_byposmut[vep, on = c("posmut" = "Uploaded_variation")]
vep[Consequence == "start_lost", Location:IMPACT]
##    Location Allele               Gene            Feature Feature_type
## 1:  MT:9459      G ENSMUSG00000064360 ENSMUST00000082411   Transcript
## 2:  MT:9461      C ENSMUSG00000064360 ENSMUST00000082411   Transcript
##    Consequence cDNA_position CDS_position Protein_position Amino_acids  Codons
## 1:  start_lost             1            1                1         M/V Att/Gtt
## 2:  start_lost             3            3                1         M/I atT/atC
##    Existing_variation IMPACT
## 1:                  -   HIGH
## 2:                  -   HIGH

## Note, VEP currently missed the alternative initial codon for mouse mito, see
## https://github.com/Ensembl/ensembl-vep/issues/1292
## As a result, 9461 should be synonymous SNV in fact. VEP mistakenly annotated it
## as nonsyn and high impact. Let's manually correct it. 
vep[Consequence == "start_lost" & posmut == "9461:T>C", c("Consequence", "IMPACT") := list("synonymous_variant", "LOW")]
vep[pos == 9461 & Feature_type == "Transcript" & Consequence == "synonymous_variant", Location:IMPACT]
##    Location Allele               Gene            Feature Feature_type
## 1:  MT:9461      A ENSMUSG00000064360 ENSMUST00000082411   Transcript
## 2:  MT:9461      C ENSMUSG00000064360 ENSMUST00000082411   Transcript
##           Consequence cDNA_position CDS_position Protein_position Amino_acids
## 1: synonymous_variant             3            3                1           M
## 2: synonymous_variant             3            3                1         M/I
##     Codons Existing_variation IMPACT
## 1: atT/atA                  -    LOW
## 2: atT/atC                  -    LOW
fwrite(vep, file = "Report/SNVs/impact/noctrl_vep.csv")
vep <- fread(file = "Report/SNVs/impact/noctrl_vep.csv")

###########################################################################
## annotate all SNVs
###########################################################################
vep <- fread(file = "Report/SNVs/impact/noctrl_vep.csv")
vep_synonymous <- vep[Feature_type == "Transcript" & Consequence == "synonymous_variant"]
vep_synonymous[, list(uniqueN(posmut), .N)]
##     V1   N
## 1: 189 189
vep_nonsynonymous <- vep[Feature_type == "Transcript" & Consequence %in% c("frameshift_variant", "incomplete_terminal_codon_variant,coding_sequence_variant", "missense_variant", "start_lost", "stop_gained")]
vep_nonsynonymous[, list(uniqueN(posmut), .N)]
##     V1   N
## 1: 440 440
vep_tRNA <- vep[Feature_type == "Transcript" & Consequence == "non_coding_transcript_exon_variant" & BIOTYPE == "Mt_tRNA"]
vep_tRNA[, list(uniqueN(posmut), .N)]
##     V1   N
## 1: 161 162
vep_tRNA[posmut == "3772:C>T"]
##      posmut SNVID  pos mut ref alt nmice ncells nmitos Location Allele
## 1: 3772:C>T  SNV5 3772 C>T   C   T     1      1      1  MT:3772      T
## 2: 3772:C>T  SNV5 3772 C>T   C   T     1      1      1  MT:3772      T
##                  Gene            Feature Feature_type
## 1: ENSMUSG00000064342 ENSMUST00000082393   Transcript
## 2: ENSMUSG00000064343 ENSMUST00000082394   Transcript
##                           Consequence cDNA_position CDS_position
## 1: non_coding_transcript_exon_variant            67            -
## 2: non_coding_transcript_exon_variant            71            -
##    Protein_position Amino_acids Codons Existing_variation   IMPACT DISTANCE
## 1:                -           -      -                  - MODIFIER       NA
## 2:                -           -      -                  - MODIFIER       NA
##    STRAND VARIANT_CLASS SYMBOL SYMBOL_SOURCE BIOTYPE CANONICAL GENE_PHENO EXON
## 1:      1           SNV  mt-Ti           MGI Mt_tRNA       YES         NA  1/1
## 2:     -1           SNV  mt-Tq           MGI Mt_tRNA       YES         NA  1/1
##                           HGVSc APPRIS ENSP SWISSPROT TREMBL UNIPARC SIFT
## 1: ENSMUST00000082393.1:n.67C>T
## 2: ENSMUST00000082394.1:n.71G>A
##    DOMAINS HGVSp HGVS_OFFSET class
## 1:                        NA  tRNA
## 2:                        NA  tRNA
## Because they have almost equal impact, we decide to keep Ti only.
vep_tRNA <- vep_tRNA[!duplicated(posmut)]
vep_rRNA <- vep[Feature_type == "Transcript" & Consequence == "non_coding_transcript_exon_variant" & BIOTYPE == "Mt_rRNA"]
##     V1   N
## 1: 147 147
vep_dloop <- vep[Feature_type == "RegulatoryFeature" & Consequence == "regulatory_region_variant" & pos > 15422][!posmut %in% Reduce(union, list(vep_nonsynonymous[, posmut], vep_synonymous[, posmut], vep_tRNA[, posmut], vep_rRNA[, posmut]))]
vep_intergenic <- vep[Feature_type == "RegulatoryFeature" & Consequence == "regulatory_region_variant" & BIOTYPE == "promoter_flanking_region"][!posmut %in% Reduce(union, list(vep_nonsynonymous[, posmut], vep_synonymous[, posmut], vep_tRNA[, posmut], vep_rRNA[, posmut], vep_dloop[, posmut]))]
vep_intergenic
##      posmut SNVID  pos mut ref alt nmice ncells nmitos Location Allele Gene
## 1: 3843:A>G  SNV5 3843 A>G   A   G     2      2      2  MT:3843      G    -
## 2: 3844:T>C  SNV5 3844 T>C   T   C     1      1      1  MT:3844      C    -
##               Feature      Feature_type               Consequence cDNA_position
## 1: ENSMUSR00000758872 RegulatoryFeature regulatory_region_variant             -
## 2: ENSMUSR00000758872 RegulatoryFeature regulatory_region_variant             -
##    CDS_position Protein_position Amino_acids Codons Existing_variation   IMPACT
## 1:            -                -           -      -                  - MODIFIER
## 2:            -                -           -      -                  - MODIFIER
##    DISTANCE STRAND VARIANT_CLASS SYMBOL SYMBOL_SOURCE                  BIOTYPE
## 1:       NA     NA           SNV                      promoter_flanking_region
## 2:       NA     NA           SNV                      promoter_flanking_region
##    CANONICAL GENE_PHENO EXON HGVSc APPRIS ENSP SWISSPROT TREMBL UNIPARC SIFT
## 1:                   NA
## 2:                   NA
##    DOMAINS HGVSp HGVS_OFFSET      class
## 1:                        NA intergenic
## 2:                        NA intergenic

vep_nonsynonymous[, class := "nonsynonymous"]
vep_synonymous[, class := "synonymous"]
vep_tRNA[, class := "tRNA"]
vep_rRNA[, class := "rRNA"]
vep_dloop[, class := "D-loop"]
vep_intergenic[, class := "intergenic"]
vep_nonsynonymous[, uniqueN(posmut)] + vep_synonymous[, uniqueN(posmut)] + vep_tRNA[, uniqueN(posmut)] + vep_rRNA[, uniqueN(posmut)] + vep_dloop[, uniqueN(posmut)] + vep_intergenic[, uniqueN(posmut)]
## [1] 1032
vep_unique <- rbindlist(list(vep_nonsynonymous, vep_synonymous, vep_tRNA, vep_rRNA, vep_dloop, vep_intergenic))
sum(vep_unique[, table(class)])
## [1] 1032
dim(vep_unique)
## [1] 1032  42
fwrite(vep_unique, file = "Report/SNVs/impact/noctrl_vep_unique.csv")
vep <- fread(file = "Report/SNVs/impact/noctrl_vep_unique.csv")

###########################################################################
## SNV occurrence by predicted functional consequence
###########################################################################
noctrl_sift <- fread("Report/20210510/SNVs/impact/noctrl_impact_added_SIFT.tsv")
pdf("Report/20210510/SNVs/impact/noctrl_impact_added_SIFT_boxplot.pdf", width = 4.5, height = 4.5)
ggplot(noctrl_sift[SIFT_confidence == "high"], aes(x = SIFT_prediction, y = `# mice`)) + geom_boxplot() + theme_classic(base_size = 12) + geom_signif(comparisons = list(c("tolerated", "deleterious"))) + xlab("High-confidence SIFT prediction") + ylab("Number of mice shared")
ggplot(noctrl_sift[SIFT_confidence == "high"], aes(x = SIFT_prediction, y = `# cells`)) + geom_boxplot() + theme_classic(base_size = 12) + geom_signif(comparisons = list(c("tolerated", "deleterious"))) + xlab("High-confidence SIFT prediction") + ylab("Number of cells shared")
ggplot(noctrl_sift[SIFT_confidence == "high"], aes(x = SIFT_prediction, y = `# mitos`)) + geom_boxplot() + theme_classic(base_size = 12) + geom_signif(comparisons = list(c("tolerated", "deleterious"))) + xlab("High-confidence SIFT prediction") + ylab("Number of mitos shared")
dev.off()

summary(glm(`# mitos` ~ SIFT_prediction, data = noctrl_sift[SIFT_confidence == "high"], family = poisson(link = "log")))
## 
## Call:
## glm(formula = `# mitos` ~ SIFT_prediction, family = poisson(link = "log"),
##     data = noctrl_sift[SIFT_confidence == "high"])
## 
## Coefficients:
##                          Estimate Std. Error z value Pr(>|z|)
## (Intercept)               1.23353    0.09853  12.519   <2e-16 ***
## SIFT_predictiontolerated  0.15276    0.12087   1.264    0.206
## ---
## Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
## 
## (Dispersion parameter for poisson family taken to be 1)
## 
##     Null deviance: 497.74  on 80  degrees of freedom
## Residual deviance: 496.11  on 79  degrees of freedom
## AIC: 715.3
## 
## Number of Fisher Scoring iterations: 6

summary(glm(`# cells` ~ SIFT_prediction, data = noctrl_sift[SIFT_confidence == "high"], family = poisson(link = "log")))
## 
## Call:
## glm(formula = `# cells` ~ SIFT_prediction, family = poisson(link = "log"),
##     data = noctrl_sift[SIFT_confidence == "high"])
## 
## Coefficients:
##                          Estimate Std. Error z value Pr(>|z|)
## (Intercept)               1.06471    0.10721   9.931   <2e-16 ***
## SIFT_predictiontolerated  0.07864    0.13320   0.590    0.555
## ---
## Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
## 
## (Dispersion parameter for poisson family taken to be 1)
## 
##     Null deviance: 301.75  on 80  degrees of freedom
## Residual deviance: 301.40  on 79  degrees of freedom
## AIC: 514.81
## 
## Number of Fisher Scoring iterations: 5

###########################################################################
## Among all SNVs, # of SNV sites per gene locus per number-of-mitos support
###########################################################################
vep <- fread(file = "Report/SNVs/impact/noctrl_vep_unique.csv")
vep[, gene := ifelse(SYMBOL == "", class, SYMBOL)]
pdf("Report/SNVs/impact/noctrl_vep_nsnvs_bygene_bynmitos.pdf", width = 6.5, height = 6)
ggplot(vep[, .N, by = c("gene", "nmitos")][order(nmitos)], aes(x = gene, y = N, fill = nmitos)) + geom_bar(stat = "identity") + scale_x_discrete(limits = vep[, .N, by = "gene"][order(-N), gene]) + theme_classic(base_size = 16) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), axis.text.y = element_text(angle = 90, hjust = 0.5, vjust = 1)) + xlab("") + ylab("# SNVs") + scale_fill_viridis_c(breaks = c(1, 10, 100, 1000), labels = c(1, 10, 100, 1000), trans = scales::log_trans()) + scale_y_continuous(breaks = seq(0, 150, by = 50))
dev.off()
pdf("Report/SNVs/impact/noctrl_vep_nsnvs_bygene_byncells.pdf", width = 6, height = 6)
ggplot(vep[, .N, by = c("gene", "ncells")][order(ncells)], aes(x = gene, y = N, fill = ncells)) + geom_bar(stat = "identity") + scale_x_discrete(limits = vep[, .N, by = "gene"][order(-N), gene]) + theme_classic(base_size = 16) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), axis.text.y = element_text(angle = 90, hjust = 0.5, vjust = 1)) + xlab("") + ylab("# SNVs") + scale_fill_viridis_c(breaks = c(1, 10, 100, 1000), labels = c(1, 10, 100, 1000), trans = scales::log_trans())
dev.off()
pdf("Report/SNVs/impact/noctrl_vep_nsnvs_bygene_bynmice.pdf", width = 6, height = 6)
ggplot(vep[, .N, by = c("gene", "nmice")][order(nmice)], aes(x = gene, y = N, fill = nmice)) + geom_bar(stat = "identity") + scale_x_discrete(limits = vep[, .N, by = "gene"][order(-N), gene]) + theme_classic(base_size = 16) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), axis.text.y = element_text(angle = 90, hjust = 0.5, vjust = 1)) + xlab("") + ylab("# SNVs") + scale_fill_viridis_c(breaks = seq(1, 13, by = 2), labels = seq(1, 13, by = 2))
dev.off()

###########################################################################
## Among all SNVs, how is the per-base SNV occurrence per region?
###########################################################################
vep <- fread(file = "Report/SNVs/impact/noctrl_vep_unique.csv")
vep[, gene := ifelse(SYMBOL == "", class, SYMBOL)]
chrmgenes <- fread("Report/artifact/chrmgenes.csv")
## get D-loop and intergenic info
chrmgenes <- rbind(chrmgenes,
    chrmbases_properties[is_in_range == "Y" & is_in_primer == "N" & SNVID == "SNV8", list(symbol = "D-loop", start = min(pos), end = max(pos), strand = "*", nbases_covered = .N)], 
    chrmbases_properties[pos %in% 3843:3844, list(symbol = "intergenic", start = min(pos), end = max(pos), strand = "*", nbases_covered = .N)]
)
npos_bygene <- vep[, list(npos = uniqueN(pos)), by = c("gene")]
npos_bygene <- chrmgenes[, c("symbol", "nbases_covered")][npos_bygene, on = c(symbol = "gene")]
npos_bygene[, npos_perbase := npos / nbases_covered]
setnames(npos_bygene, "symbol", "gene")
fwrite(npos_bygene[order(-npos_perbase)], file = "Report/SNVs/impact/noctrl_vep_npos_bygene_perbase.csv")
npos_bygene <- fread(file = "Report/SNVs/impact/noctrl_vep_npos_bygene_perbase.csv")

pdf("Report/SNVs/impact/noctrl_vep_npos_bygene_perbase.pdf", width = 5.4, height = 6)
ggplot(npos_bygene, aes(x = gene, y = npos_perbase)) + geom_bar(stat = "identity", fill = "gray") + scale_x_discrete(limits = npos_bygene[order(-npos_perbase), gene]) + theme_classic(base_size = 16) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), axis.text.y = element_text(angle = 90, hjust = 0.5, vjust = 1)) + xlab("") + ylab("# SNV sites per-base")
dev.off()

###########################################################################
## Transition probability matrix for all SNVs
###########################################################################
support_byposmut <- fread("Report/SNVs/filter/basediffperc_cutdemux_sub500k_q30_unstranded_highdepth_highaf_qcfltd_support_byposmut.csv")
transmat <- dcast.data.table(support_byposmut[, .N, by = c("ref", "alt")], ref ~ alt, value.var = "N")
transmat
##    ref  A   C   G   T del
## 1:   A NA  39 275  28  11
## 2:   C 90  NA  18 173   5
## 3:   G 89  11  NA  41  NA
## 4:   T 24 204  21  NA   3
## Because of the deletions, the sum of ti and tv will be less than 1032. 

vep <- fread(file = "Report/SNVs/impact/noctrl_vep_unique.csv")
vep_nonsynonymous <- vep[class == "nonsynonymous"]
vep_synonymous <- vep[class  == "synonymous"]
vep_tRNA <- vep[class == "tRNA"]
vep_rRNA <- vep[class == "rRNA"]
vep_dloop <- vep[class == "D-loop"]
vep_intergenic <- vep[class == "intergenic"]

vep_synonymous_posmut <- vep_synonymous[, unique(posmut)]
vep_nonsynonymous_posmut <- vep_nonsynonymous[, unique(posmut)]
vep_tRNA_posmut <- vep_tRNA[, unique(posmut)]
vep_rRNA_posmut <- vep_rRNA[, unique(posmut)]
vep_dloop_posmut <- vep_dloop[, unique(posmut)]
vep_intergenic_posmut <- vep_intergenic[, unique(posmut)]

transmat_nonsynonymous <- dcast.data.table(support_byposmut[posmut %in% vep_nonsynonymous_posmut, .N, by = c("ref", "alt")], ref ~ alt, value.var = "N")
transmat_synonymous <- dcast.data.table(support_byposmut[posmut %in% vep_synonymous_posmut, .N, by = c("ref", "alt")], ref ~ alt, value.var = "N")
transmat_tRNA <- dcast.data.table(support_byposmut[posmut %in% vep_tRNA_posmut, .N, by = c("ref", "alt")], ref ~ alt, value.var = "N")
transmat_rRNA <- dcast.data.table(support_byposmut[posmut %in% vep_rRNA_posmut, .N, by = c("ref", "alt")], ref ~ alt, value.var = "N")
transmat_dloop <- dcast.data.table(support_byposmut[posmut %in% vep_dloop_posmut, .N, by = c("ref", "alt")], ref ~ alt, value.var = "N")
transmat_intergenic <- dcast.data.table(support_byposmut[posmut %in% vep_intergenic_posmut, .N, by = c("ref", "alt")], ref ~ alt, value.var = "N")

Tools$write_xlsx(list(
    all = transmat,
    nonsynonymous = transmat_nonsynonymous, 
    synonymous = transmat_synonymous, 
    tRNA = transmat_tRNA, 
    rRNA = transmat_rRNA, 
    dloop= transmat_dloop, 
    intergenic = transmat_intergenic
), file = "Report/SNVs/impact/noctrl_transmat.xlsx", row.names = FALSE)

support_byposmut[, mut := factor(mut, levels = setdiff(paste0(rep(c("A", "C", "G", "T"), each = 5), ">", c("A", "C", "G", "T", "del")), c("A>A", "C>C", "G>G", "T>T")))]
pdf(file = "Report/SNVs/impact/noctrl_transmat_barplot.pdf", width = 5, height = 6)
par(ps = 12, lend = 2, ljoin = 1, bty = "L", mfrow = c(7, 1), mar = c(2, 2.5, 1, 0), oma = c(0, 0, 0, 0), mgp = c(1.5, 0.5, 0), cex.axis = 1)
support_byposmut[, barplot(table(mut), ylab = "# SNVs", main = "total", cex.names = 0.8, border = NA)]
support_byposmut[posmut %in% vep_nonsynonymous_posmut, barplot(table(mut), ylab = "# SNVs", main = "nonsynonymous", cex.names = 0.8, border = NA)]
support_byposmut[posmut %in% vep_synonymous_posmut, barplot(table(mut), ylab = "# SNVs", main = "synonymous", cex.names = 0.8, border = NA)]
support_byposmut[posmut %in% vep_tRNA_posmut, barplot(table(mut), ylab = "# SNVs", main = "tRNA", cex.names = 0.8, border = NA)]
support_byposmut[posmut %in% vep_rRNA_posmut, barplot(table(mut), ylab = "# SNVs", main = "rRNA", cex.names = 0.8, border = NA)]
support_byposmut[posmut %in% vep_dloop_posmut, barplot(table(mut), ylab = "# SNVs", main = "D-loop", cex.names = 0.8, border = NA)]
support_byposmut[posmut %in% vep_intergenic_posmut, barplot(table(mut), ylab = "# SNVs", main = "intergenic", cex.names = 0.8, border = NA)]
dev.off()

colmap <- structure(c(RColorBrewer::brewer.pal(5, "Set1"), "gray"), names = c("D-loop", "rRNA", "tRNA", "synonymous", "nonsynonymous", "intergenic"))[c(6, 1:5)]
nsnvs_byclass <- melt.data.table(dcast.data.table(vep, mut ~ class, fun.aggregate = function(x) uniqueN(x), value.var = "pos"), id.vars = "mut", variable.name = "class", value.name = "nsnvs")
nsnvs_byclass[, class := factor(class, levels = c("intergenic", "D-loop", "rRNA", "tRNA", "synonymous", "nonsynonymous"))]
fig <- ggplot(nsnvs_byclass, aes(x = mut, y = nsnvs, fill = class)) + geom_bar(stat = "identity") + theme_classic(12) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), axis.text.y = element_text(angle = 90, hjust = 0.5, vjust = 1)) + xlab("") + ylab("# SNVs") + scale_fill_manual(values = colmap) + scale_y_continuous(breaks = seq(0, 300, by = 50))
ggsave(fig, file = "Report/SNVs/impact/noctrl_nsnvs_byclass_bytrans_bar.pdf", width = 6, height = 4.5)

noctrl_titv <- list(
    all = Genetics$titv(Tools$dt2df(transmat)), 
    nonsynonymous = Genetics$titv(Tools$dt2df(transmat_nonsynonymous)),
    synonymous = Genetics$titv(Tools$dt2df(transmat_synonymous)),
    tRNA = Genetics$titv(Tools$dt2df(transmat_tRNA)),
    rRNA = Genetics$titv(Tools$dt2df(transmat_rRNA)),
    dloop = Genetics$titv(Tools$dt2df(transmat_dloop)), 
    intergenic = Genetics$titv(Tools$dt2df(transmat_intergenic))
)
noctrl_titv <- do.call(rbind, noctrl_titv)
noctrl_titv <- data.table(class = rownames(noctrl_titv), noctrl_titv)

fwrite(noctrl_titv, file = "Report/SNVs/impact/noctrl_titv.csv")
noctrl_titv <- fread(file = "Report/SNVs/impact/noctrl_titv.csv")

###########################################################################
## Among inherited SNVs, how is the per-base mutation occurrence per region?
###########################################################################
vep <- fread(file = "Report/SNVs/impact/noctrl_vep_unique.csv")
vep[, gene := ifelse(SYMBOL == "", class, SYMBOL)]
inherited_noctrl_vep <- vep[nmice >= 3]
chrmgenes <- fread("Report/artifact/chrmgenes.csv")
## get D-loop and intergenic info
chrmgenes <- rbind(chrmgenes,
    chrmbases_properties[is_in_range == "Y" & is_in_primer == "N" & SNVID == "SNV8", list(symbol = "D-loop", start = min(pos), end = max(pos), strand = "*", nbases_covered = .N)], 
    chrmbases_properties[pos %in% 3843:3844, list(symbol = "intergenic", start = min(pos), end = max(pos), strand = "*", nbases_covered = .N)]
)
inherited_noctrl_npos_bygene <- inherited_noctrl_vep[, list(npos = uniqueN(pos)), by = c("gene")]
inherited_noctrl_npos_bygene <- chrmgenes[, c("symbol", "nbases_covered")][inherited_noctrl_npos_bygene, on = c(symbol = "gene")]
inherited_noctrl_npos_bygene[, npos_perbase := npos / nbases_covered]
setnames(inherited_noctrl_npos_bygene, "symbol", "gene")
fwrite(inherited_noctrl_npos_bygene[order(-npos_perbase)], file = "Report/SNVs/impact/inherited_noctrl_vep_npos_bygene_perbase.csv")
inherited_noctrl_npos_bygene <- fread(file = "Report/SNVs/impact/inherited_noctrl_vep_npos_bygene_perbase.csv")

pdf("Report/SNVs/impact/inherited_noctrl_vep_npos_bygene_perbase.pdf", width = 5, height = 6)
ggplot(inherited_noctrl_npos_bygene, aes(x = gene, y = npos_perbase)) + geom_bar(stat = "identity") + scale_x_discrete(limits = inherited_noctrl_npos_bygene[order(-npos_perbase), gene]) + theme_classic(base_size = 16) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + xlab("") + ylab("# SNV sites per-base")
dev.off()

###########################################################################
## Transition probability matrix for inherited SNVs
###########################################################################
support_byposmut <- fread("Report/SNVs/filter/basediffperc_cutdemux_sub500k_q30_unstranded_highdepth_highaf_qcfltd_support_byposmut.csv")
inherited_noctrl_transmat <- dcast.data.table(support_byposmut[nmice >= 3, .N, by = c("ref", "alt")], ref ~ alt, value.var = "N")

vep <- fread(file = "Report/SNVs/impact/noctrl_vep_unique.csv")
inherited_noctrl_vep <- vep[nmice >= 3]
inherited_noctrl_vep_nonsynonymous <- inherited_noctrl_vep[class == "nonsynonymous"]
inherited_noctrl_vep_synonymous <- inherited_noctrl_vep[class  == "synonymous"]
inherited_noctrl_vep_tRNA <- inherited_noctrl_vep[class == "tRNA"]
inherited_noctrl_vep_rRNA <- inherited_noctrl_vep[class == "rRNA"]
inherited_noctrl_vep_dloop <- inherited_noctrl_vep[class == "D-loop"]

inherited_noctrl_vep_synonymous_posmut <- inherited_noctrl_vep_synonymous[, unique(posmut)]
inherited_noctrl_vep_nonsynonymous_posmut <- inherited_noctrl_vep_nonsynonymous[, unique(posmut)]
inherited_noctrl_vep_tRNA_posmut <- inherited_noctrl_vep_tRNA[, unique(posmut)]
inherited_noctrl_vep_rRNA_posmut <- inherited_noctrl_vep_rRNA[, unique(posmut)]
inherited_noctrl_vep_dloop_posmut <- inherited_noctrl_vep_dloop[, unique(posmut)]

inherited_noctrl_transmat_nonsynonymous <- dcast.data.table(support_byposmut[posmut %in% inherited_noctrl_vep_nonsynonymous_posmut, .N, by = c("ref", "alt")], ref ~ alt, value.var = "N")
inherited_noctrl_transmat_synonymous <- dcast.data.table(support_byposmut[posmut %in% inherited_noctrl_vep_synonymous_posmut, .N, by = c("ref", "alt")], ref ~ alt, value.var = "N")
inherited_noctrl_transmat_tRNA <- dcast.data.table(support_byposmut[posmut %in% inherited_noctrl_vep_tRNA_posmut, .N, by = c("ref", "alt")], ref ~ alt, value.var = "N")
inherited_noctrl_transmat_rRNA <- dcast.data.table(support_byposmut[posmut %in% inherited_noctrl_vep_rRNA_posmut, .N, by = c("ref", "alt")], ref ~ alt, value.var = "N")
inherited_noctrl_transmat_dloop <- dcast.data.table(support_byposmut[posmut %in% inherited_noctrl_vep_dloop_posmut, .N, by = c("ref", "alt")], ref ~ alt, value.var = "N")

Tools$write_xlsx(list(
    all = inherited_noctrl_transmat,
    nonsynonymous = inherited_noctrl_transmat_nonsynonymous, 
    synonymous = inherited_noctrl_transmat_synonymous, 
    tRNA = inherited_noctrl_transmat_tRNA, 
    rRNA = inherited_noctrl_transmat_rRNA, 
    dloop= inherited_noctrl_transmat_dloop
), file = "Report/SNVs/impact/inherited_noctrl_transmat.xlsx", row.names = FALSE)

inherited_noctrl_titv <- list(
    all = Genetics$titv(Tools$dt2df(inherited_noctrl_transmat)), 
    nonsynonymous = Genetics$titv(Tools$dt2df(inherited_noctrl_transmat_nonsynonymous)),
    synonymous = Genetics$titv(Tools$dt2df(inherited_noctrl_transmat_synonymous)),
    tRNA = Genetics$titv(Tools$dt2df(inherited_noctrl_transmat_tRNA)),
    rRNA = Genetics$titv(Tools$dt2df(inherited_noctrl_transmat_rRNA)),
    dloop = Genetics$titv(Tools$dt2df(inherited_noctrl_transmat_dloop))
)
inherited_noctrl_titv <- do.call(rbind, inherited_noctrl_titv)
inherited_noctrl_titv <- data.table(class = rownames(inherited_noctrl_titv), inherited_noctrl_titv)

fwrite(inherited_noctrl_titv, file = "Report/SNVs/impact/inherited_noctrl_titv.csv")
inherited_noctrl_titv <- fread(file = "Report/SNVs/impact/inherited_noctrl_titv.csv")

## Ti/Tv per gene
vep <- fread(file = "Report/SNVs/impact/noctrl_vep_unique.csv")
inherited_noctrl_vep <- vep[nmice >= 3]
inherited_noctrl_vep[, gene := SYMBOL]
inherited_noctrl_vep[, gene := ifelse(gene == "", class, gene)]
inherited_noctrl_vep[, table(gene)]
## gene
##  D-loop  mt-Co1  mt-Co2  mt-Co3 mt-Cytb  mt-Nd1  mt-Nd3  mt-Nd5  mt-Nd6 mt-Rnr2
##      32       5       3      19       3      11       6       7      30      17
##   mt-Tg  mt-Tl1   mt-Tm   mt-Tq
##      16      10       1       1
inherited_noctrl_transmat_bygene <- sapply(inherited_noctrl_vep[, sort(unique(gene))], function(g) dcast.data.table(inherited_noctrl_vep[gene == g, .N, by = c("ref", "alt")], ref ~ alt, value.var = "N"), simplify = FALSE)
inherited_noctrl_titv_bygene <- t(sapply(inherited_noctrl_transmat_bygene, function(dt) Genetics$titv(Tools$dt2df(dt))))
inherited_noctrl_titv_bygene <- data.table(gene = rownames(inherited_noctrl_titv_bygene), inherited_noctrl_titv_bygene)
fwrite(inherited_noctrl_titv_bygene, file = "Report/SNVs/impact/inherited_noctrl_titv_bygene.csv")
inherited_noctrl_titv_bygene <- fread(file = "Report/SNVs/impact/inherited_noctrl_titv_bygene.csv")

###########################################################################
## Among somatic SNVs, how is the per-base mutation occurrence per region?
###########################################################################
vep <- fread(file = "Report/SNVs/impact/noctrl_vep_unique.csv")
vep[, gene := ifelse(SYMBOL == "", class, SYMBOL)]
somatic_noctrl_vep <- vep[nmice == 1]
chrmgenes <- fread("Report/artifact/chrmgenes.csv")
## get D-loop and intergenic info
chrmgenes <- rbind(chrmgenes,
    chrmbases_properties[is_in_range == "Y" & is_in_primer == "N" & SNVID == "SNV8", list(symbol = "D-loop", start = min(pos), end = max(pos), strand = "*", nbases_covered = .N)], 
    chrmbases_properties[pos %in% 3843:3844, list(symbol = "intergenic", start = min(pos), end = max(pos), strand = "*", nbases_covered = .N)]
)
somatic_noctrl_npos_bygene <- somatic_noctrl_vep[, list(npos = uniqueN(pos)), by = c("gene")]
somatic_noctrl_npos_bygene <- chrmgenes[, c("symbol", "nbases_covered")][somatic_noctrl_npos_bygene, on = c(symbol = "gene")]
somatic_noctrl_npos_bygene[, npos_perbase := npos / nbases_covered]
setnames(somatic_noctrl_npos_bygene, "symbol", "gene")
fwrite(somatic_noctrl_npos_bygene[order(-npos_perbase)], file = "Report/SNVs/impact/somatic_noctrl_vep_npos_bygene_perbase.csv")
somatic_noctrl_npos_bygene <- fread(file = "Report/SNVs/impact/somatic_noctrl_vep_npos_bygene_perbase.csv")

pdf("Report/SNVs/impact/somatic_noctrl_vep_npos_bygene_perbase.pdf", width = 5, height = 6)
ggplot(somatic_noctrl_npos_bygene, aes(x = gene, y = npos_perbase)) + geom_bar(stat = "identity") + scale_x_discrete(limits = somatic_noctrl_npos_bygene[order(-npos_perbase), gene]) + theme_classic(base_size = 16) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + xlab("") + ylab("# SNV sites per-base")
dev.off()

###########################################################################
## Transition probability matrix for somatic SNVs
###########################################################################
support_byposmut <- fread("Report/SNVs/filter/basediffperc_cutdemux_sub500k_q30_unstranded_highdepth_highaf_qcfltd_support_byposmut.csv")
somatic_noctrl_transmat <- dcast.data.table(support_byposmut[nmice == 1, .N, by = c("ref", "alt")], ref ~ alt, value.var = "N")

vep <- fread(file = "Report/SNVs/impact/noctrl_vep_unique.csv")
somatic_noctrl_vep <- vep[nmice == 1]
somatic_noctrl_vep_nonsynonymous <- somatic_noctrl_vep[class == "nonsynonymous"]
somatic_noctrl_vep_synonymous <- somatic_noctrl_vep[class  == "synonymous"]
somatic_noctrl_vep_tRNA <- somatic_noctrl_vep[class == "tRNA"]
somatic_noctrl_vep_rRNA <- somatic_noctrl_vep[class == "rRNA"]
somatic_noctrl_vep_dloop <- somatic_noctrl_vep[class == "D-loop"]
somatic_noctrl_vep_intergenic <- somatic_noctrl_vep[class == "intergenic"]

somatic_noctrl_vep_synonymous_posmut <- somatic_noctrl_vep_synonymous[, unique(posmut)]
somatic_noctrl_vep_nonsynonymous_posmut <- somatic_noctrl_vep_nonsynonymous[, unique(posmut)]
somatic_noctrl_vep_tRNA_posmut <- somatic_noctrl_vep_tRNA[, unique(posmut)]
somatic_noctrl_vep_rRNA_posmut <- somatic_noctrl_vep_rRNA[, unique(posmut)]
somatic_noctrl_vep_dloop_posmut <- somatic_noctrl_vep_dloop[, unique(posmut)]
somatic_noctrl_vep_intergenic_posmut <- somatic_noctrl_vep_intergenic[, unique(posmut)]

somatic_noctrl_transmat_nonsynonymous <- dcast.data.table(support_byposmut[posmut %in% somatic_noctrl_vep_nonsynonymous_posmut, .N, by = c("ref", "alt")], ref ~ alt, value.var = "N")
somatic_noctrl_transmat_synonymous <- dcast.data.table(support_byposmut[posmut %in% somatic_noctrl_vep_synonymous_posmut, .N, by = c("ref", "alt")], ref ~ alt, value.var = "N")
somatic_noctrl_transmat_tRNA <- dcast.data.table(support_byposmut[posmut %in% somatic_noctrl_vep_tRNA_posmut, .N, by = c("ref", "alt")], ref ~ alt, value.var = "N")
somatic_noctrl_transmat_rRNA <- dcast.data.table(support_byposmut[posmut %in% somatic_noctrl_vep_rRNA_posmut, .N, by = c("ref", "alt")], ref ~ alt, value.var = "N")
somatic_noctrl_transmat_dloop <- dcast.data.table(support_byposmut[posmut %in% somatic_noctrl_vep_dloop_posmut, .N, by = c("ref", "alt")], ref ~ alt, value.var = "N")
somatic_noctrl_transmat_intergenic <- dcast.data.table(support_byposmut[posmut %in% somatic_noctrl_vep_intergenic_posmut, .N, by = c("ref", "alt")], ref ~ alt, value.var = "N")

Tools$write_xlsx(list(
    all = somatic_noctrl_transmat,
    nonsynonymous = somatic_noctrl_transmat_nonsynonymous, 
    synonymous = somatic_noctrl_transmat_synonymous, 
    tRNA = somatic_noctrl_transmat_tRNA, 
    rRNA = somatic_noctrl_transmat_rRNA, 
    dloop = somatic_noctrl_transmat_dloop,
    intergenic = somatic_noctrl_transmat_intergenic
), file = "Report/SNVs/impact/somatic_noctrl_transmat.xlsx", row.names = FALSE)

somatic_noctrl_titv <- list(
    all = Genetics$titv(Tools$dt2df(somatic_noctrl_transmat)), 
    nonsynonymous = Genetics$titv(Tools$dt2df(somatic_noctrl_transmat_nonsynonymous)),
    synonymous = Genetics$titv(Tools$dt2df(somatic_noctrl_transmat_synonymous)),
    tRNA = Genetics$titv(Tools$dt2df(somatic_noctrl_transmat_tRNA)),
    rRNA = Genetics$titv(Tools$dt2df(somatic_noctrl_transmat_rRNA)),
    dloop = Genetics$titv(Tools$dt2df(somatic_noctrl_transmat_dloop)), 
    intergenic = Genetics$titv(Tools$dt2df(somatic_noctrl_transmat_intergenic)) 
)
somatic_noctrl_titv <- do.call(rbind, somatic_noctrl_titv)
somatic_noctrl_titv <- data.table(class = rownames(somatic_noctrl_titv), somatic_noctrl_titv)

fwrite(somatic_noctrl_titv, file = "Report/SNVs/impact/somatic_noctrl_titv.csv")
somatic_noctrl_titv <- fread(file = "Report/SNVs/impact/somatic_noctrl_titv.csv")

## Ti/Tv per gene
vep <- fread(file = "Report/SNVs/impact/noctrl_vep_unique.csv")
somatic_noctrl_vep <- vep[nmice == 1]
somatic_noctrl_vep[, gene := SYMBOL]
somatic_noctrl_vep[, gene := ifelse(gene == "", class, gene)]
somatic_noctrl_vep[, table(gene)]
## gene
##     D-loop intergenic     mt-Co1     mt-Co2     mt-Co3    mt-Cytb     mt-Nd1
##         49          1         49         58         64         31         50
##     mt-Nd3     mt-Nd5     mt-Nd6    mt-Rnr2      mt-Tg      mt-Ti     mt-Tl1
##         14         54         68         92         32          1         30
##      mt-Tm      mt-Tq      mt-Tt
##         11         27          1
somatic_noctrl_transmat_bygene <- sapply(somatic_noctrl_vep[, sort(unique(gene))], function(g) dcast.data.table(somatic_noctrl_vep[gene == g, .N, by = c("ref", "alt")], ref ~ alt, value.var = "N"), simplify = FALSE)
somatic_noctrl_titv_bygene <- t(sapply(somatic_noctrl_transmat_bygene, function(dt) Genetics$titv(Tools$dt2df(dt))))
somatic_noctrl_titv_bygene <- data.table(gene = rownames(somatic_noctrl_titv_bygene), somatic_noctrl_titv_bygene)
fwrite(somatic_noctrl_titv_bygene, file = "Report/SNVs/impact/somatic_noctrl_titv_bygene.csv")
somatic_noctrl_titv_bygene <- fread(file = "Report/SNVs/impact/somatic_noctrl_titv_bygene.csv")

###########################################################################
## somatic vs. inherited SNVs
###########################################################################
somatic_noctrl_nsnvs_bygene <- somatic_noctrl_vep[, list(nsnvs = .N, origin = "somatic"), by = c("gene")]
somatic_noctrl_nsnvs_bygene[, psnvs := nsnvs / sum(nsnvs) * 100]
inherited_noctrl_nsnvs_bygene <- inherited_noctrl_vep[, list(nsnvs = .N, origin = "inherited"), by = c("gene")]
inherited_noctrl_nsnvs_bygene[, psnvs := nsnvs / sum(nsnvs) * 100]
comparison_noctrl_nsnvs_bygene <- rbind(somatic_noctrl_nsnvs_bygene, inherited_noctrl_nsnvs_bygene)
comparison_noctrl_nsnvs_bygene[, origin := factor(origin, levels = c("somatic", "inherited"))]

pdf("Report/SNVs/impact/comparison_noctrl_nsnvs_bygene_bar.pdf", width = 7, height = 6)
ggplot(comparison_noctrl_nsnvs_bygene, aes(x = gene, fill = origin, y = psnvs)) + geom_bar(stat = "identity", position = position_dodge2(preserve = "single")) + scale_x_discrete(limits = c(inherited_noctrl_nsnvs_bygene[order(-psnvs), gene], setdiff(somatic_noctrl_nsnvs_bygene[, gene], inherited_noctrl_nsnvs_bygene[, gene]))) + theme_classic(base_size = 16) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + xlab("") + ylab("% of SNVs from either origin") + scale_fill_brewer(palette = "Set1")
dev.off()

## SNV rate: expected SNV sites per gene, observed SNV sites per gene for somatic and for inherited SNVs
somatic_noctrl_npos_bygene[, npos_expected := somatic_noctrl_npos_bygene[, sum(npos)] * nbases_covered / sum(nbases_covered)]
inherited_noctrl_npos_bygene[, npos_expected := inherited_noctrl_npos_bygene[, sum(npos)] * nbases_covered / sum(nbases_covered)]
somatic_noctrl_npos_bygene[, ratio := npos / npos_expected]
inherited_noctrl_npos_bygene[, ratio := npos / npos_expected]
somatic_noctrl_npos_bygene[, origin := "somatic"]
inherited_noctrl_npos_bygene[, origin := "inherited"]
comparison_noctrl_npos_bygene <- rbind(somatic_noctrl_npos_bygene, inherited_noctrl_npos_bygene)

pdf("Report/SNVs/impact/comparison_noctrl_npos_bygene_bar.pdf", width = 7, height = 6)
ggplot(comparison_noctrl_npos_bygene, aes(x = gene, y = log2(ratio), fill = origin)) + geom_bar(stat = "identity", position = position_dodge2(preserve = "single")) + scale_x_discrete(limits = c(inherited_noctrl_npos_bygene[order(-ratio), gene], setdiff(somatic_noctrl_npos_bygene[, gene], inherited_noctrl_npos_bygene[, gene])[c(3, 1, 2)])) + theme_classic(base_size = 16) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + xlab("") + ylab("normalized mutation rate (log2)") + scale_fill_brewer(palette = "Set1")
dev.off()

## SNV rate by functional class (coding vs noncoding) for inherited vs somatic SNVs
comparison_noctrl_npos_bygene <- fread(file = "Report/20210510/SNVs/impact/comparison_noctrl_npos_bygene.csv")
comparison_noctrl_npos_byclass <- rbind(
    comparison_noctrl_npos_bygene[origin == "inherited" & gene %like% "Atp|Nd|Co|Cytb", list(origin = "inherited", class = "coding", nbases_covered = sum(nbases_covered), npos = sum(npos))], 
    comparison_noctrl_npos_bygene[origin == "inherited" & ! gene %like% "Atp|Nd|Co|Cytb", list(origin = "inherited", class = "noncoding", nbases_covered = sum(nbases_covered), npos = sum(npos))], 
    comparison_noctrl_npos_bygene[origin == "somatic" & gene %like% "Atp|Nd|Co|Cytb", list(origin = "somatic", class = "coding", nbases_covered = sum(nbases_covered), npos = sum(npos))], 
    comparison_noctrl_npos_bygene[origin == "somatic" & ! gene %like% "Atp|Nd|Co|Cytb", list(origin = "somatic", class = "noncoding", nbases_covered = sum(nbases_covered), npos = sum(npos))]
)

pdf("Report/20210510/SNVs/impact/comparison_noctrl_npos_byclass_bar.pdf", width = 6, height = 6)
ggplot(comparison_noctrl_npos_byclass, aes(x = origin, y = npos / nbases_covered, fill = class)) + geom_bar(stat = "identity", position = "dodge") + theme_classic(16) + xlab("") + ylab("Number of SNV site per base") + scale_fill_brewer(palette = "Set1")
dev.off()
fisher.test(as.matrix(comparison_noctrl_npos_byclass[origin == "inherited", list(nbases_covered, npos)]))
## 
##         Fisher's Exact Test for Count Data
## 
## data:  as.matrix(comparison_noctrl_npos_byclass[origin == "inherited", list(nbases_covered, npos)])
## p-value = 0.02842
## alternative hypothesis: true odds ratio is not equal to 1
## 95 percent confidence interval:
##  1.033391 2.083782
## sample estimates:
## odds ratio
##   1.468354
## 
fisher.test(as.matrix(comparison_noctrl_npos_byclass[origin == "somatic", list(nbases_covered, npos)]))
## 
##         Fisher's Exact Test for Count Data
## 
## data:  as.matrix(comparison_noctrl_npos_byclass[origin == "somatic", list(nbases_covered, npos)])
## p-value = 0.8377
## alternative hypothesis: true odds ratio is not equal to 1
## 95 percent confidence interval:
##  0.8317907 1.2553468
## sample estimates:
## odds ratio
##   1.022407


## ti/tv comparison
inherited_noctrl_titv <- fread(file = "Report/SNVs/impact/inherited_noctrl_titv.csv")
somatic_noctrl_titv <- fread(file = "Report/SNVs/impact/somatic_noctrl_titv.csv")
inherited_noctrl_titv[, origin := "inherited"]
somatic_noctrl_titv[, origin := "somatic"]
comparison_noctrl_titv <- rbind(inherited_noctrl_titv, somatic_noctrl_titv)
comparison_noctrl_titv[class == "all", class := "average"]
comparison_noctrl_titv[class == "synonymous", class := "syn"]
comparison_noctrl_titv[class == "nonsynonymous", class := "nonsyn"]
comparison_noctrl_titv[class == "dloop", class := "D-loop"]
pdf("Report/SNVs/impact/comparison_noctrl_titv_bar.pdf", width = 7, height = 6)
ggplot(comparison_noctrl_titv, aes(x = class, fill = origin, y = r)) + geom_bar(stat = "identity", position = position_dodge2(preserve = "single")) + scale_x_discrete(limits = c(comparison_noctrl_titv[origin == "inherited"][order(-r), class], "intergenic")) + theme_classic(base_size = 16)  + xlab("") + ylab("Ti/Tv") + scale_fill_brewer(palette = "Set1") + scale_y_continuous(breaks = seq(0, 12, by = 2)) ## + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
dev.off()
(sapply(c("all", "tRNA", "rRNA", "nonsynonymous", "synonymous", "dloop"), function(x) { inherited <- comparison_noctrl_titv[class == x & origin == "inherited", list(ti, tv)]; somatic <- comparison_noctrl_titv[class == x & origin == "somatic", list(ti, tv)]; y <- rbind(inherited, somatic); fisher.test(y)$p.value }, simplify = TRUE))
##           all          tRNA          rRNA nonsynonymous    synonymous
##   0.119299923   0.023629020   0.001195984   0.034943804   0.738647551
##         dloop
##   0.001920203

inherited_noctrl_titv_bygene <- fread(file = "Report/SNVs/impact/inherited_noctrl_titv_bygene.csv")
somatic_noctrl_titv_bygene <- fread(file = "Report/SNVs/impact/somatic_noctrl_titv_bygene.csv")
inherited_noctrl_titv_bygene[, origin := "inherited"]
somatic_noctrl_titv_bygene[, origin := "somatic"]
comparison_noctrl_titv_bygene <- rbind(inherited_noctrl_titv_bygene, somatic_noctrl_titv_bygene)

pdf("Report/SNVs/impact/comparison_noctrl_titv_bygene_bar.pdf", width = 7, height = 6)
ggplot(comparison_noctrl_titv_bygene, aes(x = gene, fill = origin, y = r)) + geom_bar(stat = "identity", position = position_dodge2(preserve = "single")) + scale_x_discrete(limits = setdiff(unique(c(comparison_noctrl_titv_bygene[origin == "inherited"][order(-r), gene], comparison_noctrl_titv_bygene[origin == "somatic"][order(r), gene])), c("D-loop", "intergenic"))) + theme_classic(base_size = 16)  + xlab("") + ylab("Ti/Tv") + scale_fill_brewer(palette = "Set1") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + scale_y_continuous(breaks = seq(0, 12, by = 2))
dev.off()

## Ti/Tv by functional class (coding vs noncoding) for inherited vs somatic SNVs
comparison_noctrl_titv <- fread(file = "Report/20210510/SNVs/impact/comparison_noctrl_titv.csv")
comparison_noctrl_titv_byclass <- rbind(
    comparison_noctrl_titv[class != "average" & origin == "inherited" & class %in% c("nonsyn", "syn"), list(origin = "inherited", class = "coding", ti = sum(ti), tv = sum(tv), r = sum(ti)/sum(tv))], 
    comparison_noctrl_titv[class != "average" & origin == "inherited" & !class %in% c("nonsyn", "syn"), list(origin = "inherited", class = "noncoding", ti = sum(ti), tv = sum(tv), r = sum(ti)/sum(tv))], 
    comparison_noctrl_titv[class != "average" & origin == "somatic" & class %in% c("nonsyn", "syn"), list(origin = "somatic", class = "coding", ti = sum(ti), tv = sum(tv), r = sum(ti)/sum(tv))], 
    comparison_noctrl_titv[class != "average" & origin == "somatic" & !class %in% c("nonsyn", "syn"), list(origin = "somatic", class = "noncoding", ti = sum(ti), tv = sum(tv), r = sum(ti)/sum(tv))]
)
pdf("Report/20210510/SNVs/impact/comparison_noctrl_titv_byclass_bar.pdf", width = 6, height = 6)
ggplot(comparison_noctrl_titv_byclass, aes(x = origin, fill = class, y = r)) + geom_bar(stat = "identity", position = "dodge") + theme_classic(base_size = 16)  + xlab("") + ylab("Ti/Tv") + scale_fill_brewer(palette = "Set1")
dev.off()

fisher.test(as.matrix(comparison_noctrl_titv_byclass[origin == "inherited", list(ti, tv)]))
## 
## 
##         Fisher's Exact Test for Count Data
## 
## data:  as.matrix(comparison_noctrl_titv_byclass[origin == "inherited", list(ti, tv)])
## p-value = 0.1006
## alternative hypothesis: true odds ratio is not equal to 1
## 95 percent confidence interval:
##  0.8390188 4.1757872
## sample estimates:
## odds ratio
##   1.852849
## 
fisher.test(as.matrix(comparison_noctrl_titv_byclass[origin == "somatic", list(ti, tv)]))
## 
##         Fisher's Exact Test for Count Data
## 
## data:  as.matrix(comparison_noctrl_titv_byclass[origin == "somatic", list(ti, tv)])
## p-value = 0.7905
## alternative hypothesis: true odds ratio is not equal to 1
## 95 percent confidence interval:
##  0.7390601 1.5233002
## sample estimates:
## odds ratio
##   1.062389
## 
pchisq(-2 * (log(0.01) + log(1-0.79)), df = 4, lower.tail = FALSE)
## [1] 0.01504822

###########################################################################
## How about the 2017 dataset?
###########################################################################
Morris2017_altperc_bymito_byposmut <- fread(file = "Report/Morris2017/common_altperc_bymito_byposmut.csv")
Morris2017_common_posmut <- data.table(posmut = names(Morris2017_altperc_bymito_byposmut)[-c(1:10)])
Morris2017_common_posmut[, pos := as.integer(sapply(strsplit(posmut, ":"), "[", 1))]
Morris2017_common_posmut[, ref := sapply(strsplit(sapply(strsplit(posmut, ":"), "[", 2), ">"), "[", 1)]
Morris2017_common_posmut[, alt := sapply(strsplit(sapply(strsplit(posmut, ":"), "[", 2), ">"), "[", 2)]
fwrite(Morris2017_common_posmut, file = "Report/SNVs/impact/Morris2017_common_posmut.csv")
Morris2017_common_transmat <- dcast.data.table(Morris2017_common_posmut[, .N, by = c("ref", "alt")], ref ~ alt, value.var = "N")
fwrite(Morris2017_common_transmat, file = "Report/SNVs/impact/Morris2017_common_transmat.csv")
Morris2017_common_transmat <- fread(file = "Report/SNVs/impact/Morris2017_common_transmat.csv")
setkey(Morris2017_common_transmat, ref)

Morris2017_common_titv <- Genetics$titv(Tools$dt2df(Morris2017_common_transmat))
fwrite(t(Morris2017_common_titv), file = "Report/SNVs/impact/Morris2017_common_titv.csv")
Morris2017_common_titv <- fread(file = "Report/SNVs/impact/Morris2017_common_titv.csv")

Morris2017_common_posmut[, mut := factor(mut, levels = paste0(rep(c("A", "C", "G", "T"), each = 5), ">", c("A", "C", "G", "T", "del")))]
pdf(file = "Report/SNVs/impact/Morris2017_common_transmat_barplot.pdf", width = 6, height = 1.5)
par(ps = 12, lend = 2, ljoin = 1, bty = "L", mfrow = c(1, 1), mar = c(1.5, 2.0, 1, 0), oma = c(0, 0, 0, 0), mgp = c(1.0, 0.5, 0), cex.axis = 0.8)
Morris2017_common_posmut[, barplot(table(mut), ylab = "# of SNVs", main = "Morris2017", cex.names = 0.5)]
dev.off()
