###########################################################################
## Author: Youtao Lu <luyoutao@sas.upenn.edu>
##  
## Copyright (c) 2021-2023, Youtao Lu and Junhyong Kim, Department of Biology, University of Pennsylvania
## Copyright (c) 2021-2023, James Eberwine, Perelman School of Medicine, University of Pennsylvania
## Copyright (c) 2021-2023, David Issadore, Penn Engineering, University of Pennsylvania
## 
## This Source Code is distributed under Creative Commons Attribution License 4.0 (CC BY).
###########################################################################
source("Source/release/functions.R")
library("data.table")
library("VennDiagram")
library("ggplot2")
library("pheatmap")
library("ape")

chrmproperties <- fread(file = "Report/release/artifact/chrmbases_properties.csv.gz")
support_byposmut <- fread(file = "Report/release/SNVs/filter/basediffperc_cutdemux_sub500k_q30_unstranded_highdepth_highaf_qcfltd_support_byposmut.csv")
support_bypos <- fread("Report/release/SNVs/filter/basediffperc_cutdemux_sub500k_q30_unstranded_highdepth_highaf_qcfltd_support_bypos.csv")
smito_pos <- support_byposmut[, unique(pos)]

###########################################################################
## multiple sequence alignment for 17 mouse strains
###########################################################################
strains_mclust <- read.csv("Report/release/strains/clustalo-17strains_str.tsv", sep = "\t", header = FALSE, as.is = TRUE)
bases <- lapply(strains_mclust[, 2], function(x) strsplit(x, "")[[1]])
unique(lengths(bases))
## [1] 16315
bases <- do.call(cbind, bases)
strains <- strains_mclust[, 1]
colnames(bases) <- strains
pos_C57BL_6J <- cumsum(ifelse(bases[, which(strains == "C57BL/6J")] == "-", 0, 1))
strains_mclust <- data.frame(pos_C57BL_6J = pos_C57BL_6J, bases, stringsAsFactors = FALSE, check.names = FALSE)
fwrite(strains_mclust, file = "Report/release/strains/clustalo-17strains_base.csv")
strains_mclust <- fread(file = "Report/release/strains/clustalo-17strains_base.csv")

###########################################################################
## Make a fan tree 
###########################################################################
tree_17strains <- ape::read.tree("Report/release/strains/clustalo-17strains.ph")
pdf("Report/release/strains/clustalo-17strains.pdf")
par(ps = 16, lend = 2, ljoin = 1, bty = "L", mfrow = c(1, 1), mar = c(2.5, 2.5, 1, 0), oma = c(0, 0, 0, 0), mgp = c(1.5, 0.5, 0))
plot.phylo(tree_17strains, type = "fan", use.edge.length = FALSE)
dev.off()

###########################################################################
## SMITO-evo site overlap
###########################################################################
strains_mclust <- fread(file = "Report/release/strains/clustalo-17strains_base.csv")
setDT(strains_mclust)
setkey(strains_mclust, pos_C57BL_6J)
strains_mclust_varonly <- subset(strains_mclust, alignment != "*" & `C57BL/6J` != "-")
strains_mclust_varonly <- strains_mclust_varonly[, -c("alignment")]

fig <- venn.diagram(list(SMITO = smito_pos, `all_17strains_inrange` = intersect(chrmproperties[is_in_range == "Y" & is_in_primer == "N", pos], strains_mclust_varonly[, unique(pos_C57BL_6J)])), filename = NULL, main.fontfamily = "sans", sub.fontfamily = "sans", cat.fontfamily = "sans", hyper.test = TRUE, total.population = length(chrmproperties[is_in_range == "Y" & is_in_primer == "N", pos]), lower.tail = FALSE, main.cex = 1.5, cat.cex = 1, fill = c("#984EA3", "#FFFFE0"), cex = c(1, 1, 1) * 1.5, alpha = 0.5, lwd = 1)
ggsave(fig, filename= "Report/release/strains/smito_all_vs_17strains-inrange_venn.pdf", width = 6, height = 6)

## when SMITO is stratified by the number of mouse support
smito_strains_x <- sapply(1:13, function(x) support_byposmut[nmice >= x, length(intersect(unique(pos), intersect(chrmproperties[is_in_range == "Y" & is_in_primer == "N", pos], strains_mclust_varonly[, unique(pos_C57BL_6J)])))])
smito_strains_k <- sapply(1:13, function(x) support_byposmut[nmice >= x, uniqueN(pos)])
smito_strains_m <- uniqueN(intersect(chrmproperties[is_in_range == "Y" & is_in_primer == "N", pos], strains_mclust_varonly[, unique(pos_C57BL_6J)]))
smito_strains_n <- chrmproperties[is_in_range == "Y" & is_in_primer == "N", .N] - smito_strains_m
smito_strains_pval <- phyper(q = smito_strains_x, k = smito_strains_k, n = rep(smito_strains_n, 13), m = rep(smito_strains_m, 13), lower.tail = FALSE)
smito_strains_obs <- smito_strains_x / smito_strains_k
smito_strains_bkg <- smito_strains_m / chrmproperties[is_in_range == "Y" & is_in_primer == "N", .N]
smito_strains_oddsratio <- smito_strains_obs / smito_strains_bkg

smito_strains_overlap <- data.table(
    nmice = 1:13, smito = smito_strains_k, `17strains` = smito_strains_m, 
    overlap = smito_strains_x, 
    frac_smito = smito_strains_obs, 
    frac_17strains = smito_strains_bkg, 
    oddsratio = smito_strains_oddsratio,
    pval = smito_strains_pval
)

fwrite(smito_strains_overlap, file = "Report/release/strains/smito_17strains-inrange_overlap.csv")
smito_strains_overlap <- fread(file = "Report/release/strains/smito_17strains-inrange_overlap.csv")

pdf(file = "Report/release/strains/smito_17strains-inrange_overlap_barplot.pdf", height = 6, width = 6)
par(ps = 16, lend = 2, ljoin = 1, bty = "L", mfrow = c(1, 1), mar = c(2.5, 2.5, 1, 2.5), oma = c(0, 0, 0, 0), mgp = c(1.5, 0.5, 0))
x <- barplot(smito_strains_overlap[, oddsratio], names.arg = 1:13, border = FALSE)
lines(x = x, y = smito_strains_overlap[, pval * max(oddsratio)/max(pval)])
points(x = x, y = smito_strains_overlap[, pval * max(oddsratio)/max(pval)])
axis(side = 4, at = seq(0, smito_strains_overlap[, 0.6 * max(oddsratio) / max(pval)], length.out = 4), labels = seq(0, 0.6, by = 0.2))
title(ylab = "overlap odds ratio")
mtext("hypergeometric p-value", side = 4, line = 1.5)
mtext("# mice", side = 1, line = 1.5)
dev.off()

###########################################################################
## Test whether 
## (1) the major alleles
## (2) the top minor alleles 
## from SMITO match the evolutionary data
###########################################################################
highdepth_highaf <- fread("Report/release/SNVs/filter/basediffperc_cutdemux_sub500k_q30_unstranded_highdepth_highaf_qcfltd.csv.gz")
highdepth_highaf_noctrl_17strains <- highdepth_highaf[IsCtrl == "N" & pos %in% strains_mclust_varonly[, unique(pos_C57BL_6J)]]
highdepth_highaf_noctrl_17strains <- strains_mclust_varonly[highdepth_highaf_noctrl_17strains, on = c("pos_C57BL_6J" = "pos")]
setnames(highdepth_highaf_noctrl_17strains, "pos_C57BL_6J", "pos")

highdepth_highaf_noctrl_17strains[ref == "A", A := `=`]
highdepth_highaf_noctrl_17strains[ref == "C", C := `=`]
highdepth_highaf_noctrl_17strains[ref == "G", G := `=`]
highdepth_highaf_noctrl_17strains[ref == "T", T := `=`]

smito_avgaf <- highdepth_highaf_noctrl_17strains[, { X <- colSums(t(sapply(1:nrow(.SD), function(x) { unlist(.SD[x, 1]) * unlist(.SD[x, -1]) }))); as.list(100 * X / sum(X)) }, keyby = "pos", .SDcols = c("depth", "A", "C", "G", "T", "del")]
smito_avgaf <- support_bypos[smito_avgaf, on = "pos"]

evo17strains <- strains_mclust_varonly[pos_C57BL_6J %in% highdepth_highaf_noctrl_17strains[, pos]]
setnames(evo17strains, "pos_C57BL_6J", "pos")

evo17strains_avgaf <- evo17strains[, as.list(table(factor(.SD, levels = c("A", "C", "G", "T", "-"))) / 17 * 100), keyby = "pos"]
setnames(evo17strains_avgaf, "-", "del")

smito_evo17strains_avgaf <- merge.data.table(smito_avgaf, evo17strains_avgaf, by.x = "pos", by.y = "pos")
setnames(smito_evo17strains_avgaf, 
    c("A.x", "C.x", "G.x", "T.x", "del.x", "A.y", "C.y", "G.y", "T.y", "del.y"), 
    c("smito_A", "smito_C", "smito_G", "smito_T", "smito_del", "evo_A", "evo_C", "evo_G", "evo_T", "evo_del")
)

smito_evo17strains_avgaf <- chrmproperties[smito_evo17strains_avgaf[, -"SNVID"], on = "pos"]

smito_major <- smito_evo17strains_avgaf[, apply(.SD, 1, function(x) c("A", "C", "G", "T", "del")[which.max(x)]), .SD = c("smito_A", "smito_C", "smito_G", "smito_T", "smito_del")]
evo17strains_major <- smito_evo17strains_avgaf[, apply(.SD, 1, function(x) c("A", "C", "G", "T", "del")[which.max(x)]), .SD = c("evo_A", "evo_C", "evo_G", "evo_T", "evo_del")]
smito_topminor <- smito_evo17strains_avgaf[, apply(.SD, 1, function(x) { z <- c("A", "C", "G", "T", "del"); y <- match(x[1], z); i <- which.max(x[-1][-y]); z[-y][i] }), .SDcols = c("ref", "smito_A", "smito_C", "smito_G", "smito_T", "smito_del")]
evo17strains_topminor <- smito_evo17strains_avgaf[, apply(.SD, 1, function(x) { z <- c("A", "C", "G", "T", "del"); y <- match(x[1], z); i <- which.max(x[-1][-y]); z[-y][i] }), .SDcols = c("ref", "evo_A", "evo_C", "evo_G", "evo_T", "evo_del")]
smito_evo17strains_avgaf[, smito_major := smito_major]
smito_evo17strains_avgaf[, evo_major := evo17strains_major]
smito_evo17strains_avgaf[, smito_topminor := smito_topminor]
smito_evo17strains_avgaf[, evo_topminor := evo17strains_topminor]

fwrite(smito_evo17strains_avgaf, file = "Report/release/strains/smito_17strains_avgaf.csv")
smito_evo17strains_avgaf <- fread(file = "Report/release/strains/smito_17strains_avgaf.csv")
