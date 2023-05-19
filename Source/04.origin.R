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
library("ggplot2")
library("pheatmap")
library("circlize")
library("viridisLite")

snv_info <- fread("Data/snv_loci_v2.csv")
snvIDs <- snv_info[, SNVID]
mito_barcodes <- fread("Data/mito_barcodes.csv")
mitoIDs <- mito_barcodes[, ID]

chrmbases_properties <- fread("Report/release/artifact/chrmbases_properties.csv.gz")
chrmbases <- chrmbases_properties[, ref]
nchrmbases <- length(chrmbases)
nchrmbases
## [1] 16299

MitoInfo <- fread("Report/release/metadata/MitoInfo.csv")
dim(MitoInfo)
## [1] 1717   19
MitoInfo[, ExptID := factor(ExptID)]
MitoInfo[, MitoID := factor(MitoID, levels = mitoIDs)]
MitoInfo[, CellID := factor(CellID)]

CellInfo <- fread("Report/release/metadata/CellInfo.csv")
dim(CellInfo)
## [1] 102  12
MouseInfo <- fread("Report/release/metadata/MouseInfo.csv")
dim(MouseInfo)
## [1] 13  2

MitoInfo_df <- as.data.frame(MitoInfo)
rownames(MitoInfo_df) <- MitoInfo_df[, "LibraryMitoID"]

CellInfo_df <- as.data.frame(CellInfo)
rownames(CellInfo_df) <- CellInfo_df$CellUID

MouseInfo_df <- as.data.frame(MouseInfo)
rownames(MouseInfo_df) <- MouseInfo_df$MouseID

###########################################################################
## What is the theoretical expectation of the # of SNVs shared by k mice?
## Let's do an empirical distribution.
###########################################################################
support_byposmut <- fread(file = "Report/release/SNVs/filter/basediffperc_cutdemux_sub500k_q30_unstranded_highdepth_highaf_qcfltd_support_byposmut.csv")
support_byposmut[, .N]
## [1] 1032
highdepth_qcfltd_altperc_bymito_byposmut <- fread(file = "Report/release/SNVs/filter/basediffperc_cutdemux_sub500k_q30_unstranded_highdepth_qcfltd_altperc_bymito_byposmut.csv.gz")
dim(highdepth_qcfltd_altperc_bymito_byposmut)
## [1] 1702 1051
X <- as.matrix(highdepth_qcfltd_altperc_bymito_byposmut[IsCtrl == "N", -c(1:19)])
dim(X)
## [1] 1633 1032
vaf_th <- 0.05
Y <- matrix(ifelse(X >= vaf_th * 100, 1, 0), ncol = ncol(X), dimnames = dimnames(X))
Y <- data.table(highdepth_qcfltd_altperc_bymito_byposmut[IsCtrl == "N", 1:19], Y)
Z <- Y[, lapply(.SD, function(x) ifelse(any(x > 0, na.rm = TRUE), 1, 0)), keyby = "MouseID", .SDcols = 20:1051]
nperms <- 100
res <- replicate(nperms, expr = { colSums(t(apply(Z[, -1], 1, sample))) })
dim(res)
## [1] 1032  100
nposmut_bynmice <- apply(res, 2, function(x) table(factor(x, levels = 1:support_byposmut[, max(nmice)])))
nposmut_bynmice_mean <- rowMeans(nposmut_bynmice)

pdf("Report/release/SNVs/origin/highdepth_highaf_noctrl_nposmut_bynmice_obs-vs-null.pdf", width = 6, height = 6)
par(ps = 16, lend = 2, ljoin = 1, bty = "L", mfrow = c(1, 1), mar = c(2.5, 2.5, 0.5, 0.5), oma = c(0, 0, 0, 0), mgp = c(1.5, 0.5, 0))
plot(nposmut_bynmice_mean, type = "b", xlab = "# mice sharing", ylab = "# SNVs", pch = 2, col = "black", ylim = c(0, 640))
support_byposmut[, .N, by = nmice][order(nmice)][, points(N ~ nmice, xlab = "# mice", ylab = "# SNVs", type = 'b', pch = 1, col = "red")]
abline(h = 0, lty = 2, col = "gray")
legend("topright", legend = c("expected", "observed"), pch = c(2, 1), col = c("black", "red"), box.col = NA)
dev.off()

###########################################################################
## filter total SNVs by the number of support: (1) inherited, nmice >= 3;
## (2) somatic nmice == 1
## main data
###########################################################################
## 1. inherited SNVs
highdepth_qcfltd <- fread("Report/release/SNVs/filter/basediffperc_cutdemux_sub500k_q30_unstranded_highdepth_qcfltd.csv.gz")
highdepth_qcfltd[, alt := factor(alt, levels = c("A", "C", "G", "T", "del"))]
highdepth_qcfltd[, mut := factor(paste0(ref, ">", alt), levels = paste0(rep(c("A", "C", "G", "T"), each = 5), ">", c("A", "C", "G", "T", "del")))]
randed_highdepth_qcfltd, pos, alt)
highdepth_qcfltd[, posmut := paste0(pos, ":", mut)]
highdepth_qcfltd[, posmut := factor(posmut, levels = unique(posmut))]
dim(highdepth_qcfltd)
## [1] 1021830    28

support_byposmut <- fread(file = "Report/release/SNVs/filter/basediffperc_cutdemux_sub500k_q30_unstranded_highdepth_highaf_qcfltd_support_byposmut.csv")
posmut_inherited <- support_byposmut[nmice >= 3, posmut]
pos_inherited <- support_byposmut[nmice >= 3, unique(pos)]

highdepth_inherited <- highdepth_qcfltd[as.character(posmut) %in% posmut_inherited]
dim(highdepth_inherited)
## [1] 161934     28
highdepth_inherited <- highdepth_inherited[, -c(27:28)]
fwrite(highdepth_inherited, file = "Report/release/SNVs/origin/highdepth_inherited.csv.gz")
highdepth_inherited <- fread(file = "Report/release/SNVs/origin/highdepth_inherited.csv.gz")

highdepth_qcfltd_altperc_bymito_byposmut <- fread("Report/release/SNVs/filter/basediffperc_cutdemux_sub500k_q30_unstranded_highdepth_qcfltd_altperc_bymito_byposmut.csv.gz")
highdepth_qcfltd_allele_bymito_bypos <- fread("Report/release/SNVs/filter/basediffperc_cutdemux_sub500k_q30_unstranded_highdepth_qcfltd_allele_bymito_bypos.csv.gz")

highdepth_inherited_altperc_bymito_byposmut <- highdepth_qcfltd_altperc_bymito_byposmut[, c(1:19, match(posmut_inherited, colnames(highdepth_qcfltd_altperc_bymito_byposmut))), with = FALSE]
fwrite(highdepth_inherited_altperc_bymito_byposmut, file = "Report/release/SNVs/origin/highdepth_inherited_altperc_bymito_byposmut.csv.gz")
highdepth_inherited_altperc_bymito_byposmut <- fread(file = "Report/release/SNVs/origin/highdepth_inherited_altperc_bymito_byposmut.csv.gz")

highdepth_inherited_allele_bymito_bypos <- highdepth_qcfltd_allele_bymito_bypos[, c(1:19, match(pos_inherited, colnames(highdepth_qcfltd_allele_bymito_bypos))), with = FALSE]
fwrite(highdepth_inherited_allele_bymito_bypos, file = "Report/release/SNVs/origin/highdepth_inherited_allele_bymito_bypos.csv.gz")
highdepth_inherited_allele_bymito_bypos <- fread(file = "Report/release/SNVs/origin/highdepth_inherited_allele_bymito_bypos.csv.gz")


highdepth_inherited_noctrl_altperc_bymito_byposmut <- highdepth_inherited_altperc_bymito_byposmut[IsCtrl == "N"]
highdepth_inherited_noctrl_altperc_bymito_byposmut_df <- as.data.frame(highdepth_inherited_noctrl_altperc_bymito_byposmut[, -c(1:19)])
rownames(highdepth_inherited_noctrl_altperc_bymito_byposmut_df) <- highdepth_inherited_noctrl_altperc_bymito_byposmut[, LibraryMitoID]
MitoInfo_df <- as.data.frame(MitoInfo)
rownames(MitoInfo_df) <- MitoInfo_df$LibraryMitoID
MitoInfo_df$MitoID <- factor(MitoInfo_df$MitoID, levels = mitoIDs)
MitoInfo_df$CellID <- as.integer(MitoInfo_df$CellID)
MitoInfo_df$MouseID <- factor(MitoInfo_df$MouseID, levels = sort(MouseInfo[, MouseID]))
pdf("Report/release/SNVs/origin/highdepth_inherited_noctrl_altperc_bymito_byposmut_heatmap.pdf", width = 16, height = 8)
pheatmap(highdepth_inherited_noctrl_altperc_bymito_byposmut_df[with(MitoInfo_df[rownames(highdepth_inherited_noctrl_altperc_bymito_byposmut_df), ], order(MouseID, CellType, CellID, MitoID)), ], cluster_row = FALSE, cluster_col = FALSE, na_col = "#FFFFFF", annotation_row = MitoInfo_df[, c("MouseID", "CellType")], annotation_color = list(MouseID = structure(colorRampPalette(scales::brewer_pal(palette = "Set1")(9))(length(MouseInfo[, MouseID])), names = MouseInfo[, MouseID]), CellType = c("Astrocyte" = Graphics$palette20[5], "Neuron" = Graphics$palette20[4])), fontsize_col = 6, show_rownames = FALSE, angle = 90, border_color = NA, breaks = seq(0, 100, length.out = 255), col = viridis(256))
dev.off()


## 2. somatic SNVs
posmut_somatic <- support_byposmut[nmice == 1, posmut]
pos_somatic <- support_byposmut[nmice == 1, unique(pos)]

highdepth_somatic <- highdepth_qcfltd[as.character(posmut) %in% posmut_somatic]
dim(highdepth_somatic)
## [1] 615323     28
highdepth_somatic <- highdepth_somatic[, -c(27:28)]
fwrite(highdepth_somatic, file = "Report/release/SNVs/origin/highdepth_somatic.csv.gz")
highdepth_somatic <- fread(file = "Report/release/SNVs/origin/highdepth_somatic.csv.gz")

highdepth_somatic_altperc_bymito_byposmut <- highdepth_qcfltd_altperc_bymito_byposmut[, c(1:19, match(posmut_somatic, colnames(highdepth_qcfltd_altperc_bymito_byposmut))), with = FALSE]
fwrite(highdepth_somatic_altperc_bymito_byposmut, file = "Report/release/SNVs/origin/highdepth_somatic_altperc_bymito_byposmut.csv.gz")
highdepth_somatic_altperc_bymito_byposmut <- fread(file = "Report/release/SNVs/origin/highdepth_somatic_altperc_bymito_byposmut.csv.gz")

highdepth_somatic_allele_bymito_bypos <- highdepth_qcfltd_allele_bymito_bypos[, c(1:19, match(pos_somatic, colnames(highdepth_qcfltd_allele_bymito_bypos))), with = FALSE]
fwrite(highdepth_somatic_allele_bymito_bypos, file = "Report/release/SNVs/origin/highdepth_somatic_allele_bymito_bypos.csv.gz")
highdepth_somatic_allele_bymito_bypos <- fread(file = "Report/release/SNVs/origin/highdepth_somatic_allele_bymito_bypos.csv.gz")

highdepth_somatic_noctrl_altperc_bymito_byposmut <- highdepth_somatic_altperc_bymito_byposmut[IsCtrl == "N"]
highdepth_somatic_noctrl_altperc_bymito_byposmut_df <- as.data.frame(highdepth_somatic_noctrl_altperc_bymito_byposmut[, -c(1:19)])
rownames(highdepth_somatic_noctrl_altperc_bymito_byposmut_df) <- highdepth_somatic_noctrl_altperc_bymito_byposmut[, LibraryMitoID]
MitoInfo_df <- as.data.frame(MitoInfo)
rownames(MitoInfo_df) <- MitoInfo_df$LibraryMitoID
MitoInfo_df$MitoID <- factor(MitoInfo_df$MitoID, levels = mitoIDs)
MitoInfo_df$CellID <- as.integer(MitoInfo_df$CellID)
MitoInfo_df$MouseID <- factor(MitoInfo_df$MouseID, levels = sort(MouseInfo[, MouseID]))
pdf("Report/release/SNVs/origin/highdepth_somatic_noctrl_altperc_bymito_byposmut_heatmap.pdf", width = 16, height = 8)
pheatmap(highdepth_somatic_noctrl_altperc_bymito_byposmut_df[with(MitoInfo_df[rownames(highdepth_somatic_noctrl_altperc_bymito_byposmut_df), ], order(MouseID, CellType, CellID, MitoID)), ], cluster_row = FALSE, cluster_col = FALSE, na_col = "#FFFFFF", annotation_row = MitoInfo_df[, c("MouseID", "CellType")], annotation_color = list(MouseID = structure(colorRampPalette(scales::brewer_pal(palette = "Set1")(9))(length(MouseInfo[, MouseID])), names = MouseInfo[, MouseID]), CellType = c("Astrocyte" = Graphics$palette20[5], "Neuron" = Graphics$palette20[4])), fontsize_col = 2, show_rownames = FALSE, angle = 90, border_color = NA, breaks = seq(0, 100, length.out = 255), col = viridis(256))
dev.off()

###########################################################################
## per-site QC stats 
###########################################################################
highdepth_noctrl_hasdata_bymito_bypos <- fread(file = "Report/release/SNVs/QC/highdepth_noctrl_hasdata_bymito_bypos.csv.gz", header = TRUE)
highdepth_noctrl_hasdata_bycell_bypos <- fread(file = "Report/release/SNVs/QC/highdepth_noctrl_hasdata_bycell_bypos.csv.gz", header = TRUE)
highdepth_noctrl_hasdata_bymouse_bypos <- fread(file = "Report/release/SNVs/QC/highdepth_noctrl_hasdata_bymouse_bypos.csv.gz", header = TRUE)
highdepth_noctrl_nmitoshasdata_bycell_bypos <- fread(file = "Report/release/SNVs/QC/highdepth_noctrl_nmitoshasdata_bycell_bypos.csv", header = TRUE)
highdepth_noctrl_ncellshasdata_bymouse_bypos <- fread(file = "Report/release/SNVs/QC/highdepth_noctrl_ncellshasdata_bymouse_bypos.csv", header = TRUE)
highdepth_highaf_noctrl_nvaralleles_bymito_bypos <- fread(file = "Report/release/SNVs/QC/highdepth_highaf_noctrl_nvaralleles_bymito_bypos.csv", header = TRUE)
highdepth_highaf_noctrl_nvaralleles_bycell_bypos <- fread(file = "Report/release/SNVs/QC/highdepth_highaf_noctrl_nvaralleles_bycell_bypos.csv", header = TRUE)
highdepth_highaf_noctrl_nvaralleles_bymouse_bypos <- fread(file = "Report/release/SNVs/QC/highdepth_highaf_noctrl_nvaralleles_bymouse_bypos.csv", header = TRUE)

highdepth_inherited_noctrl_hasdata_bymito_bypos <- data.table(highdepth_noctrl_hasdata_bymito_bypos[, 1:19], highdepth_noctrl_hasdata_bymito_bypos[, match(pos_inherited, names(highdepth_noctrl_hasdata_bymito_bypos)), with = FALSE])
highdepth_inherited_noctrl_hasdata_bycell_bypos <- data.table(highdepth_noctrl_hasdata_bycell_bypos[, 1:7], highdepth_noctrl_hasdata_bycell_bypos[, match(pos_inherited, names(highdepth_noctrl_hasdata_bycell_bypos)), with = FALSE])
highdepth_inherited_noctrl_hasdata_bymouse_bypos <- data.table(highdepth_noctrl_hasdata_bymouse_bypos[, 1:2], highdepth_noctrl_hasdata_bymouse_bypos[, match(pos_inherited, names(highdepth_noctrl_hasdata_bymouse_bypos)), with = FALSE])
highdepth_inherited_noctrl_nmitoshasdata_bycell_bypos <- data.table(highdepth_noctrl_nmitoshasdata_bycell_bypos[, 1:7], highdepth_noctrl_nmitoshasdata_bycell_bypos[, match(pos_inherited, names(highdepth_noctrl_nmitoshasdata_bycell_bypos)), with = FALSE])
highdepth_inherited_noctrl_ncellshasdata_bymouse_bypos <- data.table(highdepth_noctrl_ncellshasdata_bymouse_bypos[, 1:2], highdepth_noctrl_ncellshasdata_bymouse_bypos[, match(pos_inherited, names(highdepth_noctrl_ncellshasdata_bymouse_bypos)), with = FALSE])
highdepth_highaf_inherited_noctrl_nvaralleles_bymito_bypos <- data.table(highdepth_highaf_noctrl_nvaralleles_bymito_bypos[, 1:19], highdepth_highaf_noctrl_nvaralleles_bymito_bypos[, match(pos_inherited, names(highdepth_highaf_noctrl_nvaralleles_bymito_bypos)), with = FALSE])
highdepth_highaf_inherited_noctrl_nvaralleles_bycell_bypos <- data.table(highdepth_highaf_noctrl_nvaralleles_bycell_bypos[, 1:7], highdepth_highaf_noctrl_nvaralleles_bycell_bypos[, match(pos_inherited, names(highdepth_highaf_noctrl_nvaralleles_bycell_bypos)), with = FALSE])
highdepth_highaf_inherited_noctrl_nvaralleles_bymouse_bypos <- data.table(highdepth_highaf_noctrl_nvaralleles_bymouse_bypos[, 1:2], highdepth_highaf_noctrl_nvaralleles_bymouse_bypos[, match(pos_inherited, names(highdepth_highaf_noctrl_nvaralleles_bymouse_bypos)), with = FALSE])

highdepth_highaf_inherited_noctrl_nvaralleles_bymito_bypos[, -c(1:19)][, names(which(apply(.SD, 2, function(x) any(x == 2))))]
##  [1] "1224"  "6430"  "8960"  "9027"  "9461"  "13788" "16094" "16095" "16096"
## [10] "16099" "16108" "16110" "16125" "16133"
highdepth_highaf_inherited_noctrl_nvaralleles_bymito_bypos[, -c(1:19)][, names(which(apply(.SD, 2, function(x) any(x == 3))))]
## [1] "9027"  "16099"

highdepth_highaf_inherited_noctrl_nvaralleles_bycell_bypos[, -c(1:19)][, names(which(apply(.SD, 2, function(x) any(x == 2))))]
## [10] "13786" "13788" "13789" "16093" "16094" "16095" "16096" "16099" "16108"
## [19] "16110" "16120" "16125" "16133"
highdepth_highaf_inherited_noctrl_nvaralleles_bycell_bypos[, -c(1:19)][, names(which(apply(.SD, 2, function(x) any(x == 3))))]
## [1] "9027"  "16099"

highdepth_highaf_inherited_noctrl_nvaralleles_bymouse_bypos[, -c(1:19)][, names(which(apply(.SD, 2, function(x) any(x == 2))))]
## [10] "9463"  "12811" "12831" "13766" "13767" "13778" "13779" "13785" "13786"
## [19] "13788" "13789" "13796" "16090" "16093" "16094" "16096" "16099" "16108"
## [28] "16110" "16120" "16125" "16133" "16140"
highdepth_highaf_inherited_noctrl_nvaralleles_bymouse_bypos[, -c(1:19)][, names(which(apply(.SD, 2, function(x) any(x == 3))))]
## [1] "9027"  "13785" "13786" "16095" "16099"

fwrite(highdepth_inherited_noctrl_hasdata_bymito_bypos, file = "Report/release/SNVs/origin/highdepth_inherited_noctrl_hasdata_bymito_bypos.csv")
highdepth_inherited_noctrl_hasdata_bymito_bypos <- fread(file = "Report/release/SNVs/origin/highdepth_inherited_noctrl_hasdata_bymito_bypos.csv", header = TRUE)
fwrite(highdepth_inherited_noctrl_hasdata_bycell_bypos, file = "Report/release/SNVs/origin/highdepth_inherited_noctrl_hasdata_bycell_bypos.csv")
highdepth_inherited_noctrl_hasdata_bycell_bypos <- fread(file = "Report/release/SNVs/origin/highdepth_inherited_noctrl_hasdata_bycell_bypos.csv", header = TRUE)
fwrite(highdepth_inherited_noctrl_hasdata_bymouse_bypos, file = "Report/release/SNVs/origin/highdepth_inherited_noctrl_hasdata_bymouse_bypos.csv")
highdepth_inherited_noctrl_hasdata_bymouse_bypos <- fread(file = "Report/release/SNVs/origin/highdepth_inherited_noctrl_hasdata_bymouse_bypos.csv", header = TRUE)
fwrite(highdepth_inherited_noctrl_nmitoshasdata_bycell_bypos, file = "Report/release/SNVs/origin/highdepth_inherited_noctrl_nmitoshasdata_bycell_bypos.csv")
highdepth_inherited_noctrl_nmitoshasdata_bycell_bypos <- fread(file = "Report/release/SNVs/origin/highdepth_inherited_noctrl_nmitoshasdata_bycell_bypos.csv", header = TRUE)
fwrite(highdepth_inherited_noctrl_ncellshasdata_bymouse_bypos, file = "Report/release/SNVs/origin/highdepth_inherited_noctrl_ncellshasdata_bymouse_bypos.csv")
highdepth_inherited_noctrl_ncellshasdata_bymouse_bypos <- fread(file = "Report/release/SNVs/origin/highdepth_inherited_noctrl_ncellshasdata_bymouse_bypos.csv", header = TRUE)
fwrite(highdepth_highaf_inherited_noctrl_nvaralleles_bymito_bypos, file = "Report/release/SNVs/origin/highdepth_highaf_inherited_noctrl_nvaralleles_bymito_bypos.csv")
highdepth_highaf_inherited_noctrl_nvaralleles_bymito_bypos <- fread(file = "Report/release/SNVs/origin/highdepth_highaf_inherited_noctrl_nvaralleles_bymito_bypos.csv", header = TRUE)
fwrite(highdepth_highaf_inherited_noctrl_nvaralleles_bycell_bypos, file = "Report/release/SNVs/origin/highdepth_highaf_inherited_noctrl_nvaralleles_bycell_bypos.csv")
highdepth_highaf_inherited_noctrl_nvaralleles_bycell_bypos <- fread(file = "Report/release/SNVs/origin/highdepth_highaf_inherited_noctrl_nvaralleles_bycell_bypos.csv", header = TRUE)
fwrite(highdepth_highaf_inherited_noctrl_nvaralleles_bymouse_bypos, file = "Report/release/SNVs/origin/highdepth_highaf_inherited_noctrl_nvaralleles_bymouse_bypos.csv")
highdepth_highaf_inherited_noctrl_nvaralleles_bymouse_bypos <- fread(file = "Report/release/SNVs/origin/highdepth_highaf_inherited_noctrl_nvaralleles_bymouse_bypos.csv", header = TRUE)

highdepth_highaf_inherited_noctrl_nvaralleles_bymito_bypos_df <- as.data.frame(highdepth_highaf_inherited_noctrl_nvaralleles_bymito_bypos[, -c(1:19)])
highdepth_highaf_inherited_noctrl_nvaralleles_bymito_bypos_df[highdepth_highaf_inherited_noctrl_nvaralleles_bymito_bypos_df == 0] <- NA
rownames(highdepth_highaf_inherited_noctrl_nvaralleles_bymito_bypos_df) <- highdepth_highaf_inherited_noctrl_nvaralleles_bymito_bypos[, LibraryMitoID]
MitoInfo_df <- as.data.frame(MitoInfo)
rownames(MitoInfo_df) <- MitoInfo_df$LibraryMitoID
MitoInfo_df$MitoID <- factor(MitoInfo_df$MitoID, levels = mitoIDs)
MitoInfo_df$CellID <- as.integer(MitoInfo_df$CellID)
highdepth_highaf_inherited_noctrl_nvaralleles_bymito_bypos_df <- highdepth_highaf_inherited_noctrl_nvaralleles_bymito_bypos_df[highdepth_highaf_inherited_noctrl_nvaralleles_bymito_bypos[, order(MouseID, CellType, CellID, MitoID)], ]
pdf("Report/release/SNVs/origin/highdepth_highaf_inherited_noctrl_nvaralleles_bymito_bypos_heatmap.pdf", width = 16, height = 7)
pheatmap(highdepth_highaf_inherited_noctrl_nvaralleles_bymito_bypos_df, cluster_col = FALSE, cluster_row = FALSE, angle = 90, annotation_row = MitoInfo_df[rownames(highdepth_highaf_inherited_noctrl_nvaralleles_bymito_bypos_df), c("MouseID", "CellType", "CellID", "MitoID")], annotation_color = list(MouseID = structure(colorRampPalette(scales::brewer_pal(palette = "Set1")(9))(length(MouseInfo[, MouseID])), names = MouseInfo[, MouseID]), CellType = c("Astrocyte" = Graphics$palette20[5], "Neuron" = Graphics$palette20[4])), fontsize_col = 6, show_rownames = FALSE, na_col = "#FFFFFFFF")
dev.off()

highdepth_highaf_inherited_noctrl_nvaralleles_bycell_bypos_df <- as.data.frame(highdepth_highaf_inherited_noctrl_nvaralleles_bycell_bypos[, -c(1:7)])
highdepth_highaf_inherited_noctrl_nvaralleles_bycell_bypos_df[highdepth_highaf_inherited_noctrl_nvaralleles_bycell_bypos_df == 0] <- NA
rownames(highdepth_highaf_inherited_noctrl_nvaralleles_bycell_bypos_df) <- highdepth_highaf_inherited_noctrl_nvaralleles_bycell_bypos[, CellUID]
highdepth_highaf_inherited_noctrl_nvaralleles_bycell_bypos_df <- highdepth_highaf_inherited_noctrl_nvaralleles_bycell_bypos_df[highdepth_highaf_inherited_noctrl_nvaralleles_bycell_bypos[, order(MouseID, CellType, CellID)], ]
pdf("Report/release/SNVs/origin/highdepth_highaf_inherited_noctrl_nvaralleles_bycell_bypos_heatmap.pdf", width = 16, height = 4.5)
pheatmap(highdepth_highaf_inherited_noctrl_nvaralleles_bycell_bypos_df, cluster_col = FALSE, cluster_row = FALSE, angle = 90, annotation_row = CellInfo_df[rownames(highdepth_highaf_inherited_noctrl_nvaralleles_bycell_bypos_df), c("MouseID", "CellType", "CellID")], annotation_color = list(MouseID = structure(colorRampPalette(scales::brewer_pal(palette = "Set1")(9))(length(MouseInfo[, MouseID])), names = MouseInfo[, MouseID]), CellType = c("Astrocyte" = Graphics$palette20[5], "Neuron" = Graphics$palette20[4])), fontsize_col = 6, show_rownames = FALSE, na_col = "#FFFFFFFF")
dev.off()

highdepth_highaf_inherited_noctrl_nvaralleles_bymouse_bypos_df <- as.data.frame(highdepth_highaf_inherited_noctrl_nvaralleles_bymouse_bypos)
highdepth_highaf_inherited_noctrl_nvaralleles_bymouse_bypos_df[highdepth_highaf_inherited_noctrl_nvaralleles_bymouse_bypos_df == 0] <- NA
rownames(highdepth_highaf_inherited_noctrl_nvaralleles_bymouse_bypos_df) <- highdepth_highaf_inherited_noctrl_nvaralleles_bymouse_bypos[, MouseID]
highdepth_highaf_inherited_noctrl_nvaralleles_bymouse_bypos_df <- highdepth_highaf_inherited_noctrl_nvaralleles_bymouse_bypos_df[, -c(1:2)]
highdepth_highaf_inherited_noctrl_nvaralleles_bymouse_bypos_df <- highdepth_highaf_inherited_noctrl_nvaralleles_bymouse_bypos_df[highdepth_highaf_inherited_noctrl_nvaralleles_bymouse_bypos[, order(MouseID)], ]
pdf("Report/release/SNVs/origin/highdepth_highaf_inherited_noctrl_nvaralleles_bymouse_bypos_heatmap.pdf", width = 16, height = 2)
pheatmap(highdepth_highaf_inherited_noctrl_nvaralleles_bymouse_bypos_df, cluster_col = FALSE, cluster_row = FALSE, angle = 90, annotation_row = MouseInfo_df[rownames(highdepth_highaf_inherited_noctrl_nvaralleles_bymouse_bypos_df), c("RCAID"), drop = FALSE], annotation_color = list(MouseID = structure(colorRampPalette(scales::brewer_pal(palette = "Set1")(9))(length(MouseInfo[, MouseID])), names = MouseInfo[, MouseID]), CellType = c("Astrocyte" = Graphics$palette20[5], "Neuron" = Graphics$palette20[4])), fontsize_col = 6, show_rownames = TRUE, na_col = "#FFFFFFFF")
dev.off()


## 2. somatic SNVs
highdepth_somatic_noctrl_hasdata_bymito_bypos <- data.table(highdepth_noctrl_hasdata_bymito_bypos[, 1:19], highdepth_noctrl_hasdata_bymito_bypos[, match(pos_somatic, names(highdepth_noctrl_hasdata_bymito_bypos)), with = FALSE])
highdepth_somatic_noctrl_hasdata_bycell_bypos <- data.table(highdepth_noctrl_hasdata_bycell_bypos[, 1:7], highdepth_noctrl_hasdata_bycell_bypos[, match(pos_somatic, names(highdepth_noctrl_hasdata_bycell_bypos)), with = FALSE])
highdepth_somatic_noctrl_hasdata_bymouse_bypos <- data.table(highdepth_noctrl_hasdata_bymouse_bypos[, 1:2], highdepth_noctrl_hasdata_bymouse_bypos[, match(pos_somatic, names(highdepth_noctrl_hasdata_bymouse_bypos)), with = FALSE])
highdepth_somatic_noctrl_nmitoshasdata_bycell_bypos <- data.table(highdepth_noctrl_nmitoshasdata_bycell_bypos[, 1:7], highdepth_noctrl_nmitoshasdata_bycell_bypos[, match(pos_somatic, names(highdepth_noctrl_nmitoshasdata_bycell_bypos)), with = FALSE])
highdepth_somatic_noctrl_ncellshasdata_bymouse_bypos <- data.table(highdepth_noctrl_ncellshasdata_bymouse_bypos[, 1:2], highdepth_noctrl_ncellshasdata_bymouse_bypos[, match(pos_somatic, names(highdepth_noctrl_ncellshasdata_bymouse_bypos)), with = FALSE])
highdepth_highaf_somatic_noctrl_nvaralleles_bymito_bypos <- data.table(highdepth_highaf_noctrl_nvaralleles_bymito_bypos[, 1:19], highdepth_highaf_noctrl_nvaralleles_bymito_bypos[, match(pos_somatic, names(highdepth_highaf_noctrl_nvaralleles_bymito_bypos)), with = FALSE])
highdepth_highaf_somatic_noctrl_nvaralleles_bycell_bypos <- data.table(highdepth_highaf_noctrl_nvaralleles_bycell_bypos[, 1:7], highdepth_highaf_noctrl_nvaralleles_bycell_bypos[, match(pos_somatic, names(highdepth_highaf_noctrl_nvaralleles_bycell_bypos)), with = FALSE])
highdepth_highaf_somatic_noctrl_nvaralleles_bymouse_bypos <- data.table(highdepth_highaf_noctrl_nvaralleles_bymouse_bypos[, 1:2], highdepth_highaf_noctrl_nvaralleles_bymouse_bypos[, match(pos_somatic, names(highdepth_highaf_noctrl_nvaralleles_bymouse_bypos)), with = FALSE])

highdepth_highaf_somatic_noctrl_nvaralleles_bymito_bypos[, -c(1:19)][, names(which(apply(.SD, 2, function(x) any(x == 2))))]
##  [1] "1224"  "6430"  "8960"  "9027"  "9461"  "13788" "16094" "16108" "16110"
## [10] "16125" "16133"
highdepth_highaf_somatic_noctrl_nvaralleles_bymito_bypos[, -c(1:19)][, names(which(apply(.SD, 2, function(x) any(x == 3))))]
## [1] "9027"

highdepth_highaf_somatic_noctrl_nvaralleles_bycell_bypos[, -c(1:19)][, names(which(apply(.SD, 2, function(x) any(x == 2))))]
##  [1] "6430"  "6461"  "8960"  "8964"  "8990"  "9027"  "9419"  "9461"  "12831"
## [10] "13778" "13788" "16094" "16108" "16110" "16120" "16125" "16133"
highdepth_highaf_somatic_noctrl_nvaralleles_bycell_bypos[, -c(1:19)][, names(which(apply(.SD, 2, function(x) any(x == 3))))]
## [1] "9027"

highdepth_highaf_somatic_noctrl_nvaralleles_bymouse_bypos[, -c(1:19)][, names(which(apply(.SD, 2, function(x) any(x == 2))))]
##  [1] "1265"  "2659"  "2719"  "2725"  "3048"  "3106"  "3834"  "6430"  "6436"
## [10] "6461"  "7561"  "7660"  "8960"  "8964"  "8990"  "9019"  "9027"  "9393"
## [19] "9419"  "9454"  "9461"  "9463"  "12831" "13766" "13767" "13769" "13778"
## [28] "13788" "13796" "13802" "13804" "13818" "13835" "16059" "16094" "16097"
## [37] "16108" "16110" "16120" "16125" "16133" "16140"
highdepth_highaf_somatic_noctrl_nvaralleles_bymouse_bypos[, -c(1:19)][, names(which(apply(.SD, 2, function(x) any(x == 3))))]
## [1] "9027"

fwrite(highdepth_somatic_noctrl_hasdata_bymito_bypos, file = "Report/release/SNVs/origin/highdepth_somatic_noctrl_hasdata_bymito_bypos.csv")
highdepth_somatic_noctrl_hasdata_bymito_bypos <- fread(file = "Report/release/SNVs/origin/highdepth_somatic_noctrl_hasdata_bymito_bypos.csv", header = TRUE)
fwrite(highdepth_somatic_noctrl_hasdata_bycell_bypos, file = "Report/release/SNVs/origin/highdepth_somatic_noctrl_hasdata_bycell_bypos.csv")
highdepth_somatic_noctrl_hasdata_bycell_bypos <- fread(file = "Report/release/SNVs/origin/highdepth_somatic_noctrl_hasdata_bycell_bypos.csv", header = TRUE)
fwrite(highdepth_somatic_noctrl_hasdata_bymouse_bypos, file = "Report/release/SNVs/origin/highdepth_somatic_noctrl_hasdata_bymouse_bypos.csv")
highdepth_somatic_noctrl_hasdata_bymouse_bypos <- fread(file = "Report/release/SNVs/origin/highdepth_somatic_noctrl_hasdata_bymouse_bypos.csv", header = TRUE)
fwrite(highdepth_somatic_noctrl_nmitoshasdata_bycell_bypos, file = "Report/release/SNVs/origin/highdepth_somatic_noctrl_nmitoshasdata_bycell_bypos.csv")
highdepth_somatic_noctrl_nmitoshasdata_bycell_bypos <- fread(file = "Report/release/SNVs/origin/highdepth_somatic_noctrl_nmitoshasdata_bycell_bypos.csv", header = TRUE)
fwrite(highdepth_somatic_noctrl_ncellshasdata_bymouse_bypos, file = "Report/release/SNVs/origin/highdepth_somatic_noctrl_ncellshasdata_bymouse_bypos.csv")
highdepth_somatic_noctrl_ncellshasdata_bymouse_bypos <- fread(file = "Report/release/SNVs/origin/highdepth_somatic_noctrl_ncellshasdata_bymouse_bypos.csv", header = TRUE)
fwrite(highdepth_highaf_somatic_noctrl_nvaralleles_bymito_bypos, file = "Report/release/SNVs/origin/highdepth_highaf_somatic_noctrl_nvaralleles_bymito_bypos.csv")
highdepth_highaf_somatic_noctrl_nvaralleles_bymito_bypos <- fread(file = "Report/release/SNVs/origin/highdepth_highaf_somatic_noctrl_nvaralleles_bymito_bypos.csv", header = TRUE)
fwrite(highdepth_highaf_somatic_noctrl_nvaralleles_bycell_bypos, file = "Report/release/SNVs/origin/highdepth_highaf_somatic_noctrl_nvaralleles_bycell_bypos.csv")
highdepth_highaf_somatic_noctrl_nvaralleles_bycell_bypos <- fread(file = "Report/release/SNVs/origin/highdepth_highaf_somatic_noctrl_nvaralleles_bycell_bypos.csv", header = TRUE)
fwrite(highdepth_highaf_somatic_noctrl_nvaralleles_bymouse_bypos, file = "Report/release/SNVs/origin/highdepth_highaf_somatic_noctrl_nvaralleles_bymouse_bypos.csv")
highdepth_highaf_somatic_noctrl_nvaralleles_bymouse_bypos <- fread(file = "Report/release/SNVs/origin/highdepth_highaf_somatic_noctrl_nvaralleles_bymouse_bypos.csv", header = TRUE)


highdepth_highaf_somatic_noctrl_nvaralleles_bymito_bypos_df <- as.data.frame(highdepth_highaf_somatic_noctrl_nvaralleles_bymito_bypos[, -c(1:19)])
highdepth_highaf_somatic_noctrl_nvaralleles_bymito_bypos_df[highdepth_highaf_somatic_noctrl_nvaralleles_bymito_bypos_df == 0] <- NA
rownames(highdepth_highaf_somatic_noctrl_nvaralleles_bymito_bypos_df) <- highdepth_highaf_somatic_noctrl_nvaralleles_bymito_bypos[, LibraryMitoID]
MitoInfo_df <- as.data.frame(MitoInfo)
rownames(MitoInfo_df) <- MitoInfo_df$LibraryMitoID
MitoInfo_df$MitoID <- factor(MitoInfo_df$MitoID, levels = mitoIDs)
MitoInfo_df$CellID <- as.integer(MitoInfo_df$CellID)
highdepth_highaf_somatic_noctrl_nvaralleles_bymito_bypos_df <- highdepth_highaf_somatic_noctrl_nvaralleles_bymito_bypos_df[highdepth_highaf_somatic_noctrl_nvaralleles_bymito_bypos[, order(MouseID, CellType, CellID, MitoID)], ]
pdf("Report/release/SNVs/origin/highdepth_highaf_somatic_noctrl_nvaralleles_bymito_bypos_heatmap.pdf", width = 16, height = 7)
pheatmap(highdepth_highaf_somatic_noctrl_nvaralleles_bymito_bypos_df, cluster_col = FALSE, cluster_row = FALSE, angle = 90, annotation_row = MitoInfo_df[rownames(highdepth_highaf_somatic_noctrl_nvaralleles_bymito_bypos_df), c("MouseID", "CellType", "CellID", "MitoID")], annotation_color = list(MouseID = structure(colorRampPalette(scales::brewer_pal(palette = "Set1")(9))(length(MouseInfo[, MouseID])), names = MouseInfo[, MouseID]), CellType = c("Astrocyte" = Graphics$palette20[5], "Neuron" = Graphics$palette20[4])), fontsize_col = 2, show_rownames = FALSE, na_col = "#FFFFFFFF")
dev.off()

highdepth_highaf_somatic_noctrl_nvaralleles_bycell_bypos_df <- as.data.frame(highdepth_highaf_somatic_noctrl_nvaralleles_bycell_bypos[, -c(1:7)])
highdepth_highaf_somatic_noctrl_nvaralleles_bycell_bypos_df[highdepth_highaf_somatic_noctrl_nvaralleles_bycell_bypos_df == 0] <- NA
rownames(highdepth_highaf_somatic_noctrl_nvaralleles_bycell_bypos_df) <- highdepth_highaf_somatic_noctrl_nvaralleles_bycell_bypos[, CellUID]
highdepth_highaf_somatic_noctrl_nvaralleles_bycell_bypos_df <- highdepth_highaf_somatic_noctrl_nvaralleles_bycell_bypos_df[highdepth_highaf_somatic_noctrl_nvaralleles_bycell_bypos[, order(MouseID, CellType, CellID)], ]
pdf("Report/release/SNVs/origin/highdepth_highaf_somatic_noctrl_nvaralleles_bycell_bypos_heatmap.pdf", width = 16, height = 4.5)
pheatmap(highdepth_highaf_somatic_noctrl_nvaralleles_bycell_bypos_df, cluster_col = FALSE, cluster_row = FALSE, angle = 90, annotation_row = CellInfo_df[rownames(highdepth_highaf_somatic_noctrl_nvaralleles_bycell_bypos_df), c("MouseID", "CellType", "CellID")], annotation_color = list(MouseID = structure(colorRampPalette(scales::brewer_pal(palette = "Set1")(9))(length(MouseInfo[, MouseID])), names = MouseInfo[, MouseID]), CellType = c("Astrocyte" = Graphics$palette20[5], "Neuron" = Graphics$palette20[4])), fontsize_col = 2, show_rownames = FALSE, na_col = "#FFFFFFFF")
dev.off()

highdepth_highaf_somatic_noctrl_nvaralleles_bymouse_bypos_df <- as.data.frame(highdepth_highaf_somatic_noctrl_nvaralleles_bymouse_bypos)
highdepth_highaf_somatic_noctrl_nvaralleles_bymouse_bypos_df[highdepth_highaf_somatic_noctrl_nvaralleles_bymouse_bypos_df == 0] <- NA
rownames(highdepth_highaf_somatic_noctrl_nvaralleles_bymouse_bypos_df) <- highdepth_highaf_somatic_noctrl_nvaralleles_bymouse_bypos[, MouseID]
highdepth_highaf_somatic_noctrl_nvaralleles_bymouse_bypos_df <- highdepth_highaf_somatic_noctrl_nvaralleles_bymouse_bypos_df[, -c(1:2)]
highdepth_highaf_somatic_noctrl_nvaralleles_bymouse_bypos_df <- highdepth_highaf_somatic_noctrl_nvaralleles_bymouse_bypos_df[highdepth_highaf_somatic_noctrl_nvaralleles_bymouse_bypos[, order(MouseID)], ]
pdf("Report/release/SNVs/origin/highdepth_highaf_somatic_noctrl_nvaralleles_bymouse_bypos_heatmap.pdf", width = 16, height = 2)
pheatmap(highdepth_highaf_somatic_noctrl_nvaralleles_bymouse_bypos_df, cluster_col = FALSE, cluster_row = FALSE, angle = 90, annotation_row = MouseInfo_df[rownames(highdepth_highaf_somatic_noctrl_nvaralleles_bymouse_bypos_df), c("RCAID"), drop = FALSE], annotation_color = list(MouseID = structure(colorRampPalette(scales::brewer_pal(palette = "Set1")(9))(length(MouseInfo[, MouseID])), names = MouseInfo[, MouseID]), CellType = c("Astrocyte" = Graphics$palette20[5], "Neuron" = Graphics$palette20[4])), fontsize_col = 2, show_rownames = TRUE, na_col = "#FFFFFFFF")
dev.off()

###########################################################################
## per-SNV stats
###########################################################################
highdepth_noctrl_nmitoshasdata_bycell_byposmut <- fread(file = "Report/release/SNVs/QC/highdepth_noctrl_nmitoshasdata_bycell_byposmut.csv.gz")
highdepth_noctrl_ncellshasdata_bymouse_byposmut <- fread(file = "Report/release/SNVs/QC/highdepth_noctrl_ncellshasdata_bymouse_byposmut.csv.gz")
highdepth_highaf_noctrl_hassnv_bymito_byposmut <- fread(file = "Report/release/SNVs/QC/highdepth_highaf_noctrl_hassnv_bymito_byposmut.csv.gz")
highdepth_highaf_noctrl_hassnv_bycell_byposmut <- fread(file = "Report/release/SNVs/QC/highdepth_highaf_noctrl_hassnv_bycell_byposmut.csv.gz")
highdepth_highaf_noctrl_hassnv_bymouse_byposmut <- fread(file = "Report/release/SNVs/QC/highdepth_highaf_noctrl_hassnv_bymouse_byposmut.csv.gz")
highdepth_highaf_noctrl_nmitoshassnv_bycell_byposmut <- fread(file = "Report/release/SNVs/QC/highdepth_highaf_noctrl_nmitoshassnv_bycell_byposmut.csv.gz")
highdepth_highaf_noctrl_ncellshassnv_bymouse_byposmut <- fread(file = "Report/release/SNVs/QC/highdepth_highaf_noctrl_ncellshassnv_bymouse_byposmut.csv.gz")
highdepth_highaf_noctrl_fmitoshassnv_bycell_byposmut <- fread(file = "Report/release/SNVs/QC/highdepth_highaf_noctrl_fmitoshassnv_bycell_byposmut.csv.gz")
highdepth_highaf_noctrl_fcellshassnv_bymouse_byposmut <- fread(file = "Report/release/SNVs/QC/highdepth_highaf_noctrl_fcellshassnv_bymouse_byposmut.csv.gz")
highdepth_highaf_noctrl_fmitoshasdatahassnv_bycell_byposmut <- fread(file = "Report/release/SNVs/QC/highdepth_highaf_noctrl_fmitoshasdatahassnv_bycell_byposmut.csv.gz")
highdepth_highaf_noctrl_fcellshasdatahassnv_bymouse_byposmut <- fread(file = "Report/release/SNVs/QC/highdepth_highaf_noctrl_fcellshasdatahassnv_bymouse_byposmut.csv.gz")

## 1. inherited
highdepth_inherited_noctrl_nmitoshasdata_bycell_byposmut <- data.table(highdepth_noctrl_nmitoshasdata_bycell_byposmut[, 1:7], highdepth_noctrl_nmitoshasdata_bycell_byposmut[, match(posmut_inherited, names(highdepth_noctrl_nmitoshasdata_bycell_byposmut)), with = FALSE])
highdepth_inherited_noctrl_ncellshasdata_bymouse_byposmut <- data.table(highdepth_noctrl_ncellshasdata_bymouse_byposmut[, 1:2], highdepth_noctrl_ncellshasdata_bymouse_byposmut[, match(posmut_inherited, names(highdepth_noctrl_ncellshasdata_bymouse_byposmut)), with = FALSE])
highdepth_highaf_inherited_noctrl_hassnv_bymito_byposmut <- data.table(highdepth_highaf_noctrl_hassnv_bymito_byposmut[, 1:19], highdepth_highaf_noctrl_hassnv_bymito_byposmut[, match(posmut_inherited, names(highdepth_highaf_noctrl_hassnv_bymito_byposmut)), with = FALSE])
highdepth_highaf_inherited_noctrl_hassnv_bycell_byposmut <- data.table(highdepth_highaf_noctrl_hassnv_bycell_byposmut[, 1:7], highdepth_highaf_noctrl_hassnv_bycell_byposmut[, match(posmut_inherited, names(highdepth_highaf_noctrl_hassnv_bycell_byposmut)), with = FALSE])
highdepth_highaf_inherited_noctrl_hassnv_bymouse_byposmut <- data.table(highdepth_highaf_noctrl_hassnv_bymouse_byposmut[, 1:2], highdepth_highaf_noctrl_hassnv_bymouse_byposmut[, match(posmut_inherited, names(highdepth_highaf_noctrl_hassnv_bymouse_byposmut)), with = FALSE])
highdepth_highaf_inherited_noctrl_nmitoshassnv_bycell_byposmut <- data.table(highdepth_highaf_noctrl_nmitoshassnv_bycell_byposmut[, 1:7], highdepth_highaf_noctrl_nmitoshassnv_bycell_byposmut[, match(posmut_inherited, names(highdepth_highaf_noctrl_nmitoshassnv_bycell_byposmut)), with = FALSE])
highdepth_highaf_inherited_noctrl_ncellshassnv_bymouse_byposmut <- data.table(highdepth_highaf_noctrl_ncellshassnv_bymouse_byposmut[, 1:2], highdepth_highaf_noctrl_ncellshassnv_bymouse_byposmut[, match(posmut_inherited, names(highdepth_highaf_noctrl_ncellshassnv_bymouse_byposmut)), with = FALSE])
highdepth_highaf_inherited_noctrl_fmitoshassnv_bycell_byposmut <- data.table(highdepth_highaf_noctrl_fmitoshassnv_bycell_byposmut[, 1:7], highdepth_highaf_noctrl_fmitoshassnv_bycell_byposmut[, match(posmut_inherited, names(highdepth_highaf_noctrl_fmitoshassnv_bycell_byposmut)), with = FALSE])
highdepth_highaf_inherited_noctrl_fcellshassnv_bymouse_byposmut <- data.table(highdepth_highaf_noctrl_fcellshassnv_bymouse_byposmut[, 1:2], highdepth_highaf_noctrl_fcellshassnv_bymouse_byposmut[, match(posmut_inherited, names(highdepth_highaf_noctrl_fcellshassnv_bymouse_byposmut)), with = FALSE])
highdepth_highaf_inherited_noctrl_fmitoshasdatahassnv_bycell_byposmut <- data.table(highdepth_highaf_noctrl_fmitoshasdatahassnv_bycell_byposmut[, 1:7], highdepth_highaf_noctrl_fmitoshasdatahassnv_bycell_byposmut[, match(posmut_inherited, names(highdepth_highaf_noctrl_fmitoshasdatahassnv_bycell_byposmut)), with = FALSE])
highdepth_highaf_inherited_noctrl_fcellshasdatahassnv_bymouse_byposmut <- data.table(highdepth_highaf_noctrl_fcellshasdatahassnv_bymouse_byposmut[, 1:2], highdepth_highaf_noctrl_fcellshasdatahassnv_bymouse_byposmut[, match(posmut_inherited, names(highdepth_highaf_noctrl_fcellshasdatahassnv_bymouse_byposmut)), with = FALSE])

fwrite(highdepth_inherited_noctrl_nmitoshasdata_bycell_byposmut, file = "Report/release/SNVs/origin/highdepth_inherited_noctrl_nmitoshasdata_bycell_byposmut.csv")
highdepth_inherited_noctrl_nmitoshasdata_bycell_byposmut <- fread(file = "Report/release/SNVs/origin/highdepth_inherited_noctrl_nmitoshasdata_bycell_byposmut.csv")
fwrite(highdepth_inherited_noctrl_ncellshasdata_bymouse_byposmut, file = "Report/release/SNVs/origin/highdepth_inherited_noctrl_ncellshasdata_bymouse_byposmut.csv")
highdepth_inherited_noctrl_ncellshasdata_bymouse_byposmut <- fread(file = "Report/release/SNVs/origin/highdepth_inherited_noctrl_ncellshasdata_bymouse_byposmut.csv")
fwrite(highdepth_highaf_inherited_noctrl_hassnv_bymito_byposmut, file = "Report/release/SNVs/origin/highdepth_highaf_inherited_noctrl_hassnv_bymito_byposmut.csv")
highdepth_highaf_inherited_noctrl_hassnv_bymito_byposmut <- fread(file = "Report/release/SNVs/origin/highdepth_highaf_inherited_noctrl_hassnv_bymito_byposmut.csv")
fwrite(highdepth_highaf_inherited_noctrl_hassnv_bycell_byposmut, file = "Report/release/SNVs/origin/highdepth_highaf_inherited_noctrl_hassnv_bycell_byposmut.csv")
highdepth_highaf_inherited_noctrl_hassnv_bycell_byposmut <- fread(file = "Report/release/SNVs/origin/highdepth_highaf_inherited_noctrl_hassnv_bycell_byposmut.csv")
fwrite(highdepth_highaf_inherited_noctrl_hassnv_bymouse_byposmut, file = "Report/release/SNVs/origin/highdepth_highaf_inherited_noctrl_hassnv_bymouse_byposmut.csv")
highdepth_highaf_inherited_noctrl_hassnv_bymouse_byposmut <- fread(file = "Report/release/SNVs/origin/highdepth_highaf_inherited_noctrl_hassnv_bymouse_byposmut.csv")
fwrite(highdepth_highaf_inherited_noctrl_nmitoshassnv_bycell_byposmut, file = "Report/release/SNVs/origin/highdepth_highaf_inherited_noctrl_nmitoshassnv_bycell_byposmut.csv")
highdepth_highaf_inherited_noctrl_nmitoshassnv_bycell_byposmut <- fread(file = "Report/release/SNVs/origin/highdepth_highaf_inherited_noctrl_nmitoshassnv_bycell_byposmut.csv")
fwrite(highdepth_highaf_inherited_noctrl_ncellshassnv_bymouse_byposmut, file = "Report/release/SNVs/origin/highdepth_highaf_inherited_noctrl_ncellshassnv_bymouse_byposmut.csv")
highdepth_highaf_inherited_noctrl_ncellshassnv_bymouse_byposmut <- fread(file = "Report/release/SNVs/origin/highdepth_highaf_inherited_noctrl_ncellshassnv_bymouse_byposmut.csv")
fwrite(highdepth_highaf_inherited_noctrl_fmitoshassnv_bycell_byposmut, file = "Report/release/SNVs/origin/highdepth_highaf_inherited_noctrl_fmitoshassnv_bycell_byposmut.csv")
highdepth_highaf_inherited_noctrl_fmitoshassnv_bycell_byposmut <- fread(file = "Report/release/SNVs/origin/highdepth_highaf_inherited_noctrl_fmitoshassnv_bycell_byposmut.csv")
fwrite(highdepth_highaf_inherited_noctrl_fcellshassnv_bymouse_byposmut, file = "Report/release/SNVs/origin/highdepth_highaf_inherited_noctrl_fcellshassnv_bymouse_byposmut.csv")
highdepth_highaf_inherited_noctrl_fcellshassnv_bymouse_byposmut <- fread(file = "Report/release/SNVs/origin/highdepth_highaf_inherited_noctrl_fcellshassnv_bymouse_byposmut.csv")
fwrite(highdepth_highaf_inherited_noctrl_fmitoshasdatahassnv_bycell_byposmut, file = "Report/release/SNVs/origin/highdepth_highaf_inherited_noctrl_fmitoshasdatahassnv_bycell_byposmut.csv")
highdepth_highaf_inherited_noctrl_fmitoshasdatahassnv_bycell_byposmut <- fread(file = "Report/release/SNVs/origin/highdepth_highaf_inherited_noctrl_fmitoshasdatahassnv_bycell_byposmut.csv")
fwrite(highdepth_highaf_inherited_noctrl_fcellshasdatahassnv_bymouse_byposmut, file = "Report/release/SNVs/origin/highdepth_highaf_inherited_noctrl_fcellshasdatahassnv_bymouse_byposmut.csv")
highdepth_highaf_inherited_noctrl_fcellshasdatahassnv_bymouse_byposmut <- fread(file = "Report/release/SNVs/origin/highdepth_highaf_inherited_noctrl_fcellshasdatahassnv_bymouse_byposmut.csv")


## 2. somatic
highdepth_somatic_noctrl_nmitoshasdata_bycell_byposmut <- data.table(highdepth_noctrl_nmitoshasdata_bycell_byposmut[, 1:7], highdepth_noctrl_nmitoshasdata_bycell_byposmut[, match(posmut_somatic, names(highdepth_noctrl_nmitoshasdata_bycell_byposmut)), with = FALSE])
highdepth_somatic_noctrl_ncellshasdata_bymouse_byposmut <- data.table(highdepth_noctrl_ncellshasdata_bymouse_byposmut[, 1:2], highdepth_noctrl_ncellshasdata_bymouse_byposmut[, match(posmut_somatic, names(highdepth_noctrl_ncellshasdata_bymouse_byposmut)), with = FALSE])
highdepth_highaf_somatic_noctrl_hassnv_bymito_byposmut <- data.table(highdepth_highaf_noctrl_hassnv_bymito_byposmut[, 1:19], highdepth_highaf_noctrl_hassnv_bymito_byposmut[, match(posmut_somatic, names(highdepth_highaf_noctrl_hassnv_bymito_byposmut)), with = FALSE])
highdepth_highaf_somatic_noctrl_hassnv_bycell_byposmut <- data.table(highdepth_highaf_noctrl_hassnv_bycell_byposmut[, 1:7], highdepth_highaf_noctrl_hassnv_bycell_byposmut[, match(posmut_somatic, names(highdepth_highaf_noctrl_hassnv_bycell_byposmut)), with = FALSE])
highdepth_highaf_somatic_noctrl_hassnv_bymouse_byposmut <- data.table(highdepth_highaf_noctrl_hassnv_bymouse_byposmut[, 1:2], highdepth_highaf_noctrl_hassnv_bymouse_byposmut[, match(posmut_somatic, names(highdepth_highaf_noctrl_hassnv_bymouse_byposmut)), with = FALSE])
highdepth_highaf_somatic_noctrl_nmitoshassnv_bycell_byposmut <- data.table(highdepth_highaf_noctrl_nmitoshassnv_bycell_byposmut[, 1:7], highdepth_highaf_noctrl_nmitoshassnv_bycell_byposmut[, match(posmut_somatic, names(highdepth_highaf_noctrl_nmitoshassnv_bycell_byposmut)), with = FALSE])
highdepth_highaf_somatic_noctrl_ncellshassnv_bymouse_byposmut <- data.table(highdepth_highaf_noctrl_ncellshassnv_bymouse_byposmut[, 1:2], highdepth_highaf_noctrl_ncellshassnv_bymouse_byposmut[, match(posmut_somatic, names(highdepth_highaf_noctrl_ncellshassnv_bymouse_byposmut)), with = FALSE])
highdepth_highaf_somatic_noctrl_fmitoshassnv_bycell_byposmut <- data.table(highdepth_highaf_noctrl_fmitoshassnv_bycell_byposmut[, 1:7], highdepth_highaf_noctrl_fmitoshassnv_bycell_byposmut[, match(posmut_somatic, names(highdepth_highaf_noctrl_fmitoshassnv_bycell_byposmut)), with = FALSE])
highdepth_highaf_somatic_noctrl_fcellshassnv_bymouse_byposmut <- data.table(highdepth_highaf_noctrl_fcellshassnv_bymouse_byposmut[, 1:2], highdepth_highaf_noctrl_fcellshassnv_bymouse_byposmut[, match(posmut_somatic, names(highdepth_highaf_noctrl_fcellshassnv_bymouse_byposmut)), with = FALSE])
highdepth_highaf_somatic_noctrl_fmitoshasdatahassnv_bycell_byposmut <- data.table(highdepth_highaf_noctrl_fmitoshasdatahassnv_bycell_byposmut[, 1:7], highdepth_highaf_noctrl_fmitoshasdatahassnv_bycell_byposmut[, match(posmut_somatic, names(highdepth_highaf_noctrl_fmitoshasdatahassnv_bycell_byposmut)), with = FALSE])
highdepth_highaf_somatic_noctrl_fcellshasdatahassnv_bymouse_byposmut <- data.table(highdepth_highaf_noctrl_fcellshasdatahassnv_bymouse_byposmut[, 1:2], highdepth_highaf_noctrl_fcellshasdatahassnv_bymouse_byposmut[, match(posmut_somatic, names(highdepth_highaf_noctrl_fcellshasdatahassnv_bymouse_byposmut)), with = FALSE])

fwrite(highdepth_somatic_noctrl_nmitoshasdata_bycell_byposmut, file = "Report/release/SNVs/origin/highdepth_somatic_noctrl_nmitoshasdata_bycell_byposmut.csv")
highdepth_somatic_noctrl_nmitoshasdata_bycell_byposmut <- fread(file = "Report/release/SNVs/origin/highdepth_somatic_noctrl_nmitoshasdata_bycell_byposmut.csv")
fwrite(highdepth_somatic_noctrl_ncellshasdata_bymouse_byposmut, file = "Report/release/SNVs/origin/highdepth_somatic_noctrl_ncellshasdata_bymouse_byposmut.csv")
highdepth_somatic_noctrl_ncellshasdata_bymouse_byposmut <- fread(file = "Report/release/SNVs/origin/highdepth_somatic_noctrl_ncellshasdata_bymouse_byposmut.csv")
fwrite(highdepth_highaf_somatic_noctrl_hassnv_bymito_byposmut, file = "Report/release/SNVs/origin/highdepth_highaf_somatic_noctrl_hassnv_bymito_byposmut.csv")
highdepth_highaf_somatic_noctrl_hassnv_bymito_byposmut <- fread(file = "Report/release/SNVs/origin/highdepth_highaf_somatic_noctrl_hassnv_bymito_byposmut.csv")
fwrite(highdepth_highaf_somatic_noctrl_hassnv_bycell_byposmut, file = "Report/release/SNVs/origin/highdepth_highaf_somatic_noctrl_hassnv_bycell_byposmut.csv")
highdepth_highaf_somatic_noctrl_hassnv_bycell_byposmut <- fread(file = "Report/release/SNVs/origin/highdepth_highaf_somatic_noctrl_hassnv_bycell_byposmut.csv")
fwrite(highdepth_highaf_somatic_noctrl_hassnv_bymouse_byposmut, file = "Report/release/SNVs/origin/highdepth_highaf_somatic_noctrl_hassnv_bymouse_byposmut.csv")
highdepth_highaf_somatic_noctrl_hassnv_bymouse_byposmut <- fread(file = "Report/release/SNVs/origin/highdepth_highaf_somatic_noctrl_hassnv_bymouse_byposmut.csv")
fwrite(highdepth_highaf_somatic_noctrl_nmitoshassnv_bycell_byposmut, file = "Report/release/SNVs/origin/highdepth_highaf_somatic_noctrl_nmitoshassnv_bycell_byposmut.csv")
highdepth_highaf_somatic_noctrl_nmitoshassnv_bycell_byposmut <- fread(file = "Report/release/SNVs/origin/highdepth_highaf_somatic_noctrl_nmitoshassnv_bycell_byposmut.csv")
fwrite(highdepth_highaf_somatic_noctrl_ncellshassnv_bymouse_byposmut, file = "Report/release/SNVs/origin/highdepth_highaf_somatic_noctrl_ncellshassnv_bymouse_byposmut.csv")
highdepth_highaf_somatic_noctrl_ncellshassnv_bymouse_byposmut <- fread(file = "Report/release/SNVs/origin/highdepth_highaf_somatic_noctrl_ncellshassnv_bymouse_byposmut.csv")
fwrite(highdepth_highaf_somatic_noctrl_fmitoshassnv_bycell_byposmut, file = "Report/release/SNVs/origin/highdepth_highaf_somatic_noctrl_fmitoshassnv_bycell_byposmut.csv")
highdepth_highaf_somatic_noctrl_fmitoshassnv_bycell_byposmut <- fread(file = "Report/release/SNVs/origin/highdepth_highaf_somatic_noctrl_fmitoshassnv_bycell_byposmut.csv")
fwrite(highdepth_highaf_somatic_noctrl_fcellshassnv_bymouse_byposmut, file = "Report/release/SNVs/origin/highdepth_highaf_somatic_noctrl_fcellshassnv_bymouse_byposmut.csv")
highdepth_highaf_somatic_noctrl_fcellshassnv_bymouse_byposmut <- fread(file = "Report/release/SNVs/origin/highdepth_highaf_somatic_noctrl_fcellshassnv_bymouse_byposmut.csv")
fwrite(highdepth_highaf_somatic_noctrl_fmitoshasdatahassnv_bycell_byposmut, file = "Report/release/SNVs/origin/highdepth_highaf_somatic_noctrl_fmitoshasdatahassnv_bycell_byposmut.csv")
highdepth_highaf_somatic_noctrl_fmitoshasdatahassnv_bycell_byposmut <- fread(file = "Report/release/SNVs/origin/highdepth_highaf_somatic_noctrl_fmitoshasdatahassnv_bycell_byposmut.csv")
fwrite(highdepth_highaf_somatic_noctrl_fcellshasdatahassnv_bymouse_byposmut, file = "Report/release/SNVs/origin/highdepth_highaf_somatic_noctrl_fcellshasdatahassnv_bymouse_byposmut.csv")
highdepth_highaf_somatic_noctrl_fcellshasdatahassnv_bymouse_byposmut <- fread(file = "Report/release/SNVs/origin/highdepth_highaf_somatic_noctrl_fcellshasdatahassnv_bymouse_byposmut.csv")
