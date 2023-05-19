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
## Load data
###########################################################################
basedifffreq_cutdemux_q30_unstranded <- fread(file = "Report/release/SNVs/basedifffreq_cutdemux_sub500k_q30_unstranded_9027fixed.csv.gz")
basedifffreq_cutdemux_q30_unstranded[, SNVID := factor(SNVID, levels = snvIDs)]
basedifffreq_cutdemux_q30_unstranded[, MitoID := factor(MitoID, levels = mitoIDs)]
basedifffreq_cutdemux_q30_unstranded[, ExptID := factor(ExptID)]
basedifffreq_cutdemux_q30_unstranded[, CellID := factor(CellID)]
basedifffreq_cutdemux_q30_unstranded[, PlateID := factor(PlateID)]
basedifffreq_cutdemux_q30_unstranded[, ref := toupper(ref)]

###########################################################################
## How many bases have > 50 depth in each mito in each SNV region?
###########################################################################
basedifffreq_cutdemux_q30_unstranded_nbases_highdepth <- dcast(basedifffreq_cutdemux_q30_unstranded[IsCtrl == "N", .(nbases_highdepth = sum(depth >= 50, na.rm = TRUE)), by = c("LibraryMitoID", "SNVID")], LibraryMitoID ~ SNVID, value.var = "nbases_highdepth")
basedifffreq_cutdemux_q30_unstranded_nbases_highdepth_mat <- as.matrix(basedifffreq_cutdemux_q30_unstranded_nbases_highdepth[, -1])
rownames(basedifffreq_cutdemux_q30_unstranded_nbases_highdepth_mat) <- basedifffreq_cutdemux_q30_unstranded_nbases_highdepth[, LibraryMitoID]
colnames(basedifffreq_cutdemux_q30_unstranded_nbases_highdepth_mat) <- sub("SNV", "Region", colnames(basedifffreq_cutdemux_q30_unstranded_nbases_highdepth_mat))
pdf("Report/release/SNVs/QC/nbases_highdepth_heatmap.pdf", height = 9, width = 9)
pheatmap(basedifffreq_cutdemux_q30_unstranded_nbases_highdepth_mat, annotation_row = MitoInfo_df[rownames(basedifffreq_cutdemux_q30_unstranded_nbases_highdepth_mat), c("ExptID", "MitoID", "MouseID", "CellType")], cluster_rows = TRUE, cluster_cols = FALSE, show_rownames = FALSE, angle = 90)
dev.off()

###########################################################################
## How do control samples work in terms of sequencing yield? 
###########################################################################
basedifffreq_cutdemux_q30_unstranded_avgdepth <- basedifffreq_cutdemux_q30_unstranded[, .(avgdepth = sum(depth)/150), by = .(HasPrimers, HasRCA, HasMtDNA, Enzyme, CellType %in% c("Astrocyte", "Neuron"), SNVID, MitoID, LibraryMitoID, LibraryID)]
basedifffreq_cutdemux_q30_unstranded_avgdepth[, CellType := ifelse(CellType, "Cell", "mtDNA")]
basedifffreq_cutdemux_q30_unstranded_avgdepth[, CtrlType := paste(CellType, Enzyme, HasMtDNA, HasRCA, HasPrimers, sep = "_")]
table(basedifffreq_cutdemux_q30_unstranded_avgdepth[, CtrlType])
##  Cell_None_Y_Y_Y mtDNA_None_N_Y_N mtDNA_None_N_Y_Y mtDNA_None_Y_Y_N 
##            18501              144              285              166 
## mtDNA_None_Y_Y_Y 
##              249 
basedifffreq_cutdemux_q30_unstranded_avgdepth[, CtrlType:= unname(c(
    "Cell_None_Y_Y_Y" =  "single-mito mtDNA(+) RCA(+)", 
    "mtDNA_None_N_Y_N" = "pooled mtDNA(-) RCA(-)", 
    "mtDNA_None_N_Y_Y" = "pooled mtDNA(-) RCA(+)", 
    "mtDNA_None_Y_Y_N" = "pooled mtDNA(+) RCA(-)", 
    "mtDNA_None_Y_Y_Y" = "pooled mtDNA(+) RCA(+)"
)[basedifffreq_cutdemux_q30_unstranded_avgdepth[["CtrlType"]]])]
fwrite(basedifffreq_cutdemux_q30_unstranded_avgdepth, file = "Report/release/SNVs/QC/avgdepth_byCtrlType.csv.gz")
basedifffreq_cutdemux_q30_unstranded_avgdepth <- fread(file = "Report/release/SNVs/QC/avgdepth_byCtrlType.csv.gz")
figs <- sapply(snvIDs, function(snv) {
    ggplot(basedifffreq_cutdemux_q30_unstranded_avgdepth[SNVID == snv], aes(x = CtrlType, y = log10(avgdepth))) + geom_boxplot(outlier.size = 0.5) + theme_classic() + xlab("") + ylab("avg. depth (log10)") + coord_flip() + ggtitle(snv)
}, simplify = FALSE)
figs <- gridExtra::marrangeGrob(figs, ncol = 3, nrow = 4, top = "")
ggsave(figs, file = "Report/release/SNVs/QC/avgdepth_byCtrlType_boxplot.pdf", width = 12, height = 7)


png("Report/release/SNVs/QC/depth_bypos_colorbymito.png", width = 300 * 1 * length(snvIDs), height = 300 * 1 * length(mitoIDs), res = 300, type = "cairo")
par(ps = 11, lend = 2, ljoin = 1, bty = "L", mfrow = c(length(mitoIDs), length(snvIDs)), mar = c(2, 2, 1, 0.5), oma = c(0, 0, 0, 0), mgp = c(1.4, 0.5, 0))
for (mitoID in mitoIDs) {
    for (snvID in snvIDs) {
        message(mitoID, " ", snvID)
        X <- basedifffreq_cutdemux_q30_unstranded[IsCtrl == "N" & MitoID == mitoID & SNVID == snvID, .(depth, LibraryMitoID), by = "pos"]
        n <- length(unique(X[, LibraryMitoID]))
        palette(rainbow(n))
        with(X, plot(pos, depth, col = factor(LibraryMitoID), type = "p", xlim = c(snv_info[SNVID == snvID, Start], snv_info[SNVID == snvID, Start + snv_info[, max(Length)]]), log = "y", pch = '.', main = sprintf("%s %s", mitoID, sub("SNV", "Region", snvID)), xlab = "", ylab = "depth"))
        abline(v = snv_info[SNVID == snvID, Start] + 150, col = "gray", lwd = 1, lty = 2)
        abline(v = snv_info[SNVID == snvID, End], col = "black", lwd = 1, lty = 3)
    }
}
dev.off()

summary(basedifffreq_cutdemux_q30_unstranded$depth)
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##       1       8      70   12628     960  494752 
## after replacing 34 low-depth mitos
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.  
##       1       9      74   12679    1003  494752
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##       1       9      74   12560     995  494752 
## after removing SNVs inside PCR forward primers
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
##       1       9      74   12560     995  494752
## after fixing 9027
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
##       1       9      74   12582     999 2595925

###########################################################################
## more detailed depth summary
###########################################################################
basedifffreq_cutdemux_q30_unstranded_depthstats <- basedifffreq_cutdemux_q30_unstranded[, .SD, .SDcols = c("LibraryMitoID", "depth"), by = c("LibraryMitoID", "SNVID")]
basedifffreq_cutdemux_q30_unstranded_depthstats <- basedifffreq_cutdemux_q30_unstranded[, .SD, .SDcols = c("depth"), by = c("LibraryMitoID", "SNVID")]

depth_th <- 50
basedifffreq_cutdemux_q30_unstranded_depthstats <- basedifffreq_cutdemux_q30_unstranded_depthstats[, list(
    sum = sum(depth), 
    max = max(depth), 
    min = min(depth), 
    mean = mean(depth),
    median = median(depth), 
    sd = sd(depth), 
    n_highdepth = sum(depth >= depth_th)
    ), 
by = c("LibraryMitoID", "SNVID")]
basedifffreq_cutdemux_q30_unstranded_depthstats <- MitoInfo[basedifffreq_cutdemux_q30_unstranded_depthstats, on = "LibraryMitoID"]
basedifffreq_cutdemux_q30_unstranded_depthstats[, SNVID := factor(SNVID, levels = snvIDs)]

fwrite(basedifffreq_cutdemux_q30_unstranded_depthstats, file = "Report/release/SNVs/QC/depthstats.csv.gz")
basedifffreq_cutdemux_q30_unstranded_depthstats <- fread("Report/release/SNVs/QC/depthstats.csv.gz")

basedifffreq_cutdemux_q30_unstranded_depthstats_nbaseshighdepth <- dcast(basedifffreq_cutdemux_q30_unstranded_depthstats, LibraryMitoID ~ SNVID, value.var = "n_highdepth", fill = 0)
basedifffreq_cutdemux_q30_unstranded_depthstats_nbaseshighdepth[, c("sum", "max", "min", "mean", "median", "sd") := list(rowSums(.SD), apply(.SD, 1, max), apply(.SD, 1, min), apply(.SD, 1, mean), apply(.SD, 1, median), apply(.SD, 1, sd)), .SDcols = -"LibraryMitoID"]
basedifffreq_cutdemux_q30_unstranded_depthstats_nbaseshighdepth <- MitoInfo[basedifffreq_cutdemux_q30_unstranded_depthstats_nbaseshighdepth, on = "LibraryMitoID"]

fwrite(basedifffreq_cutdemux_q30_unstranded_depthstats_nbaseshighdepth, file = "Report/release/SNVs/QC/highdepthstats.csv.gz")
basedifffreq_cutdemux_q30_unstranded_highdepthstats <- fread("Report/release/SNVs/QC/highdepthstats.csv.gz")

###########################################################################
## binary table: whether nonmissing data in mito/cell/mouse X at site Y
###########################################################################
basediffperc_cutdemux_q30_unstranded_highdepth <- fread("Report/release/SNVs/filter/basediffperc_cutdemux_sub500k_q30_unstranded_highdepth_qcfltd.csv.gz")
basediffperc_cutdemux_q30_unstranded_highdepth_noctrl_hasdata_bymito_bypos <- dcast.data.table(basediffperc_cutdemux_q30_unstranded_highdepth[IsCtrl == "N"], LibraryMitoID ~ pos, value.var = "depth", fun.aggregate = function(x) ifelse(any(x), 1, 0))
basediffperc_cutdemux_q30_unstranded_highdepth_noctrl_hasdata_bymito_bypos <- merge.data.table(MitoInfo[IsCtrl == "N"], basediffperc_cutdemux_q30_unstranded_highdepth_noctrl_hasdata_bymito_bypos, by = "LibraryMitoID", all.x = TRUE)
fwrite(basediffperc_cutdemux_q30_unstranded_highdepth_noctrl_hasdata_bymito_bypos, file = "Report/release/SNVs/QC/highdepth_noctrl_hasdata_bymito_bypos.csv.gz")
basediffperc_cutdemux_q30_unstranded_highdepth_noctrl_hasdata_bymito_bypos <- fread(file = "Report/release/SNVs/QC/highdepth_noctrl_hasdata_bymito_bypos.csv.gz", header = TRUE)

basediffperc_cutdemux_q30_unstranded_highdepth_noctrl_hasdata_bycell_bypos <- dcast.data.table(basediffperc_cutdemux_q30_unstranded_highdepth[IsCtrl == "N"], CellUID ~ pos, value.var = "depth", fun.aggregate = function(x) ifelse(any(x), 1, 0))
basediffperc_cutdemux_q30_unstranded_highdepth_noctrl_hasdata_bycell_bypos <- merge.data.table(CellInfo, basediffperc_cutdemux_q30_unstranded_highdepth_noctrl_hasdata_bycell_bypos, by = "CellUID", all.x = TRUE)
fwrite(basediffperc_cutdemux_q30_unstranded_highdepth_noctrl_hasdata_bycell_bypos, file = "Report/release/SNVs/QC/highdepth_noctrl_hasdata_bycell_bypos.csv.gz")
basediffperc_cutdemux_q30_unstranded_highdepth_noctrl_hasdata_bycell_bypos <- fread(file = "Report/release/SNVs/QC/highdepth_noctrl_hasdata_bycell_bypos.csv.gz", header = TRUE)

basediffperc_cutdemux_q30_unstranded_highdepth_noctrl_hasdata_bymouse_bypos <- dcast.data.table(basediffperc_cutdemux_q30_unstranded_highdepth[IsCtrl == "N"], MouseID ~ pos, value.var = "depth", fun.aggregate = function(x) ifelse(any(x), 1, 0))
basediffperc_cutdemux_q30_unstranded_highdepth_noctrl_hasdata_bymouse_bypos <- merge.data.table(MouseInfo, basediffperc_cutdemux_q30_unstranded_highdepth_noctrl_hasdata_bymouse_bypos, by = "MouseID", all.x = TRUE)
fwrite(basediffperc_cutdemux_q30_unstranded_highdepth_noctrl_hasdata_bymouse_bypos, file = "Report/release/SNVs/QC/highdepth_noctrl_hasdata_bymouse_bypos.csv.gz")
basediffperc_cutdemux_q30_unstranded_highdepth_noctrl_hasdata_bymouse_bypos <- fread(file = "Report/release/SNVs/QC/highdepth_noctrl_hasdata_bymouse_bypos.csv.gz", header = TRUE)

###########################################################################
## For each mito/cell/mouse, how many sites have non-missing data (the sites
## are restricted to the SNV sites)?
###########################################################################
basediffperc_cutdemux_q30_unstranded_highdepth_noctrl_nposhasdata_bymito <- data.table(basediffperc_cutdemux_q30_unstranded_highdepth_noctrl_hasdata_bymito_bypos[, 1:19], nposhasdata = rowSums(as.matrix(basediffperc_cutdemux_q30_unstranded_highdepth_noctrl_hasdata_bymito_bypos[, -c(1:19)])))
fwrite(basediffperc_cutdemux_q30_unstranded_highdepth_noctrl_nposhasdata_bymito, file = "Report/release/SNVs/QC/highdepth_noctrl_nposhasdata_bymito.csv")
basediffperc_cutdemux_q30_unstranded_highdepth_noctrl_nposhasdata_bymito <- fread(file = "Report/release/SNVs/QC/highdepth_noctrl_nposhasdata_bymito.csv")

basediffperc_cutdemux_q30_unstranded_highdepth_noctrl_nposhasdata_bycell <- data.table(basediffperc_cutdemux_q30_unstranded_highdepth_noctrl_hasdata_bycell_bypos[, 1:7], nposhasdata = rowSums(as.matrix(basediffperc_cutdemux_q30_unstranded_highdepth_noctrl_hasdata_bycell_bypos[, -c(1:7)])))
fwrite(basediffperc_cutdemux_q30_unstranded_highdepth_noctrl_nposhasdata_bycell, file = "Report/release/SNVs/QC/highdepth_noctrl_nposhasdata_bycell.csv")
basediffperc_cutdemux_q30_unstranded_highdepth_noctrl_nposhasdata_bycell <- fread(file = "Report/release/SNVs/QC/highdepth_noctrl_nposhasdata_bycell.csv")

basediffperc_cutdemux_q30_unstranded_highdepth_noctrl_nposhasdata_bymouse <- data.table(basediffperc_cutdemux_q30_unstranded_highdepth_noctrl_hasdata_bymouse_bypos[, 1:2], nposhasdata = rowSums(as.matrix(basediffperc_cutdemux_q30_unstranded_highdepth_noctrl_hasdata_bymouse_bypos[, -c(1:2)])))
fwrite(basediffperc_cutdemux_q30_unstranded_highdepth_noctrl_nposhasdata_bymouse, file = "Report/release/SNVs/QC/highdepth_noctrl_nposhasdata_bymouse.csv")
basediffperc_cutdemux_q30_unstranded_highdepth_noctrl_nposhasdata_bymouse <- fread(file = "Report/release/SNVs/QC/highdepth_noctrl_nposhasdata_bymouse.csv")

###########################################################################
## For each site, how many mitos/cells/mice have non-missing data?
###########################################################################
## make sure positions in column names are sorted. 
all(diff(as.integer(names(basediffperc_cutdemux_q30_unstranded_highdepth_noctrl_hasdata_bymito_bypos[, -c(1:19)]))) > 0)
## [1] TRUE
length(as.integer(names(basediffperc_cutdemux_q30_unstranded_highdepth_noctrl_hasdata_bymito_bypos[, -c(1:19)])))
## 838

basediffperc_cutdemux_q30_unstranded_highdepth_noctrl_hasdata_bymito_bypos <- fread(file = "Report/release/SNVs/QC/highdepth_noctrl_hasdata_bymito_bypos.csv.gz", header = TRUE)
basediffperc_cutdemux_q30_unstranded_highdepth_noctrl_hasdata_bycell_bypos <- fread(file = "Report/release/SNVs/QC/highdepth_noctrl_hasdata_bycell_bypos.csv.gz", header = TRUE)
basediffperc_cutdemux_q30_unstranded_highdepth_noctrl_hasdata_bymouse_bypos <- fread(file = "Report/release/SNVs/QC/highdepth_noctrl_hasdata_bymouse_bypos.csv.gz", header = TRUE)
basediffperc_cutdemux_q30_unstranded_highdepth_noctrl_nsampleshasdata_bypos <- data.table(
    pos = as.integer(names(basediffperc_cutdemux_q30_unstranded_highdepth_noctrl_hasdata_bymito_bypos[, -c(1:19)])), 
    nmitoshasdata = colSums(basediffperc_cutdemux_q30_unstranded_highdepth_noctrl_hasdata_bymito_bypos[, names(basediffperc_cutdemux_q30_unstranded_highdepth_noctrl_hasdata_bymito_bypos[, -c(1:19)]), with = FALSE]), 
    ncellshasdata = colSums(basediffperc_cutdemux_q30_unstranded_highdepth_noctrl_hasdata_bycell_bypos[, names(basediffperc_cutdemux_q30_unstranded_highdepth_noctrl_hasdata_bymito_bypos[, -c(1:19)]), with = FALSE]), 
    nmicehasdata = colSums(basediffperc_cutdemux_q30_unstranded_highdepth_noctrl_hasdata_bymouse_bypos[, names(basediffperc_cutdemux_q30_unstranded_highdepth_noctrl_hasdata_bymito_bypos[, -c(1:19)]), with = FALSE])
)
fwrite(basediffperc_cutdemux_q30_unstranded_highdepth_noctrl_nsampleshasdata_bypos, file = "Report/release/SNVs/QC/highdepth_noctrl_nsampleshasdata_bypos.csv")
basediffperc_cutdemux_q30_unstranded_highdepth_noctrl_nsampleshasdata_bypos <- fread(file = "Report/release/SNVs/QC/highdepth_noctrl_nsampleshasdata_bypos.csv", header = TRUE)

pdf("Report/release/SNVs/QC/highdepth_noctrl_nsampleshasdata_bypos_barplot.pdf", width = 12, height = 4)
ggplot(basediffperc_cutdemux_q30_unstranded_highdepth_noctrl_nsampleshasdata_bypos, aes(y = nmitoshasdata, x = pos)) + theme_classic(base_size = 16) + theme(axis.line.x = element_blank(), axis.ticks.x = element_blank(), panel.grid.major.y = element_line(linetype = 2)) + geom_bar(stat = "identity", color = "grey", fill = NA) + ylab("# mitos without missing data") + xlab("")
ggplot(basediffperc_cutdemux_q30_unstranded_highdepth_noctrl_nsampleshasdata_bypos, aes(y = ncellshasdata, x = pos)) + theme_classic(base_size = 16) + theme(axis.line.x = element_blank(), axis.ticks.x = element_blank(), panel.grid.major.y = element_line(linetype = 2)) + geom_bar(stat = "identity", color = "grey", fill = NA) + ylab("# cells without missing data") + xlab("") 
ggplot(basediffperc_cutdemux_q30_unstranded_highdepth_noctrl_nsampleshasdata_bypos, aes(y = nmicehasdata, x = pos)) + theme_classic(base_size = 16) + theme(axis.line.x = element_blank(), axis.ticks.x = element_blank(), panel.grid.major.y = element_line(linetype = 2)) + geom_bar(stat = "identity", color = "grey", fill = NA) + ylab("# mice without missing data") + xlab("") 
dev.off()

###########################################################################
## how many non-missing-data mitos per cell and how many non-missing cells per mouse
## at each site
###########################################################################
basediffperc_cutdemux_q30_unstranded_highdepth_noctrl_nmitoshasdata_bycell_bypos <- dcast.data.table(basediffperc_cutdemux_q30_unstranded_highdepth[IsCtrl == "N"], CellUID ~ pos, value.var = "LibraryMitoID", fun.aggregate = uniqueN)
basediffperc_cutdemux_q30_unstranded_highdepth_noctrl_nmitoshasdata_bycell_bypos <- CellInfo[basediffperc_cutdemux_q30_unstranded_highdepth_noctrl_nmitoshasdata_bycell_bypos, on = "CellUID"]
fwrite(basediffperc_cutdemux_q30_unstranded_highdepth_noctrl_nmitoshasdata_bycell_bypos, file = "Report/release/SNVs/QC/highdepth_noctrl_nmitoshasdata_bycell_bypos.csv")
basediffperc_cutdemux_q30_unstranded_highdepth_noctrl_nmitoshasdata_bycell_bypos <- fread(file = "Report/release/SNVs/QC/highdepth_noctrl_nmitoshasdata_bycell_bypos.csv", header = TRUE)

basediffperc_cutdemux_q30_unstranded_highdepth_noctrl_ncellshasdata_bymouse_bypos <- dcast.data.table(basediffperc_cutdemux_q30_unstranded_highdepth[IsCtrl == "N"], MouseID ~ pos, value.var = "CellUID", fun.aggregate = uniqueN)
basediffperc_cutdemux_q30_unstranded_highdepth_noctrl_ncellshasdata_bymouse_bypos <- MouseInfo[basediffperc_cutdemux_q30_unstranded_highdepth_noctrl_ncellshasdata_bymouse_bypos, on = "MouseID"]
fwrite(basediffperc_cutdemux_q30_unstranded_highdepth_noctrl_ncellshasdata_bymouse_bypos, file = "Report/release/SNVs/QC/highdepth_noctrl_ncellshasdata_bymouse_bypos.csv")
basediffperc_cutdemux_q30_unstranded_highdepth_noctrl_ncellshasdata_bymouse_bypos <- fread(file = "Report/release/SNVs/QC/highdepth_noctrl_ncellshasdata_bymouse_bypos.csv", header = TRUE)

basediffperc_cutdemux_q30_unstranded_highdepth_noctrl_nmitoshasdata_bycell_byposmut_df <- as.data.frame(basediffperc_cutdemux_q30_unstranded_highdepth_noctrl_nmitoshasdata_bycell_byposmut)
rownames(basediffperc_cutdemux_q30_unstranded_highdepth_noctrl_nmitoshasdata_bycell_byposmut_df) <- basediffperc_cutdemux_q30_unstranded_highdepth_noctrl_nmitoshasdata_bycell_byposmut_df[, 1]
basediffperc_cutdemux_q30_unstranded_highdepth_noctrl_nmitoshasdata_bycell_byposmut_df <- basediffperc_cutdemux_q30_unstranded_highdepth_noctrl_nmitoshasdata_bycell_byposmut_df[, -c(1:7)]

pdf("Report/release/SNVs/QC/highdepth_noctrl_nmitoshasdata_bycell_byposmut_pheatmap.pdf", width = 18, height = 7)
pheatmap(basediffperc_cutdemux_q30_unstranded_highdepth_noctrl_nmitoshasdata_bycell_byposmut_df, cluster_col = FALSE, cluster_row = TRUE, angle = 90, annotation_row = CellInfo_df[rownames(basediffperc_cutdemux_q30_unstranded_highdepth_noctrl_nmitoshasdata_bycell_byposmut_df), c("MouseID", "CellType", "CellID")], annotation_color = list(MouseID = structure(colorRampPalette(scales::brewer_pal(palette = "Set1")(9))(length(MouseInfo[, MouseID])), names = MouseInfo[, MouseID]), CellType = c("Astrocyte" = Graphics$palette20[5], "Neuron" = Graphics$palette20[4])), fontsize_col = 1, fontsize_row = 4, show_rownames = TRUE, border_color = NA)
dev.off()

basediffperc_cutdemux_q30_unstranded_highdepth_noctrl_ncellshasdata_bymouse_byposmut_df <- as.data.frame(basediffperc_cutdemux_q30_unstranded_highdepth_noctrl_ncellshasdata_bymouse_byposmut)
rownames(basediffperc_cutdemux_q30_unstranded_highdepth_noctrl_ncellshasdata_bymouse_byposmut_df) <- basediffperc_cutdemux_q30_unstranded_highdepth_noctrl_ncellshasdata_bymouse_byposmut_df[, 1]
basediffperc_cutdemux_q30_unstranded_highdepth_noctrl_ncellshasdata_bymouse_byposmut_df <- basediffperc_cutdemux_q30_unstranded_highdepth_noctrl_ncellshasdata_bymouse_byposmut_df[, -c(1:2)]

pdf("Report/release/SNVs/QC/highdepth_noctrl_ncellshasdata_bymouse_byposmut_pheatmap.pdf", width = 18, height = 2)
pheatmap(basediffperc_cutdemux_q30_unstranded_highdepth_noctrl_ncellshasdata_bymouse_byposmut_df, cluster_col = FALSE, cluster_row = TRUE, angle = 90, annotation_row = MouseInfo_df[rownames(basediffperc_cutdemux_q30_unstranded_highdepth_noctrl_ncellshasdata_bymouse_byposmut_df), c("RCAID"), drop = FALSE], fontsize_col = 1, fontsize_row = 8, show_rownames = TRUE, border_color = NA)
dev.off()

###########################################################################
## Load VAF filtered data, compute SNV stats per SNV site
###########################################################################
basediffperc_cutdemux_q30_unstranded_highdepth_highaf <- fread("Report/release/SNVs/filter/basediffperc_cutdemux_sub500k_q30_unstranded_highdepth_highaf_qcfltd.csv.gz")
support_byposmut <- fread("Report/release/SNVs/filter/basediffperc_cutdemux_sub500k_q30_unstranded_highdepth_highaf_qcfltd_support_byposmut.csv")
dim(basediffperc_cutdemux_q30_unstranded_highdepth_highaf)
## [1] 6705   29
## [1] 6354   29
## [1] 7175   29 after fixing 9027
dim(basediffperc_cutdemux_q30_unstranded_highdepth_highaf[IsCtrl == "N"])
## [1] 6314   29
## [1] 5988   29
## [1] 6763   29 after fixing 9027

basediffperc_cutdemux_q30_unstranded_highdepth_highaf[IsCtrl == "N", list(N = uniqueN(pos)), by = "LibraryMitoID"][, summary(N)]
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
##   1.000   2.000   3.000   4.452   6.000  17.000
basediffperc_cutdemux_q30_unstranded_highdepth_highaf[IsCtrl == "N", list(N = uniqueN(pos)), by = "LibraryMitoID"][, sd(N)]
## [1] 3.464509
basediffperc_cutdemux_q30_unstranded_highdepth_highaf[IsCtrl == "N", list(N = uniqueN(pos)), by = "CellUID"][, summary(N)]
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
##    2.00   15.00   24.50   27.46   36.75   97.00
basediffperc_cutdemux_q30_unstranded_highdepth_highaf[IsCtrl == "N", list(N = uniqueN(pos)), by = "CellUID"][, sd(N)]
## [1] 17.9232
basediffperc_cutdemux_q30_unstranded_highdepth_highaf[IsCtrl == "N", list(N = uniqueN(pos)), by = "MouseID"][, summary(N)]
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
##    26.0    59.0   151.0   135.8   165.0   263.0
basediffperc_cutdemux_q30_unstranded_highdepth_highaf[IsCtrl == "N", list(N = uniqueN(pos)), by = "MouseID"][, sd(N)]
## [1] 76.49309

pdf("Report/release/SNVs/QC/highdepth_highaf_noctrl_vaf_hist.pdf", width = 6, height = 6)
par(ps = 16, lend = 2, ljoin = 1, bty = "L", mfrow = c(1, 1), mar = c(2.5, 2.5, 0.5, 0.5), oma = c(0, 0, 0, 0), mgp = c(1.5, 0.5, 0))
basediffperc_cutdemux_q30_unstranded_highdepth_highaf[IsCtrl == "N" , hist(100 - `=`, xlab = "VAF (%)", ylab = "# non-control mito-sites", main = "", freq = TRUE, col = "gray", border = NA)]
dev.off()

basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl <- basediffperc_cutdemux_q30_unstranded_highdepth_highaf[IsCtrl == "N"][order(MouseID, CellType, CellID, LibraryMitoID)]
n_samples <- nrow(basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl)
basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl[, i := seq(n_samples)]
fig <- ggplot(basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl, aes(y = i, x = pos)) + geom_tile(color = "red") + theme_classic(base_size = 11) + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), plot.margin = unit(c(0.3, 0.05, 0.05, 0.05), "in")) + xlab("") + ylab("single mito") + coord_cartesian(clip = "off") + 
    annotate(geom = "rect", xmin = 16800, xmax = 16900, ymin = 1:n_samples, ymax = 1:n_samples, col = colorRampPalette(scales::brewer_pal(palette = "Set1")(9))(length(unique(basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl$MouseID)))[factor(basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl$MouseID)]) + 
    annotate(geom = "rect", xmin = 17000, xmax = 17100, ymin = 1:n_samples, ymax = 1:n_samples, col = c("Astrocyte" =  Graphics$palette20[5], "Neuron" =  Graphics$palette20[4])[basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl$CellType]) + 
    annotate(geom = "rect", xmin = 17200, xmax = 17300, ymin = 1:n_samples, ymax = 1:n_samples, col = colorRampPalette(scales::brewer_pal(palette = "Set3")(9))(length(unique(basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl$CellID)))[factor(basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl$CellID)]) + 
    annotate(geom = "rect", xmin = 17400, xmax = 17500, ymin = 1:n_samples, ymax = 1:n_samples, col = colorRampPalette(rainbow(9))(length(unique(basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl$LibraryID)))[factor(basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl$LibraryID)]) + 
    annotate(geom = "rect", xmin = 17600, xmax = 17700, ymin = 1:n_samples, ymax = 1:n_samples, col = colorRampPalette(rainbow(9))(length(unique(basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl$MitoID)))[factor(basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl$MitoID)]) + 
    annotate(geom = "rect", xmin = snv_info$Start, xmax = snv_info$End, ymin = n_samples + 30, ymax = n_samples + 15, fill = "#407FA0") + 
    annotate(geom = "text", x = c(16800, 17000, 17200, 17400, 17600), y = 0, label = c("MouseID", "CellType", "CellID", "LibraryID", "MitoID"), size = 2, angle = 60, vjust = 1, hjust = 1) + 
    annotate(geom = "text", x = 0.5 * (snv_info$Start + snv_info$End), y = n_samples + 50, label = snv_info$SNVID, size = 2, angle = 90, hjust = 0)
ggsave(fig, filename = "Report/release/SNVs/QC/highdepth_highaf_noctrl_hassnv_per_pos_per_mito.pdf", width = 10, height = 4)

###########################################################################
## Circos plot for the support per base position
###########################################################################
chrmgenes <- fread("Data/mm10.mito/chrM.gtf")
chrmgenes[, symbol := sapply(regmatches(V9, regexec('gene_name "(.*)"', V9, perl = TRUE)), '[', 2)]
chrmgenes <- chrmgenes[, list(symbol, V4, V5, V7)]
colpal <- scales::hue_pal()(3)
pdf("Report/release/SNVs/QC/highdepth_highaf_noctrl_support_circos.pdf", height = 7, width = 7)
circos.clear()
circos.par(start.degree = 90)
circos.initialize(sectors = rep("chrM", 16299), x = 1:16299, ring = TRUE)
circos.track(sectors = rep("chrM", nrow(chrmgenes)), ylim = c(0, 1), panel.fun = function(x, y) { 
    circos.rect(xleft = chrmgenes[, V4], xright = chrmgenes[, V5], ytop = 1, ybottom = 0.8, border = NA, col = ifelse(chrmgenes[, V7] == "+",adjustcolor("grey50", alpha.f = 0.4), adjustcolor("grey50", alpha.f = 0.4)))
    circos.text(x = chrmgenes[, 0.5 * (V4 + V5)], y = 0.25, labels = chrmgenes[, symbol], facing = "clockwise", col = ifelse(chrmgenes[, V7] == "+", Graphics$palette20[15], Graphics$palette20[19]), cex = 0.8)}, 
bg.border = NA, track.height = 0.15)
circos.track(sectors = rep("chrM", nrow(support_byposmut)), x = support_byposmut[, pos], y = support_byposmut[, log10(nmitos)], panel.fun = function(x, y) { circos.barplot(pos = x, value = y, border = colpal[1], col = NA, lwd = 1); circos.lines(y = rep(log10(1), 16299), x = 1:16299, col = "gray50", lty = 1, lwd = 1) }, bg.border = NA)
circos.track(sectors = rep("chrM", nrow(support_byposmut)), x = support_byposmut[, pos], y = support_byposmut[, ncells], panel.fun = function(x, y) { circos.barplot(pos = x, value = y, border = colpal[2], col = NA, lwd = 1); circos.lines(y = rep(0, 16299), x = 1:16299, col = "#FFFFFFFF", lty = 1, lwd = 1) }, bg.border = NA)
circos.track(sectors = rep("chrM", nrow(support_byposmut)), x = support_byposmut[, pos], y = support_byposmut[, nmice], panel.fun = function(x, y) { circos.barplot(pos = x, value = y, border = colpal[3], col = NA, lwd = 1); circos.lines(y = rep(0, 16299), x = 1:16299, col = "#FFFFFFFF", lty = 1, lwd = 1) }, bg.border = NA)
circos.xaxis(major.at = seq(0, 16299, by = 1000), h = 60)
circos.yaxis(at = c(0, 1, 2, 3), side = "right", track.index = 2, labels = c("1", "10", "100", "1000"))
circos.yaxis(at = c(0, 50, 100), side = "right", track.index = 3, labels = c("0", "50", "100"))
circos.yaxis(at = c(0, 5, 10), side = "right", track.index = 4, labels = c("0", "5", "10"))
dev.off()

###########################################################################
## # of SNV sites per mito
###########################################################################
basediffperc_cutdemux_q30_unstranded_highdepth_highaf <- fread("Report/release/SNVs/filter/basediffperc_cutdemux_sub500k_q30_unstranded_highdepth_highaf_qcfltd.csv.gz", header = TRUE)
basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl_npos_bymito <- basediffperc_cutdemux_q30_unstranded_highdepth_highaf[IsCtrl == "N", .(npos = uniqueN(pos)), by = "LibraryMitoID"]
basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl_npos_bymito <- merge.data.table(MitoInfo[IsCtrl == "N"], basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl_npos_bymito, by = "LibraryMitoID", all.x = TRUE)
basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl_npos_bymito[, MitoID := factor(MitoID, levels = mitoIDs)]
basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl_npos_bymito[, CellID := factor(CellID)]
basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl_npos_bymito[, ExptID := factor(ExptID)]
basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl_npos_bymito[, npos := ifelse(is.na(npos), 0, npos)]
fwrite(basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl_npos_bymito, file = "Report/release/SNVs/QC/highdepth_highaf_noctrl_npos_bymito.csv")
basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl_npos_bymito <- fread(file = "Report/release/SNVs/QC/highdepth_highaf_noctrl_npos_bymito.csv")

basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl_npos_bycell <- basediffperc_cutdemux_q30_unstranded_highdepth_highaf[IsCtrl == "N", .(npos = uniqueN(pos)), by = "CellUID"]
basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl_npos_bycell <- merge.data.table(CellInfo, basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl_npos_bycell, by = "CellUID", all.x = TRUE)
fwrite(basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl_npos_bycell, file = "Report/release/SNVs/QC/highdepth_highaf_noctrl_npos_bycell.csv")
basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl_npos_bycell <- fread(file = "Report/release/SNVs/QC/highdepth_highaf_noctrl_npos_bycell.csv")

basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl_npos_bymouse <- basediffperc_cutdemux_q30_unstranded_highdepth_highaf[IsCtrl == "N", .(npos = uniqueN(pos)), by = "MouseID"]
basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl_npos_bymouse <- merge.data.table(MouseInfo, basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl_npos_bymouse, by = "MouseID", all.x = TRUE)
fwrite(basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl_npos_bymouse, file = "Report/release/SNVs/QC/highdepth_highaf_noctrl_npos_bymouse.csv")
basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl_npos_bymouse <- fread(file = "Report/release/SNVs/QC/highdepth_highaf_noctrl_npos_bymouse.csv")

fig <- ggplot(basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl_npos_bymito[IsCtrl == "N"], aes(x = CellUID, y = npos)) + geom_boxplot(size = 0.6, outlier.color = NA) + geom_jitter(size = 0.5, aes(color = MitoID)) + theme_classic(base_size = 11) + theme(plot.margin = unit(c(0.0, 0.00, 0.00, 0.00), "in"), axis.text.x = element_text(size = 8, angle = 90, vjust = 0.5, hjust = 1)) + guides(color = guide_legend(override.aes = list(size = 1))) + xlab("") + ylab("# SNV sites")
ggsave(fig, filename = "Report/release/SNVs/QC/highdepth_highaf_noctrl_npos_bycell_boxplot.pdf", width = 12, height = 6)

###########################################################################
## binary table: whether nonmissing data in mito/cell/mouse X at SNV Y
###########################################################################
support_byposmut <- fread("Report/release/SNVs/filter/basediffperc_cutdemux_sub500k_q30_unstranded_highdepth_highaf_qcfltd_support_byposmut.csv")
basediffperc_cutdemux_q30_unstranded_highdepth_noctrl_hasdata_bymito_bypos <- fread(file = "Report/release/SNVs/QC/highdepth_noctrl_hasdata_bymito_bypos.csv.gz", header = TRUE)
basediffperc_cutdemux_q30_unstranded_highdepth_noctrl_hasdata_bycell_bypos <- fread(file = "Report/release/SNVs/QC/highdepth_noctrl_hasdata_bycell_bypos.csv.gz", header = TRUE)
basediffperc_cutdemux_q30_unstranded_highdepth_noctrl_hasdata_bymouse_bypos <- fread(file = "Report/release/SNVs/QC/highdepth_noctrl_hasdata_bymouse_bypos.csv.gz", header = TRUE)

basediffperc_cutdemux_q30_unstranded_highdepth_noctrl_hasdata_bymito_byposmut <- data.table(
    basediffperc_cutdemux_q30_unstranded_highdepth_noctrl_hasdata_bymito_bypos[, 1:19], 
    basediffperc_cutdemux_q30_unstranded_highdepth_noctrl_hasdata_bymito_bypos[, match(support_byposmut[, pos], names(basediffperc_cutdemux_q30_unstranded_highdepth_noctrl_hasdata_bymito_bypos)), with = FALSE]
)
names(basediffperc_cutdemux_q30_unstranded_highdepth_noctrl_hasdata_bymito_byposmut)[-c(1:19)] <- support_byposmut[, posmut]
fwrite(basediffperc_cutdemux_q30_unstranded_highdepth_noctrl_hasdata_bymito_byposmut, file = "Report/release/SNVs/QC/highdepth_noctrl_hasdata_bymito_byposmut.csv.gz")
basediffperc_cutdemux_q30_unstranded_highdepth_noctrl_hasdata_bymito_byposmut <- fread(file = "Report/release/SNVs/QC/highdepth_noctrl_hasdata_bymito_byposmut.csv.gz", header = TRUE)

basediffperc_cutdemux_q30_unstranded_highdepth_noctrl_hasdata_bycell_byposmut <- data.table(
    basediffperc_cutdemux_q30_unstranded_highdepth_noctrl_hasdata_bycell_bypos[, 1:7], 
    basediffperc_cutdemux_q30_unstranded_highdepth_noctrl_hasdata_bycell_bypos[, match(support_byposmut[, pos], names(basediffperc_cutdemux_q30_unstranded_highdepth_noctrl_hasdata_bycell_bypos)), with = FALSE]
)
names(basediffperc_cutdemux_q30_unstranded_highdepth_noctrl_hasdata_bycell_byposmut)[-c(1:7)] <- support_byposmut[, posmut]
fwrite(basediffperc_cutdemux_q30_unstranded_highdepth_noctrl_hasdata_bycell_byposmut, file = "Report/release/SNVs/QC/highdepth_noctrl_hasdata_bycell_byposmut.csv.gz")
basediffperc_cutdemux_q30_unstranded_highdepth_noctrl_hasdata_bycell_byposmut <- fread(file = "Report/release/SNVs/QC/highdepth_noctrl_hasdata_bycell_byposmut.csv.gz", header = TRUE)

basediffperc_cutdemux_q30_unstranded_highdepth_noctrl_hasdata_bymouse_byposmut <- data.table(
    basediffperc_cutdemux_q30_unstranded_highdepth_noctrl_hasdata_bymouse_bypos[, 1:2], 
    basediffperc_cutdemux_q30_unstranded_highdepth_noctrl_hasdata_bymouse_bypos[, match(support_byposmut[, pos], names(basediffperc_cutdemux_q30_unstranded_highdepth_noctrl_hasdata_bymouse_bypos)), with = FALSE]
)
names(basediffperc_cutdemux_q30_unstranded_highdepth_noctrl_hasdata_bymouse_byposmut)[-c(1:2)] <- support_byposmut[, posmut]
fwrite(basediffperc_cutdemux_q30_unstranded_highdepth_noctrl_hasdata_bymouse_byposmut, file = "Report/release/SNVs/QC/highdepth_noctrl_hasdata_bymouse_byposmut.csv.gz")
basediffperc_cutdemux_q30_unstranded_highdepth_noctrl_hasdata_bymouse_byposmut <- fread(file = "Report/release/SNVs/QC/highdepth_noctrl_hasdata_bymouse_byposmut.csv.gz", header = TRUE)

###########################################################################
## Per SNV, how many mitos/cells/mice have non-missing data?
###########################################################################
basediffperc_cutdemux_q30_unstranded_highdepth_noctrl_hasdata_bymito_byposmut <- fread(file = "Report/release/SNVs/QC/highdepth_noctrl_hasdata_bymito_byposmut.csv.gz", header = TRUE)
basediffperc_cutdemux_q30_unstranded_highdepth_noctrl_hasdata_bycell_byposmut <- fread(file = "Report/release/SNVs/QC/highdepth_noctrl_hasdata_bycell_byposmut.csv.gz", header = TRUE)
basediffperc_cutdemux_q30_unstranded_highdepth_noctrl_hasdata_bymouse_byposmut <- fread(file = "Report/release/SNVs/QC/highdepth_noctrl_hasdata_bymouse_byposmut.csv.gz", header = TRUE)

basediffperc_cutdemux_q30_unstranded_highdepth_noctrl_nsampleshasdata_byposmut <- data.table(
    posmut = posmut, 
    nmitoshasdata = colSums(basediffperc_cutdemux_q30_unstranded_highdepth_noctrl_hasdata_bymito_byposmut[, support_byposmut[, posmut], with = FALSE]), 
    ncellshasdata = colSums(basediffperc_cutdemux_q30_unstranded_highdepth_noctrl_hasdata_bycell_byposmut[, support_byposmut[, posmut], with = FALSE]), 
    nmicehasdata = colSums(basediffperc_cutdemux_q30_unstranded_highdepth_noctrl_hasdata_bymouse_byposmut[, support_byposmut[, posmut], with = FALSE])
)
fwrite(basediffperc_cutdemux_q30_unstranded_highdepth_noctrl_nsampleshasdata_byposmut, file = "Report/release/SNVs/QC/highdepth_noctrl_nsampleshasdata_byposmut.csv")
basediffperc_cutdemux_q30_unstranded_highdepth_noctrl_nsampleshasdata_byposmut <- fread(file = "Report/release/SNVs/QC/highdepth_noctrl_nsampleshasdata_byposmut.csv", header = TRUE)

###########################################################################
## binary table: whether mito/cell/mouse X carries SNV Y?
## For each mito/cell/mouse, how mamy total SNVs does it have?
###########################################################################
support_byposmut <- fread(file = "Report/release/SNVs/filter/basediffperc_cutdemux_sub500k_q30_unstranded_highdepth_highaf_qcfltd_support_byposmut.csv")
basediffperc_cutdemux_q30_unstranded_highdepth_highaf <- fread("Report/release/SNVs/filter/basediffperc_cutdemux_sub500k_q30_unstranded_highdepth_highaf_qcfltd.csv.gz")
basediffperc_cutdemux_q30_unstranded_highdepth_highaf[ref == "A", A := `=`]
basediffperc_cutdemux_q30_unstranded_highdepth_highaf[ref == "C", C := `=`]
basediffperc_cutdemux_q30_unstranded_highdepth_highaf[ref == "G", G := `=`]
basediffperc_cutdemux_q30_unstranded_highdepth_highaf[ref == "T", T := `=`]
basediffperc_cutdemux_q30_unstranded_highdepth_highaf <- melt.data.table(basediffperc_cutdemux_q30_unstranded_highdepth_highaf, measure.vars = c("A", "C", "G", "T", "del"), variable.name = "alt", value.name = "altperc")
setnames(basediffperc_cutdemux_q30_unstranded_highdepth_highaf_qcfltd_noctrl, "=", "refperc")
basediffperc_cutdemux_q30_unstranded_highdepth_highaf[, alt := factor(alt, levels = c("A", "C", "G", "T", "del"))]
basediffperc_cutdemux_q30_unstranded_highdepth_highaf[, mut := factor(paste0(ref, ">", alt), levels = paste0(rep(c("A", "C", "G", "T"), each = 5), ">", c("A", "C", "G", "T", "del")))]
setkey(basediffperc_cutdemux_q30_unstranded_highdepth_highaf, pos, alt)
basediffperc_cutdemux_q30_unstranded_highdepth_highaf[, posmut := paste0(pos, ":", mut)]
basediffperc_cutdemux_q30_unstranded_highdepth_highaf[, posmut := factor(posmut, levels = unique(posmut))]
dim(basediffperc_cutdemux_q30_unstranded_highdepth_highaf)
## [1] 35785    28

vaf_th <- 0.05
basediffperc_cutdemux_q30_unstranded_highdepth_highaf <- basediffperc_cutdemux_q30_unstranded_highdepth_highaf[altperc >= vaf_th * 100]
basediffperc_cutdemux_q30_unstranded_highdepth_highaf <- basediffperc_cutdemux_q30_unstranded_highdepth_highaf[ref != alt]
dim(basediffperc_cutdemux_q30_unstranded_highdepth_highaf)
## [1] 7270     28
(basediffperc_cutdemux_q30_unstranded_highdepth_highaf[, setdiff(posmut, support_byposmut[, posmut])])
## [1] "1267:C>T"  "1273:C>G"  "3846:G>A"  "15197:C>G"
(basediffperc_cutdemux_q30_unstranded_highdepth_highaf[IsCtrl == "N", setdiff(posmut, support_byposmut[, posmut])])
## character(0)
## There extra SNVs are unique to the control samples, hence they don't appear in the 1032 total list.

basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl_hassnv_bymito_byposmut <- dcast.data.table(basediffperc_cutdemux_q30_unstranded_highdepth_highaf[IsCtrl == "N"],  LibraryMitoID ~ posmut, value.var = "alt", fun.aggregate = uniqueN)
basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl_hassnv_bymito_byposmut <- MitoInfo[basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl_hassnv_bymito_byposmut, on = "LibraryMitoID"]
fwrite(basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl_hassnv_bymito_byposmut, file = "Report/release/SNVs/QC/highdepth_highaf_noctrl_hassnv_bymito_byposmut.csv.gz")
basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl_hassnv_bymito_byposmut <- fread(file = "Report/release/SNVs/QC/highdepth_highaf_noctrl_hassnv_bymito_byposmut.csv.gz")

basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl_hassnv_bycell_byposmut <- dcast.data.table(basediffperc_cutdemux_q30_unstranded_highdepth_highaf[IsCtrl == "N"], CellUID ~ posmut, value.var = "alt", fun.aggregate = uniqueN)
basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl_hassnv_bycell_byposmut <- CellInfo[basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl_hassnv_bycell_byposmut, on = "CellUID"]
fwrite(basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl_hassnv_bycell_byposmut, file = "Report/release/SNVs/QC/highdepth_highaf_noctrl_hassnv_bycell_byposmut.csv.gz")
basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl_hassnv_bycell_byposmut <- fread(file = "Report/release/SNVs/QC/highdepth_highaf_noctrl_hassnv_bycell_byposmut.csv.gz")

basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl_hassnv_bymouse_byposmut <- dcast.data.table(basediffperc_cutdemux_q30_unstranded_highdepth_highaf[IsCtrl == "N"], MouseID ~ posmut, value.var = "alt", fun.aggregate = uniqueN)
basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl_hassnv_bymouse_byposmut <- MouseInfo[basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl_hassnv_bymouse_byposmut, on = "MouseID"]
fwrite(basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl_hassnv_bymouse_byposmut, file = "Report/release/SNVs/QC/highdepth_highaf_noctrl_hassnv_bymouse_byposmut.csv.gz")
basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl_hassnv_bymouse_byposmut <- fread(file = "Report/release/SNVs/QC/highdepth_highaf_noctrl_hassnv_bymouse_byposmut.csv.gz")

###########################################################################
## For each SNV, how many mitos carry it per cell, how many cells carry it per mouse?
###########################################################################
basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl_nmitoshassnv_bycell_byposmut <- dcast.data.table(basediffperc_cutdemux_q30_unstranded_highdepth_highaf[IsCtrl == "N"], CellUID ~ posmut, value.var = "LibraryMitoID", fun.aggregate = uniqueN)
basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl_nmitoshassnv_bycell_byposmut <- CellInfo[basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl_nmitoshassnv_bycell_byposmut, on = "CellUID"]
fwrite(basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl_nmitoshassnv_bycell_byposmut, file = "Report/release/SNVs/QC/highdepth_highaf_noctrl_nmitoshassnv_bycell_byposmut.csv.gz")
basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl_nmitoshassnv_bycell_byposmut <- fread(file = "Report/release/SNVs/QC/highdepth_highaf_noctrl_nmitoshassnv_bycell_byposmut.csv.gz")

basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl_ncellshassnv_bymouse_byposmut <- dcast.data.table(basediffperc_cutdemux_q30_unstranded_highdepth_highaf[IsCtrl == "N"], MouseID ~ posmut, value.var = "CellUID", fun.aggregate = uniqueN)
basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl_ncellshassnv_bymouse_byposmut <- MouseInfo[basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl_ncellshassnv_bymouse_byposmut, on = "MouseID"]
fwrite(basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl_ncellshassnv_bymouse_byposmut, file = "Report/release/SNVs/QC/highdepth_highaf_noctrl_ncellshassnv_bymouse_byposmut.csv.gz")
basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl_ncellshassnv_bymouse_byposmut <- fread(file = "Report/release/SNVs/QC/highdepth_highaf_noctrl_ncellshassnv_bymouse_byposmut.csv.gz")


basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl_nmitoshassnv_bycell_byposmut_df <- as.data.frame(basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl_nmitoshassnv_bycell_byposmut)
rownames(basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl_nmitoshassnv_bycell_byposmut_df) <- basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl_nmitoshassnv_bycell_byposmut_df[, 1]
basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl_nmitoshassnv_bycell_byposmut_df <- basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl_nmitoshassnv_bycell_byposmut_df[, -c(1:7)]

pdf("Report/release/SNVs/QC/highdepth_highaf_noctrl_nmitoshassnv_bycell_byposmut_pheatmap.pdf", width = 15, height = 7)
pheatmap(basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl_nmitoshassnv_bycell_byposmut_df, cluster_col = FALSE, cluster_row = TRUE, angle = 90, annotation_row = CellInfo_df[rownames(basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl_nmitoshassnv_bycell_byposmut_df), c("MouseID", "CellType", "CellID")], annotation_color = list(MouseID = structure(colorRampPalette(scales::brewer_pal(palette = "Set1")(9))(length(MouseInfo[, MouseID])), names = MouseInfo[, MouseID]), CellType = c("Astrocyte" = Graphics$palette20[5], "Neuron" = Graphics$palette20[4])), fontsize_col = 1, fontsize_row = 4, show_rownames = TRUE, border_color = NA)
dev.off()

basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl_ncellshassnv_bymouse_byposmut_df <- as.data.frame(basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl_ncellshassnv_bymouse_byposmut)
rownames(basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl_ncellshassnv_bymouse_byposmut_df) <- basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl_ncellshassnv_bymouse_byposmut_df[, 1]
basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl_ncellshassnv_bymouse_byposmut_df <- basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl_ncellshassnv_bymouse_byposmut_df[, -c(1:2)]

pdf("Report/release/SNVs/QC/highdepth_highaf_noctrl_ncellshassnv_bymouse_byposmut_pheatmap.pdf", width = 15, height = 2)
pheatmap(basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl_ncellshassnv_bymouse_byposmut_df, cluster_col = FALSE, cluster_row = TRUE, angle = 90, annotation_row = MouseInfo_df[rownames(basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl_ncellshassnv_bymouse_byposmut_df), c("RCAID"), drop = FALSE], fontsize_col = 1, fontsize_row = 8, show_rownames = TRUE, border_color = NA)
dev.off()

###########################################################################
## For each SNV, how many mitos/cells/mice carry the SNV?
###########################################################################
basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl_hassnv_bymito_byposmut <- fread(file = "Report/release/SNVs/QC/highdepth_highaf_noctrl_hassnv_bymito_byposmut.csv.gz")
basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl_hassnv_bycell_byposmut <- fread(file = "Report/release/SNVs/QC/highdepth_highaf_noctrl_hassnv_bycell_byposmut.csv.gz")
basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl_hassnv_bymouse_byposmut <- fread(file = "Report/release/SNVs/QC/highdepth_highaf_noctrl_hassnv_bymouse_byposmut.csv.gz")
basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl_nsampleshassnv_byposmut <- data.table(
    posmut = support_byposmut[, posmut], 
    nmitoshassnv = colSums(basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl_hassnv_bymito_byposmut[, support_byposmut[, posmut], with = FALSE]), 
    ncellshassnv = colSums(basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl_hassnv_bycell_byposmut[, support_byposmut[, posmut], with = FALSE]), 
    nmicehassnv = colSums(basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl_hassnv_bymouse_byposmut[, support_byposmut[, posmut], with = FALSE])
)
fwrite(basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl_nsampleshassnv_byposmut, file = "Report/release/SNVs/QC/highdepth_highaf_noctrl_nsampleshassnv_byposmut.csv")
basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl_nsampleshassnv_byposmut <- fread(file = "Report/release/SNVs/QC/highdepth_highaf_noctrl_nsampleshassnv_byposmut.csv")

###########################################################################
## Get general nsamples stats for each SNV 
###########################################################################
basediffperc_cutdemux_q30_unstranded_highdepth_noctrl_nsampleshasdata_byposmut <- fread(file = "Report/release/SNVs/QC/highdepth_noctrl_nsampleshasdata_byposmut.csv")
setkey(basediffperc_cutdemux_q30_unstranded_highdepth_noctrl_nsampleshasdata_byposmut, posmut)
basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl_nsampleshassnv_byposmut <- fread(file = "Report/release/SNVs/QC/highdepth_highaf_noctrl_nsampleshassnv_byposmut.csv")
setkey(basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl_nsampleshassnv_byposmut, posmut)

basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl_nsamples_byposmut <- data.table(
    support_byposmut, 
    basediffperc_cutdemux_q30_unstranded_highdepth_noctrl_nsampleshasdata_byposmut[support_byposmut[, posmut], -1, on = "posmut"], 
    basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl_nsampleshassnv_byposmut[support_byposmut[, posmut], -1, on = "posmut"]
)
basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl_nsamples_byposmut[, c("fmitoshasdata", "fcellshasdata", "fmicehasdata") := list(nmitoshasdata / MitoInfo[IsCtrl == "N", uniqueN(LibraryMitoID)], ncellshasdata / MitoInfo[IsCtrl == "N", uniqueN(CellUID)], nmicehasdata / MitoInfo[IsCtrl == "N", uniqueN(MouseID)])]
basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl_nsamples_byposmut[, c("fmitoshassnv", "fcellshassnv", "fmicehassnv") := list(nmitoshassnv / nmitoshasdata, ncellshassnv / ncellshasdata, nmicehassnv / nmicehasdata)]
basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl_nsamples_byposmut[, c("fmitoshasdatahassnv", "fcellshasdatahassnv", "fmicehasdatahassnv") := list(nmitoshassnv / MitoInfo[IsCtrl == "N", uniqueN(LibraryMitoID)], ncellshassnv / MitoInfo[IsCtrl == "N", uniqueN(CellUID)], nmicehassnv / MitoInfo[IsCtrl == "N", uniqueN(MouseID)])]

fwrite(basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl_nsamples_byposmut, file = "Report/release/SNVs/QC/highdepth_highaf_noctrl_nsamples_byposmut.csv")
basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl_nsamples_byposmut <- fread(file = "Report/release/SNVs/QC/highdepth_highaf_noctrl_nsamples_byposmut.csv")

###########################################################################
## # of variant alleles per mito/cell/mouse per SNV site
###########################################################################
support_byposmut <- fread(file = "Report/release/SNVs/filter/basediffperc_cutdemux_sub500k_q30_unstranded_highdepth_highaf_qcfltd_support_byposmut.csv")
basediffperc_cutdemux_q30_unstranded_highdepth_highaf <- fread("Report/release/SNVs/filter/basediffperc_cutdemux_sub500k_q30_unstranded_highdepth_highaf_qcfltd.csv.gz")
basediffperc_cutdemux_q30_unstranded_highdepth_highaf[ref == "A", A := `=`]
basediffperc_cutdemux_q30_unstranded_highdepth_highaf[ref == "C", C := `=`]
basediffperc_cutdemux_q30_unstranded_highdepth_highaf[ref == "G", G := `=`]
basediffperc_cutdemux_q30_unstranded_highdepth_highaf[ref == "T", T := `=`]
basediffperc_cutdemux_q30_unstranded_highdepth_highaf <- melt.data.table(basediffperc_cutdemux_q30_unstranded_highdepth_highaf, measure.vars = c("A", "C", "G", "T", "del"), variable.name = "alt", value.name = "altperc")
setnames(basediffperc_cutdemux_q30_unstranded_highdepth_highaf_qcfltd_noctrl, "=", "refperc")
basediffperc_cutdemux_q30_unstranded_highdepth_highaf[, alt := factor(alt, levels = c("A", "C", "G", "T", "del"))]
basediffperc_cutdemux_q30_unstranded_highdepth_highaf[, mut := factor(paste0(ref, ">", alt), levels = paste0(rep(c("A", "C", "G", "T"), each = 5), ">", c("A", "C", "G", "T", "del")))]
setkey(basediffperc_cutdemux_q30_unstranded_highdepth_highaf, pos, alt)
basediffperc_cutdemux_q30_unstranded_highdepth_highaf[, posmut := paste0(pos, ":", mut)]
basediffperc_cutdemux_q30_unstranded_highdepth_highaf[, posmut := factor(posmut, levels = unique(posmut))]
dim(basediffperc_cutdemux_q30_unstranded_highdepth_highaf)

vaf_th <- 0.05
basediffperc_cutdemux_q30_unstranded_highdepth_highaf <- basediffperc_cutdemux_q30_unstranded_highdepth_highaf[altperc >= vaf_th * 100]
basediffperc_cutdemux_q30_unstranded_highdepth_highaf <- basediffperc_cutdemux_q30_unstranded_highdepth_highaf[ref != alt]

basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl_nvaralleles_bymito_bypos <- dcast.data.table(basediffperc_cutdemux_q30_unstranded_highdepth_highaf[IsCtrl == "N"], LibraryMitoID ~ pos, value.var = "alt", fun.aggregate = uniqueN)
basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl_nvaralleles_bymito_bypos <- MitoInfo[basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl_nvaralleles_bymito_bypos, on = "LibraryMitoID"]
dim(basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl_nvaralleles_bymito_bypos)
## [1] 1519  857
fwrite(basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl_nvaralleles_bymito_bypos, file = "Report/release/SNVs/QC/highdepth_highaf_noctrl_nvaralleles_bymito_bypos.csv")
basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl_nvaralleles_bymito_bypos <- fread(file = "Report/release/SNVs/QC/highdepth_highaf_noctrl_nvaralleles_bymito_bypos.csv", header = TRUE)

basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl_nvaralleles_bycell_bypos <- dcast.data.table(basediffperc_cutdemux_q30_unstranded_highdepth_highaf[IsCtrl == "N"], CellUID ~ pos, value.var = "alt", fun.aggregate = uniqueN)
basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl_nvaralleles_bycell_bypos <- CellInfo[basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl_nvaralleles_bycell_bypos, on = "CellUID"]
dim(basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl_nvaralleles_bycell_bypos)
## [1] 102 845 
fwrite(basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl_nvaralleles_bycell_bypos, file = "Report/release/SNVs/QC/highdepth_highaf_noctrl_nvaralleles_bycell_bypos.csv")
basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl_nvaralleles_bycell_bypos <- fread(file = "Report/release/SNVs/QC/highdepth_highaf_noctrl_nvaralleles_bycell_bypos.csv", header = TRUE)

basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl_nvaralleles_bymouse_bypos <- dcast.data.table(basediffperc_cutdemux_q30_unstranded_highdepth_highaf[IsCtrl == "N"], MouseID ~ pos, value.var = "alt", fun.aggregate = uniqueN)
basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl_nvaralleles_bymouse_bypos <- MouseInfo[basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl_nvaralleles_bymouse_bypos, on = "MouseID"]
dim(basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl_nvaralleles_bymouse_bypos)
## [1] 13  840
fwrite(basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl_nvaralleles_bymouse_bypos, file = "Report/release/SNVs/QC/highdepth_highaf_noctrl_nvaralleles_bymouse_bypos.csv")
basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl_nvaralleles_bymouse_bypos <- fread(file = "Report/release/SNVs/QC/highdepth_highaf_noctrl_nvaralleles_bymouse_bypos.csv", header = TRUE)

basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl_nvaralleles_bymito_bypos_df <- as.data.frame(basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl_nvaralleles_bymito_bypos[, -c(1:19)])
rownames(basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl_nvaralleles_bymito_bypos_df) <- basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl_nvaralleles_bymito_bypos[, LibraryMitoID]
MitoInfo_df <- as.data.frame(MitoInfo)
rownames(MitoInfo_df) <- MitoInfo_df$LibraryMitoID
MitoInfo_df$MitoID <- factor(MitoInfo_df$MitoID, levels = mitoIDs)
MitoInfo_df$CellID <- as.integer(MitoInfo_df$CellID)
MitoInfo_df$MouseID <- factor(MitoInfo_df$MouseID, levels = sort(MouseInfo[, MouseID]))
basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl_nvaralleles_bymito_bypos_df <- basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl_nvaralleles_bymito_bypos_df[with(MitoInfo_df[rownames(basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl_nvaralleles_bymito_bypos_df), ], order(MouseID, CellType, CellID, MitoID)), ]

pdf("Report/release/SNVs/QC/highdepth_highaf_noctrl_nvaralleles_bymito_bypos_heatmap.pdf", width = 12, height = 7)
pheatmap(basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl_nvaralleles_bymito_bypos_df, cluster_col = FALSE, cluster_row = FALSE, angle = 90, annotation_row = MitoInfo_df[rownames(basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl_nvaralleles_bymito_bypos_df), c("MouseID", "CellType", "CellID", "MitoID")], annotation_color = list(MouseID = structure(colorRampPalette(scales::brewer_pal(palette = "Set1")(9))(length(MouseInfo[, MouseID])), names = MouseInfo[, MouseID]), CellType = c("Astrocyte" = Graphics$palette20[5], "Neuron" = Graphics$palette20[4])), fontsize_col = 1, show_rownames = FALSE, na_col = "#FFFFFFFF")
dev.off()

basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl_nvaralleles_bycell_bypos_df <- as.data.frame(basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl_nvaralleles_bycell_bypos[, -c(1:7)])
rownames(basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl_nvaralleles_bycell_bypos_df) <- basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl_nvaralleles_bycell_bypos[, CellUID]
CellInfo_df <- as.data.frame(CellInfo)
rownames(CellInfo_df) <- CellInfo_df$CellUID
CellInfo_df$CellID <- as.integer(CellInfo_df$CellID)
CellInfo_df$MouseID <- factor(CellInfo_df$MouseID, levels = sort(MouseInfo[, MouseID]))
basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl_nvaralleles_bycell_bypos_df <- basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl_nvaralleles_bycell_bypos_df[with(CellInfo_df[rownames(basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl_nvaralleles_bycell_bypos_df), ], order(MouseID, CellType, CellID)), ]

pdf("Report/release/SNVs/QC/highdepth_highaf_noctrl_nvaralleles_bycell_bypos_heatmap.pdf", width = 12, height = 4.5)
pheatmap(basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl_nvaralleles_bycell_bypos_df, cluster_col = FALSE, cluster_row = FALSE, angle = 90, annotation_row = CellInfo_df[rownames(basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl_nvaralleles_bycell_bypos_df), c("MouseID", "CellType", "CellID")], annotation_color = list(MouseID = structure(colorRampPalette(scales::brewer_pal(palette = "Set1")(9))(length(MouseInfo[, MouseID])), names = MouseInfo[, MouseID]), CellType = c("Astrocyte" = Graphics$palette20[5], "Neuron" = Graphics$palette20[4])), fontsize_col = 1, show_rownames = FALSE, na_col = "#FFFFFFFF")
dev.off()

basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl_nvaralleles_bymouse_bypos_df <- as.data.frame(basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl_nvaralleles_bymouse_bypos[, -c(1:2)])
rownames(basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl_nvaralleles_bymouse_bypos_df) <- basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl_nvaralleles_bymouse_bypos[, MouseID]
MouseInfo_df <- as.data.frame(MouseInfo)
rownames(MouseInfo_df) <- MouseInfo_df$MouseID
MouseInfo_df$MouseID <- factor(MouseInfo_df$MouseID, levels = sort(MouseInfo[, MouseID]))
MouseInfo_df$RCAID <- factor(MouseInfo_df$RCAID, levels = sort(MouseInfo[, RCAID]))
basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl_nvaralleles_bymouse_bypos_df <- basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl_nvaralleles_bymouse_bypos_df[with(MouseInfo_df[rownames(basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl_nvaralleles_bymouse_bypos_df), ], order(MouseID, RCAID)), ]

pdf("Report/release/SNVs/QC/highdepth_highaf_noctrl_nvaralleles_bymouse_bypos_heatmap.pdf", width = 12, height = 2)
pheatmap(basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl_nvaralleles_bymouse_bypos_df, cluster_col = FALSE, cluster_row = FALSE, angle = 90, annotation_row = MouseInfo_df[rownames(basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl_nvaralleles_bymouse_bypos_df), c("RCAID"), drop = FALSE], annotation_color = list(MouseID = structure(colorRampPalette(scales::brewer_pal(palette = "Set1")(9))(length(MouseInfo[, MouseID])), names = MouseInfo[, MouseID]), CellType = c("Astrocyte" = Graphics$palette20[5], "Neuron" = Graphics$palette20[4])), fontsize_col = 1, show_rownames = TRUE, na_col = "#FFFFFFFF")
dev.off()

###########################################################################
## show that VAF is correlated to the number of mitos/cells/mice sharing
###########################################################################
support_byposmut <- fread(file = "Report/release/SNVs/filter/basediffperc_cutdemux_sub500k_q30_unstranded_highdepth_highaf_qcfltd_support_byposmut.csv")
basediffperc_cutdemux_q30_unstranded_highdepth <- fread("Report/release/SNVs/filter/basediffperc_cutdemux_sub500k_q30_unstranded_highdepth.csv.gz")
basediffperc_cutdemux_q30_unstranded_highdepth[ref == "A", A := `=`]
basediffperc_cutdemux_q30_unstranded_highdepth[ref == "C", C := `=`]
basediffperc_cutdemux_q30_unstranded_highdepth[ref == "G", G := `=`]
basediffperc_cutdemux_q30_unstranded_highdepth[ref == "T", T := `=`]
basediffperc_cutdemux_q30_unstranded_highdepth <- melt(basediffperc_cutdemux_q30_unstranded_highdepth, measure.vars = c("A", "C", "G", "T", "del"), variable.name = "alt", value.name = "altperc")
setnames(basediffperc_cutdemux_q30_unstranded_highdepth, "=", "refperc")
basediffperc_cutdemux_q30_unstranded_highdepth[, alt := factor(alt, levels = c("A", "C", "G", "T", "del"))]
basediffperc_cutdemux_q30_unstranded_highdepth[, mut := factor(paste0(ref, ">", alt), levels = paste0(rep(c("A", "C", "G", "T"), each = 5), ">", c("A", "C", "G", "T", "del")))]
setkey(basediffperc_cutdemux_q30_unstranded_highdepth, pos, alt)
basediffperc_cutdemux_q30_unstranded_highdepth[, posmut := paste0(pos, ":", mut)]
basediffperc_cutdemux_q30_unstranded_highdepth[, posmut := factor(posmut, levels = unique(posmut))]
dim(basediffperc_cutdemux_q30_unstranded_highdepth)
## [1] 8512540      28
## after removing SNVs inside PCR forward primer.
## [1] 8512540      28 
## Nothing changed and this is because highdepth is before QC filter.
## [1] 8515765      28 after fixing 9027
basediffperc_cutdemux_q30_unstranded_highdepth <- basediffperc_cutdemux_q30_unstranded_highdepth[as.character(posmut) %in% support_byposmut[, unique(posmut)]]
dim(basediffperc_cutdemux_q30_unstranded_highdepth)
## [1] 1045370      28
## [1] 1018393      28
## [1] 1021830      28 after fixing 9027
basediffperc_cutdemux_q30_unstranded_highdepth_noctrl <- basediffperc_cutdemux_q30_unstranded_highdepth[IsCtrl == "N"] 
dim(basediffperc_cutdemux_q30_unstranded_highdepth_noctrl)
## [1] 995513     28
## [1] 969813     28
## [1] 973133     28 after fixing 9027

basediffperc_cutdemux_q30_unstranded_highdepth_avgvaf_support_byposmut <- merge.data.table(basediffperc_cutdemux_q30_unstranded_highdepth_noctrl[, .(avgvaf = mean(altperc)), by = "posmut"], support_byposmut, by = "posmut")
basediffperc_cutdemux_q30_unstranded_highdepth_avgvaf_support_byposmut[, ncells1 := cut(ncells, breaks = c(0, 10, 20, 30, 40, 50, 60, 90, 100), include.lowest = FALSE)]
basediffperc_cutdemux_q30_unstranded_highdepth_avgvaf_support_byposmut[, nmitos1 := cut(nmitos, breaks = c(0, 10, 20, 30, 40, 50, 60, 70, 80, 100, 120, 140, 160, 200, 400, 600, 800, 1000, 1400), include.lowest = FALSE, dig.lab = 5)]

pdf("Report/release/SNVs/QC/highdepth_noctrl_avgvaf_support_byposmut_boxplot.pdf", width = 6, height = 6)
par(ps = 16, lend = 2, ljoin = 1, bty = "L", mar = c(5.5, 2.5, 1, 0.5), oma = c(0, 0, 0, 0), mgp = c(1.5, 0.5, 0))
boxplot(avgvaf ~ nmitos1, xlab = "", main = "# mitos", ylab = "mean VAF (%)", horizontal = FALSE, data = basediffperc_cutdemux_q30_unstranded_highdepth_avgvaf_support_byposmut, pch = 16, col = "gray", cex = 0.5, las = 2, drop = TRUE)
par(ps = 16, lend = 2, ljoin = 1, bty = "L", mar = c(4.5, 2.5, 1, 0.5), oma = c(0, 0, 0, 0), mgp = c(1.5, 0.5, 0))
boxplot(avgvaf ~ ncells1, xlab = "", main = "# cells", ylab = "mean VAF (%)", horizontal = FALSE, data = basediffperc_cutdemux_q30_unstranded_highdepth_avgvaf_support_byposmut, pch = 16, col = "gray", cex = 0.5, las = 2, drop = TRUE)
par(ps = 16, lend = 2, ljoin = 1, bty = "L", mar = c(2.5, 2.5, 0.5, 0.5), oma = c(0, 0, 0, 0), mgp = c(1.5, 0.5, 0))
boxplot(avgvaf ~ nmice, xlab = "", main = "# mice", ylab = "mean VAF (%)", horizontal = FALSE, data = basediffperc_cutdemux_q30_unstranded_highdepth_avgvaf_support_byposmut, pch = 16, col = "gray")
dev.off()

pdf("Report/release/SNVs/QC/highdepth_highaf_noctrl_nposmut_bysupport.pdf", width = 6, height = 6)
par(ps = 16, lend = 2, ljoin = 1, bty = "L", mfrow = c(1, 1), mar = c(2.5, 2.5, 0.5, 0.5), oma = c(0, 0, 0, 0), mgp = c(1.5, 0.5, 0))
support_byposmut[, .N, by = nmice][order(nmice)][, plot(N ~ nmice,  xlab = "# mice", ylab = "# SNVs", type = 'b')]
support_byposmut[, .N, by = ncells][order(ncells)][, plot(N ~ ncells, log = "xy", xlab = "# cells", ylab = "# SNVs", type = "b")]
support_byposmut[, .N, by = nmitos][order(nmitos)][, plot(N ~ nmitos, log = "xy", xlab = "# mitos", ylab = "# SNVs", type = "b")]
dev.off()

###########################################################################
## How is the # of mice/cells/mitos at each SNV site correlated to the 
## # of variant alleles at each SNV site? 
## Here we use all 838 sites.
###########################################################################
support_bypos <- fread(file = "Report/release/SNVs/filter/basediffperc_cutdemux_sub500k_q30_unstranded_highdepth_highaf_qcfltd_support_bypos.csv")
nvaralleles_bypos <- support_byposmut[, list(nvaralleles = uniqueN(posmut)), by = 'pos']
support_nvaralleles_bypos <- support_bypos[nvaralleles_bypos, on = "pos"]

pdf("Report/release/SNVs/QC/highdepth_highaf_noctrl_support_bynvaralleles.pdf", width = 6, height = 6)
ggplot(support_nvaralleles_bypos, aes(x = nvaralleles, y = nmitos, group = nvaralleles)) + theme_classic(base_size = 16) + geom_boxplot(outlier.color = NA) + geom_jitter(size = 0.5) + xlab("# variant alleles") + ylab("# mitos sharing")
ggplot(support_nvaralleles_bypos, aes(x = nvaralleles, y = ncells, group = nvaralleles)) + theme_classic(base_size = 16) + geom_boxplot(outlier.color = NA) + geom_jitter(size = 0.5) + xlab("# variant alleles") + ylab("# cells sharing")
ggplot(support_nvaralleles_bypos, aes(x = nvaralleles, y = nmice, group = nvaralleles)) + theme_classic(base_size = 16) + geom_boxplot(outlier.color = NA) + geom_jitter(size = 0.5) + xlab("# variant alleles") + ylab("# mice sharing")
dev.off()
