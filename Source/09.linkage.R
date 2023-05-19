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
library("circlize")

snv_info <- fread("Data/snv_loci_v2.csv")
snvIDs <- snv_info[, SNVID]
mito_barcodes <- read.csv("Data/mito_barcodes.csv", as.is = TRUE)
mitoIDs <- mito_barcodes[, "ID"] 

chrmbases_properties <- fread("Report/artifact/chrmbases_properties.csv.gz")
chrmbases <- chrmbases_properties[, ref]
nchrmbases <- length(chrmbases)
nchrmbases
## [1] 16299

chrmgenes <- fread("Data/mm10.mito/chrM.gtf")
chrmgenes[, symbol := sapply(regmatches(V9, regexec('gene_name "(.*)"', V9, perl = TRUE)), '[', 2)]
chrmgenes <- chrmgenes[, list(symbol, V4, V5, V7)]

bases <- c("A", "G", "C", "T", "-")

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

MitoInfo_df <- as.data.frame(MitoInfo)
rownames(MitoInfo_df) <- MitoInfo_df[, "LibraryMitoID"]

CellInfo_df <- as.data.frame(CellInfo)
rownames(CellInfo_df) <- CellInfo_df$CellUID

MouseInfo_df <- as.data.frame(MouseInfo)
rownames(MouseInfo_df) <- MouseInfo_df$MouseID

############################################################################
## Because we saw a lot of linkages heppening to the "suboptimal" reads, we went to 
## filter these reads out and did the SNV AF call again. 
## Here, we decided to work on the strictly filtered data for fears that alignment
## artifacts could impact linkage analysis the most.
############################################################################
loose_support_byposmut <- fread(file = "Report/SNVs/filter/basediffperc_cutdemux_sub500k_q30_unstranded_highdepth_highaf_qcfltd_support_byposmut.csv")
strict_support_byposmut <- fread(file = "Report/SNVs/filter/basediffperc_cutdemux_strict_sub500k_q30_unstranded_highdepth_highaf_qcfltd_support_byposmut.csv")
support_byposmut <- strict_support_byposmut[pos %in% loose_support_byposmut[, unique(`pos`)]]
inherited_support_byposmut <- support_byposmut[nmice >= 3]
inherited_pos <- inherited_support_byposmut[, unique(pos)]

############################################################################
## Linkage per cell
############################################################################
inherited_allele_bymito_bypos <- fread("Report/SNVs/filter/basediffperc_cutdemux_strict_sub500k_q30_unstranded_highdepth_inherited_allele_bymito_bypos.csv.gz")
inherited_allele_bymito_bypos <- inherited_allele_bymito_bypos[IsCtrl == "N", c(1:19, match(inherited_pos, names(inherited_allele_bymito_bypos))), with = FALSE]
dim(inherited_allele_bymito_bypos)
## [1] 1610  155

mouseIDs <- MouseInfo[, sort(MouseID)] ##  != "Mouse16&17"
cellUIDs <- sort(unique(inherited_allele_bymito_bypos[, CellUID]))
libraryMitoIDs <- unique(inherited_allele_bymito_bypos[, LibraryMitoID])

nmitosbycell <- inherited_allele_bymito_bypos[LibraryMitoID %in% libraryMitoIDs, uniqueN(LibraryMitoID), keyby = "CellUID"]

## We want each cell to have >= 20 mitos
cellUIDs <- nmitosbycell[V1 >= 20, CellUID]
cellUIDs
##  [1] "Mouse03_Neuron_2"     "Mouse03_Neuron_7"     "Mouse04_Neuron_11"
##  [4] "Mouse04_Neuron_2"     "Mouse05_Neuron_4"     "Mouse05_Neuron_6"
##  [7] "Mouse06_Neuron_2"     "Mouse06_Neuron_3"     "Mouse07_Astrocyte_1"
## [10] "Mouse07_Astrocyte_2"  "Mouse07_Astrocyte_4"  "Mouse07_Astrocyte_5"
## [13] "Mouse07_Astrocyte_6"  "Mouse07_Astrocyte_8"  "Mouse08_Astrocyte_21"
## [16] "Mouse09_Astrocyte_1"  "Mouse09_Astrocyte_5"  "Mouse09_Neuron_1"
## [19] "Mouse09_Neuron_3"     "Mouse09_Neuron_5"     "Mouse12_Astrocyte_1"
## [22] "Mouse12_Astrocyte_2"  "Mouse12_Astrocyte_3"  "Mouse12_Neuron_2"
## [25] "Mouse12_Neuron_3"     "Mouse16&17_Neuron_4"

inherited_allele_bymito_bypos_bycell <- sapply(cellUIDs, function(x) inherited_allele_bymito_bypos[CellUID == x], simplify = FALSE)
Tools$write_xlsx(inherited_allele_bymito_bypos_bycell, file = "Report/SNVs/linkage/inherited_allele_bymito_bypos_bycell.xlsx", row.names = FALSE)

## further filter out sites that have < 10 non-missing mitos
inherited_allele_nafltd_bymito_bypos_bycell <- lapply(inherited_allele_bymito_bypos_bycell, function(Y) { 
    X <- as.matrix(Y[, !1:19])
    rownames(X) <- Y[, LibraryMitoID]
    nonna <- apply(X, 2, function(x) sum(x != ""))
    Y <- X[, nonna >= 10]
})
sapply(inherited_allele_nafltd_bymito_bypos_bycell, dim)
##      Mouse03_Neuron_2 Mouse03_Neuron_7 Mouse04_Neuron_11 Mouse04_Neuron_2
## [1,]               20               24                22               23
## [2,]              107              105                73              104
##      Mouse05_Neuron_4 Mouse05_Neuron_6 Mouse06_Neuron_2 Mouse06_Neuron_3
## [1,]               28               31               25               25
## [2,]              120              120               91              104
##      Mouse07_Astrocyte_1 Mouse07_Astrocyte_2 Mouse07_Astrocyte_4
## [1,]                  53                  20                  49
## [2,]                 120                 120                 121
##      Mouse07_Astrocyte_5 Mouse07_Astrocyte_6 Mouse07_Astrocyte_8
## [1,]                  33                  34                  24
## [2,]                 120                 120                 112
##      Mouse08_Astrocyte_21 Mouse09_Astrocyte_1 Mouse09_Astrocyte_5
## [1,]                   28                  33                  49
## [2,]                  120                 119                 126
##      Mouse09_Neuron_1 Mouse09_Neuron_3 Mouse09_Neuron_5 Mouse12_Astrocyte_1
## [1,]               32               50               47                  36
## [2,]              117              123              105                 128
##      Mouse12_Astrocyte_2 Mouse12_Astrocyte_3 Mouse12_Neuron_2 Mouse12_Neuron_3
## [1,]                  34                  26               35               31
## [2,]                 122                 101              130               91
##      Mouse16&17_Neuron_4
## [1,]                  31
## [2,]                 100

## filter for sites that have at least one minor allele with >= 5% frequency out of all mitos
## This is because we do not want such case as T T T ... A T where the major allele T is 99% and the minor allele A is 1%. In this case, we do not have enough variation. 
bases <- c("A", "G", "C", "T", "-")
inherited_allele_nafltd_margin_bycell <- lapply(inherited_allele_nafltd_bymito_bypos_bycell, function(X) {
        apply(X, 2, function(x) sapply(bases, function(b) mean(grepl(x[x != "" & !is.na(x)], pattern = b))))
})

inherited_allele_nafltd_rarefltd_bycell <- mapply(function(X, Y) { y <- apply(Y, 2, function(k) { sum(k >= 0.05) >= 2 }); X[, y, drop = FALSE] }, X = inherited_allele_nafltd_bymito_bypos_bycell, Y = inherited_allele_nafltd_margin_bycell)
sapply(inherited_allele_nafltd_rarefltd_bycell, ncol)
##     Mouse03_Neuron_2     Mouse03_Neuron_7    Mouse04_Neuron_11
##                   14                   14                    6
##     Mouse04_Neuron_2     Mouse05_Neuron_4     Mouse05_Neuron_6
##                    3                    6                   10
##     Mouse06_Neuron_2     Mouse06_Neuron_3  Mouse07_Astrocyte_1
##                    6                    3                    1
##  Mouse07_Astrocyte_2  Mouse07_Astrocyte_4  Mouse07_Astrocyte_5
##                    6                    1                    3
##  Mouse07_Astrocyte_6  Mouse07_Astrocyte_8 Mouse08_Astrocyte_21
##                    5                    8                   14
##  Mouse09_Astrocyte_1  Mouse09_Astrocyte_5     Mouse09_Neuron_1
##                   10                    9                    6
##     Mouse09_Neuron_3     Mouse09_Neuron_5  Mouse12_Astrocyte_1
##                    8                    5                   10
##  Mouse12_Astrocyte_2  Mouse12_Astrocyte_3     Mouse12_Neuron_2
##                    8                    4                   13
##     Mouse12_Neuron_3  Mouse16&17_Neuron_4
##                    4                   11
Tools$write_xlsx(inherited_allele_nafltd_rarefltd_bycell[sort(cellUIDs)], file = "Report/SNVs/linkage/inherited_allele_nafltd_rarefltd_bycell.xlsx", row.names = TRUE)

bases <- c("A", "G", "C", "T", "-")
inherited_allele_nafltd_rarefltd_contab_bycell <- sapply(names(inherited_allele_nafltd_rarefltd_bycell), function(cellUID) { 
    X <- inherited_allele_nafltd_rarefltd_bycell[[cellUID]]
    Genetics$link_bin_contab(X, lab = cellUID)
}, simplify = FALSE)

inherited_allele_nafltd_rarefltd_fet_bycell <- sapply(names(inherited_allele_nafltd_rarefltd_bycell), function(cellUID) {
    X <- inherited_allele_nafltd_rarefltd_contab_bycell[[cellUID]]
    Genetics$link_bin_fet(X, lab = cellUID)
}, simplify = FALSE)

inherited_allele_nafltd_rarefltd_oddsratio_bycell <- do.call(rbind, lapply(names(inherited_allele_nafltd_rarefltd_contab_bycell), function(cellUID) { 
    X <- inherited_allele_nafltd_rarefltd_contab_bycell[[cellUID]]
    pos <- names(X)
    E <- do.call(rbind, lapply(pos, function(j) {
        do.call(rbind, lapply(pos, function(i) {
            Z <- sapply(bases, function(b) { 
                sapply(bases, function(a) {
                    message(cellUID, " ", i, " ", j, " ", a, " ", b)
                    Y <- X[[j]][[i]][[b]][[a]]
                    unname(inherited_allele_nafltd_rarefltd_fet_bycell[[cellUID]][[j]][[i]][[b]][[a]]["oddsratio"])
                }, simplify = TRUE)
            }, simplify = TRUE)
            K <- reshape2::melt(Z)
            colnames(K) <- c("allele1", "allele2", "oddsratio")
            data.table(pos1 = as.integer(i), pos2 = as.integer(j), K)
        }))
    }))
    data.table(CellUID = cellUID, E)
}))

inherited_allele_nafltd_rarefltd_pval_bycell <- do.call(rbind, lapply(names(inherited_allele_nafltd_rarefltd_contab_bycell), function(cellUID) { 
    X <- inherited_allele_nafltd_rarefltd_contab_bycell[[cellUID]]
    pos <- names(X)
    E <- do.call(rbind, lapply(pos, function(j) {
        do.call(rbind, lapply(pos, function(i) {
            Z <- sapply(bases, function(b) { 
                sapply(bases, function(a) {
                    message(cellUID, " ", i, " ", j, " ", a, " ", b)
                    Y <- X[[j]][[i]][[b]][[a]]
                    unname(inherited_allele_nafltd_rarefltd_fet_bycell[[cellUID]][[j]][[i]][[b]][[a]]["pval"])
                }, simplify = TRUE)
            }, simplify = TRUE)
            K <- reshape2::melt(Z)
            colnames(K) <- c("allele1", "allele2", "pval")
            data.table(pos1 = as.integer(i), pos2 = as.integer(j), K)
        }))
    }))
    data.table(CellUID = cellUID, E)
}))

## Note, we did a lot of unnecessary comparisons for the sake of simple computation 
## because not all allele pairs exist in our data. For example, we often have such 
## contingency tables:
##              | allele 2 == A | allele 2 != A
## allele1 == T |             0 |             0
## allele1 != T |             0 |             5
## or 
##              | allele 2 == A | allele 2 != A
## allele1 == T |             0 |             3
## allele1 != T |             0 |             5
## or
##              | allele 2 == A | allele 2 != A
## allele1 == T |             0 |             0
## allele1 != T |             4 |             5
## So we prefilter them before multiple correction. 

inherited_allele_nafltd_rarefltd_oddsratio_bycell <- inherited_allele_nafltd_rarefltd_oddsratio_bycell[pos1 < pos2]
inherited_allele_nafltd_rarefltd_pval_bycell <- inherited_allele_nafltd_rarefltd_pval_bycell[pos1 < pos2]
inherited_allele_nafltd_rarefltd_fetstat_bycell <- merge.data.table(inherited_allele_nafltd_rarefltd_oddsratio_bycell, inherited_allele_nafltd_rarefltd_pval_bycell, on = c("CellUID", "pos1", "pos2", "allele1", "allele2"))

inherited_allele_nafltd_rarefltd_fetstat_bycell <- inherited_allele_nafltd_rarefltd_fetstat_bycell[oddsratio > 0]
inherited_allele_nafltd_rarefltd_fetstat_bycell[, padj := p.adjust(pval), by = "CellUID"]

inherited_allele_nafltd_rarefltd_fetstat_bycell[, `:=`(
    ref1 = chrmbases[as.integer(pos1)], 
    ref2 = chrmbases[as.integer(pos2)])
]
inherited_allele_nafltd_rarefltd_fetstat_bycell[, `:=`(
    mut1 = factor(paste0(ref1, ">", allele1), levels = paste0(rep(c("A", "C", "G", "T"), each = 5), ">", c("A", "C", "G", "T", "-"))), 
    mut2 = factor(paste0(ref2, ">", allele2), levels = paste0(rep(c("A", "C", "G", "T"), each = 5), ">", c("A", "C", "G", "T", "-"))))
]
inherited_allele_nafltd_rarefltd_fetstat_bycell[, `:=`(
    posmut1 = paste0(pos1, ":", mut1), 
    posmut2 = paste0(pos2, ":", mut2))
]
setcolorder(inherited_allele_nafltd_rarefltd_fetstat_bycell, names(inherited_allele_nafltd_rarefltd_fetstat_bycell)[c(1:3, 9:14, 4:8)])

fwrite(inherited_allele_nafltd_rarefltd_fetstat_bycell, file = "Report/SNVs/linkage/inherited_allele_nafltd_rarefltd_fetstat_bycell.csv")
inherited_allele_nafltd_rarefltd_fetstat_bycell <- fread(file = "Report/SNVs/linkage/inherited_allele_nafltd_rarefltd_fetstat_bycell.csv")

pdf("Report/SNVs/linkage/inherited_allele_nafltd_rarefltd_fetstat_bycell_chord_padj0.1.pdf")
for (cellUID in unique(inherited_allele_nafltd_rarefltd_fetstat_bycell[, CellUID])) {
    X <- inherited_allele_nafltd_rarefltd_fetstat_bycell[padj < 0.1 & CellUID == cellUID]
    if (nrow(X) > 0) {
        Y <- unique(X[, data.table(from = pos1, to = pos2, padj = padj)])
        circos.clear()
        circos.par(start.degree = 90)
        circos.initialize(sectors = rep("chrM", 16299), x = 1:16299, ring = TRUE)
        circos.track(sectors = rep("chrM", nrow(chrmgenes)), ylim = c(0, 1), panel.fun = function(x, y) { 
            circos.rect(xleft = 1, xright = 16299, ytop = 1, ybottom = 0.8, border = NA, col = adjustcolor("grey50", alpha.f = 0.4)) }, 
        bg.border = NA, track.height = 0.15)
        for (i in 1:nrow(Y)) {
            Y[i, circos.link(sector.index1 = "chrM", point1 = from, sector.index2 = "chrM", point2 = to, h.ratio = 0.5, lty = 1, col = adjustcolor("red", alpha = 1))]
        }
        circos.xaxis(major.at = seq(0, 16299, by = 1000))
        title(main = cellUID)
        circos.clear()
    }
}
dev.off()

###########################################################################
## All mice together
###########################################################################
loose_support_byposmut <- fread(file = "Report/SNVs/filter/basediffperc_cutdemux_sub500k_q30_unstranded_highdepth_highaf_qcfltd_support_byposmut.csv")
strict_support_byposmut <- fread(file = "Report/SNVs/filter/basediffperc_cutdemux_strict_sub500k_q30_unstranded_highdepth_highaf_qcfltd_support_byposmut.csv")
support_byposmut <- strict_support_byposmut[pos %in% loose_support_byposmut[, unique(`pos`)]]
inherited_support_byposmut <- support_byposmut[nmice >= 3]
inherited_pos <- inherited_support_byposmut[, unique(pos)]

inherited_allele_bymito_bypos <- fread("Report/SNVs/filter/basediffperc_cutdemux_strict_sub500k_q30_unstranded_highdepth_inherited_allele_bymito_bypos.csv.gz")
inherited_allele_bymito_bypos <- inherited_allele_bymito_bypos[IsCtrl == "N", c(1:19, match(inherited_pos, names(inherited_allele_bymito_bypos))), with = FALSE]
dim(inherited_allele_bymito_bypos)
## [1] 1609  228

libraryMitoIDs <- inherited_allele_bymito_bypos[, LibraryMitoID]

## further filter out sites that have < 50 non-missing mitos
inherited_allele_bymito_bypos_mat <- as.matrix(inherited_allele_bymito_bypos[, !1:19])
rownames(inherited_allele_bymito_bypos_mat) <- inherited_allele_bymito_bypos[, LibraryMitoID]

inherited_allele_nafltd_bymito_bypos <- inherited_allele_bymito_bypos_mat[, apply(inherited_allele_bymito_bypos_mat, 2, function(x) sum(x != "" & !is.na(x))) >= 50]
dim(inherited_allele_nafltd_bymito_bypos)
## [1] 1609  209

## filter for sites that have at least one minor allele with >= 5% frequency out of all mitos
## This is because we do not want such case as T T T ... A T where the major allele T is 99% and the minor allele A is 1%. In this case, we do not have enough variation. 
bases <- c("A", "G", "C", "T", "-")
inherited_allele_nafltd_margin <- apply(inherited_allele_nafltd_bymito_bypos, 2, function(x) sapply(bases, function(b) mean(grepl(x[x != "" & !is.na(x)], pattern = b))))

inherited_allele_nafltd_rarefltd <- inherited_allele_nafltd_bymito_bypos[, apply(inherited_allele_nafltd_margin, 2, function(k) { sum(k >= 0.05) >= 2 }), drop = FALSE]
dim(inherited_allele_nafltd_rarefltd)
## [1] 1609   17

inherited_allele_nafltd_rarefltd_dt <- data.table(LibraryMitoID = rownames(inherited_allele_nafltd_rarefltd), inherited_allele_nafltd_rarefltd)
inherited_allele_nafltd_rarefltd_dt <- MitoInfo[inherited_allele_nafltd_rarefltd_dt, on = "LibraryMitoID"]
fwrite(inherited_allele_nafltd_rarefltd_dt, file = "Report/SNVs/linkage/inherited_allele_nafltd_rarefltd.csv")
inherited_allele_nafltd_rarefltd_dt <- fread(file = "Report/SNVs/linkage/inherited_allele_nafltd_rarefltd.csv")

bases <- c("A", "G", "C", "T", "-")
inherited_allele_nafltd_rarefltd_contab <- local({
    X <- inherited_allele_nafltd_rarefltd
    Genetics$link_bin_contab(X, lab = "all")
})

inherited_allele_nafltd_rarefltd_fet <- local({
    X <- inherited_allele_nafltd_rarefltd_contab
    Genetics$link_bin_fet(X, lab = "all")
})

inherited_allele_nafltd_rarefltd_oddsratio <- local({ 
    X <- inherited_allele_nafltd_rarefltd_fet
    pos <- names(X)
    E <- do.call(rbind, lapply(pos, function(j) {
        do.call(rbind, lapply(pos, function(i) {
            Z <- sapply(bases, function(b) { 
                sapply(bases, function(a) {
                    message(i, " ", j, " ", a, " ", b)
                    unname(X[[j]][[i]][[b]][[a]][["oddsratio"]])
                }, simplify = TRUE)
            }, simplify = TRUE)
            K <- reshape2::melt(Z)
            colnames(K) <- c("allele1", "allele2", "oddsratio")
            data.table(pos1 = as.integer(i), pos2 = as.integer(j), K)
        }))
    }))
    data.table(E)
})

inherited_allele_nafltd_rarefltd_pval <- local({
    X <- inherited_allele_nafltd_rarefltd_fet
    pos <- names(X)
    E <- do.call(rbind, lapply(pos, function(j) {
        do.call(rbind, lapply(pos, function(i) {
            Z <- sapply(bases, function(b) { 
                sapply(bases, function(a) {
                    message(i, " ", j, " ", a, " ", b)
                    unname(X[[j]][[i]][[b]][[a]][["pval"]])
                }, simplify = TRUE)
            }, simplify = TRUE)
            K <- reshape2::melt(Z)
            colnames(K) <- c("allele1", "allele2", "pval")
            data.table(pos1 = as.integer(i), pos2 = as.integer(j), K)
        }))
    }))
    data.table(E)
})

inherited_allele_nafltd_rarefltd_oddsratio <- inherited_allele_nafltd_rarefltd_oddsratio[pos1 < pos2]
inherited_allele_nafltd_rarefltd_pval <- inherited_allele_nafltd_rarefltd_pval[pos1 < pos2]
inherited_allele_nafltd_rarefltd_fetstat <- merge.data.table(inherited_allele_nafltd_rarefltd_oddsratio, inherited_allele_nafltd_rarefltd_pval, on = c("pos1", "pos2", "allele1", "allele2"))

inherited_allele_nafltd_rarefltd_fetstat <- inherited_allele_nafltd_rarefltd_fetstat[oddsratio > 0]
inherited_allele_nafltd_rarefltd_fetstat[, padj := p.adjust(pval)]
inherited_allele_nafltd_rarefltd_fetstat[, `:=`(
    ref1 = chrmbases[as.integer(pos1)], 
    ref2 = chrmbases[as.integer(pos2)])
]
inherited_allele_nafltd_rarefltd_fetstat[, `:=`(
    mut1 = factor(paste0(ref1, ">", allele1), levels = paste0(rep(c("A", "C", "G", "T"), each = 5), ">", c("A", "C", "G", "T", "-"))), 
    mut2 = factor(paste0(ref2, ">", allele2), levels = paste0(rep(c("A", "C", "G", "T"), each = 5), ">", c("A", "C", "G", "T", "-"))))
]
inherited_allele_nafltd_rarefltd_fetstat[, `:=`(
    posmut1 = paste0(pos1, ":", mut1), 
    posmut2 = paste0(pos2, ":", mut2))
]
setcolorder(inherited_allele_nafltd_rarefltd_fetstat, names(inherited_allele_nafltd_rarefltd_fetstat)[c(1:2, 8:13, 3:7)])

fwrite(inherited_allele_nafltd_rarefltd_fetstat, file = "Report/SNVs/linkage/inherited_allele_nafltd_rarefltd_fetstat.csv")
inherited_allele_nafltd_rarefltd_fetstat <- fread(file = "Report/SNVs/linkage/inherited_allele_nafltd_rarefltd_fetstat.csv")

inherited_allele_nafltd_rarefltd_fetstat_sig <- inherited_allele_nafltd_rarefltd_fetstat[padj < 0.1]
pdf("Report/SNVs/linkage/inherited_allele_nafltd_rarefltd_locuspair_chord_padj0.1.pdf", width = 6, height = 6)
chordDiagram(unique(inherited_allele_nafltd_rarefltd_fetstat_sig[, list(pos1, pos2)]))
title(main = "All mitos (padj < 0.1)")
circos.clear()
dev.off()
