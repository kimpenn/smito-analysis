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
library("reshape2")
library("pheatmap")
library("scales")
library("ggplot2")
library("ggrepel")
library("ggsignif")
library("parallel")

snv_info <- fread("Data/snv_loci_v2.csv")
snvIDs <- snv_info[,SNVID]
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
## Load the data
###########################################################################
inherited_altperc_bymito_byposmut <- fread(file = "Report/release/SNVs/filter/basediffperc_cutdemux_sub500k_q30_unstranded_highdepth_inherited_altperc_bymito_byposmut.csv.gz")

inherited_noctrl_altperc_bymito_byposmut <- inherited_altperc_bymito_byposmut[IsCtrl == "N"]
inherited_noctrl_altperc_bymito_byposmut[, CellTypeID := paste0(CellType, "_", CellID)]
dim(inherited_noctrl_altperc_bymito_byposmut)
## 1631 286
## 1631 270 after removing 16 inside-primer SNVs
## 1631 270 after fixing 9027

support_byposmut <- fread("Report/release/SNVs/filter/basediffperc_cutdemux_sub500k_q30_unstranded_highdepth_highaf_qcfltd_support_byposmut.csv")
inherited_posmuts <- support_byposmut[nmice >= 3, posmut]


mito_metacnames <- names(inherited_noctrl_altperc_bymito_byposmut)[!names(inherited_noctrl_altperc_bymito_byposmut) %in% inherited_posmuts]

## arcsin VAF per mito
inherited_noctrl_altasin_bymito_byposmut <- inherited_noctrl_altperc_bymito_byposmut[, data.table(inherited_noctrl_altperc_bymito_byposmut[, mito_metacnames, with = FALSE], apply(.SD, 2, Stats$parcsin)), .SDcols = inherited_posmuts]

## arcsin VAF per cell
inherited_noctrl_altasin_bycell_byposmut <- inherited_noctrl_altasin_bymito_byposmut[, { x <- colMeans(.SD[, inherited_posmuts, with = FALSE], na.rm = TRUE); x[is.nan(x)] <- NA; data.table(t(x)) }, by = "CellUID"]
inherited_noctrl_altasin_bycell_byposmut <- CellInfo[inherited_noctrl_altasin_bycell_byposmut, on = "CellUID"]
cell_metacnames <- names(inherited_noctrl_altasin_bycell_byposmut)[!names(inherited_noctrl_altasin_bycell_byposmut) %in% inherited_posmuts]

## arcsin VAF per mouse
inherited_noctrl_altasin_bymouse_byposmut <- inherited_noctrl_altasin_bycell_byposmut[, { x <- colMeans(.SD[, inherited_posmuts, with = FALSE], na.rm = TRUE); x[is.nan(x)] <- NA; data.table(t(x)) }, by = "MouseID"]
inherited_noctrl_altasin_bymouse_byposmut <- MouseInfo[inherited_noctrl_altasin_bymouse_byposmut, on = "MouseID"]
mouse_metacnames <- names(inherited_noctrl_altasin_bymouse_byposmut)[!names(inherited_noctrl_altasin_bymouse_byposmut) %in% inherited_posmuts]

fwrite(inherited_noctrl_altasin_bymito_byposmut, file = "Report/release/SNVs/hierarchy/inherited_noctrl_altasin_bymito_byposmut.csv")
fwrite(inherited_noctrl_altasin_bycell_byposmut, file = "Report/release/SNVs/hierarchy/inherited_noctrl_altasin_bycell_byposmut.csv")
fwrite(inherited_noctrl_altasin_bymouse_byposmut, file = "Report/release/SNVs/hierarchy/inherited_noctrl_altasin_bymouse_byposmut.csv")
inherited_noctrl_altasin_bymito_byposmut <- fread(file = "Report/release/SNVs/hierarchy/inherited_noctrl_altasin_bymito_byposmut.csv")
inherited_noctrl_altasin_bycell_byposmut <- fread(file = "Report/release/SNVs/hierarchy/inherited_noctrl_altasin_bycell_byposmut.csv")
inherited_noctrl_altasin_bymouse_byposmut <- fread(file = "Report/release/SNVs/hierarchy/inherited_noctrl_altasin_bymouse_byposmut.csv")

inherited_noctrl_altasinctrd_bymito_byposmut <- sapply(inherited_posmuts, function(posmut) {
    X <- inherited_noctrl_altasin_bymito_byposmut[, .(permito = get(posmut), LibraryMitoID, CellUID)]
    Y <- inherited_noctrl_altasin_bycell_byposmut[, .(percell = get(posmut), CellUID)]
    Z <- merge.data.table(X, Y, by = "CellUID", all.x = TRUE)
    setkey(Z, LibraryMitoID)
    Z[J(inherited_noctrl_altasin_bymito_byposmut[, LibraryMitoID]), permito - percell]
})
inherited_noctrl_altasinctrd_bymito_byposmut <- data.table(inherited_noctrl_altasin_bymito_byposmut[, mito_metacnames, with= FALSE], inherited_noctrl_altasinctrd_bymito_byposmut)

inherited_noctrl_altasinctrd_bycell_byposmut <- sapply(inherited_posmuts, function(posmut) {
    X <- inherited_noctrl_altasin_bycell_byposmut[, .(percell = get(posmut), CellUID, MouseID)]
    Y <- inherited_noctrl_altasin_bymouse_byposmut[, .(permouse = get(posmut), MouseID)]
    Z <- merge.data.table(X, Y, by = "MouseID", all.x = TRUE)
    setkey(Z, CellUID)
    Z[J(inherited_noctrl_altasin_bycell_byposmut[, CellUID]), percell - permouse]
})
inherited_noctrl_altasinctrd_bycell_byposmut <- data.table(inherited_noctrl_altasin_bycell_byposmut[, cell_metacnames, with= FALSE], inherited_noctrl_altasinctrd_bycell_byposmut)

inherited_noctrl_altasinctrd_bymouse_byposmut <- sapply(inherited_posmuts, function(posmut) {
    inherited_noctrl_altasin_bymouse_byposmut[, get(posmut)] - mean(inherited_noctrl_altasin_bymouse_byposmut[, get(posmut)], na.rm = TRUE)
})
inherited_noctrl_altasinctrd_bymouse_byposmut <- data.table(inherited_noctrl_altasin_bymouse_byposmut[, mouse_metacnames, with = FALSE], inherited_noctrl_altasinctrd_bymouse_byposmut)

fwrite(inherited_noctrl_altasinctrd_bymito_byposmut, file = "Report/release/SNVs/hierarchy/inherited_noctrl_altasinctrd_bymito_byposmut.csv")
fwrite(inherited_noctrl_altasinctrd_bycell_byposmut, file = "Report/release/SNVs/hierarchy/inherited_noctrl_altasinctrd_bycell_byposmut.csv")
fwrite(inherited_noctrl_altasinctrd_bymouse_byposmut, file = "Report/release/SNVs/hierarchy/inherited_noctrl_altasinctrd_bymouse_byposmut.csv")
inherited_noctrl_altasinctrd_bymito_byposmut <- fread(file = "Report/release/SNVs/hierarchy/inherited_noctrl_altasinctrd_bymito_byposmut.csv")
inherited_noctrl_altasinctrd_bycell_byposmut <- fread(file = "Report/release/SNVs/hierarchy/inherited_noctrl_altasinctrd_bycell_byposmut.csv")
inherited_noctrl_altasinctrd_bymouse_byposmut <- fread(file = "Report/release/SNVs/hierarchy/inherited_noctrl_altasinctrd_bymouse_byposmut.csv")

###########################################################################
## using ANOVA to decompose the variance
###########################################################################
inherited_noctrl_altasin_nestedaov_bylevel_byposmut <- sapply(inherited_posmuts, function(posmut) {
    X <- inherited_noctrl_altasin_bymito_byposmut[MouseID != "Mouse16&17", c(posmut, "MouseID", "CellUID"), with = FALSE]
    names(X)[1] <- "VAF"
    mod <- aov(VAF ~ MouseID + Error(CellUID), data = X)
    summary(mod)
}, simplify = FALSE)
## Note, here we don't have to use car::Anova for type 3 SS because there is no order issue. 

inherited_noctrl_altasin_nestedaov_df_bylevel_byposmut <- sapply(inherited_noctrl_altasin_nestedaov_bylevel_byposmut, function(X) { structure(c(X[[1]][[1]][, "Df"], X[[2]][[1]][, "Df"]), names = c("mousewise", "cellwise", "mitowise")) })
inherited_noctrl_altasin_nestedaov_ms_bylevel_byposmut <- sapply(inherited_noctrl_altasin_nestedaov_bylevel_byposmut, function(X) { structure(c(X[[1]][[1]][, "Mean Sq"], X[[2]][[1]][, "Mean Sq"]), names = c("mousewise", "cellwise", "mitowise")) })

inherited_noctrl_altasin_nestedaov_ftest_mousetocell_byposmut <- sapply(inherited_posmuts, function(x) {
    df <- inherited_noctrl_altasin_nestedaov_df_bylevel_byposmut[c("mousewise", "cellwise"), x]
    ms <- inherited_noctrl_altasin_nestedaov_ms_bylevel_byposmut[c("mousewise", "cellwise"), x]
    f <- ms["mousewise"] / ms["cellwise"]
    p <- pf(f, df1 = df["mousewise"], df2 = df["cellwise"], lower.tail = FALSE)
    c(fstat = unname(f), pval = unname(p))
})

inherited_noctrl_altasin_nestedaov_ftest_celltomito_byposmut <- sapply(inherited_posmuts, function(x) {
    df <- inherited_noctrl_altasin_nestedaov_df_bylevel_byposmut[c("cellwise", "mitowise"), x]
    ms <- inherited_noctrl_altasin_nestedaov_ms_bylevel_byposmut[c("cellwise", "mitowise"), x]
    f <- ms["cellwise"] / ms["mitowise"]
    p <- pf(f, df1 = df["cellwise"], df2 = df["mitowise"], lower.tail = FALSE)
    c(fstat = unname(f), pval = unname(p))
})

inherited_noctrl_altasin_nestedaov_msprop_bylevel_byposmut <- t(t(inherited_noctrl_altasin_nestedaov_ms_bylevel_byposmut) / colSums(inherited_noctrl_altasin_nestedaov_ms_bylevel_byposmut))

inherited_noctrl_altasin_nestedaov_ftest_bylevel_byposmut <- t(rbind(
   local({rownames(inherited_noctrl_altasin_nestedaov_df_bylevel_byposmut) <- paste0("df_", rownames(inherited_noctrl_altasin_nestedaov_df_bylevel_byposmut)); inherited_noctrl_altasin_nestedaov_df_bylevel_byposmut }), 
   local({rownames(inherited_noctrl_altasin_nestedaov_ms_bylevel_byposmut) <- paste0("ms_", rownames(inherited_noctrl_altasin_nestedaov_ms_bylevel_byposmut)); inherited_noctrl_altasin_nestedaov_ms_bylevel_byposmut }), 
   local({rownames(inherited_noctrl_altasin_nestedaov_msprop_bylevel_byposmut) <- paste0("msprop_", rownames(inherited_noctrl_altasin_nestedaov_msprop_bylevel_byposmut)); inherited_noctrl_altasin_nestedaov_msprop_bylevel_byposmut }),
   local({rownames(inherited_noctrl_altasin_nestedaov_ftest_celltomito_byposmut) <- paste0(rownames(inherited_noctrl_altasin_nestedaov_ftest_celltomito_byposmut), "_celltomito"); inherited_noctrl_altasin_nestedaov_ftest_celltomito_byposmut }) , 
   local({rownames(inherited_noctrl_altasin_nestedaov_ftest_mousetocell_byposmut) <- paste0(rownames(inherited_noctrl_altasin_nestedaov_ftest_mousetocell_byposmut), "_mousetocell"); inherited_noctrl_altasin_nestedaov_ftest_mousetocell_byposmut }) 
))

inherited_noctrl_altasin_nestedaov_ftest_bylevel_byposmut <- data.table(
    posmut = inherited_posmuts, 
    inherited_noctrl_altasin_nestedaov_ftest_bylevel_byposmut
)
fwrite(inherited_noctrl_altasin_nestedaov_ftest_bylevel_byposmut, file = "Report/release/SNVs/hierarchy/inherited_noctrl_altasin_nestedaov_ftest_bylevel_byposmut.csv")
inherited_noctrl_altasin_nestedaov_ftest_bylevel_byposmut <- fread(file = "Report/release/SNVs/hierarchy/inherited_noctrl_altasin_nestedaov_ftest_bylevel_byposmut.csv")

inherited_noctrl_altasin_nestedaov_df_bylevel_byposmut_df <- reshape2::melt(inherited_noctrl_altasin_nestedaov_df_bylevel_byposmut)
inherited_noctrl_altasin_nestedaov_ms_bylevel_byposmut_df <- reshape2::melt(inherited_noctrl_altasin_nestedaov_ms_bylevel_byposmut)
inherited_noctrl_altasin_nestedaov_msprop_bylevel_byposmut_df <- reshape2::melt(inherited_noctrl_altasin_nestedaov_msprop_bylevel_byposmut)

inherited_noctrl_altasin_nestedaov_bylevel_byposmut_dt <- data.table(
    posmut = inherited_noctrl_altasin_nestedaov_df_bylevel_byposmut_df[, 2],
    level = inherited_noctrl_altasin_nestedaov_df_bylevel_byposmut_df[, 1],
    df = inherited_noctrl_altasin_nestedaov_df_bylevel_byposmut_df[, 3], 
    ms = inherited_noctrl_altasin_nestedaov_ms_bylevel_byposmut_df[, 3],
    ms_prop = inherited_noctrl_altasin_nestedaov_msprop_bylevel_byposmut_df[, 3]
)

fig <- ggplot(inherited_noctrl_altasin_nestedaov_bylevel_byposmut_dt, aes(x = posmut, fill = level, y = 100 * ms_prop)) + geom_bar(stat = "identity") + scale_x_discrete(limits = inherited_noctrl_altasin_nestedaov_bylevel_byposmut_dt[level == "mousewise"][order(-ms_prop)][, posmut]) + theme_classic(base_size = 12) + theme(axis.text.x = element_text(size = 4, angle = 90, vjust = 0.5, hjust = 1)) + xlab("") + ylab("% mean squares")
ggsave(fig, file = "Report/release/SNVs/hierarchy/inherited_noctrl_altasin_nestedaov_ms_bylevel_byposmut_bar.pdf", width = 12, height = 4)

pdf("Report/release/SNVs/hierarchy/inherited_noctrl_altasin_nestedaov_msprop_bylevel_hist.pdf", width = 12, height = 6)
par(ps = 16, lend = 2, ljoin = 1, bty = "L", mfrow = c(1, 2), mar = c(2.5, 2.5, 1, 0), oma = c(0, 0, 0, 0), mgp = c(1.5, 0.5, 0))
hist(log2(inherited_noctrl_altasin_nestedaov_msprop_bylevel_byposmut["mousewise", ] / inherited_noctrl_altasin_nestedaov_msprop_bylevel_byposmut["cellwise", ]), main = "between-mice to between-cells", xlab = "log2 ratios of mean squares", nclass = 100, border = FALSE)
hist(log2(inherited_noctrl_altasin_nestedaov_msprop_bylevel_byposmut["cellwise", ] / inherited_noctrl_altasin_nestedaov_msprop_bylevel_byposmut["mitowise", ]), main = "between-cells to between-mitos", xlab = "log2 ratios of mean squares", nclass = 100, border = FALSE)
dev.off()

pdf(file = "Report/release/SNVs/hierarchy/inherited_noctrl_altasin_bycell_byposmut_boxplot.pdf", width = 11, height = 4)
for (posmut in inherited_noctrl_altasin_nestedaov_bylevel_byposmut_dt[level == "mitowise"][order(ms_prop)][, posmut]) {
    message(posmut)
    print(ggplot(inherited_noctrl_altasin_bymito_byposmut[MouseID != "Mouse16&17", c("MouseID", "CellUID", posmut), with = FALSE], aes(y = get(posmut), x = CellUID)) + geom_boxplot(outlier.color = NA) + geom_jitter(aes(col = MouseID), size = 0.5) + theme_classic(base_size = 12) + theme(axis.text.x = element_text(size = 8, angle = 90, vjust = 0.5, hjust = 1)) + xlab("") + ylab("arcsin AF") + ggtitle(posmut))
}
dev.off()

pdf("Report/release/SNVs/hierarchy/inherited_noctrl_altasin_bymouse_byposmut_boxplot.pdf", width = 4.5, height = 4.5, useDingbats = FALSE)
for (posmut in inherited_noctrl_altasin_nestedaov_bylevel_byposmut_dt[level == "cellwise"][order(ms_prop)][, posmut]) {
    message(posmut)
    print(ggplot(inherited_noctrl_altasin_bycell_byposmut[MouseID != "Mouse16&17"], aes(x = MouseID, y = get(posmut), fill = MouseID)) + geom_boxplot(outlier.color = NA) + geom_jitter(size = 0.6) + theme_classic(base_size = 12) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + xlab("") + ylab("arcsin AF") + ggtitle(posmut))
}
dev.off()

pdf("Report/release/SNVs/hierarchy/inherited_noctrl_altasinctrd_bycell_byposmut_boxplot.pdf", width = 11, height = 4, useDingbats = FALSE)
for (posmut in inherited_noctrl_altasin_nestedaov_bylevel_byposmut_dt[level == "mitowise"][order(ms_prop)][, posmut]) {
    message(posmut)
    print(ggplot(inherited_noctrl_altasinctrd_bymito_byposmut[MouseID != "Mouse16&17", c("MouseID", "CellUID", posmut), with = FALSE], aes(y = get(posmut), x = CellUID)) + geom_boxplot(outlier.color = NA) + geom_jitter(aes(col = MouseID), size = 0.5) + theme_classic(base_size = 12) + theme(axis.text.x = element_text(size = 8, angle = 90, vjust = 0.5, hjust = 1)) + xlab("") + ylab("arcsin centered AF") + ggtitle(posmut))
}
dev.off()

pdf("Report/release/SNVs/hierarchy/inherited_noctrl_altasinctrd_bymouse_byposmut_boxplot.pdf", width = 4.5, height = 4.5, useDingbats = FALSE)
for (posmut in inherited_noctrl_altasin_nestedaov_bylevel_byposmut_dt[level == "cellwise"][order(ms_prop)][, posmut]) {
    message(posmut)
    print(ggplot(inherited_noctrl_altasinctrd_bycell_byposmut[MouseID != "Mouse16&17"], aes(x = MouseID, y = get(posmut), fill = MouseID)) + geom_boxplot(outlier.color = NA) + geom_jitter(size = 0.6) + theme_classic(base_size=12) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + xlab("") + ylab("arcsin centered AF") + ggtitle(posmut))
}
dev.off()

pdf("Report/release/SNVs/hierarchy/inherited_noctrl_altasinctrd_allmice_byposmut_boxplot.pdf", width = 2.5, height = 4.5, useDingbats = FALSE)
for (posmut in inherited_noctrl_altasin_nestedaov_bylevel_byposmut_dt[level == "cellwise"][order(ms_prop)][, posmut]) {
    message(posmut)
    print(ggplot(inherited_noctrl_altasinctrd_bymouse_byposmut[MouseID != "Mouse16&17"], aes(x = 1, y = get(posmut))) + geom_boxplot(outlier.color = NA) + geom_jitter(size = 0.6, aes(col = MouseID)) + theme_classic(base_size=12) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + xlab("") + ylab("arcsin centered AF") + ggtitle(posmut))
}
dev.off()

###########################################################################
## In order to rule out the possibility that the variation proportion is
## driven by the <5% AF noise, we went to shuffle the data and assess the 
## null distribution of the mouse-wise, cell-wise and mito-wise variability
## proportion. The outcome is indeed 1/3, 1/3, 1/3. 
###########################################################################
## Let's run 1000 times permutation and see if the null distribution is 1/3, 1/3, 1/3.
set.seed(2022)
rnd_altasin_bymito_byposmut <- replicate(1000, {
    rnd <- copy(inherited_noctrl_altasin_bymito_byposmut)
    idx <- sample(seq(nrow(rnd)))
    rnd[, MouseID := MouseID[idx]]
    rnd[, CellUID := CellUID[idx]]
    rnd
}, simplify = FALSE)

rnd_altasin_nestedaov_bylevel_byposmut <- mclapply(1:1000, function(i) {
    message(i)
    sapply(inherited_posmuts, function(posmut) {
        X <- rnd_altasin_bymito_byposmut[[i]][MouseID != "Mouse16&17", c(posmut, "MouseID", "CellUID"), with = FALSE]
        names(X)[1] <- "VAF"
        mod <- aov(VAF ~ MouseID + Error(CellUID), data = X)
        summary(mod)
    }, simplify = FALSE)
}, mc.cores = 50)

rnd_altasin_nestedaov_df_bylevel_byposmut <- lapply(1:1000, function(i) { sapply(rnd_altasin_nestedaov_bylevel_byposmut[[i]], function(X) { structure(c(X[[1]][[1]][, "Df"], X[[2]][[1]][, "Df"]), names = c("mousewise", "cellwise", "mitowise")) }) })
rnd_altasin_nestedaov_ms_bylevel_byposmut <- lapply(1:1000, function(i) { sapply(rnd_altasin_nestedaov_bylevel_byposmut[[i]], function(X) { structure(c(X[[1]][[1]][, "Mean Sq"], X[[2]][[1]][, "Mean Sq"]), names = c("mousewise", "cellwise", "mitowise")) }) })

rnd_altasin_nestedaov_ftest_mousetocell_byposmut <- lapply(1:1000, function(i) { 
    sapply(inherited_posmuts, function(x) {
        df <- rnd_altasin_nestedaov_df_bylevel_byposmut[[i]][c("mousewise", "cellwise"), x]
        ms <- rnd_altasin_nestedaov_ms_bylevel_byposmut[[i]][c("mousewise", "cellwise"), x]
        f <- ms["mousewise"] / ms["cellwise"]
        p <- pf(f, df1 = df["mousewise"], df2 = df["cellwise"], lower.tail = FALSE)
        c(fstat = unname(f), pval = unname(p))
    })
})
rnd_altasin_nestedaov_ftest_celltomito_byposmut <- lapply(1:1000, function(i) {
    sapply(inherited_posmuts, function(x) {
        df <- rnd_altasin_nestedaov_df_bylevel_byposmut[[i]][c("cellwise", "mitowise"), x]
        ms <- rnd_altasin_nestedaov_ms_bylevel_byposmut[[i]][c("cellwise", "mitowise"), x]
        f <- ms["cellwise"] / ms["mitowise"]
        p <- pf(f, df1 = df["cellwise"], df2 = df["mitowise"], lower.tail = FALSE)
        c(fstat = unname(f), pval = unname(p))
    })
})

rnd_altasin_nestedaov_msprop_bylevel_byposmut <- lapply(1:1000, function(i) { t(t(rnd_altasin_nestedaov_ms_bylevel_byposmut[[i]]) / colSums(rnd_altasin_nestedaov_ms_bylevel_byposmut[[i]])) })

rnd_altasin_nestedaov_ftest_bylevel_byposmut <- lapply(1:1000, function(i) {
    t(rbind(
        local({rownames(rnd_altasin_nestedaov_df_bylevel_byposmut[[i]]) <- paste0("df_", rownames(rnd_altasin_nestedaov_df_bylevel_byposmut[[i]])); rnd_altasin_nestedaov_df_bylevel_byposmut[[i]] }), 
        local({rownames(rnd_altasin_nestedaov_ms_bylevel_byposmut[[i]]) <- paste0("ms_", rownames(rnd_altasin_nestedaov_ms_bylevel_byposmut[[i]])); rnd_altasin_nestedaov_ms_bylevel_byposmut[[i]] }), 
        local({rownames(rnd_altasin_nestedaov_msprop_bylevel_byposmut[[i]]) <- paste0("msprop_", rownames(rnd_altasin_nestedaov_msprop_bylevel_byposmut[[i]])); rnd_altasin_nestedaov_msprop_bylevel_byposmut[[i]] }),
        local({rownames(rnd_altasin_nestedaov_ftest_celltomito_byposmut[[i]]) <- paste0(rownames(rnd_altasin_nestedaov_ftest_celltomito_byposmut[[i]]), "_celltomito"); rnd_altasin_nestedaov_ftest_celltomito_byposmut[[i]] }) , 
        local({rownames(rnd_altasin_nestedaov_ftest_mousetocell_byposmut[[i]]) <- paste0(rownames(rnd_altasin_nestedaov_ftest_mousetocell_byposmut[[i]]), "_mousetocell"); rnd_altasin_nestedaov_ftest_mousetocell_byposmut[[i]]}) 
    ))
})

rnd_altasin_nestedaov_ftest_bylevel_byposmut <- lapply(1:1000, function(i) {
    data.table(
        perm = i,
        posmut = inherited_posmuts, 
        rnd_altasin_nestedaov_ftest_bylevel_byposmut[[i]]
    )
})
rnd_altasin_nestedaov_ftest_bylevel_byposmut <- rbindlist(rnd_altasin_nestedaov_ftest_bylevel_byposmut)
fwrite(rnd_altasin_nestedaov_ftest_bylevel_byposmut, file = "Report/release/SNVs/hierarchy/rnd_altasin_nestedaov_ftest_bylevel_byposmut.csv.gz")
rnd_altasin_nestedaov_ftest_bylevel_byposmut <- fread(file = "Report/release/SNVs/hierarchy/rnd_altasin_nestedaov_ftest_bylevel_byposmut.csv.gz")

pdf(file = "Report/release/SNVs/hierarchy/rnd_altasin_nestedaov_ftest_bylevel_density.pdf", width = 10, height = 4)
par(ps = 16, lend = 2, ljoin = 1, bty = "L", mfrow = c(1, 3), mar = c(2.5, 2.5, 1, 0), oma = c(0, 0, 0, 0), mgp = c(1.5, 0.5, 0))
inherited_noctrl_altasin_nestedaov_ftest_bylevel_byposmut[, plot(density(100 * msprop_mousewise), col = scales::hue_pal()(3)[1], lwd = 3, xlim = c(-0.1, 1) * 100, ylim = c(0, 7)/100, main = "between mice", xlab = "mean squares proportion (%)", ylab = "density")]
rnd_altasin_nestedaov_ftest_bylevel_byposmut[, lines(density(100 * msprop_mousewise), col= "#00000001"), keyby = "perm"]
inherited_noctrl_altasin_nestedaov_ftest_bylevel_byposmut[, lines(density(100 * msprop_mousewise), col = scales::hue_pal()(3)[1], lwd = 3)]

inherited_noctrl_altasin_nestedaov_ftest_bylevel_byposmut[, plot(density(100 * msprop_cellwise), col = scales::hue_pal()(3)[2], lwd = 3, xlim = c(-0.1, 1) * 100, ylim = c(0, 7)/100, main = "between cells", xlab = "mean squares proportion (%)", ylab = "density")]
rnd_altasin_nestedaov_ftest_bylevel_byposmut[, lines(density(100 * msprop_cellwise), col = "#00000001"), keyby = "perm"]
inherited_noctrl_altasin_nestedaov_ftest_bylevel_byposmut[, lines(density(100 * msprop_cellwise), col = scales::hue_pal()(3)[2], lwd = 3)]

inherited_noctrl_altasin_nestedaov_ftest_bylevel_byposmut[, plot(density(100 * msprop_mitowise), col = scales::hue_pal()(3)[3], lwd = 3, xlim = c(-0.1, 1) * 100, ylim = c(0, 7)/100, main = "between mitos", xlab = "mean squares proportion (%)", ylab = "density")]
rnd_altasin_nestedaov_ftest_bylevel_byposmut[, lines(density(100 * msprop_mitowise), col= "#00000001"), keyby = "perm"]
inherited_noctrl_altasin_nestedaov_ftest_bylevel_byposmut[, lines(density(100 * msprop_mitowise), col = scales::hue_pal()(3)[3], lwd = 3)]
dev.off()
