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
library("uwot")
library("RColorBrewer")

inherited_altasin_bymito_byposmut <- fread(file = "Report/SNVs/hierarchy/inherited_noctrl_altasin_bymito_byposmut.csv")
inherited_altasin_bymito_byposmut <- inherited_altasin_bymito_byposmut[IsCtrl == "N" & MouseID != "Mouse16&17"]


inherited_altasin_mat <- inherited_altasin_bymito_byposmut[, -(1:20)][, apply(.SD, 2, as.numeric)]
rownames(inherited_altasin_mat) <- inherited_altasin_bymito_byposmut[, LibraryMitoID]
colnames(inherited_altasin_mat) <- colnames(inherited_altasin_bymito_byposmut)[-(1:20)]
inherited_altasin_mat[is.na(inherited_altasin_mat)] <- 0
inherited_altasin_umap <- umap(prcomp(inherited_altasin_mat)$x[, 1:5], metric = "euclidean", n_components = 2, n_neighbors = 15, min_dist = 0.5, spread = 1, verbose = TRUE, seed = 42, n_epochs = 500)
pdf("Report/SNVs/umap/inherited_altasin_umap.pdf", width = 12, height = 6)
par(ps = 12, lend = 2, ljoin = 1, bty = "L", mfrow = c(1, 2), mar = c(2.5, 2.5, 1, 0), oma = c(0, 0, 0, 0), mgp = c(1.5, 0.5, 0))
plot(inherited_altasin_umap, col = colorRampPalette(brewer.pal(n = 9, "Set1"))(inherited_altasin_bymito_byposmut[, uniqueN(MouseID)])[factor(inherited_altasin_bymito_byposmut[, MouseID])], pch = 19, cex = 1, xlab = "UMAP 1", ylab = "UMAP 2")
legend("topright", col = colorRampPalette(brewer.pal(n = 9, "Set1"))(inherited_altasin_bymito_byposmut[, uniqueN(MouseID)]), legend = levels(factor(inherited_altasin_bymito_byposmut[, MouseID])), pch = 19, bty = "n", cex = 0.8)
plot(inherited_altasin_umap, col = colorRampPalette(rainbow(9))(inherited_altasin_bymito_byposmut[, uniqueN(CellUID)])[factor(inherited_altasin_bymito_byposmut[, CellUID])], pch = 19, cex = 1, xlab = "UMAP 1", ylab = "UMAP 2")
legend("topright", col = colorRampPalette(rainbow(9))(inherited_altasin_bymito_byposmut[, uniqueN(CellUID)]), legend = levels(factor(inherited_altasin_bymito_byposmut[, CellUID])), pch = 19, bty = "n", cex = 0.3, ncol = 5)
dev.off()


inherited_altasin_mat <- inherited_altasin_bymito_byposmut[, -(1:20)][, apply(.SD, 2, as.numeric)]
rownames(inherited_altasin_mat) <- inherited_altasin_bymito_byposmut[, LibraryMitoID]
colnames(inherited_altasin_mat) <- colnames(inherited_altasin_bymito_byposmut)[-(1:20)]
inherited_altasin_mat <- apply(inherited_altasin_mat, 2, function(x) replace(x, is.na(x), mean(x, na.rm = TRUE)))
inherited_altasin_umap <- umap(prcomp(inherited_altasin_mat)$x[, 1:5], metric = "euclidean", n_components = 2, n_neighbors = 15, min_dist = 0.5, spread = 1, verbose = TRUE, seed = 42, n_epochs = 500)
pdf("Report/SNVs/umap/inherited_altasin_imputed-mean_umap.pdf", width = 12, height = 6)
par(ps = 12, lend = 2, ljoin = 1, bty = "L", mfrow = c(1, 2), mar = c(2.5, 2.5, 1, 0), oma = c(0, 0, 0, 0), mgp = c(1.5, 0.5, 0))
plot(inherited_altasin_umap, col = colorRampPalette(brewer.pal(n = 9, "Set1"))(inherited_altasin_bymito_byposmut[, uniqueN(MouseID)])[factor(inherited_altasin_bymito_byposmut[, MouseID])], pch = 19, cex = 1, xlab = "UMAP 1", ylab = "UMAP 2")
legend("topright", col = colorRampPalette(brewer.pal(n = 9, "Set1"))(inherited_altasin_bymito_byposmut[, uniqueN(MouseID)]), legend = levels(factor(inherited_altasin_bymito_byposmut[, MouseID])), pch = 19, bty = "n", cex = 0.9)
plot(inherited_altasin_umap, col = colorRampPalette(rainbow(9))(inherited_altasin_bymito_byposmut[, uniqueN(CellUID)])[factor(inherited_altasin_bymito_byposmut[, CellUID])], pch = 19, cex = 1, xlab = "UMAP 1", ylab = "UMAP 2")
legend("topleft", col = colorRampPalette(rainbow(9))(inherited_altasin_bymito_byposmut[, uniqueN(CellUID)]), legend = levels(factor(inherited_altasin_bymito_byposmut[, CellUID])), pch = 19, bty = "n", cex = 0.3, ncol = 2)
dev.off()


inherited_altasin_mat <- inherited_altasin_bymito_byposmut[, -(1:20)][, apply(.SD, 2, as.numeric)]
inherited_altasin_meanbymouse <- do.call(rbind, by(inherited_altasin_mat, INDICES = inherited_altasin_bymito_byposmut[, MouseID], FUN = function(X) colMeans(X, na.rm = TRUE)))
inherited_altasin_meanbymouse <- apply(inherited_altasin_meanbymouse, 2, function(x) replace(x, is.na(x), mean(x, na.rm = TRUE)))
inherited_altasin_mat <- inherited_altasin_bymito_byposmut[, -(1:20)][, apply(.SD, 2, as.numeric)]
rownames(inherited_altasin_mat) <- inherited_altasin_bymito_byposmut[, LibraryMitoID]
colnames(inherited_altasin_mat) <- colnames(inherited_altasin_bymito_byposmut)[-(1:20)]
inherited_altasin_mat <- sapply(colnames(inherited_altasin_mat), function(s) {
    message(s)
    sapply(rownames(inherited_altasin_mat), function(m) {
        x <- inherited_altasin_mat[m, s]
        if (is.na(x)) {
            k <- inherited_altasin_bymito_byposmut[LibraryMitoID == m, MouseID]
            inherited_altasin_meanbymouse[k, s]
        } else {
            x
        }
    })
})
inherited_altasin_umap <- umap(prcomp(inherited_altasin_mat)$x[, 1:5], metric = "euclidean", n_components = 2, n_neighbors = 15, min_dist = 0.4, spread = 1, verbose = TRUE, seed = 42, n_epochs = 500)
pdf("Report/SNVs/umap/inherited_altasin_imputed-meanbymouse_umap.pdf", width = 12, height = 6)
par(ps = 12, lend = 2, ljoin = 1, bty = "L", mfrow = c(1, 2), mar = c(2.5, 2.5, 1, 0), oma = c(0, 0, 0, 0), mgp = c(1.5, 0.5, 0))
plot(inherited_altasin_umap, col = colorRampPalette(brewer.pal(n = 9, "Set1"))(inherited_altasin_bymito_byposmut[, uniqueN(MouseID)])[factor(inherited_altasin_bymito_byposmut[, MouseID])], pch = 19, cex = 1, xlab = "UMAP 1", ylab = "UMAP 2")
legend("topleft", col = colorRampPalette(brewer.pal(n = 9, "Set1"))(inherited_altasin_bymito_byposmut[, uniqueN(MouseID)]), legend = levels(factor(inherited_altasin_bymito_byposmut[, MouseID])), pch = 19, bty = "n", cex = 0.9)
plot(inherited_altasin_umap, col = colorRampPalette(rainbow(9))(inherited_altasin_bymito_byposmut[, uniqueN(CellUID)])[factor(inherited_altasin_bymito_byposmut[, CellUID])], pch = 19, cex = 1, xlab = "UMAP 1", ylab = "UMAP 2")
legend("topleft", col = colorRampPalette(rainbow(9))(inherited_altasin_bymito_byposmut[, uniqueN(CellUID)]), legend = levels(factor(inherited_altasin_bymito_byposmut[, CellUID])), pch = 19, bty = "n", cex = 0.3, ncol = 4)
dev.off()
