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
library("ggrepel")
library("lme4")
library("car")

mito_barcodes <- read.csv("Data/mito_barcodes.csv", as.is = TRUE)
mitoIDs <- mito_barcodes[, "ID"] 

snv_info <- read.csv("Data/snv_loci_v2.csv", as.is = TRUE)
snvIDs <- snv_info[,"SNVID"]

chrmbases <- readLines("Data/mm10.mito/chrM.fa")[-1]
chrmbases <- paste0(chrmbases, collapse = "")
chrmbases <- strsplit(chrmbases, "")[[1]]
nchrmbases <- length(chrmbases)

MitoInfo <- fread(file = "Report/metadata/MitoInfo.csv")
MitoInfo[, ExptID := factor(ExptID)]
MitoInfo[, MitoID := factor(MitoID, levels = mitoIDs)]

CellInfo <- fread("Report/metadata/CellInfo.csv")
dim(CellInfo)
## [1] 102  12
MouseInfo <- fread("Report/metadata/MouseInfo.csv")
dim(MouseInfo)

support_byposmut <- fread("Report/SNVs/filter/basediffperc_cutdemux_sub500k_q30_unstranded_highdepth_highaf_qcfltd_support_byposmut.csv")
inherited_posmuts <- support_byposmut[nmice >= 3, posmut]
somatic_posmuts <- support_byposmut[nmice == 1, posmut]

###########################################################################
## ANOVA to test whether cell type contribute to the # of SNV sites variance
###########################################################################
noctrl_nposhasdata_bymito <- fread("Report/SNVs/QC/basediffperc_cutdemux_q30_unstranded_highdepth_noctrl_nposhasdata_bymito.csv")
noctrl_nposhasdata_bymito[is.na(nposhasdata), nposhasdata := 0]
noctrl_npos_bymito <- fread("Report/SNVs/QC/basediffperc_cutdemux_q30_unstranded_highdepth_highaf_noctrl_npos_bymito.csv")
noctrl_npos_bymito <- merge.data.table(noctrl_npos_bymito, noctrl_nposhasdata_bymito[, c("LibraryMitoID", "nposhasdata")], on = "LibraryMitoID")

## Are total # of sites differential?
## ANOVA with MouseID as random effects, using # bases with non-missing data instead of the barcode
## We first tried to use CellUID as the other random effect i.e. (1 | MouseID/CellUID) as the random term in
## the formula but that resulted in "boundary (singular) fit: see help('isSingular')" for some cases, possibly
## becaused of the small sample size for some cells. 
car::Anova(mod <- lme4::lmer(log1p(npos) ~ CellType + (1 | MouseID) + log1p(nposhasdata), data = noctrl_npos_bymito[MouseID != "Mouse16&17"]), type = 3, test.statistic = "Chisq")
## Analysis of Deviance Table (Type III Wald chisquare tests)
## 
## Response: log1p(npos)
##                        Chisq Df Pr(>Chisq)
## (Intercept)         232.2444  1     <2e-16 ***
## CellType              1.5848  1     0.2081
## log1p(nposhasdata) 1161.0929  1     <2e-16 ***
## ---
## Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
 
###########################################################################
## Let's skip ANOVA but just use the mixed-effect linear model coefficients.
###########################################################################
mod <- lme4::lmer(log1p(npos) ~ CellType + (1 | MouseID) + log1p(nposhasdata), data = noctrl_npos_bymito[MouseID != "Mouse16&17"])
summary(mod)[[10]]["CellTypeNeuron", "Estimate"]
## [1] -0.04986401

set.seed(2021)
rnd_npos_ls <- replicate(1000, { X <- copy(noctrl_npos_bymito[MouseID != "Mouse16&17"]); X[, CellType := sample(CellType)]; X }, simplify = FALSE)
rnd_npos_lme <- lapply(seq_along(rnd_npos_ls), function(i) {
    message(i)
    lme4::lmer(log1p(npos) ~ CellType + (1 | MouseID) + log1p(nposhasdata), data = rnd_npos_ls[[i]])
})
rnd_npos_lme_celltypestats <- sapply(rnd_npos_lme, function(X) summary(X)[[10]]["CellTypeNeuron", "Estimate"])

pdf("Report/SNVs/diff/noctrl_rndcelltype_npos_lme_celltypestats_hist.pdf", width = 6, height = 6)
par(ps = 16, lend = 2, ljoin = 1, bty = "L", mfrow = c(1, 1), mar = c(2.5, 2.5, 1, 0.5), oma = c(0, 0, 0, 0), mgp = c(1.5, 0.5, 0))
hist(rnd_npos_lme_celltypestats, main = sprintf("p = %.3f", if (summary(mod)[[10]]["CellTypeNeuron", "Estimate"] < 0) { mean(rnd_npos_lme_celltypestats < summary(mod)[[10]]["CellTypeNeuron", "Estimate"]) } else { mean(rnd_npos_lme_celltypestats > summary(mod)[[10]]["CellTypeNeuron", "Estimate"]) }), xlab = "linear mixed-effects model cell-type coef.", ylab = "# permutations", border = "grey70", nclass = 1000)
abline(v = summary(mod)[[10]]["CellTypeNeuron", "Estimate"], lty = 2, lwd = 3, col = "red")
dev.off()

###########################################################################
## Use `hassnv` as a binary indicator and do a logistic regression to see
## if one cell type is more likely to have more # of SNVs.
###########################################################################
inherited_noctrl_hassnv_bymito_byposmut <- fread("Report/SNVs/origin/highdepth_highaf_inherited_noctrl_hassnv_bymito_byposmut.csv")
inherited_noctrl_hasdata_bymito_byposmut <- fread("Report/SNVs/origin/highdepth_inherited_noctrl_hasdata_bymito_byposmut.csv")
inherited_noctrl_hassnv_bymito_byposmut <- inherited_noctrl_hassnv_bymito_byposmut[MitoInfo[IsCtrl == "N", c("LibraryMitoID")], on = "LibraryMitoID"]
inherited_noctrl_hasdata_bymito_byposmut <- inherited_noctrl_hasdata_bymito_byposmut[MitoInfo[IsCtrl == "N", c("LibraryMitoID")], on = "LibraryMitoID"]
inherited_noctrl_hassnv_bymito_byposmut <- data.table(MitoInfo[IsCtrl == "N", 1:19], inherited_noctrl_hassnv_bymito_byposmut[, -c(1:19)])
inherited_noctrl_hasdata_bymito_byposmut <- data.table(MitoInfo[IsCtrl == "N", 1:19], inherited_noctrl_hasdata_bymito_byposmut[, -c(1:19)])
setnafill(inherited_noctrl_hassnv_bymito_byposmut, fill = 0, cols = names(inherited_noctrl_hassnv_bymito_byposmut)[-c(1:19)])
setnafill(inherited_noctrl_hasdata_bymito_byposmut, fill = 0, cols = names(inherited_noctrl_hasdata_bymito_byposmut)[-c(1:19)])

## Since we now condition the model on the non-missing data only, we can ignore the barcode effect, which was shown to impact the missing data rate. 
inherited_noctrl_hassnv_logit_byposmut <- sapply(inherited_posmuts, function(posmut) {
    message(posmut)
    X <- inherited_noctrl_hassnv_bymito_byposmut[MouseID != "Mouse16&17", c(posmut, "MouseID", "CellType"), with = FALSE] 
    names(X)[1] <- "hassnv"
    Y <- inherited_noctrl_hasdata_bymito_byposmut[MouseID != "Mouse16&17", c(posmut), with = FALSE]
    names(Y)[1] <- "hasdata"
    idx <- Y[, hasdata] == 1
    X <- X[idx]
    tryCatch( {
        mod <- lme4::glmer(hassnv ~ CellType + (1 | MouseID), data = X, family = "binomial")
        summary(mod)[[10]]["CellTypeNeuron", ]
    }, error = function(e) {
        structure(rep(NA, 4), names = c("Estimate",  "Std. Error",  "z value",  "Pr(>|z|)"))
        ## Watch out for nonconverging estimations! To be safe, we simply don't trust any of them.
    }, warning = function(w) {
        structure(rep(NA, 4), names = c("Estimate",  "Std. Error",  "z value",  "Pr(>|z|)"))
    })
}, simplify = FALSE)

inherited_noctrl_hassnv_logit_byposmut <- do.call(rbind, inherited_noctrl_hassnv_logit_byposmut)
inherited_noctrl_hassnv_logit_byposmut <- data.table(posmut = inherited_posmuts, inherited_noctrl_hassnv_logit_byposmut)
inherited_noctrl_hassnv_logit_byposmut[, padj := p.adjust(`Pr(>|z|)`)]
inherited_noctrl_hassnv_logit_byposmut <- inherited_noctrl_hassnv_logit_byposmut[order(`Pr(>|z|)`)]

fwrite(inherited_noctrl_hassnv_logit_byposmut, file = "Report/SNVs/diff/inherited_noctrl_hassnv_logit_byposmut.csv")
inherited_noctrl_hassnv_logit_byposmut <- fread(file = "Report/SNVs/diff/inherited_noctrl_hassnv_logit_byposmut.csv")

###########################################################################
## Comparing the fraction of mitos with SNV per mito between astrocytes and
## neurons. We set a limit to include cells with enough number of mitos. 
###########################################################################
inherited_noctrl_fmitoshassnv_bycell_byposmut <- fread("Report/SNVs/origin/highdepth_highaf_inherited_noctrl_fmitoshassnv_bycell_byposmut.csv")
inherited_noctrl_nmitoshasdata_bycell_byposmut <- fread(file = "Report/SNVs/origin/highdepth_inherited_noctrl_nmitoshasdata_bycell_byposmut.csv")
inherited_noctrl_fmitoshassnv_bycell_byposmut <- inherited_noctrl_fmitoshassnv_bycell_byposmut[CellInfo[, "CellUID"], on = "CellUID"]
inherited_noctrl_nmitoshasdata_bycell_byposmut <- inherited_noctrl_nmitoshasdata_bycell_byposmut[CellInfo[, "CellUID"], on = "CellUID"]

nmitoshasdata_th <- 3
inherited_noctrl_fmitoshassnv_nolowmito_logit_byposmut <- sapply(inherited_posmuts, function(posmut) {
    X <- inherited_noctrl_nmitoshasdata_bycell_byposmut[MouseID != "Mouse16&17"]
    Y <- inherited_noctrl_fmitoshassnv_bycell_byposmut[MouseID != "Mouse16&17"]
    x <- X[, get(posmut)] 
    idx <- x >= nmitoshasdata_th
    Y <- Y[idx, c(posmut, "MouseID", "CellType", "CellUID"), with = FALSE]
    x <- x[idx]
    names(Y)[1] <- "fmitoshassnv"
    Y[, nmitoshasdata := x]
    tryCatch( {
    ## https://stats.stackexchange.com/questions/233366/how-to-fit-a-mixed-model-with-response-variable-between-0-and-1
    ## We need to watch out for failed SNVs. Just set them to NA. 
        res <- lme4::glmer(fmitoshassnv ~ CellType + (1 | MouseID), data = Y, family = "binomial", weight = nmitoshasdata)
        summary(res)[[10]]["CellTypeNeuron", ]
    }, error = function(e) { 
        structure(rep(NA, 4), names = c("Estimate",  "Std. Error",  "z value",  "Pr(>|z|)")) 
    }, warning = function(w) {
        structure(rep(NA, 4), names = c("Estimate",  "Std. Error",  "z value",  "Pr(>|z|)")) 
    })
}, simplify = FALSE)
inherited_noctrl_fmitoshassnv_nolowmito_logit_byposmut <- do.call(rbind, inherited_noctrl_fmitoshassnv_nolowmito_logit_byposmut)
inherited_noctrl_fmitoshassnv_nolowmito_logit_byposmut <- data.table(posmut = rownames(inherited_noctrl_fmitoshassnv_nolowmito_logit_byposmut), inherited_noctrl_fmitoshassnv_nolowmito_logit_byposmut)
inherited_noctrl_fmitoshassnv_nolowmito_logit_byposmut[, padj := p.adjust(`Pr(>|z|)`)]
inherited_noctrl_fmitoshassnv_nolowmito_logit_byposmut <- inherited_noctrl_fmitoshassnv_nolowmito_logit_byposmut[order(`Pr(>|z|)`)]

fwrite(inherited_noctrl_fmitoshassnv_nolowmito_logit_byposmut, file = "Report/SNVs/diff/inherited_noctrl_fmitoshassnv_nolowmito_logit_byposmut.csv")
inherited_noctrl_fmitoshassnv_nolowmito_logit_byposmut <- fread(file = "Report/SNVs/diff/inherited_noctrl_fmitoshassnv_nolowmito_logit_byposmut.csv")

###########################################################################
## make a figure pooling all mice, using the stats from fmitoshassnv logistic 
## regression, excluding the low-mito cells
###########################################################################
inherited_noctrl_fmitoshassnv_bycell_byposmut <- fread("Report/SNVs/origin/highdepth_highaf_inherited_noctrl_fmitoshassnv_bycell_byposmut.csv")
inherited_noctrl_nmitoshasdata_bycell_byposmut <- fread(file = "Report/SNVs/origin/highdepth_inherited_noctrl_nmitoshasdata_bycell_byposmut.csv")
inherited_noctrl_fmitoshassnv_bycell_byposmut <- inherited_noctrl_fmitoshassnv_bycell_byposmut[CellInfo[, "CellUID"], on = "CellUID"]
inherited_noctrl_nmitoshasdata_bycell_byposmut <- inherited_noctrl_nmitoshasdata_bycell_byposmut[CellInfo[, "CellUID"], on = "CellUID"]

inherited_noctrl_fmitoshassnv_nolowmito_logit_byposmut <- fread(file = "Report/SNVs/diff/inherited_noctrl_fmitoshassnv_nolowmito_logit_byposmut.csv")
nmitoshasdata_th <- 3
pdf(file = "Report/SNVs/diff/inherited_noctrl_fmitoshassnv_nolowmito_logit_byposmut_micepooled_violin.pdf", width = 4, height = 6)
for (posmut in inherited_noctrl_fmitoshassnv_nolowmito_logit_byposmut[!is.na(`Pr(>|z|)`)][, posmut]) {
    message(posmut)
    pval <- inherited_noctrl_fmitoshassnv_nolowmito_logit_byposmut[posmut, `Pr(>|z|)`, on = "posmut"]
    padj <- inherited_noctrl_fmitoshassnv_nolowmito_logit_byposmut[posmut, `padj`, on = "posmut"]
    X <- copy(inherited_noctrl_fmitoshassnv_bycell_byposmut[MouseID != "Mouse16&17"])
    X[, MouseID := factor(MouseID, levels = sort(unique(MouseID)))]
    Y <- inherited_noctrl_nmitoshasdata_bycell_byposmut[MouseID != "Mouse16&17"]
    y <- Y[, get(posmut)]
    idx <- y >= nmitoshasdata_th
    print(ggplot(X[idx], aes(x = CellType, y = 100 * get(posmut))) + theme_classic(base_size = 16) + geom_violin() + geom_jitter(alpha = 1, aes(col = MouseID)) + xlab("") + ylab("% mitos with SNV per cell") + ggtitle(sprintf("%s, pval = %.3g, padj = %.3g", posmut, pval, padj)) + scale_x_discrete(labels = c("Astrocyte" = "A", "Neuron" = "N")))
}
dev.off()

###########################################################################
## Finally, we went to compare the variance between astrocytes and neurons.
###########################################################################
inherited_noctrl_altasin_bymito_byposmut <- fread(file = "Report/SNVs/hierarchy/inherited_noctrl_altasin_bymito_byposmut.csv")
inherited_noctrl_altasin_astro_nestedaov_byposmut <- sapply(inherited_posmuts, function(posmut) {
    X <- inherited_noctrl_altasin_bymito_byposmut[MouseID != "Mouse16&17" & CellType == "Astrocyte", c(posmut, "MouseID", "CellUID", "CellType", "MitoID"), with = FALSE]
    names(X)[1] <- "VAF"
    mod <- aov(VAF ~ MouseID + Error(CellUID), data = X)
    X <- rbind(as.data.frame(summary(mod)[[1]][[1]]), as.data.frame(summary(mod)[[2]][[1]]))
    rownames(X) <- c("mouse", "cell", "mito")
    X
}, simplify = FALSE)
inherited_noctrl_altasin_neuron_nestedaov_byposmut <- sapply(inherited_posmuts, function(posmut) {
    X <- inherited_noctrl_altasin_bymito_byposmut[MouseID != "Mouse16&17" & CellType == "Neuron", c(posmut, "MouseID", "CellUID", "CellType", "MitoID"), with = FALSE]
    names(X)[1] <- "VAF"
    mod <- aov(VAF ~ MouseID + Error(CellUID), data = X)
    X <- rbind(as.data.frame(summary(mod)[[1]][[1]]), as.data.frame(summary(mod)[[2]][[1]]))
    rownames(X) <- c("mouse", "cell", "mito")
    X
}, simplify = FALSE)

inherited_noctrl_altasin_vartest_byposmut <- t(mapply(
    function(X, Y) {
        dfx <- X["mito", "Df"]; dfy <- Y["mito", "Df"]
        msx <- X["mito", "Mean Sq"]; msy <- Y["mito", "Mean Sq"]
        fx <- msx / msy
        pval <- ifelse(
            fx > 1, 
            pf(fx, df1 = dfx, df2 = dfy, lower.tail = FALSE) + pf(1/fx, df1 = dfx, df2 = dfy, lower.tail = TRUE), 
            ifelse(
                fx < 1, 
                pf(fx, df1 = dfx, df2 = dfy, lower.tail = TRUE) + pf(1/fx, df1 = dfx, df2 = dfy, lower.tail = FALSE), 
                NA
            )
        )
        c(df_astro = dfx, 
          df_neuron = dfy,
          ms_astro = msx, 
          ms_neuron = msy, 
          fstat_astro_to_neuron = fx, 
          fstat_neuron_to_astro = 1/fx, 
          pval = pval
        )
    }, X = inherited_noctrl_altasin_astro_nestedaov_byposmut, Y = inherited_noctrl_altasin_neuron_nestedaov_byposmut
))
inherited_noctrl_altasin_vartest_byposmut <- data.table(posmut = inherited_posmuts, inherited_noctrl_altasin_vartest_byposmut)
inherited_noctrl_altasin_vartest_byposmut[, padj := p.adjust(pval)]
inherited_noctrl_altasin_vartest_byposmut <- inherited_noctrl_altasin_vartest_byposmut[order(pval)]

fwrite(inherited_noctrl_altasin_vartest_byposmut, file = "Report/SNVs/diff/inherited_noctrl_altasin_vartest_byposmut.csv")
inherited_noctrl_altasin_vartest_byposmut <- fread(file = "Report/SNVs/diff/inherited_noctrl_altasin_vartest_byposmut.csv")

pdf(file = "Report/SNVs/diff/inherited_noctrl_altasin_vartest_nolowvar_scatter.pdf", width = 8, height = 6, useDingbats = FALSE)
ggplot(inherited_noctrl_altasin_vartest_byposmut[(ms_astro > 0.01 | ms_neuron > 0.01) & padj < 0.05], aes(x = (ms_astro), y = (ms_neuron))) + geom_point(aes(color = -log10(padj), size = log2(pmax(fstat_astro_to_neuron, fstat_neuron_to_astro)))) + geom_text_repel(aes(x = (ms_astro), y = (ms_neuron), label = posmut), size = 4, max.overlaps = 10, data = inherited_noctrl_altasin_vartest_byposmut[(ms_astro > 0.01 | ms_neuron > 0.01) & padj < 0.05]) + geom_abline(slope = 1, intercept = 0, color = "red", linetype = 2) + theme_classic(base_size = 16) + xlab("between-mitos variance  in astrocytes") + ylab("between-mitos variance in neurons") + scale_x_continuous(limits = c(0, 0.13)) + scale_y_continuous(limits = c(0, 0.13)) + guides(size = guide_legend("var ratio (log2)")) + scale_color_viridis_c()
dev.off()
