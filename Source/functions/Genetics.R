###########################################################################
## Author: Youtao Lu <luyoutao@sas.upenn.edu>
##  
## Copyright (c) 2021-2023, Youtao Lu and Junhyong Kim, Department of Biology, University of Pennsylvania
## Copyright (c) 2021-2023, James Eberwine, Perelman School of Medicine, University of Pennsylvania
## Copyright (c) 2021-2023, David Issadore, Penn Engineering, University of Pennsylvania
## 
## This Source Code is distributed under Creative Commons Attribution License 4.0 (CC BY).
###########################################################################

if (!exists("Genetics") || is.environment(Genetics)) Genetics <- new.env(parent = emptyenv()) 
local({
    .VERSION = "0.1"

    titv <- function(transmat) {
        ti <- sum(c(transmat["A", "G"], transmat["G", "A"], transmat["C", "T"], transmat["T", "C"]), na.rm = TRUE)
        tv <- sum(c(transmat["A", "C"], transmat["A", "T"], transmat["C", "A"], transmat["C", "G"], 
                    transmat["G", "C"], transmat["G", "T"], transmat["T", "A"], transmat["T", "G"]), na.rm = TRUE)
        c(ti = ti, tv = tv, r = ti/tv)
    }

    link_bin_contab <- function(df, lab = "") {
        pos <- colnames(df)
        sapply(pos, function(j) {
            sapply(pos, function(i) {
                Y <- df[, c(i, j)]
                y1 <- Y[, 1]
                y2 <- Y[, 2]
                idx <- y1 != "" & y2 != ""
                Z <- Y[idx, ]
                y1 <- y1[idx]
                y2 <- y2[idx]
                sapply(bases, function(b) { 
                    sapply(bases, function(a) {
                        message(lab, " ", i, " ", j, " ", a, " ", b)
                        a1 <- grepl(y1, pattern = a)
                        b2 <- grepl(y2, pattern = b)
                        x11 <- sum(a1 & b2)
                        x21 <- sum((!a1) & b2)
                        x12 <- sum(a1 & !b2)
                        x22 <- sum((!a1) & (!b2))
                        contab <- matrix(c(x11, x21, x12, x22), ncol = 2, nrow = 2, dimnames = list(c(paste0("allele1==", a), paste0("allele1!=", a)), c(paste0("allele2==", b), paste0("allele2!=", b))))
                    }, simplify = FALSE)
                }, simplify = FALSE)
            }, simplify = FALSE)
        }, simplify = FALSE)
    }

    link_bin_fet <- function(contab, lab = "") {
        pos <- names(contab)
        sapply(pos, function(j) {
            sapply(pos, function(i) {
                sapply(bases, function(b) { 
                    sapply(bases, function(a) {
                        message(lab, " ", i, " ", j, " ", a, " ", b)
                        Y <- contab[[j]][[i]][[b]][[a]]
                        c(oddsratio = unname(fisher.test(Y, alternative = "greater")[["estimate"]]), 
                          pval      = unname(fisher.test(Y, alternative = "greater")[["p.value"]]))
                    }, simplify = FALSE)
                }, simplify = FALSE)
            }, simplify = FALSE)
        }, simplify = FALSE)
    }

    for (obj in ls()) 
        assign(obj, get(obj), envir = Genetics)
})
