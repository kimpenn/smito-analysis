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

    .dt_to_df <- function(dt) {
        k <- key(dt)
        if (is.null(k)) {
            warning("No key found! Using the first column as the key...")
            k <- names(dt)[1]
        } 
        df <- as.data.frame(dt)
        if (any(duplicated(dt[, get(k)]))) {
            df
        } else {
            rownames(df) <- dt[, get(k)]
            df[, -match(k, names(dt)), drop = FALSE]
        }
    }

    mutate_from_transmut <- function(pos, ref, size, transmut) {
        ref <- toupper(ref)
        transmut <- as.matrix(.dt_to_df(transmut))
        transmut[is.na(transmut)] <- 0
        n <- sum(transmut)
        transmut <- transmut / n
        margin <- rowSums(transmut)
        if (missing(size)) {
            size <- n
        }
        ns <- round(size * margin)
        refbases <- c('A', 'C', 'G', 'T')
        ms <- structure(as.vector(table(factor(ref, levels = refbases))), names = refbases)
        nm <- ms / ns
        i <- which.min(nm)

        if (nm[i] < 1) {
            ns <- round(ns * nm[i])
            warning(sprintf("Some base (%s) has fewer reference sites than can be sampled, shrinking `size` to %d...", refbases[i], sum(ns)))
        }
        if (any(ns > ms)) {
            stop("Cannot afford the sample size even after shrinking!")
        }

        res <- lapply(refbases, function(refbase) {
            p <- pos[ref == refbase]
            i <- which(refbases == refbase)
            K <- ns[i]
            r <- unlist(transmut[refbase, ])
            r <- r / sum(r)
            m <- round(r * K)
            f <- rep(factor(seq(length(r))), times = m)
            K <- sum(m)
            q <- sample(p, size = K, replace = FALSE)
            t <- split(q, f = f)
            s <- paste0(refbase, '>', names(r))
            R <- mapply(FUN = function(x, y) {
                if (length(x) > 0) {
                    data.table::data.table(pos = x, mut = y)
                } else {
                    data.table::data.table(pos = integer(0), mut = character(0))
                }
            }, x = t, y = s, SIMPLIFY = FALSE)
            rbindlist(R)
        })
        res <- rbindlist(res)
        res[, ref := sapply(strsplit(mut, '>'), '[', 1)]
        res[, alt := sapply(strsplit(mut, '>'), '[', 2)]
        res[, c("pos", "ref", "alt", "mut")]
        res[order(pos)]
    }

    KaKs <- function(total_pos, var_pos, codon_table) {
        nonsyn_sites <- vertmtcodon[pos %in% var_pos, list(syn = mean(S), nonsyn = mean(N)), by = c("symbol", "strand", "codon_index")][, sum(nonsyn)]
        nonsyn_total <- vertmtcodon[pos %in% total_pos, list(syn = mean(S), nonsyn = mean(N)), by = c("symbol", "strand", "codon_index")][, sum(nonsyn)]
        syn_sites <- vertmtcodon[pos %in% var_pos, list(syn = mean(S), nonsyn = mean(N)), by = c("symbol", "strand", "codon_index")][, sum(syn)]
        syn_total <- vertmtcodon[pos %in% total_pos, list(syn = mean(S), nonsyn = mean(N)), by = c("symbol", "strand", "codon_index")][, sum(syn)]
        (nonsyn_sites / nonsyn_total) / (syn_sites / syn_total)
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
