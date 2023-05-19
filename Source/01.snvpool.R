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

mito_barcodes <- read.csv("Data/mito_barcodes.csv", as.is = TRUE)
mitoIDs <- mito_barcodes[, "ID"] 

snv_info <- read.csv("Data/snv_loci_v2.csv", as.is = TRUE)
snvIDs <- snv_info[,"ID"]

chrmbases <- readLines("Data/mm10.mito/chrM.fa")[-1]
chrmbases <- paste0(chrmbases, collapse = "")
chrmbases <- strsplit(chrmbases, "")[[1]]
nchrmbases <- length(chrmbases)

MitoInfo <- fread(file = "Report/metadata/MitoInfo.csv")

basedifffreq_cutdemux_rerun_q30_unstranded <- read.csv("Report/SNVs/basedifffreq_cutdemux_rerun_sub500k_q30_unstranded.csv.gz", as.is = TRUE, check.names = FALSE)
dim(basedifffreq_cutdemux_rerun_q30_unstranded)
## [1] 1738811      13
basedifffreq_cutdemux_E749_q30_unstranded <- read.csv("Report/SNVs/basedifffreq_cutdemux_E749_sub500k_q30_unstranded.csv.gz", as.is = TRUE, check.names = FALSE)
dim(basedifffreq_cutdemux_E749_q30_unstranded)
## [1] 436025     13
basedifffreq_cutdemux_E754_q30_unstranded <- read.csv("Report/SNVs/basedifffreq_cutdemux_E754_sub500k_q30_unstranded.csv.gz", as.is = TRUE, check.names = FALSE)
dim(basedifffreq_cutdemux_E754_q30_unstranded)
## [1] 352826     13
basedifffreq_cutdemux_E757_q30_unstranded <- read.csv("Report/SNVs/basedifffreq_cutdemux_E757_sub500k_q30_unstranded.csv.gz", as.is = TRUE, check.names = FALSE)
dim(basedifffreq_cutdemux_E757_q30_unstranded)
## [1] 160858     13
basedifffreq_cutdemux_E771_q30_unstranded <- read.csv("Report/SNVs/basedifffreq_cutdemux_E771_sub500k_q30_unstranded.csv.gz", as.is = TRUE, check.names = FALSE)
dim(basedifffreq_cutdemux_E771_q30_unstranded)
## [1] 148146     13
basedifffreq_cutdemux_E777_q30_unstranded <- read.csv("Report/SNVs/basedifffreq_cutdemux_E777_sub500k_q30_unstranded.csv.gz", as.is = TRUE, check.names = FALSE)
dim(basedifffreq_cutdemux_E777_q30_unstranded)
## [1] 313974     13
basedifffreq_cutdemux_E799_q30_unstranded <- read.csv("Report/SNVs/basedifffreq_cutdemux_E799_sub500k_q30_unstranded.csv.gz", as.is = TRUE, check.names = FALSE)
dim(basedifffreq_cutdemux_E799_q30_unstranded)
## [1] 67242    13

basedifffreq_cutdemux_q30_unstranded <- rbind(
    basedifffreq_cutdemux_rerun_q30_unstranded, 
    basedifffreq_cutdemux_E749_q30_unstranded, 
    basedifffreq_cutdemux_E754_q30_unstranded,
    basedifffreq_cutdemux_E757_q30_unstranded,
    basedifffreq_cutdemux_E771_q30_unstranded,
    basedifffreq_cutdemux_E777_q30_unstranded,
    basedifffreq_cutdemux_E799_q30_unstranded
)
setDT(basedifffreq_cutdemux_q30_unstranded)
basedifffreq_cutdemux_q30_unstranded[, LibraryMitoID := paste0(LibraryID, "_", MitoID)]
dim(basedifffreq_cutdemux_q30_unstranded)
## [1] 3217882      14
setcolorder(basedifffreq_cutdemux_q30_unstranded, "LibraryMitoID")
setnames(basedifffreq_cutdemux_q30_unstranded, "D", "del")

## remove the old 34 mitos as well as those don't exist in the metadata
basedifffreq_cutdemux_q30_unstranded <- basedifffreq_cutdemux_q30_unstranded[LibraryMitoID %in% MitoInfo[, LibraryMitoID]]
dim(basedifffreq_cutdemux_q30_unstranded)
## [1] 3117733      14

basedifffreq_cutdemux_q30_unstranded <- MitoInfo[basedifffreq_cutdemux_q30_unstranded, on = "LibraryMitoID", nomatch = NULL]
basedifffreq_cutdemux_q30_unstranded[, c("i.ExptID", "i.LibraryID", "i.MitoID") := NULL]
dim(basedifffreq_cutdemux_q30_unstranded)
## [1] 3117733      29

fwrite(basedifffreq_cutdemux_q30_unstranded, file = "Report/SNVs/basedifffreq_cutdemux_sub500k_q30_unstranded.csv.gz")
basedifffreq_cutdemux_q30_unstranded <- fread("Report/SNVs/basedifffreq_cutdemux_sub500k_q30_unstranded.csv.gz")

###########################################################################
## gather all insertion tables
###########################################################################
MitoInfo <- fread(file = "Report/metadata/MitoInfo.csv")
dim(MitoInfo)
## [1] 1717   19

mpileups_ins_cutdemux_rerun <- read.csv("Report/SNVs/ins/mpileups_ins_cutdemux_rerun_sub500k.csv.gz", as.is = TRUE)
mpileups_ins_cutdemux_E749 <- read.csv("Report/SNVs/ins/mpileups_ins_cutdemux_E749_sub500k.csv.gz", as.is = TRUE)
mpileups_ins_cutdemux_E754 <- read.csv("Report/SNVs/ins/mpileups_ins_cutdemux_E754_sub500k.csv.gz", as.is = TRUE)
mpileups_ins_cutdemux_E757 <- read.csv("Report/SNVs/ins/mpileups_ins_cutdemux_E757_sub500k.csv.gz", as.is = TRUE)
mpileups_ins_cutdemux_E771 <- read.csv("Report/SNVs/ins/mpileups_ins_cutdemux_E771_sub500k.csv.gz", as.is = TRUE)
mpileups_ins_cutdemux_E777 <- read.csv("Report/SNVs/ins/mpileups_ins_cutdemux_E777_sub500k.csv.gz", as.is = TRUE)
mpileups_ins_cutdemux_E799 <- read.csv("Report/SNVs/ins/mpileups_ins_cutdemux_E799_sub500k.csv.gz", as.is = TRUE)
mpileups_ins_cutdemux <- rbind(mpileups_ins_cutdemux_rerun, 
                               mpileups_ins_cutdemux_E749, 
                               mpileups_ins_cutdemux_E754, 
                               mpileups_ins_cutdemux_E757,
                               mpileups_ins_cutdemux_E771,
                               mpileups_ins_cutdemux_E777,
                               mpileups_ins_cutdemux_E799
)
setDT(mpileups_ins_cutdemux)
dim(mpileups_ins_cutdemux)
## [1] 84528    10
## [1] 86094    10

## remove the old 34 mitos as well as those don't exist in the metadata
## remove technical replicates
## remove restriction enzyme digested controls
mpileups_ins_cutdemux <- mpileups_ins_cutdemux[paste0(LibraryID, "_", MitoID) %in% MitoInfo[, LibraryMitoID]]
dim(mpileups_ins_cutdemux)
## [1] 84603    10

mpileups_ins_cutdemux[, LibraryMitoID := paste0(LibraryID, "_", MitoID)]
mpileups_ins_cutdemux[, altperctot := altfreqtot / depth * 100]
setdiff(mpileups_ins_cutdemux[, LibraryMitoID], MitoInfo[, LibraryMitoID])
## character(0)

mpileups_ins_cutdemux <- MitoInfo[mpileups_ins_cutdemux, on = "LibraryMitoID", nomatch = NULL]
mpileups_ins_cutdemux[, c("i.ExptID", "i.LibraryID", "i.MitoID") := NULL]
dim(mpileups_ins_cutdemux)
## [1] 83945    27
## [1] 85490    27
## [1] 84603    27

chrmproperties <- fread(file = "Report/artifact/chrmbases_properties.csv.gz")
setkey(chrmproperties, "pos")

mpileups_ins_cutdemux <- chrmproperties[mpileups_ins_cutdemux, on = "pos"]
mpileups_ins_cutdemux[, i.ref := NULL]
setcolorder(mpileups_ins_cutdemux, colnames(mpileups_ins_cutdemux)[c(10:29, 1:9, 30:34)])

fwrite(mpileups_ins_cutdemux, file = "Report/SNVs/ins/mpileups_ins_cutdemux_sub500k.csv.gz")
mpileups_ins_cutdemux <- fread(file = "Report/SNVs/ins/mpileups_ins_cutdemux_sub500k.csv.gz")
