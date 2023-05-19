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

dir.create("Report/SNVs/filter", FALSE, TRUE)

snv_info <- fread("Data/snv_loci_v2.csv")
snvIDs <- snv_info[, SNVID]

chrmproperties <- fread("Report/artifact/chrmbases_properties.csv.gz")
setkey(chrmproperties, "pos")

MitoInfo <- fread(file = "Report/metadata/MitoInfo.csv")
dim(MitoInfo)
## [1] 1717   19
CellInfo <- fread("Report/metadata/CellInfo.csv")
dim(CellInfo)
## [1] 102  12
MouseInfo <- fread("Report/metadata/MouseInfo.csv")
dim(MouseInfo)
## [1] 13  2

###########################################################################
## load data, the original version that had missing data at 9027 because of
## the 150/151 mixed cycles
###########################################################################
basedifffreq_cutdemux_q30_unstranded <- fread("Report/SNVs/basedifffreq_cutdemux_sub500k_q30_unstranded.csv.gz")
## make sure the lower-case is converted to upper-case
basedifffreq_cutdemux_q30_unstranded[, ref := toupper(ref)]
basedifffreq_cutdemux_q30_unstranded[pos == 9027, list(LibraryMitoID, `=`, `A`, `C`, `G`, `T`)]
##       LibraryMitoID  =   A C G T
##    1:   L1R24P1C_M1  1  28 0 0 0
##    2:   L1R24P1C_M2 37  44 0 0 0
##    3:   L1R24P1C_M3 53  46 0 0 0
##    4:   L1R24P1C_M4 12 167 0 0 0
##    5:   L1R24P1C_M5  0   1 0 0 0
##   ---
## 1655:  DepSeq-L5_M2 34 157 0 0 0
## 1656:  DepSeq-L5_M3 11 142 1 0 0
## 1657:  DepSeq-L5_M4 11 286 0 0 0
## 1658:  DepSeq-L5_M6 32   0 0 0 0
## 1659:  DepSeq-L5_M7  1  38 0 0 0

## load 9027 matched and clipped counts
basedifffreq_end_q30_bysnv <- fread(file = "Report/artifact/basedifffreq_end_q30_bysnv.csv.gz")
setkey(basedifffreq_end_q30_bysnv, LibraryMitoID)
basedifffreq_9027_q30 <- basedifffreq_end_q30_bysnv[pos == 9027 & IsClip == "Y"][J(basedifffreq_cutdemux_q30_unstranded[pos == 9027, LibraryMitoID]), list(LibraryMitoID, `A`, `C`, `G`, `T`)]
basedifffreq_9027_q30[is.na(`A`), `A` := 0]
basedifffreq_9027_q30[is.na(`C`), `C` := 0]
basedifffreq_9027_q30[is.na(`G`), `G` := 0]
basedifffreq_9027_q30[is.na(`T`), `T` := 0]
basedifffreq_9027_q30
##       LibraryMitoID      A  C G  T
##    1:   L1R24P1C_M1  12255  0 0  0
##    2:   L1R24P1C_M2 255117  7 0 28
##    3:   L1R24P1C_M3 290043  4 0 31
##    4:   L1R24P1C_M4 127776  2 0 12
##    5:   L1R24P1C_M5    110  0 0  0
##   ---
## 1655:  DepSeq-L5_M2   2783  1 0  7
## 1656:  DepSeq-L5_M3 199386 23 0 37
## 1657:  DepSeq-L5_M4  59197  9 0 11
## 1658:  DepSeq-L5_M6      5  0 0  0
## 1659:  DepSeq-L5_M7     31  0 0  0

## compensate the softclipped bases to 9027
basedifffreq_cutdemux_q30_unstranded[pos == 9027, `A` := `A` + basedifffreq_9027_q30[, `A`]]
basedifffreq_cutdemux_q30_unstranded[pos == 9027, `C` := `C` + basedifffreq_9027_q30[, `C`]]
basedifffreq_cutdemux_q30_unstranded[pos == 9027, `=` := `=` + basedifffreq_9027_q30[, `G`]] ## must be zero G
basedifffreq_cutdemux_q30_unstranded[pos == 9027, `T` := `T` + basedifffreq_9027_q30[, `T`]]
basedifffreq_cutdemux_q30_unstranded[pos == 9027, list(LibraryMitoID, `=`, `A`, `C`, `G`, `T`)]
##       LibraryMitoID  =      A  C G  T
##    1:   L1R24P1C_M1  1  12283  0 0  0
##    2:   L1R24P1C_M2 37 255161  7 0 28
##    3:   L1R24P1C_M3 53 290089  4 0 31
##    4:   L1R24P1C_M4 12 127943  2 0 12
##    5:   L1R24P1C_M5  0    111  0 0  0
##   ---
## 1655:  DepSeq-L5_M2 34   2940  1 0  7
## 1656:  DepSeq-L5_M3 11 199528 24 0 37
## 1657:  DepSeq-L5_M4 11  59483  9 0 11
## 1658:  DepSeq-L5_M6 32      5  0 0  0
## 1659:  DepSeq-L5_M7  1     69  0 0  0
basedifffreq_cutdemux_q30_unstranded[pos == 9027, `depth` := `=` + `A` + `C` + `G` + `T` + `del`]
basedifffreq_cutdemux_q30_unstranded[pos == 9027, list(LibraryMitoID, depth, `=`, `A`, `C`, `G`, `T`, del)]
##       LibraryMitoID  depth  =      A  C G  T del
##    1:   L1R24P1C_M1  12284  1  12283  0 0  0   0
##    2:   L1R24P1C_M2 255233 37 255161  7 0 28   0
##    3:   L1R24P1C_M3 290177 53 290089  4 0 31   0
##    4:   L1R24P1C_M4 127969 12 127943  2 0 12   0
##    5:   L1R24P1C_M5    111  0    111  0 0  0   0
##   ---
## 1655:  DepSeq-L5_M2   2982 34   2940  1 0  7   0
## 1656:  DepSeq-L5_M3 199600 11 199528 24 0 37   0
## 1657:  DepSeq-L5_M4  59514 11  59483  9 0 11   0
## 1658:  DepSeq-L5_M6     37 32      5  0 0  0   0
## 1659:  DepSeq-L5_M7     70  1     69  0 0  0   0
fwrite(basedifffreq_cutdemux_q30_unstranded, file = "Report/SNVs/basedifffreq_cutdemux_sub500k_q30_unstranded_9027fixed.csv.gz")
basedifffreq_cutdemux_q30_unstranded <- fread(file = "Report/SNVs/basedifffreq_cutdemux_sub500k_q30_unstranded_9027fixed.csv.gz")

basediffperc_cutdemux_q30_unstranded <- copy(basedifffreq_cutdemux_q30_unstranded)
basediffperc_cutdemux_q30_unstranded[, c("=", "A", "C", "G", "T", "del") := lapply(.SD, function(x) x / depth * 100),.SDcols = c("=", "A", "C", "G", "T", "del")]
dim(basediffperc_cutdemux_q30_unstranded)
## [1] 3117733      29

depth_th <- 50
basediffperc_cutdemux_q30_unstranded_highdepth <- basediffperc_cutdemux_q30_unstranded[depth >= depth_th]
dim(basediffperc_cutdemux_q30_unstranded_highdepth)
## [1] 1433039      36 when depth_th = 100
## [1] 1675309      36
## [1] 1717206      36 when deep-seq-ed 34 mitos replaced the old data
## [1] 1702508      29 when enzyme digested mitos and replicate mito are removed
## [1] 1703153      29 after patching 9027 softclipped A/C/T bases
fwrite(basediffperc_cutdemux_q30_unstranded_highdepth, file = "Report/SNVs/filter/basediffperc_cutdemux_sub500k_q30_unstranded_highdepth.csv.gz")
basediffperc_cutdemux_q30_unstranded_highdepth <- fread("Report/SNVs/filter/basediffperc_cutdemux_sub500k_q30_unstranded_highdepth.csv.gz")

###########################################################################
## SNV calls
###########################################################################
vaf_th <- 0.05
basediffperc_cutdemux_q30_unstranded_highdepth_highaf <- basediffperc_cutdemux_q30_unstranded_highdepth[`A` >= vaf_th * 100 | `C` >= vaf_th * 100 | `G` >= vaf_th * 100 | `T` >= vaf_th * 100 | `del` >= vaf_th * 100]
dim(basediffperc_cutdemux_q30_unstranded_highdepth_highaf)
## [1] 5240   36 when depth_th = 100
## [1] 7030   36
## [1] 7219   36 when deep-seq-ed 34 mitos replaced the old data
## [1] 7164   29 when enzyme digested mitos and replicate mito are removed
## [1] 7967   29 after patching 9027
fwrite(basediffperc_cutdemux_q30_unstranded_highdepth_highaf, file = "Report/SNVs/filter/basediffperc_cutdemux_sub500k_q30_unstranded_highdepth_highaf.csv.gz")
basediffperc_cutdemux_q30_unstranded_highdepth_highaf <- fread("Report/SNVs/filter/basediffperc_cutdemux_sub500k_q30_unstranded_highdepth_highaf.csv.gz")

###########################################################################
## BWA vs. STAR, take the overlapping SNV sites
## Note, is_by_bwa is conditioned on the SNVs called by BWA at 5%, while 
## is_in_range and is_in_primer are irrespective to SNV calls. 
###########################################################################
pos_noctrl_starvsbwa <- fread("Report/bwa/snv_pos_noctrl.csv")
pos_noctrl_starvsbwa <- na.omit(pos_noctrl_starvsbwa[, intersect(all_star, all_bwa)])
pos_noctrl_starvsbwa <- data.table(pos = basediffperc_cutdemux_q30_unstranded_highdepth_highaf[IsCtrl == "N", sort(unique(pos))], is_by_bwa = ifelse(basediffperc_cutdemux_q30_unstranded_highdepth_highaf[IsCtrl == "N", sort(unique(pos))] %in% pos_noctrl_starvsbwa, "Y", "N"))
pos_noctrl_starvsbwa[, table(is_by_bwa)]
## is_by_bwa
##   N   Y 
##  54 896
pos_noctrl_qc <- merge.data.table(pos_noctrl_starvsbwa, chrmproperties[, c("pos", "is_in_range", "is_in_primer")], by = "pos")
pos_noctrl_qc[is_by_bwa == "Y" & is_in_range == "Y" & is_in_primer == "N", .N]
## [1] 838
fwrite(pos_noctrl_qc, file = "Report/SNVs/filter/basediffperc_cutdemux_sub500k_q30_unstranded_highdepth_highaf_qc.csv")
pos_noctrl_qc <- fread(file = "Report/SNVs/filter/basediffperc_cutdemux_sub500k_q30_unstranded_highdepth_highaf_qc.csv")

###########################################################################
## filter for the BWA-STAR overlapping, in-range, out-of-primers SNV sites
###########################################################################
basediffperc_cutdemux_q30_unstranded_highdepth_highaf_qcfltd <- basediffperc_cutdemux_q30_unstranded_highdepth_highaf[pos %in% pos_noctrl_qc[is_by_bwa == "Y" & is_in_range == "Y" & is_in_primer == "N", pos]]
dim(basediffperc_cutdemux_q30_unstranded_highdepth_highaf_qcfltd)
## [1] 6705   29
## [1] 6354   29 after filtering out SNVs inside PCR forward primers
## [1] 7157   29 after fixing 9027

fwrite(basediffperc_cutdemux_q30_unstranded_highdepth_highaf_qcfltd, file = "Report/SNVs/filter/basediffperc_cutdemux_sub500k_q30_unstranded_highdepth_highaf_qcfltd.csv.gz")
basediffperc_cutdemux_q30_unstranded_highdepth_highaf_qcfltd <- fread("Report/SNVs/filter/basediffperc_cutdemux_sub500k_q30_unstranded_highdepth_highaf_qcfltd.csv.gz")

###########################################################################
## How many mitos/cells/mice support each SNV after QC filtering?
###########################################################################
basediffperc_cutdemux_q30_unstranded_highdepth_highaf_qcfltd_noctrl <- basediffperc_cutdemux_q30_unstranded_highdepth_highaf_qcfltd[IsCtrl == "N"]
dim(basediffperc_cutdemux_q30_unstranded_highdepth_highaf_qcfltd_noctrl)
## [1] 6314   29
## [1] 5988   29
## [1] 6763   29 after fixing 9027
basediffperc_cutdemux_q30_unstranded_highdepth_highaf_qcfltd_noctrl[ref == "A", A := `=`]
basediffperc_cutdemux_q30_unstranded_highdepth_highaf_qcfltd_noctrl[ref == "C", C := `=`]
basediffperc_cutdemux_q30_unstranded_highdepth_highaf_qcfltd_noctrl[ref == "G", G := `=`]
basediffperc_cutdemux_q30_unstranded_highdepth_highaf_qcfltd_noctrl[ref == "T", T := `=`]

basediffperc_cutdemux_q30_unstranded_highdepth_highaf_qcfltd_noctrl <- melt.data.table(basediffperc_cutdemux_q30_unstranded_highdepth_highaf_qcfltd_noctrl, measure.vars = c("A", "C", "G", "T", "del"), variable.name = "alt", value.name = "altperc")
setnames(basediffperc_cutdemux_q30_unstranded_highdepth_highaf_qcfltd_noctrl, "=", "refperc")
dim(basediffperc_cutdemux_q30_unstranded_highdepth_highaf_qcfltd_noctrl)
## [1] 31570    26
## [1] 29940    26
## [1] 33815    26 after fixing 9027

## keep only variant alleles
basediffperc_cutdemux_q30_unstranded_highdepth_highaf_qcfltd_noctrl <- basediffperc_cutdemux_q30_unstranded_highdepth_highaf_qcfltd_noctrl[ref != alt]
dim(basediffperc_cutdemux_q30_unstranded_highdepth_highaf_qcfltd_noctrl)
## [1] 25256    26 after digested and replicate mitos are removed
## [1] 23952    26
## [1] 27052    26

basediffperc_cutdemux_q30_unstranded_highdepth_highaf_qcfltd_noctrl[, summary(altperc)]
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
##   0.000   0.000   0.000   8.111   5.073 100.000
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
##   0.000   0.000   0.000   8.436   5.085 100.000
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.   after fixing 9027
##   0.000   0.000   0.000  10.857   5.085 100.000

## keep only variant alleles whose AF are >= 5% 
vaf_th <- 0.05
basediffperc_cutdemux_q30_unstranded_highdepth_highaf_qcfltd_noctrl <- basediffperc_cutdemux_q30_unstranded_highdepth_highaf_qcfltd_noctrl[altperc >= vaf_th * 100]
dim(basediffperc_cutdemux_q30_unstranded_highdepth_highaf_qcfltd_noctrl)
## [1] 6411   26 after digested and replicate mitos are removed
## [1] 6085   26
## [1] 6863   26 after fixing 9027

basediffperc_cutdemux_q30_unstranded_highdepth_highaf_qcfltd_noctrl[, alt := factor(alt, levels = c("A", "C", "G", "T", "del"))]
basediffperc_cutdemux_q30_unstranded_highdepth_highaf_qcfltd_noctrl[, mut := factor(paste0(ref, ">", alt), levels = paste0(rep(c("A", "C", "G", "T"), each = 5), ">", c("A", "C", "G", "T", "del")))]
setkey(basediffperc_cutdemux_q30_unstranded_highdepth_highaf_qcfltd_noctrl, pos, alt)
basediffperc_cutdemux_q30_unstranded_highdepth_highaf_qcfltd_noctrl[, posmut := paste0(pos, ":", mut)]
basediffperc_cutdemux_q30_unstranded_highdepth_highaf_qcfltd_noctrl[, posmut := factor(posmut, levels = unique(posmut))]
dim(basediffperc_cutdemux_q30_unstranded_highdepth_highaf_qcfltd_noctrl)
## [1] 6085   28
## [1] 6863   26 after fixing 9027
basediffperc_cutdemux_q30_unstranded_highdepth_highaf_qcfltd_noctrl[, summary(altperc)]
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
##   5.000   6.835  11.250  31.464  45.535 100.000
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
##   5.000   6.897  11.765  32.696  50.704 100.000
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. after fixing 9027
##   5.000   7.317  15.855  42.322  96.926 100.000

## number of sample support by variant alleles
basediffperc_cutdemux_q30_unstranded_highdepth_highaf_qcfltd_support <- basediffperc_cutdemux_q30_unstranded_highdepth_highaf_qcfltd_noctrl[, list(SNVID = unique(SNVID), nmice = uniqueN(MouseID), ncells = uniqueN(CellUID), nmitos = uniqueN(LibraryMitoID)), by = "posmut"]
dim(basediffperc_cutdemux_q30_unstranded_highdepth_highaf_qcfltd_support)
## [1] 1056    5
## [1] 1030    5
## [1] 1032    5 after fixing 9027
basediffperc_cutdemux_q30_unstranded_highdepth_highaf_qcfltd_support[, pos := as.integer(sapply(strsplit(as.character(posmut), ":"), '[', 1))]
basediffperc_cutdemux_q30_unstranded_highdepth_highaf_qcfltd_support[, mut := sapply(strsplit(as.character(posmut), ":"), '[', 2)]
basediffperc_cutdemux_q30_unstranded_highdepth_highaf_qcfltd_support[, ref := sapply(strsplit(as.character(mut), ">"), '[', 1)]
basediffperc_cutdemux_q30_unstranded_highdepth_highaf_qcfltd_support[, alt := sapply(strsplit(as.character(mut), ">"), '[', 2)]
setcolorder(basediffperc_cutdemux_q30_unstranded_highdepth_highaf_qcfltd_support, c("posmut", "SNVID", "pos", "mut", "ref", "alt", "nmice", "ncells", "nmitos"))
fwrite(basediffperc_cutdemux_q30_unstranded_highdepth_highaf_qcfltd_support, file = "Report/SNVs/filter/basediffperc_cutdemux_sub500k_q30_unstranded_highdepth_highaf_qcfltd_support_byposmut.csv")
basediffperc_cutdemux_q30_unstranded_highdepth_highaf_qcfltd_support_byposmut <- fread(file = "Report/SNVs/filter/basediffperc_cutdemux_sub500k_q30_unstranded_highdepth_highaf_qcfltd_support_byposmut.csv")

## number of sample support by SNV site
basediffperc_cutdemux_q30_unstranded_highdepth_highaf_qcfltd_support <- basediffperc_cutdemux_q30_unstranded_highdepth_highaf_qcfltd[IsCtrl == "N", list(SNVID = unique(SNVID), nmice = uniqueN(MouseID), ncells = uniqueN(CellUID), nmitos = uniqueN(LibraryMitoID)), by = "pos"]
dim(basediffperc_cutdemux_q30_unstranded_highdepth_highaf_qcfltd_support)
## [1] 862   5 after digested and replicate mitos are removed
## [1] 838   5
## [1] 838   5 after fixing 9027
fwrite(basediffperc_cutdemux_q30_unstranded_highdepth_highaf_qcfltd_support, file = "Report/SNVs/filter/basediffperc_cutdemux_sub500k_q30_unstranded_highdepth_highaf_qcfltd_support_bypos.csv")
basediffperc_cutdemux_q30_unstranded_highdepth_highaf_qcfltd_support_bypos <- fread(file = "Report/SNVs/filter/basediffperc_cutdemux_sub500k_q30_unstranded_highdepth_highaf_qcfltd_support_bypos.csv")

###########################################################################
## don't filter for nmitos/ncells/nmice >= k SNVs, first output all of them
## make a data table AF per SNV
###########################################################################
basediffperc_cutdemux_q30_unstranded_highdepth <- fread("Report/SNVs/filter/basediffperc_cutdemux_sub500k_q30_unstranded_highdepth.csv.gz")
basediffperc_cutdemux_q30_unstranded_highdepth_highaf_qcfltd_support_byposmut <- fread(file = "Report/SNVs/filter/basediffperc_cutdemux_sub500k_q30_unstranded_highdepth_highaf_qcfltd_support_byposmut.csv")
basediffperc_cutdemux_q30_unstranded_highdepth_qcfltd <- basediffperc_cutdemux_q30_unstranded_highdepth[pos %in% basediffperc_cutdemux_q30_unstranded_highdepth_highaf_qcfltd_support_byposmut[, pos]]
dim(basediffperc_cutdemux_q30_unstranded_highdepth_qcfltd)
## [1] 832872     29 ## after fixing 9027

basediffperc_cutdemux_q30_unstranded_highdepth_qcfltd[ref == "A", A := `=`]
basediffperc_cutdemux_q30_unstranded_highdepth_qcfltd[ref == "C", C := `=`]
basediffperc_cutdemux_q30_unstranded_highdepth_qcfltd[ref == "G", G := `=`]
basediffperc_cutdemux_q30_unstranded_highdepth_qcfltd[ref == "T", T := `=`]

basediffperc_cutdemux_q30_unstranded_highdepth_qcfltd <- melt.data.table(basediffperc_cutdemux_q30_unstranded_highdepth_qcfltd, measure.vars = c("A", "C", "G", "T", "del"), variable.name = "alt", value.name = "altperc")
setnames(basediffperc_cutdemux_q30_unstranded_highdepth_qcfltd, "=", "refperc")
dim(basediffperc_cutdemux_q30_unstranded_highdepth_qcfltd)
## [1] 4164360      26

basediffperc_cutdemux_q30_unstranded_highdepth_qcfltd <- basediffperc_cutdemux_q30_unstranded_highdepth_qcfltd[ref != alt]
dim(basediffperc_cutdemux_q30_unstranded_highdepth_qcfltd)
## [1] 3331488      26

basediffperc_cutdemux_q30_unstranded_highdepth_qcfltd[, alt := factor(alt, levels = c("A", "C", "G", "T", "del"))]
basediffperc_cutdemux_q30_unstranded_highdepth_qcfltd[, mut := factor(paste0(ref, ">", alt), levels = paste0(rep(c("A", "C", "G", "T"), each = 5), ">", c("A", "C", "G", "T", "del")))]
setkey(basediffperc_cutdemux_q30_unstranded_highdepth_qcfltd, pos, alt)
basediffperc_cutdemux_q30_unstranded_highdepth_qcfltd[, posmut := paste0(pos, ":", mut)]
basediffperc_cutdemux_q30_unstranded_highdepth_qcfltd[, posmut := factor(posmut, levels = unique(posmut))]
dim(basediffperc_cutdemux_q30_unstranded_highdepth_qcfltd)
## [1] 3331488      28

basediffperc_cutdemux_q30_unstranded_highdepth_qcfltd <- basediffperc_cutdemux_q30_unstranded_highdepth_qcfltd[as.character(posmut) %in% as.character(basediffperc_cutdemux_q30_unstranded_highdepth_highaf_qcfltd_support_byposmut[, posmut])]
dim(basediffperc_cutdemux_q30_unstranded_highdepth_qcfltd)
## [1] 1021830      28

basediffperc_cutdemux_q30_unstranded_highdepth_qcfltd <- basediffperc_cutdemux_q30_unstranded_highdepth_qcfltd[, -(27:28)]
fwrite(basediffperc_cutdemux_q30_unstranded_highdepth_qcfltd, file = "Report/SNVs/filter/basediffperc_cutdemux_sub500k_q30_unstranded_highdepth_qcfltd.csv.gz")
basediffperc_cutdemux_q30_unstranded_highdepth_qcfltd <- fread(file = "Report/SNVs/filter/basediffperc_cutdemux_sub500k_q30_unstranded_highdepth_qcfltd.csv.gz")

###########################################################################
## for all QC filtered SNVs, get the VAF table: mito X SNV
###########################################################################
basediffperc_cutdemux_q30_unstranded_highdepth_qcfltd <- fread(file = "Report/SNVs/filter/basediffperc_cutdemux_sub500k_q30_unstranded_highdepth_qcfltd.csv.gz")
basediffperc_cutdemux_q30_unstranded_highdepth_qcfltd[, alt := factor(alt, levels = c("A", "C", "G", "T", "del"))]
basediffperc_cutdemux_q30_unstranded_highdepth_qcfltd[, mut := factor(paste0(ref, ">", alt), levels = paste0(rep(c("A", "C", "G", "T"), each = 5), ">", c("A", "C", "G", "T", "del")))]
setkey(basediffperc_cutdemux_q30_unstranded_highdepth_qcfltd, pos, alt)
basediffperc_cutdemux_q30_unstranded_highdepth_qcfltd[, posmut := paste0(pos, ":", mut)]
basediffperc_cutdemux_q30_unstranded_highdepth_qcfltd[, posmut := factor(posmut, levels = unique(posmut))]
dim(basediffperc_cutdemux_q30_unstranded_highdepth_qcfltd)
## [1] 1021830      28

basediffperc_cutdemux_q30_unstranded_highdepth_qcfltd_altperc_bymito_byposmut <- dcast.data.table(basediffperc_cutdemux_q30_unstranded_highdepth_qcfltd, LibraryMitoID ~ posmut, value.var = "altperc")
dim(basediffperc_cutdemux_q30_unstranded_highdepth_qcfltd_altperc_bymito_byposmut)
## [1] 1702 1033

basediffperc_cutdemux_q30_unstranded_highdepth_qcfltd_altperc_bymito_byposmut <- merge.data.table(MitoInfo, basediffperc_cutdemux_q30_unstranded_highdepth_qcfltd_altperc_bymito_byposmut, by = "LibraryMitoID", all.x = FALSE, all.y = TRUE)
dim(basediffperc_cutdemux_q30_unstranded_highdepth_qcfltd_altperc_bymito_byposmut)
## [1] 1702 1051

fwrite(basediffperc_cutdemux_q30_unstranded_highdepth_qcfltd_altperc_bymito_byposmut, file = "Report/SNVs/filter/basediffperc_cutdemux_sub500k_q30_unstranded_highdepth_qcfltd_altperc_bymito_byposmut.csv.gz")
basediffperc_cutdemux_q30_unstranded_highdepth_qcfltd_altperc_bymito_byposmut <- fread(file = "Report/SNVs/filter/basediffperc_cutdemux_sub500k_q30_unstranded_highdepth_qcfltd_altperc_bymito_byposmut.csv.gz")

###########################################################################
## for all QC filtered SNV sites, get the allele table: mitos X sites 
## (A, AG, etc. whatever bases >= 5%)
###########################################################################
basediffperc_cutdemux_q30_unstranded_highdepth <- fread("Report/SNVs/filter/basediffperc_cutdemux_sub500k_q30_unstranded_highdepth.csv.gz")
basediffperc_cutdemux_q30_unstranded_highdepth_highaf_qcfltd_support_byposmut <- fread(file = "Report/SNVs/filter/basediffperc_cutdemux_sub500k_q30_unstranded_highdepth_highaf_qcfltd_support_byposmut.csv")
basediffperc_cutdemux_q30_unstranded_highdepth_qcfltd <- basediffperc_cutdemux_q30_unstranded_highdepth[pos %in% basediffperc_cutdemux_q30_unstranded_highdepth_highaf_qcfltd_support_byposmut[, pos]]
basediffperc_cutdemux_q30_unstranded_highdepth_qcfltd[ref == "A", A := `=`]
basediffperc_cutdemux_q30_unstranded_highdepth_qcfltd[ref == "C", C := `=`]
basediffperc_cutdemux_q30_unstranded_highdepth_qcfltd[ref == "G", G := `=`]
basediffperc_cutdemux_q30_unstranded_highdepth_qcfltd[ref == "T", T := `=`]

vaf_th <- 0.05
basediffperc_cutdemux_q30_unstranded_highdepth_qcfltd[, allele := apply(.SD, 1, function(x) paste0(c("A", "C", "G", "T", "-")[x >= vaf_th * 100], collapse = "")), .SDcols = c("A", "C", "G", "T", "del")]
basediffperc_cutdemux_q30_unstranded_highdepth_qcfltd[, table(allele)]
## allele
##      A     A-     AC    AC-    ACG   ACG-   ACGT    ACT     AG    AGT     AT
## 298972     49    855     57     15      2      1      2   1439     25    505
##      C     C-     CG    CGT     CT      G     GT      T     T-
## 204063     28    339      8   1015 100040    853 224601      3
basediffperc_cutdemux_q30_unstranded_highdepth_qcfltd_allele_bymito_bypos <- dcast.data.table(basediffperc_cutdemux_q30_unstranded_highdepth_qcfltd, LibraryMitoID ~ pos, value.var = "allele")
dim(basediffperc_cutdemux_q30_unstranded_highdepth_qcfltd_allele_bymito_bypos)
## [1] 1702  839

basediffperc_cutdemux_q30_unstranded_highdepth_qcfltd_allele_bymito_bypos <- merge.data.table(MitoInfo, basediffperc_cutdemux_q30_unstranded_highdepth_qcfltd_allele_bymito_bypos, by = "LibraryMitoID", all.x = FALSE, all.y = TRUE)
dim(basediffperc_cutdemux_q30_unstranded_highdepth_qcfltd_allele_bymito_bypos)
## [1] 1702  857

fwrite(basediffperc_cutdemux_q30_unstranded_highdepth_qcfltd_allele_bymito_bypos, file = "Report/SNVs/filter/basediffperc_cutdemux_sub500k_q30_unstranded_highdepth_qcfltd_allele_bymito_bypos.csv.gz")
basediffperc_cutdemux_q30_unstranded_highdepth_qcfltd_allele_bymito_bypos <- fread(file = "Report/SNVs/filter/basediffperc_cutdemux_sub500k_q30_unstranded_highdepth_qcfltd_allele_bymito_bypos.csv.gz")

###########################################################################
## for all QC filtered SNVs, get the depth table: mito X SNV site
## Note, we need to recover the low-depth sites, so we have to go back to
## the original table before depth filtering. 
###########################################################################
basedifffreq_cutdemux_q30_unstranded <- fread(file = "Report/SNVs/basedifffreq_cutdemux_sub500k_q30_unstranded_9027fixed.csv.gz")
basediffperc_cutdemux_q30_unstranded_highdepth_highaf_qcfltd_support_byposmut <- fread(file = "Report/SNVs/filter/basediffperc_cutdemux_sub500k_q30_unstranded_highdepth_highaf_qcfltd_support_byposmut.csv")
basedifffreq_cutdemux_q30_unstranded_qcfltd <- basedifffreq_cutdemux_q30_unstranded[pos %in% basediffperc_cutdemux_q30_unstranded_highdepth_highaf_qcfltd_support_byposmut[, unique(pos)]]
basedifffreq_cutdemux_q30_unstranded_qcfltd_depth_bymito_bypos <- dcast.data.table(basedifffreq_cutdemux_q30_unstranded_qcfltd, LibraryMitoID ~ pos, value.var = "depth", fill = 0)
basedifffreq_cutdemux_q30_unstranded_qcfltd_depth_bymito_bypos <- MitoInfo[basedifffreq_cutdemux_q30_unstranded_qcfltd_depth_bymito_bypos, on = "LibraryMitoID"]
fwrite(basedifffreq_cutdemux_q30_unstranded_qcfltd_depth_bymito_bypos, file = "Report/SNVs/filter/basedifffreq_cutdemux_sub500k_q30_unstranded_qcfltd_depth_bymito_bypos.csv.gz")
basedifffreq_cutdemux_q30_unstranded_qcfltd_depth_bymito_bypos <- fread(file = "Report/SNVs/filter/basedifffreq_cutdemux_sub500k_q30_unstranded_qcfltd_depth_bymito_bypos.csv.gz", header = TRUE)

###########################################################################
## filter insersions
###########################################################################
depth_th <- 50
vaf_th <- 0.05

mpileups_ins_cutdemux <- read.csv("Report/SNVs/ins/mpileups_ins_cutdemux_sub500k.csv.gz", as.is = TRUE, check.names = FALSE)
setDT(mpileups_ins_cutdemux)
mpileups_ins_cutdemux_highdepth_highaf <- mpileups_ins_cutdemux[depth >= depth_th & altperctot >= vaf_th * 100]
dim(mpileups_ins_cutdemux_highdepth_highaf)
## [1] 123  34
## [1] 134  34 when E.799 replaced the old mitos
## [1] 134  34 after filtering out technical and digested controls
fwrite(mpileups_ins_cutdemux_highdepth_highaf, file = "Report/SNVs/filter/mpileups_ins_cutdemux_sub500k_highdepth_highaf.csv.gz", na = "NA")
mpileups_ins_cutdemux_highdepth_highaf <- fread("Report/SNVs/filter/mpileups_ins_cutdemux_sub500k_highdepth_highaf.csv.gz")

chrmproperties <- fread("Report/artifact/chrmbases_properties.csv.gz")
mpileups_ins_cutdemux_highdepth_highaf_qcfltd <- mpileups_ins_cutdemux_highdepth_highaf[pos %in% chrmproperties[is_in_range == "Y" & is_in_primer == "N", pos]]
fwrite(mpileups_ins_cutdemux_highdepth_highaf_qcfltd, file = "Report/SNVs/filter/mpileups_ins_cutdemux_sub500k_highdepth_highaf_qcfltd.csv.gz")
