###########################################################################
## Author: Youtao Lu <luyoutao@sas.upenn.edu>
##  
## Copyright (c) 2021-2023, Youtao Lu and Junhyong Kim, Department of Biology, University of Pennsylvania
## Copyright (c) 2021-2023, James Eberwine, Perelman School of Medicine, University of Pennsylvania
## Copyright (c) 2021-2023, David Issadore, Penn Engineering, University of Pennsylvania
## 
## This Source Code is distributed under Creative Commons Attribution License 4.0 (CC BY).
###########################################################################

if (!exists("Tools") || is.environment(Tools)) Tools <- new.env(parent = emptyenv())
local({
    .VERSION = "0.3"

    write_xlsx <- function(dflist, filename, sheetnames = NULL, row.names = TRUE) {
        requireNamespace("openxlsx")
        wb <- openxlsx::createWorkbook()
        k <- 0
        if (is.null(sheetnames)) {
            if (is.null(names(dflist))) { stop("No names in data.frame list nor sheetnames provided!") }
            sheetnames <- names(dflist)
        }
        for (x in sheetnames) {
            k <- k + 1
            openxlsx::addWorksheet(wb, sheetName = substr(x, 1, 31))
            df <- dflist[[k]]
            if (!is.data.frame(df)) { df <- as.data.frame(df) }
            openxlsx::writeDataTable(wb, sheet = k, x = df, keepNA = TRUE, headerStyle = NULL, tableStyle = "none", bandedRows = FALSE, bandedCols = FALSE, withFilter = TRUE, rowNames = row.names)
        }
        openxlsx::saveWorkbook(wb, file = filename, overwrite = TRUE)
    }

    read_xlsx <- function(filename, sheetnames, ...) {
        requireNamespace("openxlsx")
        dfs <- lapply(sheetnames, function(sheet) { openxlsx::read.xlsx(xlsxFile = filename, sheet = sheet, ...) })
        names(dfs) <- sheetnames
        dfs
    }

    for (.obj in ls()) {
        assign(.obj, get(.obj), envir = Tools)
    }
})
