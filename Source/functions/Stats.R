###########################################################################
## Author: Youtao Lu <luyoutao@sas.upenn.edu>
##  
## Copyright (c) 2021-2023, Youtao Lu and Junhyong Kim, Department of Biology, University of Pennsylvania
## Copyright (c) 2021-2023, James Eberwine, Perelman School of Medicine, University of Pennsylvania
## Copyright (c) 2021-2023, David Issadore, Penn Engineering, University of Pennsylvania
## 
## This Source Code is distributed under Creative Commons Attribution License 4.0 (CC BY).
###########################################################################

if (!exists("Stats") || is.environment(Stats)) Stats <- new.env(parent = emptyenv()) 
local({
    .VERSION = "0.1"

    parcsin <- function(x) {
        if (any(!is.na(x) & x < 0) || any(!is.na(x) & x > 100)) { 
            stop("Input is expected to be a percentage!")
        }
        asin(sqrt(x / 100))
    }

    for (obj in ls()) 
        assign(obj, get(obj), envir = Stats)
})
