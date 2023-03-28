#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom data.table :=
#' @importFrom data.table .BY
#' @importFrom data.table .EACHI
#' @importFrom data.table .GRP
#' @importFrom data.table .I
#' @importFrom data.table .N
#' @importFrom data.table .NGRP
#' @importFrom data.table .SD
#' @importFrom data.table data.table
#' @importFrom data.table melt.data.table
#' @importFrom data.table dcast.data.table
#' @importFrom data.table as.data.table
#' @importFrom data.table merge.data.table
#' @importFrom data.table rbindlist
#' @importFrom magrittr %>%
#' @importFrom stats approx
#' @importFrom stats loess
#' @importFrom stats setNames
#' @importFrom stats predict
#' @importFrom stats rnorm
#' @importFrom stats runif
#' @importFrom stats na.omit
## usethis namespace: end
NULL

utils::globalVariables(
  c(
    ".",
    "SampleID",
    "abd",
    "ad",
    "ag",
    "amd",
    "ap",
    "aph",
    "bbd",
    "bbph",
    "bd",
    "bph",
    "comp",
    "variable",
    "value"
  )
)
