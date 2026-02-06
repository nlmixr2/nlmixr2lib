#' Derive one-hot encoded race indicator columns for nirsevimab model
#'
#' Creates RACE_BLK_OTH and RACE_ASN_AMI_MUL binary columns from a character
#' race column, following the groupings in Clegg et al. (2024).
#'
#' @param data A data.frame containing a race column.
#' @param race_col Character string giving the name of the source race column.
#' @return The input data.frame with two additional integer columns.
#' @export
Clegg_2024_nirsevimab_derive_race_indicators <- function(data, race_col = "RACE") {
  data$RACE_BLK_OTH <-
    as.integer(data[[race_col]] %in% c("Black or African American", "Other"))
  data$RACE_ASN_AMI_MUL <-
    as.integer(data[[race_col]] %in% c("Asian", "American Indian or Alaskan Native", "Multiple"))
  data
}
