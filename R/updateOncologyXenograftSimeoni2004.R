#' Update an oncology xenograft model based on Simeoni 2004
#'
#' @inheritParams nlmixr2est::nlmixr2
#' @param ncmt The desired number of damaged cell compartments
#' @param damagedCmtName,undamagedCmtName,tumorVolName character string names
#'   for the compartments for damaged cells, undamaged cells, and the calcualted
#'   tumor volume (the sum of undamaged and damaged cells)
#' @param drugEffectName,transitRateName character string names of the drug effect and transit rate (as used in the model block)
#' @return An updated model with the new number of compartments
#' @examples
#' readModelDb("oncology_xenograft_simeoni_2004") %>%
#'   updateOncologyXenograftSimeoni2004(ncmt = 5)
#' @export
updateOncologyXenograftSimeoni2004 <- function(object, ncmt, damagedCmtName = "damagedCells", drugEffectName = "drugEffectCyclingCells", undamagedCmtName = "cyclingCells", tumorVolName = "tumorVol", transitRateName = "damageTransit") {
  checkmate::assert_integerish(ncmt, lower = 1, upper = 100, any.missing = FALSE, len = 1, null.ok = FALSE)
  if (is.function(object)) {
    # Convert a function to something able to be queried
    object <- rxode2::rxode(object)
  }
  # what compartments exist now?
  currentStates <- object$state
  # The tumor calculation will be calculated by addition.  Ensure that it is not
  # a compartment name.
  stopifnot(!(tumorVolName %in% currentStates))
  # The undamaged compartment will not be modified, but it will be an input to
  # the damaged compartments.  Ensure that it exists.
  stopifnot(undamagedCmtName %in% currentStates)
  newStates <- paste0(damagedCmtName, seq_len(ncmt))
  stopifnot(newStates[1] %in% currentStates)

  # Generate the new model lines
  tumorLine <- sprintf("%s <- %s", tumorVolName, paste(c(undamagedCmtName, newStates), collapse = " + "))
  damagedLines <-
    sprintf(
      "d/dt(%s) <- %s*%s - %s*%s",
      newStates[1],
      drugEffectName, undamagedCmtName,
      transitRateName, newStates[1]
    )
  if (ncmt > 1) {
    damagedLines <-
      c(
        damagedLines,
        sprintf(
          "d/dt(%s) <- %s*(%s - %s)",
          newStates[-1], transitRateName,
          newStates[-length(newStates)], newStates[-1]
        )
      )
  }
  dropOldLines <- sprintf("-d/dt(%s)", currentStates[startsWith(currentStates, damagedCmtName)])
  objectNoOld <- do.call(rxode2::model, append(list(object), lapply(X = dropOldLines, FUN = str2lang)))
  objectTumor <- do.call(rxode2::model, append(list(objectNoOld), lapply(X = tumorLine, FUN = str2lang)))
  objectDamage <- do.call(rxode2::model, append(list(objectNoOld, append = TRUE), lapply(X = damagedLines, FUN = str2lang)))
  objectDamage
}
