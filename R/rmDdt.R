#' Remove all the d/dt(cmts) in modelLines
#'
#' @param modelLines list of model lines to modify
#' @param cmts compartment names to remove
#' @return modelLines with compartments removed
#' @noRd
#' @author Matthew L. Fidler
.rmDdt <- function(modelLines, cmts) {
  .w <- which(vapply(seq_along(modelLines),
                     function(i) {
                       .cur <- modelLines[[i]]
                       any(vapply(cmts,
                                  function(cmt) {
                                    .ddtCentral1 <- str2lang(paste0("d/dt(",
                                                                    cmt, ") <- ."))
                                    .ddtCentral2 <- str2lang(paste0("d/dt(",
                                                                    cmt, ") = ."))
                                    rxode2::.matchesLangTemplate(.cur, .ddtCentral1) ||
                                      rxode2::.matchesLangTemplate(.cur, .ddtCentral2)
                                  }, logical(1), USE.NAMES = FALSE))

                     }, logical(1), USE.NAMES = FALSE))
  lapply(seq_along(modelLines)[-.w],
         function(i) {
           modelLines[[i]]
         })
}
