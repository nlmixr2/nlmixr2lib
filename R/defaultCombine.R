.combineEnv <- new.env(parent = emptyenv())
.combineEnv$default <- "camel"

#' Case combine two strings depending on .combineEnv$defualt
#'
#' @param v1 string 1
#' @param v2 string 2
#' @return defaultCase combined string
#' @noRd
#' @author Matthew L. Fidler
.defaultCombine2 <- function(v1, v2) {
  if (checkmate::testCharacter(v1, min.len=2, any.missing = FALSE)) {
    v1 <- do.call(defaultCombine, as.list(v1))
  }
  if (checkmate::testCharacter(v1, min.len=2, any.missing = FALSE)) {
    v2 <- do.call(defaultCombine, as.list(v2))
  }
  checkmate::assertCharacter(v1, len=1, any.missing = FALSE)
  checkmate::assertCharacter(v2, len=1, any.missing = FALSE)
  .n1 <-  nchar(v2)
  if (.n1 == 0L) {
    v1
  } else if (.n1 == 1L) {
    switch(.combineEnv$default,
           camel = paste0(v1, toupper(v2)),
           snake = paste0(v1, "_", v2),
           dot   = paste0(v1, ".", v2),
           blank = paste0(v1, v2))
  } else {
    switch(.combineEnv$default,
           camel = paste0(v1, toupper(substr(v2, 1, 1)), substr(v2, 2, .n1)),
           snake = paste0(v1, "_", v2),
           dot   = paste0(v1, ".", v2),
           blank = paste0(v1, v2))
  }
}
#' Default combine strings
#'
#' @param ... uses default to combine strings
#' @param type type of combine to use (set with `setCombineType()`)
#' @return combined strings
#' @export
#' @author Matthew L. Fidler
#' @examples
#'
#' # default combine
#'
#' defaultCombine("f", "depot")
#'
#' defaultCombine(list(c("a", "funny", "c")))
#'
#' defaultCombine(c("a", "funny", "c"))
#'
#' # snake combine
#'
#' snakeCombine("f", "depot")
#'
#' snakeCombine(list(c("a", "funny", "c")))
#'
#' snakeCombine(c("a", "funny", "c"))
#'
#' # dot combine
#'
#' dotCombine("f", "depot")
#'
#' dotCombine(list(c("a", "funny", "c")))
#'
#' dotCombine(c("a", "funny", "c"))
#'
#' # blank combine
#'
#' blankCombine("f", "depot")
#'
#' blankCombine(list(c("a", "funny", "c")))
#'
#' blankCombine(c("a", "funny", "c"))
#'
defaultCombine <- function(...) {
  .args <- list(...)
  .n <- length(.args)
  if (.n == 0L) {
    stop("no arguments provided",
         call.=FALSE)
  } else if (.n == 1L) {
    .ret <- .args[[1]]
    if (checkmate::testCharacter(.ret, len=1, any.missing=FALSE)) {
      .ret
    } else if (checkmate::testCharacter(.ret, min.len=2, any.missing=FALSE)) {
      do.call(defaultCombine, as.list(.ret))
    } else if (is.list(.ret)) {
      do.call(defaultCombine, .ret)
    } else {
      stop("invalid argument",
           call.=FALSE)
    }
  } else {
    Reduce(.defaultCombine2, .args)
  }
}
#' @describeIn defaultCombine
#' @export
snakeCombine <- function(...) {
  .combineEnv$old <- .combineEnv$default
  .combineEnv$default <- "snake"
  on.exit({.combineEnv$default <- .combineEnv$old})
  defaultCombine(...)
}
#' @describeIn defaultCombine
#' @export
camelCombine <- function(...) {
  .combineEnv$old <- .combineEnv$default
  .combineEnv$default <- "camel"
  on.exit({.combineEnv$default <- .combineEnv$old})
  defaultCombine(...)
}

#' @describeIn defaultCombine
#' @export
dotCombine <- function(...) {
  .combineEnv$old <- .combineEnv$default
  .combineEnv$default <- "dot"
  on.exit({.combineEnv$default <- .combineEnv$old})
  defaultCombine(...)
}

#' @describeIn defaultCombine
#' @export
blankCombine <- function(...) {
  .combineEnv$old <- .combineEnv$default
  .combineEnv$default <- "blank"
  on.exit({.combineEnv$default <- .combineEnv$old})
  defaultCombine(...)
}

#' @describeIn defaultCombine
#' @export
setCombineType <- function(type=c("snake", "camel", "dot", "blank")) {
  .combineEnv$default <- type
}