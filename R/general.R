#' Search within a model to replace part of the model
#' 
#' @inheritParams nlmixr2::nlmixr2
#' @return \code{object} with \code{find} replaced with \code{replace}
#' @export
searchReplace <- function(object, find, replace) {
  # ensure that it is an nlmixr2 object
  object <- nlmixr2::nlmixr2(object)
  if (is.character(find)) {
    find <- str2lang(find)
  }
  if (is.character(replace)) {
    replace <- str2lang(replace)
  }
  nlmixr2::nlmixr2(
    searchReplaceHelper(object = object$fun, find = find, replace = replace)
  )
}

#' @export
searchReplaceHelper <- function(object, find, replace) {
  UseMethod("searchReplaceHelper")
}

#' @export
searchReplaceHelper.function <- function(object, find, replace) {
  functionBody(object) <- searchReplaceHelper(object = functionBody(object), find = find, replace = replace)
  object
}

#' @export
searchReplaceHelper.call <- function(object, find, replace) {
  if (identical(object[[1]], as.name("ini"))) {
    # no replacement within ini()
    return(object)
  } else {
    if (identical(object, find)) {
      return(replace)
    } else {
      for (idx in seq_along(object)) {
        object[[idx]] <- searchReplaceHelper(object[[idx]], find, replace)
      }
    }
  }
  object
}

#' @export
searchReplaceHelper.default <- function(object, find, replace) {
  if (identical(object, find)) {
    message("found")
    object <- replace
  } else if (length(object) != 1) {
    for (idx in seq_along(object)) {
      object[[idx]] <- searchReplaceHelper(object[[idx]], find = find, replace = replace)
    }
  }
  # if the length is 1 and it is not identical, then do nothing
  object
}
