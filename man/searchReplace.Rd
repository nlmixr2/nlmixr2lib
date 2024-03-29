% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/searchReplace.R
\name{searchReplace}
\alias{searchReplace}
\alias{searchReplaceHelper}
\title{Search within a model to replace part of the model}
\usage{
searchReplace(object, find, replace)

searchReplaceHelper(object, find, replace)
}
\arguments{
\item{object}{function specifying the nlmixr2 model}

\item{find, replace}{Character scalars of parts of the model to replace}
}
\value{
\code{object} with \code{find} replaced with \code{replace}
}
\description{
Search within a model to replace part of the model
}
\section{Functions}{
\itemize{
\item \code{searchReplaceHelper()}: A helper function for searchReplace (not intended
for users to use directly)

}}
\keyword{Internal}
