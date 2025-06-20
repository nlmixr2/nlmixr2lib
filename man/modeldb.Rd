% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/modeldb.R
\docType{data}
\name{modeldb}
\alias{modeldb}
\title{Model library for nlmixr2}
\format{
A data frame with 53 rows and 9 columns
\describe{
  \item{name}{Model name that can be used to extract the model from the model library}
  \item{description}{Model description in free from text; in model itself}
  \item{parameters}{A comma separated string listing either the parameter in the model defined by population/individual effects or a population effect parameter}
  \item{DV}{The definition of the dependent variable(s)}
  \item{linCmt}{Logical flag indicating if solved models are used (TRUE) or not (FALSE)}
  \item{algebraic}{Logical flag indicating if the model is purely algebraic: TRUE no linCmt() and no ODEs; FALSE otherwise}
  \item{dosing}{A comma separated string of identified dosing compartments}
  \item{depends}{A comma separated string of objects the model depends on}
  \item{filename}{Filename of the model.  By default these are installed in the model library and read on demand}
}
}
\usage{
modeldb
}
\description{
This is a data frame of the available models in nlmixr2lib, it is generated
with the package.  Custom modeldb may be used.
}
\keyword{datasets}
