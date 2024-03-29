% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/addTransit.R
\name{addTransit}
\alias{addTransit}
\title{To add transit compartments to the model}
\usage{
addTransit(
  model,
  transit,
  central = "central",
  depot = "depot",
  transitComp = "transit",
  ktr = "ktr"
)
}
\arguments{
\item{model}{The model as a function}

\item{transit}{the number of transit compartments to be added}

\item{central}{a character vector representing the central
compartment}

\item{depot}{a character vector representing the depot compartment}

\item{transitComp}{the transit compartment prefix}

\item{ktr}{the parameter name for the transit compartment rate}
}
\value{
a model with transit compartment added
}
\description{
To add transit compartments to the model
}
\examples{
readModelDb("PK_1cmt_des") |>
  addTransit(3)
}
