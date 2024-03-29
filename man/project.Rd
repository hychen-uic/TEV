% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{project}
\alias{project}
\title{Project x to the orthogonal complement of c.}
\usage{
project(cfd, x)
}
\arguments{
\item{cfd}{a matrix of nxq dimensional confounders}

\item{x}{a matrix of nxp dimensional exposures}
}
\value{
projtocorth, the projected x.
}
\description{
This function project exposures x to the orthogonal complement of confounders c.
}
\details{
This function uses c to form projection P=I-c(c'c)^(-)c'
and then compute projection z=Px
}
