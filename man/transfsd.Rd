% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{transfsd}
\alias{transfsd}
\title{Decorrelation}
\usage{
transfsd(x, X)
}
\arguments{
\item{x}{a matrix of nxp dimension. Input matrix to be transformed.}

\item{X}{a matrix of Nxp dimension. Supplementary data for covariance matrix estimation.}
}
\value{
Decorrelated data matrix.
}
\description{
Decorrelation
}
\details{
This function use both the input data and the supplementary data to estimate
the covariance matrix, and then use the estimated covariate matrix to decorrelate
the input data matrix. require (n+N)>p.
}
