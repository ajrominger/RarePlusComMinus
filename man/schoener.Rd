% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/schoener.R
\name{schoener}
\alias{schoener}
\title{Schoener similarity}
\usage{
schoener(x)
}
\arguments{
\item{x}{site by species matrix (sites as rows, species as columns, abundances in cells)}
}
\value{
an \code{ncol(x)} by \code{ncol(x)} matrix of species-species similarities
}
\description{
Calculate Schoener similarity metric
}
\details{
Faster implementation than what is availible in \code{spaa::niche.overlap}.
}
\author{
Andy Rominger <ajrominger@gmail.com>
}
