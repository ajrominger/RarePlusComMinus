% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plusMinusIndSwap.R
\name{plusMinusIndSwap}
\alias{plusMinusIndSwap}
\title{Duplicates \code{plusMinus} but using the independent swap algorithm
from \{picante\}}
\usage{
plusMinusIndSwap(x, alpha = 0.05, B = 999)
}
\arguments{
\item{x}{site by species matrix (sites as rows, species as columns, abundances in cells)}

\item{alpha}{the significance level}

\item{B}{number of replicates for the null model permutations}
}
\value{
an \code{ncol(x)} by \code{ncol(x)} matrix of species-species similarities
}
\description{
Duplicates \code{plusMinus} but using the independent swap algorithm
from \{picante\}
}
\details{
All details are the same as \code{plusMinus} but this function uses
the independent swap algorithm from \{picante\} instead of the fixed-fixed
null algorithm. This function is only really intended to be used with
\code{indSwapTest}.
}
\author{
Andy Rominger <ajrominger@gmail.com>
}
