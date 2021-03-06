% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{rnorm_truc}
\alias{rnorm_truc}
\title{Sample from truncated normal distribution}
\usage{
rnorm_truc(n, mean, sd, a, b)
}
\arguments{
\item{n}{sample size}

\item{mean}{the mean of the normal distribution}

\item{sd}{the standard deviation of the normal distribution}

\item{a}{the minimum value allowed}

\item{b}{the maximum value allowed}
}
\value{
a vector of size n
}
\description{
This function is borrowed from \code{\link[splatter]{splatter}}.
}
