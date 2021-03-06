% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/alpha_diversity.R
\name{hill_number}
\alias{hill_number}
\title{hill_number}
\usage{
hill_number(physeq, c(1:3))
}
\arguments{
\item{phyloseq}{,  A \code{phyloseq} object containing merged information of abundance,
sample data including the measured variables and categorical information of the samples.}

\item{q}{, A vevtor eg c(0:1), c(0:3)}

\item{samID}{, A character indicate the sample ID is in row or col; include "row" or "col"}
}
\value{
data.frame
}
\description{
calculate the hill number
}
\details{
10/01/2020  ShenZhen China
}
\examples{
data(physeq_data)
hill_number(physeq_data, c(1:3))
}
\author{
(s)   Huahui Ren
}
