% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MCD.R
\name{MCD}
\alias{MCD}
\title{Capture \emph{MCD} value}
\usage{
MCD(x, maxgap = 8)
}
\arguments{
\item{maxgap}{max gap for calculate ratio when getting \emph{MCD} value, default is 8.}
}
\description{
return value of Month of Cyclical Dominance(MCD).
}
\examples{
x =seq(80, to = 120, length.out = 100)
MCD(x = x, maxgap = 20)


}
