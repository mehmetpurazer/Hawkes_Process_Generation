% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/NHPP_Inv.R
\name{NHPP.inversion}
\alias{NHPP.inversion}
\title{Non-Homogeneous Possion Process Simulation with Inversion (Transformation)}
\usage{
NHPP.inversion(T, lam)
}
\arguments{
\item{lam}{The intensity function of NHPP}
}
\description{
This function simulates the event times of Non-Homogeneous Possion Process (NHHP).

It uses Inversion (Transformation) Algorithm
}
\examples{
NHPP.inversion(T=20,lam = function(x) {0*x+1})
NHPP.inversion(T=20,lam = function(x) {3.9*exp(-4*x)})
NHPP.inversion(T=20,lam = function(x) {2000*(x+87)^(-2.33)})

}
