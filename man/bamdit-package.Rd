% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bamdit-package.R
\docType{package}
\name{bamdit-package}
\alias{bamdit-package}
\alias{bamdit}
\title{Bayesian Meta-Analysis of Diagnostic Test Data}
\description{
Bayesian meta-analysis of diagnostic test data based on a scale mixtures
bivariate random-effects model.
This package was developed with the aim of simplifying the use of meta-analysis
models that up to now have demanded great statistical expertise in Bayesian meta-analysis.
The package implements a series of innovative statistical techniques including:
the BSROC (Bayesian Summary ROC) curve, the BAUC (Bayesian AUC), predictive surfaces,
the use of prior distributions  that avoid boundary estimation problems of component
of variance and correlation parameters, analysis of conflict of evidence and robust
estimation of model parameters. In addition, the package comes with several published
examples of meta-analysis that can be used for illustration or further research in
this area.
}
\details{
\tabular{ll}{
Package: \tab bamdit    \cr
Type: \tab Package      \cr
Version: \tab 3.4.0     \cr
Date: \tab 2022-04-04  \cr
License: \tab GPL (>= 2)\cr
LazyLoad: \tab yes\cr }
}
\references{
Verde P. E. (2010). Meta-analysis of diagnostic test data: A
bivariate Bayesian modeling approach. Statistics in Medicine. 29(30):3088-102.
doi: 10.1002/sim.4055.

Verde P. E. (2018). bamdit: An R Package for Bayesian Meta-Analysis
of Diagnostic Test Data. Journal of Statistical Software. Volume 86, issue 10, pages 1--32.
}
\author{
PD Dr. Pablo Emilio Verde \email{pabloemilio.verde@hhu.de}
}
\keyword{AUC}
\keyword{JAGS}
\keyword{MCMC}
\keyword{Meta-Analysis}
\keyword{Outliers}
\keyword{ROC}
\keyword{Sensitivity}
\keyword{Specificity}
