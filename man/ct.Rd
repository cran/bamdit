\name{ct}
\docType{data}
\alias{ct}
\title{Diagnosis of appendicities with computer tomography scans}
\description{
This data frame corresponds to 51 clinical studies reporting the accuracy of
computer tomography (CT) scans for the diagnosis of appendicities.

}
\usage{ct}
\format{A matrix with 51 rows and 4 columns. Each row represents
study results, the columns are:
\describe{
	\item{\code{tp}}{true positives}
	\item{\code{n1}}{number of patients with disease}
	\item{\code{fp}}{false positives}
	\item{\code{n2}}{number of patients without disease}
}
}
\source{
Ohmann C, Verde PE, Gilbers T, Franke C, Fuerst G, Sauerland S, Boehner H.
(2006) Systematic review of CT investigation in suspected
acute appendicitis. \emph{Final Report; Coordination Centre for
Clinical Trials, Heinrich-Heine University}.
Moorenstr. 5, D-40225 Duesseldorf
Germany.
}
\references{
Verde P. E. (2010). Meta-analysis of diagnostic test data: A bivariate
Bayesian modeling approach. \emph{Statistics in Medicine}. \bold{29},
3088-3102.}
\keyword{datasets}
