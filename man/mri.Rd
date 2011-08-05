\name{mri}
\docType{data}
\alias{mri}
\title{Diagnosis of lymph node metastasis with magnetic resonance imaging}
\description{
This data frame corresponds to 10 clinical studies reporting the accuracy of magnetic resonance
imaging (MRI) for the diagnosis of lymph  node metastasis in women with cervical cancer
(Table III of Scheidler et al. 1997).
}

\usage{mri}
\format{A matrix with 10 rows and 4 columns. Each row represents study results, the columns
are:
\describe{
	\item{\code{tp}}{true positives}
	\item{\code{n1}}{number of patients with disease}
	\item{\code{fp}}{false positives}
	\item{\code{n2}}{number of patients without disease}
}
}
\source{
The data were obtainded from

Scheidler J, Hricak H, Yu KK, Subak L, Segal MR. (1997)
Radiological evaluation of lymph node metastases in patients with cervical cancer: a
meta-analysis.
\emph{The Journal of the American Medical Association};
\bold{278}:1096-1101.}

\references{
Verde P. E. (2010). Meta-analysis of diagnostic test data: A bivariate
Bayesian modeling approach. \emph{Statistics in Medicine}. \bold{29},
3088-3102.}

\keyword{datasets}
