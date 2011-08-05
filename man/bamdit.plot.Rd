\name{bamdit.plot}
\alias{bamdit.plot}
\title{The Bamdit Plot Function}
\description{
This function plots a scatterplot of posterior predictive sensitivity and specificity
with their marginal histograms. The original rates are circled blue.
}
\usage{
bamdit.plot(model,data=NULL,S=200)
}
\arguments{
	\item{data}{The data use for analysis which the model is generated.}
	\item{model}{The model created from metadiag().}
	\item{S}{The number of predictive posterior rates display in the plot.}
}
\seealso{
  \code{\link{metadiag}}.
}
\examples{
## execute analysis
data(mri)
bamdit.plot(model=metadiag(mri),data=mri)
}
\keyword{file}