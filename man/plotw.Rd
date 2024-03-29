% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotw.R
\name{plotw}
\alias{plotw}
\title{Plot for the conflict of evidence parameters w1 and w2}
\usage{
plotw(
  m,
  group = NULL,
  title = "Posterior quantiles (25\%, 50\%, 75\%)",
  group.colors = c("blue", "red")
)
}
\arguments{
\item{m}{The object generated by metadiag. The model object must be fitted with the options: re = "sm" and split.w = TRUE.}

\item{group}{An optional argument which is a variable name indicating a group factor.
If set, then the plot is colored by groups.}

\item{title}{The title of the plot.}

\item{group.colors}{A character vector with two color names.}
}
\description{
Conflict of evidence plot: this plot displays the posterior distribution of the study's weights w1 and w1.
These weights indicate potential conflict of evidence of the studies. The weight w1 indicates deviations
with respect to the specificity and w2 to the sensitivity.
}
\examples{

## execute analysis
\dontrun{
data(ep)
ep$design = factor(ep$d1,labels = c("prospective", "retrospective"))
m.ep <- metadiag(ep, re = "sm", re.model = "SeSp",
                split.w = TRUE,
                df.estimate = TRUE)

plotw(m.ep)
#Relationship between conflict and study design
plotw(m.ep, group = "design")

}

}
\seealso{
\code{\link{metadiag}}.
}
\keyword{file}
