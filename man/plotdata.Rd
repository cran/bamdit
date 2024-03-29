% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotdata.R
\name{plotdata}
\alias{plotdata}
\title{Basic function to plot the data of meta-analysis of diagnostic test}
\usage{
plotdata(
  data,
  two.by.two = FALSE,
  group = NULL,
  x.lo = 0,
  x.up = 1,
  y.lo = 0,
  y.up = 1,
  alpha.p = 0.7,
  max.size = 15
)
}
\arguments{
\item{data}{Either a data frame with at least 4 columns containing the true positives (tp),
number of patients with disease (n1), false positives (fp), number of patients without
disease (n2), or for two.by.two = TRUE a data frame where each line contains the
diagnostic results as a two by two table, where the column names are:
TP, FP, TN, FN.}

\item{two.by.two}{If TRUE indicates that the diagnostic results are given as: TP, FP, TN, FN.}

\item{group}{a variable name indicating a group factor}

\item{x.lo}{lower limit of the x-axis}

\item{x.up}{upper limit of the x-axis}

\item{y.lo}{lower limit of the y-axis}

\item{y.up}{upper limit of the y-axis}

\item{alpha.p}{transparency of the points}

\item{max.size}{scale parameter of the maximum size}
}
\description{
This function plots the true positive rates vs the false positive rates of each study included
in the meta-analysis. Study results are displayed by circles, the diameter of each circle is proportional
to the sample size of the study (or table). If subgroups are displayed each group is represented by
different colours. This function use the package \emph{ggplot2}.
}
\examples{

## execute analysis
\dontrun{

data(ct)
ct$design <- with(ct, factor(design,
             labels = c("Prospective", "Retrospective")))

plotdata(ct,              # Data frame
        group = "design", # Groupping variable
        y.lo = 0.75,      # Lower limit of y-axis
        x.up = 0.75,      # Upper limit of x-axis
        alpha.p = 0.5,    # Transparency of the balls
        max.size = 5)     # Scale the circles
}


}
