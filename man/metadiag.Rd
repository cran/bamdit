\name{metadiag}
\alias{metadiag}
\title{Bayesian Meta-Analysis of diagnostic test data}
\description{
This function performers a Bayesian meta-analysis of diagnostic test
data by fitting a bivariate random effects model. The number of true
positives and false positives are modeled with two conditional Binomial
distributions and the random-effects are based on a bivariate scale
mixture of Normals. Computations are done by calling (Just Another Gibbs Sampler) to perform
MCMC (Markov Chain Monte Carlo) sampling and returning an object of
the class \emph{mcmc.list}. }
\usage{
metadiag(data, re="normal", link="logit", R = matrix(c(1, 0, 0, 1),
byrow=TRUE, nrow=2), m.0 = c(0,0), pre.mu = c(0.25, 0.25), k = 3, 
nu.0 = 1, n.chains = 3, n.iter = 40000, n.burnin = 20000, 
n.thin = 10, verbose = TRUE)
}
\arguments{
 \item{data}{A data frame with 4 columns containing the true positives,
  number of patients with disease, false positives, number of patients
  without disease}
	\item{re}{Random effects distribution for the resulting model. Possible values are \emph{normal} and \emph{scalemix}}
	\item{link}{The link function used in the model. Possible values are \emph{logit} and \emph{cloglog}. }
  \item{R}{ A 2x2 matrix parameters for the prior of Lambda, default value is matrix(c(1, 0, 0, 1), byrow=TRUE, nrow=2)}
	\item{m.0}{A two dimensional vector for the priors of mu, where the default values 
	are c(0, 0)}
	\item{pre.mu}{A two dimensional vector for the precision of mu, where the defaults 
	are c(0.25, 0.25)}
	\item{k}{The degrees of freedom for the Wishart prior of Lambda. Default: 3}
	\item{nu.0}{The parameter for the prior of nu, default value is nu.0=1}
	\item{n.chains}{Number of chains for the Models. default: 3}
	\item{n.iter}{Number of Iterations. iter-burnin = obtained values. default: 40k}
	\item{n.burnin}{Number of Iterations that are skipped at the beginning of the simulation. Default: 20k}
	\item{n.thin}{Default: 10}
	\item{verbose}{Default: TRUE}
}

\examples{
## execute analysis
data(mri)
mri
md <- metadiag(mri, re = "scalemix")
summary(md)
}
\keyword{file}