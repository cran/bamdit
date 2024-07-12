#' Bayesian Meta-Analysis of diagnostic test data
#'
#' This function performers a Bayesian meta-analysis of diagnostic test data by
#' fitting a bivariate random effects model. The number of true positives and
#' false positives are modeled with two conditional Binomial distributions and
#' the random-effects are based on a bivariate scale mixture of Normals.
#' Computations are done by calling JAGS (Just Another Gibbs Sampler) to perform
#' MCMC (Markov Chain Monte Carlo) sampling and returning an object of the
#' class \emph{mcmc.list}.
#'
#'
#' Installation of JAGS: It is important to note that R 3.3.0 introduced a major change in the
#' use of toolchain for Windows. This new toolchain is incompatible with older packages written in C++.
#' As a consequence, if the installed version of JAGS does not match the R installation, then the rjags
#' package will spontaneously crash. Therefore, if a user works with R version >= 3.3.0, then JAGS must
#' be installed with the installation program JAGS-4.2.0-Rtools33.exe. For users who continue using R 3.2.4 or
#' an earlier version, the installation program for JAGS is the default installer JAGS-4.2.0.exe.
#'
#'
#'
#' @param data            Either a data frame with at least 4 columns containing the true positives (tp),
#'                        number of patients with disease (n1), false positives (fp), number of patients without
#'                        disease (n2), or for two.by.two = TRUE a data frame where each line contains the
#'                        diagnostic results as a two by two table, where the column names are:
#'                        TP, FP, TN, FN.
#'
#' @param two.by.two      If TRUE indicates that the diagnostic results are given as: TP, FP, TN, FN.
#'
#' @param re              Random effects distribution for the resulting model. Possible
#'                        values are \emph{normal} for bivariate random effects and \emph{sm} for scale mixtures.
#'
#' @param link            The link function used in the model. Possible values are
#'                        \emph{logit}, \emph{cloglog} \emph{probit}.
#'
#' @param re.model        If re.model = "DS" indicates that the sum and differences of TPR and FPR are modeled as random effects and re.model = "SeSp" indicates that the Sensitivity and Specificity are modeled as ranodm effects.
#'                        The default value is re.model = "DS".
#'
#' @param mean.mu.D prior Mean of D, default value is 0.
#'
#' @param mean.mu.S prior Mean of S, default value is 0.
#'
#' @param sd.mu.D prior   Standard deviation of D, default value is 1 (the prior of mu.D is a logistic distribution).
#'
#' @param sd.mu.S prior   Standard deviation of S, default value is 1 (the prior of mu.S is a logistic distribution).
#'
#' @param sigma.D.upper   Upper bound of the uniform prior of sigma.S, default value is 10.
#'
#' @param sigma.S.upper   Upper bound of the uniform prior of sigma.S, default value is 10.
#'
#' @param mean.Fisher.rho Mean of rho in the Fisher scale default value is 0.
#'
#' @param sd.Fisher.rho   Standard deviation of rho in the Fisher scale, default value is 1/sqrt(2).
#'
#' @param df              If de.estimate = FALSE, then df is the degrees of freedom for the scale mixture distribution, default value is 4.
#'
#' @param df.estimate     Estimate the posterior of df. The defualt value is FALSE.
#'
#' @param df.lower        Lower bound of the prior of df. The defulat value is 3.
#'
#' @param df.upper        Upper bound of the prior of df. The defulat value is 30.
#'
#' @param split.w         Split the w parameter in two independent weights one for each random effect. The default value is FALSE.
#'
#' @param n.1.new         Number of patients with disease in a predictive study default is 50.
#'
#' @param n.2.new         Number of patients with non-disease in a predictive study default is 50.
#'
#' @param nr.chains       Number of chains for the MCMC computations, default 5.
#'
#' @param nr.iterations   Number of iterations after adapting the MCMC, default is 10000. Some models may need more iterations.
#'
#' @param nr.adapt        Number of iterations in the adaptation process, default is 1000. Some models may need more iterations during adaptation.
#'
#' @param nr.burnin       Number of iteration discard for burn-in period, default is 1000. Some models may need a longer burnin period.
#'
#' @param nr.thin         Thinning rate, it must be a positive integer, the default value 1.
#'
#' @param be.quiet        Do not print warning message if the model does not adapt default value is FALSE. If you are not sure about the adaptation period choose be.quiet=TRUE.
#'
#' @param r2jags          Which interface is used to link R to JAGS (rjags and R2jags) default value is R2Jags TRUE.
#'
#'
#'
#'
#' @return This function returns an object of the class metadiag. This object contains the MCMC output of
#' each parameter and hyper-parameter in the model, the data frame used for fitting the model, the link function,
#' type of random effects distribution and the splitting information for conflict of evidence analysis.
#'
#' The results of the object of the class metadiag can be extracted with R2jags or with rjags. In addition
#' a summary, a print and a plot functions are implemented for this type of object.
#'
#'
#' @references Verde P. E. (2010). Meta-analysis of diagnostic test data: A
#' bivariate Bayesian modeling approach. Statistics in Medicine. 29(30):3088-102.
#' doi: 10.1002/sim.4055.
#'
#'
#' @references Verde P. E. (2018). bamdit: An R Package for Bayesian Meta-Analysis
#' of Diagnostic Test Data. Journal of Statisticsl Software. Volume 86, issue 10, pages 1--32.
#'
#'
#' @examples
#'
#'
#' \dontrun{
#'
#' # Example: data from Glas et al. (2003).....................................
#' library(bamdit)
#' data("glas")
#' glas.t <- glas[glas$marker == "Telomerase", 1:4]
#'
#' glas.t <- glas[glas$marker == "Telomerase", 1:4]
#'
#' # Simple visualization ...
#'
#' plotdata(glas.t,                # Data frame
#'          two.by.two = FALSE     # Data is given as: (tp, n1, fp, n2)
#'          )
#'
#' glas.m1 <- metadiag(glas.t,                # Data frame
#'                     two.by.two = FALSE,    # Data is given as: (tp, n1, fp, n2)
#'                     re = "normal",         # Random effects distribution
#'                     re.model = "DS",       # Random effects on D and S
#'                     link = "logit",        # Link function
#'                     sd.Fisher.rho   = 1.7, # Prior standard deviation of correlation
#'                     nr.burnin = 1000,      # Iterations for burnin
#'                     nr.iterations = 10000, # Total iterations
#'                     nr.chains = 2,         # Number of chains
#'                     r2jags = TRUE)         # Use r2jags as interface to jags
#'
#'
#'  summary(glas.m1, digit=3)
#'
#'  plot(glas.m1,                    # Fitted model
#'       level = c(0.5, 0.75, 0.95), # Credibility levels
#'       parametric.smooth = TRUE)   # Parametric curve
#'
#'
#'# Plot results: based on a non-parametric smoother of the posterior predictive rates .......
#'
#' plot(glas.m1,                    # Fitted model
#'      level = c(0.5, 0.75, 0.95), # Credibility levels
#'      parametric.smooth = FALSE)  # Non-parametric curve
#'
#'
#'# Using the pipe command in the package dplyr ...............................................
#'
#'library(dplyr)
#'
#'glas.t %>%
#'   metadiag(re = "normal", re.model ="SeSp") %>%
#'   plot(parametric.smooth = FALSE, color.pred.points = "red")
#'
#'
#'
#'# Visualization of posteriors of hyper-parameters .........................................
#'library(ggplot2)
#'library(GGally)
#'library(R2jags)
#'attach.jags(glas.m1)
#'hyper.post <- data.frame(mu.D, mu.S, sigma.D, sigma.S, rho)
#'ggpairs(hyper.post,                  # Data frame
#'        title = "Hyper-Posteriors",          # title of the graph
#'        lower = list(continuous = "density") # contour plots
#'        )
#'
#'
#' #............................................................................
#'
#'# List of different statistical models:
#'#    1) Different link functions: logit, cloglog and probit
#'
#'#    2) Different parametrization of random effects in the link scale:
#'#         DS = "differences of TPR and FPR"
#'#         SeSp = "Sensitivity and Specificity"
#'
#'#    3) Different random effects distributions:
#'#       "normal" or "sm = scale mixtures".
#'
#'#    4) For the scale mixture random effects:
#'#       split.w = TRUE => "split the weights".
#'
#'#    5) For the scale mixture random effects:
#'#       df.estimate = TRUE => "estimate the degrees of freedom".
#'
#'#    6) For the scale mixture random effects:
#'#       df.estimate = TRUE => "estimate the degrees of freedom".
#'
#'#    7) For the scale mixture random effects:
#'#       df = 4 => "fix the degrees of freedom to a particual value".
#'#       Note that df = 1 fits a Cauchy bivariate distribution to the random effects.
#'
#'# logit-normal-DS
#' m <- metadiag(glas.t, re = "normal", re.model = "DS", link = "logit")
#' summary(m)
#' plot(m)
#'
#'# cloglog-normal-DS
#' summary(metadiag(glas.t, re = "normal", re.model = "DS", link = "cloglog"))
#'
#'# probit-normal-DS
#' summary(metadiag(glas.t, re = "normal", re.model = "DS", link = "probit"))
#'# logit-normal-SeSp
#' summary(metadiag(glas.t, re = "normal", re.model = "SeSp", link = "logit"))
#'
#'# cloglog-normal-SeSp
#' summary(metadiag(glas.t, re = "normal", re.model = "SeSp", link = "cloglog"))
#'# probit-normal-SeSp
#' summary(metadiag(glas.t, re = "normal", re.model = "SeSp", link = "probit"))
#'
#'# logit-sm-DS
#' summary(metadiag(glas.t, re = "sm", re.model = "DS", link = "logit", df = 1))
#'
#'# cloglog-sm-DS
#' summary(m<-metadiag(glas.t, re = "sm", re.model = "DS", link = "cloglog", df = 1))
#' plot(m, parametric.smooth = FALSE)
#'
#'# probit-sm-DS
#' summary(m<-metadiag(glas.t, re = "sm", re.model = "DS", link = "probit", df = 1))
#' plot(m, parametric.smooth = FALSE)
#'
#'# logit-sm-SeSp
#' summary(m<-metadiag(glas.t, re = "sm", re.model = "SeSp", link = "logit", df = 1))
#' plot(m, parametric.smooth = FALSE, level = c(0.5, 0.9))
#'
#'# cloglog-sm-SeSp
#' summary(m<-metadiag(glas.t, re = "sm", re.model = "SeSp", link = "cloglog", df = 1))
#' plot(m, parametric.smooth = FALSE, level = c(0.5, 0.9))
#'
#'# probit-sm-SeSp
#' summary(m<-metadiag(glas.t, re = "sm", re.model = "SeSp", link = "probit", df = 1))
#' plot(m, parametric.smooth = FALSE, level = c(0.5, 0.9))
#'
#'# logit-sm-DS-df
#' summary(m<-metadiag(glas.t, re = "sm", re.model = "DS", link = "logit",
#'  df.estimate = TRUE))
#' plot(m, parametric.smooth = FALSE, level = c(0.5, 0.9))
#'
#'# cloglog-sm-DS-df
#' summary(m<-metadiag(glas.t, re = "sm", re.model = "DS", link = "cloglog",
#' df.estimate = TRUE))
#' plot(m, parametric.smooth = FALSE, level = c(0.5, 0.9))
#'
#'# probit-sm-DS-df
#' summary(m<-metadiag(glas.t, re = "sm", re.model = "DS", link = "probit",
#' df.estimate = TRUE))
#' plot(m, parametric.smooth = FALSE, level = c(0.5, 0.9))
#'
#'# logit-sm-SeSp-df
#' summary(m<-metadiag(glas.t, re = "sm", re.model = "SeSp", link = "probit",
#' df.estimate = TRUE))
#' plot(m, parametric.smooth = FALSE, level = c(0.5, 0.9))
#'
#'# cloglog-sm-SeSp-df
#' summary(m<-metadiag(glas.t, re = "sm", re.model = "SeSp", link = "cloglog",
#' df.estimate = TRUE))
#' plot(m, parametric.smooth = FALSE, level = c(0.5, 0.9))
#'
#'# probit-sm-SeSp-df
#' summary(m<-metadiag(glas.t, re = "sm", re.model = "SeSp", link = "probit",
#' df.estimate = TRUE))
#' plot(m, parametric.smooth = FALSE, level = c(0.5, 0.9))
#'
#'# split.w ...................................................................
#'
#'# logit-sm-DS
#' summary(m <- metadiag(glas.t, re = "sm", re.model = "DS", link = "logit", split.w = TRUE, df = 10))
#' plot(m)
#'
#'# cloglog-sm-DS
#' summary(m<-metadiag(glas.t, re = "sm", re.model = "DS", link = "cloglog", split.w = TRUE, df = 4))
#' plot(m)
#'
#'# probit-sm-DS
#' summary(m<-metadiag(glas.t, re = "sm", re.model = "DS", link = "probit", split.w = TRUE, df = 4))
#' plot(m, parametric.smooth = FALSE)
#'
#'# logit-sm-SeSp
#' summary(m<-metadiag(glas.t, re = "sm", re.model = "SeSp", link = "logit", split.w = TRUE, df = 1))
#' plot(m, parametric.smooth = FALSE, level = c(0.5, 0.9))
#' plotw(m)
#'
#'# cloglog-sm-SeSp
#' summary(m<-metadiag(glas.t, re = "sm", re.model = "SeSp", link = "cloglog", split.w = TRUE, df = 1))
#' plot(m, parametric.smooth = FALSE, level = c(0.5, 0.9))
#' plotw(m)
#'
#'# probit-sm-SeSp
#' summary(m<-metadiag(glas.t, re = "sm", re.model = "SeSp", link = "probit", split.w = TRUE, df = 1))
#' plot(m, parametric.smooth = FALSE, level = c(0.5, 0.9))
#' plotw(m)
#'
#'
#'# logit-sm-DS-df
#' summary(m<-metadiag(glas.t, re = "sm", re.model = "DS", link = "logit", split.w = TRUE,
#'  df.estimate = TRUE))
#' plot(m, parametric.smooth = FALSE, level = c(0.5, 0.9))
#' plotw(m)
#'
#'# cloglog-sm-DS-df
#' summary(m<-metadiag(glas.t, re = "sm", re.model = "DS", link = "cloglog", split.w = TRUE,
#' df.estimate = TRUE))
#' plot(m, parametric.smooth = FALSE, level = c(0.5, 0.9))
#' plotw(m)
#'
#'# probit-sm-DS-df
#' summary(m<-metadiag(glas.t, re = "sm", re.model = "DS", link = "probit", split.w = TRUE,
#' df.estimate = TRUE))
#' plot(m, parametric.smooth = FALSE, level = c(0.5, 0.9))
#' plotw(m)
#'
#'# logit-sm-SeSp-df
#' summary(m<-metadiag(glas.t, re = "sm", re.model = "SeSp", link = "probit", split.w = TRUE,
#' df.estimate = TRUE))
#' plot(m, parametric.smooth = FALSE, level = c(0.5, 0.9))
#' plotw(m)
#'
#'# cloglog-sm-SeSp-df
#' summary(m<-metadiag(glas.t, re = "sm", re.model = "SeSp", link = "cloglog", split.w = TRUE,
#' df.estimate = TRUE))
#' plot(m, parametric.smooth = FALSE, level = c(0.5, 0.9))
#' plotw(m)
#'
#'# probit-sm-SeSp-df
#' summary(m<-metadiag(glas.t, re = "sm", re.model = "SeSp", link = "probit", split.w = TRUE,
#' df.estimate = TRUE))
#' plot(m, parametric.smooth = FALSE, level = c(0.5, 0.9))
#' plotw(m)
#'
#' }
#'
#' @import R2jags
#' @import rjags
#'
#' @export
#'
metadiag <- function(data,
                     two.by.two = FALSE,
                     # Arguments for the model:
                     re              = "normal",
                     re.model        = "DS",
                     link            = "logit",
                     # Hyperpriors parameters............................................
                     mean.mu.D       = 0,
                     mean.mu.S       = 0,
                     sd.mu.D         = 1,
                     sd.mu.S         = 1,
                     sigma.D.upper   = 10,
                     sigma.S.upper   = 10,
                     mean.Fisher.rho = 0,
                     sd.Fisher.rho   = 1/sqrt(2),
                     df              = 4,
                     df.estimate     = FALSE,
                     df.lower        = 3,
                     df.upper        = 20,

                     # Split weights
                     split.w         = FALSE,

                     # Predictions
                     n.1.new         = 50,
                     n.2.new         = 50,

                     # MCMC setup........................................................
                     nr.chains       = 2,
                     nr.iterations   = 10000,
                     nr.adapt        = 1000,
                     nr.burnin       = 1000,
                     nr.thin         = 1,

                     # Further options to link jags and R ...............................
                     be.quiet        = FALSE,
                     r2jags          = TRUE
                     )UseMethod("metadiag")

#' @export
metadiag.default <- function(
          # Data
           data,
           two.by.two = FALSE,
          # Arguments for the model:
          re              = "normal",

          # Parametrization...................................................
          re.model        = "SeSp",
          link            = "logit",

					# Hyperpriors parameters............................................
					mean.mu.D       = 0,
					mean.mu.S       = 0,
					sd.mu.D         = 1,
					sd.mu.S         = 1,
					sigma.D.upper   = 10,
					sigma.S.upper   = 10,
					mean.Fisher.rho = 0,
					sd.Fisher.rho   = 1/sqrt(2),
					df              = 4,

					df.estimate     = FALSE,
					df.lower        = 3,
					df.upper        = 20,

          # Split weights
          split.w         = FALSE,

          # Predictions
          n.1.new         = 50,
          n.2.new         = 50,

          # MCMC setup........................................................
					nr.chains       = 2,
					nr.iterations   = 10000,
					nr.adapt        = 1000,
					nr.burnin       = 1000,
          nr.thin         = 1,

          # Further options to link jags and R ...............................
					be.quiet        = FALSE,
          r2jags          = TRUE
          )
{


# Model errors checking-----

if(re=="normal" & split.w==TRUE)stop("Normal random effects and splitting weights are not compatible options")

re.test <- re %in% c("normal", "sm")

if(!re.test)stop("This random effects distribution is not implemented")

re.model.test <- re.model %in% c("SeSp", "DS")

if(!re.model.test)stop("This random effects parametrization is not implemented")

link.test <- link %in% c("logit", "cloglog", "probit")

if(!link.test)stop("This link function is not implemented")

	# Setting up hyperparameters ...
	      pre.mu.D <- 1/(sd.mu.D*sd.mu.D)
	      pre.mu.S <- 1/(sd.mu.S*sd.mu.S)
	pre.Fisher.rho <- 1/(sd.Fisher.rho * sd.Fisher.rho)


  # Setting up data nodes ...

	if(two.by.two == FALSE)
	{
	  tp <- data[,1]
	  n1 <- data[,2]
	  fp <- data[,3]
	  n2 <- data[,4]
	} else
	{
	  tp <- data$TP
	  fp <- data$FP
	  fn <- data$FN
	  tn <- data$TN
	  n1 <- tp + fn
	  n2 <- fp + tn
	}

	N <- length(n1)

  # Data errors
  if(any(tp>n1) || any(fp>n2))stop("the data is inconsistent")

  if(missing(data))stop("NAs are not alow in this function")

	# Statistical warning!
	if(N<5)warning("\n ************************************** \n
	                The number of studies is lower than the number of fixed parameters!\n
	                Posteriors have been estimated correctly, but you should check the
	                setup of the priors.")


#-----
# Data, initial values and parameters .....................................................
	data.model <-
            list(N = N,
	              tp = tp,
	              fp = fp,
	              n1 = n1,
	              n2 = n2,
	       mean.mu.D = mean.mu.D,
	       mean.mu.S = mean.mu.S,
	        pre.mu.D = pre.mu.D,
	        pre.mu.S = pre.mu.S,
	   sigma.D.upper = sigma.D.upper,
	   sigma.S.upper = sigma.S.upper,
	 mean.Fisher.rho = mean.Fisher.rho,
	  pre.Fisher.rho = pre.Fisher.rho,
	 # Predictions
	 n.1.new         = 50,
	 n.2.new         = 50)

	if(re == "sm" & df.estimate == TRUE){data.model$df.lower <- df.lower
	                                    data.model$df.upper <- df.upper}
	if(re == "sm" & df.estimate == FALSE){data.model$df <- df}

  # Parameters to monitor ....................................................................
	parameters.model.DS <-
    c("se.pool",
      "sp.pool",
       "se.new",
       "sp.new",
       "tp.new",
       "fp.new",
         "mu.D",
         "mu.S",
      "sigma.D",
      "sigma.S",
            "rho",
      "se",
      "sp")

	parameters.model.SeSp <-
	  c("se.pool",
	    "sp.pool",
	    "se.new",
	    "sp.new",
	    "tp.new",
	    "fp.new",
	    "mu.Se",
	    "mu.Sp",
	    "sigma.Se",
	    "sigma.Sp",
	    "rho",
	    "se",
	    "sp")

# Choose the list of parameters according to the parametrization ...
	if(re.model=="DS") parameters.model <- parameters.model.DS
	  else parameters.model <- parameters.model.SeSp

# This take the weights for the scale mixture random effects model, etc...

        if(re == "sm" & split.w == TRUE & df.estimate == TRUE)
          parameters.model <- c(parameters.model[], "w1", "w2", "p.w1", "p.w2", "df")
        else
          if(re == "sm" & split.w == TRUE & df.estimate == FALSE)
            parameters.model <- c(parameters.model[], "w1", "w2", "p.w1", "p.w2")
	      else
          if(re=="sm" & split.w == FALSE & df.estimate == TRUE)
            parameters.model <- c(parameters.model, "w", "p.w", "df")


# Model construction

inits.model <- list(mu.D = 0,
                    mu.S = 0,
                 sigma.D = 1,
                 sigma.S = 1,
                 z       = 0)

# Model BUGS script

#----
blueprint <- function(link = "logit", re = "normal", re.model = "DS", split.w = FALSE, df.estimate = FALSE)
{

   if(split.w == FALSE & df.estimate == TRUE ) re <- "sm.df"
   if(split.w == TRUE  & df.estimate == FALSE) re <- "sm.split"
   if(split.w == TRUE  & df.estimate == TRUE ) re <- "sm.split.df"

  #----
  # Block for data model ......................................................................
  dm <-
    "
  model
{
  for(i in 1:N)
{
  tp[i] ~ dbin(TPR[i], n1[i])
  fp[i] ~ dbin(FPR[i], n2[i])

  se[i] = TPR[i]
  sp[i] = 1 - FPR[i]

  "

  #----


  #----
  # Block for the link function................................................................
  link.logit.DS <-
    "
  logit(TPR[i]) <- (D[i] + S[i])/2
  logit(FPR[i]) <- (S[i] - D[i])/2
  "

  link.cloglog.DS <-
    "
  cloglog(TPR[i]) <- (D[i] + S[i])/2
  cloglog(FPR[i]) <- (S[i] - D[i])/2
  "

  link.probit.DS <-
    "
  probit(TPR[i]) <- (D[i] + S[i])/2
  probit(FPR[i]) <- (S[i] - D[i])/2
  "

  link.logit.SeSp <-
    "
  logit(TPR[i]) <-    l.se[i]
  logit(FPR[i]) <- -1*l.sp[i]
  "


  link.cloglog.SeSp <-
    "
  cloglog(TPR[i]) <- l.se[i]
  cloglog(FPR[i]) <- -1*l.sp[i]
  "

  link.probit.SeSp <-
    "
  probit(TPR[i]) <- l.se[i]
  probit(FPR[i]) <- -1*l.sp[i]
  "

  # Choose the link according to parametrization ...
  if(re.model == "DS"){  link.logit <- link.logit.DS
                       link.cloglog <- link.cloglog.DS
                        link.probit <- link.probit.DS}
     else {
       link.logit <- link.logit.SeSp
     link.cloglog <- link.cloglog.SeSp
      link.probit <- link.probit.SeSp}


  #----

  # Block for structural distribution .........................................................
  re.normal.DS <-
    "
  S[i] ~ dnorm(mu.S, pre.S)
  D[i] ~ dnorm(mu.D.S[i], pre.D.S)
  mu.D.S[i] <- mu.D + rho * sigma.D / sigma.S * (S[i] - mu.S)

}

  # Hyper priors
  mu.D ~  dlogis(mean.mu.D, pre.mu.D)
  mu.S ~  dlogis(mean.mu.S, pre.mu.S)

  # Dispersion parameters
  sigma.D ~ dunif(0, sigma.D.upper)
  sigma.S ~ dunif(0, sigma.S.upper)
  pre.D <- 1/(sigma.D*sigma.D)
  pre.S <- 1/(sigma.S*sigma.S)

  #Conditional precision
  pre.D.S <- pre.D/(1-rho*rho)

  # Correlation
  z ~ dnorm(mean.Fisher.rho, pre.Fisher.rho)
  rho <- 2*exp(z)/(1+exp(z)) - 1

  # Predictions ...
  mu.hat[1] <- mu.D
  mu.hat[2] <- mu.S

  Sigma.hat[1, 1] <- pow(sigma.D, 2)
  Sigma.hat[2, 2] <- pow(sigma.S, 2)
  Sigma.hat[1, 2] <- rho * sigma.D * sigma.S
  Sigma.hat[2, 1] <- Sigma.hat[1, 2]

  Omega.hat[1:2, 1:2] <- inverse(Sigma.hat[1:2, 1:2])
  DSnew[1:2] ~ dmnorm(mu.hat[1:2], Omega.hat[1:2 ,1:2])
  "

  re.normal.SeSp <-
    "
  l.sp[i] ~ dnorm(mu.Sp, pre.Sp)
  l.se[i] ~ dnorm(mu.Se.Sp[i], pre.Se.Sp)
  mu.Se.Sp[i] <- mu.Se + rho * sigma.Se / sigma.Sp * (l.sp[i] - mu.Sp)
}

  # Hyper priors
  mu.Se ~  dlogis(mean.mu.D, pre.mu.D)
  mu.Sp ~  dlogis(mean.mu.S, pre.mu.S)

  # Dispersion parameters
  sigma.Se ~ dunif(0, sigma.D.upper)
  sigma.Sp ~ dunif(0, sigma.S.upper)
  pre.Se <- 1/(sigma.Se*sigma.Se)
  pre.Sp <- 1/(sigma.Sp*sigma.Sp)

  #Conditional precision
  pre.Se.Sp <- pre.Se/(1-rho*rho)

  # Correlation
  z ~ dnorm(mean.Fisher.rho, pre.Fisher.rho)
  rho <- 2*exp(z)/(1+exp(z)) - 1

  # Predictions ...
  mu.hat[1] <- mu.Se
  mu.hat[2] <- mu.Sp

  Sigma.hat[1, 1] <- pow(sigma.Se, 2)
  Sigma.hat[2, 2] <- pow(sigma.Sp, 2)
  Sigma.hat[1, 2] <- rho * sigma.Se * sigma.Sp
  Sigma.hat[2, 1] <- Sigma.hat[1, 2]

  Omega.hat[1:2, 1:2] <- inverse(Sigma.hat[1:2, 1:2])
  l.sesp.new[1:2] ~ dmnorm(mu.hat[1:2], Omega.hat[1:2 ,1:2])
  "

#Choose the model according to parametrization ...

    if(re.model=="DS"){re.normal <- re.normal.DS}
     else {re.normal <- re.normal.SeSp}


  #----
  re.sm.DS <-
    "
  S[i] ~ dnorm(mu.S, pre.w.S[i])
  D[i] ~ dnorm(mu.D.S[i], pre.w.D.S[i])
  mu.D.S[i] <- mu.D + rho * sigma.D / sigma.S * (S[i] - mu.S)

  lambda[i] ~ dchisqr(df)
  w[i] <- df/lambda[i]

  pre.w.D[i] <- pre.D / w[i]
  pre.w.S[i] <- pre.S / w[i]

  #Conditional precision
  pre.w.D.S[i] <- pre.w.D[i] / (1-rho*rho)
}


  # Hyper priors
  mu.D ~  dlogis(mean.mu.D, pre.mu.D)
  mu.S ~  dlogis(mean.mu.S, pre.mu.S)

  # Dispersion parameters
  sigma.D ~ dunif(0, sigma.D.upper)
  sigma.S ~ dunif(0, sigma.S.upper)
  pre.D <- 1/(sigma.D*sigma.D)
  pre.S <- 1/(sigma.S*sigma.S)


  # Correlation
  z ~ dnorm(mean.Fisher.rho, pre.Fisher.rho)
  rho <- 2*exp(z)/(1+exp(z)) - 1

  # Predictions ...
  mu.hat[1] <- mu.D
  mu.hat[2] <- mu.S

  Sigma.hat[1, 1] <- pow(sigma.D, 2)
  Sigma.hat[2, 2] <- pow(sigma.S, 2)
  Sigma.hat[1, 2] <- rho * sigma.D * sigma.S
  Sigma.hat[2, 1] <- Sigma.hat[1, 2]

  lambda.new ~ dchisqr(df)
  w.new <- df/lambda.new
  Omega.hat[1:2,1:2] <- inverse(Sigma.hat[1:2, 1:2]) / w.new

  DSnew[1:2] ~ dmnorm(mu.hat[1:2], Omega.hat[1:2 ,1:2])

  # Probability of outlier
  for( i in 1:N)
  {
  p.w[i] <- step(w[i]-0.99)
  }

  "

#----
re.sm.SeSp <-
  "
  l.sp[i] ~ dnorm(mu.Sp, pre.w.Sp[i])
  l.se[i] ~ dnorm(mu.Se.Sp[i], pre.w.Se.Sp[i])
mu.Se.Sp[i] <- mu.Se + rho * sigma.Se / sigma.Sp * (l.sp[i] - mu.Sp)

lambda[i] ~ dchisqr(df)
w[i] <- df/lambda[i]

pre.w.Se[i] <- pre.Se / w[i]
pre.w.Sp[i] <- pre.Sp / w[i]

#Conditional precision

pre.w.Se.Sp[i] <- pre.w.Se[i] / (1-rho*rho)
}


# Hyper priors
  mu.Se ~ dlogis(mean.mu.D, pre.mu.D)
  mu.Sp ~ dlogis(mean.mu.S, pre.mu.S)

# Dispersion parameters
sigma.Se ~ dunif(0, sigma.D.upper)
sigma.Sp ~ dunif(0, sigma.S.upper)
pre.Se <- 1/(sigma.Se*sigma.Se)
pre.Sp <- 1/(sigma.Sp*sigma.Sp)

# Correlation
z ~ dnorm(mean.Fisher.rho, pre.Fisher.rho)
rho <- 2*exp(z)/(1+exp(z)) - 1

# Predictions ...
mu.hat[1] <- mu.Se
mu.hat[2] <- mu.Sp

Sigma.hat[1, 1] <- pow(sigma.Se, 2)
Sigma.hat[2, 2] <- pow(sigma.Sp, 2)
Sigma.hat[1, 2] <- rho * sigma.Se * sigma.Sp
Sigma.hat[2, 1] <- Sigma.hat[1, 2]

lambda.new ~ dchisqr(df)
w.new <- df/lambda.new
Omega.hat[1:2,1:2] <- inverse(Sigma.hat[1:2, 1:2]) / w.new

l.sesp.new[1:2] ~ dmnorm(mu.hat[1:2], Omega.hat[1:2 ,1:2])

# Probability of outlier
for( i in 1:N)
{
  p.w[i] <- step(w[i]-0.99)
}

"


#Choose the model according to parametrization ...

if(re.model=="DS"){re.sm <- re.sm.DS}
   else {re.sm <- re.sm.SeSp}

#--------------------------------------------------

  re.sm.split.DS <-
    "
  S[i] ~ dnorm(mu.S, pre.w.S[i])
  D[i] ~ dnorm(mu.D.S[i], pre.w.D.S[i])
  mu.D.S[i] <- mu.D + rho * sigma.D / sigma.S * (S[i] - mu.S)

  lambda1[i] ~ dchisqr(df)
  w1[i] <- df/lambda1[i]

  lambda2[i] ~ dchisqr(df)
  w2[i] <- df/lambda2[i]

  pre.w.D[i] <- pre.D / w1[i]
  pre.w.S[i] <- pre.S / w2[i]

  #Conditional precision
  pre.w.D.S[i] <- pre.w.D[i] / (1-rho*rho)
}

  # Hyper priors
  mu.D ~  dnorm(mean.mu.D, pre.mu.D)
  mu.S ~  dnorm(mean.mu.S, pre.mu.S)

  # Dispersion parameters
  sigma.D ~ dunif(0, sigma.D.upper)
  sigma.S ~ dunif(0, sigma.S.upper)
  pre.D <- 1/(sigma.D*sigma.D)
  pre.S <- 1/(sigma.S*sigma.S)


  # Correlation
  z ~ dnorm(mean.Fisher.rho, pre.Fisher.rho)
  rho <- 2*exp(z)/(1+exp(z)) - 1

  # Predictions ............................................................
  mu.hat[1] <- mu.D
  mu.hat[2] <- mu.S

  lambda1.new ~ dchisqr(df)
  lambda2.new ~ dchisqr(df)

  w1.new <- df/lambda1.new
  w2.new <- df/lambda2.new

  Sigma.hat[1, 1] <- pow(sigma.D, 2)/w1.new
  Sigma.hat[2, 2] <- pow(sigma.S, 2)/w2.new
  Sigma.hat[1, 2] <- rho * sigma.D * sigma.S / sqrt(w1.new*w2.new)
  Sigma.hat[2, 1] <- Sigma.hat[1, 2]

  Omega.hat[1:2,1:2] <- inverse(Sigma.hat[1:2, 1:2])

  DSnew[1:2] ~ dmnorm(mu.hat[1:2], Omega.hat[1:2 ,1:2])


  # Probability of outlier
  for( i in 1:N)
  {
  p.w1[i] <- step(w1[i]-0.99)
  p.w2[i] <- step(w2[i]-0.99)
  }
  "


re.sm.split.SeSp <-
  "
l.sp[i] ~ dnorm(mu.Sp, pre.w.Sp[i])
l.se[i] ~ dnorm(mu.Se.Sp[i], pre.w.Se.Sp[i])
mu.Se.Sp[i] <- mu.Se + rho * sigma.Se / sigma.Sp * (l.sp[i] - mu.Sp)

lambda1[i] ~ dchisqr(df)
w1[i] <- df/lambda1[i]

lambda2[i] ~ dchisqr(df)
w2[i] <- df/lambda2[i]

pre.w.Se[i] <- pre.Se / w1[i]
pre.w.Sp[i] <- pre.Sp / w2[i]

#Conditional precision
pre.w.Se.Sp[i] <- pre.w.Se[i] / (1-rho*rho)
}


# Hyper priors
mu.Se ~ dlogis(mean.mu.D, pre.mu.D)
mu.Sp ~ dlogis(mean.mu.S, pre.mu.S)

# Dispersion parameters
sigma.Se ~ dunif(0, sigma.D.upper)
sigma.Sp ~ dunif(0, sigma.S.upper)
pre.Se <- 1/(sigma.Se*sigma.Se)
pre.Sp <- 1/(sigma.Sp*sigma.Sp)

# Correlation
z ~ dnorm(mean.Fisher.rho, pre.Fisher.rho)
rho <- 2*exp(z)/(1+exp(z)) - 1

# Predictions ............................................................

mu.hat[1] <- mu.Se
mu.hat[2] <- mu.Sp

  lambda1.new ~ dchisqr(df)
  lambda2.new ~ dchisqr(df)

w1.new <- df/lambda1.new
w2.new <- df/lambda2.new

  Sigma.hat[1, 1] <- pow(sigma.Se, 2)/w1.new
  Sigma.hat[2, 2] <- pow(sigma.Sp, 2)/w2.new
  Sigma.hat[1, 2] <- rho * sigma.Se * sigma.Sp / sqrt(w1.new*w2.new)
  Sigma.hat[2, 1] <- Sigma.hat[1, 2]

Omega.hat[1:2,1:2] <- inverse(Sigma.hat[1:2, 1:2])

l.sesp.new[1:2] ~ dmnorm(mu.hat[1:2], Omega.hat[1:2 ,1:2])

# Probability of outlier
  for( i in 1:N)
{
  p.w1[i] <- step(w1[i]-0.99)
  p.w2[i] <- step(w2[i]-0.99)
}
"


#Choose the model according to parametrization ...

if(re.model=="DS"){re.sm.split <- re.sm.split.DS}
else {re.sm.split <- re.sm.split.SeSp}


#------------------------------------------------------------------------

re.sm.df.DS <-
  "
S[i] ~ dnorm(mu.S, pre.w.S[i])
D[i] ~ dnorm(mu.D.S[i], pre.w.D.S[i])
mu.D.S[i] <- mu.D + rho * sigma.D / sigma.S * (S[i] - mu.S)

lambda[i] ~ dchisqr(df)
w[i] <- df/lambda[i]

pre.w.D[i] <- pre.D / w[i]
pre.w.S[i] <- pre.S / w[i]

#Conditional precision
pre.w.D.S[i] <- pre.w.D[i] / (1-rho*rho)
}


# Hyper priors
mu.D ~  dlogis(mean.mu.D, pre.mu.D)
mu.S ~  dlogis(mean.mu.S, pre.mu.S)

# Dispersion parameters
sigma.D ~ dunif(0, sigma.D.upper)
sigma.S ~ dunif(0, sigma.S.upper)
pre.D <- 1/(sigma.D*sigma.D)
pre.S <- 1/(sigma.S*sigma.S)

# Correlation
z ~ dnorm(mean.Fisher.rho, pre.Fisher.rho)
rho <- 2*exp(z)/(1+exp(z)) - 1

# Degrees of freedom
  a <- 1/df.upper
  b <- 1/df.lower
 df <- 1/U
  U ~ dunif(a, b)

# Predictions ...
mu.hat[1] <- mu.D
mu.hat[2] <- mu.S

Sigma.hat[1, 1] <- pow(sigma.D, 2)
Sigma.hat[2, 2] <- pow(sigma.S, 2)
Sigma.hat[1, 2] <- rho * sigma.D * sigma.S
Sigma.hat[2, 1] <- Sigma.hat[1, 2]

lambda.new ~ dchisqr(df)
w.new <- df/lambda.new
Omega.hat[1:2,1:2] <- inverse(Sigma.hat[1:2, 1:2]) / w.new

DSnew[1:2] ~ dmnorm(mu.hat[1:2], Omega.hat[1:2 ,1:2])

# Probability of outlier
for( i in 1:N)
{
  p.w[i] <- step(w[i]-0.99)
}

"



re.sm.df.SeSp <-
  "
  l.sp[i] ~ dnorm(mu.Sp, pre.w.Sp[i])
l.se[i] ~ dnorm(mu.Se.Sp[i], pre.w.Se.Sp[i])
mu.Se.Sp[i] <- mu.Se + rho * sigma.Se / sigma.Sp * (l.sp[i] - mu.Sp)

lambda[i] ~ dchisqr(df)
w[i] <- df/lambda[i]

pre.w.Se[i] <- pre.Se / w[i]
pre.w.Sp[i] <- pre.Sp / w[i]

#Conditional precision
pre.w.Se.Sp[i] <- pre.w.Se[i] / (1-rho*rho)
}


# Hyper priors
mu.Se ~ dlogis(mean.mu.D, pre.mu.D)
mu.Sp ~ dlogis(mean.mu.S, pre.mu.S)

# Dispersion parameters
sigma.Se ~ dunif(0, sigma.D.upper)
sigma.Sp ~ dunif(0, sigma.S.upper)
pre.Se <- 1/(sigma.Se*sigma.Se)
pre.Sp <- 1/(sigma.Sp*sigma.Sp)

# Correlation
z ~ dnorm(mean.Fisher.rho, pre.Fisher.rho)
rho <- 2*exp(z)/(1+exp(z)) - 1


# Degrees of freedom
a <- 1/df.upper
b <- 1/df.lower
df <- 1/U
U ~ dunif(a, b)


# Predictions ...
mu.hat[1] <- mu.Se
mu.hat[2] <- mu.Sp

Sigma.hat[1, 1] <- pow(sigma.Se, 2)
Sigma.hat[2, 2] <- pow(sigma.Sp, 2)
Sigma.hat[1, 2] <- rho * sigma.Se * sigma.Sp
Sigma.hat[2, 1] <- Sigma.hat[1, 2]

lambda.new ~ dchisqr(df)
w.new <- df/lambda.new
Omega.hat[1:2,1:2] <- inverse(Sigma.hat[1:2, 1:2]) / w.new

l.sesp.new[1:2] ~ dmnorm(mu.hat[1:2], Omega.hat[1:2 ,1:2])

# Probability of outlier
for( i in 1:N)
{
  p.w[i] <- step(w[i]-0.99)
}
"

#Choose the model according to parametrization ...
if(re.model=="DS"){re.sm.df <- re.sm.df.DS}
else {re.sm.df <- re.sm.df.SeSp}

#------------------------------------------------------------------------

re.sm.split.df.DS <-
  "
S[i] ~ dnorm(mu.S, pre.w.S[i])
D[i] ~ dnorm(mu.D.S[i], pre.w.D.S[i])
mu.D.S[i] <- mu.D + rho * sigma.D / sigma.S * (S[i] - mu.S)

lambda1[i] ~ dchisqr(df)
w1[i] <- df/lambda1[i]

lambda2[i] ~ dchisqr(df)
w2[i] <- df/lambda2[i]

pre.w.D[i] <- pre.D / w1[i]
pre.w.S[i] <- pre.S / w2[i]

#Conditional precision
pre.w.D.S[i] <- pre.w.D[i] / (1-rho*rho)
}

# Hyper priors
mu.D ~  dnorm(mean.mu.D, pre.mu.D)
mu.S ~  dnorm(mean.mu.S, pre.mu.S)

# Dispersion parameters
sigma.D ~ dunif(0, sigma.D.upper)
sigma.S ~ dunif(0, sigma.S.upper)
pre.D <- 1/(sigma.D*sigma.D)
pre.S <- 1/(sigma.S*sigma.S)


# Correlation
z ~ dnorm(mean.Fisher.rho, pre.Fisher.rho)
rho <- 2*exp(z)/(1+exp(z)) - 1

# Degrees of freedom
  a <- 1/df.upper
b <- 1/df.lower
df <- 1/U
U ~ dunif(a, b)

# Predictions ............................................................

mu.hat[1] <- mu.D
mu.hat[2] <- mu.S

lambda1.new ~ dchisqr(df)
lambda2.new ~ dchisqr(df)

w1.new <- df/lambda1.new
w2.new <- df/lambda2.new

Sigma.hat[1, 1] <- pow(sigma.D, 2)/w1.new
Sigma.hat[2, 2] <- pow(sigma.S, 2)/w2.new
Sigma.hat[1, 2] <- rho * sigma.D * sigma.S / sqrt(w1.new*w2.new)
Sigma.hat[2, 1] <- Sigma.hat[1, 2]

Omega.hat[1:2,1:2] <- inverse(Sigma.hat[1:2, 1:2])

DSnew[1:2] ~ dmnorm(mu.hat[1:2], Omega.hat[1:2 ,1:2])


# Probability of outlier

for( i in 1:N)
{
  p.w1[i] <- step(w1[i]-0.99)
  p.w2[i] <- step(w2[i]-0.99)
}
"


re.sm.split.df.SeSp <-
  "
l.sp[i] ~ dnorm(mu.Sp, pre.w.Sp[i])
l.se[i] ~ dnorm(mu.Se.Sp[i], pre.w.Se.Sp[i])
mu.Se.Sp[i] <- mu.Se + rho * sigma.Se / sigma.Sp * (l.sp[i] - mu.Sp)

lambda1[i] ~ dchisqr(df)
w1[i] <- df/lambda1[i]

lambda2[i] ~ dchisqr(df)
w2[i] <- df/lambda2[i]

pre.w.Se[i] <- pre.Se / w1[i]
pre.w.Sp[i] <- pre.Sp / w2[i]

#Conditional precision
pre.w.Se.Sp[i] <- pre.w.Se[i] / (1-rho*rho)
}


# Hyper priors
mu.Se ~ dlogis(mean.mu.D, pre.mu.D)
mu.Sp ~ dlogis(mean.mu.S, pre.mu.S)

# Dispersion parameters
sigma.Se ~ dunif(0, sigma.D.upper)
sigma.Sp ~ dunif(0, sigma.S.upper)
pre.Se <- 1/(sigma.Se*sigma.Se)
pre.Sp <- 1/(sigma.Sp*sigma.Sp)

# Correlation
z ~ dnorm(mean.Fisher.rho, pre.Fisher.rho)
rho <- 2*exp(z)/(1+exp(z)) - 1

# Degrees of freedom
  a <- 1/df.upper
  b <- 1/df.lower
 df <- 1/U
   U ~ dunif(a, b)


# Predictions ............................................................

mu.hat[1] <- mu.Se
mu.hat[2] <- mu.Sp

  lambda1.new ~ dchisqr(df)
  lambda2.new ~ dchisqr(df)

  w1.new <- df/lambda1.new
  w2.new <- df/lambda2.new

Sigma.hat[1, 1] <- pow(sigma.Se, 2)/w1.new
Sigma.hat[2, 2] <- pow(sigma.Sp, 2)/w2.new
Sigma.hat[1, 2] <- rho * sigma.Se * sigma.Sp / sqrt(w1.new*w2.new)
Sigma.hat[2, 1] <- Sigma.hat[1, 2]

Omega.hat[1:2,1:2] <- inverse(Sigma.hat[1:2, 1:2])

l.sesp.new[1:2] ~ dmnorm(mu.hat[1:2], Omega.hat[1:2 ,1:2])

# Probability of outlier
for( i in 1:N)
{
  p.w1[i] <- step(w1[i]-0.99)
  p.w2[i] <- step(w2[i]-0.99)
}
"

#Choose the model according to parametrization ...
if(re.model=="DS"){re.sm.split.df <- re.sm.split.df.DS}
   else {re.sm.split.df <- re.sm.split.df.SeSp}


  # Block of parameters of interest depending on the links .................................
  par.logit.DS <- "
  # Parameters of interest
  # Pooled summaries ...
  x <- (mu.D + mu.S)/2
  y <- (mu.S - mu.D)/2
  se.pool <- ilogit(x)
  sp.pool <- 1 - ilogit(y)

  # Predictive summaries ...
  x.new <- (DSnew[1] + DSnew[2])/2
  y.new <- (DSnew[2] - DSnew[1])/2
  se.new <- ilogit(x.new)
  fpr.new <- ilogit(y.new)
  sp.new <- 1 - fpr.new

  tp.new ~ dbin( se.new, n.1.new)
  fp.new ~ dbin(fpr.new, n.2.new)
}"

  par.logit.SeSp <- "
  # Parameters of interest
  # Pooled summaries ...
  se.pool <- ilogit(mu.Se)
  sp.pool <- ilogit(mu.Sp)

  # Predictive summaries ...
  se.new <- ilogit(l.sesp.new[1])
  sp.new <- ilogit(l.sesp.new[2])
  fpr.new <- 1- sp.new
  tp.new ~ dbin( se.new, n.1.new)
  fp.new ~ dbin(fpr.new, n.2.new)
}"


  #Choose the link according to parametrization ...
  if(re.model=="DS"){par.logit <- par.logit.DS}
  else {par.logit <- par.logit.SeSp}


  par.cloglog.DS <- "
  # Parameters of interest
  # Pooled summaries ...
  x <- (mu.D + mu.S)/2
  y <- (mu.S - mu.D)/2
  se.pool <- icloglog(x)
  sp.pool <- 1 - icloglog(y)

  # Predictive summaries ...
  x.new <- (DSnew[1] + DSnew[2])/2
  y.new <- (DSnew[2] - DSnew[1])/2
  se.new <- icloglog(x.new)
  fpr.new <- icloglog(y.new)
  sp.new <- 1 - fpr.new

  tp.new ~ dbin( se.new, n.1.new)
  fp.new ~ dbin(fpr.new, n.2.new)
}"


  par.cloglog.SeSp <- "
 # Parameters of interest
 # Pooled summaries ...
  se.pool <- icloglog(mu.Se)
  sp.pool <- icloglog(mu.Sp)

  # Predictive summaries ...
  se.new <- icloglog(l.sesp.new[1])
  sp.new <- icloglog(l.sesp.new[2])
  fpr.new <- 1- sp.new
  tp.new ~ dbin( se.new, n.1.new)
  fp.new ~ dbin(fpr.new, n.2.new)
  }
  "

  #Choose the link according to parametrization ...
  if(re.model=="DS"){par.cloglog <- par.cloglog.DS}
     else {par.cloglog <- par.cloglog.SeSp}


  par.probit.DS <-
  "
  # Parameters of interest
  # Pooled summaries ...
  x <- (mu.D + mu.S)/2
  y <- (mu.S - mu.D)/2
  se.pool <- phi(x)
  sp.pool <- 1 - phi(y)

  # Predictive summaries ...
  x.new <- (DSnew[1] + DSnew[2])/2
  y.new <- (DSnew[2] - DSnew[1])/2
  se.new <- phi(x.new)
  fpr.new <- phi(y.new)
  sp.new <- 1 - fpr.new

  tp.new ~ dbin(se.new, n.1.new)
  fp.new ~ dbin(fpr.new, n.2.new)
  }
  "

  par.probit.SeSp <- "
 # Parameters of interest
 # Pooled summaries ...
  se.pool <- phi(mu.Se)
  sp.pool <- phi(mu.Sp)

  # Predictive summaries ...
  se.new <- phi(l.sesp.new[1])
  sp.new <- phi(l.sesp.new[2])
  fpr.new <- 1- sp.new
  tp.new ~ dbin( se.new, n.1.new)
  fp.new ~ dbin(fpr.new, n.2.new)
}
  "

  #Choose the link according to parametrization ...
  if(re.model=="DS"){par.probit <- par.probit.DS}
  else {par.probit <- par.probit.SeSp}



#----
# possible models
#normal random effects
m1 <- paste(dm, link.logit,   re.normal, par.logit)
m2 <- paste(dm, link.cloglog, re.normal, par.cloglog)
m3 <- paste(dm, link.probit,  re.normal, par.probit)

# sm random effects
m4 <- paste(dm, link.logit,   re.sm, par.logit)
m5 <- paste(dm, link.cloglog, re.sm, par.cloglog)
m6 <- paste(dm, link.probit,  re.sm, par.probit)

# sm random effects with two t-distributions one for each random effect
m7 <- paste(dm, link.logit,   re.sm.split,  par.logit)
m8 <- paste(dm, link.cloglog, re.sm.split,  par.cloglog)
m9 <- paste(dm, link.probit,  re.sm.split,  par.probit)

# sm random effects
m10 <- paste(dm, link.logit,   re.sm.df, par.logit)
m11 <- paste(dm, link.cloglog, re.sm.df, par.cloglog)
m12 <- paste(dm, link.probit,  re.sm.df, par.probit)

# sm random effects with two t-distributions one for each random effect
m13 <- paste(dm, link.logit,   re.sm.split.df,  par.logit)
m14 <- paste(dm, link.cloglog, re.sm.split.df,  par.cloglog)
m15 <- paste(dm, link.probit,  re.sm.split.df,  par.probit)


#----
  switch(re,
         normal = switch(link,
                           logit = return(m1),
                         cloglog = return(m2),
                          probit = return(m3),
                         stop("The model you requested is not implemented.")
                         ),
             sm = switch(link,
                           logit = return(m4),
                         cloglog = return(m5),
                          probit = return(m6),
                         stop("The model you requested is not implemented.")
                         ),
         sm.split = switch(link,
                           logit = return(m7),
                         cloglog = return(m8),
                          probit = return(m9),
                         stop("The model you requested is not implemented.")
                         ),
           sm.df = switch(link,
                           logit = return(m10),
                         cloglog = return(m11),
                          probit = return(m12),
                         stop("The model you requested is not implemented.")
         ),
         sm.split.df = switch(link,
                           logit = return(m13),
                         cloglog = return(m14),
                          probit = return(m15),
                         stop("The model you requested is not implemented.")
         ),
        stop("The model you requested is not implemented.")
        )
}

model.bugs <- blueprint(link, re, re.model, split.w, df.estimate)

if(r2jags == TRUE){
  # Use R2jags as interface for JAGS ...
  results <- jags(              data = data.model,
                  parameters.to.save = parameters.model,
                            #   inits = inits.model,
                          model.file = textConnection(model.bugs),
                            n.chains = nr.chains,
                              n.iter = nr.iterations,
                            n.burnin = nr.burnin,
                              n.thin = nr.thin
                       )
  }
  else {
  # Use rjags as interface for JAGS ...
  # Send the model to JAGS, check syntax, run ...
	jm <- jags.model(file     = textConnection(model.bugs),
	                 data     = data.model,
                   inits    = inits.model,
	                 n.chains = nr.chains,
	                 n.adapt  = nr.adapt,
	                 quiet    = be.quiet)

	results <- coda.samples(jm,
	                        variable.names = parameters.model,
	                        n.iter         = nr.iterations)
  }

if(r2jags == FALSE)
  {cat("You are using the package rjags as interface to JAGS.", "\n")
   cat("The plot functions for output analysis are not implemented in this bamdit version", "\n")
}

# Extra outputs that are linked with other functions

results$link <- link
results$re <- re
results$re.model <- re.model
results$data <- data
results$two.by.two <- two.by.two
#results$r2jags <- r2jags
results$split.w <- split.w
#results$df.estimate <- df.estimate

class(results) <- c("metadiag")

return(results)
}










