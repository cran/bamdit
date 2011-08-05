metadiag <- function(	data,
# arguments for the model:
          re="normal",
          link="logit",
# arguments for the priors:
					R = matrix(c(1, 0, 0, 1), byrow=TRUE, nrow=2), 
					m.0 = c(0,0),
					pre.mu = c(0.25, 0.25),
					k = 3,
					nu.0 = 1,
# arguments for the MCMC:			
					n.chains = 3, 
					n.iter = 40000, 
					n.burnin = 20000, 
					n.thin = 10,
					verbose = TRUE)
{

	# test the data an throw an error otherwise
	# the dimension is returned
	n <- .test.data(data)
	
	# get the data together
	tp <- data[,1]
	n1 <- data[,2]
	fp <- data[,3]
	n2 <- data[,4]
	
	# setup data vector for jags
	#data <- list("R", "tp", "n1", "fp", "n2", "n", "m.0", "pre.mu", "k", "nu.0")
	data <- list(tp=tp, n1=n1, fp=fp, n2=n2, n=n, k=k, R=R, m.0=m.0, nu.0=nu.0,
  pre.mu=pre.mu)
	
	# say the user that we writing the model
	if (verbose) .println("Writing model...")
	
	# write the bugs model:
	.write.bugs.model(re,link)
	if(!.check.params(re,link))
    .throw("The model you requested (",re,", ",link,") was not implemented yet.")
	
	# say the user that we are checking the model
	if (verbose) .println("Compiling model...")
	
	# Check the Model ...
	m <- jags.model(.in.tmp("bamdit-model.bug"),
	                 data = data,
	                 n.chains = 3)
	
	# Nodes to monitor ...
	modelPar <- c("pool.se", "pool.sp",
	               "new.se", "new.sp",
	               "mu", "sigmaD", "sigmaS",
	               "rhoDS")

	# Say that we are doing the burn-in
	if (verbose) .println(paste("Running ",n.burnin," MCMC iterations to get convergency...",sep=""))

	# Run MCMC itereations to get convergence...
	x <- jags.samples(m, modelPar, n.iter = n.burnin, thin=n.thin)

	# Say something smart
	if (verbose) .println("Obtaining statistics of ",n.iter - n.burnin," samples of monitored nodes...")

	# Update the MCMC iterations and coerces the output to a single mcmc.list object 
	x.out <- coda.samples(m, modelPar, n.iter = (n.iter - n.burnin), thin=(n.thin/2))
	
	# clean up											
	rm(data,R,n,tp,n1,fp,n2,m,modelPar,x)
	
	return(x.out)
}