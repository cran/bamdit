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
#	if (verbose) .println("Generating model blueprint ...")
	


	# It's done
#	if (verbose) .println("Done.")
	
	# Fit the Model the Model ...
	# ... or die!
	# set tries
	tries <- 3
	for (i in 1:tries){
		# say the user that we are checking the model
		if (verbose) .println("Trying to fit model in ",n.burnin," iterations ...")
		# try to make the model
		# get the handler
		blueprint = .model.handler(re,link)
		model <- jags.model(blueprint, n.adapt = n.burnin,
	                 	data = data,
	                 	n.chains = 3)
		close(blueprint)
		# if the model doesn't complete adaption
		if (.Call("is_adapting",model$ptr(),PACKAGE="rjags")){
			# if we are out of tries
			if (i < tries){
				if (verbose) .println("The model failed to adapt. I'll try to fix this myself, please stay put ...")
			} else {
				# bring this unnessery tradgedy to its end
				.throw("The model failed to adapt three times. I'm not able to fix this myself. Please try to change hyper-parameters values or to set n.iter and n.burnin parameters to higher values and retry yourself.\n")
			}
		} else break
	}
	
	# Nodes to monitor ...
	
	modelPar <- c("pool.se", "pool.sp",
	               "new.se", "new.sp",
	               "mu", "sigmaD", "sigmaS",
	               "rhoDS")
 
 # This take the weights too as parameters to return.
  if(re=="scalemix") modelPar <- c(modelPar, "w")
  
  
	# Say something smart
	if (verbose) .println("Completed fitting. Obtaining statistics of ",n.iter - n.burnin," samples of monitored nodes...")

	# Update the MCMC iterations and coerces the output to a single mcmc.list object 
	x.out <- coda.samples(model, modelPar, n.iter = (n.iter - n.burnin), thin=n.thin)
	
	# clean up											
	rm(data,R,n,tp,n1,fp,n2,model,modelPar)
	
	return(x.out)
}