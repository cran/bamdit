# Gives back a textConnection to the model-blueprint

.model.handler <- function(re="normal",link="logit"){
	
	normal.pre.link =		"model
		{
			for( i in 1 : n ) {
				tp[i] ~ dbin(tpr[i], n1[i])
				fp[i] ~ dbin(fpr[i], n2[i])
				"

	normal.post.link = "			    
						m[i,1:2] ~ dmnorm(mu[], sigma.inv[1:2 ,1:2] )
				  }
					# Priors ...
					mu[1] ~ dnorm(m.0[1], pre.mu[1])
					mu[2] ~ dnorm(m.0[2], pre.mu[2])
					sigma.inv[1:2,1:2] ~ dwish(R[1:2,1:2], k)

					# Summary statistics ...
					# Pooled summaries
					x <- (mu[1]+mu[2])/2
					y <- (mu[2]-mu[1])/2
					pool.se <- exp(x) / ( 1 + exp(x) )      # with logit link
					pool.sp <- 1 - exp(y) / ( 1 + exp(y) )  # with logit link

					# Predictive summaries ...
					m.star[1:2] ~ dmnorm(mu[], sigma.inv[1:2 ,1:2] )
					x.star <- (m.star[1]+m.star[2])/2
					y.star <- (m.star[2]-m.star[1])/2
					new.se <- exp(x.star)/(1 + exp(x.star))       # with logit link
					new.sp <- 1 - exp(y.star) /(1 + exp(y.star))  # with logit link

					# Variance covariance matrix for random effects
					sigma[1:2, 1:2] <- inverse(sigma.inv[ , ])


					# Variance covariance matrix for random effects
					sigmaD <- sigma[1,1]
					sigmaS <- sigma[2, 2]
				 	rhoDS  <- sigma[1,2]/(pow(sigmaD, 0.5) * pow(sigmaS, 0.5))

					# just done to suppress a warning
					nu ~ dexp(nu.0)

				}"			

	sm.pre.link = "model
				{
					for( i in 1 : n ) {
						tp[i] ~ dbin(tpr[i], n1[i])  
						fp[i] ~ dbin(fpr[i], n2[i])

				 	  m[i,1:2] ~ dmnorm(mu.0[1:2 ], sigma.inv[1:2, 1:2])

				#		Line changed for JAGS 3 compatibility
				 	  w[i] ~ dgamma(nu.2, nu.2) T(0.1, 3)
				#  	w[i] ~ dgamma(nu.2, nu.2) I(0.1, 3)

				  	y[i, 1] <- mu[1] + m[i, 1] / sqrt(w[i])
						y[i, 2] <- mu[2] + m[i, 2] / sqrt(w[i])
						"

	sm.post.link = "
						se[i] <- tpr[i]
				    sp[i] <- 1-fpr[i]

				 }

					# Priors ...
					mu[1] ~ dnorm(m.0[1], pre.mu[1])
					mu[2] ~ dnorm(m.0[2], pre.mu[2])
					mu.0[1] <- 0
					mu.0[2] <- 0

					# Weights distribution
					nu.2 <- nu / 2
					nu ~ dexp(nu.0) # prior for df

					sigma.inv[1:2,1:2] ~ dwish(R[1:2,1:2], k)
					sigma[1:2, 1:2] <- inverse(sigma.inv[1:2, 1:2])

					# Pooled summaries
					x.pool <- (mu[1]+mu[2])/2
					y.pool <- (mu[2]-mu[1])/2
					pool.se <- exp(x.pool) / ( 1 + exp(x.pool) )      # with logit link
					pool.sp <- 1 - exp(y.pool) / ( 1 + exp(y.pool) )  # with logit link

					# Predictive summaries...

				  m.new[1:2]~dmnorm(mu.0[], sigma.inv[1:2 ,1:2])

				#		Line changed for JAGS 3 compatibility 
				#	w.new ~ dgamma(nu.2, nu.2) I(0.01, 5)
					w.new ~ dgamma(nu.2, nu.2) T(0.01, 5)

					y.new[1] <- mu[1] + m.new[1] / sqrt(w.new)
				  y.new[2] <- mu[2] + m.new[2] / sqrt(w.new)

				  new.x  <- (y.new[1] + y.new[2])/2
				  new.y  <- (y.new[2] - y.new[1])/2
				  new.se <- exp(new.x) / (1 + exp(new.x) )
				  new.sp <- 1 - exp(new.y) / (1 + exp(new.y))

					# Variance covariance matrix for random effects
					sigmaD <- sigma[1, 1]
					sigmaS <- sigma[2, 2]
					rhoDS <- sigma[1, 2]/(pow(sigmaD, 0.5) * pow(sigmaS, 0.5))
				}"

	sm.logit.link = "	    
				logit(tpr[i]) <- (y[i, 1] + y[i, 2])/2
		    logit(fpr[i]) <- (y[i, 2] - y[i, 1])/2
				"

	normal.logit.link =	"
				    logit(tpr[i]) <- m[i,1]/2 + m[i,2]/2      # (Di + Si)/2
				    logit(fpr[i]) <- m[i,2]/2 - m[i,1]/2      # (Di - Si)/2
				   	"	

	normal.cloglog.link =	"
						    cloglog(tpr[i]) <- m[i,1]/2 + m[i,2]/2      # (Di + Si)/2
						    cloglog(fpr[i]) <- m[i,2]/2 - m[i,1]/2      # (Di - Si)/2
						  "
		
		# the expression in textConnection(expression) can't be longer
		# then 60 signs... that very odd but it is like it is....

		switch(re,
			normal = 
			switch(link,
				logit = return( textConnection( paste(normal.pre.link
																								,normal.logit.link,
																								normal.post.link,sep=""))),
				cloglog = return( textConnection( paste(normal.pre.link,
																									normal.cloglog.link,
																									normal.post.link,sep=""))) 
				),
			scalemix = 
			switch(link,
				logit = return( textConnection( paste(sm.pre.link, sm.logit.link,sm.post.link,sep=""))),
				cloglog = .throw("Under Construction")
				),
			.throw("The model you requested does not exists. Please contact the developer and we'll what we can do for you.")
		)

	}