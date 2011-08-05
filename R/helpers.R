###################################################################
#	INTERNAL HELPERS!!!

# println
.println <- function(...){
	cat(paste(list(...),"\n",sep=""))
}

# throw an error to users face!
.throw <- function(...){
	stop(...,call.=F)
}

# gives back the path as if the file is in the 
# systems tmp-directory
.in.tmp <- function(file){
	paste(tempdir(),file,sep=.Platform$file.sep)
}

# tests data for consistency
.test.data <- function(data){
	# Get the data length
	n <- dim(data)[1]
	# Columns should be:
	# true-pos, real-pos, false-pos, real-neg
	if (dim(data)[2] > 4)
		.throw("Data matrix has more than 4 columns. Please read the help.")
	if (dim(data)[2] < 4)
		.throw("Data matrix has less than 4 columns. Please read the help.")
	# Test for cosistency	
	for (i in 1:n){
		if (data[i,1] > data [i,2])
			.throw("tp > n1 in row ", i, "! Data inconsistent.")
		if (data[i,3] > data [i,4])
			.throw("fp > n2 in row ", i, "! Data inconsistent.")
	}
	return(n)
}

# check model parameters
.check.params <- function(re,link){
  if (re=="normal"){
    if (link=="logit" || link=="cloglog")return(T)
  } else if (re=="scalemix"){
    if (link=="logit")return(T)
  }
  return(F) 
}