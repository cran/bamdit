read.RevMan <- function(file){
	tmp <- read.csv(file, header = T, sep = ",")
	ret <- matrix(data = NA, nrow = dim(tmp)[1], ncol = 4)
	colnames(ret)<-c("tp","n1","fp","n2")
	for (i in 1:dim(tmp)[1]) {
		ret[i,1] = tmp$TP[i]
		ret[i,2] = tmp$TP[i] + tmp$FN[i]
		ret[i,3] = tmp$FP[i]
		ret[i,4] = tmp$FP[i] + tmp$TN[i]
	}
	return(ret)
}