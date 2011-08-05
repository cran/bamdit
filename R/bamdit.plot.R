#####################################################################
# Special plot function for bamdit
#####################################################################
bamdit.plot <- function(model, data=NULL, S=200)
{
# plot parametrization
nf <- layout(matrix(c(2, 0, 1, 3), 2,2, byrow = TRUE), 
c(3,1), c(1,3), TRUE)
layout.show(nf)

if(!missing(data)){
  .test.data(data)
  # Observed rates...
  tpr <-  data[,1] / data[,2]
  fpr <-  data[,3] / data[,4]
}
   
# Marginally predicted rates...
x <- 1 - as.mcmc(as.matrix(model))[ ,"new.sp"]
y <- as.mcmc(as.matrix(model))[ ,"new.se"]

xhist <- hist(x, plot = F, breaks = 50)
yhist <- hist(y, plot = F, breaks = 50)

top <- max(xhist$counts, yhist$counts)
xrange <- c(min(x), max(x))
yrange <- c(min(y), max(y))
Ind <- sample(1:S)

par(mar=c(3+2.5,3+2.5,1,1))# here I make room for the text
plot(x[Ind], y[Ind], xlim =xrange, ylim=yrange, ylab="TPR (Sensitivity)",
        xlab="FPR (1- Specificity)",
        cex.axis=1.5,
        cex.lab=1.5,
        cex=0.85)
        
if (!missing(data)) points(fpr, tpr, col="blue", cex=3, lwd=2)

par(mar = c(0, 3+2.5, 1, 1))
barplot(xhist$counts, axes = FALSE, ylim = c(0, top), space =0)

par(mar = c(3+2.5, 0, 1, 1))
barplot(yhist$counts, axes = FALSE, xlim = c(0, top), space = 0 , hori =TRUE)
  legend(x="bottomright",legend=c("Observed rates", "Predicted rates"), pch=c("o", "o"),
		col=c("blue", "black"),bty = "n", cex=0.8)
}
#####################################################################
