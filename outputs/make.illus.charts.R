# =================================================================
# Codes generating illustrative charts.
# =================================================================


FILE = c("figures/formula.pdf")
pdf(file=FILE,pointsize=10,width=5,height=3)

min.x <- -.3
max.x <- .3
x <- seq(min.x,max.x,by=.01)
alpha <- c(.5,1,2)

nb.Q <- 4 # number of quarters

par(mfrow=c(1,1))
par(plt=c(.15,.95,.25,.97))

plot(-x,proba.def.Q.quarters.fction.of.surplus(x,alpha[1],nb.Q),lwd=3,type="l",las=1,ylim=c(0,1),
     xlab=expression(paste("Surplus gap (",s,"*-",s[t],")",sep="")),
     ylab="Probability of default")

polygon(c(x[1]-1,0,0,x[1]-1),c(-1,-1,2,2),
        col=rgb(.85,.85,.85),border = NA)
lines(-x,proba.def.Q.quarters.fction.of.surplus(x,alpha[1],nb.Q),col="black",lwd=2)
lines(-x,proba.def.Q.quarters.fction.of.surplus(x,alpha[2],nb.Q),col="dark grey",lwd=2)
lines(-x,proba.def.Q.quarters.fction.of.surplus(x,alpha[3],nb.Q),col="black",lwd=2,lty=2)

text(x = min.x,y=.05,labels = "Negative surplus gap",pos=4)
text(x = min.x+.6*(max.x-min.x),y=.4,labels = "Positive surplus gap",pos=4)

legend("topright",
       c(expression(paste(alpha," = 0.5",sep="")),
         expression(paste(alpha," = 1.0",sep="")),
         expression(paste(alpha," = 2.0",sep=""))),
       lty=c(1,1,2), # gives the legend appropriate symbols (lines)       
       lwd=2, # line width
       col=c("black","dark grey","black"),
       #title="For dates:",
       seg.len=3,
       cex=1,
       bg="white")


dev.off()

