

maturity.in.Q <- 40

all.yds <- list()
all.CDS <- list()
all.dates <- list()
all.names <- list()

count <- 0
for(ctry in list.of.ctries){
  count <- count + 1
  
  # load dataset:
  source(paste('load.data.files/load.data.',ctry,'.R',sep=""))
  
  indic.mat.yds <- which(DATA$maturity.yields==maturity.in.Q)
  indic.mat.CDS <- which(DATA$maturity.CDS==maturity.in.Q)
  
  all.yds[[count]]   <- DATA$yields[,indic.mat.yds] - DATA$CDS[,indic.mat.CDS]
  all.CDS[[count]]   <- DATA$CDS[,indic.mat.CDS]
  all.dates[[count]] <- DATA$DATE
  all.names[[count]] <- DATA$name.of.country$long.name
}



FILE = "figures/figure_correlYDS_CDS.pdf"

pdf(file=FILE,pointsize=11,width=8, height=7)

par(mfrow=c(2,2))
par(plt=c(.1,.9,.15,.85))
count <- 0
for(ctry in list.of.ctries){
  count <- count + 1
  
  yds <- 400*all.yds[[count]]
  cds <- 40000*all.CDS[[count]]
  
  # yds <- yds - cds/100
  
  corrl <- cor(yds,cds)
  
  y.lim <- c(min(0,min(yds)),1.3*max(yds))
  
  plot(all.dates[[count]],yds,type="l",ylim=y.lim,
       lwd=2,
       main=paste(all.names[[count]],sep=""),
       # main=paste(all.names[[count]],
       #            " (correlation = ",toString(round(100*corrl,0)),"%)",sep=""),
       las=1,
       xlab="",ylab="")
  par(new=TRUE)
  y.lim <- c(0,max(cds))
  plot(all.dates[[count]],cds,type="l",col="red",lwd=2,lty=2,
       xaxt="n",yaxt="n",xlab="",ylim=y.lim,ylab="")
  axis(4, ylim=y.lim, col="red",col.axis="red",las=1)
  
  if(count==1){
    legend("topright", 
           c(paste(toString(maturity.in.Q/4),"-yr yield, percent",sep=""),
             paste(toString(maturity.in.Q/4),"-yr CDS, in bps",sep="")),
           lty=c(1,2), # gives the legend appropriate symbols (lines)       
           lwd=c(2,2), # line width
           col=c("black","red"),
           bg="white",
           seg.len = 2,
           cex=parcexleg)
  }
}

dev.off()
