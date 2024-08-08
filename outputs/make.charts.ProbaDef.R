# ==============================================================================
# Produces figure with default probabilities
# ==============================================================================

nb.of.CDS.maturities <- length(Horizons)

ColS <- ifelse(indic.colors,rgb(0,0,1),rgb(0.5,0.5,0.5))

################################################################################
ALLDATAPD <- data.frame()
start.min.date <- c("2004-01-01")
end.max.date   <- c("2021-04-01")

vec.dates.df <- seq(as.Date(start.min.date), as.Date(end.max.date), by = "quarter")

TIME <- data.frame(TIME=vec.dates.df)

DATAPD <- data.frame()
################################################################################

# First step, construct names of files

# compute number of files:
nb.ctries <- length(list.of.ctries)
nb.files  <- nb.ctries/max.nb.rows.per.figure
nb.files <- ifelse(nb.files==round(nb.files),nb.files,round(nb.files)+1)

list.names.of.file <- NULL
list.nb.ctries <- list()
ctries.already.done <- 0
for(i in 1:nb.files){
  ctries.on.chart.i <- (ctries.already.done+1):(ctries.already.done+min(nb.ctries-ctries.already.done,max.nb.rows.per.figure))
  ctries.already.done <- ctries.already.done + length(ctries.on.chart.i)
  list.nb.ctries[[i]] <- ctries.on.chart.i
  ctries.string <- list.of.ctries[ctries.on.chart.i[1]]
  for(j in 2:length(ctries.on.chart.i)){
    ctries.string <- paste(ctries.string,list.of.ctries[ctries.on.chart.i[j]],  sep="-")
  }
  list.names.of.file <- c(list.names.of.file,paste("FigureProbaDef", ctries.string,  sep="_"))
}


for(j in 1:length(list.names.of.file)){
  
  list_of_countries <- list.nb.ctries[[j]]
  nb.ctries.on.chart.j <- length(list_of_countries)
  
  FILE = paste("figures/",list.names.of.file[j],suffix,
               ".pdf",sep="")
  
  pdf(file=FILE,pointsize=11,width=9, height=2*nb.ctries.on.chart.j)
  
  par(mfrow=c(nb.ctries.on.chart.j,3))
  
  par(plt=c(.1,.95,.25,.77))
  
  counter.ctry <- 0
  
  for(ctry in list.of.ctries[list_of_countries]){
    counter.ctry <- counter.ctry + 1
    
    path.results <- paste("results/res_estim_",ctry,suffix,".Rdat",sep="")
    load(file=path.results)
    
    vec.dates <- DATA$DATE
    
    # Construct estimated model:
    model.est <- Theta2Model(full_theta_est,model_ini_sol,bounds)
    model.sol.est <- solve_model(model.est,indic_compute_beta = 1)
    
    n_w <- dim(model.sol.est$Phi)[1]
    n_x <- dim(model.sol.est$Phi_x)[1]
    
    resEKF_est <- EKF_smoother(model.sol.est,list.stdv,DATA)
    
    T     <- dim(DATA$CDS)[1]
    all.X <- matrix(0,T,n_x)
    all.X[,1:n_w] <- t(resEKF_est$W.smoothed)
    all.X[,n_w+1] <- DATA$r
    all.X[,n_w+2] <- DATA$d
    all.X[,n_w+3] <- c(DATA$d[1],DATA$d[1:(T-1)])
    
    DATA$maturity.CDS <- Horizons
    
    Proba.def.Q <- PDef(model.sol.est,all.X,
                        max.h = max(DATA$maturity.CDS),
                        indic.P = 0)
    Proba.def.P <- PDef(model.sol.est,all.X,
                        max.h = max(DATA$maturity.CDS),
                        indic.P = 1)

    CRP <- Proba.def.Q - Proba.def.P
    cor.P.CRP <- diag(cor(Proba.def.P, CRP))
    
    COUNTRY <- rep(ctry, length(TIME))
    
    eval(parse(text = gsub(" ","",paste("DATAPD.i <- data.frame(TIME=vec.dates,
                                      PD.P = Proba.def.P,
                                      PD.Q = Proba.def.Q,
                                      CRP = CRP)", sep=""))))
    
    DATAPDct <- data.frame(TIME,COUNTRY)
    DATAPD_by_country <- merge(DATAPDct,DATAPD.i,by="TIME",all = TRUE)
    
    DATAPD <- rbind(DATAPD,DATAPD_by_country)
    
    for(i.maturity in nb.of.CDS.maturities:1){
      
      count.measure <- 0
      for(measure in c("P","Q","diff")){
        
        count.measure <- count.measure + 1
        
        AUX.all.P <- Proba.def.P[,DATA$maturity.CDS]
        AUX.P     <- Proba.def.P[,DATA$maturity.CDS[i.maturity]]
        AUX.all.Q <- Proba.def.Q[,DATA$maturity.CDS]
        AUX.Q     <- Proba.def.Q[,DATA$maturity.CDS[i.maturity]]
        
        P.def.P <- 100*AUX.P
        P.def.Q <- 100*AUX.Q
        AUX.all.Q <- 100 * AUX.all.Q
        
        min.y <- 0
        max.y <- 1.0*max(AUX.all.Q)
        if(ctry == "CA"){
          max.y <- min(10,max.y)
        }
        
        par(mfg=c(counter.ctry,count.measure))
        
        Title <- paste(DATA$name.of.country$long.name,", under ",measure,sep="")
        if(measure == "diff"){
          Title <- "Risk premium"
        }
        if(measure == "P"){
          y <- P.def.P
        }else if(measure == "Q"){
          y <- P.def.Q
        }else{
          y <- P.def.Q - P.def.P
        }
        
        if(i.maturity==nb.of.CDS.maturities){
          plot(vec.dates,y,type="l",
               las=1,lwd=2,lty=1,
               ylim=c(min.y,max.y),
               xlab="",ylab="",col="dark grey",
               main = Title,
               cex.main=parcex, cex.lab=parcex, cex.axis=parcex)
        }else{
          lines(vec.dates,y,col="black",lwd=2,lty=i.maturity) 
        }
        
        #abline(h=0,col="grey",lty=1)
        
        if(i.maturity==1){
          vec.maturities <- paste(toString(DATA$maturity.CDS[1]/4),ifelse(DATA$maturity.CDS[1]/4<=1," year"," years"),sep="")
          for(jj in 2:nb.of.CDS.maturities){
            vec.maturities <- c(paste(toString(DATA$maturity.CDS[jj]/4)," years",sep=""),
                                vec.maturities)
          }
          
          if((counter.ctry==1)*(i.maturity==1)){
            legend("topright", # places a legend at the appropriate place c("Health","Defense"), # puts text in the legend 
                   vec.maturities,
                   lty=1, # gives the legend appropriate symbols (lines)       
                   lwd=c(2,2), # line width
                   col=c("dark grey","black"),
                   bg="white",
                   pch=c(NaN,NaN,NaN),
                   #pt.bg=c(NaN,"grey"),
                   #pt.cex = c(NaN,2),
                   seg.len = 2,
                   ncol=1,
                   # xpd=TRUE, inset=c(0, -0.1),
                   cex=parcexleg)
          }
        }
      }
    }
  }
  
  dev.off()
  
  
}

write.csv(DATAPD,file=paste("results/PD_panel.csv",sep=""))

