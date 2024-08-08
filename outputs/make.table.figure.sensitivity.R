# =================================================================
# Compute counterfactual exercises (to study spread sensitivity)
# =================================================================

CR <- read.csv2("data/CreditRatings/creditratings.csv")
scale.Moodys <- c("Aaa",
                  "Aa1","Aa2","Aa3",
                  "A1","A2","A3",
                  "Baa1","Baa2","Baa3",
                  "Ba1","Ba2","Ba3",
                  "B1","B2","B3")
scale.SP <- c("AAA",
              "AA+","AA","AA-",
              "A+","A","A-",
              "BBB+","BBB","BBB-",
              "BB+","BB","BB-",
              "B+","B","B-")
scale.SP.4.charts <- c("AAA",
                       "AA+","AA","AA-",
                       "A+","A","A-",
                       "BBB+","BBB","BBB-",
                       "BB+","BB","BB-",
                       "B+","B","B-")
scale.Fitch <- scale.SP


par(plt=c(.1,.9,.1,.9))

Horizon <- 4 # horizon, in quarters.

all.types.date                     <- c("last","min FS")
all.values.of.surpluses.in.percent <- seq(0,20,by=1)
deficit.in.table                   <- c(1,5,10)
all.maturities                     <- c(8,20,40) # in quarters


all.tightest.FS.dates <- NULL

# LOOP ON MATURITIES
for(considered.maturity.of.CDS in all.maturities){
  
  print("====================================")
  print(paste("Processing ",toString(considered.maturity.of.CDS/4),"-year CDS",sep=""))
  
  all.unit.effects  <- array(NaN,
                             c(length(list.of.ctries),
                               length(deficit.in.table),
                               length(all.types.date)))
  all.total.effects <- array(NaN,
                             c(length(list.of.ctries),
                               length(deficit.in.table),
                               length(all.types.date)))
  
  all.Moodys.ratings <- NULL
  all.SP.ratings     <- NULL
  all.Fitch.ratings  <- NULL
  
  all.long.names <- NULL
  
  par(mfrow=c(4,3))
  
  count.ctry <- 0
  
  for(ctry in list.of.ctries){
    
    count.ctry <- count.ctry + 1
    
    indic.in.rating.table <- which(CR[,1]==ctry)
    Moodys.rating <- CR[indic.in.rating.table,2]
    SP.rating     <- CR[indic.in.rating.table,3]
    Fitch.rating  <- CR[indic.in.rating.table,4]
    all.Moodys.ratings <- c(all.Moodys.ratings,which(Moodys.rating==scale.Moodys))
    all.SP.ratings     <- c(all.SP.ratings,    which(SP.rating    ==scale.SP))
    all.Fitch.ratings  <- c(all.Fitch.ratings, which(Fitch.rating ==scale.Fitch))
    
    path.results <- paste("results/res_estim_",ctry,suffix,".Rdat",sep="")
    load(file=path.results)
    
    vec.dates <- DATA$DATE
    
    # Construct estimated model:
    model.est     <- Theta2Model(full_theta_est,model_ini_sol,bounds)
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
    
    X0 <- all.X
    
    all.long.names <- c(all.long.names,
                        DATA$name.of.country$long.name)
    
    # Compute fiscal space:
    elld <- DATA$d/4
    ell1 <- compute.debt.limits(model.sol.est,all.X,Proba=.01,h=4)/model.sol.est$freq
    FS <- ell1 - elld

    count.dates <- 0
    
    for(index.type.date in 1:length(all.types.date)){
      
      count.dates <-   count.dates + 1
      
      if(all.types.date[index.type.date]=="last"){
        T <- dim(X0)[1]
      }else if(all.types.date[index.type.date]=="min FS"){
        T <- which(FS==min(FS))
        if(length(all.tightest.FS.dates)<length(list.of.ctries)){
          all.tightest.FS.dates <- c(all.tightest.FS.dates,DATA$DATE[T])
        }
      }else{
        print("Problem")
      }
      
      modeled.CDS.Q.0    <- compute.CDS(model.sol.est,
                                        X0,40,
                                        indic.cpp = 1)$CDS.spds
      
      all.CDS.chges <- NULL
      
      for(surplus.in.percent in all.values.of.surpluses.in.percent){
        
        # Compute CDSs under scenario:
        X1 <- X0
        X1[T,n_w+2] <- X1[T,n_w+2] + surplus.in.percent*4/100

        modeled.CDS.Q.1    <- compute.CDS(model.sol.est,
                                          X1,40,
                                          indic.cpp = 1)$CDS.spds
        
        modeled.cds.q <- 40000*(modeled.CDS.Q.1[T,considered.maturity.of.CDS] -
                                  modeled.CDS.Q.0[T,considered.maturity.of.CDS])
        
        all.CDS.chges <- c(all.CDS.chges,modeled.cds.q)
      }
      
      x.axis <- all.values.of.surpluses.in.percent[2:length(all.values.of.surpluses.in.percent)]
      effect.unit.increase <- all.CDS.chges[2:length(all.values.of.surpluses.in.percent)] -
        all.CDS.chges[1:(length(all.values.of.surpluses.in.percent)-1)]
      effect.total <- all.CDS.chges[2:length(all.values.of.surpluses.in.percent)]
      
      all.unit.effects[count.ctry,,count.dates]  <- effect.unit.increase[deficit.in.table]
      all.total.effects[count.ctry,,count.dates] <- effect.total[deficit.in.table]
      
      if(index.type.date==1){
        plot(x.axis,effect.unit.increase,type="l",main=ctry,
             ylim=c(0,2*max(effect.unit.increase)),lwd=2)
      }else{
        lines(x.axis,effect.unit.increase,
              #lty=,
              col="grey",
              lwd=ifelse(index.type.date==2,1,3))
      }
    }
    
  }
  
  
  # --------------------------------
  # Prepare Latex table
  # --------------------------------
  
  name.of.file <- paste("tables/table_sensitivity_",
                        toString(considered.maturity.of.CDS/4),"yr",sep="")
  
  Top.of.Table <- rbind(
    "\\begin{table}[H]",
    paste("\\caption{",toString(considered.maturity.of.CDS/4),"-year CDS sensitivity to deficits}",sep=""),
    paste("\\label{tab:sensitivity_",toString(considered.maturity.of.CDS/4),"yr}",sep=""),
    "\\begin{footnotesize}",
    "\\begin{tabular*}{\\textwidth}{l@{\\extracolsep{\\fill}}lrrrrrrrrrr}",
    "")
  
  latex.table <- NULL
  
  count.dates <- 0
  
  for(index.type.date in all.types.date){
    
    count.dates <- count.dates + 1
    
    all.lines.4.types.date <- NULL
    
    count.ctry <- 0
    
    for(ctry in list.of.ctries){
      
      count.ctry <-   count.ctry + 1
      
      this.line <- paste(all.long.names[count.ctry],sep="")
      
      first.line.with.deficits <- "Fiscal shock:"
      
      for(count.deficit in 1:length(deficit.in.table)){
        
        tot.effect <- all.total.effects[count.ctry,
                                        count.deficit,
                                        count.dates]
        unit.effect <- all.unit.effects[count.ctry,
                                        count.deficit,
                                        count.dates]
        
        first.line.with.deficits <- paste(first.line.with.deficits,"&&",
                                          "\\multicolumn{2}{c}{+",toString(deficit.in.table[count.deficit])," p.p. of GDP}",
                                          sep="")
        
        if(tot.effect>1000){
          this.line <- c(paste(this.line,"&&",
                               round.fixed.length(tot.effect/1000,1),"$\\times10^3$",
                               "&[\\textit{",
                               round.fixed.length(unit.effect/1000,1),"$\\times10^3$",
                               "}]",
                               sep=""))
        }else{
          this.line <- c(paste(this.line,"&&",
                               round.fixed.length(tot.effect,1),
                               "&[\\textit{",
                               round.fixed.length(unit.effect,1),
                               "}]",
                               sep=""))
        }
      }
      
      all.lines.4.types.date <- rbind(all.lines.4.types.date,
                                      paste(this.line,"\\\\ ", sep=""))
      
    }
    
    if(all.types.date[count.dates]=="last"){
      T <- dim(X0)[1]
    }else if(all.types.date[count.dates]=="min FS"){
      T <- which(FS==min(FS))
    }
    
    title.line <- paste("\\multicolumn{10}{l}{{\\bf Panel ",LETTERS[count.dates]," - Date: ",
                        DATA$DATE[T],"}}",sep="")
    
    latex.table <- rbind(latex.table,
                         "\\hline",
                         title.line,
                         "\\\\",
                         paste(first.line.with.deficits,"\\\\",sep=""),
                         "\\cmidrule(lr){3-4}\\cmidrule(lr){6-7}\\cmidrule(lr){9-10}",
                         all.lines.4.types.date,
                         "\\\\")
  }
  
  
  latex.table <- rbind(Top.of.Table,
                       "\\hline",
                       latex.table)
  
  latex.table <- rbind(latex.table,
                       "\\hline",
                       "\\hline",
                       "\\end{tabular*}",
                       "\\end{footnotesize}",
                       "\\begin{footnotesize}",
                       "\\parbox{\\linewidth}{",
                       paste("Note: This table documents the sensitivity of the ",
                             toString(considered.maturity.of.CDS/4),"-year CDS spreads to fiscal conditions.
                       We consider three sizes of fiscal shocks (increases in primary deficits by 1\\%, 5\\% and 10\\% of GDP).
                       The reported figures are in basis points.
                       The number in square brackets correspond to the marginal influence of an additional unit increase in the deficit.
                       Panel A reports the results for the last quarter of the estimation sample; Panel B corresponds to the quarter featuring the smallest fiscal space.
                       ",sep=""),
                       "}",
                       "\\end{footnotesize}",
                       "\\end{table}")
  
  latex.file <- paste(name.of.file,suffix,".txt", sep="")
  write(latex.table, latex.file)
  
  
  
  # ============================================================================
  # Produce chart
  # ============================================================================
  
  pch.agency <- c(1,4,19)
  cex.agency <- c(2,1.8,1)
  
  all.ratings <- rbind(all.Moodys.ratings,
                       all.SP.ratings,
                       all.Fitch.ratings)
  
  FILE = paste("figures/Figure_CDS_Sensitivity_",
               toString(considered.maturity.of.CDS/4),"yr",suffix,".pdf",sep="")
  pdf(file=FILE,pointsize=10,width=9, height=9)
  
  par(mfrow=c(2,2))
  
  par(plt=c(.15,.95,.2,.8))
  
  select.choice.of.date.type     <- c(1,2)
  select.choice.of.deficit.type  <- c(1,3)
  
  for(choice.of.date.type in select.choice.of.date.type){
    for(choice.of.deficit.type in select.choice.of.deficit.type){
      
      plot(1:length(scale.Moodys),rep(0,length(scale.Moodys)),
           ylim=c(-2,log(50)),col="white",
           xlim=c(0,14),
           xaxt="n",yaxt="n",
           xlab="",
           cex.lab=1.3,
           #ylab="CDS sensitivity (in bps)",
           ylab="")
      axis(side=1,
           at=1:length(scale.Moodys),
           labels=scale.SP.4.charts,
           #pos=,
           lty=1,
           col="black",
           cex.axis=1.3,
           las=2
           #tck=
      )
      tick.4.grid <- c(.1,1,2,5,10,20,50,100,200,500,1000,2000,5000)
      tick.4.grid.string <- NULL
      for(i in 1:length(tick.4.grid)){
        tick.4.grid.string <- c(tick.4.grid.string,toString(tick.4.grid[i]))
        abline(h=log(tick.4.grid[i]),col="light grey",lty=3)
      }
      axis(side=2,
           at=log(tick.4.grid),
           labels=tick.4.grid.string,
           #pos=,
           lty=1,
           col="black",
           cex.axis=1.3,
           las=2
           #tck=
      )
      
      count.ctry <- 0
      for(ctry in list.of.ctries){
        count.ctry <- count.ctry + 1
        
        for(agency in 1:3){
          y.values <- log(all.unit.effects[,choice.of.deficit.type,
                                           choice.of.date.type])
          if(agency==1){
            text(x=max(all.ratings[,count.ctry])+.5,y.values[count.ctry],ctry)
          }
          points(all.ratings[agency,count.ctry],y.values[count.ctry],
                 xlab="Ratings",ylab="CDS sensitivity",
                 pch=pch.agency[agency],
                 lwd=2,
                 cex=cex.agency[agency])
        }
        
        ratings <- sort(all.ratings[,count.ctry])
        # detect different ratings:
        marg <- .15
        different.ratings <- sort(as.numeric(levels(as.factor(ratings))))
        if(length(different.ratings)==2){
          lines(c(different.ratings[1]+marg,different.ratings[2]-marg),
                rep(y.values[count.ctry],2))
        }
        if(length(different.ratings)==3){
          lines(c(different.ratings[1]+marg,different.ratings[2]-marg),
                rep(y.values[count.ctry],2))
          lines(c(different.ratings[2]+marg,different.ratings[3]-marg),
                rep(y.values[count.ctry],2))
        }
        
      }
      
      if((choice.of.date.type==1)&(choice.of.deficit.type==1)){
        legend("topleft", 
               c("Moody's","S&P","Fitch"),
               lty=c(NaN,NaN,NaN), # gives the legend appropriate symbols (lines)
               lwd=2, # line width
               col=c("black"),
               bg="white",
               pch=pch.agency,
               seg.len = 1,
               ncol=1,
               title="Rating agency:",
               cex=1.4,
               pt.cex=cex.agency)
      }
    }
  }
  
  count.dates.ini <- 0
  for(choice.of.date.type in select.choice.of.date.type){
    count.dates.ini <- count.dates.ini + 1
    par(mfrow=c(2,1))
    par(plt=c(.05,.95,.05,.86))
    par(mfg=c(count.dates.ini,1))
    title(
      paste("Panel ",LETTERS[count.dates.ini]," - ",
            ifelse(all.types.date[choice.of.date.type]=="last",
                   "Last date in sample",
                   "Date with lowest fiscal space"),
            sep="")
    )
    count.deficit.shock <- 0
    for(choice.of.deficit.type in select.choice.of.deficit.type){
      count.deficit.shock <-   count.deficit.shock + 1
      par(mfrow=c(2,2))
      par(plt=c(.05,.95,.05,.8))
      par(mfg=c(count.dates.ini,count.deficit.shock))
      title(
        paste(LETTERS[count.dates.ini],".",toString(count.deficit.shock),
              " Deficit shock: +",toString(deficit.in.table[choice.of.deficit.type])," p.p. of GDP",sep="")
      )
    }
  }
  
  dev.off()
  
  
}