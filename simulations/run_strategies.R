# ==============================================================================
# Analysis of public debt management strategies
# ==============================================================================


outputs <- c("mean_d","stdv_d","DaR95","mean_rr","stdv_rr",
             "stdv_Delta_d","avg_PD[maxH]","avg_spreads[maxH]")
# Convert the names of output into "expressions" for charts labels & titles:
outputs4charts <- c(expression(paste(E(d),sep="")),
                    expression(paste(sqrt(V(d)),sep="")),
                    expression(paste(q[95](d),sep="")),
                    expression(paste(E(r),sep="")),
                    expression(paste(sqrt(V(r)),sep="")),
                    expression(paste(sqrt(Delta(d)),sep="")),
                    expression(paste(E(PD[10]),sep="")),
                    expression(paste(E(spd[10]),sep="")))

grids <- make_grid(nb_grid = nb_grid,
                   min_d  = min_d,
                   max_d  = max_d,
                   min_rr = min_rr,
                   max_rr = max_rr,
                   sigma_eps = Model$sigma_eps,
                   all_quantiles_eps = c(-2,-1,1,2))

M <- NULL
parameters <- NULL

for(chi in values_of_chi){
  for(kappa_pi in values_of_kappa_pi){
    for(kappa_y in values_of_kappa_y){
      
      print(paste("chi = ",chi," , kappa_pi = ",kappa_pi,
                  " , kappa_y = ",kappa_y,sep=""))
      Model_i <- Model
      Model_i$kappa_pi <- kappa_pi
      Model_i$kappa_y  <- kappa_y
      Model_i$chi      <- chi
      
      # ============================================
      # To model liquidity risks:
      Model_i$delta <- Model$delta - .0 * kappa_pi
      # ============================================

      Model_solved_i <- solve_ToyModel(Model_i,grids,
                                       nb_iter = nb_iter,
                                       nb_iter_sdf = nb_iter_sdf)
      
      # -----------------------------------
      p <- compute_uncond_distri(Model_solved_i$indicators_x,
                                 Model_solved_i$Probas,1000)
      distri_d  <- compute_distri_x(grids$all_d,Model_solved_i$d,p)
      plot(grids$all_d,distri_d,type="l")
      distri_rr  <- compute_distri_x(grids$all_rr,Model_solved_i$rr,p)
      plot(grids$all_rr,distri_rr,type="l")
      # -----------------------------------
      
      strat_i <- run_strategy(Model_solved_i,maxH=10,
                              nb_iter_sdf = nb_iter_sdf)
      
      thisLine <- NULL
      for(output in outputs){
        eval(parse(text = gsub(" ","",paste("aux <- strat_i$",output,sep=""))))
        thisLine <- c(thisLine,aux)
      }
      
      M <- rbind(M,thisLine)
      parameters <- rbind(parameters,
                          c(chi,kappa_pi,kappa_y))
    }
  }
}

colnames(parameters) <- c("chi","kappa_pi","kappa_y")
colnames(M)          <- outputs

plot(M[,"mean_d"],M[,"stdv_d"])
plot(M[,"mean_rr"],M[,"avg_PD[maxH]"])
plot(M[,"mean_rr"],M[,"DaR95"])

plot(M[(parameters[,"kappa_y"]==0),"mean_d"],M[(parameters[,"kappa_y"]==0),"stdv_d"])
plot(M[(parameters[,"kappa_y"]==0)&(parameters[,"chi"]==0.9),"mean_rr"],
     M[(parameters[,"kappa_y"]==0)&(parameters[,"chi"]==0.9),"avg_PD[maxH]"])




# ==============================================================================
# Prepare charts
# ==============================================================================


# Prepare colors:
nb_kappa_pi <- length(values_of_kappa_pi)
nb_kappa_y  <- length(values_of_kappa_y)
nb_chi      <- length(values_of_chi)

pal1 <- colorRampPalette(c("white", rgb(1,0,0)), space = "rgb")
col1 <- val2col(values_of_kappa_pi, col=pal1(nb_kappa_pi))
pal2 <- colorRampPalette(c("white", rgb(0,0,1)), space = "rgb")
col2 <- val2col(values_of_kappa_y, col=pal2(nb_kappa_y))
col3 <- NA*seq(nrow(parameters))
pch3 <- NA*seq(nrow(parameters))
for(i in seq(nrow(parameters))){
  xpos <- which(parameters[i,"kappa_pi"]==values_of_kappa_pi)
  ypos <- which(parameters[i,"kappa_y"]==values_of_kappa_y)
  coltmp <- (col2rgb(col1[xpos])/2) + (col2rgb(col2[ypos])/2)
  col3[i] <- rgb(coltmp[1], coltmp[2], coltmp[3], maxColorValue = 255)
  pch3[i] <- 20 + which(parameters[i,"chi"]==values_of_chi)
}

names_legend_4_chi <- NULL
for(i in 1:nb_chi){
  names_legend_4_chi <- c(names_legend_4_chi,as.character(values_of_chi[i]))
}


FILE = paste("figures/Figure_strategies_perf.pdf",sep="")
pdf(file=FILE, pointsize=10, width=6, height=6)

#plot
layout(t(matrix(c(1,2,3,4,5,6,7,7,7,7), 5, 2)),
       widths=c(4,1,1,1,1), heights=c(4,4), respect=T)

count_chart <- 0

for(quadrant in 1:4){
  
  if(quadrant == 2){
    par(mar=c(4,0,3,5))
    plot.new()
    
    image(x=1, y=values_of_kappa_pi, z=t(as.matrix(values_of_kappa_pi)),
          col=pal1(nb_kappa_pi), xaxt="n", yaxt="n", xlab="", ylab="")
    box()
    axis(4,las=1,at=values_of_kappa_pi)
    mtext(expression(paste("Inflation indexation ",kappa[pi],sep="")), side=4, line=3, cex=0.7)
    image(x=1, y=values_of_kappa_y, z=t(as.matrix(values_of_kappa_y)),
          col=pal2(nb_kappa_pi), xaxt="n", yaxt="n", xlab="", ylab="")
    box()
    axis(4,las=1,at=values_of_kappa_y)
    mtext(expression(paste("Real GDP indexation ",kappa[y],sep="")), side=4, line=3, cex=0.7)
    
    plot.new()
  }else{
    
    count_chart <- count_chart + 1
    
    # Define x and y axes:
    xvariable <- outputs4chart[count_chart,1]
    yvariable <- outputs4chart[count_chart,2]
    
    x <- M[,xvariable]
    y <- M[,yvariable]
    
    xlim <- c(min(x),1.1*max(x))
    ylim <- c(min(y),max(y))
    
    par(mar=c(4,5,3,2))
    plot(0,0,xlim=xlim,ylim=ylim,
         xlab=outputs4charts[which(xvariable==outputs)],
         ylab=outputs4charts[which(yvariable==outputs)],las=1,
         main = paste("(",letters[count_chart],")",sep=""))
    for(i in seq(nrow(parameters))){
      points(x[i],y[i],bg=col3[i],pch=pch3[i],col="dark grey",cex=2)
    }
    if(count_chart==1){
      legend("topright",
             names_legend_4_chi,
             lty=rep(NaN,nb_chi), # gives the legend appropriate symbols (lines)
             lwd=1, # line width
             col=c("black"),
             bg="white",
             pch=seq(21,21+nb_chi-1),
             seg.len = 2,
             cex=1,
             pt.cex = 2,
             title=expression(paste("Coupon decay rate ",chi,sep="")))
    }
  }
}

dev.off()

