
if(!exists("outputs")){
  outputs <- c("mean_d","stdv_d","DaR95","mean_rr","stdv_rr",
               "stdv_Delta_d","avg_PD[maxH]","avg_spreads[maxH]")
}
if(!exists("outputs4charts")){
  outputs4charts <- c(expression(paste(E(d),sep="")),
                      expression(paste(sqrt(V(d)),sep="")),
                      expression(paste(q[95](d),sep="")),
                      expression(paste(E(r),sep="")),
                      expression(paste(sqrt(V(r)),sep="")),
                      expression(paste(sqrt(Delta(d)),sep="")),
                      expression(paste(E(PD[10]),sep="")),
                      expression(paste(E(spd[10]),sep="")))
}
if(!exists("outputs4chart")){
  outputs4chart <- matrix(NaN,3,2)
  outputs4chart[,1] <- "mean_d"
  outputs4chart[,2] <- c("stdv_d","DaR95","avg_PD[maxH]")
}

matrix_values_4arrows <- NULL
matrix_values_4arrows <- rbind(matrix_values_4arrows,
                               c(.9,1,.3))
matrix_values_4arrows <- rbind(matrix_values_4arrows,
                               c(.9,0,.3))
matrix_values_4arrows <- rbind(matrix_values_4arrows,
                               c(.9,1,0))
matrix_values_4arrows <- rbind(matrix_values_4arrows,
                               parameters[which.min(M[,"DaR95"]),])
matrix_values_4arrows <- rbind(matrix_values_4arrows,
                               parameters[which.min(M[,"stdv_d"]),])
matrix_values_4arrows <- rbind(matrix_values_4arrows,
                               parameters[which.min(M[,"avg_PD[maxH]"]),])
colnames(matrix_values_4arrows) <- c("chi","kappa_pi","kappa_y")


# Identify the single strategy approximating the effective US debt portfolio.
tol_strategy <- 1e-8
is_us_strategy <- rep(FALSE, nrow(parameters))
if(exists("additional_strategies")){
  us_parameters <- as.numeric(additional_strategies[1, c("chi", "kappa_pi", "kappa_y")])
}else if(exists("Model")){
  us_parameters <- c(Model$chi, 1935/25734, 0)
}else{
  us_parameters <- c(NA_real_, NA_real_, NA_real_)
}
if(all(is.finite(us_parameters))){
  is_us_strategy <- apply(parameters, 1, function(x){
    all(abs(x - us_parameters) < tol_strategy)
  })
}
regular_strategy <- !is_us_strategy

# Prepare colors.  The US strategy is not used to define coupon-decay symbols,
# because it is shown separately as a large cross.
values_of_chi      <- sort(unique(parameters[regular_strategy,"chi"]))
values_of_kappa_pi <- sort(unique(parameters[regular_strategy,"kappa_pi"]))
values_of_kappa_y  <- sort(unique(parameters[regular_strategy,"kappa_y"]))

nb_kappa_pi <- length(values_of_kappa_pi)
nb_kappa_y  <- length(values_of_kappa_y)
nb_chi      <- length(values_of_chi)
kappa_pi_ticks <- seq(min(values_of_kappa_pi), max(values_of_kappa_pi), by=.2)
kappa_y_ticks <- values_of_kappa_y

pal1 <- colorRampPalette(c("white", rgb(1,0,0)), space = "rgb")
col1 <- val2col(parameters[,"kappa_pi"],
                zlim=range(values_of_kappa_pi),
                col=pal1(nb_kappa_pi))
pal2 <- colorRampPalette(c("white", rgb(0,0,1)), space = "rgb")
col2 <- val2col(parameters[,"kappa_y"],
                zlim=range(values_of_kappa_y),
                col=pal2(nb_kappa_y))
col3 <- NA*seq(nrow(parameters))
pch3 <- NA*seq(nrow(parameters))
for(i in seq(nrow(parameters))){
  coltmp <- (col2rgb(col1[i])/2) + (col2rgb(col2[i])/2)
  col3[i] <- rgb(coltmp[1], coltmp[2], coltmp[3], maxColorValue = 255)
  if(is_us_strategy[i]){
    pch3[i] <- 4
  }else{
    pch3[i] <- 20 + which(abs(parameters[i,"chi"] - values_of_chi) < tol_strategy)
  }
}

names_legend_4_chi <- NULL
for(i in 1:nb_chi){
  names_legend_4_chi <- c(names_legend_4_chi,sprintf("%.3f", values_of_chi[i]))
}

find_strategy_rows <- function(targets){
  rows <- integer(0)
  for(i in seq_len(nrow(targets))){
    rows <- c(rows, which(apply(parameters, 1, function(x){
      all(abs(x - targets[i,]) < tol_strategy)
    })))
  }
  unique(rows)
}

draw_annotations <- function(x, y, annotation_rows, is_us_annotation,
                             range_x, range_y, mid_x, xlim, ylim){
  if(length(annotation_rows) == 0){
    return(invisible(NULL))
  }
  ordered <- order(y[annotation_rows], decreasing = TRUE)
  annotation_rows <- annotation_rows[ordered]
  is_us_annotation <- is_us_annotation[ordered]
  n_annotations <- length(annotation_rows)
  x_label <- rep(max(x) + .12 * range_x, n_annotations)
  y_label <- seq(ylim[2] - .10 * range_y,
                 ylim[1] + .10 * range_y,
                 length.out = n_annotations)
  
  for(j in seq_along(annotation_rows)){
    row_id <- annotation_rows[j]
    arrows(x0 = x_label[j] - .025 * range_x,
           y0 = y_label[j],
           x1 = x[row_id],
           y1 = y[row_id],
           length = .07,
           angle = 12,
           lwd = 1,
           col = "grey35")
    if(!is_us_annotation[j]){
      label_pi <- substitute(kappa[pi] == KP,
                             list(KP = round(parameters[row_id,"kappa_pi"], 2)))
      label_y <- substitute(kappa[y] == KY,
                            list(KY = round(parameters[row_id,"kappa_y"], 2)))
      text(x_label[j], y_label[j],
           labels = label_pi,
           adj = c(0, .5),
           cex = 1.08)
      text(x_label[j], y_label[j] - .036 * range_y,
           labels = label_y,
           adj = c(0, .5),
           cex = 1.08)
    }
  }
  invisible(NULL)
}


FILE = paste("figures/Figure_strategies_perf.pdf",sep="")
pdf(file=FILE, pointsize=10, width=8, height=8)

#plot
layout(t(matrix(c(1,2,3,4,5,6,7,7,7,7), 5, 2)),
       widths=c(4,1,1,1,1), heights=c(4,4), respect=T)

count_chart <- 0

for(quadrant in 1:4){
  
  if(quadrant == 2){
    par(mar=c(4,0,3,5))
    plot.new()
    
    image(x=1, y=values_of_kappa_pi, z=t(as.matrix(seq_along(values_of_kappa_pi))),
          col=pal1(nb_kappa_pi), xaxt="n", yaxt="n", xlab="", ylab="")
    box()
    axis(4,las=1,at=kappa_pi_ticks,labels=sprintf("%.1f", kappa_pi_ticks))
    mtext(expression(paste("Inflation indexation ",kappa[pi],sep="")), side=4, line=3, cex=0.7)
    image(x=1, y=values_of_kappa_y, z=t(as.matrix(seq_along(values_of_kappa_y))),
          col=pal2(nb_kappa_y), xaxt="n", yaxt="n", xlab="", ylab="")
    box()
    axis(4,las=1,at=kappa_y_ticks,labels=sprintf("%.1f", kappa_y_ticks))
    mtext(expression(paste("Real GDP indexation ",kappa[y],sep="")), side=4, line=3, cex=0.7)
    
    plot.new()
  }else{
    
    count_chart <- count_chart + 1
    
    # Define x and y axes:
    xvariable <- outputs4chart[count_chart,1]
    yvariable <- outputs4chart[count_chart,2]
    
    x <- M[,xvariable]
    y <- M[,yvariable]
    
    range_x <- max(x) - min(x)
    range_y <- max(y) - min(y)
    mid_x <- (max(x) + min(x))/2
    xlim <- c(min(x) - .18*range_x,max(x) + .72*range_x)
    ylim <- c(min(y),max(y))
    if(count_chart == 1){
      ylim[2] <- ylim[2] + .45 * range_y
    }
    
    par(mar=c(4,5,3,2))
    plot(0,0,xlim=xlim,ylim=ylim,
         xlab=outputs4charts[which(xvariable==outputs)],
         ylab=outputs4charts[which(yvariable==outputs)],las=1,
         main = paste("(",letters[count_chart],")",sep=""))
    grid()
    for(i in which(!is_us_strategy)){
      points(x[i], y[i], bg = col3[i], pch = pch3[i],
             col = "dark grey", cex = 2)
    }
    
    annotation_rows <- find_strategy_rows(matrix_values_4arrows)
    annotation_rows <- setdiff(annotation_rows, which(is_us_strategy))
    is_us_annotation <- is_us_strategy[annotation_rows]
    draw_annotations(x, y, annotation_rows, is_us_annotation,
                     range_x, range_y, mid_x, xlim, ylim)
    
    if(any(is_us_strategy)){
      i_us <- which(is_us_strategy)[1]
      points(x[i_us], y[i_us], pch = 4, col = "black",
             cex = 2.8, lwd = 3)
    }
    
    if(count_chart==1){
      legend("topleft",
             names_legend_4_chi,
             lty = rep(NaN, nb_chi), # gives the legend appropriate symbols
             lwd = 1,
             col = "black",
             bg = "white",
             pch = seq(21,21+nb_chi-1),
             pt.lwd = rep(1, nb_chi),
             seg.len = 2,
             cex = 1,
             pt.cex = rep(2, nb_chi),
             title = expression(paste("Coupon decay rate ",chi,sep="")))
    }
  }
}

dev.off()
