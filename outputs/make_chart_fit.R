
maxH <- 10

res_SS <- make_StateSpace(Model)
M      <- res_SS$M
N      <- res_SS$N
F      <- res_SS$F
dates  <- res_SS$dates

res_KH <- KH_filter(F,M,N,Model$Omega)
fitted <- res_KH$ksi_matrix %*% t(M)

par(mfrow=c(2,3))
for(i in 1:dim(M)[1]){
  plot(F[,i],type="l")
  lines(fitted[,i],col="red")
}

nb_m <- dim(Model$Omega)[1] # number of macro regimes


RES <- prepare_returns_yds(Model,maxH)

plot(c(RES$avg_nominal_yds),type="l",ylim=c(0,.08))
lines(c(RES$avg_ILB_yds),col="red")

plot(c(RES$avg_nominal_returns),type="l")
lines(c(RES$avg_ILB_returns),col="red")
lines(c(RES$avg_GDPLB_returns),col="black",lty=2)



# ==============================================================================
# Prepare figure illustrating the model fit
# ==============================================================================

fitted_yds_nom  <- res_KH$ksi_matrix %*% RES$res_LTnominal_prices$all_LT_rth
fitted_yds_real <- res_KH$ksi_matrix %*% RES$res_LTILB_prices$all_LT_rth

FILE = paste("figures/Figure_fit.pdf",sep="")
pdf(file=FILE, pointsize=10, width=6, height=6)

par(mfrow=c(2,2))
par(plt=c(.1,.95,.1,.85))

plot(dates,F[,1],type="l",lwd=2,main="(a) Inflation",
     xlab="",ylab="")
lines(dates,fitted[,1],lty=2,lwd=2,col="dark grey")

legend("topright", # places a legend at the appropriate place c("Health","Defense"), # puts text in the legend
       c("Data","Model"),
       lty=c(1,2), # gives the legend appropriate symbols (lines)
       lwd=c(2,2), # line width
       col=c("black","dark grey"),
       bg="white",
       #pch=c(3,NaN,NaN),
       seg.len = 2,
       cex=1)

plot(dates,F[,2],type="l",lwd=2,main="(b) GDP growth",
     xlab="",ylab="")
lines(dates,fitted[,2],lty=2,lwd=2,col="dark grey")

plot(DATA$date,DATA$SVENY01/100,type="l",lwd=2,main="(c) 1-yr nominal rate",
     xlab="",ylab="")
lines(dates,fitted_yds_nom[,1],lty=2,lwd=2,col="dark grey")

plot(DATA$date,DATA$SVENY10/100,type="l",lwd=2,main="(d) 10-yr nominal rate",
     xlab="",ylab="")
lines(dates,fitted_yds_nom[,10],lty=2,lwd=2,col="dark grey")

dev.off()


# ==============================================================================
# Prepare figure showing average yield curves
# ==============================================================================

FILE = paste("figures/Figure_avg_yc.pdf",sep="")
pdf(file=FILE, pointsize=10, width=6, height=6)

par(mfrow=c(1,1))
par(plt=c(.1,.95,.1,.85))

plot(1:maxH,rep(0,maxH),
     xlab="",ylab="",type="l",ylim=c(-.01,.10),lwd=2,col="white")

lower_bound <- RES$avg_nominal_yds - 1.0 * RES$Std_nominal_yds
upper_bound <- RES$avg_nominal_yds + 1.0 * RES$Std_nominal_yds
polygon(c(1:maxH,rev(1:maxH)),c(lower_bound,rev(upper_bound)),border = NaN,
        col="#66AAAA44")
lines(1:maxH,RES$avg_nominal_yds,col="blue",lwd=2)

lower_bound <- RES$avg_ILB_yds - 1.0 * RES$Std_ILB_yds
upper_bound <- RES$avg_ILB_yds + 1.0 * RES$Std_ILB_yds
polygon(c(1:maxH,rev(1:maxH)),c(lower_bound,rev(upper_bound)),border = NaN,
        col="#AA66AA44")
lines(1:maxH,RES$avg_ILB_yds,col="red",lwd=2)

grid()

dev.off()


# ==============================================================================
# Prepare figure showing expected (nominal) annual returns
# ==============================================================================

FILE = paste("figures/Figure_expected_returns.pdf",sep="")
pdf(file=FILE, pointsize=10, width=6, height=6)

par(mfrow=c(1,1))
par(plt=c(.1,.95,.1,.85))

plot(1:maxH,rep(0,maxH),
     xlab="",ylab="",type="l",ylim=c(0.02,.12),lwd=2,col="white")

lower_bound <- RES$avg_annual_nominal_returns - 1.0 * RES$Std_annual_nominal_returns
upper_bound <- RES$avg_annual_nominal_returns + 1.0 * RES$Std_annual_nominal_returns
polygon(c(1:maxH,rev(1:maxH)),c(lower_bound,rev(upper_bound)),border = NaN,
        col="#66AAAA44")
lines(1:maxH,c(RES$avg_annual_nominal_returns),col="blue",lwd=2)

# lower_bound <- avg_real_yds - 1.0 * Std_real_yds
# upper_bound <- avg_real_yds + 1.0 * Std_real_yds
# polygon(c(1:maxH,rev(1:maxH)),c(lower_bound,rev(upper_bound)),border = NaN,
#         col="#AA66AA44")
lines(1:maxH,RES$avg_annual_ILB_returns,col="red",lwd=2)

lower_bound <- RES$avg_annual_GDPLB_returns - 1.0 * RES$Std_annual_GDPLB_returns
upper_bound <- RES$avg_annual_GDPLB_returns + 1.0 * RES$Std_annual_GDPLB_returns
polygon(c(1:maxH,rev(1:maxH)),c(lower_bound,rev(upper_bound)),border = NaN,
        col="#AAAAAA44")
lines(1:maxH,RES$avg_annual_GDPLB_returns,col="black",lwd=2,lty=2)

grid()

dev.off()





