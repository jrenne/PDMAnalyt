
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

RES <- prepare_returns_yds(Model,maxH)

plot(c(RES$avg_nominal_yds),type="l",ylim=c(0,.08))
lines(c(RES$avg_ILB_yds),col="red")

plot(c(RES$avg_nominal_returns),type="l")
lines(c(RES$avg_ILB_returns),col="red")
lines(c(RES$avg_GDPLB_returns),col="black",lty=2)

