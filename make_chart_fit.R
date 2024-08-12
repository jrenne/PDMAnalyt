# Compute (LT) nominal bond prices:
Model_nominal <- Model
Model_nominal$kappa_pi <- 0
Model_nominal$kappa_y  <- 0
res_LTnominal_prices <- compute_LTRF_bond_prices(Model_nominal,maxH=maxH)

# Compute (LT) real bond prices:
Model_TIPS <- Model
Model_TIPS$kappa_pi <- 1
Model_TIPS$kappa_y  <- 0
res_LTTIPS_prices <- compute_LTRF_bond_prices(Model_TIPS,maxH=maxH)

# Compute (LT) TIPS bond prices:
Model_GDPLB <- Model
Model_GDPLB$kappa_pi <- 1
Model_GDPLB$kappa_y  <- 1
res_LTGDPLB_prices   <- compute_LTRF_bond_prices(Model_GDPLB,maxH=maxH)


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

stat_distri <- compute_stat_distri(Model)
avg_nom_yds  <- t(stat_distri) %*% res_LTnominal_prices$all_LT_rth
avg_TIPS_yds <- t(stat_distri) %*% res_LTTIPS_prices$all_LT_rth

avg_nom_returns    <- t(stat_distri) %*% res_LTnominal_prices$all_LT_ExpReturn_th
avg_TIPS_returns   <- t(stat_distri) %*% res_LTTIPS_prices$all_LT_ExpReturn_th
avg_GDPLB_returns  <- t(stat_distri) %*% res_LTGDPLB_prices$all_LT_ExpReturn_th

plot(c(avg_nom_yds),type="l",ylim=c(0,.08))
lines(c(avg_TIPS_yds),col="red")

plot(c(avg_nom_returns),type="l")
lines(c(avg_TIPS_returns),col="red")
lines(c(avg_GDPLB_returns),col="black",lty=2)


# ==============================================================================
# Prepare table showing model parameterization
# ==============================================================================

make.entry <- function(x,format.nb){
  output <- paste("$",sprintf(format.nb,x),"$",sep="")
  return(output)
}

format.nb0 <- paste("%.",0,"f",sep="")
format.nb1 <- paste("%.",1,"f",sep="")
format.nb2 <- paste("%.",2,"f",sep="")
format.nb3 <- paste("%.",3,"f",sep="")


columns <- NULL
for(i in 1:5){
  columns <- paste(columns,"r",sep="")
}
for(i in 1:nb_m){
  columns <- paste(columns,"c",sep="")
}

latex_table <- rbind("\\begin{table}[ph!]",
                     "\\caption{Model parameterization}",
                     "\\label{tab:param}",
                     paste("\\begin{tabular*}{\\textwidth}{c@{\\extracolsep{\\fill}}",columns,"}",sep=""),
                     "\\hline",
                     paste("Regime&$\\mu_\\pi$&&$\\mu_y$&&\\multicolumn{",nb_m,"}{c}{$\\Omega$}\\\\",sep=""),
                     "\\hline")

for(i in 1:nb_m){
  pi_i <- NULL
  for(j in 1:nb_m){
    pi_i <- paste(pi_i,"&",make.entry(Model$Omega[i,j],format.nb3),sep="")
  }
  this_line <- paste(i,"&",make.entry(Model$mu_pi[i],format.nb3),
                     "&&",make.entry(Model$mu_y[i],format.nb3),
                     "&&",pi_i,
                     "\\\\",sep="")
  latex_table <- rbind(latex_table,
                       this_line)
}

latex_table <- rbind(latex_table,
                     "\\hline",
                     "\\end{tabular*}",
                     "\\begin{footnotesize}",
                     "\\parbox{\\linewidth}{\\textit{Notes}: This table shows the model parameterization of the baseline model.}",
                     "\\end{footnotesize}",
                     "\\end{table}")

name.of.file <- "table_param"
latex.file <- paste(name.of.file,".txt", sep="")
write(latex_table, paste("tables/",latex.file,sep=""))



# ==============================================================================
# Prepare figure illustrating the model fit
# ==============================================================================

fitted_yds_nom  <- res_KH$ksi_matrix %*% res_LTnominal_prices$all_LT_rth
fitted_yds_real <- res_KH$ksi_matrix %*% res_LTreal_prices$all_LT_rth

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

Var_nom_yds  <- t(stat_distri) %*% res_LTnominal_prices$all_LT_rth^2 - avg_nom_yds^2
Var_real_yds <- t(stat_distri) %*% res_LTreal_prices$all_LT_rth^2    - avg_real_yds^2

Std_nom_yds  <- sqrt(Var_nom_yds)
Std_real_yds <- sqrt(Var_real_yds)

max_mat <- dim(res_LTnominal_prices$all_LT_rth)[2]

FILE = paste("figures/Figure_avg_yc.pdf",sep="")
pdf(file=FILE, pointsize=10, width=6, height=6)

par(mfrow=c(1,1))
par(plt=c(.1,.95,.1,.85))

plot(1:max_mat,rep(0,max_mat),
     xlab="",ylab="",type="l",ylim=c(-.01,.10),lwd=2,col="white")

lower_bound <- avg_nom_yds - 1.0 * Std_nom_yds
upper_bound <- avg_nom_yds + 1.0 * Std_nom_yds
polygon(c(1:max_mat,rev(1:max_mat)),c(lower_bound,rev(upper_bound)),border = NaN,
        col="#66AAAA44")
lines(1:max_mat,avg_nom_yds,col="blue",lwd=2)

lower_bound <- avg_real_yds - 1.0 * Std_real_yds
upper_bound <- avg_real_yds + 1.0 * Std_real_yds
polygon(c(1:max_mat,rev(1:max_mat)),c(lower_bound,rev(upper_bound)),border = NaN,
        col="#AA66AA44")
lines(1:max_mat,avg_real_yds,col="red",lwd=2)

grid()

dev.off()


# ==============================================================================
# Prepare figure showing expected (nominal) annual returns
# ==============================================================================

exp_nominal_returns <- log(res_LTnominal_prices$all_LT_ExpReturn_th)/t(matrix(1:maxH,maxH,nb_m))
exp_TIPS_returns    <- log(res_LTTIPS_prices$all_LT_ExpReturn_th)/t(matrix(1:maxH,maxH,nb_m))
exp_GDPLB_returns   <- log(res_LTGDPLB_prices$all_LT_ExpReturn_th)/t(matrix(1:maxH,maxH,nb_m))

avg_nominal_returns <- c(t(stat_distri) %*% exp_nominal_returns)
avg_TIPS_returns    <- c(t(stat_distri) %*% exp_TIPS_returns)
avg_GDPLB_returns   <- c(t(stat_distri) %*% exp_GDPLB_returns)

Var_nominal_returns <- t(stat_distri) %*% exp_nominal_returns^2 - avg_nominal_returns^2
Var_TIPS_returns    <- t(stat_distri) %*% exp_TIPS_returns^2    - avg_TIPS_returns^2
Var_GDPLB_returns   <- t(stat_distri) %*% exp_GDPLB_returns^2   - avg_GDPLB_returns^2

Std_nominal_returns <- sqrt(Var_nominal_returns)
Std_TIPS_returns    <- sqrt(Var_TIPS_returns)
Std_GDPLB_returns   <- sqrt(Var_GDPLB_returns)

FILE = paste("figures/Figure_expected_returns.pdf",sep="")
pdf(file=FILE, pointsize=10, width=6, height=6)

par(mfrow=c(1,1))
par(plt=c(.1,.95,.1,.85))

plot(1:max_mat,rep(0,max_mat),
     xlab="",ylab="",type="l",ylim=c(-.01,.10),lwd=2,col="white")

lower_bound <- avg_nominal_returns - 1.0 * Std_nominal_returns
upper_bound <- avg_nominal_returns + 1.0 * Std_nominal_returns
polygon(c(1:max_mat,rev(1:max_mat)),c(lower_bound,rev(upper_bound)),border = NaN,
        col="#66AAAA44")
lines(1:max_mat,avg_nominal_returns,col="blue",lwd=2)

# lower_bound <- avg_real_yds - 1.0 * Std_real_yds
# upper_bound <- avg_real_yds + 1.0 * Std_real_yds
# polygon(c(1:max_mat,rev(1:max_mat)),c(lower_bound,rev(upper_bound)),border = NaN,
#         col="#AA66AA44")
lines(1:max_mat,avg_TIPS_returns,col="red",lwd=2)

lower_bound <- avg_GDPLB_returns - 1.0 * Std_GDPLB_returns
upper_bound <- avg_GDPLB_returns + 1.0 * Std_GDPLB_returns
polygon(c(1:max_mat,rev(1:max_mat)),c(lower_bound,rev(upper_bound)),border = NaN,
        col="#AAAAAA44")
lines(1:max_mat,avg_GDPLB_returns,col="black",lwd=2,lty=2)

grid()

dev.off()





