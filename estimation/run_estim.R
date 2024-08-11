

load(file="Data/data.Rda")

start_year <- 1970

# Resize dataset:
indic_fst <- which(format(DATA$date,"%Y")==start_year)
DATA <- DATA[indic_fst:dim(DATA)[1],]

# Number of regimes:
nb_m <- 5

# Optimization set up:
maxit.nlminb <- 100
maxit.NlMd   <- 5000
nb_loops     <- 2

# Standard deviations of measurement errors (for Hamilton filter):
std_Pi <- .01
std_Dy <- .01
std_nom_yd <- .01

min_Pi <- -.02
max_Pi <- +.20
min_Dy <- -.10
max_Dy <- +.06
min_gamma <- 1
max_gamma <- 10

max_abs_param <- 8

# Targetted moments:
targets <- list(
  target_10_nom = mean(DATA$SVENY10,na.rm=TRUE)/100,
  target_01_nom = mean(DATA$SVENY01,na.rm=TRUE)/100,
  target_10_rea = mean(DATA$TIPSY10,na.rm=TRUE)/100,
  target_02_rea = mean(DATA$TIPSY02,na.rm=TRUE)/100,
  target_slop_nom = mean(DATA$SVENY10,na.rm=TRUE)/100 - 
    mean(DATA$SVENY01,na.rm=TRUE)/100,
  target_slop_rea = mean(DATA$TIPSY10,na.rm=TRUE)/100 -
    mean(DATA$TIPSY02,na.rm=TRUE)/100,
  target_avg_Pi = mean(DATA$pi,na.rm=TRUE),
  target_avg_Dy = mean(DATA$dy,na.rm=TRUE),
  target_std_10_nom = sd(DATA$SVENY10,na.rm=TRUE)/100,
  target_std_10_rea = sd(DATA$TIPSY10,na.rm=TRUE)/100)

# targets <- list(
#   target_10_nom = .06,
#   target_01_nom = .045,
#   target_10_rea = .015,
#   target_01_rea = .01,
#   target_slop_nom = .015,
#   target_slop_rea = .005,
#   target_avg_Pi = .03,
#   target_avg_Dy = .015)

# Set initial Model ------------------------------------------------------------
gamma     <- 5       # risk aversion
delta     <- .99     # preference for present
beta      <- .3      # sensitivity of surplus to debt level
d_star    <- .6      # targeted debt level
sigma_eps <- .04     # std dev of surplus shocks
alpha     <- 1     # elasticity of default proba
s_star    <- .05     # max surplus
nu_y      <- -.2 * 0
nu_pi     <- .05 * 0

# Perpetuity specification:
kappa_pi <- 0 # Inflation indexation (1 = TIPS)
kappa_y  <- 0 # Real GDP indexation
RR       <- .3 # Recovery rate
chi      <- .7 # coupon decay

rho <- .5
Omega <- rho * diag(nb_m)
Omega[1:(nb_m-1),2:nb_m] <- Omega[1:(nb_m-1),2:nb_m] + (1 - rho)/2 * diag(nb_m-1)
Omega[2:nb_m,1:(nb_m-1)] <- Omega[2:nb_m,1:(nb_m-1)] + (1 - rho)/2 * diag(nb_m-1)
Omega[1,2] <- 1 - Omega[1,1]
Omega[nb_m,nb_m-1] <- 1 - Omega[nb_m,nb_m]
Omega <- Omega + .01 * matrix(1,nb_m,nb_m)
Omega <- Omega / (matrix(apply(Omega,1,sum),ncol=1) %*% matrix(1,1,nb_m))

mu_y  <- matrix(c(min(DATA$dy,na.rm=TRUE),rep(mean(DATA$dy,na.rm=TRUE),nb_m-2),.05),nb_m,1)
mu_pi <- matrix(c(rep(mean(DATA$pi,na.rm=TRUE),nb_m-1),min(DATA$pi,na.rm = TRUE)),nb_m,1)

mu_eta <- .5 * mu_y

Model <- list(gamma=1.2,
              delta     = delta,
              chi       = chi,
              beta      = beta,
              d_star    = d_star,
              sigma_eps = sigma_eps,
              alpha     = alpha,
              s_star    = s_star,
              RR        = RR,
              Omega     = Omega,
              mu_pi     = mu_pi,
              mu_y      = mu_y,
              mu_eta    = mu_eta,
              nu_pi     = nu_pi,
              nu_y      = nu_y,
              kappa_pi  = kappa_pi,
              kappa_y   = kappa_y)

param <- model2param(Model)
param <- param + .6*rnorm(param)

#param <- param + .1*rnorm(param)

algo <- "Nelder-Mead"
#algo <- "nlminb"
for(ii in 1:nb_loops){
  for(algo in c("Nelder-Mead","nlminb")){
    res.estim <- optimx(param,
                        #compute_distance,
                        #compute_loglik,
                        compute_total_distance,
                        targets = targets,
                        method = algo,
                        control=
                          list(trace=1,
                               maxit=ifelse(algo=="Nelder-Mead",
                                            maxit.NlMd,maxit.nlminb),
                               kkt = FALSE))
    param     <- c(as.matrix(res.estim)[1:length(param)])
  }
}
#save(Model,file="results/res_772024.Rdat")
#save(Model,file="results/res_882024.Rdat")
#save(Model,file="results/res_11082024.Rdat")

compute_total_distance(param,targets)

Model <- make_model(param)

# Compute (LT) nominal bond prices:
Model_nominal <- Model
Model_nominal$kappa_pi <- 0
Model_nominal$kappa_y  <- 0
res_LTnominal_prices <- compute_LTRF_bond_prices(Model_nominal,maxH=10)

# Compute (LT) real bond prices:
Model_real <- Model
Model_real$kappa_pi <- 1
Model_real$kappa_y  <- 0
res_LTreal_prices <- compute_LTRF_bond_prices(Model_real,maxH=10)

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
avg_real_yds <- t(stat_distri) %*% res_LTreal_prices$all_LT_rth

plot(c(avg_nom_yds),type="l",ylim=c(0,.08))
lines(c(avg_real_yds),col="red")


# est_yds_real    <- res_KH$ksi_matrix %*% res_LTreal_prices$all_LT_rth
# est_yds_nominal <- res_KH$ksi_matrix %*% res_LTnominal_prices$all_LT_rth
# plot(est_yds_real[,2],type="l",ylim=c(-.02,.04))
# lines(DATA$TIPSY02/100,lty=2)
# plot(est_yds_real[,10],type="l",ylim=c(-.02,.04))
# lines(DATA$TIPSY10/100,lty=2)

# T <- 30
# res_simul <- simul_RS(Model,maxH=10,T)
# par(mfrow=c(2,1))
# plot(res_simul$Pi,type="l")
# plot(res_simul$Dy,type="l")
# 
# # Check KH filter formula:
# noise_Pi <- 3*c(.01,.02,.01)
# noise_Dy <- 3*c(.01,.01,.02)
# N <- rbind(noise_Pi,
#            noise_Dy)
# M <- rbind(Model$mu_pi,Model$mu_y)
# noisy_Pi <- res_simul$Pi + res_simul$states %*% noise_Pi * rnorm(T)
# noisy_Dy <- res_simul$Dy + res_simul$states %*% noise_Dy * rnorm(T)
# plot(noisy_Pi,type="l")
# plot(noisy_Dy,type="l")
# 
# F <- cbind(noisy_Pi,noisy_Dy)
# 
# res_KH <- KH_filter(F,M,N,Model$Omega)
# 
# par(mfrow = c(1,nb_m))
# for(i in 1:nb_m){
#   plot(res_simul$states[,i])
#   lines(res_KH$ksi_matrix[,i],col="red")
# }
# 


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

