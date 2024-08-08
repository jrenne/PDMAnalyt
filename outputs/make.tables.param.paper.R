

threshold.corner <- log(.999/(1-.999)) # threshold to detemrine whether a parameter is

nb.ctries <- length(list.of.ctries)

name.of.file <- c("param_table")


# ------------------------------------------------------------------------------
# Main table ===================================================================
# ------------------------------------------------------------------------------

parameters.2.show <- c("alpha","d_star","gamma","delta",
                       "Phi_11","Phi_22","Phi_33","Phi_44",
                       "sigma_y_1","sigma_y_5",
                       "sigma_pi_2","sigma_pi_6",
                       "sigma_s_3","sigma_s_7",
                       "sigma_star_4",
                       "mu_y","mu_pi","mu_star")

all.param.values       <- NULL
all.param.upper.corner <- NULL
all.param.lower.corner <- NULL
all.param.stdv         <- NULL

counter.ctry <- 0

for(ctry in list.of.ctries){
  counter.ctry <- counter.ctry + 1
  
  print(paste("===========",ctry,"==========="))
  
  path.results <- paste("results/res_estim_",ctry,suffix,".Rdat",sep="")
  load(file=path.results)
  
  param.and.stuff <- show.param(full_theta_est,
                                model_ini_sol,
                                parameters.2.show,
                                bounds,
                                threshold.corner,
                                names.4.latex)
  
  if(indic.comput.stdv.param==1){
    res.stdv <- compute.stdv.param.2.show.in.table(full_theta_est,
                                                   max_abs_val_param4CovMat = threshold.corner,
                                                   List.param.estim = list.param.estim[[1]],
                                                   parameters.2.show = parameters.2.show,
                                                   model_ini_sol,
                                                   bounds,
                                                   list.stdv,
                                                   DATA,
                                                   threshold.corner,
                                                   names.4.latex)
    # Save results (can be used in the Hamilton 1986 approach):
    save(res.stdv,
         file=paste("results/res_MatCov_",ctry,suffix,".Rdat",sep=""))
    
    stdv <- res.stdv$stdv.param
  }else{
    stdv <- rep(0,length(param.and.stuff$param.values))
  }

  all.param.values       <- cbind(all.param.values,param.and.stuff$param.values)
  all.param.lower.corner <- cbind(all.param.lower.corner,
                                  param.and.stuff$param.lower.corner)
  all.param.upper.corner <- cbind(all.param.upper.corner,
                                  param.and.stuff$param.upper.corner)
  all.param.stdv         <- cbind(all.param.stdv,stdv)
  
}


# Parameters that are the same for all countries:
multipl.factor <- param.and.stuff$multipl.factor
param.notation <- param.and.stuff$param.notation
param.constr   <- param.and.stuff$param.constr

# ====  Change param.constr for non-binding contraints:
# param.constr["alpha"] <- 
#   paste("$\\le ",toString(round(bounds["alpha",2],3)),"$",sep="")
# param.constr["mu_star"] <- 
#   paste("$\\ge ",toString(round(100*bounds["mu_star",1],1)),"\\%$",sep="")
param.constr["alpha"]   <- ""
param.constr["mu_star"] <- ""


# remove stdv for Chi:
all.param.stdv[param.notation=="\\chi",] <- 0


# ------------------------------------------------------------------------------
# Prepare Latex table
# ------------------------------------------------------------------------------

latex.table <- rbind(
  "\\begin{table}[H]",
  "\\caption{Models' parameterization}",
  "\\label{tab:param}",
  "\\begin{footnotesize}",
  "\\begin{tabular*}{\\textwidth}{l@{\\extracolsep{\\fill}}crrrrrrrrrrrr}",
  "")

first.line <- " & Constr. & Mult."
for(ctry in list.of.ctries){
  
  path.results <- paste("results/res_estim_",ctry,".Rdat",sep="")
  
  load(file = path.results)
  
  first.line <- paste(first.line,"&",
                      DATA$name.of.country$short.name,
                      sep="")
}
first.line <- paste(first.line,"\\\\",sep="")

latex.table <- rbind(latex.table,
                     "\\hline","\\hline",
                     first.line,
                     "\\hline")

for(i in 1:length(param.notation)){# Loop on table's lines 
  
  if((sum(all.param.lower.corner[i,])!=0)|(sum(all.param.upper.corner[i,])!=0)|
     !str_detect(param.constr[i],'\\[')|
     !str_detect(param.constr[i],'\\]')){
    constr.2.show <- param.constr[i]
  }else{
    constr.2.show <- ""
  }
  
  this.line <- paste(
    "$",param.notation[i],"$&",constr.2.show,"&",
    ifelse(multipl.factor[i,1]!=0,
           paste("$\\times 10^{",toString(multipl.factor[i,1]),"}$",sep=""),
           ""),sep="")
  for(j in 1:nb.ctries){
    if(all.param.stdv[i,j]==0){
      this.line <-  paste(this.line,"&$",
                          paste(ifelse(all.param.values[i,j]==0,
                                       "-",
                                       round.fixed.length(all.param.values[i,j]*10^multipl.factor[i,1],
                                                          multipl.factor[i,2])),sep=""),
                          ifelse(all.param.lower.corner[i,j]==1,"^\\dagger",
                                 ifelse(all.param.upper.corner[i,j]==1,"^\\ddagger","")),"$",sep="")
    }else{
      this.line <-  paste(this.line,
                          "&$\\underset{(",
                          round.fixed.length(all.param.stdv[i,j]*10^multipl.factor[i,1],
                                             multipl.factor[i,2]),")}{",
                          round.fixed.length(all.param.values[i,j]*10^multipl.factor[i,1],
                                             multipl.factor[i,2]),
                          "}$",sep="")
    }
  }
  this.line <- paste(this.line,"\\\\")
  if(multipl.factor[i,3]>0){
    latex.table <- rbind(
      latex.table,
      ifelse(multipl.factor[i,3]==2,"\\hline",""),
      this.line
    )
  }
  
}

latex.table <- rbind(latex.table,
                     "\\hline",
                     "\\hline",
                     "\\end{tabular*}",
                     "\\end{footnotesize}",
                     "\\vspace{.1cm}",
                     "\\begin{footnotesize}",
                     "\\parbox{\\linewidth}{",
                     "Note: This table presents the models' parameterizations. Asymptotic standard deviations (based on the outer product of the log-likelihood gradient) are reported in parentheses. (These standard deviations are calculated for parameters that are sufficiently far from the imposed bounds, if any.)
  The second column of the table indicates the constraints that have been imposed during the 
  maximizaiton of the likelihood function. $\\dagger$ (respectively $\\ddagger$) indicates that the lower (resp. upper) bound is binding. ``Data'' means the parameters is calibrated from sample moments; for $d^*$ (see Eq.\\,\\ref{eq:sd}), ``Data $\\pm 10\\,p.p.$'' means that the parameter cannot be smaller (larger) than the sample mean minus (plus) 10 percentage points; ``pre-covid sd'' is the standard deviation of the budget surplus, as measured over the pre-covid period. We impose $\\sigma^*_{4}=\\sigma^*_{8}$. 
  $H$ is the duration of the perpetuity, expressed in years (see supplementary meterial for sources). $\\chi$ is the decay rate of perpetuities' coupons. 
  $RR$ denotes the recovery rate. $\\gamma$ stands for the coefficient of relative risk aversion (Eq.\\,\\ref{eq:sdf}). $b_y$ is the output fall upon default (Eq.\\,\\ref{eq:cypis_affine}).
  $b_\\pi$ is the inflation increase upon default (Eq.\\,\\ref{eq:cypis_affine}). $\\alpha$ is the elasticity of the probability of default to the surplus gap ($s_t-s_{t}^*$, see Eq.\\,\\ref{eq:probadef}).
  $\\mu_i$ and $\\sigma_i$ determine the specification of the macroeconomic variables, with $i=y,\\,\\pi,\\,s,\\,*$ (Eqs.\\,\\ref{eq:cypis_affine}, \\ref{eq:sd}, and \\ref{eq:s_star}), $\\mu_y$ and $\\mu_\\pi$ are quarterly growth rates.
  In the state-space model, the standard deviations of the measurement errors associated with yields and CDS spreads are respectively set to 20\\% of the sample standard deviations of yields and CDS spreads.
  The standard deviations of the measurement errors associated with cumulated growth, cumulated inflation, and debt forecasts are increasing with the horizon ($h=$ 1, 2, 3, 4, and 5 years); for growth and inflation and horizon $h$, the standard deviation is set to $0.2 \\times h \\times sd(x)/5$, where $sd(x)$ denotes the sample standard deviation of growth and inflation; for debt forecasts, it is set to $0.5 \\times h \\times sd(x)/5$.}",
                     "\\end{footnotesize}",
                     "\\end{table}")


vec.dates <- DATA$DATES

latex.file <- paste(name.of.file,suffix,".txt", sep="")
write(latex.table, paste("tables/",latex.file,sep=""))





# ------------------------------------------------------------------------------
# Sigma table ==================================================================
# ------------------------------------------------------------------------------

all.param.values <- NULL

counter.ctry <- 0

for(ctry in list.of.ctries){
  counter.ctry <- counter.ctry + 1
  
  path.results <- paste("results/res_estim_",ctry,suffix,".Rdat",sep="")
  load(file=path.results)
  model.est <- Theta2Model(full_theta_est,model_ini_sol,bounds)
  model.sol.est <- solve_model(model.est,indic_compute_beta = 1)
  Sigma_persist <- model.sol.est$Sigma.w[1:4,1:4]
  Sigma_volatil <- model.sol.est$Sigma.w[5:8,5:8]
  param.values     <- c(Sigma_persist[!upper.tri(Sigma_persist)],
                        Sigma_volatil[!upper.tri(Sigma_volatil)])
  all.param.values <- cbind(all.param.values,param.values)
}

indic.param <- apply(matrix(parameters.2.show,ncol=1),1,
                     function(x){which(x==rownames(bounds))})

param.notation <- NULL
for(i in 1:4){
  for(j in 1:i){
    param.notation <- c(param.notation,
                        paste("\\Sigma_{w,",i,",",j,"}",sep=""))
  }
}
for(i in 5:8){
  for(j in 5:i){
    param.notation <- c(param.notation,
                        paste("\\Sigma_{w,",i,",",j,"}",sep=""))
  }
}

multipl.factor <- matrix(0,length(param.notation),3)
multipl.factor[,2] <- 3
multipl.factor[,3] <- 1
multipl.factor[11,3] <- 2 # put a hline there

# ------------------------------------------------------------------------------
# Prepare Latex table
# ------------------------------------------------------------------------------


latex.table <- rbind(
  "\\begin{table}[H]",
  "\\caption{Models' parameterization ($\\Sigma_w$)}",
  "\\label{tab:paramSigma}",
  "\\begin{footnotesize}",
  "\\begin{tabular*}{\\textwidth}{l@{\\extracolsep{\\fill}}rrrrrrrrrrr}",
  "")

first.line <- ""

for(ctry in list.of.ctries){
  
  path.results <- paste("results/res_estim_",ctry,".Rdat",sep="")
  
  load(file = path.results)
  
  first.line <- paste(first.line,"&",
                      DATA$name.of.country$short.name,
                      sep="")
}
first.line <- paste(first.line,"\\\\",sep="")

latex.table <- rbind(latex.table,
                     "\\hline","\\hline",
                     first.line,
                     "\\hline")


for(i in 1:length(param.notation)){# Loop on table's lines 
  this.line <- paste("$",param.notation[i],"$",sep="")
  for(j in 1:nb.ctries){
    this.line <-  paste(this.line,"&$",
                        paste(ifelse(all.param.values[i,j]==0,
                                     "-",
                                     round.fixed.length(all.param.values[i,j]*10^multipl.factor[i,1],
                                                        multipl.factor[i,2])),sep=""),"$",sep="")
  }
  this.line <- paste(this.line,"\\\\")
  if(multipl.factor[i,3]>0){
    latex.table <- rbind(
      latex.table,
      ifelse(multipl.factor[i,3]==2,"\\hline",""),
      this.line
    )
  }
  
}

latex.table <- rbind(latex.table,
                     "\\hline",
                     "\\hline",
                     "\\end{tabular*}",
                     "\\end{footnotesize}",
                     "\\vspace{.1cm}",
                     "\\begin{footnotesize}",
                     "\\parbox{\\linewidth}{",
                     "Note: This table reports the estimated parameterization of $\\Sigma_w$. Given Eq. \\eqref{eq:w}, we have that $\\Sigma_w\\Sigma_w'$ is the conditional covariance matrix of $w_{t+1}$ (as of date $t$). ",
                     "This matrix is block diagonal. The $4 \\times 4$ upper-left block (respectively lower-right block) is lower triangular and corresponds to the persistent components (resp. volatile components) of $w_t$, its specification is given in the upper part of the table (resp. in the lower part of the table). ",
                     "The parameterization is such that, for the sake of identification, the unconditional variance of each of the $w_{i,t}$'s is equal to one.}",
                     "\\end{footnotesize}",
                     "\\end{table}")

vec.dates <- DATA$DATES

name.of.file <- "param_table_Sigma"

latex.file <- paste(name.of.file,suffix,".txt", sep="")
write(latex.table, paste("tables/",latex.file,sep=""))

