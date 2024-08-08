#################################################################################################
library(stargazer)

dpath <- getwd()
save.path <- paste(dpath,"/tables", sep="")



PD <- read.csv(paste(dpath,"/data/panel/data/EPU_FS/EPU_FSL.csv", sep=""), dec=".")

#install.packages("lmtest", dependencies = TRUE)
library(lmtest)
# install.packages("plm", dependencies = TRUE)
library(plm)

names(PD)[1] <- "TIME"
PDdf <- pdata.frame(PD, index = c("COUNTRY", "TIME"), drop.index = FALSE, row.names = TRUE,
                    stringsAsFactors = default.stringsAsFactors(),
                    replace.non.finite = FALSE, drop.NA.series = FALSE,
                    drop.const.series = FALSE, drop.unused.levels = FALSE)


# Country FE

Panel.REG.FS <- plm(FS ~ lag(FS,1) + EPU, data = PDdf, index=c("COUNTRY", "TIME"), effect =c("individual"), model = "within")
Panel.REG.FS.w <- coeftest(Panel.REG.FS, vcov = vcovHC(Panel.REG.FS, method = "white2", type = "HC3"))
summary(Panel.REG.FS)

Panel.REG.FL <- plm(FL ~ lag(FL,1) + EPU, data = PDdf, index=c("COUNTRY", "TIME"), effect =c("individual"), model = "within")
Panel.REG.FL.w <- coeftest(Panel.REG.FL, vcov = vcovHC(Panel.REG.FL, method = "white2", type = "HC3"))
summary(Panel.REG.FL)

# Time (quarter) FE

Panel.REG.FS.t <- plm(FS ~ lag(FS,1) + EPU, data = PDdf, index=c("COUNTRY", "TIME"), effect =c("time"), model = "within")
Panel.REG.FS.t.w <- coeftest(Panel.REG.FS.t, vcov = vcovHC(Panel.REG.FS.t, method = "white2", type = "HC3"))
summary(Panel.REG.FS.t)

Panel.REG.FL.t <- plm(FL ~ lag(FL,1) + EPU, data = PDdf, index=c("COUNTRY", "TIME"), effect =c("time"), model = "within")
Panel.REG.FL.t.w <- coeftest(Panel.REG.FL.t, vcov = vcovHC(Panel.REG.FL.t, method = "white2", type = "HC3"))
summary(Panel.REG.FL.t)

# Both Country and Time (quarter) FE

Panel.REG.FS.2 <- plm(FS ~ lag(FS,1) + EPU, data = PDdf, index=c("COUNTRY", "TIME"), effect =c("twoways"), model = "within")
Panel.REG.FS.2.w <- coeftest(Panel.REG.FS.2, vcov = vcovHC(Panel.REG.FS.2, method = "white2", type = "HC3"))
summary(Panel.REG.FS.2)

Panel.REG.FL.2 <- plm(FL ~ lag(FL,1) + EPU, data = PDdf, index=c("COUNTRY", "TIME"), effect =c("twoways"), model = "within")
Panel.REG.FL.2.w <- coeftest(Panel.REG.FL.2, vcov = vcovHC(Panel.REG.FL.2, method = "white2", type = "HC3"))
summary(Panel.REG.FL.2)

# Pooling

Panel.REG.FSp.2 <- plm(FS ~ lag(FS,1) + EPU, data = PDdf, index=c("COUNTRY", "TIME"), effect =c("twoways"), model = "pooling")
Panel.REG.FSp.2.w <- coeftest(Panel.REG.FSp.2, vcov = vcovHC(Panel.REG.FSp.2, method = "white2", type = "HC3"))
summary(Panel.REG.FSp.2)

Panel.REG.FLp.2 <- plm(FL ~ lag(FL,1) + EPU, data = PDdf, index=c("COUNTRY", "TIME"), effect =c("twoways"), model = "pooling")
Panel.REG.FLp.2.w <- coeftest(Panel.REG.FLp.2, vcov = vcovHC(Panel.REG.FLp.2, method = "white2", type = "HC3"))
summary(Panel.REG.FLp.2)

####################################################################################

####################################################################################

# AE cluster #

####################################################################################
PDae <- read.csv(paste(dpath,"/data/panel/data/EPU_FS/EPU_FSL_advanced.csv", sep=""), dec=".")
names(PDae)[1] <- "TIME"


PDdf <- pdata.frame(PDae, index = c("COUNTRY", "TIME"), drop.index = FALSE, row.names = TRUE,
                    stringsAsFactors = default.stringsAsFactors(),
                    replace.non.finite = FALSE, drop.NA.series = FALSE,
                    drop.const.series = FALSE, drop.unused.levels = FALSE)


# Country FE

Panel.REG.FS.ae <- plm(FS ~ lag(FS,1) + EPU, data = PDdf, index=c("COUNTRY", "TIME"), effect =c("individual"), model = "within")
Panel.REG.FS.ae.w <- coeftest(Panel.REG.FS.ae, vcov = vcovHC(Panel.REG.FS.ae, method = "white2", type = "HC3"))
summary(Panel.REG.FS.ae)

Panel.REG.FL.ae <- plm(FL ~ lag(FL,1) + EPU, data = PDdf, index=c("COUNTRY", "TIME"), effect =c("individual"), model = "within")
Panel.REG.FL.ae.w <- coeftest(Panel.REG.FL.ae, vcov = vcovHC(Panel.REG.FL.ae, method = "white2", type = "HC3"))
summary(Panel.REG.FL.ae)

# Time (quarter) FE

Panel.REG.FS.ae.t <- plm(FS ~ lag(FS,1) + EPU, data = PDdf, index=c("COUNTRY", "TIME"), effect =c("time"), model = "within")
Panel.REG.FS.ae.t.w <- coeftest(Panel.REG.FS.ae.t, vcov = vcovHC(Panel.REG.FS.ae.t, method = "white2", type = "HC3"))
summary(Panel.REG.FS.ae.t)

Panel.REG.FL.ae.t <- plm(FL ~ lag(FL,1) + EPU, data = PDdf, index=c("COUNTRY", "TIME"), effect =c("time"), model = "within")
Panel.REG.FL.ae.t.w <- coeftest(Panel.REG.FL.ae.t, vcov = vcovHC(Panel.REG.FL.ae.t, method = "white2", type = "HC3"))
summary(Panel.REG.FL.ae.t)

# Both Country and Time (quarter) FE

Panel.REG.FS.ae.2 <- plm(FS ~ lag(FS,1) + EPU, data = PDdf, index=c("COUNTRY", "TIME"), effect =c("twoways"), model = "within")
Panel.REG.FS.ae.2.w <- coeftest(Panel.REG.FS.ae.2, vcov = vcovHC(Panel.REG.FS.ae.2, method = "white2", type = "HC3"))
summary(Panel.REG.FS.ae.2)

Panel.REG.FL.ae.2 <- plm(FL ~ lag(FL,1) + EPU, data = PDdf, index=c("COUNTRY", "TIME"), effect =c("twoways"), model = "within")
Panel.REG.FL.ae.2.w <- coeftest(Panel.REG.FL.ae.2, vcov = vcovHC(Panel.REG.FL.ae.2, method = "white2", type = "HC3"))
summary(Panel.REG.FL.ae.2)

# Pooling

Panel.REG.FS.aep.2 <- plm(FS ~ lag(FS,1) + EPU, data = PDdf, index=c("COUNTRY", "TIME"), effect =c("twoways"), model = "pooling")
Panel.REG.FS.aep.2.w <- coeftest(Panel.REG.FS.aep.2, vcov = vcovHC(Panel.REG.FS.aep.2, method = "white2", type = "HC3"))
summary(Panel.REG.FS.aep.2)

Panel.REG.FL.aep.2 <- plm(FL ~ lag(FL,1) + EPU, data = PDdf, index=c("COUNTRY", "TIME"), effect =c("twoways"), model = "pooling")
Panel.REG.FL.aep.2.w <- coeftest(Panel.REG.FL.aep.2, vcov = vcovHC(Panel.REG.FL.aep.2, method = "white2", type = "HC3"))
summary(Panel.REG.FL.aep.2)

####################################################################################################################

# EME cluster

####################################################################################################################

PD.eme <- read.csv(paste(dpath,"/data/panel/data/EPU_FS/EPU_FSL_emerging.csv", sep=""), dec=".")
names(PD.eme)[1] <- "TIME"

PDdf <- pdata.frame(PD.eme, index = c("COUNTRY", "TIME"), drop.index = FALSE, row.names = TRUE,
                    stringsAsFactors = default.stringsAsFactors(),
                    replace.non.finite = FALSE, drop.NA.series = FALSE,
                    drop.const.series = FALSE, drop.unused.levels = FALSE)


# Country FE

Panel.REG.FS.eme <- plm(FS ~ lag(FS,1) + EPU, data = PDdf, index=c("COUNTRY", "TIME"), effect =c("individual"), model = "within")
Panel.REG.FS.eme.w <- coeftest(Panel.REG.FS.eme, vcov = vcovHC(Panel.REG.FS.eme, method = "white2", type = "HC3"))
summary(Panel.REG.FS.eme)

Panel.REG.FL.eme <- plm(FL ~ lag(FL,1) + EPU, data = PDdf, index=c("COUNTRY", "TIME"), effect =c("individual"), model = "within")
Panel.REG.FL.eme.w <- coeftest(Panel.REG.FL.eme, vcov = vcovHC(Panel.REG.FL.eme, method = "white2", type = "HC3"))
summary(Panel.REG.FL.eme)

# Time (quarter) FE

Panel.REG.FS.eme.t <- plm(FS ~ lag(FS,1) + EPU, data = PDdf, index=c("COUNTRY", "TIME"), effect =c("time"), model = "within")
Panel.REG.FS.eme.t.w <- coeftest(Panel.REG.FS.eme.t, vcov = vcovHC(Panel.REG.FS.eme.t, method = "white2", type = "HC3"))
summary(Panel.REG.FS.eme.t)

Panel.REG.FL.eme.t <- plm(FL ~ lag(FL,1) + EPU, data = PDdf, index=c("COUNTRY", "TIME"), effect =c("time"), model = "within")
Panel.REG.FL.eme.t.w <- coeftest(Panel.REG.FL.eme.t, vcov = vcovHC(Panel.REG.FL.eme.t, method = "white2", type = "HC3"))
summary(Panel.REG.FL.eme.t)

# Both Country and Time (quarter) FE

Panel.REG.FS.eme.2 <- plm(FS ~ lag(FS,1) + EPU, data = PDdf, index=c("COUNTRY", "TIME"), effect =c("twoways"), model = "within")
Panel.REG.FS.eme.2.w <- coeftest(Panel.REG.FS.eme.2, vcov = vcovHC(Panel.REG.FS.eme.2, method = "white2", type = "HC3"))
summary(Panel.REG.FS.eme.2)

Panel.REG.FL.eme.2 <- plm(FL ~ lag(FL,1) + EPU, data = PDdf, index=c("COUNTRY", "TIME"), effect =c("twoways"), model = "within")
Panel.REG.FL.eme.2.w <- coeftest(Panel.REG.FL.eme.2, vcov = vcovHC(Panel.REG.FL.eme.2, method = "white2", type = "HC3"))
summary(Panel.REG.FL.eme.2)

# Pooling

Panel.REG.FS.emep.2 <- plm(FS ~ lag(FS,1) + EPU, data = PDdf, index=c("COUNTRY", "TIME"), effect =c("twoways"), model = "pooling")
Panel.REG.FS.emep.2.w <- coeftest(Panel.REG.FS.emep.2, vcov = vcovHC(Panel.REG.FS.emep.2, method = "white2", type = "HC3"))
summary(Panel.REG.FS.emep.2)

Panel.REG.FL.emep.2 <- plm(FL ~ lag(FL,1) + EPU, data = PDdf, index=c("COUNTRY", "TIME"), effect =c("twoways"), model = "pooling")
Panel.REG.FL.emep.2.w <- coeftest(Panel.REG.FL.emep.2, vcov = vcovHC(Panel.REG.FL.emep.2, method = "white2", type = "HC3"))
summary(Panel.REG.FL.emep.2)



####################################################################################

stargazer(Panel.REG.FS.2.w, Panel.REG.FS.w,
          Panel.REG.FS.t.w,
          title="All countries - Panel Regression Results")

stargazer(Panel.REG.FS.ae.2.w, Panel.REG.FS.ae.w,
          Panel.REG.FS.ae.t.w,
          title="Advanced Economies - Panel Regression Results")

stargazer(Panel.REG.FS.eme.2.w, Panel.REG.FS.eme.w,
          Panel.REG.FS.eme.t.w,
          title="Emerging Economies - Panel Regression Results")


####################################################################################

ALLreg <- stargazer(Panel.REG.FS.2.w, Panel.REG.FS.w,
                    Panel.REG.FS.t.w,
                    covariate.labels=c("$FS_{t-1}$","$EPU_{t}$"),
                    title="All countries - Panel Regression Results")

AEreg <- stargazer(Panel.REG.FS.ae.2.w, Panel.REG.FS.ae.w,
                   Panel.REG.FS.ae.t.w,
                   covariate.labels=c("$FS_{t-1}$","$EPU_{t}$"),
                   title="Advanced Economies - Panel Regression Results")

EMEreg <- stargazer(Panel.REG.FS.eme.2.w, Panel.REG.FS.eme.w,
                    Panel.REG.FS.eme.t.w,
                    covariate.labels=c("$FS_{t-1}$","$EPU_{t}$"),
                    title="Emerging Economies - Panel Regression Results")


stargazer(Panel.REG.FS.2.w, Panel.REG.FS.w,
          Panel.REG.FS.t.w,
          Panel.REG.FS.ae.2.w, Panel.REG.FS.ae.w,
          Panel.REG.FS.ae.t.w,
          Panel.REG.FS.eme.2.w, Panel.REG.FS.eme.w,
          Panel.REG.FS.eme.t.w, title="Fiscal space and economic policy uncertainty - Panel regression results",
          align=TRUE, dep.var.labels=c("$FS_{t}$","$FS_{t}$"),
          covariate.labels=c("$FS_{t-1}$","$EPU_{t}$"),
          omit.stat=c("LL","ser","f"), no.space=TRUE)





# --------------------------------
# Prepare Latex table
# --------------------------------

latex.table <- c(
  "\\begin{table}[!htbp] \\centering 
  \\caption{Fiscal space and economic policy uncertainty - Panel regression results} 
   \\label{tab:prFS}  
\\begin{tabular*}{\\textwidth}{l@{\\extracolsep{\\fill}}ccc}
%{@{\\extracolsep{5pt}}lccccccccc} 
\\hline 
\\hline 
\\\\
 & \\multicolumn{3}{c}{{\\bf Panel A} - \\textit{All countries}} \\\\ 
\\cline{2-4} 
& \\multicolumn{3}{c}{ } \\\\ 
& Country \\& Time FE & Country FE & Time FE\\\\ 
& $FS_{t}$ & $FS_{t}$ & $FS_{t}$\\\\
\\hline ")

allreg.lines <- rbind(ALLreg[15],ALLreg[16],ALLreg[18],ALLreg[19])

fill.line.1 <- c("\\hline \\\\
 & \\multicolumn{3}{c}{{\\bf Panel B} - \\textit{Advanced Economies (US, UK, JP, CA)}} \\\\ 
\\cline{2-4} 
& \\multicolumn{3}{c}{ } \\\\ 
& Country \\& Time FE & Country FE & Time FE\\\\ 
& $FS_{t}$ & $FS_{t}$ & $FS_{t}$\\\\
\\hline ") 

aereg.lines <- rbind(AEreg[15],AEreg[16],AEreg[18],AEreg[19])

fill.line.2 <- c("\\hline \\\\
 & \\multicolumn{3}{c}{{\\bf Panel C} - \\textit{Emerging Economies (BR, CN, IN, RU)}} \\\\ 
\\cline{2-4} 
& \\multicolumn{3}{c}{ } \\\\ 
& Country \\& Time FE & Country FE & Time FE\\\\ 
& $FS_{t}$ & $FS_{t}$ & $FS_{t}$\\\\
\\hline ")

emereg.lines <- rbind(EMEreg[15],EMEreg[16],EMEreg[18],EMEreg[19])



latex.table <- rbind(latex.table,
                     allreg.lines,
                     fill.line.1,
                     aereg.lines,
                     fill.line.2,
                     emereg.lines,
                     "\\hline 
                      \\hline 
                      \\end{tabular*} 
                      \\begin{footnotesize}
                      \\parbox{\\linewidth}{
                      Note: This table reports the results of panel regressions of fiscal space ($FS$) estimates on the Economic Policy Uncertainty ($EPU$) indices.  The estimation sample goes from $2004Q1$ to $2022Q3$. See text for more details. FE stands for Fixed Effects. $^{*}$p$<$0.1; $^{**}$p$<$0.05; $^{***}$p$<$0.01}.
                      \\end{footnotesize}
                     \\end{table} ")

name.of.file <- c("PanelRegFS")

latex.file <- paste(name.of.file, ".tex", sep="")
write(latex.table, paste(save.path,"/",latex.file,sep=""))


################################################################################
stop()
# Granger Test - FS and EPU 

first.date <- as.Date("01/04/2004","%d/%m/%Y")
last.date <- as.Date("01/01/2021","%d/%m/%Y")

All_FL_FS <-  read.csv(paste(dpath,"/data/panel/data/all_FL_FS.csv", sep=""), dec=".")
All_FL_FS$date <- as.Date(All_FL_FS$date)
EPU_q <- read.csv(paste(dpath,"/data/panel/data/EPU_FS/EPU_DB_quarterly.csv",sep=""), dec=".")
EPU_q$DATE <- as.Date(EPU_q$DATE)

indic.first.epu <- which(EPU_q$DATE == first.date)
indic.last.epu <-  which(EPU_q$DATE == last.date)
indic.first.fs <-  which(All_FL_FS$date == first.date)
indic.last.fs <- which(All_FL_FS$date == last.date)

All_FL_FS <- All_FL_FS[indic.first.fs:indic.last.fs,]
EPU_q <- EPU_q[indic.first.epu:indic.last.epu,]

library(lmtest)

GT.US.EPU <- grangertest(EPU_q$US~All_FL_FS$FS.US, order = 4)
GT.US.FS <- grangertest(All_FL_FS$FS.US~EPU_q$US, order = 4)

GT.UK.EPU <- grangertest(EPU_q$UK~All_FL_FS$FS.UK, order = 4)
GT.UK.FS <- grangertest(All_FL_FS$FS.UK~EPU_q$UK, order = 4)

GT.JP.EPU <- grangertest(EPU_q$JP~All_FL_FS$FS.JP, order = 4)
GT.JP.FS <- grangertest(All_FL_FS$FS.JP~EPU_q$JP, order = 4)

GT.CA.EPU <- grangertest(EPU_q$CA~All_FL_FS$FS.CA, order = 4)
GT.CA.FS <- grangertest(All_FL_FS$FS.CA~EPU_q$CA, order = 4)

GT.BR.EPU <- grangertest(EPU_q$BR~All_FL_FS$FS.BR, order = 4)
GT.BR.FS <- grangertest(All_FL_FS$FS.BR~EPU_q$BR, order = 4)

GT.CN.EPU <- grangertest(EPU_q$CN~All_FL_FS$FS.CN, order = 4)
GT.CN.FS <- grangertest(All_FL_FS$FS.CN~EPU_q$CN, order = 4)

GT.IN.EPU <- grangertest(EPU_q$IN~All_FL_FS$FS.IN, order = 4)
GT.IN.FS <- grangertest(All_FL_FS$FS.IN~EPU_q$IN, order = 4)

GT.RU.EPU <- grangertest(EPU_q$RU~All_FL_FS$FS.RU, order = 4)
GT.RU.FS <- grangertest(All_FL_FS$FS.RU~EPU_q$RU, order = 4)

Granger.tests <- list(
  GT.US.EPU=GT.US.EPU,
  GT.US.FS=GT.US.FS,
  
  GT.UK.EPU=GT.UK.EPU,
  GT.UK.FS=GT.UK.FS,
  
  GT.JP.EPU=GT.JP.EPU,
  GT.JP.FS=GT.JP.FS,
  
  GT.CA.EPU=GT.CA.EPU,
  GT.CA.FS=GT.CA.FS,
  
  GT.BR.EPU=GT.BR.EPU,
  GT.BR.FS=GT.BR.FS,
  
  GT.CN.EPU=GT.CN.EPU,
  GT.CN.FS=GT.CN.FS,
  
  GT.IN.EPU=GT.IN.EPU,
  GT.IN.FS=GT.IN.FS,
  
  GT.RU.EPU=GT.RU.EPU,
  GT.RU.FS=GT.RU.FS
)