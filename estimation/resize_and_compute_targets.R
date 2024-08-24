
# Resize dataset:
indic_fst <- which(format(DATA$date,"%Y")==start_year)
DATA <- DATA[indic_fst:dim(DATA)[1],]

# Targeted moments:
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
  target_std_10_rea = sd(DATA$TIPSY10,na.rm=TRUE)/100,
  IRP10 = 0 # 10-year inflation risk premium
)
