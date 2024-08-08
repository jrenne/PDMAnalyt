

nb_d  <- 7
nb_rr <- 7

indic_dd <- which(res0$d ==all_d[nb_d])
indic_rr <- which(res0$rr==all_rr[nb_rr])

indic <- which(indic_dd %in% indic_rr)

plot(res0$q[indic])


