# ==============================================================================
# Prepare chart showing average average redemption schedule
# ==============================================================================

# ==============================================================================
make_chart_issuances_data <- function(auctions,GDP,
                                      first_year,first_quarter,
                                      last_year=first_year,last_quarter,
                                      tips_filter="No",
                                      main.t="",
                                      indic_yearly=TRUE){
  
  nb_quarters <- 4*(last_year - first_year) + last_quarter - first_quarter + 1
  
  max_maturity   <- max(as.numeric(levels(as.factor(auctions$Maturity.in.quarters))))
  all_maturities <-  1:max_maturity
  
  all_Issuances_GDP <- matrix(NaN,max_maturity,nb_quarters)
  
  count_quarter <- 0
  for(year_issuances in first_year:last_year){
    if(year_issuances==first_year){
      fst_Q <- first_quarter
    }else{fst_Q <- 1}
    if(year_issuances==last_year){
      lst_Q <- last_quarter
    }else{lst_Q <- 4}
    for(Q in fst_Q:lst_Q){
      count_quarter <- count_quarter + 1
      
      DATE     <- as.Date(paste(year_issuances,"-",3*Q-2,"-01",sep=""))
      gdp_index  <- which(GDP$DATE==DATE)
      GDP_year <- GDP$GDP[gdp_index] * 10^9
      
      # Select auctions taking place in the relevant quarter:
      auctions.red <- 
        subset(auctions,
               (Auction.year==year_issuances)&(Auction.quarter==Q)&(TIPS==tips_filter))
      Issuances <- matrix(0,max_maturity,1)
      # count <- 0
      # for(m in all_maturities){
      #   count <- count + 1
      #   tempo <- subset(auctions.red,Maturity.in.quarters==m)
      #   Issuances[count] <- sum(tempo$Total.Accepted)
      # }
      if(dim(auctions.red)[1]>0){
        for(i in 1:dim(auctions.red)[1]){
          
          # Treat only if repaid at least next year:
          issue_year <- format(auctions.red$Auction.Date[i],"%Y")
          matur_year <- format(auctions.red$Maturity.Date[i],"%Y")
          mat.in.q <- auctions.red$Maturity.in.quarters[i]
          mat.in.y <- auctions.red$Maturity.in.years[i]
          #if(matur_year>issue_year){
          if(mat.in.q>=4){
            Issuances[mat.in.q] <- Issuances[mat.in.q] + auctions.red$Total.Accepted[i]
            # Add interest payments:
            Issuances[1:mat.in.q] <- Issuances[1:mat.in.q] + 
              auctions.red$Interest.Rate[i]/4 * auctions.red$Total.Accepted[i]
          }
        }
      }
      Issuances_GDP <- Issuances/GDP_year
      all_Issuances_GDP[,count_quarter] <- Issuances_GDP
    }
  }
  
  if(isTRUE(indic_yearly)){
    # check that max_maturity multiple of 4: 
    if(4*trunc(max_maturity/4)!=max_maturity){
      new_max_maturity <- 4*(1 + trunc(max_maturity/4))
      addit_q <- new_max_maturity - max_maturity
      all_Issuances_GDP <- 
        rbind(all_Issuances_GDP,matrix(0,addit_q,dim(all_Issuances_GDP)[2]))
    }
    max_mat_in_years <- dim(all_Issuances_GDP)[1]/4
    all_Issuances_GDP <- (diag(max_mat_in_years) %x% matrix(1,1,4)) %*% all_Issuances_GDP
    all_maturities <-  1:dim(all_Issuances_GDP)[1]
  }
  
  print(all_Issuances_GDP)
  
  avg_Issuances_GDP <- apply(all_Issuances_GDP,1,function(x){mean(x,na.rm = TRUE)})
  if(last_year>first_year){
    std_Issuances_GDP <- apply(all_Issuances_GDP,1,function(x){sd(x,na.rm = TRUE)})
    # upper_bound <- avg_Issuances_GDP+std_Issuances_GDP
    # lower_bound <- avg_Issuances_GDP-std_Issuances_GDP
    lower_bound <- apply(all_Issuances_GDP,1,
                         function(x){quantile(x,.25,na.rm = TRUE)})
    upper_bound <- apply(all_Issuances_GDP,1,
                         function(x){quantile(x,.75,na.rm = TRUE)})
  }else{
    upper_bound <- avg_Issuances_GDP
  }
  
  # Correct for the fact that quarterly data:
  avg_Issuances_GDP <- 4 * avg_Issuances_GDP
  lower_bound <- 4 * lower_bound
  upper_bound <- 4 * upper_bound
  
  fig_bar <- barplot(avg_Issuances_GDP ~ all_maturities,
                     xlab="Maturities, in years",
                     ylab="Amount issued, percent of GDP",
                     main=main.t,las=1,
                     ylim=c(0,1.1*max(upper_bound)))
  if(last_year>first_year){
    points(x = fig_bar,lower_bound,pch=3,lwd=2)
    points(x = fig_bar,upper_bound,pch=3,lwd=2)
    for(i in 1:length(fig_bar)){
      lines(c(fig_bar[i],fig_bar[i]),
            c(lower_bound[i],upper_bound[i]),lty=3)
    }
  }
  return(list(avg_Issuances_GDP=avg_Issuances_GDP,
              lower_bound = lower_bound,
              upper_bound = upper_bound,
              all_Issuances_GDP = all_Issuances_GDP))
}

make_chart_schedule_data <- function(auctions,GDP,
                                     first_year,first_quarter,
                                     last_year,last_quarter,
                                     tips_filter="No",
                                     main.t = "",
                                     indic_yearly=TRUE){
  
  nb_quarters <- 4*(last_year - first_year) + last_quarter - first_quarter + 1
  
  max_maturity   <- max(as.numeric(levels(as.factor(auctions$Maturity.in.quarters))))
  all_maturities <-  1:max_maturity
  
  redemptions <- matrix(0,max_maturity,nb_quarters)
  
  count_quarter <- 0
  for(year_issuances in first_year:last_year){
    if(year_issuances==first_year){fst_Q <- first_quarter
    }else{fst_Q <- 1}
    if(year_issuances==last_year){lst_Q <- last_quarter
    }else{lst_Q <- 4}
    for(Q in fst_Q:lst_Q){
      count_quarter <- count_quarter + 1
      
      DATE     <- as.Date(paste(year_issuances,"-",3*Q-2,"-01",sep=""))
      gdp_index  <- which(GDP$DATE==DATE)
      GDP_year <- GDP$GDP[gdp_index] * 10^9
      
      DATE_year    <- as.numeric(format(DATE,"%Y"))
      DATE_month   <- as.numeric(format(DATE,"%m"))
      DATE_quarter <- 1 + trunc((DATE_month-1)/3)
      
      auctions$alive <- auctions$Maturity.Date > DATE
      auctions.red   <- subset(auctions,
                               (alive==1)&(Issue.Date<DATE)&(TIPS==tips_filter))
      auctions.red$residual_maturity_in_quarters <-
        4*(auctions.red$Maturing.year - DATE_year) +
        auctions.red$Maturing.quarter - DATE_quarter
      
      if(dim(auctions.red)[1]>0){
        for(i in 1:dim(auctions.red)[1]){
          mat.in.q <- auctions.red$residual_maturity_in_quarters[i]
          
          redemptions[mat.in.q,count_quarter] <- redemptions[mat.in.q,count_quarter] +
            auctions.red$Total.Accepted[i]/GDP_year
          # Add interest payments:
          redemptions[1:mat.in.q,count_quarter] <- redemptions[1:mat.in.q,count_quarter] + 
            auctions.red$Interest.Rate[i]/4 * auctions.red$Total.Accepted[i]/GDP_year
        }
      }
    }
  }
  
  if(isTRUE(indic_yearly)){
    # check that max_maturity multiple of 4: 
    if(4*trunc(max_maturity/4)!=max_maturity){
      new_max_maturity <- 4*(1 + trunc(max_maturity/4))
      addit_q <- new_max_maturity - max_maturity
      redemptions <- 
        rbind(redemptions,matrix(0,addit_q,dim(redemptions)[2]))
    }
    max_mat_in_years <- dim(redemptions)[1]/4
    redemptions <- (diag(max_mat_in_years) %x% matrix(1,1,4)) %*% redemptions
    all_maturities <-  1:dim(redemptions)[1]
  }
  
  avg_redemptions <- apply(redemptions,1,function(x){mean(x,na.rm = TRUE)})
  # ==============================================
  # ==============================================
  # ==============================================
  #if(length(ddates)>1){
  # ==============================================
  # ==============================================
  # ==============================================
  if(nb_quarters>1){
    std_redemptions <- apply(redemptions,1,function(x){sd(x,na.rm = TRUE)})
    # upper_bound <- avg_redemptions+std_redemptions
    # lower_bound <- avg_redemptions-std_redemptions
    lower_bound <- apply(redemptions,1,
                         function(x){quantile(x,.25,na.rm = TRUE)})
    upper_bound <- apply(redemptions,1,
                         function(x){quantile(x,.75,na.rm = TRUE)})
  }else{
    upper_bound <- avg_redemptions
  }
  
  fig_bar <- barplot(avg_redemptions ~ all_maturities,
                     xlab="Maturities, in years",
                     ylab="Debt outstanding, percent of GDP",
                     main=main.t,las=1,
                     ylim=c(0,1.1*max(upper_bound)))
  if(nb_quarters>1){
    points(x = fig_bar,lower_bound,pch=3,lwd=2)
    points(x = fig_bar,upper_bound,pch=3,lwd=2)
    for(i in 1:length(fig_bar)){
      lines(c(fig_bar[i],fig_bar[i]),
            c(lower_bound[i],upper_bound[i]),lty=3)
    }
  }
  
  # barplot(avg_redemptions ~ Maturities,
  #         main=main.t)
  
  return(list(
    avg_redemptions = avg_redemptions,
    lower_bound = lower_bound,
    upper_bound = upper_bound,
    redemptions = redemptions))
}
# ==============================================================================




# Load GDP
GDP <- read.csv("Data/US Debt Portfolio/GDP.csv")
GDP$DATE <- as.Date(GDP$DATE)
plot(GDP$DATE,GDP$GDP,type="l")

for(i in 1:11){
  File <- paste("Data/US Debt Portfolio/Securities-",i,".csv",sep="")
  Securities <- read.csv(File)
  if(i==1){
    auctions <- Securities
  }else{
    auctions <- rbind(auctions,Securities)
  }
}

auctions$Total.Accepted <- gsub(",","",auctions$Total.Accepted)
auctions$Total.Accepted <- gsub("\\$","",auctions$Total.Accepted)
auctions$Total.Accepted <- as.numeric(auctions$Total.Accepted)

auctions$Offering.Amount <- gsub(",","",auctions$Offering.Amount)
auctions$Offering.Amount <- gsub("\\$","",auctions$Offering.Amount)
auctions$Offering.Amount <- as.numeric(auctions$Offering.Amount)

auctions$Interest.Rate <- as.numeric(gsub("%","",auctions$Interest.Rate))/100
auctions$Interest.Rate[is.na(auctions$Interest.Rate)] <- 0

auctions$Maturity.Date <- as.Date(auctions$Maturity.Date,"%m/%d/%Y")
auctions$Issue.Date    <- as.Date(auctions$Issue.Date,"%m/%d/%Y")
auctions$Auction.Date  <- as.Date(auctions$Auction.Date,"%m/%d/%Y")

auctions$Maturing.year    <- as.numeric(format(auctions$Maturity.Date,"%Y"))
auctions$Maturing.month   <- as.numeric(format(auctions$Maturity.Date,"%m"))
auctions$Maturing.quarter <- 1+trunc((auctions$Maturing.month-1)/3)

auctions$Auction.year    <- as.numeric(format(auctions$Auction.Date,"%Y"))
auctions$Auction.month   <- as.numeric(format(auctions$Auction.Date,"%m"))
auctions$Auction.quarter <- 1+trunc((auctions$Auction.month-1)/3)

# auctions$Maturity.in.years <- as.numeric(format(auctions$Maturity.Date,"%Y")) -
#   as.numeric(format(auctions$Issue.Date,"%Y"))
auctions$Maturity.in.quarters <- 4*(auctions$Maturing.year - auctions$Auction.year) +
  auctions$Maturing.quarter - auctions$Auction.quarter
auctions$Maturity.in.years    <- auctions$Maturing.year - auctions$Auction.year

# # Compute average issuance schedule (over GDP): --------------------------------
# 
# res_Issuances_data_nom <- make_chart_issuances_data(auctions,GDP,
#                                                     first_year=2000,
#                                                     first_quarter=1,
#                                                     last_year = 2023,
#                                                     last_quarter = 4,
#                                                     tips_filter="No",
#                                                     main.t = "Issuances - nominal bonds",
#                                                     indic_yearly=TRUE)
# res_Issuances_data_ILB <- make_chart_issuances_data(auctions,GDP,
#                                                     first_year=2000,
#                                                     first_quarter=1,
#                                                     last_year = 2023,
#                                                     last_quarter = 4,
#                                                     tips_filter="Yes",
#                                                     main.t = "Issuances - ILBs",
#                                                     indic_yearly=TRUE)

# Compute average repayment schedule (over GDP): -------------------------------

FILE = paste("figures/Figure_avg_US_redempt.pdf",sep="")
pdf(file=FILE, pointsize=10, width=6, height=4)

par(plt=c(.15,.95,.2,.95))
res_Schedule_data_nom <- make_chart_schedule_data(auctions,GDP,
                                                  first_year=2000,
                                                  first_quarter=1,
                                                  last_year = 2023,
                                                  last_quarter = 4,
                                                  tips_filter = "No",
                                                  main.t = "",
                                                  indic_yearly=TRUE)
lines(.10 * .84^(0:30),lwd=2,col="red")

dev.off()

# res_Schedule_data_ILB <- make_chart_schedule_data(auctions,GDP,
#                                                   first_year=2000,
#                                                   first_quarter=1,
#                                                   last_year = 2023,
#                                                   last_quarter = 4,
#                                                   tips_filter = "Yes",
#                                                   main.t = "Repayment schedule - ILBs",
#                                                   indic_yearly=TRUE)
