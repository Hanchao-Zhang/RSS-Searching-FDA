library(chron)
library(refund)
library(mgcv)
library(fda)
library(chron)
library(plotly)
library(tidyverse)
library(parallel)
library(foreach)
library(doParallel)
library(doSNOW)
library(dplyr)




err_cal <- function(data = dat_co2, co2 = T, i = 1, k = 10, n = 3, lno = 1, d = 5, proplno = T, consecutive = T){

  id <- unique(data$StudyID)[i]
  dati <- data[data$StudyID == id,] ## data need to be changed
  diff <- unique(  as.numeric(chron::times(gsub(" PM", "", dati$Time)))  )[2]- unique(  as.numeric(chron::times(gsub(" PM", "", dati$Time)))  )[1]
  if( co2 == T){
    dati$Time2 =  (as.numeric(chron::times(gsub(" PM", "", dati$Time))) - as.numeric(chron::times(gsub(" PM", "", dati$Time)[dim(dati)[1]]))) / diff
  }else{
    dati$Time2 = 4*   (as.numeric(chron::times(gsub(" PM", "", dati$Time))) - as.numeric(chron::times(gsub(" PM", "", dati$Time)[dim(dati)[1]]))) / diff
  }

  obs_all <- dati$StudyID

  if(proplno == T){
    lno <- ceiling(dim(dati)[1]/lno)
  }else{
    #lno <- ceiling(dim(dati)[1]/sum(diff(dati$y) < 6)) + lno
    lno <- lno
  }

  if( lno == dim(dati)[1] ){
    stop('prop lno cannot equal to 1')
  }

  ## fit a gam function with all variables, jump to the next patients if the data cannot fit the gam function
  fit.gam_all <- tryCatch(
    {
      gam(as.numeric(y) ~ s(Time2, bs="cs", k=k), method="REML" , data = dati)
    }, error = function(x) NA
  )
  if( is.na(fit.gam_all)[1] ){
    print( (paste('patients',i,'cannot fit gam function'))  )
    return(NA)
    next
  }

  res_all <- fit.gam_all$residuals
  pred0 <- predict(fit.gam_all, se.fit = T)
  upp <- pred0$fit + n * pred0$se.fit
  low <- pred0$fit - n * pred0$se.fit
  indx3sdfit <- dati$y > low & dati$y < upp

  ## delete consecutive points
  if(consecutive == T){
    diff_0 <- which(diff(dati$y) == 0)
    diff_0_0 <- diff(diff_0)
    leadit<-function(x) x!=lead(x, default="what")
    rows <- leadit(dati$y)
    dati <- dati[rows, ]
  }

  ## Create a indicator variable for detecting oulier by the window method
  out1 <- rep(T, dim(dati)[1])
  diff_fit_pred <- rep(0, dim(dati)[1])

  for ( in_1 in ( 1: (dim(dati)[1]-lno) ) ){
    dati[!out1,] <- NA
    i_1 <- in_1:(in_1+lno)
    fit.gam <- tryCatch(
      {
        fit.gam <-gam(as.numeric(y) ~ s(Time2, bs="cs", k=k), method="REML" , data = dati[-i_1,])
      }, error = function(x) NA
    )

    if( is.na(fit.gam)[1] ){
      print( (paste('patients',i,'cannot fit gam function in the window process'))  )
      return(NA)
      next
    }

    pred1 <- predict(fit.gam, newdata = dati[i_1,], se.fit = T)
    diff_pred_y <- pred1$fit - dati[i_1,]$y
    diff_fit <- (fitted(fit.gam) - na.omit(dati$y[-i_1]) )
    if( in_1 <  ( dim(dati)[1]-  2 * lno) ){
      out1[i_1][c(1:d, (lno-d):lno )] <- (abs(diff_pred_y) < n * pred1$se.fit  )[c(1:d, (lno-d):lno )]
      input <- diff_pred_y[c(1:d, (lno-d):lno )]
      diff_fit_pred[i_1][  (c(1:d, (lno-d):lno)) [!is.na(input)]  ] <-  input[!is.na(input)]
      out1[is.na(out1)] <- F
    }else{
      if( 2*d > lno){
        out1[i_1] <- (abs(diff_pred_y) < n * pred1$se.fit  )
        input <- diff_pred_y[c(1:d, (lno-d):lno )]
        diff_fit_pred[i_1][!is.na(input)]  <-  input[!is.na(input)]
        out1[is.na(out1)] <- F
      }else{
        out1[i_1][c(1:d, (lno-d):lno )] <- (abs(diff_pred_y) < n * pred1$se.fit  )[c(1:d, (lno-d):lno )]
        input <- diff_pred_y[c(1:d, (lno-d):lno )]
        diff_fit_pred[i_1][  (c(1:d, (lno-d):lno)) [!is.na(input)]  ] <-  input[!is.na(input)]
        out1[is.na(out1)] <- F
      }
    }
  }


  dati <- na.omit(dati[out1,])

  indgam3 <- obs_all %in% dati$StudyID

  msep1 <- ( mean(res_all[!indgam3 & !indx3sdfit]^2, na.rm = T) + mean( diff_fit_pred[!indgam3 & !indx3sdfit]^2, na.rm = T) )/2
  msep1 <- ifelse(is.na(msep1), yes = 0, no = msep1)
  msep2 <- mean(diff_fit_pred[indgam3 & !indx3sdfit]^2, na.rm = T)
  msep2 <- ifelse(is.na(msep2), yes = 0, no = msep2)
  msep3 <- mean(res_all[!indgam3 & indx3sdfit]^2, na.rm = T)
  msep3 <- ifelse(is.na(msep3), yes = 0, no = msep3)

  return( msep1 -  msep2 - msep3 )

}



dat_o2 <- read.csv("datall.csv", header=T)########



plt_outlier_removal <- function(data = dat_o2, co2 = F, i = 10, k = 10, n = 3, lno = 46, d = 5, proplno = T, consecutive = T){

  id <- unique(data$StudyID)[i]
  dati <- data[data$StudyID == id,] ## data need to be changed
  diff <- unique(  as.numeric(chron::times(gsub(" PM", "", dati$Time)))  )[2]- unique(  as.numeric(chron::times(gsub(" PM", "", dati$Time)))  )[1]
  if( co2 == T){
    dati$Time2 =  (as.numeric(chron::times(gsub(" PM", "", dati$Time))) - as.numeric(chron::times(gsub(" PM", "", dati$Time)[dim(dati)[1]]))) / diff
  }else{
    dati$Time2 = 4*   (as.numeric(chron::times(gsub(" PM", "", dati$Time))) - as.numeric(chron::times(gsub(" PM", "", dati$Time)[dim(dati)[1]]))) / diff
  }

  obs_all <- dati$StudyID

  if(proplno == T){
    lno <- ceiling(dim(dati)[1]/lno)
  }else{
    #lno <- ceiling(dim(dati)[1]/sum(diff(dati$y) < 6)) + lno
    lno <- lno
  }

  if( lno == dim(dati)[1] ){
    print('prop lno cannot equal to 1')
    next
  }

  ## fit a gam function with all variables, jump to the next patients if the data cannot fit the gam function
  fit.gam_all <- tryCatch(
    {
      gam(as.numeric(y) ~ s(Time2, bs="cs", k=k), method="REML" , data = dati)
    }, error = function(x) NA
  )
  if( is.na(fit.gam_all)[1] ){
    print( (paste('patients',i,'cannot fit gam function'))  )
    return(NA)
    next
  }
  par(mfrow = c(1,2))

  plot(dati$Time2, dati$y, type="n",
       xlab = "Duration of CPR Monitoring (Time in Seconds Before End of CPR)",
       ylab="rSO2 %",
       # xlim=c(-1600,0), ylim=c(0,100),
       xlim=c(min(dati$Time2),0), ylim=c(0,100),
       main="CPR: rSO2 versus Time")
  points(dati$Time2, dati$y, pch=19, col=2, cex=.3)
  lines(dati$Time2, fitted(fit.gam_all), lwd=4, col=4)
  #dev.off()
  text(min(dati$Time2)/2, 80, id)


  res_all <- fit.gam_all$residuals
  pred0 <- predict(fit.gam_all, se.fit = T)
  upp <- pred0$fit + n * pred0$se.fit
  low <- pred0$fit - n * pred0$se.fit
  indx3sdfit <- dati$y > low & dati$y < upp

  ## delete consecutive points
  if(consecutive == T){
    diff_0 <- which(diff(dati$y) == 0)
    diff_0_0 <- diff(diff_0)
    leadit<-function(x) x!=lead(x, default="what")
    rows <- leadit(dati$y)
    dati <- dati[rows, ]
  }

  ## Create a indicator variable for detecting oulier by the window method
  out1 <- rep(T, dim(dati)[1])
  diff_fit_pred <- rep(0, dim(dati)[1])

  for ( in_1 in ( 1: (dim(dati)[1]-lno) ) ){
    dati[!out1,] <- NA
    i_1 <- in_1:(in_1+lno)
    fit.gam <- tryCatch(
      {
        fit.gam <-gam(as.numeric(y) ~ s(Time2, bs="cs", k=k), method="REML" , data = dati[-i_1,])
      }, error = function(x) NA
    )

    if( is.na(fit.gam)[1] ){
      print( (paste('patients',i,'cannot fit gam function in the window process'))  )
      return(NA)
      next
    }

    pred1 <- predict(fit.gam, newdata = dati[i_1,], se.fit = T)
    diff_pred_y <- pred1$fit - dati[i_1,]$y
    diff_fit <- (fitted(fit.gam) - na.omit(dati$y[-i_1]) )
    if( in_1 <  ( dim(dati)[1]-  2 * lno) ){
      out1[i_1][c(1:d, (lno-d):lno )] <- (abs(diff_pred_y) < n * pred1$se.fit  )[c(1:d, (lno-d):lno )]
      input <- diff_pred_y[c(1:d, (lno-d):lno )]
      diff_fit_pred[i_1][  (c(1:d, (lno-d):lno)) [!is.na(input)]  ] <-  input[!is.na(input)]
      out1[is.na(out1)] <- F
    }else{
      if( 2*d > lno){
        out1[i_1] <- (abs(diff_pred_y) < n * pred1$se.fit  )
        input <- diff_pred_y[c(1:d, (lno-d):lno )]
        diff_fit_pred[i_1][!is.na(input)]  <-  input[!is.na(input)]
        out1[is.na(out1)] <- F
      }else{
        out1[i_1][c(1:d, (lno-d):lno )] <- (abs(diff_pred_y) < n * pred1$se.fit  )[c(1:d, (lno-d):lno )]
        input <- diff_pred_y[c(1:d, (lno-d):lno )]
        diff_fit_pred[i_1][  (c(1:d, (lno-d):lno)) [!is.na(input)]  ] <-  input[!is.na(input)]
        out1[is.na(out1)] <- F
      }
    }
  }


  dati <- na.omit(dati[out1,])
  try(
   {
     fit.gam <-gam(as.numeric(y) ~ s(Time2, bs="cs", k=k), method="REML" , data = dati)
     plot(dati$Time2, dati$y, type="n",
          xlab = "Duration of CPR Monitoring (Time in Seconds Before End of CPR)",
          ylab="rSO2 %",
          # xlim=c(-1600,0), ylim=c(0,100),
          xlim=c(min(dati$Time2),0), ylim=c(0,100),
          main="CPR: rSO2 versus Time")
     points(dati$Time2, dati$y, pch=19, col=2, cex=.3)
     lines(na.omit(dati$Time2), fitted(fit.gam), lwd=4, col=4)
     #dev.off()
     text(min(dati$Time2)/2, 80, id)
   }
 )

return(dati)
  #indgam3 <- obs_all %in% dati$StudyID

  #msep1 <- ( mean(res_all[!indgam3 & !indx3sdfit]^2, na.rm = T) + mean( diff_fit_pred[!indgam3 & !indx3sdfit]^2, na.rm = T) )/2
  #msep1 <- ifelse(is.na(msep1), yes = 0, no = msep1)
  #msep2 <- mean(diff_fit_pred[indgam3 & !indx3sdfit]^2, na.rm = T)
  #msep2 <- ifelse(is.na(msep2), yes = 0, no = msep2)
  #msep3 <- mean(res_all[!indgam3 & indx3sdfit]^2, na.rm = T)
  #msep3 <- ifelse(is.na(msep3), yes = 0, no = msep3)

  #return( msep1 -  msep2 - msep3 )

}
