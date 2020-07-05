library(chron)
library(refund)
library(mgcv)
library(fda)
library(tidyverse)
library(parallel)
library(foreach)
library(doParallel)
library(doSNOW)
library(dplyr)

#dat_o2 <- read.csv("datall.csv", header=T)########

#i = 4

#data = dat_o2
#co2 = F
#k = 10
#n = 3
#lno = 10
#d = 'all'
#proplno = T
#consecutive = T
#imputation = 'na'


plt_outlier_removal <- function(data = dat_o2, co2 = F, i = 3, k = 10, n = 3, lno = 50, d = 2, proplno = T, consecutive = T, imputation = 'na', trim = 3){

  id <- unique(data$StudyID)[i]
  dati <- data[data$StudyID == id,] ## data need to be changed
  diff <- unique(  as.numeric(chron::times(gsub(" PM", "", dati$Time)))  )[2]- unique(  as.numeric(chron::times(gsub(" PM", "", dati$Time)))  )[1]
  if( co2 == T){
    dati$Time2 =  (as.numeric(chron::times(gsub(" PM", "", dati$Time))) - as.numeric(chron::times(gsub(" PM", "", dati$Time)[dim(dati)[1]]))) / diff
  }else{
    dati$Time2 = 4 *   (as.numeric(chron::times(gsub(" PM", "", dati$Time))) - as.numeric(chron::times(gsub(" PM", "", dati$Time)[dim(dati)[1]]))) / diff
  }


  if(proplno == T){
    lno <- ceiling(dim(dati)[1]/lno)
  }else{
    #lno <- ceiling(dim(dati)[1]/sum(diff(dati$y) < 6)) + lno
    lno <- lno
  }

  if( lno == dim(dati)[1] ){
    stop('prop lno cannot equal to 1')

  }




  dati <- dati[ dati$y > 25 & dati$y < 75, ]
    #dati$y[ dati$y < 25 | dati$y > 75 ] <- NA

  ind_nonmissing <- which(is.na(dati$y) == F)

  if( is.numeric(trim) ){

      trim <- ceiling(trim * 60 / 4)

      ind_begin <- ind_nonmissing[trim]

      ind_end <- ind_nonmissing[ length(ind_nonmissing) - trim ]

      if( ind_begin <  ind_end){
          dati <- dati[ind_begin:ind_end, ]
      }else{
          next
      }

  }


  datii <- dati
  obs_all <- dati$X


  ## fit a gam function with all variables, jump to the next patients if the data cannot fit the gam function
  fit.gam_all <- tryCatch(
    {
      gam(as.numeric(y) ~ s(Time2, bs="cs", k=k), method="REML" , data = dati)
    }, error = function(x) NA
  )
  if( is.na(fit.gam_all)[1] ){
    #print( (paste('patients',i,'cannot fit gam function'))  )
    return(datii)
    next
  }
  par(mfrow = c(1,2))

  plot(dati$Time2, dati$y, type="n",
       xlab = "Duration of CPR Monitoring (Time in Seconds Before End of CPR)",
       ylab="rSO2 %",
       # xlim=c(-1600,0), ylim=c(0,100),
       xlim=c(min(dati$Time2),0), ylim=c(0,100),
       main="CPR: rSO2 versus Time")
  points(dati$Time2, dati$y, pch=19, col='blue', cex=.1)
  lines(dati$Time2, fitted(fit.gam_all), lwd=3, col=alpha(4, 0.5))
  #dev.off()
  text(min(dati$Time2)/2, 80, id)



  fit.gam_all <- tryCatch(
    {
      gam(as.numeric(y) ~ s(Time2, bs="cs", k=k), method="REML" , data = dati)
    }, error = function(x) NA
  )
  if( is.na(fit.gam_all)[1] ){
    #print( (paste('patients',i,'cannot fit gam function'))  )
    #return((paste('patients',i,'cannot fit gam function')) )
    next
  }

  res_all <- fit.gam_all$residuals
  pred0 <- predict(fit.gam_all, newdata = dati,se.fit = T)

  upp <- pred0$fit + n * pred0$se.fit
  low <- pred0$fit - n * pred0$se.fit
  indx3sdfit <- dati$y > low & dati$y < upp

  ## delete consecutive points
  if(consecutive == T){
    diff_0 <- which(diff(dati$y) == 0)
    diff_0_0 <- diff(diff_0)
    leadit<-function(x) x!=lead(x, default="what")
    rows <- leadit(dati$y)
    dati$y[!rows] <- NA
    #dati <- dati[rows, ]
  }

  ## Create a indicator variable for detecting oulier by the window method
  out1 <- rep(T, dim(dati)[1])
  diff_fit_pred <- rep(0, dim(dati)[1])
  imputed_vec <- rep(NA, dim(dati)[1])



  for ( in_1 in ( 1: (dim(dati)[1]-lno) ) ){

    i_1 <- in_1:(in_1+lno)
    fit.gam <- tryCatch(
      {
        fit.gam <-gam(as.numeric(y) ~ s(Time2, bs="cs", k=k), method="REML" , data = dati[-i_1,] )
      }, error = function(x) NA
    )

    if( is.na(fit.gam)[1] ){
      #print( (paste('patients',i,'cannot fit gam function in the window process', 'in the', in_1, 'loop'))  )
      #return((paste('patients',i,'cannot fit gam function in the window process', 'in the', in_1, 'loop')) )
      next
    }

    pred1 <- predict(fit.gam, newdata = dati[i_1,], se.fit = T)
    y_imputed <- dati[i_1,]$y
    y_imputed[is.na(y_imputed)] <- fitted(fit.gam)[i_1][is.na(y_imputed)]

    diff_pred_y <- pred1$fit - y_imputed
    diff_fit <- (fitted(fit.gam) - na.omit(dati$y[-i_1]) )


    if( d == 'all' ){
        d = 10000
    }



    if( imputation == 'weighted'){
        if( 2*d < lno ){
            out1[i_1][c(1:d, (lno-d):lno )] <- (abs(diff_pred_y) < n * pred1$se.fit  )[c(1:d, (lno-d):lno )]
            ind_outliers <- ( abs(diff_pred_y) > n * pred1$se.fit  )[c(1:d, (lno-d):lno )]
            ind_outliers <- ind_outliers[!is.na(ind_outliers)]
            dati$y[i_1][c(1:d, (lno-d):lno )][ind_outliers] <- 0.5 * (pred1$fit[c(1:d, (lno-d):lno )][ind_outliers] + fitted(fit.gam_all)[i_1][c(1:d, (lno-d):lno )][ind_outliers])
            imputed_vec[i_1][c(1:d, (lno-d):lno )][ind_outliers] <- 0.5 * (pred1$fit[c(1:d, (lno-d):lno )][ind_outliers] + fitted(fit.gam_all)[i_1][c(1:d, (lno-d):lno )][ind_outliers])
        }else{
            out1[i_1] <- (abs(diff_pred_y) < n * pred1$se.fit  )
            ind_outliers <- ( abs(diff_pred_y) > n * pred1$se.fit  )
            ind_outliers <- ind_outliers[!is.na(ind_outliers)]
            dati$y[i_1][ind_outliers] <- 0.5 * (pred1$fit[ind_outliers] + fitted(fit.gam_all)[i_1][ind_outliers])
            imputed_vec[i_1][ind_outliers] <- 0.5 * (pred1$fit[ind_outliers] + fitted(fit.gam_all)[i_1][ind_outliers])
        }

    }else if( imputation == 'fitted'){
        if( 2*d < lno ){
            out1[i_1][c(1:d, (lno-d):lno )] <- (abs(diff_pred_y) < n * pred1$se.fit  )[c(1:d, (lno-d):lno )]
            ind_outliers <- ( abs(diff_pred_y) > n * pred1$se.fit  )[c(1:d, (lno-d):lno )]
            ind_outliers <- ind_outliers[!is.na(ind_outliers)]
            dati$y[i_1][c(1:d, (lno-d):lno )][ind_outliers] <- fitted(fit.gam_all)[i_1][c(1:d, (lno-d):lno )][ind_outliers]
            imputed_vec[i_1][c(1:d, (lno-d):lno )][ind_outliers] <- fitted(fit.gam_all)[i_1][c(1:d, (lno-d):lno )][ind_outliers]
        }else{
            ind_outliers <- ( abs(diff_pred_y) > n * pred1$se.fit  )
            ind_outliers <- ind_outliers[!is.na(ind_outliers)]
            dati$y[i_1][ind_outliers] <- fitted(fit.gam_all)[i_1][ind_outliers]
            imputed_vec[i_1][ind_outliers] <- fitted(fit.gam_all)[i_1][ind_outliers]
        }

    }else if( imputation == 'predicted' ){
        if( 2*d < lno ){
            out1[i_1][c(1:d, (lno-d):lno )] <- (abs(diff_pred_y) < n * pred1$se.fit  )[c(1:d, (lno-d):lno )]
            ind_outliers <- ( abs(diff_pred_y) > n * pred1$se.fit  )[c(1:d, (lno-d):lno )]
            ind_outliers <- ind_outliers[!is.na(ind_outliers)]
            dati$y[i_1][c(1:d, (lno-d):lno )][ind_outliers] <- pred1$fit[c(1:d, (lno-d):lno )][ind_outliers]
            imputed_vec[i_1][c(1:d, (lno-d):lno )][ind_outliers] <- pred1$fit[c(1:d, (lno-d):lno )][ind_outliers]
        }else{
            out1[i_1] <- (abs(diff_pred_y) < n * pred1$se.fit  )
            ind_outliers <- ( abs(diff_pred_y) > n * pred1$se.fit  )
            ind_outliers <- ind_outliers[!is.na(ind_outliers)]
            dati$y[i_1][ind_outliers] <- pred1$fit[ind_outliers]
            imputed_vec[i_1][ind_outliers] <- pred1$fit[ind_outliers]
        }
    }else if( imputation == 'na' | imputation ==  'na.na' ){
        if( 2*d < lno ){
            out1[i_1][c(1:d, (lno-d):lno )] <- (abs(diff_pred_y) < n * pred1$se.fit  )[c(1:d, (lno-d):lno )]
            ind_outliers <- ( abs(diff_pred_y) > n * pred1$se.fit  )[c(1:d, (lno-d):lno )]
            ind_outliers <- ind_outliers[!is.na(ind_outliers)]
            dati$y[i_1][c(1:d, (lno-d):lno )][ind_outliers] <- NA
        }else{
            out1[i_1] <- (abs(diff_pred_y) < n * pred1$se.fit  )
            ind_outliers <- ( abs(diff_pred_y) > n * pred1$se.fit  )
            ind_outliers <- ind_outliers[!is.na(ind_outliers)]
            dati$y[i_1][ind_outliers] <- NA
        }
    }else if( imputation == 'na.weighted' ){
        if( 2*d < lno ){
            out1[i_1][c(1:d, (lno-d):lno )] <- (abs(diff_pred_y) < n * pred1$se.fit  )[c(1:d, (lno-d):lno )]
            ind_outliers <- ( abs(diff_pred_y) > n * pred1$se.fit  )[c(1:d, (lno-d):lno )]
            ind_outliers <- ind_outliers[!is.na(ind_outliers)]
            dati$y[i_1][c(1:d, (lno-d):lno )][ind_outliers] <- NA
            imputed_vec[i_1][c(1:d, (lno-d):lno )][ind_outliers] <- 0.5 * (pred1$fit[c(1:d, (lno-d):lno )][ind_outliers] + fitted(fit.gam_all)[i_1][c(1:d, (lno-d):lno )][ind_outliers])
        }else{
            out1[i_1] <- (abs(diff_pred_y) < n * pred1$se.fit  )
            ind_outliers <- ( abs(diff_pred_y) > n * pred1$se.fit  )
            ind_outliers <- ind_outliers[!is.na(ind_outliers)]
            dati$y[i_1][ind_outliers] <- NA
            imputed_vec[i_1][ind_outliers] <- 0.5 * (pred1$fit[ind_outliers] + fitted(fit.gam_all)[i_1][ind_outliers])
        }

    }else if( imputation == 'na.fitted' ){
        if( 2*d < lno ){
            out1[i_1][c(1:d, (lno-d):lno )] <- (abs(diff_pred_y) < n * pred1$se.fit  )[c(1:d, (lno-d):lno )]
            ind_outliers <- ( abs(diff_pred_y) > n * pred1$se.fit  )[c(1:d, (lno-d):lno )]
            ind_outliers <- ind_outliers[!is.na(ind_outliers)]
            imputed_vec[i_1][c(1:d, (lno-d):lno )][ind_outliers] <- fitted(fit.gam_all)[i_1][c(1:d, (lno-d):lno )][ind_outliers]
            dati$y[i_1][c(1:d, (lno-d):lno )][ind_outliers] <- NA
        }else{
            out1[i_1] <- (abs(diff_pred_y) < n * pred1$se.fit  )
            ind_outliers <- ( abs(diff_pred_y) > n * pred1$se.fit  )
            ind_outliers <- ind_outliers[!is.na(ind_outliers)]
            imputed_vec[i_1][ind_outliers] <- fitted(fit.gam_all)[i_1][ind_outliers]
            dati$y[i_1][ind_outliers] <- NA
        }
    }else if( imputation == 'na.predicted' ){
        if( 2*d < lno ){
            out1[i_1][c(1:d, (lno-d):lno )] <- (abs(diff_pred_y) < n * pred1$se.fit  )[c(1:d, (lno-d):lno )]
            ind_outliers <- ( abs(diff_pred_y) > n * pred1$se.fit  )[c(1:d, (lno-d):lno )]
            ind_outliers <- ind_outliers[!is.na(ind_outliers)]
            imputed_vec[i_1][c(1:d, (lno-d):lno )][ind_outliers] <- pred1$fit[c(1:d, (lno-d):lno )][ind_outliers]
            dati$y[i_1][c(1:d, (lno-d):lno )][ind_outliers] <- NA
        }else{
            out1[i_1] <- (abs(diff_pred_y) < n * pred1$se.fit  )
            ind_outliers <- ( abs(diff_pred_y) > n * pred1$se.fit  )
            ind_outliers <- ind_outliers[!is.na(ind_outliers)]
            imputed_vec[i_1][ind_outliers] <- pred1$fit[ind_outliers]
            dati$y[i_1][ind_outliers] <- NA
        }

    }

  }


  if( imputation == 'na.weighted' |  imputation == 'na.fitted' | imputation == 'na.predicted'){
      dati$y[is.na(dati$y)] <- imputed_vec[is.na(dati$y)]
  }



  fit.gam.smoothed <- tryCatch({
    gam(as.numeric(y) ~ s(Time2, bs="cs", k=k), method="REML" , data = dati)
  }, error = function(x)NA)

  if( is.na(fit.gam.smoothed)[1] ){
    #print( (paste('patients',i,'cannot fit gam function'))  )
    return(datii)
    next
  }

  pred2 <- predict(fit.gam.smoothed, newdata = dati, se.fit = T)
  if( imputation == 'na' ){
      #dati$y[ which((!out1)[is.na(dati$y)] == T) ] <- pred2$fit[  which((!out1)[is.na(dati$y)] == T)  ]
      #dati$y[ which(is.na(dati$y)[!out1] == T) ] <- pred2$fit[ which(is.na(dati$y)[!out1] == T) ]
      dati$y[is.na(dati$y)] <- tryCatch({
          pred2$fit[is.na(dati$y)]
      }, error = function(x)NA)
      #dati$y[ which((!out1)[is.na(dati$y)] == T) ] <- (pred2$fit[  which((!out1)[is.na(dati$y)] == T)  ] + imputed_vec[  which((!out1)[is.na(dati$y)] == T)  ]) * 0.5
 }

 if( imputation == 'na.na' ){
     #dati$y[ which((!out1)[is.na(dati$y)] == T) ] <- pred2$fit[  which((!out1)[is.na(dati$y)] == T)  ]
     #dati$y[ which(is.na(dati$y)[!out1] == T) ] <- pred2$fit[ which(is.na(dati$y)[!out1] == T) ]
     dati <- dati[!is.na(dati$y),]
 }




  fit.gam <- tryCatch({
    gam(as.numeric(y) ~ s(Time2, bs="cs", k=k), method="REML" , data = dati)
  }, error = function(x)NA)
  #dati <- na.omit(dati[out1,])

  pred3 <- predict(fit.gam, newdata = dati,se.fit = T)
  #rmind <- abs(pred3$fit - datii$y) > n * pred3$se.fit







  try(
   {
     plot(dati$Time2, dati$y, type="n",
          xlab = "Duration of CPR Monitoring (Time in Seconds Before End of CPR)",
          ylab="rSO2 %",
          # xlim=c(-1600,0), ylim=c(0,100),
          xlim=c(min(dati$Time2),0), ylim=c(0,100),
          main="CPR: rSO2 versus Time")
     points(dati$Time2, dati$y, pch=19, col=alpha(4, .5), cex=.1)
     #points(datii$Time2[rmind], datii$y[rmind], pch=19, col=alpha(2, .7), cex=.35)
     lines(na.omit(dati$Time2), predict(fit.gam, newdata = dati), lwd=3, col=alpha(4, .5) )
     legend('topright', paste0('k=',k, ', n=', n ))
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






err_cal <- function(data = dat_o2, co2 = T, i = 3, k = 10, n = 3, lno = 1, d = 5, proplno = T, consecutive = T, imputation = 'na'){

  id <- unique(data$StudyID)[i]
  dati <- data[data$StudyID == id,] ## data need to be changed

  dati$y[ dati$y < 15 | dati$y > 80 ] <- NA

  #dati$y[dati$y < 15 | dati$y > 80] <- NA

  diff <- unique(  as.numeric(chron::times(gsub(" PM", "", dati$Time)))  )[2]- unique(  as.numeric(chron::times(gsub(" PM", "", dati$Time)))  )[1]
  if( co2 == T){
    dati$Time2 =  (as.numeric(chron::times(gsub(" PM", "", dati$Time))) - as.numeric(chron::times(gsub(" PM", "", dati$Time)[dim(dati)[1]]))) / diff
  }else{
    dati$Time2 = 4*   (as.numeric(chron::times(gsub(" PM", "", dati$Time))) - as.numeric(chron::times(gsub(" PM", "", dati$Time)[dim(dati)[1]]))) / diff
  }



  obs_all <- dati$X

  if(proplno == T){
    lno <- ceiling(dim(dati)[1]/lno)
  }else{
    #lno <- ceiling(dim(dati)[1]/sum(diff(dati$y) < 6)) + lno
    lno <- lno
  }

  if( lno == dim(dati)[1] ){
    next('prop lno cannot equal to 1')
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
  pred0 <- predict(fit.gam_all, newdata = dati,se.fit = T)
  upp <- pred0$fit + n * pred0$se.fit
  low <- pred0$fit - n * pred0$se.fit
  indx3sdfit <- dati$y > low & dati$y < upp

  ## delete consecutive points
  if(consecutive == T){
    diff_0 <- which( diff(dati$y) == 0 )
    diff_0_0 <- diff(diff_0)
    leadit<-function(x) x!=lead(x, default="what")
    rows <- leadit(dati$y)
    dati$y[!rows] <- NA
  }

  ## Create a indicator variable for detecting oulier by the window method
  out1 <- rep(T, dim(dati)[1])
  diff_fit_pred <- rep(0, dim(dati)[1])



    for ( in_1 in ( 1: (dim(dati)[1]-lno) ) ){

      i_1 <- in_1:(in_1+lno)
      fit.gam <- tryCatch(
        {
          fit.gam <-gam(as.numeric(y) ~ s(Time2, bs="cs", k=k), method="REML" , data = dati[-i_1,] )
        }, error = function(x) NA
      )

      if( is.na(fit.gam)[1] ){
        #print( (paste('patients',i,'cannot fit gam function in the window process', 'in the', in_1, 'loop'))  )
        #return((paste('patients',i,'cannot fit gam function in the window process', 'in the', in_1, 'loop')) )
        next
      }

      pred1 <- predict(fit.gam, newdata = dati[i_1,], se.fit = T)
      y_imputed <- dati[i_1,]$y
      y_imputed[is.na(y_imputed)] <- fitted(fit.gam)[i_1][is.na(y_imputed)]

      diff_pred_y <- pred1$fit - y_imputed
      diff_fit <- (fitted(fit.gam) - na.omit(dati$y[-i_1]) )


      if( d == 'all' ){
          d = 10000
      }



      if( imputation == 'weighted'){
          if( 2*d < lno ){
              out1[i_1][c(1:d, (lno-d):lno )] <- (abs(diff_pred_y) < n * pred1$se.fit  )[c(1:d, (lno-d):lno )]
              ind_outliers <- ( abs(diff_pred_y) > n * pred1$se.fit  )[c(1:d, (lno-d):lno )]
              ind_outliers <- ind_outliers[!is.na(ind_outliers)]
              dati$y[i_1][c(1:d, (lno-d):lno )][ind_outliers] <- 0.5 * (pred1$fit[c(1:d, (lno-d):lno )][ind_outliers] + fitted(fit.gam_all)[i_1][c(1:d, (lno-d):lno )][ind_outliers])
              imputed_vec[i_1][c(1:d, (lno-d):lno )][ind_outliers] <- 0.5 * (pred1$fit[c(1:d, (lno-d):lno )][ind_outliers] + fitted(fit.gam_all)[i_1][c(1:d, (lno-d):lno )][ind_outliers])
          }else{
              out1[i_1] <- (abs(diff_pred_y) < n * pred1$se.fit  )
              ind_outliers <- ( abs(diff_pred_y) > n * pred1$se.fit  )
              ind_outliers <- ind_outliers[!is.na(ind_outliers)]
              dati$y[i_1][ind_outliers] <- 0.5 * (pred1$fit[ind_outliers] + fitted(fit.gam_all)[i_1][ind_outliers])
              imputed_vec[i_1][ind_outliers] <- 0.5 * (pred1$fit[ind_outliers] + fitted(fit.gam_all)[i_1][ind_outliers])
          }

      }else if( imputation == 'fitted'){
          if( 2*d < lno ){
              out1[i_1][c(1:d, (lno-d):lno )] <- (abs(diff_pred_y) < n * pred1$se.fit  )[c(1:d, (lno-d):lno )]
              ind_outliers <- ( abs(diff_pred_y) > n * pred1$se.fit  )[c(1:d, (lno-d):lno )]
              ind_outliers <- ind_outliers[!is.na(ind_outliers)]
              dati$y[i_1][c(1:d, (lno-d):lno )][ind_outliers] <- fitted(fit.gam_all)[i_1][c(1:d, (lno-d):lno )][ind_outliers]
              imputed_vec[i_1][c(1:d, (lno-d):lno )][ind_outliers] <- fitted(fit.gam_all)[i_1][c(1:d, (lno-d):lno )][ind_outliers]
          }else{
              ind_outliers <- ( abs(diff_pred_y) > n * pred1$se.fit  )
              ind_outliers <- ind_outliers[!is.na(ind_outliers)]
              dati$y[i_1][ind_outliers] <- fitted(fit.gam_all)[i_1][ind_outliers]
              imputed_vec[i_1][ind_outliers] <- fitted(fit.gam_all)[i_1][ind_outliers]
          }

      }else if( imputation == 'predicted' ){
          if( 2*d < lno ){
              out1[i_1][c(1:d, (lno-d):lno )] <- (abs(diff_pred_y) < n * pred1$se.fit  )[c(1:d, (lno-d):lno )]
              ind_outliers <- ( abs(diff_pred_y) > n * pred1$se.fit  )[c(1:d, (lno-d):lno )]
              ind_outliers <- ind_outliers[!is.na(ind_outliers)]
              dati$y[i_1][c(1:d, (lno-d):lno )][ind_outliers] <- pred1$fit[c(1:d, (lno-d):lno )][ind_outliers]
              imputed_vec[i_1][c(1:d, (lno-d):lno )][ind_outliers] <- pred1$fit[c(1:d, (lno-d):lno )][ind_outliers]
          }else{
              out1[i_1] <- (abs(diff_pred_y) < n * pred1$se.fit  )
              ind_outliers <- ( abs(diff_pred_y) > n * pred1$se.fit  )
              ind_outliers <- ind_outliers[!is.na(ind_outliers)]
              dati$y[i_1][ind_outliers] <- pred1$fit[ind_outliers]
              imputed_vec[i_1][ind_outliers] <- pred1$fit[ind_outliers]
          }
      }else if( imputation == 'na' | imputation ==  'na.na' ){
          if( 2*d < lno ){
              out1[i_1][c(1:d, (lno-d):lno )] <- (abs(diff_pred_y) < n * pred1$se.fit  )[c(1:d, (lno-d):lno )]
              ind_outliers <- ( abs(diff_pred_y) > n * pred1$se.fit  )[c(1:d, (lno-d):lno )]
              ind_outliers <- ind_outliers[!is.na(ind_outliers)]
              dati$y[i_1][c(1:d, (lno-d):lno )][ind_outliers] <- NA
          }else{
              out1[i_1] <- (abs(diff_pred_y) < n * pred1$se.fit  )
              ind_outliers <- ( abs(diff_pred_y) > n * pred1$se.fit  )
              ind_outliers <- ind_outliers[!is.na(ind_outliers)]
              dati$y[i_1][ind_outliers] <- NA
          }
      }else if( imputation == 'na.weighted' ){
          if( 2*d < lno ){
              out1[i_1][c(1:d, (lno-d):lno )] <- (abs(diff_pred_y) < n * pred1$se.fit  )[c(1:d, (lno-d):lno )]
              ind_outliers <- ( abs(diff_pred_y) > n * pred1$se.fit  )[c(1:d, (lno-d):lno )]
              ind_outliers <- ind_outliers[!is.na(ind_outliers)]
              dati$y[i_1][c(1:d, (lno-d):lno )][ind_outliers] <- NA
              imputed_vec[i_1][c(1:d, (lno-d):lno )][ind_outliers] <- 0.5 * (pred1$fit[c(1:d, (lno-d):lno )][ind_outliers] + fitted(fit.gam_all)[i_1][c(1:d, (lno-d):lno )][ind_outliers])
          }else{
              out1[i_1] <- (abs(diff_pred_y) < n * pred1$se.fit  )
              ind_outliers <- ( abs(diff_pred_y) > n * pred1$se.fit  )
              ind_outliers <- ind_outliers[!is.na(ind_outliers)]
              dati$y[i_1][ind_outliers] <- NA
              imputed_vec[i_1][ind_outliers] <- 0.5 * (pred1$fit[ind_outliers] + fitted(fit.gam_all)[i_1][ind_outliers])
          }

      }else if( imputation == 'na.fitted' ){
          if( 2*d < lno ){
              out1[i_1][c(1:d, (lno-d):lno )] <- (abs(diff_pred_y) < n * pred1$se.fit  )[c(1:d, (lno-d):lno )]
              ind_outliers <- ( abs(diff_pred_y) > n * pred1$se.fit  )[c(1:d, (lno-d):lno )]
              ind_outliers <- ind_outliers[!is.na(ind_outliers)]
              imputed_vec[i_1][c(1:d, (lno-d):lno )][ind_outliers] <- fitted(fit.gam_all)[i_1][c(1:d, (lno-d):lno )][ind_outliers]
              dati$y[i_1][c(1:d, (lno-d):lno )][ind_outliers] <- NA
          }else{
              out1[i_1] <- (abs(diff_pred_y) < n * pred1$se.fit  )
              ind_outliers <- ( abs(diff_pred_y) > n * pred1$se.fit  )
              ind_outliers <- ind_outliers[!is.na(ind_outliers)]
              imputed_vec[i_1][ind_outliers] <- fitted(fit.gam_all)[i_1][ind_outliers]
              dati$y[i_1][ind_outliers] <- NA
          }
      }else if( imputation == 'na.predicted' ){
          if( 2*d < lno ){
              out1[i_1][c(1:d, (lno-d):lno )] <- (abs(diff_pred_y) < n * pred1$se.fit  )[c(1:d, (lno-d):lno )]
              ind_outliers <- ( abs(diff_pred_y) > n * pred1$se.fit  )[c(1:d, (lno-d):lno )]
              ind_outliers <- ind_outliers[!is.na(ind_outliers)]
              imputed_vec[i_1][c(1:d, (lno-d):lno )][ind_outliers] <- pred1$fit[c(1:d, (lno-d):lno )][ind_outliers]
              dati$y[i_1][c(1:d, (lno-d):lno )][ind_outliers] <- NA
          }else{
              out1[i_1] <- (abs(diff_pred_y) < n * pred1$se.fit  )
              ind_outliers <- ( abs(diff_pred_y) > n * pred1$se.fit  )
              ind_outliers <- ind_outliers[!is.na(ind_outliers)]
              imputed_vec[i_1][ind_outliers] <- pred1$fit[ind_outliers]
              dati$y[i_1][ind_outliers] <- NA
          }

      }




    }


    if( imputation == 'na.weighted' |  imputation == 'na.fitted' | imputation == 'na.predicted'){
        dati$y[is.na(dati$y)] <- imputed_vec[is.na(dati$y)]
    }



  dati <- na.omit(dati[na.omit(out1),])


  indgam3 <- obs_all %in% dati$X


  msep1 <- ( mean(res_all[!indgam3 & !indx3sdfit]^2, na.rm = T) + mean( diff_fit_pred[!indgam3 & !indx3sdfit]^2, na.rm = T) )/2
  #msep1 <- mean( diff_fit_pred[!indgam3 & !indx3sdfit]^2, na.rm = T)
  msep1 <- ifelse(is.na(msep1), yes = 0, no = msep1)
  msep2 <- mean(diff_fit_pred[indgam3 & !indx3sdfit]^2, na.rm = T)
  #msep2 <- mean(res_all[indgam3 & !indx3sdfit]^2, na.rm = T)
  msep2 <- ifelse(is.na(msep2), yes = 0, no = msep2)
  msep3 <- mean(res_all[!indgam3 & indx3sdfit]^2, na.rm = T)
  #msep3 <- mean(diff_fit_pred[!indgam3 & indx3sdfit]^2, na.rm = T)
  msep3 <- ifelse(is.na(msep3), yes = 0, no = msep3)

  return( msep1 -  msep2 - msep3 )

}













plt_outlier_removal_onplt <- function(data = dat_o2, co2 = F, i = 3, k = 10, n = 3, lno = 50, d = 2, proplno = T, consecutive = T, imputation = 'na', trim = 3){

  id <- unique(data$StudyID)[i]
  dati <- data[data$StudyID == id,] ## data need to be changed
  diff <- unique(  as.numeric(chron::times(gsub(" PM", "", dati$Time)))  )[2]- unique(  as.numeric(chron::times(gsub(" PM", "", dati$Time)))  )[1]
  if( co2 == T){
    dati$Time2 =  (as.numeric(chron::times(gsub(" PM", "", dati$Time))) - as.numeric(chron::times(gsub(" PM", "", dati$Time)[dim(dati)[1]]))) / diff
  }else{
    dati$Time2 = 4 *   (as.numeric(chron::times(gsub(" PM", "", dati$Time))) - as.numeric(chron::times(gsub(" PM", "", dati$Time)[dim(dati)[1]]))) / diff
  }


  if(proplno == T){
    lno <- ceiling(dim(dati)[1]/lno)
  }else{
    #lno <- ceiling(dim(dati)[1]/sum(diff(dati$y) < 6)) + lno
    lno <- lno
  }

  if( lno == dim(dati)[1] ){
    stop('prop lno cannot equal to 1')

  }




  dati <- dati[ dati$y > 25 & dati$y < 75, ]
    #dati$y[ dati$y < 25 | dati$y > 75 ] <- NA

  ind_nonmissing <- which(is.na(dati$y) == F)

  if( is.numeric(trim) ){

      trim <- ceiling(trim * 60 / 4)

      ind_begin <- ind_nonmissing[trim]

      ind_end <- ind_nonmissing[ length(ind_nonmissing) - trim ]

      if( ind_begin <  ind_end){
          dati <- dati[ind_begin:ind_end, ]
      }else{
          next
      }

  }



  datii <- dati
  obs_all <- dati$X


  ## fit a gam function with all variables, jump to the next patients if the data cannot fit the gam function
  fit.gam_all <- tryCatch(
    {
      gam(as.numeric(y) ~ s(Time2, bs="cs", k=k), method="REML" , data = dati)
    }, error = function(x) NA
  )
  if( is.na(fit.gam_all)[1] ){
    #print( (paste('patients',i,'cannot fit gam function'))  )
    return(datii)
    next
  }
  #par(mfrow = c(1,2))

  plot(dati$Time2, dati$y, type="n",
       xlab = "Duration of CPR Monitoring (Time in Seconds Before End of CPR)",
       ylab="rSO2 %",
       # xlim=c(-1600,0), ylim=c(0,100),
       xlim=c(min(dati$Time2),0), ylim=c(0,100),
       main="CPR: rSO2 versus Time")
  points(dati$Time2, dati$y, pch=19, col='blue', cex=.1)
  #lines(dati$Time2, fitted(fit.gam_all), lwd=3, col=alpha(4, 0.5))
  #dev.off()
  #text(min(dati$Time2)/2, 80, id)



  fit.gam_all <- tryCatch(
    {
      gam(as.numeric(y) ~ s(Time2, bs="cs", k=k), method="REML" , data = dati)
    }, error = function(x) NA
  )
  if( is.na(fit.gam_all)[1] ){
    #print( (paste('patients',i,'cannot fit gam function'))  )
    #return((paste('patients',i,'cannot fit gam function')) )
    next
  }

  res_all <- fit.gam_all$residuals
  pred0 <- predict(fit.gam_all, newdata = dati,se.fit = T)

  upp <- pred0$fit + n * pred0$se.fit
  low <- pred0$fit - n * pred0$se.fit
  indx3sdfit <- dati$y > low & dati$y < upp

  ## delete consecutive points
  if(consecutive == T){
    diff_0 <- which(diff(dati$y) == 0)
    diff_0_0 <- diff(diff_0)
    leadit<-function(x) x!=lead(x, default="what")
    rows <- leadit(dati$y)
    dati$y[!rows] <- NA
    #dati <- dati[rows, ]
  }

  ## Create a indicator variable for detecting oulier by the window method
  out1 <- rep(T, dim(dati)[1])
  diff_fit_pred <- rep(0, dim(dati)[1])
  imputed_vec <- rep(NA, dim(dati)[1])



  for ( in_1 in ( 1: (dim(dati)[1]-lno) ) ){

    i_1 <- in_1:(in_1+lno)
    fit.gam <- tryCatch(
      {
        fit.gam <-gam(as.numeric(y) ~ s(Time2, bs="cs", k=k), method="REML" , data = dati[-i_1,] )
      }, error = function(x) NA
    )

    if( is.na(fit.gam)[1] ){
      #print( (paste('patients',i,'cannot fit gam function in the window process', 'in the', in_1, 'loop'))  )
      #return((paste('patients',i,'cannot fit gam function in the window process', 'in the', in_1, 'loop')) )
      next
    }

    pred1 <- predict(fit.gam, newdata = dati[i_1,], se.fit = T)
    y_imputed <- dati[i_1,]$y
    y_imputed[is.na(y_imputed)] <- fitted(fit.gam)[i_1][is.na(y_imputed)]

    diff_pred_y <- pred1$fit - y_imputed
    diff_fit <- (fitted(fit.gam) - na.omit(dati$y[-i_1]) )


    if( d == 'all' ){
        d = 10000
    }



    if( imputation == 'weighted'){
        if( 2*d < lno ){
            out1[i_1][c(1:d, (lno-d):lno )] <- (abs(diff_pred_y) < n * pred1$se.fit  )[c(1:d, (lno-d):lno )]
            ind_outliers <- ( abs(diff_pred_y) > n * pred1$se.fit  )[c(1:d, (lno-d):lno )]
            ind_outliers <- ind_outliers[!is.na(ind_outliers)]
            dati$y[i_1][c(1:d, (lno-d):lno )][ind_outliers] <- 0.5 * (pred1$fit[c(1:d, (lno-d):lno )][ind_outliers] + fitted(fit.gam_all)[i_1][c(1:d, (lno-d):lno )][ind_outliers])
            imputed_vec[i_1][c(1:d, (lno-d):lno )][ind_outliers] <- 0.5 * (pred1$fit[c(1:d, (lno-d):lno )][ind_outliers] + fitted(fit.gam_all)[i_1][c(1:d, (lno-d):lno )][ind_outliers])
        }else{
            out1[i_1] <- (abs(diff_pred_y) < n * pred1$se.fit  )
            ind_outliers <- ( abs(diff_pred_y) > n * pred1$se.fit  )
            ind_outliers <- ind_outliers[!is.na(ind_outliers)]
            dati$y[i_1][ind_outliers] <- 0.5 * (pred1$fit[ind_outliers] + fitted(fit.gam_all)[i_1][ind_outliers])
            imputed_vec[i_1][ind_outliers] <- 0.5 * (pred1$fit[ind_outliers] + fitted(fit.gam_all)[i_1][ind_outliers])
        }

    }else if( imputation == 'fitted'){
        if( 2*d < lno ){
            out1[i_1][c(1:d, (lno-d):lno )] <- (abs(diff_pred_y) < n * pred1$se.fit  )[c(1:d, (lno-d):lno )]
            ind_outliers <- ( abs(diff_pred_y) > n * pred1$se.fit  )[c(1:d, (lno-d):lno )]
            ind_outliers <- ind_outliers[!is.na(ind_outliers)]
            dati$y[i_1][c(1:d, (lno-d):lno )][ind_outliers] <- fitted(fit.gam_all)[i_1][c(1:d, (lno-d):lno )][ind_outliers]
            imputed_vec[i_1][c(1:d, (lno-d):lno )][ind_outliers] <- fitted(fit.gam_all)[i_1][c(1:d, (lno-d):lno )][ind_outliers]
        }else{
            ind_outliers <- ( abs(diff_pred_y) > n * pred1$se.fit  )
            ind_outliers <- ind_outliers[!is.na(ind_outliers)]
            dati$y[i_1][ind_outliers] <- fitted(fit.gam_all)[i_1][ind_outliers]
            imputed_vec[i_1][ind_outliers] <- fitted(fit.gam_all)[i_1][ind_outliers]
        }

    }else if( imputation == 'predicted' ){
        if( 2*d < lno ){
            out1[i_1][c(1:d, (lno-d):lno )] <- (abs(diff_pred_y) < n * pred1$se.fit  )[c(1:d, (lno-d):lno )]
            ind_outliers <- ( abs(diff_pred_y) > n * pred1$se.fit  )[c(1:d, (lno-d):lno )]
            ind_outliers <- ind_outliers[!is.na(ind_outliers)]
            dati$y[i_1][c(1:d, (lno-d):lno )][ind_outliers] <- pred1$fit[c(1:d, (lno-d):lno )][ind_outliers]
            imputed_vec[i_1][c(1:d, (lno-d):lno )][ind_outliers] <- pred1$fit[c(1:d, (lno-d):lno )][ind_outliers]
        }else{
            out1[i_1] <- (abs(diff_pred_y) < n * pred1$se.fit  )
            ind_outliers <- ( abs(diff_pred_y) > n * pred1$se.fit  )
            ind_outliers <- ind_outliers[!is.na(ind_outliers)]
            dati$y[i_1][ind_outliers] <- pred1$fit[ind_outliers]
            imputed_vec[i_1][ind_outliers] <- pred1$fit[ind_outliers]
        }
    }else if( imputation == 'na'  | imputation ==  'na.na'){
        if( 2*d < lno ){
            out1[i_1][c(1:d, (lno-d):lno )] <- (abs(diff_pred_y) < n * pred1$se.fit  )[c(1:d, (lno-d):lno )]
            ind_outliers <- ( abs(diff_pred_y) > n * pred1$se.fit  )[c(1:d, (lno-d):lno )]
            ind_outliers <- ind_outliers[!is.na(ind_outliers)]
            dati$y[i_1][c(1:d, (lno-d):lno )][ind_outliers] <- NA
        }else{
            out1[i_1] <- (abs(diff_pred_y) < n * pred1$se.fit  )
            ind_outliers <- ( abs(diff_pred_y) > n * pred1$se.fit  )
            ind_outliers <- ind_outliers[!is.na(ind_outliers)]
            dati$y[i_1][ind_outliers] <- NA
        }
    }else if( imputation == 'na.weighted' ){
        if( 2*d < lno ){
            out1[i_1][c(1:d, (lno-d):lno )] <- (abs(diff_pred_y) < n * pred1$se.fit  )[c(1:d, (lno-d):lno )]
            ind_outliers <- ( abs(diff_pred_y) > n * pred1$se.fit  )[c(1:d, (lno-d):lno )]
            ind_outliers <- ind_outliers[!is.na(ind_outliers)]
            dati$y[i_1][c(1:d, (lno-d):lno )][ind_outliers] <- NA
            imputed_vec[i_1][c(1:d, (lno-d):lno )][ind_outliers] <- 0.5 * (pred1$fit[c(1:d, (lno-d):lno )][ind_outliers] + fitted(fit.gam_all)[i_1][c(1:d, (lno-d):lno )][ind_outliers])
        }else{
            out1[i_1] <- (abs(diff_pred_y) < n * pred1$se.fit  )
            ind_outliers <- ( abs(diff_pred_y) > n * pred1$se.fit  )
            ind_outliers <- ind_outliers[!is.na(ind_outliers)]
            dati$y[i_1][ind_outliers] <- NA
            imputed_vec[i_1][ind_outliers] <- 0.5 * (pred1$fit[ind_outliers] + fitted(fit.gam_all)[i_1][ind_outliers])
        }

    }else if( imputation == 'na.fitted' ){
        if( 2*d < lno ){
            out1[i_1][c(1:d, (lno-d):lno )] <- (abs(diff_pred_y) < n * pred1$se.fit  )[c(1:d, (lno-d):lno )]
            ind_outliers <- ( abs(diff_pred_y) > n * pred1$se.fit  )[c(1:d, (lno-d):lno )]
            ind_outliers <- ind_outliers[!is.na(ind_outliers)]
            imputed_vec[i_1][c(1:d, (lno-d):lno )][ind_outliers] <- fitted(fit.gam_all)[i_1][c(1:d, (lno-d):lno )][ind_outliers]
            dati$y[i_1][c(1:d, (lno-d):lno )][ind_outliers] <- NA
        }else{
            out1[i_1] <- (abs(diff_pred_y) < n * pred1$se.fit  )
            ind_outliers <- ( abs(diff_pred_y) > n * pred1$se.fit  )
            ind_outliers <- ind_outliers[!is.na(ind_outliers)]
            imputed_vec[i_1][ind_outliers] <- fitted(fit.gam_all)[i_1][ind_outliers]
            dati$y[i_1][ind_outliers] <- NA
        }
    }else if( imputation == 'na.predicted' ){
        if( 2*d < lno ){
            out1[i_1][c(1:d, (lno-d):lno )] <- (abs(diff_pred_y) < n * pred1$se.fit  )[c(1:d, (lno-d):lno )]
            ind_outliers <- ( abs(diff_pred_y) > n * pred1$se.fit  )[c(1:d, (lno-d):lno )]
            ind_outliers <- ind_outliers[!is.na(ind_outliers)]
            imputed_vec[i_1][c(1:d, (lno-d):lno )][ind_outliers] <- pred1$fit[c(1:d, (lno-d):lno )][ind_outliers]
            dati$y[i_1][c(1:d, (lno-d):lno )][ind_outliers] <- NA
        }else{
            out1[i_1] <- (abs(diff_pred_y) < n * pred1$se.fit  )
            ind_outliers <- ( abs(diff_pred_y) > n * pred1$se.fit  )
            ind_outliers <- ind_outliers[!is.na(ind_outliers)]
            imputed_vec[i_1][ind_outliers] <- pred1$fit[ind_outliers]
            dati$y[i_1][ind_outliers] <- NA
        }

    }

  }


  if( imputation == 'na.weighted' |  imputation == 'na.fitted' | imputation == 'na.predicted'){
      dati$y[is.na(dati$y)] <- imputed_vec[is.na(dati$y)]
  }



  fit.gam.smoothed <- tryCatch({
    gam(as.numeric(y) ~ s(Time2, bs="cs", k=k), method="REML" , data = dati)
  }, error = function(x)NA)

  if( is.na(fit.gam.smoothed)[1] ){
    #print( (paste('patients',i,'cannot fit gam function'))  )
    return(datii)
    next
  }

  pred2 <- predict(fit.gam.smoothed, newdata = dati, se.fit = T)
  if( imputation == 'na' ){
      #dati$y[ which((!out1)[is.na(dati$y)] == T) ] <- pred2$fit[  which((!out1)[is.na(dati$y)] == T)  ]
      #dati$y[ which(is.na(dati$y)[!out1] == T) ] <- pred2$fit[ which(is.na(dati$y)[!out1] == T) ]
      dati$y[is.na(dati$y)] <- tryCatch({
          pred2$fit[is.na(dati$y)]
      }, error = function(x)NA)
      #dati$y[ which((!out1)[is.na(dati$y)] == T) ] <- (pred2$fit[  which((!out1)[is.na(dati$y)] == T)  ] + imputed_vec[  which((!out1)[is.na(dati$y)] == T)  ]) * 0.5
 }

 if( imputation == 'na.na' ){
     #dati$y[ which((!out1)[is.na(dati$y)] == T) ] <- pred2$fit[  which((!out1)[is.na(dati$y)] == T)  ]
     #dati$y[ which(is.na(dati$y)[!out1] == T) ] <- pred2$fit[ which(is.na(dati$y)[!out1] == T) ]
     dati <- dati[!is.na(dati$y),]
 }




  fit.gam <- tryCatch({
    gam(as.numeric(y) ~ s(Time2, bs="cs", k=k), method="REML" , data = dati)
  }, error = function(x)NA)
  #dati <- na.omit(dati[out1,])

  pred3 <- predict(fit.gam, newdata = dati,se.fit = T)
  #rmind <- abs(pred3$fit - datii$y) > n * pred3$se.fit







  try(
   {
     #plot(dati$Time2, dati$y, type="n",
    #      xlab = "Duration of CPR Monitoring (Time in Seconds Before End of CPR)",
    #      ylab="rSO2 %",
          # xlim=c(-1600,0), ylim=c(0,100),
    #      xlim=c(min(dati$Time2),0), ylim=c(0,100),
    #      main="CPR: rSO2 versus Time")
    #points(dati$Time2, dati$y, pch=19, col=alpha(4, .5), cex=.1)
     #points(datii$Time2[rmind], datii$y[rmind], pch=19, col=alpha(2, .7), cex=.35)
     lines(na.omit(dati$Time2), predict(fit.gam, newdata = dati), lwd=3, col=alpha('red', .5) )
     legend('topright', paste0('k=',k, ', n=', n ))
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
