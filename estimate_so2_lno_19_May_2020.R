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

dat_co2 <- read.csv("ETCO2_2_No-0-127_cum_2019-09-05_ROSCadded_2019-09-24.csv", header=T)########

#dat_co2 = dat_co2[(1 + 10*(i-1)):(10*i), ]
#rownmames(dat_co2)= NULL

#length(unique(dat_co2$StudyID))

dat_co2$y <- dat_co2$ETCO2
dat_co2$X <- 1:nrow(dat_co2)

source('Functions.R')


clnum <- detectCores()
cl <- makeCluster(4)
registerDoSNOW(cl)
pb <- txtProgressBar(max = 28, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
t1<-Sys.time()

pid <- length(unique(dat_co2$StudyID))

errlist <- list()

i  = 1
for( i in 1:pid){
    id <- unique(dat_co2$StudyID)[i]
    dati <- dat_co2[dat_co2$StudyID == id,]
    a <- ceiling(nrow(dati)/120)
    b <- nrow(dati)/4
    a <- ifelse( a == 1, yes = 10, no = a)
    b <- ifelse( a < b, yes = b, no = a)
    lno_grid <- seq(a, b, by = ceiling(  (nrow(dati)/4)/ 40 ))
    loop <- 1:length(lno_grid)
    err <- foreach(x=lno_grid,.packages=c('chron', 'refund', 'mgcv', 'fda', 'tidyverse'), .options.snow = opts, .combine = rbind ) %dopar% err_cal(data = dat_co2, co2 = T, i = i, k = 10, n = 3, lno = x, d = 2, proplno = T, consecutive = T)
    errlist[[i]] <- data.frame(lno_grid = lno_grid, err = err)
    print(i)
}


close(pb)
stopCluster(cl)
t2<-Sys.time()
t2-t1

save(errlist, file = 'errlist_so2_k10.RData')


errlist
