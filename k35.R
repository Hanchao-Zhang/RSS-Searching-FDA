ppkgs <- c('chron','doSNOW', 'refund', 'mgcv', 'fda', 'plotly', 'tidyverse', 'parallel', 'foreach', 'doParallel', 'doSNOW','dplyr')
install.packages( ppkgs[!ppkgs %in% rownames(installed.packages())], repos = "http://cran.us.r-project.org", dependencies = T)

library(chron)
library(refund)
library(mgcv)
library(fda)
library(plotly)
library(tidyverse)
library(parallel)
library(foreach)
library(doParallel)
library(doSNOW)
library(dplyr)

source('Functions.R')

n = 3
k = 35

dat_o2 <- read.csv("datall.csv", header=T)########

#dat_co2 = dat_co2[(1 + 10*(i-1)):(10*i), ]
#rownmames(dat_co2)= NULL

#length(unique(dat_co2$StudyID))

#dat_co2$y <- dat_co2$ETCO2
#dat_co2$X <- 1:nrow(dat_co2)


clnum <- detectCores()
clnum <- 4
cl <- makeCluster(clnum)
registerDoSNOW(cl)
pb <- txtProgressBar(max = 28, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
t1<-Sys.time()

pid <- 1:length(unique(dat_o2$StudyID))

errlist <- list()


for( i in pid){
    id <- unique(dat_o2$StudyID)[i]
    dati <- dat_o2[dat_o2$StudyID == id,]
    
    lno_grid<- 1: ceiling(nrow(dati)/6)
    loop <- 1:length(lno_grid)
    
    err <- foreach(x=lno_grid,.packages=c('chron', 'refund', 'mgcv', 'fda', 'tidyverse'), .options.snow = opts, .combine = rbind ) %dopar% {
        tryCatch({
            err_cal(data = dat_o2, co2 = F, i = i, k = k, n = n, lno = x, d = 'all', proplno = F, consecutive = F)
        }, error = function(x)NA)
    }
    errlist[[i]] <- data.frame(lno_grid = lno_grid, err = err)
    print(i)
}




close(pb)
stopCluster(cl)
t2<-Sys.time()
t2-t1



n = 3
k = 35

save(errlist, file = paste0('errlist_o2_k',k, '_n_', n ,'.Rdata')  )
