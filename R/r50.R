# read in data and X50$X7 is the daily data 
library(readr)
i2005 <- read_csv("aiqtt_data/i2005.txt", col_names = FALSE)
i2009 <- read_csv("aiqtt_data/i2009.txt", col_names = FALSE)
i2005ret <- diff(log(i2005$X7))
i2009ret <- diff(log(i2009$X7))
source('~/Projects/R/aiquant/R/autogarch.R')
i2005ret_fit = garchAuto(i2005ret, cores=2, trace=TRUE)
i2009ret_fit = garchAuto(i2009ret, cores=2, trace=TRUE)