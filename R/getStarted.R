X50 <- readRDS("~/Projects/R/aiquant/aiqtt_data/X50.rds")
ret50 <- readRDS("~/Projects/R/aiquant/aiqtt_data/ret50.rds")
fitX50 = garchAuto(ret50, cores=2, trace=TRUE)
plot(fitX50)
