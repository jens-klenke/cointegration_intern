#Packages
library(dplyr)
#Data Example
data("mts-examples", package="MTS")



## Variables 
trend <- c('none', 'constant', "trend")
lags <- 5
crit <- 0.05






plot(qgdp$us, type = 'l' , ylim = c(0,max(c(qgdp$us, qgdp$uk))))
lines(qgdp$uk)


##### Johansen

j_lags <- lags + 1 # why? is it a stat thing? 

?urca::ca.jo()

johansen <- urca::ca.jo(qgdp%>%dplyr::select(us,uk), ecdet= 'none', K = 2)
summary(johansen)






