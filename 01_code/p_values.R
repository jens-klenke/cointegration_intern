### Parameters zum Versuchen

stat <- c(-2.4529, 21.8242, -2.1006, 10.0941)
pval <- rep(0, length(stat))

trendtype <- 2
nvar <- 3
N <- nrow(null_dist)        # nrow(df) ### Numbers of observations



### Obtain p-values

load(here::here('00_data/null_dist.rda'))
# Erstellen des Vectors  
basecase <- 44 * (trendtype - 1) + 4 * (nvar -2) 

### Calulating p-values

for (i in 1:4) {
    case = basecase + i
    if (i %in% c(1, 3)) {
        n <- sum(stat[i] > null_dist[, case])
        pval[i] <-  (n/N) + .000000000001
    }
    else {
        if (i %in% c(2, 4)) {
            n <- sum(stat[i] < null_dist[, case])
            pval[i] <-  (n/N) + .000000000001
        }
    }
}
pval


######
#load('null_dist.rda')

N <- nrow(null_dist)

basecase <- 44 * (trendtype - 1) + 4 * (nvar - 2)

for (i in 1:4) {
    case = basecase + i
    if (i %in% c(1, 3)) {
        n <- sum(pval.stat[i] > null_dist[, case])
        pval.stat[i] <-  (n/N) + .000000000001
    } else {
        if (i %in% c(2, 4)) {
            n <- sum(pval.stat[i] < null_dist[, case])
            pval.stat[i] <-  (n/N) + .000000000001
        }
    }
}


