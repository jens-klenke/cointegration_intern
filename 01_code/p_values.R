### Parameters zum Versuchen

trendtype <- 2
nvar <- 3
N <- 10000        # nrow(df) ### Numbers of observations

stat <- c(-2, 21.8242, 3, 4)

### Obtain p-values

load(here::here('/null_dist.rda'))
pval <- stat # Erstellen des Vectors  
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


