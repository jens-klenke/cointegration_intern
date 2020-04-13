library(bayerhanck)

bayerhanck()

library(haven)

hendry <- read.table(here::here('06_christoph/Gregory data/clements_hendry95.dat'))

str(hendry)

#DATE        RN          V
hendryhendry%>%
    dplyr::mutate(V2 = as.numeric(V2),
                  V3 = as.numeric(V3))

A <- bayerhanck(V2 ~ V3, data = hendry, lags = 4)

summary(A)

# DATE       FOODPC       NDPC         NDSPC         WAGEFOOD     WAGEND       WAGENDS

cooley_ogaki <- read.table(here::here('06_christoph/Gregory data/cooley_ogaki96.dat'))


b <- bayerhanck(V2 ~ V3 + V4, data = cooley_ogaki)

summary(b)

####      FL           SL
martens <- read.table(here::here(
    '06_christoph/Gregory data/martens_kofman_vorst_1_98.dat'))


c <- bayerhanck(V1 ~ V2, data = martens, lags = 100)

summary(c)

##############


d <- bayerhanck(sp ~ibm, data = ibmspko, lags = 45, trend = 'trend')

summary(d)


###########Lütkephol 


l <- bayerhanck(linvestment ~ lincome + lconsumption, data = Lutkephol, lags = 2, trend = 'trend')
summary(l)


