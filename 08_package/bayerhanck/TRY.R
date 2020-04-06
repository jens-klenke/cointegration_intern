library(devtools)
load_all()


A <- bayerhanck(linvestment ~ lincome + lconsumption, data = Lutkephol, lags = 1, trend = "const", test = "all", crit = 0.05)
A$bh.crit.val


