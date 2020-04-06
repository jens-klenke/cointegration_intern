library(devtools)
load_all()


A <- bayerhanck(linvestment ~ lincome + lconsumption, data = Lutkephol, lags = 1, trend = "const", test = "all", crit = 0.05)
A$bh.crit.val


test1 <- bayerhanck(linvestment ~ lincome + lconsumption, Lutkephol, lags = 3, trend = "none", test = "eg-j", crit = 0.10)
plot(test1)


banerjee(linvestment ~ lincome + lconsumption, Lutkephol, lags = 1, trend = "none") #differ: stata(-2.1739) R(-1.3212)

banerjee(linvestment ~ lincome + lconsumption, Lutkephol, lags = 1, trend = "const")
