
<!-- README.md is generated from README.Rmd. Please edit that file -->

# pvurt

<!-- badges: start -->

<!-- badges: end -->

This package has been created as part of my master thesis. The package
has only one user facing function: <tt>computePVale()</tt>, which
approximates p-value (cumulative probability) for augmented
Dickey-Fuller unit root test based on the distribution of the test
statistics. The approximation is based on pre-trained generalized
additive logistic model or logistic regression with polynomial terms.

The function can take either a numeric value or the output from the four
functions of the following three packages: <tt>adfTest()</tt> and
<tt>unitrootTest()</tt> from <tt>fUnitRoots</tt>, <tt>adf.test()</tt>
from <tt>tseries</tt> and, lastly, <tt>ur.df()</tt> from <tt>urca</tt>.
If output from any of the function is passed on, <tt>computePVale()</tt>
appends its approximation result to the original output.

## Installation

The package has not been released in [CRAN](https://CRAN.R-project.org).
To install use:

``` r
devtools::install_github("mlincon/pvurt")
```

Since the package must be complied, ensure that <tt>Rtools.exe</tt> is
installed beforehand.

## Example

``` r
library(pvurt)

y <- arima.sim(model = list(order = c(0, 1, 0)), n = 100)

# Test type: with drift and trend
# package: fUnitRoots
library(fUnitRoots)
#> Loading required package: timeDate
#> Loading required package: timeSeries
#> Loading required package: fBasics
computePValue(adfTest(y, lags = 3, type = "ct"))
#> 
#> Title:
#>  Augmented Dickey-Fuller Test
#> 
#> Test Results:
#>   PARAMETER:
#>     Lag Order: 3
#>   STATISTIC:
#>     Dickey-Fuller: -1.1702
#>   P VALUE:
#>     t: 0.9078 
#>     pvurt: 0.9099 
#> 
#> Description:
#>  Fri Sep 13 11:18:55 2019 by user: User

computePValue(unitrootTest(y, lags = 3, type = "ct"))
#> 
#> Title:
#>  Augmented Dickey-Fuller Test
#> 
#> Test Results:
#>   PARAMETER:
#>     Lag Order: 3
#>   STATISTIC:
#>     DF: -1.1702
#>   P VALUE:
#>     t: 0.9107 
#>     n: 0.9851 
#>     pvurt: 0.9099 
#> 
#> Description:
#>  Fri Sep 13 11:18:55 2019 by user: User

# package: urca
library(urca)
#> 
#> Attaching package: 'urca'
#> The following objects are masked from 'package:fUnitRoots':
#> 
#>     punitroot, qunitroot, unitrootTable
computePValue(ur.df(y, lags = 3, type = "trend"))
#> 
#> ############################################################### 
#> # Augmented Dickey-Fuller Test Unit Root / Cointegration Test # 
#> ############################################################### 
#> 
#> The value of the test statistic is: -1.1702 0.5717 0.8238 
#> Approximated P value is (pvurt):    0.9099 NA NA
# print summary
summary(computePValue(ur.df(y, lags = 3, type = "trend")))
#> 
#> ############################################### 
#> # Augmented Dickey-Fuller Test Unit Root Test # 
#> ############################################### 
#> 
#> Test regression trend 
#> 
#> 
#> Call:
#> lm(formula = z.diff ~ z.lag.1 + 1 + tt + z.diff.lag)
#> 
#> Residuals:
#>      Min       1Q   Median       3Q      Max 
#> -2.51378 -0.61868  0.01927  0.54949  2.15283 
#> 
#> Coefficients:
#>              Estimate Std. Error t value Pr(>|t|)
#> (Intercept)  0.017396   0.205271   0.085    0.933
#> z.lag.1     -0.036016   0.030778  -1.170    0.245
#> tt           0.001004   0.003939   0.255    0.799
#> z.diff.lag1  0.021370   0.105427   0.203    0.840
#> z.diff.lag2  0.051329   0.105211   0.488    0.627
#> z.diff.lag3 -0.013270   0.105714  -0.126    0.900
#> 
#> Residual standard error: 0.879 on 91 degrees of freedom
#> Multiple R-squared:  0.01995,    Adjusted R-squared:  -0.0339 
#> F-statistic: 0.3705 on 5 and 91 DF,  p-value: 0.8677
#> 
#> 
#> Value of test-statistic is:     -1.1702 0.5717 0.8238 
#> Approximated P value is (pvurt): 0.9099 NA NA 
#> 
#> Critical values for test statistics: 
#>       1pct  5pct 10pct
#> tau3 -3.99 -3.43 -3.13
#> phi2  6.22  4.75  4.07
#> phi3  8.43  6.49  5.47

# package: tseries
library(tseries)
computePValue(adf.test(y, alternative = "stationary", k = 3))
#> 
#>  Augmented Dickey-Fuller Test
#> 
#> data:  y
#> Dickey-Fuller = -1.17, Lag order = 3, p-value = 0.9078, pvurt =
#> 0.9099
#> alternative hypothesis: stationary

# no packages
tStat <- -2.239
sampleSize <- 100
computePValue(tStat, n = sampleSize, model = "gam", type = "ct")
#> [1] 0.4626582
```
