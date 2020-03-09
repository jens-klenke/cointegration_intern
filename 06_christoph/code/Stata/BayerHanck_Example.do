*-------------------------------------------------------------------------
* Usage example of Bayer-Hanck (2009) Test.
* Copy bayerhanck.ado and NullDistr.dta in 
* the PERSONAL Folder set in:
 
sysdir set PERSONAL "C:\Documents and Settings\Bayer\My Documents\Stata"
*sysdir set PERSONAL "C:\Users\bayer\Documents\Stata"

* then type 
*
* bayerhanck LHSvar, rhs(RHSvarlist) [trend(none | constant | trend) lags(num) crit(1 | 5 | 10)]
*
* for the test. The test provides P-Values for underlying tests and the 
* test statistics for the combined Fisher type test.
* For further documentation of trend options see help file for dfuller.ado
* Defaults are: trend=constant, lags=1, crit=5
*
* Version 0.9
*-------------------------------------------------------------------------
clear all
set more off
set matsize 10000

local rep=2
mat def testPower=J(`rep',2,999)
mat def testSize=J(`rep',2,999)

forv z=1/`rep' {
	clear
	
	*qui{
		set obs 700	
		gen dx=rnormal()
		gen x=sum(dx)+50
		forv j=1/5 {
			gen dz`j'=rnormal()
			gen z`j'=sum(dz`j')
		}
		gen T=_n
		tsset T
		gen u=rnormal()
		replace u=u+0.9*l.u if T>1
		gen y=x+u
		drop if T<200
		replace x=x+T
		
	

		di "Power Example"
		bayerhanck x, rhs(y) trend(trend) lags(1)
		mat testPower[`z',1]=`e(EJ)'
		mat testPower[`z',2]=`e(BECREJ)'
		di "Size Example"
		bayerhanck x, rhs(z*) trend(trend) lags(1) crit(10)
		mat testSize[`z',1]=`e(EJ)'
		mat testSize[`z',2]=`e(BECREJ)'
	*}
}
ereturn list
svmat testPower
svmat testSize
keep test*
drop if testPower1==.
count if testSize1>e(CRIT_EJ)
count if testSize2>e(CRIT_BECREJ)
count if testPower1>e(CRIT_EJ)
count if testPower2>e(CRIT_BECREJ)
