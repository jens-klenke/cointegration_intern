The folder contains all MatLab Files used for the large and small Sample Simulations.
The subfolders are edited (to speed up calculations) versions of the Packages from LeSage.

main_asy.m: 	Batch that performs the simulations for the asymptotic test
main_boot.m: 	Batch that performs simulations for the bootstrap test
boottest.m: 	Bootstrap Cointegration test
asytest.m: 	Asyptotic Cointegration test
asy.m: 		Generate Asymptotic distributions.
gen_ar1.m: 	Subfunction

Matlab Workspaces:
critical_values.mat:	Array of Critical Values
NullDistributions.mat:	Distributions of underlying tests under the Null to calculate p-values

Note that the simulations in the paper are based on 100,000 replications / draws from the Null Distribution.
The Matrix NullDistribution.mat here contains only 10,000 elements (i.e. every 10th element) to save memory. 
For practical purposes this should be sufficient.

Note also that running the code requires the `distributed computing' package. If you do not have access to this 
package replace all instances of `parfor' in the code by `for', i.e. line 88 in boottest.m and line 196 in asy.m 

Also, comment out lines 14-18 in main_boot.m and 39-43 in asy.m i.e.

pool=matlabpool('size')
if pool>0
    matlabpool close
end
matlabpool open 4

To run boottests.m or asytests.m as stand-alone functions (e.g., to conduct an empirical analysis) 
the paths to the subfolders and the main folder need to be added using the addpath command, as in main_asy.m and 
main_boot.m

-------------------------------------------------------




