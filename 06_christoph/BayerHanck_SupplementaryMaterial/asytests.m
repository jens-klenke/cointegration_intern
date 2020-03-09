function [pEngleGranger, pJohansen, pErrCorr, pBoswijk, FisherBECR, FisherBJ, ...
    FisherBE, FisherEJ, FisherBJE, FisherBECRJ, FisherBECRJE, MinBECR, MinEJ]...
    =asytests(X, nlag, nlag_r,NullDistr)
% ASYTESTS performs the test based on asymptotic 
% critical values for the underlying tests. Critical values are supplied in 
% Critical_values.mat. 
% -------------------------------------------------------------------------
% USAGE: [pEngleGranger, pJohansen, pErrCorr, pBoswijk, FisherBECR, FisherBJ, ...
%    FisherBE, FisherEJ, FisherBJE, FisherBECRJ, FisherBECRJE, MinBECR, MinEJ]...
%    =asytests(X, nlag, nlag_r)
%--------------------------------------------------------------------------
% It requires as input variables
% X:        Matrix of T x p observations, where the first column is the left hand
%           variable in the Engle Granger test.
% nlag:     number of lags used in Johansen, Boswijk, and Banerjee test
% nlag_r:   number of lags used in the Engle Granger Test
%--------------------------------------------------------------------------
% OPTIONAL usage:  ASYTESTS(X, nlag, nlag_r, NullDistr)
% NullDistr: Array containing Simulated Distributions under the null
%--------------------------------------------------------------------------
% The function returns P-values for the underlying tests (first 4 output 
% arguments) and then test statistics in the following order:
%             1:  EngleGranger (p-value)         
%             2:  Johansen l-max (p-value)         
%             3:  ErrCorr (p-value)         
%             4:  Boswijk (p-value)
%             5:  FisherBECR;
%             6:  FisherBJ
%             7:  FisherBE
%             8:  FisherEJ
%             9:  FisherBJE
%             10: FisherBECRJ
%             11: FisherBECRJE
%             12: MinBECR
%             13: MinEJ
%
% The appropriate critical values for 5-13 are given in Tables 1 and A.1
% in the extended working paper version.
%
% Date: 02 Apr 2009
%--------------------------------------------------------------------------
% If the Array NullDistr is supplied, it must be organized as follows 
% (default: use the one supplied with the test codes)
% NullDistr cols  1: EngleGranger,      CV dim 3: nvars-1
%                 2: Johansen,          CV dim 4: deterministics (trend+1)
%                 3: ErrCorr,
%                 4: Boswijk
%--------------------------------------------------------------------------
%%
if nargin==3
    load NullDistributions
end
%%
trend=1; % Note that here trend has a different notation: -1 no trend, 0 demean, 1 linear trend
r=0; %Number of co-integration vectors under the null
nvars=size(X,2); %Number of variables
T=size(X,1)-nlag-1; %Time-series length after lags

%% Step 1: Run underlying tests and obtain p-values
%% Johansen Test
jres = johansen(X,trend,nlag); %Johansen Test
lmax=jres.lr2(1+r);
pJohansen=mean(lmax<NullDistr(:,2,nvars-1,trend+2));

%% Engle-Granger Test
results = cadf(X(:,1),X(:,2:end),trend,nlag_r,1); 
cad=results.adf;
pEngleGranger=mean(cad>NullDistr(:,1,nvars-1,trend+2));

%% Boswijk Test & Banerjee Test
[WaldStat,ErrCorrStat]=boswijk(X(:,1),X(:,2:end),trend,nlag);
pBoswijk=mean(WaldStat<NullDistr(:,4,nvars-1,trend+2));
pErrCorr=mean(ErrCorrStat>NullDistr(:,3,nvars-1,trend+2));

%% Calculate Vector of Log P-Values
logP=log([pEngleGranger, pJohansen, pErrCorr, pBoswijk]);
%% Define Fisher stats
FisherBECR=-2*sum(logP([3 4]));
FisherBJ=-2*sum(logP([2 4]));
FisherBE=-2*sum(logP([1 4]));
FisherEJ=-2*sum(logP([1 2]));
FisherBJE=-2*sum(logP([1 2 4]));
FisherBECRJ=-2*sum(logP([2 3 4]));
FisherBECRJE=-2*sum(logP([1 2 3 4]));

MinBECR=min(pErrCorr, pBoswijk);
MinEJ=min(pEngleGranger, pJohansen);