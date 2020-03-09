function [Wstat,tstat]=boswijk(y,z,p,nlag)

% PURPOSE: Compute Boswijk (1994, JECM)
% ----------------------------------------------------------------
% USAGE: Wstat = boswijk(y,x,p,nlag)
% where: y = dependent variable time-series vector
%        z = explanatory variables matrix
%             p = order of time polynomial in the null-hypothesis
%                 p = -1, no deterministic part
%                 p =  0, for constant term
%                 p =  1, for constant plus time-trend
%                 p >  1, for higher order polynomial
%     nlag = # of lagged changes of the residuals to include in regression
% ----------------------------------------------------------------
% RETURNS: Wald statistic from the Boswijk test
%       and t-Statistics from Banerjee t-Test (ECR)
% ---------------------------------------------------------------- 
% References: Peter Boswijk (1994) 'Testing for an Unstable Root 
% in Conditional and Unconditional Error Correction Models', 
% Biometrika, Volume 63, pp. 37-60.
%
% written by:
% Christian Bayer and Christoph Hanck
% Universität Bonn 
% and Universiteit Maastricht
% christian.bayer@unibocconi.it
% c.hanck@ke.unimaas.nl
%% Demean
y = detrend(y,p);
z = detrend(z,p);



%% Generate Regressor Matrix
aux=[y, z];
Xlag = lag(aux,1); % Lagged Levels
Xlag = trimr(Xlag,1,0) ;
Ydif = diff(y);

aux = diff(Xlag);
W   = diff(z);
if nlag > 0
    for runlag=1:nlag
        W=[W(2:end,:),aux];
        aux = lag(aux,1); % Lagged FD
        aux = trimr(aux,1,0);
    end
end
Xlag = trimr(Xlag,nlag,0);
%%
Ydif = trimr(Ydif,nlag,0);
[nobs, nvar]=size(Xlag);

%M = eye(nobs)-W*inv(W'*W)*W';
resX = Xlag-W*(W\Xlag);
resY = Ydif-W*(W\Ydif);


delta = (resX)\(resY);

u = resY-resX*delta;
omega=(nobs-1)/(nobs-(nlag+2)*nvar+1)* var(u);
Vinv=resX'*resX/omega;
Wstat=delta'*Vinv*delta;

Vinv2=inv(Vinv);
tstat=delta(1)/sqrt(Vinv2(1,1));





