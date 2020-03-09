function [pEngleGranger, pJohansen, pErrCorr, pBoswijk, pFisherBECR, pFisherBJ,...
    pFisherBE, pFisherEJ, pFisherBJE, pFisherBECRJ, pFisherBECRJE, pMinBECR, pMinEJ]...
    =boottests(X, nlag, nlag_r,B)
% BOOTTESTS performs the bootstrap version of the test.
% -------------------------------------------------------------------------
% USAGE: [pEngleGranger, pJohansen, pErrCorr, pBoswijk, pFisherBECR, pFisherBJ,...
%    pFisherBE, pFisherEJ, pFisherBJE, pFisherBECRJ, pFisherBECRJE, pMinBECR, pMinEJ]...
%    =boottests(X, nlag, nlag_r)
% -------------------------------------------------------------------------
% It requires as input variables
% X:        Matrix of T x p observations, where the first column is the left hand
%           variable in the Engle Granger test.
% nlag:     number of lags used in Johansen, Boswijk, and Banerjee test
% nlag_r:   number of lags used in the Engle Granger Test
% -------------------------------------------------------------------------
% Optional USAGE: boottests(X, nlag, nlag_r,B)
% B:        number of Bootstrap Replications
%--------------------------------------------------------------------------
% The function returns bootstraped P-values for the various test in the 
% following order:
%             1:  EngleGranger          
%             2:  Johansen l-max          
%             3:  ErrCorr          
%             4:  Boswijk 
%             5:  FisherBECR;
%             6:  FisherBJ
%             7:  FisherBE
%             8:  FisherEJ
%             9:  FisherBJE
%             10: FisherBECRJ
%             11: FisherBECRJE
%             12: MinBECR
%             13: MinEJ
% NOTE: The test returns NaN if explosive root found.
% Date: 02 Apr 2009
%--------------------------------------------------------------------------
%%
if nargin==3
    B=10000;
end
%%
trend=1; % Note that here trend has a different notation: -1 no trend, 0 demean, 1 linear trend
r=0; %Number of co-integration vectors under the null
nvars=size(X,2); %Number of variables
T=size(X,1)-nlag-1; %Time-series length after lags

%% Swenson Algorithm to Bootstrap an ECM
% Step 1: Estimation as a VAR
results2=vare(X, nlag+1);
for j=1:nvars
    eps(:,j)=results2(j).resid;
    Coeff2(:,j)=results2(j).beta;
end
Coeff3=zeros(nlag*nvars+1,nvars);
for j=1:nvars
    for lag=nlag:-1:1
        Coeff3((j-1)*(nlag)+lag,:)=Coeff2((j-1)*(nlag+1)+lag+1,:)+Coeff3((j-1)*(nlag)+lag+1,:);
    end
end
Coeff=-Coeff3;
Coeff(end,:)=Coeff2(end,:);

    

%% Step 2: Check for explosive root

for j=1:nlag
    G(1:nvars,1+(j-1)*nvars:j*nvars)=Coeff(nlag-j+1+(0:nvars-1).*nlag,:)';
end
test=2;
if nlag >0
    Temp=polymroot([G'; eye(nvars)]);
    test=min(abs(Temp));
end
if test>1
    %% Johansen Test
    jres = johansen(X,trend,nlag); %Johansen Test
    lmax=jres.lr2(1+r);
    ecvectors = jres.evec;
    %% Engle-Granger Test
    results = cadf(X(:,1),X(:,2:end),trend,nlag_r,1); % Need to check 'trend' % last parameter is 1 if crit vals are calculated
    cad=results.adf;
        
    %% Boswijk & Banerjee Test
    [WaldStat,ErrCorrStat]=boswijk(X(:,1),X(:,2:end),trend,nlag);
    
    
    %% Generate Bootstrap replication
    %initialize test statistics vectors
    st_lmax=zeros(B,1);
    st_cad=zeros(B,1);
    st_Wald=zeros(B,1);
    st_ECR=zeros(B,1);
    
    ZZ=X(1:nlag+1,:);
    warning('off','stats:regress:NoConst')
    parfor b=1:B
        %% Step 3: Pseudo-Observations
        index=ceil(rand(length(eps),1)*length(eps)); % Draw index numbers
        s_err=eps(index,:); % Generate pseudo-shock vector 
        Y=DGP(ZZ,Coeff, nlag, s_err); % generate pseudo-observations from Shocks
        %% Step 4: Calculate various test statistics
        % Johansen test
        jres = johansen(Y,trend,nlag);
        st_lmax(b)=jres.lr2(1+r);
    
        % Engle-Granger
        results = cadf(Y(:,1),Y(:,2:end),trend,nlag_r,1); % last parameter is 1 if crit vals are calculated
        st_cad(b)=results.adf;
    
        % boswijk & Banerjee Test
        [st_Wald(b),st_ECR(b)]=boswijk(Y(:,1),Y(:,2:end),trend,nlag);
    
    end
    warning('off','stats:regress:NoConst')
    %% Step 5: Generate pseudo p-values
    ST=[st_lmax st_cad st_Wald st_ECR];
    
    ST=sortrows(ST,1); %Johansen Test
    ST=[ST, ((1:-1/B:1/B))'];
    
    ST=sortrows(ST,2); %Engle Granger Test
    ST=[ST, ((1/B:1/B:1))'];
    
    ST=sortrows(ST,3); % Boswijk Test
    ST=[ST ((1:-1/B:1/B))'];
    
    ST=sortrows(ST,4); %Banerjee Test
    ST=[ST ((1/B:1/B:1))'];
    
    % Read P-value vectors from ST
    pB_J=ST(:,5);
    pB_EG=ST(:,6);
    pB_Wald=ST(:,7);
    pB_ECR=ST(:,8);
    
    % Pseudo P-values for the individual tests
    pJohansen=mean(st_lmax>=lmax)+.0000000000000000000000000000001;
    pEngleGranger=mean(st_cad<=cad)+.00000000000000000000000000001;
    pErrCorr=mean(st_ECR<=ErrCorrStat)+.00000000000000000000000000001;
    pBoswijk=mean(st_Wald>=WaldStat)+.00000000000000000000000000001;
    
    
    %% Step 6: Distribution of Combination Test statistics
    % Fisher Tests
    logP=log([pEngleGranger, pJohansen, pErrCorr, pBoswijk]);
    FisherBECR=-2*sum(logP([3 4]));
    FisherBJ=-2*sum(logP([2 4]));
    FisherBE=-2*sum(logP([1 4]));
    FisherEJ=-2*sum(logP([1 2]));
    FisherBJE=-2*sum(logP([1 2 4]));
    FisherBECRJ=-2*sum(logP([2 3 4]));
    FisherBECRJE=-2*sum(logP([1 2 3 4]));
    
    % Min Tests
    MinBECR=min(pErrCorr, pBoswijk);
    MinEJ=min(pEngleGranger, pJohansen);
    
    %% Step 7: Calculate Combination Test statistics for actual tests (on real data) 
    logP_B=log([pB_EG, pB_J, pB_ECR, pB_Wald]);
    st_FisherBECR=-2*sum(logP_B(:,[3 4]),2);
    st_FisherBJ=-2*sum(logP_B(:,[2 4]),2);
    st_FisherBE=-2*sum(logP_B(:,[1 4]),2);
    st_FisherEJ=-2*sum(logP_B(:,[1 2]),2);
    st_FisherBJE=-2*sum(logP_B(:,[1 2 4]),2);
    st_FisherBECRJ=-2*sum(logP_B(:,[2 3 4]),2);
    st_FisherBECRJE=-2*sum(logP_B(:,[1 2 3 4]),2);
    
    st_MinBECR=min(pB_ECR, pB_Wald);
    st_MinEJ=min(pB_EG, pB_J);
    
    %% Step 8: Obtain P-Values comparing statistics from Step 7 with distributions from step 6
    pFisherBECR=mean(st_FisherBECR>=FisherBECR);
    pFisherBJ=mean(st_FisherBJ>=FisherBJ);
    pFisherBE=mean(st_FisherBE>=FisherBE);
    pFisherEJ=mean(st_FisherEJ>=FisherEJ);
    pFisherBJE=mean(st_FisherBJE>=FisherBJE);
    pFisherBECRJ=mean(st_FisherBECRJ>=FisherBECRJ);
    pFisherBECRJE=mean(st_FisherBECRJE>=FisherBECRJE);
    
    pMinBECR=mean(st_MinBECR<=MinBECR);
    pMinEJ=mean(st_MinEJ<=MinEJ);
else
    %% Return NaN if explosive root
    pFisherBECR=NaN;
    pFisherBJ=NaN;
    pFisherBE=NaN;
    pFisherEJ=NaN;
    pFisherBJE=NaN;
    pFisherBECRJ=NaN;
    pFisherBECRJE=NaN;
    
    pMinBECR=NaN;
    pMinEJ=NaN;
end

function Y=DGP(START,Coeff, nlag, err)
% This subfunction to boottest.m calculates The realisations of an
% Error-Correction Process from the Coefficient Matrix COEFF, the error
% Matrix ERR, the Starting values START, and number of lags NLAG
%--------------------------------------------------------------------------
% Coefficients must be ordered as: differences variable1 lag1 ...
% variable1 lag_n, variable2 lag_1, ..., variable_J lag_n, levels variable1
% ...variable_J, intercept; THE DIFFERENT VARIABLES ARE IN COLUMNS
%--------------------------------------------------------------------------
T=size(err,1);
nvar=size(err,2);
Y=zeros(T+nlag+1,nvar);
Y(1:nlag+1,:)=START;
dY=Y(2:nlag+1,:)-Y(1:nlag,:);
dx=zeros(T,nvar);
dY=dY(end:-1:1,:); % sort in reverse order, newest on top
Fix=Coeff(end,:);
Coeff(end,:)=[];
if size(Coeff,1)>=1
    for t=1:T
        dx(t,:)= dY(:)'*Coeff+Fix+err(t,:);
        dY(2:end,:)=dY(1:end-1,:);
        dY(1,:)=dx(t,:);
    end
    dx(1,:)=dx(1,:)+Y(nlag+1,:);
    X=cumsum(dx);
    Y=[Y(1:nlag+1,:);X];
else
    for j=1:nvar
        err(:,j)=err(:,j)+Fix(j);
    end
    err(1,:)=err(1,:)+Y(nlag+1,:);
    X=cumsum(err);
    Y=[Y(1:nlag+1,:);X];
end