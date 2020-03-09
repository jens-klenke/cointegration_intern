%-----------------------------------------------------------------------
% This is the Batch File that runs the small sample simulations for the
% asymptotic Co-Integration test
% Date: 03 Apr 2009
%-----------------------------------------------------------------------
%% Make use of the LeSage Toolbox for some Co-Integration tests
addpath('coint')
addpath('util')

%% Set Program Parameters,
TL=[50 75 100 150 200]% Number of Years
p=2; % Number of Variables
B=10000; % Number of bootstrap replications
H=5000; %Number of monte-carlo experiments
critlist=[0.01 0.05 0.10];
crittype=2
h=1;
discard=0;
EX=6
%% Define DGP
%DGP (A)
cor=0; % Short-run correlation of shocks
s2=1; %Variance of x2
%DGP(B)
Gamma_coint=[-0.1 0.1;0.1 -0.1]; % Specification from Swenson (just that he uses 6 variables w/ 2 cointegrated)
[TEMP1,TEMP2]=eig([1 cor;cor s2]);
Sig=TEMP1*sqrt(TEMP2); % Spectral decomposition
Gamm=[0.2 0; 0 0.2]; % Short run VAR Coeff
%DGP(C)
beta=1; % Beta-factor
alpha=[0 0]; % intercept
endog=0.5; % endogeneity
C=[0, 15]; %Local to Asymptotic Parameter

tic

for j=1:length(TL) % Different time-length
    T=TL(j);
    clear v1 v2 v3 v4
    
    % pre-define Matrices to collect p-values
    pEngleGranger=ones(H,length(C),EX);
    pJohansen=ones(H,length(C),EX);
    pErrCorr=ones(H,length(C),EX);
    pBoswijk=ones(H,length(C),EX);
    pFisherBECR=ones(H,length(C),EX);
    pFisherBJ=ones(H,length(C),EX);
    pFisherBE=ones(H,length(C),EX);
    pFisherEJ=ones(H,length(C),EX);
    pFisherBJE=ones(H,length(C),EX);
    pFisherBECRJ=ones(H,length(C),EX);
    pFisherBECRJE=ones(H,length(C),EX);
    pMinBECR=ones(H,length(C),EX);
    pMinEJ=ones(H,length(C),EX);
    
    % Load Asymptotic Distributions under the Null
    load NullDistributions
    for h=1:H %Replications of the Monte Carlo experiments
         % Set lag length to its true value
        nlag=0; %Lag length for all Error-Correction Tests
        nlag_r=0; % Lag Length for Engle-Granger Test
        u=normrnd(0,1,T+30,p); % Draw Shocks 
        load critical_values % Load Simulated Critical Values for asymptotic tests
        for s=1:length(C)  % Loop over local to asymptotic parameters
            for ex=1:EX-2 % Loop over various DGP(A) specs                
                delta=0.25*(ex-1); %Pesavento R2
                Omega=[1 sqrt(delta); sqrt(delta) 1]; %Correlation Matrix
                % Generate correlated shocks
                v=u*chol(Omega);
                z1=cumsum(v(:,1)); % first I(1) variable
                if s==1 %Size
                    z2=cumsum(v(:,2)); % second Variable
                else %Power
                    z2=gen_ar1(T+30,1-C(s)/T,v(:,2)); % second Variable
                end
                X=[beta*z1+z2,z1]; %Possibly cointegrated system
                X(1:30,:)=[]; % Drop initial 30 observation
                % Run Test
                [pEngleGranger(h,s,ex), pJohansen(h,s,ex), pErrCorr(h,s,ex), pBoswijk(h,s,ex), FisherBECR(h,s,ex), FisherBJ(h,s,ex), ...
                    FisherBE(h,s,ex), FisherEJ(h,s,ex), FisherBJE(h,s,ex), FisherBECRJ(h,s,ex), FisherBECRJE(h,s,ex), ...
                    MinBECR(h,s,ex), MinEJ(h,s,ex)]=asytests(X, nlag, nlag_r,NullDistr);    
            end
        end
        
        %% DGP (B)
        ex=EX-1; %DGP(B)
        
        nlag=1;%Lag length for all Error-Correction Tests
        nlag_r=1;% Lag Length for Engle-Granger Test
        X=zeros(T+30,p); %Predefine Matrix to hold obs
        for s=1:2; %Size /Power
            X(1,:)=u(1,:); % Initialize ECR
            dx=zeros(1,2); % initialize first differences
            if s==1
                for t=2:T+30 %Generate simulated observations, Size
                    X(t,:)=X(t-1,:)+dx*Gamm+u(t,:)*Sig;
                    dx=X(t,:)-X(t-1,:);
                end
            else
                for t=2:T+30 %Generate simulated observations, Power
                    X(t,:)=X(t-1,:)*(eye(2)+Gamma_coint)+dx*Gamm+u(t,:)*Sig;
                    dx=X(t,:)-X(t-1,:);
                end
            end
            X(1:30,:)=[]; %Drop initial Observations
            % Run test
            [pEngleGranger(h,s,ex), pJohansen(h,s,ex), pErrCorr(h,s,ex), pBoswijk(h,s,ex), FisherBECR(h,s,ex), FisherBJ(h,s,ex), ...
                FisherBE(h,s,ex), FisherEJ(h,s,ex), FisherBJE(h,s,ex), FisherBECRJ(h,s,ex), FisherBECRJE(h,s,ex), ...
                MinBECR(h,s,ex), MinEJ(h,s,ex)]=asytests(X, nlag, nlag_r,NullDistr);
            
        end
        
        %% DGP (C)
        nlag=0; %Lag length for all Error-Correction Tests
        nlag_r=0; % Lag Length for Engle-Granger Test
        v1=cumsum(u(:,1)); %Generate I(1) Variable
        ex=EX; %DGP(C)
        for s=1:length(C); %Loop over Local to Asymptotic
            v3=gen_ar1(T+30,1-C(s)/T,u(:,2)); % Generate second shock
            X=[alpha(1)+beta*v1+v3, alpha(2)+v1+endog*v3]; %Possibly cointegrated system
            X(1:30,:)=[]; %Drop initial Observations
            [pEngleGranger(h,s,ex), pJohansen(h,s,ex), pErrCorr(h,s,ex), pBoswijk(h,s,ex), FisherBECR(h,s,ex), FisherBJ(h,s,ex), ...
                FisherBE(h,s,ex), FisherEJ(h,s,ex), FisherBJE(h,s,ex), FisherBECRJ(h,s,ex), FisherBECRJE(h,s,ex), ...
                MinBECR(h,s,ex), MinEJ(h,s,ex)]=asytests(X, nlag, nlag_r,NullDistr);
        end

    end
    trend=2; % Variable to select case (iii) from the asymptotic critical values
    nvars=size(X,2); % number of variables for the asymptotic critical values
    crit=critlist(crittype); % Nominal level
    clc
    toc
    
    head={' EG ', ' J','ErrCorr', 'Wald',' Fisher BECR ',' Fisher BJ ',' Fisher BE ',' Fisher EJ ',' Fisher BJE ',' Fisher BECRJ ',' Fisher BECRJE ', 'Min BECR', 'Min EJ','naive EJ', 'naive BECR'};
    disp(['rejection rates at ' num2str(100*crit) '% lvl'])
            disp('-------------------------------------------')
        disp(['DGP (A)'])
        disp('-------------------------------------------')
    for s=1:length(C)
        if s==1
            disp('-------------------------------------------')
            disp(['Size'])
            disp('-------------------------------------------')
        else
            disp('-------------------------------------------')
            disp(['Power'])
            disp('-------------------------------------------')
        end
        for ex=1:EX-2
            disp('-------------------------------------------')
            disp(['delta = ' num2str(0.25*(ex-1))])
            disp('-------------------------------------------')
            % Critical values loaded from critval.mat, Array CV is
            % constructed as:
            %  CV: rows:  1:  Boswijk           CV cols : nvars-1
            %             2:  ErrCorr           CV dim 3: deterministics (trend+1)
            %             3:  Johansen          CV dim 4: LvL 1% 5% 10%
            %             4:  EngleGranger
            %             5:  FisherBECR;
            %             6:  FisherBJ
            %             7:  FisherBE
            %             8:  FisherEJ
            %             9:  FisherBJE
            %             10: FisherBECRJ
            %             11: FisherBECRJE
            %             12: InvNormBJE
            %             13: MinBECR
            %             14: MinEJ
            
            %% Count number of rejections
            aux=[pEngleGranger(:,s,ex)<crit, pJohansen(:,s,ex)<crit, pErrCorr(:,s,ex)<crit, pBoswijk(:,s,ex)<crit,...
                FisherBECR(:,s,ex)>CV(5,nvars-1,trend+1,crittype), FisherBJ(:,s,ex)>CV(6,nvars-1,trend+1,crittype), ...
                FisherBE(:,s,ex)>CV(7,nvars-1,trend+1,crittype), FisherEJ(:,s,ex)>CV(8,nvars-1,trend+1,crittype),...
                FisherBJE(:,s,ex)>CV(9,nvars-1,trend+1,crittype), FisherBECRJ(:,s,ex)>CV(10,nvars-1,trend+1,crittype),...
                FisherBECRJE(:,s,ex)>CV(11,nvars-1,trend+1,crittype), MinBECR(:,s,ex)<CV(13,nvars-1,trend+1,crittype), ...
                MinEJ(:,s,ex)<CV(14,nvars-1,trend+1,crittype),min(pEngleGranger(:,s,ex), pJohansen(:,s,ex))<crit, ...
                min(pErrCorr(:,s,ex), pBoswijk(:,s,ex))<crit];
            table=head;
            base=(sum(aux)');
            for ttt=1:length(base)
                table{2,ttt}=base(ttt)/h;
            end
            Collect(j+(ex-1)*length(TL),:,s)=[T base'/h];
            disp(table)
            
        end
        
        ex=EX-1;
        disp('-------------------------------------------')
        disp(['DGP (B)'])
        disp('-------------------------------------------')
        %% Count number of rejections
        aux=[pEngleGranger(:,s,ex)<crit, pJohansen(:,s,ex)<crit, pErrCorr(:,s,ex)<crit, pBoswijk(:,s,ex)<crit,...
            FisherBECR(:,s,ex)>CV(5,nvars-1,trend+1,crittype), FisherBJ(:,s,ex)>CV(6,nvars-1,trend+1,crittype), ...
            FisherBE(:,s,ex)>CV(7,nvars-1,trend+1,crittype), FisherEJ(:,s,ex)>CV(8,nvars-1,trend+1,crittype),...
            FisherBJE(:,s,ex)>CV(9,nvars-1,trend+1,crittype), FisherBECRJ(:,s,ex)>CV(10,nvars-1,trend+1,crittype),...
            FisherBECRJE(:,s,ex)>CV(11,nvars-1,trend+1,crittype), MinBECR(:,s,ex)<CV(13,nvars-1,trend+1,crittype), ...
            MinEJ(:,s,ex)<CV(14,nvars-1,trend+1,crittype),min(pEngleGranger(:,s,ex), pJohansen(:,s,ex))<crit, ...
            min(pErrCorr(:,s,ex), pBoswijk(:,s,ex))<crit];
        table=head;
        base=(sum(aux)');
        for ttt=1:length(base)
            table{2,ttt}=base(ttt)/h;
        end
        disp(table)
        Collect(j+(ex-1)*length(TL),:,s)=[T base'/h];
        ex=EX;
        disp('-------------------------------------------')
        disp(['DGP(C)'])
        disp('-------------------------------------------')
        %% Count number of rejections
        aux=[pEngleGranger(:,s,ex)<crit, pJohansen(:,s,ex)<crit, pErrCorr(:,s,ex)<crit, pBoswijk(:,s,ex)<crit,...
            FisherBECR(:,s,ex)>CV(5,nvars-1,trend+1,crittype), FisherBJ(:,s,ex)>CV(6,nvars-1,trend+1,crittype), ...
            FisherBE(:,s,ex)>CV(7,nvars-1,trend+1,crittype), FisherEJ(:,s,ex)>CV(8,nvars-1,trend+1,crittype),...
            FisherBJE(:,s,ex)>CV(9,nvars-1,trend+1,crittype), FisherBECRJ(:,s,ex)>CV(10,nvars-1,trend+1,crittype),...
            FisherBECRJE(:,s,ex)>CV(11,nvars-1,trend+1,crittype), MinBECR(:,s,ex)<CV(13,nvars-1,trend+1,crittype), ...
            MinEJ(:,s,ex)<CV(14,nvars-1,trend+1,crittype),min(pEngleGranger(:,s,ex), pJohansen(:,s,ex))<crit, ...
            min(pErrCorr(:,s,ex), pBoswijk(:,s,ex))<crit];
        base=(sum(aux)');
        table=head;
        for ttt=1:length(base)
            table{2,ttt}=base(ttt)/h;
        end
        disp(table)
        Collect(j+(ex-1)*length(TL),:,s)=[T base'/h];
    end
    %% Save Results
    savename=['results_SMALL_ASYTESTS_T_all_tests_varlag_3=' num2str(T)];
    save(savename)
end
%% Generate Table of all results
head={'T', ' EG ', ' J','ErrCorr', 'Wald',' Fisher BECR ',' Fisher BJ ',' Fisher BE ',' Fisher EJ ',' Fisher BJE ',' Fisher BECRJ ',' Fisher BECRJE ', 'Min BECR', 'Min EJ','naive EJ', 'naive BECR'};
matrix2latex(Collect(:,:,1),'SmallSample_Size.tex','columnlabels', head)
for s=2:length(C)
    matrix2latex(Collect(:,:,s),['SmallSample_Power_C=', num2str(C(s)), '.tex'],'columnlabels', head)
end


%%
save results2
