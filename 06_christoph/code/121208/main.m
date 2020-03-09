addpath(genpath('../packages'))
pool=matlabpool('size')
if pool>0
    matlabpool close
end
matlabpool open

% Define DGP
clear
clc

cor=0; % Short-run correlation of shocks
s2=1; %Variance of x2
Gamma_coint=[-0.15 0.15;0 0]; % Specification from Swenson (just that he uses 6 variables w/ 2 cointegrated)
[TEMP1,TEMP2]=eig([1 cor;cor s2]);
Sig=TEMP1*sqrt(TEMP2); % Spectral decomposition
Gamm=[0.2 0; 0 0.2]; % Short run VAR Coeff

TL=[150]% Number of Years
p=2; % Number of Variables 
B=10000; % Number of bootstrap replications
H=2000; %Number of monte-carlo experiments
nl=2;    % Number of lags in  the test
beta=1; % Beta-factor 
alpha=[0 0]; % intercept
endog=0.5;

rho=[1 0.85 1/3];
h=1;
discard=0;
for j=1:length(TL)
    T=TL(j);
    h=1;
    clear v1 v2 v3 v4
    p_hut1=ones(H,2,3);
    p_hut2=ones(H,2,3);
    p_hut3=ones(H,2,3);
    p_hut4=ones(H,2,3);
    p_hut5=ones(H,2,3);
    p_hut6=ones(H,2,3);
    p_hut7=ones(H,2,3);
    p_hut8=ones(H,2,3);
    p_hut9=ones(H,2,3);
    p_hut10=ones(H,2,3);
    pseudo_plmax=ones(H,2,3);
    pseudo_pcadf=ones(H,2,3);
    p1=ones(H,2,3);
    p2=ones(H,2,3);
    p3=ones(H,2,3);
    p4=ones(H,2,3);
    while h<=H
        %% Monte Carlo runs for size
        tic
        %     nlag=1;
        %     nlag_r=floor(4*(T/100)^.25);
        %
        nlag=2;
        nlag_r=2;
        u=normrnd(0,1,T+30,p);
        v1=cumsum(u(:,1));
        v2=cumsum(u(:,2));

        % Johansen DGP mit white noise error
        X=zeros(T+30,p);
        X(1,:)=u(1,:);
        dx=zeros(1,2);
        for t=2:T+30
            X(t,:)=X(t-1,:)+dx*Gamm+u(t,:)*Sig;
            
            dx=X(t,:)-X(t-1,:);
        end
        X(1:30,:)=[];
        
        ex=1;
        s=1;
        [Xstar, st_lmax, st_cad, p_lmax, p_cad, chi_star,p1(h,s,ex),p2(h,s,ex), p_hut1(h,s,ex),p_hut2(h,s,ex), ...
            p_hut3(h,s,ex), p_hut4(h,s,ex), p_hut5(h,s,ex), p_hut6(h,s,ex), p_hut7(h,s,ex),p_hut8(h,s,ex), p_hut9(h,s,ex),p_hut10(h,s,ex),...
            p3(h,s,ex),p4(h,s,ex), pseudo_plmax(h,s,ex), pseudo_pcadf(h,s,ex)]=replicate(B, X, nlag, nlag_r);

        X(1,:)=u(1,:);
        dx=zeros(1,2);
        for t=2:T+30
            X(t,:)=X(t-1,:)*(eye(2)+Gamma_coint)+dx*Gamm+u(t,:)*Sig;
            dx=X(t,:)-X(t-1,:);
        end
        X(1:30,:)=[];
        
        s=2;
        [Xstar, st_lmax, st_cad, p_lmax, p_cad, chi_star,p1(h,s,ex),p2(h,s,ex), p_hut1(h,s,ex),p_hut2(h,s,ex), ...
            p_hut3(h,s,ex), p_hut4(h,s,ex), p_hut5(h,s,ex), p_hut6(h,s,ex), p_hut7(h,s,ex),p_hut8(h,s,ex), p_hut9(h,s,ex),p_hut10(h,s,ex),...
            p3(h,s,ex),p4(h,s,ex), pseudo_plmax(h,s,ex), pseudo_pcadf(h,s,ex)]=replicate(B, X, nlag, nlag_r);
        %nlag=floor(4*(T/100)^.25);
        %     nlag=2;
        %     nlag_r=2;

        %% EG DGP
        %size
        X=[alpha(1)+beta*v1+v2, alpha(2)+v1+endog*v2];
        X(1:30,:)=[];
        ex=2;
        s=1;
        [Xstar, st_lmax, st_cad, p_lmax, p_cad, chi_star,p1(h,s,ex),p2(h,s,ex), p_hut1(h,s,ex),p_hut2(h,s,ex), ...
            p_hut3(h,s,ex), p_hut4(h,s,ex), p_hut5(h,s,ex), p_hut6(h,s,ex), p_hut7(h,s,ex),p_hut8(h,s,ex), p_hut9(h,s,ex),p_hut10(h,s,ex),...
            p3(h,s,ex),p4(h,s,ex), pseudo_plmax(h,s,ex), pseudo_pcadf(h,s,ex)]=replicate(B, X, nlag, nlag_r);                                                            
        % Power experiment
        v3=gen_ar1(T+30,rho(2),u(:,2));
        X=[alpha(1)+beta*v1+v3, alpha(2)+v1+endog*v3];
        X(1:30,:)=[];
       
        s=2;
        [Xstar, st_lmax, st_cad, p_lmax, p_cad, chi_star,p1(h,s,ex),p2(h,s,ex), p_hut1(h,s,ex),p_hut2(h,s,ex), ...
            p_hut3(h,s,ex), p_hut4(h,s,ex), p_hut5(h,s,ex), p_hut6(h,s,ex), p_hut7(h,s,ex),p_hut8(h,s,ex), p_hut9(h,s,ex),p_hut10(h,s,ex),...
            p3(h,s,ex),p4(h,s,ex), pseudo_plmax(h,s,ex), pseudo_pcadf(h,s,ex)]=replicate(B, X, nlag, nlag_r);


        %% Johansen DGP mit AR innovation
        %     nlag=floor(4*(T/100)^.25);
        %     nlag_r=floor(4*(T/100)^.25);
        %
        nlag=2;
        nlag_r=2;
        v4(:,1)=gen_ar1(T+30,rho(3),u(:,1)); % To be used as one alternative specification
        v4(:,2)=gen_ar1(T+30,rho(3),u(:,2));

        X(1,:)=u(1,:);
        dx=zeros(1,2);
        for t=2:T+30
            X(t,:)=X(t-1,:)*(eye(2)+Gamma_coint)+dx*Gamm+v4(t,:)*Sig;
            dx=X(t,:)-X(t-1,:);
        end
        X(1:30,:)=[];

        ex=3;
        s=2;
        [Xstar, st_lmax, st_cad, p_lmax, p_cad, chi_star,p1(h,s,ex),p2(h,s,ex), p_hut1(h,s,ex),p_hut2(h,s,ex), ...
            p_hut3(h,s,ex), p_hut4(h,s,ex), p_hut5(h,s,ex), p_hut6(h,s,ex), p_hut7(h,s,ex),p_hut8(h,s,ex), p_hut9(h,s,ex),p_hut10(h,s,ex),...
            p3(h,s,ex),p4(h,s,ex), pseudo_plmax(h,s,ex), pseudo_pcadf(h,s,ex)]=replicate(B, X, nlag, nlag_r);
        X(1,:)=u(1,:);
        dx=zeros(1,2);
        for t=2:T+30
            X(t,:)=X(t-1,:)+dx*Gamm+v4(t,:)*Sig;
            dx=X(t,:)-X(t-1,:);
        end
        X(1:30,:)=[];
        
        s=1;
        [Xstar, st_lmax, st_cad, p_lmax, p_cad, chi_star,p1(h,s,ex),p2(h,s,ex), p_hut1(h,s,ex),p_hut2(h,s,ex), ...
            p_hut3(h,s,ex), p_hut4(h,s,ex), p_hut5(h,s,ex), p_hut6(h,s,ex), p_hut7(h,s,ex),p_hut8(h,s,ex), p_hut9(h,s,ex),p_hut10(h,s,ex),...
            p3(h,s,ex),p4(h,s,ex), pseudo_plmax(h,s,ex), pseudo_pcadf(h,s,ex)]=replicate(B, X, nlag, nlag_r);


        time_new=toc;
        if h>1
            time=mean([time,time_new]);
        else
            time=time_new;
        end
        test=isnan(p_hut1(h,1,1)) | isnan(p_hut1(h,1,2))| isnan(p_hut1(h,1,3)) | isnan(p_hut1(h,2,1))|isnan(p_hut1(h,2,2)) | isnan(p_hut1(h,2,3));
        if test~=1
            h=h+1;
        else
            discard=discard+1;
        end

        crit=0.05;
        clc
        head={' Johansen ', ' Engle-Granger', ' Johansen* ', ' Engle-Granger*', 'additive', 'log', 'pesavento1 ', 'pesavento2', 'pesavento3', ' normal ', ' min ','max','Clay','Translog', 'naive'};
        disp('Percentage of Experiments done:')
        disp([num2str(round(h/H*1000)/10) '%'])
        disp('expected time until done:')
        disp( [num2str(round(time*(H-h)/6)/10) ' minutes'])
        disp(' ')
        disp(' ')
        disp('-------------------------------------------')
        disp('             RESULTS EG-Process')
        disp('-------------------------------------------')
        disp('-------------------------------------------')
        disp('Size Test')
        disp(['rejection rates at ' num2str(100*crit) '% lvl'])
        %disp('BH basis 2 | BH Pesavento |BH 4 tests |naive min | Hartung | Johansen | Engle-Granger| Johansen* | Engle-Granger*')
        s=1;
        ex=2;
        base=(sum([p1(:,s,ex),p2(:,s,ex), pseudo_plmax(:,s,ex), pseudo_pcadf(:,s,ex), p_hut1(:,s,ex), p_hut2(:,s,ex), ...
            p_hut3(:,s,ex), p_hut4(:,s,ex), p_hut5(:,s,ex), p_hut10(:,s,ex), p_hut6(:,s,ex),p_hut7(:,s,ex),p_hut8(:,s,ex), p_hut9(:,s,ex), min(p1(:,s,ex),p2(:,s,ex))]<crit)');
        table=head;
        for ttt=1:length(base)
            table{2,ttt}=base(ttt)/h;
        end
        disp(table)
        disp('-------------------------------------------')
        disp('Power Test')
        disp(['rejection rates at ' num2str(100*crit) '% lvl'])
        s=2;
        base=(sum([p1(:,s,ex),p2(:,s,ex), pseudo_plmax(:,s,ex), pseudo_pcadf(:,s,ex), p_hut1(:,s,ex), p_hut2(:,s,ex), ...
            p_hut3(:,s,ex), p_hut4(:,s,ex), p_hut5(:,s,ex),p_hut10(:,s,ex), p_hut6(:,s,ex),p_hut7(:,s,ex),p_hut8(:,s,ex), p_hut9(:,s,ex), min(p1(:,s,ex),p2(:,s,ex))]<crit)');
        table=head;
        for ttt=1:length(base)
            table{2,ttt}=base(ttt)/h;
        end
        disp(table)
        disp('-------------------------------------------')
        disp('             RESULTS Johansen-Process')
        disp('-------------------------------------------')
        disp('-------------------------------------------')
        disp('Size Test')
        disp(['rejection rates at ' num2str(100*crit) '% lvl'])
        s=1;
        ex=1;
        base=(sum([p1(:,s,ex),p2(:,s,ex), pseudo_plmax(:,s,ex), pseudo_pcadf(:,s,ex), p_hut1(:,s,ex), p_hut2(:,s,ex), ...
            p_hut3(:,s,ex), p_hut4(:,s,ex), p_hut5(:,s,ex), p_hut10(:,s,ex),p_hut6(:,s,ex),p_hut7(:,s,ex), p_hut8(:,s,ex), p_hut9(:,s,ex),min(p1(:,s,ex),p2(:,s,ex))]<crit)');
        table=head;
        for ttt=1:length(base)
            table{2,ttt}=base(ttt)/h;
        end
        disp(table)
        disp('-------------------------------------------')
        disp('Power Test')
        disp(['rejection rates at ' num2str(100*crit) '% lvl'])
        s=2;
        ex=1;
       base=(sum([p1(:,s,ex),p2(:,s,ex), pseudo_plmax(:,s,ex), pseudo_pcadf(:,s,ex), p_hut1(:,s,ex), p_hut2(:,s,ex), ...
            p_hut3(:,s,ex), p_hut4(:,s,ex), p_hut5(:,s,ex),p_hut10(:,s,ex), p_hut6(:,s,ex),p_hut7(:,s,ex), p_hut8(:,s,ex), p_hut9(:,s,ex),min(p1(:,s,ex),p2(:,s,ex))]<crit)');
        table=head;
        for ttt=1:length(base)
            table{2,ttt}=base(ttt)/h;
        end
        disp(table)
        disp('-------------------------------------------')
        disp('             RESULTS mixed-Process')
        disp('-------------------------------------------')
        disp('-------------------------------------------')
        disp('Size Test')
        disp(['rejection rates at ' num2str(100*crit) '% lvl'])
        s=1;
        ex=3;
        base=(sum([p1(:,s,ex),p2(:,s,ex), pseudo_plmax(:,s,ex), pseudo_pcadf(:,s,ex), p_hut1(:,s,ex), p_hut2(:,s,ex), ...
            p_hut3(:,s,ex), p_hut4(:,s,ex), p_hut5(:,s,ex),p_hut10(:,s,ex), p_hut6(:,s,ex),p_hut7(:,s,ex), p_hut8(:,s,ex), p_hut9(:,s,ex),min(p1(:,s,ex),p2(:,s,ex))]<crit)');
        table=head;
        for ttt=1:length(base)
            table{2,ttt}=base(ttt)/h;
        end
        disp(table)
        disp('-------------------------------------------')
        disp('Power Test')
        disp(['rejection rates at ' num2str(100*crit) '% lvl'])
        s=2;
        ex=3;
        base=(sum([p1(:,s,ex),p2(:,s,ex), pseudo_plmax(:,s,ex), pseudo_pcadf(:,s,ex), p_hut1(:,s,ex), p_hut2(:,s,ex), ...
            p_hut3(:,s,ex), p_hut4(:,s,ex), p_hut5(:,s,ex),p_hut10(:,s,ex), p_hut6(:,s,ex),p_hut7(:,s,ex),p_hut8(:,s,ex), p_hut9(:,s,ex), min(p1(:,s,ex),p2(:,s,ex))]<crit)');
        able=head;
        for ttt=1:length(base)
            table{2,ttt}=base(ttt)/h;
        end
        disp(table)
    end
    savename=['results_T_all_tests_varlag_2=' num2str(T)];
    save(savename)
end


%%
save results2
%% Graphical and regression analysis
return
close all
figure
plot(norminv(p1_size),norminv(p2_size),'.')
axis([-4 4 -4 4])
b=regress(min(max(norminv(p1_size)',-4),4),[min(max(norminv(p2_size)',-4),4) ones(H,1)])
corr(min(max(norminv(p1_size)',-4),4),min(max(norminv(p2_size)',-4),4))


hold on
plot(-4:4,b(1)*(-4:4)+b(2),'c')
%figure
plot(norminv(p1_power),norminv(p2_power),'r.')
axis([-4 4 -4 4])
b=regress(min(max(norminv(p2_power)',-4),4),[min(max(norminv(p1_power)',-4),4) ones(H,1)])
b=regress(min(max(norminv(p1_power)',-4),4),[min(max(norminv(p2_power)',-4),4) ones(H,1)])
hold on
plot(-4:4,b(1)*(-4:4)+b(2),'y')

figure
plot((p1_size),(p2_size),'.')
axis([0 1 0 1])
b=regress(min(max(norminv(p1_size)',-4),4),[min(max(norminv(p2_size)',-4),4) ones(H,1)])
corr(min(max(norminv(p1_size)',-4),4),min(max(norminv(p2_size)',-4),4))


figure
test=min(max(norminv(p1_size)',-4),4)-.55*min(max(norminv(p2_size)',-4),4);
hist(test,30)
figure
test=min(max(norminv(p1_power)',-4),4)-.55*min(max(norminv(p2_power)',-4),4);
hist(test,30)

figure
plot(norminv(p21_size),norminv(p22_size),'.')
axis([-4 4 -4 4])
b=regress(min(max(norminv(p21_size)',-4),4),[min(max(norminv(p22_size)',-4),4) ones(H,1)])
corr(min(max(norminv(p21_size)',-4),4),min(max(norminv(p22_size)',-4),4))
hold on
plot(-4:4,b(1)*(-4:4)+b(2),'c')
%figure
plot(norminv(p21_power),norminv(p22_power),'r.')
axis([-4 4 -4 4])
b=regress(min(max(norminv(p21_power)',-4),4),[min(max(norminv(p22_power)',-4),4) ones(H,1)])
hold on
plot(-4:4,b(1)*(-4:4)+b(2),'y')

figure
test=min(max(norminv(p21_size)',-4),4)-.55*min(max(norminv(p22_size)',-4),4);
hist(test,30)

figure
test=min(max(norminv(p21_power)',-4),4)-.55*min(max(norminv(p22_power)',-4),4);
hist(test,30)

figure
plot(norminv(p31_size),norminv(p32_size),'.')
axis([-4 4 -4 4])
b=regress(min(max(norminv(p31_size)',-4),4),[min(max(norminv(p32_size)',-4),4) ones(H,1)])
corr(min(max(norminv(p31_size)',-4),4),min(max(norminv(p32_size)',-4),4))
hold on
plot(-4:4,b(1)*(-4:4)+b(2),'c')



%figure
plot(norminv(p31_power),norminv(p32_power),'r.')
axis([-4 4 -4 4])
b=regress(min(max(norminv(p31_power)',-4),4),[min(max(norminv(p32_power)',-4),4) ones(H,1)])
hold on
plot(-4:4,b(1)*(-4:4)+b(2),'y')