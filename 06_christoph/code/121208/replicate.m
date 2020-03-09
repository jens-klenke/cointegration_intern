function [Xstar, st_lmax, st_cad, p_lmax, p_cad, chi_star2,p1,p2, p_hut1,p_hut2, p_hut3, p_hut4, p_hut5, p_hut6, p_hut7, p_hut8, p_hut9, p_hut10,p3,p4, pseudo_plmax, pseudo_pcadf ]=replicate(B, X, nlag, nlag_r)
trend=1;
r=0;
nvar=size(X,2);
T=size(X,1)-nlag-1;
% Swenson Algorithm to Bootstrap an ECM
% Step (i)
results=ecm(X, nlag, nvar); 
% CH 11.12.08:
% replace with results=ecm(X, nlag, 0); (VAR in 1st diffs) or result = vare(X,nlag) (VAR in levels)
% (wir haben also bisher r=p CI-Beziehungen erzwungen. also eher das
% gegenteil von dem was wir dem hier zufolge machen sollten...vermutlich
% weil er "estimate (1) with r=p" schreibt, was ich verwirrend finde.
% hoffen wir, dass es jetzt keine probleme macht.)
% I am not sure whether the Gamma matrices still work out correctly 
% - does vare output like ecm?
% - do we need to do something about lag length if we estimate in levels
% (also: notation changes in the paper)?
eps=zeros(T,nvar);
for j=1:nvar
    eps(:,j)=results(j).resid;
    Coeff(:,j)=results(j).beta;
end
Coint_M=Coeff(end-nvar:end-1,:);
Coeff(end-nvar:end-1,:)=[];
% CH 11.12.08: d.h. in Coeff sind die Gamma-Matrizen aber nicht die
% Pi-Matrix drin? dann w�rde es unten in polymroot mE passen

% CB 17.03.08: Would it make sense to estimate a VAR in first differences
% instead?
% CH 11.12.08: ja, das sieht so aus---siehe unsere fn 2!

% Columns are the different variables, lines the time
%% Check for explosive root

for j=1:nlag
    G(1:nvar,1+(j-1)*nvar:j*nvar)=Coeff(nlag-j+1+(0:nvar-1).*nlag,:)';
end

Temp=polymroot([G'; eye(nvar)]);
test=min(abs(Temp));

if test>1
% CH 11.12.08: da wir hier nur nach roots streng gr��er 1 gucken spricht das
% sehr daf�r, dass wir das schon so wie in meiner neuen fassung des
% algorithmus beschrieben machen, da sonst ja immer eine root=1 dabei w�re
% und die bedingung oben nie erf�llt w�re.

% Step (ii)
%results=ecm(X, nlag, r);
%% Johansen Test
jres = johansen(X,trend,nlag); %Johansen Test
ecvectors = jres.evec; % recover error correction vectors
lmax=jres.lr2(1+r);
ltr=jres.lr1(1+r);

crit_valm=jres.cvm;
test_level=sum(lmax>crit_valm(1,:))+1;
temp=[0.5 0.09 0.04 0.009]; 
pseudo_plmax=temp(test_level);

% To be done:
% only generate bootstrap replication if (iii) in Algorithm 1 / 2 of Swenson
% holds
% Think about trend specification

%% Engle-Granger Test
results = cadf(X(:,1),X(:,2:end),trend,nlag_r,1); % Need to check 'trend' % last parameter is 1 if crit vals are calculated
%ols_res=ols(detrend(X(:,1),trend),detrend(X(:,2:end),trend));
R2_t=results.R2; 
cad=results.adf;

crit_cadf=results.crit;
test_level=sum(cad<crit_cadf)+1;
temp=[0.999 .98 .94 .89 0.09 0.04 0.009]; 
pseudo_pcadf=temp(test_level);
%% Phillips Perron Test
results = phillips(X(:,1),X(:,2:end),trend,nlag_r);
php=results.pstat;

%% Generate Bootstrap replication
Xstar=zeros(size(X,1),size(X,2),B);
st_lmax=zeros(B,1);
st_cad=zeros(B,1);
st_ltr=zeros(B,1);
st_php=zeros(B,1);
R2_boot=zeros(B,1);
%p_lmax=zeros(B,1);
%p_ltr=zeros(B,1);
%p_cad=zeros(B,1);
%p_php=zeros(B,1);
ZZ=X(1:nlag+1,:);
warning('off','stats:regress:NoConst')
parfor b=1:B
    index=ceil(rand(length(eps),1)*length(eps));
    s_err=eps(index,:);
    Y=DGP(ZZ,Coeff, nlag, s_err, ecvectors, r);
    Xstar(:,:,b)=Y;
    result = johansen(Y,trend,nlag);
    st_lmax(b)=result.lr2(1+r); % Do Johansen test on generated data
    st_ltr(b)=result.lr1(1+r);
    results = cadf(Y(:,1),Y(:,2:end),trend,nlag,0); %Do EG-Test on generated data
    st_cad(b)=results.adf;
    R2_boot(b)=results.R2;
    results2 = phillips(Y(:,1),Y(:,2:end),trend,nlag);
    st_php(b)=results2.pstat;
end
warning('off','stats:regress:NoConst')
R2=mean(R2_boot);
% for b=1:B
%     p_lmax(b)=mean(st_lmax>=st_lmax(b)); % give each boostrap replication an orderinmg number on the lmax stat
%     p_cad(b)=mean(st_cad<=st_cad(b));
% end
ST=[st_lmax st_ltr st_cad st_php R2_boot];
ST=sortrows(ST,1);
ST=[ST ((1:-1/B:1/B))'];
ST=sortrows(ST,2);
ST=[ST ((1:-1/B:1/B))'];
ST=sortrows(ST,3);
ST=[ST ((1/B:1/B:1))'];
ST=sortrows(ST,4);
ST=[ST ((1/B:1/B:1))'];
p_lmax=ST(:,6);
p_ltr=ST(:,7);
p_php=ST(:,9);
p_cad=ST(:,8);
R2_boot=ST(:,5);
%corr([p_lmax p_ltr p_php p_cad])


p1=mean(st_lmax>=lmax)+.0000000000000000000000000000001;
p2=mean(st_cad<=cad)+.00000000000000000000000000001;
p4=mean(st_php<=php)+.00000000000000000000000000001;
p3=mean(st_ltr>=ltr)+.00000000000000000000000000001;


%% Direct Combination of p-values
chi_star1=-2*((p_lmax)+(p_cad));
stat_hut1=-2*sum(([p1 p2]));
p_hut1=mean(stat_hut1<=chi_star1);


%% Log Combination P-Values
% Define t-stats on the basis of probits
%T_ST=norminv(max(min(ST(:,5:8),.9999),0.0001));
chi_star2=-2*(log(p_lmax)+log(p_cad));
stat_hut2=-2*sum(log([p1 p2]));
p_hut2=mean(stat_hut2<=chi_star2);
%D=inv(C);
%statistic=zeros(1,B);
% for b=1:B
%     statistic(b)=-2*log(prod(normcdf([t_lmax(b),t_cad(b)]*D))); %Correlation corrected stats
% end
%statistic=-2*sum(log(normcdf([t_lmax,t_cad]*D)),2); %Correlation corrected stats


%% Approach based on pesavento

chi_star3=-2*((1-R2)*log(p_lmax)+R2*log(p_cad));
stat_hut3=-2*((1-R2)*log(p1)+R2*log(p2));
p_hut3=mean(stat_hut3<=chi_star3);

%% Alternative 
chi_star4=-2*((R2)*log(p_lmax)+(1-R2)*log(p_cad));
stat_hut4=-2*((R2)*log(p1)+(1-R2)*log(p2));
p_hut4=mean(stat_hut4<=chi_star4);



%% Normal transformation

%t_stat=(t_lmax+t_cad)/sqrt(2+2*rho(2,1));
stat_hut5=-(norminv(p1)+norminv(p2));
%p_hart_hut=mean(t_hart_hut>t_stat);
chi_star5=-(norminv(p_lmax)+norminv(p_cad));
p_hut5=mean(stat_hut5<=chi_star5);
%% Minimum
chi_star6=min(p_lmax,p_cad);
stat_hut6=min(p1,p2);
p_hut6=mean(stat_hut6>=chi_star6);

%% Maximum
chi_star7=max(p_lmax,p_cad);
stat_hut7=max(p1,p2);
p_hut7=mean(stat_hut7>=chi_star7);

%% Clay Copula (lambda=-1)
chi_star8=(p_lmax.^-1+p_cad.^-1-1).^-1;
stat_hut8=(p1.^-1+p2.^-1-1).^-1;
p_hut8=mean(stat_hut8>=chi_star8);
%% Translog
chi_star9=-2*(log(p_lmax)+log(p_cad)+.3*(log(p_lmax).*log(p_cad)));
stat_hut9=-2*(log(p1)+log(p2)+.3*(log(p1).*log(p2)));
p_hut9=mean(stat_hut9<=chi_star9);
%% Translog
chi_star10=-2*((1-R2_boot).*log(p_lmax)+R2_boot.*log(p_cad));
stat_hut10=-2*((1-R2_t)*log(p1)+R2_t*log(p2));
p_hut10=mean(stat_hut10<=chi_star10);
else
    Xstar=NaN;
    st_lmax=NaN;
    st_cad=NaN;
    p_lmax=NaN;
    p_cad=NaN;
    chi_star2=NaN;
    p1=NaN;
    p2=NaN;
    p_hut1=NaN;
    p_hut2=NaN;
    p_hut3=NaN;
    p_hut4=NaN;
    p_hut5=NaN;
    p_hut6=NaN;
    p_hut7=NaN;
    p_hut8=NaN;
    p_hut9=NaN;
    p_hut10=NaN;
    p3=NaN;
    p4=NaN;
    pseudo_plmax=NaN;
    pseudo_pcadf=NaN;
end
    
    
    