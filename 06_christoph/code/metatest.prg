new; /* a joint test for cointegration. some Monte Carlo */
cls;
library fisher coint panel pgraph;

dat_start = date;
mth_st = ftos(dat_start[2],"%*.*lf",1,0);
day_st = ftos(dat_start[3],"%*.*lf",1,0);
hr_st = trunc(dat_start[4]/(100*60*60));
min_st = dat_start[4]/(100*60) - hr_st*60;
min_st = ftos(min_st,"%*.*lf",1,0);
hr_st = ftos(hr_st,"%*.*lf",1,0);

dir = ChangeDir("c:/dokumente und einstellungen/christoph und steffi/eigene dateien/promotion/research/joint test for cointegration/results");
/*dir = ChangeDir("d:/benutzer/hanck/research/panel cointegration/results/cross sec dep");*/
{theta,thetaspec} = respsurfcoef(3 | 4 | 5 | 221 | 5);
begloop = time;

VAR_bootstrap = 0; /* set to 1 if you want a VAR bootstrap a la Chang et al (2006) JECM, too */

/* DGP: */
K = 2;
N = 1; /* # cross sections */
T = 200; /* # time series observations */
discarded = 30; /* to be discarded later to remove initial conditions effect */
a1 = 0;
a2 = -1;

psi = 0.5;
sigma = 1;
V = 1~psi*sigma|psi*sigma~sigma^2;  
L = chol(V);

rho = 1|.8;
gamma_ma = 0; /* if non-zero we have MA errors */

rep = 500; 
smpl_ini = 30; /* # of obs. to be discarded in the w_hat_star series to get swung-in series. 
                   this number must not exceed the minimum (effective!) time series length!  */
B = 500; 
inc = 2;
noq = 2; /* # of q's you try */

for gam (1,rows(gamma_ma),1);
    size_mat = arrayinit(rows(T)|rows(N)|rows(rho)|noq|rows(a1)|2|3,4);
    rej_mat_bootstrap = arrayinit(rows(T)|rows(N)|rows(rho)|noq|rows(a1)|rep|2,4);
   
    for rh (1,rows(rho),1);
    for tt (1,rows(T),1);
    for endog (1,rows(a1),1);
        /*q = seqa(2,inc,trunc(sqrt(T[tt]))/inc);*/ /* start w 1 not yet possible */
        q = 2|trunc(4*(T[tt]/100)^(.25));
        if rows(q) /= noq;
        print "Inappropriate noq / rows(q) pair!";
        stop;
        endif;
        for qq (1,rows(q),1);
        for nn (1,rows(N),1);
            alfa = 5*rndu(N[nn],1);
            beta = 1+rndu(1,N[nn]);
            for jj (1,rep,1);
                state = -1;
                {R,state} = rndKMn((T[tt]+discarded)*N[nn],K,state);
                R2 = R*L;
                u1 = (reshape(R2[.,1],N[nn],T[tt]+discarded))';
            
                u2 = zeros(1,N[nn])|(reshape(R2[.,2],N[nn],T[tt]+discarded))';
                u2 = u2[2:rows(u2),.] + gamma_ma[gam]*u2[1:(rows(u2)-1),.];
               
                z = recserar(u1,zeros(1,N[nn]),rho[rh]*ones(1,N[nn])); /* with rho<1 you get a power study */
                w = recserar(u2,zeros(1,N[nn]),ones(1,N[nn]));  /* w has a unit root. */
              
                x = ((a2-a1[endog]*beta)^-1).*(a1[endog]*reshape(alfa,rows(z),N[nn])+(a1[endog]*z-w)); /* figure out whether constant is appropriate */
                y = ((a2-a1[endog].*beta)^-1).*(a2*reshape(alfa,rows(z),N[nn])+(a2*z-beta.*w));
                
                /* now, discard */
                x = x[discarded+1:rows(x),.];
                y = y[discarded+1:rows(y),.];
                
                /* ---------- AR sieve bootstrap ------------- */
                all_eps_res = zeros(rows(x)-1,N[nn]);
                big_A_mat = zeros(N[nn],q[qq]);
            
                q_opt = zeros(N[nn],1);

                {alph,stat_aeg,crits,res,b_ols} = cadfw_lagch(y[.,1],x[.,1],0);
                {tr,lam} = trace_test(y[.,1]~x[.,1],-1,1); /* Gutierrez */
                stat_joh = tr[1];      
                /* ADF and Johansen version of the combination tests */
                {chisquarestat, margvalsc,invn,margvalsi,logitteststar,margvalsl,pvalues_adf} = fisher(&mackinn,stat_aeg|stat_joh,2|2,1|4,2|1,T[tt]|0,theta,thetaspec,0|0);
                 
                /*{q_opt,aic} = var_aic(diff(res,1)~diff(x[.,nn],1));
                 q_opt = var_lags(diff(res,1)~diff(x[.,nn],1),trunc(T^(1/3))); SBC by Gutierrez 
                 {q_opt[nn,.],lr_stats} = var_lag_lr(diff(res,1)~diff(x[.,nn],1));*/ 
              
                q_opt = q[qq];
            
                /* now with an AR instead: */
                /*{k,bic,pic} = arord(diff(res,1),trunc(12*(T/100)^(1/5)),-1);
                k=minindc(k[.,1]);
                k=12; */ 

                {phihat,sighat} = yule(diff(res,1),q_opt[1,.]);
                /* residual series */
                phi_times_dyk = zeros(rows(diff(res,1))-q_opt[1,.],q_opt[1,.]);
                all_eps_res[q_opt[1,.]+1:rows(all_eps_res),1] = trimr(diff(res,1),q_opt[1,.],0);
                for j (1,q_opt[1,.],1);
                    phi_times_dyk[.,j] = phihat[j]*trimr(diff(res,1),q_opt[1,.]-j,j);
                    all_eps_res[q_opt[1,.]+1:rows(all_eps_res),1] = all_eps_res[q_opt[1,.]+1:rows(all_eps_res),1] - phi_times_dyk[.,j];
                endfor;
                                           
                big_A_mat[1,1:rows(phihat)] = phihat';
                /* throw away the number of observations that are zero at the beginning */ 
                throw_away = maxc(q_opt); 
                all_eps_res_cut = all_eps_res[throw_away+1:rows(all_eps_res),.];
                
                /* centered residuals*/
                all_eps_res_cen = all_eps_res_cut - reshape(meanc(all_eps_res_cut)',rows(all_eps_res_cut),N[nn]);
                
                /* resample---bootstrap loop */
                aeg_star = zeros(B,1);
                joh_star = zeros(B,1);
                                
                for bb (1,B,1);
                    /* AR resamples */
                    all_eps_res_star = resample_mat(all_eps_res_cen,smpl_ini); 
                    
                    /* AR sieve bootstrap dataset */
                    what_star = recserar(all_eps_res_star,zeros(q[qq],N[nn]),big_A_mat');
                    ustar = cumsumc(what_star[smpl_ini-q[qq]:rows(all_eps_res_star),.]); 
                    
                    /* retrieving the constant---like this only for k=2? */
                    a_ols = reshape(meanc(y)'-(b_ols'.*meanc(x))',T[tt],N[nn]);    
        
                    ystar = a_ols + b_ols.*x + ustar;
                                      
                    {alph,aeg_star[bb,1],crits,res_star,b_star} = cadfw_lagch(ystar[.,1],x[.,1],0);
                    {tr,lam} = trace_test(ystar[.,1]~x[.,1],-1,1); /* Gutierrez */
                    joh_star[bb,1] = tr[1];    
                                      
                endfor; /* end B loop */
                p_aeg_star = rankindx(aeg_star,1)/B;
                p_joh_star = (N+1-rankindx(joh_star,1))/B;
                
                fisherstat = -2*ln(p_aeg_star.*p_joh_star);

                bootstrap_pval = sumc(fisherstat .> chisquarestat*ones(B,1))/B;

                /* record rejections here */
                rej_mat_bootstrap[tt,nn,rh,qq,endog,jj,1] = (bootstrap_pval lt 0.05);
                                                        
            endfor; /* end rep loop */
            format 8,4;
            print "-------------------------------------------------------------------";
            print "Case: rho=" rho[rh] " gamma=" gamma_ma[gam] " q=" q[qq] " T=" T[tt] "  N=" N[nn] " a1=" a1[endog]; 
            print "-------------------------------------------------------------------";
            size_mat[tt,nn,rh,qq,endog,1,1] = meanc(rej_mat_bootstrap[tt,nn,rh,qq,endog,.,1]);
            print "P_chi^2 bootstrap rejection rate: " size_mat[tt,nn,rh,qq,endog,1,1];
            print "";

            /* save output into file */
            /* TBD: Write this directly into EXCEL with getmatrix and xlswritem */
            dat = date;
            
            mth = ftos(dat[2],"%*.*lf",1,0);
            day = ftos(dat[3],"%*.*lf",1,0);
            hr = trunc(dat[4]/(100*60*60));
            min = dat[4]/(100*60) - hr*60;
            min = ftos(min,"%*.*lf",1,0);
            hr = ftos(hr,"%*.*lf",1,0);
            curr_T = ftos(T[tt],"%*.*lf",1,0);
            curr_N = ftos(N[nn],"%*.*lf",1,0);
            curr_q = ftos(q[qq],"%*.*lf",1,0);
            curr_gam = ftos(10*gamma_ma[gam],"%*.*lf",1,0);
            curr_rho = ftos(100*rho[rh],"%*.*lf",1,0);
            curr_a1 = ftos(a1[endog],"%*.*lf",1,0);
            filename = "st" $+ day_st $+ "_" $+ mth_st $+ "_" $+ hr_st $+ "h" $+ min_st $+ "cpT" $+ curr_T $+ "cpN" $+ curr_N $+ "rho" $+ curr_rho $+ "q" $+ curr_q $+ "gam" $+ curr_gam $+ "a1_" $+ curr_a1 $+ "d" $+ day $+ "_" $+ mth $+ "_" $+ hr $+ "h" $+ min $+ ".asc";
            
            format 8,3;
            output file = ^filename reset;
            screen off;
            print "replications: " rep;
            print "bootstraps: " B;        
            print " ";
            print "a1:" a1';
            print "mean beta:" meanc(beta');
            
            print "psi:" psi;
            print "serial correlation: " gamma_ma';
            print "rho:" rho';
            print "relevant q for current T:" q';
            print "T:" T';
            print "N:" N'; 
            format 8,4;
            print "------------------------------------------------------------";
            print "bootstrap, asymptotic and VAR bootstrap rejection rates:";
            print "------------------------------------------------------------";
            print "read: T,N,rho,fitted lag length,endogeneity,cross-sectional dependency";
            print size_mat;
            output off;
            screen on;
        
        endfor; /* end N loop */   
        endfor; /* end q loop */   
    endfor; /* end a1 loop */
    endfor; /* end T loop */
    endfor; /* end rho loop */
endfor; /* end gam_ma loop */
endloop = time;
print "time required for the loop: " endloop-begloop; 
/* incidentally, GAUSS does not want a blank in between the times, whysoever */

end;

/* resamples a matrix of length N+smpl_inilength (for discarding initials) from an Nx1 vector x*/
/* this is the required version for preserving cross-sectional dependency */
proc resample_mat(x,smpl_inilength);
local n;
    n = rows(x);
retp(x[floor(n*rndu(n,1)+1),.]|x[floor(smpl_inilength*rndu(smpl_inilength,1)+1),.]);
endp;

/* MODIFIED: Certain checks switched off. */ 
proc (5) = cadfw_lagch(y,x,p) ; 
     local beta,dep,k,z,res,so,var_cov,xx,r;
     local timep,t,m,xmat,nobs,laglgth,b ;
     /*if (p < -1);
        "Error: p cannot be < -1";
         retp(0,0,zeros(6,1));
     endif ;*/
     nobs    = rows(x);
     /*if (nobs - (2*l) + 1 < 1) ;
        "Error: l is too large; negative degrees of freedom.";
         retp(0,0,zeros(6,1));
     endif ;*/
     y       = detrend(y,p);
     x       = detrend(x,p);
     b       = y/x;   
     r       = y - x*b;
     {laglgth} = lag_choice(r,0,1);
     dep     = diff(r,1);
     k       = 0     ;
     z       = trimr(lagn(r,1),1,0) ;
     if (laglgth > 0) ;
        do until k >= laglgth ;
           k = k + 1 ;
       z = z~lagn(dep,k) ;
        endo ;
     endif;
     z       = trimr(z,laglgth,0) ;
     dep     = trimr(dep,laglgth,0) ;
     beta    = detrend(dep,0)/detrend(z,0) ;
     res     = dep - z*beta ;
     so      = (res'res)/(rows(dep)-cols(z));
     var_cov = so*inv(z'z) ;
     retp(beta[1,1]+1.0,beta[1,1]/sqrt(var_cov[1,1]),rztcrit(nobs,cols(x),p),r,b) ;
endp ;

proc (2) = trace_test(y,co,lags);
    local dy, z, p, r0t, rkt, sig,a,c,b,d,e,lr1, zz, i, obs,mat ;
    local s00, sk0, skk, ly ;

    dy = trimr(y,1,0)-trimr(lag(y),1,0) ;
    z = lag(dy) ;
    p = 2;
    do while p lt lags;
        z = z~lagn(dy,p);
        p = p + 1 ;
    endo ;

    if co eq -1 ;
        z = trimr(z,lags,0);
        dy = trimr(dy,lags,0);
        r0t = dy - z*(dy/z);
        ly = trimr(lagn(y,lags),lags+1,0);
        rkt = ly - z*(ly/z);

    elseif co eq 0 ;

        z = trimr(z,lags,0);
        dy = trimr(dy,lags,0);
        z = z-meanc(z)';
        dy = dy-meanc(dy)';
        r0t = dy - z*(dy/z);
        ly = trimr(lagn(y,lags),lags+1,0);
        ly = ly -meanc(ly)';
        rkt = ly - z*(ly/z);

    elseif co eq 1 ;

        z = trimr(z,lags,0);
        dy = trimr(dy,lags,0);
        obs = rows(z);
        mat = ones(obs,1)~seqa(0,1,obs)./obs;
        z = z - mat*(z/mat);
        dy = dy - mat*(dy/mat);
        r0t = dy - z*(dy/z);
        ly = trimr(lagn(y,lags),lags+1,0);
        ly = ly - mat*(ly/mat);
        rkt = ly - z*(ly/z);
    endif ;

    skk = rkt'rkt/rows(rkt) ;
    sk0 = rkt'r0t/rows(rkt) ;
    s00 = r0t'r0t/rows(r0t) ;
    sig = sk0*inv(s00)*(sk0');

    { a,b,d,c } = eigrg2(inv(skk)*sig);

    d = d*inv(chol(d'skk*d)) ;

    zz = (-a)~(d') ;
    zz = sortc(zz,1) ;
    a = -zz[.,1] ;
    d = zz[.,2:cols(zz)]' ;

    lr1 = zeros(rows(a),1);
    i = 1 ;
    do while (i <= rows(a));
        lr1[i,1] = -rows(rkt)*sumc(trimr(ln(1-a),i-1,0));
        i = i + 1 ;
    endo ;

    retp(lr1,a);
endp ;

/* yule-walker of an Nx1 vector x of order k */
proc (2) = yule(x,k);
local xbar,gamm,j,gamma_k,phihat,sighat;
    xbar = meanc(x);
    gamm = zeros(k+1,1);
    for j (0,k,1);
        gamm[j+1] = ((trimr(x,j,0)-reshape(xbar,rows(x)-j,1))'(trimr(x,0,j)-reshape(xbar,rows(x)-j,1)))/(rows(x)-j);
    endfor;
    if k>1;
        gamma_k = toeplitz(gamm[1:k]);
    else;
        gamma_k = gamm[1];
    endif;
    /* coef vector acc. Schlittgen/Streitberg p.254 */
    phihat = inv(gamma_k)*gamm[2:k+1];
    /* ATTENTION:  if it turns out that non-invertible gamma_k arise also with inv instead
                   of invpd, replace this procedure with the one from Schlittgen (check arimafit.src, 
                   delete all the aic and variance-covariance stuff, you just need "alpha" */
    /* residual variance acc. Schlittgen/Streitberg p.257, above eq.6.1.1.3 solved
       for sig^2 and with consistent estimates  */
    sighat = gamm[1] - phihat'gamm[2:k+1];
retp(phihat,sighat);
endp;