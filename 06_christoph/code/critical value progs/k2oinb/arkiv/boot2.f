      PROGRAM boot 
C Program for bootstrapping simulated VAR modell.
C Model with linear term with restrictions. Use of OLS
C estimates for short run parameters. Computes 
C estimates of rank based on
C trace statistics
C Estimates rank by bootstrapping
C
C Calls the following routines:
C dpgsim, errsim, dgpini, rempmat, rrest2, rrlsq1, 
C berrsim, outcnt.
C 
C In addition the following routines from
C Press, W.H., S.A. Teukolski, W.T. Vetterling and 
C B.P. Flannery (1992): Numerical recipes in Fortran (2nd edn.)
C Cambridge University Press are needed
C ludcmp, lubksb, choldc, gasdev, ran2, dchldc, eigsrt,
C jacobi, svbksb svdcmp, pythag, balanc, elmhes, hpq
      INTEGER NNOBS,PD,KP,MD      
      PARAMETER(NNOBS=300,PD=7,KP=3,MD=3)
      INTEGER nobs,idum,i,j,l,l1,k,d,m,p,r,nsim,nboot,
     +ilamax,itrce,i1trce,i1lamax,i2lamax,i2trce,i3lamax,i3trce,
     +dim,eri,esala,esatr,rboot,esbla,esbtr,tab1t(PD,PD),
     +cnt,esim
      REAL endog(PD,NNOBS),exog(MD,NNOBS),ini(PD,KP),eerr(PD,NNOBS),
     +s00(PD,PD),s01(PD,PD),s10(PD,PD),s11(PD,PD),
     +beta(PD,PD),lambda(PD),gaen(PD,PD*KP),
     +alpha(PD,PD),mu(PD),gaex(PD,MD),res(PD,NNOBS),bootres(PD,NNOBS),
     +sbeta(PD,PD),slambda(PD),sgaen(PD,PD*KP),
     +salpha(PD,PD),smu(PD),sgaex(PD,MD),
     +simob(PD,NNOBS),pi(PD,PD),bs00(PD,PD),bs01(PD,PD),
     +bs10(PD,PD),bs11(PD,PD),bbeta(PD,PD),blam(PD),
     +balph(PD,PD),blamax,btrce,lamax(PD),trce(PD),
     +nomlev(3),pr,const,coverr(PD,PD),cov2(PD,PD),
     +crtrce(PD),crlamx(PD)
      nobs=100
      nsim=100
      esim=nsim
      nboot=500
      m=1    
      do i=1,nobs
        exog(1,i)=1.0
      enddo

      idum=-3
      nomlev(1)=0.025
      nomlev(2)=0.05
      nomlev(3)=0.10
      ilamax=0
      itrce=0
      i1trce=0
      i1lamax=0
      i2lamax=0
      i2trce=0
      i3lamax=0
      i3trce=0
C Parameters of DGP VAR 
      p=5
      k=1
      r=1
      dim=r+p*k
      do i=1,p
        do j=1,r
          beta(i,j)=0.0
          alpha(i,j)=0.0
        enddo
        do j=1,p
          endog(i,j)=0.0
        enddo
        mu(i)=0.0
      enddo
      do i=1,r
        beta(p+1,i)=0.0      
      enddo
      alpha(1,1)=-0.4
      alpha(2,1)=-0.4
      beta(1,1)=1.0
      do i=1,p
        do j=1,p
          coverr(i,j)=0.0      
        enddo
      coverr(i,i)=1.0
      enddo
C Type of residuals, resind=1: OLS,resind=2: restricted ?
      resind=1
      do i=1,p+1
        do j=1,p+1
          tab1t(i,j)=0
        enddo
      enddo 
C Type of residuals, resind=1: OLS,resind=2: restricted ?
      resind=1
C Any mixture of errors?
C er1=1: Standard normal, eri=2:Mixture of normals
      eri=1
      pr=0.9
      const=10.0 
      do i=1,p
        do j=1,p
          cov2(i,j)=(pr+const*const*(1-pr))*coverr(i,j)
        enddo
      enddo
      if (eri.EQ.2) then
        do i=1,p
          do j=1,p         
            coverr(i,j)=cov2(i,j)
          enddo
        enddo
      endif

C Simulates NOSIM realizations of DGP
      esatr=0
      esala=0
      esbtr=0
      esbla=0
      do i=1,nsim
        call dgpini(k,KP,p,PD,r,dim,
     +beta,alpha,gaen,mu,ini,idum,cov2)
        do j=1,p
        enddo
        call errsim(nobs,NNOBS,k,KP,p,PD,eri,pr,const,idum,eerr,coverr) 
        call dgpsim(nobs,NNOBS,m,MD,k,KP,p,PD,r,
     +beta,alpha,gaen,mu,exog,mu,ini,eerr,idum,endog)
C Computes estimates of parameters in simulated realizations 
        call rempmat2(endog,exog,nobs,NNOBS,k,KP,p,PD,m,MD,s00,s01,
     +s10,s11)
        call rrest2(s00,s01,s11,p,PD,lambda,sbeta,salpha)
        if (resind.EQ.1) then
          call rrlsq1(endog,exog,nobs,NNOBS,k,KP,p,PD,m,mD,p,sbeta,
     +    sgaen,salpha,sgaex,res)
        elseif (resind.EQ.2) then
          call rrlsq1(endog,exog,nobs,NNOBS,k,KP,p,PD,m,mD,r,sbeta,
     +    sgaen,salpha,sgaex,res)
        else
           continue
        endif
        do j=1,p 
          cnt=0
          call outcnt(k,KP,p,PD,r,
     +    sbeta,salpha,sgaen,cnt)
          if (cnt.GT.0) then 
            write(*,*)'cnt',j,cnt
            esim=esim-1
            goto 10
          endif 
        enddo
        do l=1,p
          lamax(l)=-(nobs-k)*log(1-lambda(l))
        enddo
        trce(p)=lamax(p)
        do l=1,p-1        
          trce(p-l)= trce(p-l+1)-(nobs-k)*log(1-lambda(p-l))
        enddo
        crlamx(1)=12.25
        crlamx(2)=18.96
        crlamx(3)=25.54
        crlamx(4)=31.46
        crlamx(5)=37.52
        crtrce(1)=12.25
        crtrce(2)=25.32
        crtrce(3)=42.44
        crtrce(4)=62.99
        crtrce(5)=87.31
        do j=1,p
        if (trce(j).LE.crtrce(p+1-j)) then
        esatr=j-1
        goto 11
        endif
        enddo
        esatr=p
11      continue
C     ******  Integer r determines rank of 
C     impactmatrix used for constructing pseudoobservations.
        do rboot=1,p
          itrce=0
          ilamax=0 
          do j=1,nboot
            call berrsim(res,nobs,NNOBS,k,KP,p,PD,idum,bootres)
            call dgpsim(nobs,NNOBS,m,MD,k,KP,p,PD,rboot-1,
     +sbeta,salpha,sgaen,sgaex,exog,sgaex,ini,bootres,idum,simob)
            call rempmat2(simob,exog,nobs,NNOBS,k,KP,p,PD,m,MD,     
     +bs00,bs01,bs10,bs11)
            call rrest2(bs00,bs01,bs11,p,PD,blam,bbeta,balph)
            btrce=0.0
            do l=rboot,p
              btrce= btrce-(nobs-k)*log(1.0-blam(l))
            enddo
            if(btrce.GT.trce(rboot))itrce=itrce+1 
            if(-(nobs-k)*log(1.0-blam(rboot)).GT.lamax(rboot))
     +      ilamax=ilamax+1 
          enddo
C end of bootstrap loop
          if(REAL(itrce)/REAL(nboot).GT.0.05)then
            esbtr=rboot-1
            goto 12
          endif
          esbtr=p        
        enddo
C end of "rank estimation" loop
12     continue 
        tab1t(esatr+1,esbtr+1)=tab1t(esatr+1,esbtr+1)+1
10    continue 
      enddo
C end of simulation loop

      write(*,*)"Number of simulations",nsim
      write(*,*)"Number of effective simulations",esim      
      write(*,*)"Number of bootstraps",nboot
      write(*,*)"Number of observations",nobs
      write(*,*)"Dimension",p
      write(*,*)"Bootstapped rank",r
      write(*,*)"laglength",k
      write(*,*)
      write(*,*)"beta"
      do i=1,p+1
        write(*,*)(beta(i,j),j=1,r)
      enddo
      write(*,*)
      write(*,*)"alpha"
      do i=1,p
          write(*,*)(alpha(i,j),j=1,r)
      enddo
      write(*,*)"coverr"
      do i=1,p
          write(*,*)(coverr(i,j),j=1,p)
      enddo
      if (eri.eq.1) then
        write(*,*)'No mixture'
      elseif (eri.EQ.2) then
        write(*,*)'prob', pr
        write(*,*)'const', const
      else
      continue 
      endif
      write(*,*)  
      write(*,*)
      do i=1,p+1
        write(*,*)(REAL(tab1t(i,j))/REAL(esim),j=1,p+1)       
      enddo

      END
   




   















