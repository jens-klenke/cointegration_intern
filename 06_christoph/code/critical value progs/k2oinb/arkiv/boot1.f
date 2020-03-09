      PROGRAM boot 
C Program for bootstrapping simulated VAR modell.
C Model with linear term with restrictions. Use of OLS
C estimates for short run parameters. Computes
C p-values of trace statistics and lambda_max statistics
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
     +ilamax,itrce,iala(6),iatr(6),iitrce(6),iilama(6),
     +dim,eri,cnt,esim,resind 
      REAL endog(PD,NNOBS),exog(MD,NNOBS),ini(PD,KP),eerr(PD,NNOBS),
     +s00(PD,PD),s01(PD,PD),s10(PD,PD),s11(PD,PD),
     +beta(PD,PD),lambda(PD),gaen(PD,PD*KP),
     +alpha(PD,PD),mu(PD),gaex(PD,MD),res(PD,NNOBS),bootres(PD,NNOBS),
     +sbeta(PD,PD),slambda(PD),sgaen(PD,PD*KP),
     +salpha(PD,PD),smu(PD),sgaex(PD,MD),
     +simob(PD,NNOBS),pi(PD,PD),bs00(PD,PD),bs01(PD,PD),
     +bs10(PD,PD),bs11(PD,PD),bbeta(PD,PD),blam(PD),
     +balph(PD,PD),lamax,trce,btrce,
     +nomlev(6),pr,const,coverr(PD,PD),asqla(6,4),asqtr(6,4),cov2(PD,PD)
      nobs=50
       nsim=100
       esim=nsim
      nboot=100
      m=1    
      do i=1,nobs
        exog(1,i)=1.0
      enddo
      idum=-3
      nomlev(1)=0.5
      nomlev(2)=0.2
      nomlev(3)=0.1
      nomlev(4)=0.05
      nomlev(5)=0.025
      nomlev(6)=0.01
      ilamax=0
      itrce=0
      asqla(1,1)=5.55
      asqla(2,1)=8.65
      asqla(3,1)=10.49
      asqla(4,1)=12.25
      asqla(5,1)=14.21
      asqla(6,1)=16.26
      asqla(1,2)=10.90
      asqla(2,2)=14.70
      asqla(3,2)=16.85
      asqla(4,2)=18.96
      asqla(5,2)=21.14
      asqla(6,2)=23.65
      asqla(1,3)=16.24
      asqla(2,3)=20.45
      asqla(3,3)=23.11
      asqla(4,3)=25.54
      asqla(5,3)=27.68
      asqla(6,3)=30.34
      asqla(1,4)=21.50
      asqla(2,4)=26.30
      asqla(3,4)=29.12
      asqla(4,4)=31.46
      asqla(5,4)=33.6
      asqla(6,4)=36.65
      asqtr(1,1)=5.55
      asqtr(2,1)=8.65
      asqtr(3,1)=10.49
      asqtr(4,1)=12.25
      asqtr(5,1)=14.21
      asqtr(6,1)=16.26
      asqtr(1,2)=15.59
      asqtr(2,2)=20.19
      asqtr(3,2)=22.76
      asqtr(4,2)=25.32
      asqtr(5,2)=27.75
      asqtr(6,2)=30.45
      asqtr(1,3)=29.53
      asqtr(2,3)=35.56
      asqtr(3,3)=39.06
      asqtr(4,3)=42.44
      asqtr(5,3)=45.42
      asqtr(6,3)=48.45
      asqtr(1,4)=47.17
      asqtr(2,4)=54.80
      asqtr(3,4)=59.14
      asqtr(4,4)=62.99
      asqtr(5,4)=66.25
      asqtr(6,4)=70.05
      do i=1,6
        iitrce(i)=0
        iilama(i)=0
        iala(i)=0
        iatr(i)=0
      enddo

C Parameters of DGP VAR 2
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
      alpha(1,1)=-0.1
      alpha(2,1)=-0.1
      beta(1,1)=1.0
      do i=1,p
        do j=1,p
          coverr(i,j)=0.0      
        enddo
      coverr(i,i)=1.0
      enddo
C Type of residuals, resind=1: OLS,resind=2: restricted ?
      resind=1
C Any mixture of errors?
C er1=1: Standard normal, eri=2:Mixture of normals
      eri=2
      pr=0.9
      const=10.0 
      do i=1,p
        do j=1,p
          cov2(i,j)=(pr+const*const*(1-pr))*coverr(i,j)
        enddo
      enddo
      dim=r+p*k

      if (eri.EQ.2) then
        do i=1,p
          do j=1,p         
            coverr(i,j)=cov2(i,j)
          enddo
        enddo
      endif

C Simulates NOSIM realizations of DGP
      do i=1,nsim
        call dgpini(k,KP,p,PD,r,dim,
     +beta,alpha,gaen,mu,ini,idum,coverr)
      do j=1,p
      enddo 
        call errsim(nobs,NNOBS,k,KP,p,PD,eri,pr,const,idum,eerr,coverr)
        do j=1,p
        enddo
        call dgpsim(nobs,NNOBS,m,MD,k,KP,p,PD,r,
     +beta,alpha,gaen,mu,exog,mu,ini,eerr,idum,endog)
C Computes estimates of parameters in simulated realizations 
        call rempmat2(endog,exog,nobs,NNOBS,k,KP,p,PD,m,MD,s00,s01,
     +s10,s11)
        call rrest2(s00,s01,s11,p,PD,lambda,sbeta,salpha)
         if (resind.EQ.1) then
           call rrlsq1(endog,exog,nobs,NNOBS,k,KP,p,PD,m,mD,p,sbeta,
     +  sgaen,salpha,sgaex,res)
         elseif (resind.EQ.2) then
           call rrlsq1(endog,exog,nobs,NNOBS,k,KP,p,PD,m,mD,r,sbeta,
     +  sgaen,salpha,sgaex,res)
         else
           continue
         endif
        do l1=1,p+1
        enddo
         do l1=1,p
         enddo

        if(r.ge.p)then
        else
          lamax=-(nobs-k)*log(1-lambda(r+1))
          trce=lamax
          do j=r+2,p
            trce= trce-(nobs-k)*log(1-lambda(j))
          enddo
        endif
        do j=1,6
          if(lamax.GT.asqla(j,4))iala(j)=iala(j)+1
          if(trce.GT.asqtr(j,4))iatr(j)=iatr(j)+1
        enddo
         do l1=1,p
          enddo

         cnt=0
         call outcnt(k,KP,p,PD,r,
     +sbeta,salpha,sgaen,cnt)
         if (cnt.GT.0) then 
           esim=esim-1
           goto 10
         endif        
C     ******  Integer r determines rank of 
C     impactmatrix used for constructing pseudoobsercvations.
      itrce=0
      ilamax=0
        do j=1,nboot
        call berrsim(res,nobs,NNOBS,k,KP,p,PD,idum,bootres)
          do l1=1,p
          enddo
        call dgpsim(nobs,NNOBS,m,MD,k,KP,p,PD,r,
     +sbeta,salpha,sgaen,sgaex,exog,sgaex,ini,bootres,idum,simob)
          do l1=1,p
          enddo         
          call rempmat2(simob,exog,nobs,NNOBS,k,KP,p,PD,m,MD,     
     +bs00,bs01,bs10,bs11)
          call rrest2(bs00,bs01,bs11,p,PD,blam,bbeta,balph)
           btrce=0.0
          do l=r+1,p
            btrce= btrce-(nobs-k)*log(1.0-blam(l))
          enddo
        if(btrce.GT.trce)itrce=itrce+1 
        if(-(nobs-k)*log(1.0-blam(r+1)).GT.lamax)ilamax=ilamax+1 
        enddo
C end of bootstrap loop
        do j=1,6
          if(REAL(itrce)/REAL(nboot).LE.nomlev(j))iitrce(j)=iitrce(j)+1
          if(REAL(ilamax)/REAL(nboot).LE.nomlev(j))iilama(j)=iilama(j)+1
        enddo
   10 continue    
      enddo
C end of simulation loop
      if (resind.EQ.1) then
        write(*,*)'Ordinary residuals bootstrapped'
      elseif (resind.EQ.2) then
        write(*,*)'Restricted  residuals bootstrapped, rank:',r
      else
        continue
      endif
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
      write(*,*)"Nom. lev and emp. level (boot. and as.tr. test)"
      do i=1,6
        write(*,'(1X,3F10.2)')nomlev(i),REAL(iitrce(i))/REAL(esim),
     +  REAL(iatr(i))/REAL(nsim)
      enddo
      write(*,*)"Nom. lev. and emp. level (boot. and as. l-m. test)"
      do i=1,6
        write(*,'(1X,3F10.2)')nomlev(i),REAL(iilama(i))/REAL(esim),
     +  REAL(iala(i))/REAL(nsim)
      enddo

      END
   




   















