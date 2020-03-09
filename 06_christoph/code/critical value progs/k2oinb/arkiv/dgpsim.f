      SUBROUTINE dgpsim(n,np,m,mp,k,kp,p,pp,r,
     +beta,alpha,gaen,mu,exog,gaex,iniobs,err,idum,endog)
      INTEGER n,np,m,mp,p,pp,k,kp,r,idum
      REAL err(pp,np),endog(pp,np),exog(mp,np),beta(pp,pp),
     +alpha(pp,pp),gaen(pp,(kp-1)*pp),mu(pp),chdia(r),
     +gaex(pp,mp),iniobs(pp,kp),bepbe(r,r),ibepbe(r,r),sum,
     +bebar(p,r) 
C     Simulates observations from a VAR GDP with restrictions
C     on linear term
C     k initial observations given
C     Output: Simulated values are returned in endog  
C     Input:   
C        p = dimension of VAR
C        k = order of VAR
C        m = number of exogenous variables (includes constant!!)
C        n = length of series
C        pi,ga,mu=parameters in DGP
C        ini = pxk matrix with initial observations
C        exog = mxn matrix containing exogenous variables
C        gaex = pxm matrix of coefficients of  exogenous variables
C        err = px(n-k) matrix of errors
      INTEGER i,j,l
      REAL z(pp*(kp-1)+mp+4),
     +mu1(p),pi(p,p)
      do i=1,p
        do j=1,p
          pi(i,j)=0.0
          do l=1,r
            pi(i,j)=pi(i,j)+alpha(i,l)*beta(j,l)
          enddo
        enddo
      enddo
      do i=1,p
        mu1(i)=0
        do l=1,r
          mu1(i)=mu1(i)+alpha(i,l)*beta(p+1,l)
        enddo
      enddo
      do i=1,k     
        do j=1,p         
          endog(j,i)=iniobs(j,i)
        enddo
      enddo   

C Starts simulating
      do i=k+1,n
        do j=1,p
          endog(j,i)=0.0
        enddo      
        do j=1,k-1 
          do l=1,p
            z((j-1)*p+l)=(endog(l,i-j)-endog(l,i-j-1))
          enddo
        enddo
        do j=1,p
          endog(j,i)=endog(j,i-1)
          endog(j,i)= endog(j,i) +(i-k)* mu1(j)         
        enddo 
        do j=1,p
          do l=1,p
            endog(j,i)=endog(j,i)+ pi(j,l)*endog(l,i-1)
          enddo
        enddo  
        do j=1,p 
          do l=1,p*(k-1)
            endog(j,i)= endog(j,i)+gaen(j,l)*z(l)
          enddo
        enddo
        do j=1,p
          do l=1,m
            endog(j,i)= endog(j,i)+gaex(j,l)*exog(l,i)
          enddo
          endog(j,i)= endog(j,i) + err(j,i)
        enddo
      enddo 

      return

      END  

     







