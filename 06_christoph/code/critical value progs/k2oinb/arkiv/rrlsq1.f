      SUBROUTINE rrlsq1(endog,exog,n,np,k,kp,p,pp,m,md,r,beta,gaen,
     +alpha,gaex,res)
      INTEGER n,np,p,pp,k,kp,m,md,r
      REAL endog(pp,np),exog(md,np),beta(pp,pp),alpha(pp,pp),
     +gaen(pp,(kp-1)*pp),gaex(pp,md),res(pp,np) 
C     Computes the least squars esimators forestimation of  parameters 
C     of a reduced rank VAR model with restricted linear term. 
C     beta ((p+1)xp) contains the cointegration 
C     vectors
C        p=dimension of VAR
C        k=order of VAR
C        m=number of exogenous variables
C        r=rank impact matrix
C        c= 0 if no constant term, 1 if unrestrictred,
C        d= 0 if no quarterly seasonal dummies, 1 otherwise
C        n=length of series
C        endo= nxp matrix containing endogenous variables
C        exo= nxm matrix containing exogenous variables
      INTEGER i,j,l,l1,dim,nrot
      REAL z(np,pp*kp+md+1),
     +v(pp*kp+md+1,pp*kp+md+1),w(pp*kp+md+1),pi(p,p+1),
     +b(np),x(pp*kp+md+1),sum,par(pp*kp+md+1,pp),thres,TOL,wmax,
     +a(np,pp*kp+md+1)
      PARAMETER(TOL=0.00001)   
      do i=1,p
        do j=1,p+1
          pi(i,j)=0.0
          do l=1,r
            pi(i,j)= pi(i,j)+alpha(i,l)*beta(j,l)
          enddo
        enddo
      enddo
C Starts summing over observations
      dim=p*(k-1)+m
      do i=k+1,n
        do j=1,k-1
          do l=1,p
            z(i-k,l+p*(j-1))= endog(l,i-j)- endog(l,i-j-1)
          enddo
        enddo
        do j=1,m
          z(i-k,j+p*(k-1))=exog(j,i)
        enddo
      enddo
      do i=1,n-k
        do j=1,dim
          a(i,j)=z(i,j)
        enddo
      enddo
C End of construction of design matrix.
C
C
C Singular value decomposition
      call svdcmp(z,n-k,dim,np,pp*kp+md+1,w,v)
C Delets small singular values
      wmax=0.0
      do i=1,dim  
        if(w(i).gt.wmax)wmax=w(i)   
      enddo
      thres=TOL*wmax
      do i=1,dim  
        if(w(i).lt.thres)w(i)=0.0   
      enddo
C Find estimated coefficient

      do i=1,p
        do j=k+1,n
          b(j-k)=endog(i,j)-endog(i,j-1)
          do l=1,p
            b(j-k)=b(j-k)-pi(i,l)*endog(l,j-1)
          enddo
          b(j-k)=b(j-k)-pi(i,p+1)*j     
        enddo
      call svbksb(z,w,v,n-k,dim,np,pp*kp+md+1,b,x)
        do j=1,dim
          par(j,i)=x(j)
        enddo
        do j=1,n-k
          sum=0.0
          do l=1,dim
            sum=sum+a(j,l)*x(l)
          enddo
          res(i,j)=b(j)-sum
        enddo  
      enddo
      do i=1,p
        do j=1,p*(k-1)
          gaen(i,j)=par(j,i)
        enddo
      enddo
      do i=1,p
        do j=p*(k-1)+1,p*(k-1)+m
          gaex(i,j-p*(k-1))=par(j,i)
        enddo
      enddo     

      return

      END  

     




