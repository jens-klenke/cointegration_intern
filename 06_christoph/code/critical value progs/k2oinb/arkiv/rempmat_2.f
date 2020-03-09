      SUBROUTINE rempmat2(endog,exog,n,np,k,kp,p,pp,m,md,s00,s01,
     +s10,s11)
      INTEGER n,np,p,pp,k,kp,m,md
      REAL endog(pp,np),exog(md,np),s00(pp,pp),s01(pp,pp),
     +s10(pp,pp),s11(pp,pp) 
C     Computes empirical moment matrices for estimation of  parameters 
C     of a reduced rank VAR model with restrictions on linear term,
C     which excludes a quadratic trend.
C        p=dimension of VAR
C        k=order of VAR
C        r=reduced rank
C        m=number of exogenous variables
C        c= 0 if no constant term, 1 if unrestrictred
C        d= 0 if no quarterly seasonal dummies, 1 otherwise
C        n=length of series
C        endo= nxp matrix containing endogenous variables
C        exo= nxm matrix containing exogenous variables
      INTEGER i,j,l,l1,dim,nrot
      DOUBLE PRECISION sdum1(n+4),sdum2(n+4),sdum3(n+4),m00(p,p),
     +m22(k*p+m,k*p+m),m20(k*p+m,p),m02(p,k*p+m),
     +m11(p+1,p+1),
     +m10(p+1,p),m01(p,p+1),m12(p+1,k*p+m),m21(k*p+m,p+1),
     +chdia(k*p+m),z(p*(k-1)+m),dendog(p,n),
     +sum,m22inv(k*p+m,k*p+m),ds00(p,p),ds01(p,p+1),
     +ds10(p+1,p),ds11(p+1,p+1)
      do i=1,n
        do j=1,p
          dendog(j,i)=dble(endog(j,i))
        enddo
      enddo
      do i=1,p
        do j=1,p
          m00(i,j)=0.0
          m10(i,j)=0.0
          m01(i,j)=0.0
          m11(i,j)=0.0
          ds00(i,j)=0.0
          ds10(i,j)=0.0
          ds01(i,j)=0.0
          ds11(i,j)=0.0
        enddo
        m11(i,p+1)=0.0
        m01(i,p+1)=0.0
        m10(p+1,i)=0.0
        ds11(i,p+1)=0.0
        ds11(p+1,1)=0.0
        ds01(i,p+1)=0.0
        ds10(p+1,i)=0.0
      enddo
      m11(p+1,p+1)=0.0
      ds11(p+1,p+1)=0.0
      do i=1,k*p+m
        do j=1,k*p+m
          m22(i,j)=0.0
          m22inv(i,j)=0.0
        enddo
      enddo
      do i=1,k*p+m
        do j=1,p
          m20(i,j)=0.0
          m21(i,j)=0.0
        enddo
        m21(i,p+1)=0.0
      enddo
      do i=1,p
        do j=1,k*p+m
          m02(i,j)=0.0
          m12(i,j)=0.0
        enddo
      enddo
      do i=1,k*p+m
        m12(p+1,i)=0.0
      enddo
      dim=p*(k-1)+m
C Starts summing over observations
      do i=k+1,n
        do j=1,k-1
          do l=1,p
            z(l+p*(j-1))= endog(l,i-j)- endog(l,i-j-1)
          enddo
        enddo
        do j=1,m
             z(p*(k-1)+j)= exog(j,i)
        enddo
C Computes the empirical moment matrices. Notation as in Johansen (1996)
        do j=1,dim
          do l=1,dim
          m22(j,l)= m22(j,l)+ z(j)*z(l)/(n-k)
          enddo
        enddo
        do j=1,p
          do l=1,p
          m00(j,l)= m00(j,l)+ (endog(l,i)- endog(l,i-1))*
     +    (endog(j,i)- endog(j,i-1))/(n-k) 
          enddo
        enddo
        do j=1,p
          do l=1,dim
            m02(j,l)= m02(j,l) + (endog(j,i)- endog(j,i-1))*z(l)/(n-k)
            m20(l,j)= m02(j,l)
          enddo
        enddo    
        do j=1,p
          do l=1,p
            m11(j,l)= m11(j,l) + endog(l,i-1)*endog(j,i-1)/(n-k)
          enddo
          m11(j,p+1)=  m11(j,p+1)+ (i)*endog(j,i-1)/(n-k)
          m11(p+1,j)=  m11(p+1,j)+ (i)*endog(j,i-1)/(n-k)
        enddo
        m11(p+1,p+1)=  m11(p+1,p+1)+
     +    (dfloat((i)*(i)))/dfloat(n-k)
        do j=1,p
          do l=1,p
          m01(j,l)= m01(j,l) + (endog(j,i)- endog(j,i-1))*
     +    endog(l,i-1)/(n-k)
          m10(l,j)=  m01(j,l)
          enddo
        m01(j,p+1)= m01(j,p+1) + (endog(j,i)- endog(j,i-1))*
     +  (i)/(n-k)
        m10(p+1,j)=  m01(j,p+1)
        enddo
        do j=1,p
          do l=1,dim
            m12(j,l)= m12(j,l)+ (endog(j,i-1))*z(l)/(n-k)
            m21(l,j)= m12(j,l)
          enddo
        enddo
        do l=1,dim
            m12(p+1,l)= m12(p+1,l)+ (i)*z(l)/(n-k)
            m21(l,p+1)= m12(j,l)
          enddo
      enddo
C End of construction of moment matrices

C Inversion of M22
      call dchldc(m22,dim,k*p+m,chdia)
      do i=1,dim
        m22(i,i)=1.0/chdia(i)
        do j=i+1,dim
          sum=0.0
          do l=i,j-1
            sum=sum-m22(j,l)*m22(l,i)
          enddo
        m22(j,i)=sum/chdia(j)
        enddo
      enddo 
      do i=1,dim
        do j=i+1,dim
          m22(i,j)=0.0
        enddo
      enddo 
      do i=1,dim
        do j=1,dim
          do l=1,dim
            m22inv(i,j)=m22inv(i,j)+m22(l,i)*m22(l,j)
          enddo
        enddo
      enddo 
C end of inversion of M22

C computation of S00, S01, S11.
      do i=1,p
        do j=1,p
        ds00(i,j)=m00(i,j)
        ds11(i,j)=m11(i,j)
        ds01(i,j)=m01(i,j)
        ds10(i,j)=m10(i,j)
          do l=1,dim
            do l1=1,dim         
              ds00(i,j)=ds00(i,j)-m02(i,l)*m22inv(l,l1)*m20(l1,j)
              ds11(i,j)=ds11(i,j)-m12(i,l)*m22inv(l,l1)*m21(l1,j)
              ds01(i,j)=ds01(i,j)-m02(i,l)*m22inv(l,l1)*m21(l1,j)
              ds10(i,j)=ds10(i,j)-m12(i,l)*m22inv(l,l1)*m20(l1,j)
            enddo
          enddo
          s00(i,j)=real(ds00(i,j))
          s11(i,j)=real(ds11(i,j))
          s01(i,j)=real(ds01(i,j))
          s10(i,j)=real(ds10(i,j))
        enddo
      enddo 
      do i=1,p
        ds11(i,p+1)=m11(i,p+1)
        ds01(i,p+1)=m01(i,p+1)
        ds10(p+1,i)=m10(i,p+1)
        do l=1,dim
          do l1=1,dim         
            ds11(i,p+1)=ds11(i,p+1)-m12(i,l)*m22inv(l,l1)*m21(l1,p+1)
            ds01(i,p+1)=ds01(i,p+1)-m02(i,l)*m22inv(l,l1)*m21(l1,p+1)
          enddo
        enddo
          ds10(p+1,i)=ds01(i,p+1)
          ds11(p+1,i)=ds11(i,p+1)
          s00(i,p+1)=real(ds00(i,p+1))
          s11(i,p+1)=real(ds11(i,p+1))
          s11(p+1,i)=s11(i,p+1)
          s01(i,p+1)=real(ds01(i,p+1))
          s10(p+1,i)=real(ds10(p+1,i))
      enddo
      ds11(p+1,p+1)=m11(p+1,p+1)
      do l=1,dim
        do l1=1,dim         
          ds11(p+1,p+1)=ds11(p+1,p+1)-m12(p+1,l)*m22inv(l,l1)*
     +m21(l1,p+1)
        enddo
      enddo
      s11(p+1,p+1)=real(ds11(p+1,p+1))

      return 

      END  

     




