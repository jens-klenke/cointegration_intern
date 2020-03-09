      SUBROUTINE rrest2(s00,s01,s11,p,pp,lambda,beta,alpha)
      INTEGER p,pp
      REAL s00(pp,pp),s01(pp,pp),s11(pp,pp),lambda(pp),beta(pp,pp),
     +alpha(pp,pp)
C     Solves a generalised eigenvalue problem for a reduced
C     rank VAR-modell WITH restrictions on the constant term,
C     or WITH restrictions on a linear term and no  
C     restrictions on a constant term. 

      INTEGER i,j,k,l,nrot1,nrot2
      REAL chdia(p),sum,s00inv(p+1,p+1),a0(pp,pp),b(pp,pp),diag1(pp),
     +diag2(pp),w(pp,pp),bisqr(p+1,p+1),a1(pp,pp),v(pp,pp),v2(pp,pp)

C Inversion of s00
      call choldc(s00,p,pp,chdia)
      do i=1,p
        s00(i,i)=1.0/chdia(i)
        do j=i+1,p
          sum=0.0
          do l=i,j-1
            sum=sum-s00(j,l)*s00(l,i)
          enddo
        s00(j,i)=sum/chdia(j)
        enddo
      enddo
      do i=1,p
        do j=i+1,p
          s00(i,j)=0.0
        enddo
      enddo
      do i=1,p
        do j=1,p
          s00inv(i,j)=0.0
        enddo
      enddo 
      do i=1,p
        do j=1,p
          do l=1,p
            s00inv(i,j)=s00inv(i,j)+s00(l,i)*s00(l,j)
          enddo
        enddo
      enddo 

C Computes a0=s10*s00inv*s01
      do i=1,p+1
        do j=1,p+1
          a0(i,j)=0.0
        enddo
      enddo
      do i=1,p+1
        do j=1,p+1
          do k=1,p
            do l=1,p
              a0(i,j)=a0(i,j)+s01(k,i)*s00inv(k,l)*s01(l,j)
            enddo
          enddo
        enddo
      enddo 

C Computes s11=w*diag-1/2*w'
      do i=1,p+1
        do j=1,p+1
          b(i,j)=s11(i,j)
        enddo
      enddo
      call jacobi(b,p+1,pp,diag1,w,nrot1)
      call eigsrt(diag1,w,p+1,pp)
      do i=1,p+1
        do j=1,p+1
          bisqr(i,j)=0.0
        enddo
      enddo
      do i=1,p+1
        do  j=1,p+1
          do  k=1,p+1
            bisqr(i,j)=bisqr(i,j)+w(i,k)*(1/sqrt(diag1(k)))*w(j,k)
          enddo
        enddo
      enddo
      do i=1,p+1
        do j=1,p+1
          a1(i,j)=0.0
        enddo
      enddo 
      do i=1,p+1
        do j=1,p+1
          do k=1,p+1
            do l=1,p+1
              a1(i,j)=a1(i,j)+bisqr(i,k)*a0(k,l)*bisqr(l,j)
            enddo
          enddo
        enddo
      enddo 

      call jacobi(a1,p+1,pp,lambda,v,nrot2)
      call eigsrt(lambda,v,p+1,pp)

C Computes beta=(bsqrt)*v
      do i=1,p+1
        do j=1,p
          beta(i,j)=0.0
        enddo
      enddo 
      do i=1,p+1
        do j=1,p
          do k=1,p+1
            beta(i,j)=beta(i,j)+ bisqr(i,k)*v(k,j)
          enddo
        enddo
      enddo

C Computes alpha=s01*beta
      do i=1,p
        do j=1,p
          alpha(i,j)=0.0
        enddo
      enddo 
      do i=1,p
        do j=1,p
          do k=1,p+1
            alpha(i,j)=alpha(i,j)+ s01(i,k)*beta(k,j)
          enddo
        enddo
      enddo


      return 

      END  

     




