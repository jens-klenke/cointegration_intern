      SUBROUTINE dgpini(k,kp,p,pp,r,dim,
     +beta,alpha,ga,mu,endog,idum,coverr)
      INTEGER p,pp,k,kp,r,idum,dim
      REAL beta(pp,pp),endog(pp,kp),coverr(pp,pp),
     +alpha(pp,pp),ga(pp,(kp-1)*pp),mu(pp)
C     Computes k initial observations 
C     for simulation fram a VAR DGP with constant
C     and restricted linear term.
C     Output: k p-vectors in endog  
C     Input:   
C        p = dimension of VAR
C        k = order of VAR
C        dim=r+p*(k-1)
C        alpha,beta,ga,mu=parameters in DGP
      INTEGER i,j,l,m,dim2,dim3
 
      REAL sum,amat(dim,dim),amati(dim,dim),ini(dim),
     +aamat(dim*dim,dim*dim),chdia(dim*dim),inicon(dim),
     +aamati(dim*dim,dim*dim),chdia2(dim),indx(dim*dim),inimn(dim),
     +eps(dim),veci(dim*dim),vecx(dim*dim),cov(dim,dim),d,d2,
     +indx2(dim),chdia3(dim),pre(p+r,p),
     +temp1(dim*dim,dim*dim),temp2(dim*dim,dim*dim),
     +temp3(dim*dim,dim*dim),temp4(dim,dim),
     +covsh(dim,dim),bepbe(r,r),ibepbe(r,r),bebar(p,r) 
      dim2=dim*dim
C Initiating
      if (k.EQ.1) then   
        dim3=r+k*p
      else
        dim3=r+(k-1)*P
      endif 
      do i=1,dim3
        do j=1,dim3
          amat(i,j)=0.0
          aamat(i,j)=0.0          
        enddo         
      enddo

C first r rows
      do i=1,r
        do j=1,r
          do l=1,p
            amat(i,j)=amat(i,j)+beta(l,i)*alpha(l,j)
          enddo 
        if (i.eq.j) amat(i,j)=amat(i,j)+ 1
        enddo         
      enddo
      do m=1,k-1
        do i=1,r
          do j=1,p     
            do l=1,p    
              amat(i,r+(m-1)*p+j)=amat(i,r+(m-1)*p+j)+
     +        beta(l,i)*ga(l,(m-1)*p+j)
            enddo
          enddo
        enddo
      enddo
   
C next p  rows
      do i=r+1,r+p
        do j=1,r
          amat(i,j)=alpha(i-r,j)
        enddo
        do m=1,k-1
          do j=1,p
            amat(i,r+(m-1)*p+j)=ga(i-r,(m-1)*p+j)
          enddo
        enddo 
      enddo
C rest of matrix
        do i=r+p,dim3
          do j=r+1,r+(k-2)*p
            if (i.eq.j+p)amat(i,j)=1
          enddo
        enddo
 
C     computes I-Kron.prod(A,A)
      do i=1,dim3   
        do j=1,dim3
          do l=1,dim3
            do m=1,dim3
              aamat(i+dim3*(l-1),j+dim3*(m-1))=amat(i,j)*amat(l,m)
            enddo
          enddo
        enddo
      enddo
      do i=1,dim2
        do j=1,dim2
           temp1(i,j)=aamat(i,j)
        enddo
      enddo
      do i=1,dim2
        do j=1,dim2
          aamat(i,j)=-aamat(i,j)
          if (i.eq.j) aamat(i,j)=aamat(i,j)+1.0
          temp1(i,j)=aamat(i,j) 
        enddo   
      enddo

C     computes covariance matrix
C     
      do i=1,p+r
        do j=1,p
          pre(i,j)=0.0
        enddo
      enddo
      do i=1,r
        do j=1,p
          pre(i,j)=beta(j,i)
        enddo
      enddo
      do i=r+1,r+p
        pre(i,i-r)=1.0
      enddo
      do i=1,dim3
        do j=1,dim3
          covsh(i,j)=0.0
        enddo
      enddo
      do i=1,r+p
        do j=1,r+p
          do l=1,p
            do m=1,p
              covsh(i,j)=covsh(i,j)+pre(i,l)*coverr(l,m)*pre(j,m)
            enddo        
          enddo 
        enddo
      enddo   
      do i=1,dim3
        do j=1,dim3
          veci(i+dim3*(j-1))=covsh(i,j)
        enddo
      enddo
      call ludcmp(aamat,dim2,dim2,indx,d)
      call lubksb(aamat,dim2,dim2,indx,veci)
      do i=1,dim3
        do j=1,dim3
          cov(j,i)=veci(j+(i-1)*dim3)
        enddo
      enddo

C computes constant in autoregression
      do i=1,dim3
        inicon(i)=0.0
        inimn(i)=0.0
      enddo 

      do i=1,r
        do j=1,p
          inicon(i)=inicon(i)+beta(j,i)*mu(j)
        enddo
        inicon(i)=inicon(i)+beta(p+1,i)
      enddo  
      do i=r+1,r+p
        inicon(i)=mu(i-r)
      enddo



C     computes stationary expectation
C     inverts (I-A)
      do i=1,dim3
        do j=1,dim3
          amat(i,j)=-amat(i,j)
          if (i.eq.j) amat(i,j)=amat(i,j)+1.0
        enddo
      enddo
      do i=1,dim3
        do j=1,dim3
          temp3(i,j)=amat(i,j)
        enddo
      enddo
      do i=1,dim3
        do j=1,dim3
          amati(i,j)=0.0
        enddo
        amati(i,i)=1.0  
      enddo
      call ludcmp(amat,dim3,dim3,indx2,d2)
      call lubksb(amat,dim3,dim3,indx2,inicon)
      do i=1,dim3   
        do j=1,dim3
         temp4(i,j)=0.0
          do l=1,dim3
            temp4(i,j)=temp4(i,j)+temp3(i,l)*amati(l,j)
          enddo
        enddo
      enddo 

C Starts computing the initial distribution
      do i=1,dim3
        chdia2(i)=0.0
      enddo
      do i=1,dim3
        do j=1,dim3
          amat(i,j)=0.0
        enddo         
      enddo
      call choldc(cov,dim3,dim3,chdia2)
      do i=1,dim3
        cov(i,i)=chdia2(i)
        eps(i)=gasdev(idum)
        ini(i)=0.0
      enddo
      do i=1,dim3
        do j=i+1,dim3
          cov(i,j)=0.0
        enddo
      enddo
      do i=1,dim3
        do j=1,dim3
          ini(i)=ini(i)+cov(i,j)*eps(j)
        enddo
        ini(i)=ini(i)+inicon(i)
      enddo
      do i=1,r
        do j=1,r
          bepbe(i,j)=0.0
          ibepbe(i,j)=0.0
          do l=1,p
            bepbe(i,j)=bepbe(i,j)+beta(l,i)*beta(l,j)
          enddo
        enddo
      enddo

C Computes initial observations for VAR
C Inversion of bepbe
      call choldc(bepbe,r,r,chdia3)
      do i=1,r
        bepbe(i,i)=1.0/chdia3(i)
        do j=i+1,r
          sum=0.0
          do l=i,j-1
            sum=sum-bepbe(j,l)*bepbe(l,i)
          enddo
          bepbe(j,i)=sum/chdia3(j)
        enddo
      enddo 
      do i=1,r
        do j=1,i-1
          bepbe(j,i)=0.0
        enddo
      enddo           
      do i=1,r
        do j=1,r
          do l=1,r
            ibepbe(i,j)=ibepbe(i,j)+bepbe(l,i)*bepbe(l,j)
          enddo
        enddo
      enddo
      do i=1,p
        do j=1,r
          bebar(i,j)=0.0
          do l=1,r
            bebar(i,j)=bebar(i,j)+beta(i,l)*ibepbe(l,j)
          enddo
        enddo
      enddo 
C Initiating
      do i=1,p
        endog(i,k)=0.0
        do j=1,r
          endog(i,k)=endog(i,k)+bebar(i,j)*ini(j)
        enddo
      enddo
      do i=1,k-1
        do j=1,p
          endog(j,k-i)=endog(j,k-i+1)-ini(r+(i-1)*p+j)
        enddo
      enddo

      return

      END  

     







