      SUBROUTINE outcnt(k,kp,p,pp,r,
     +beta,alpha,ga,cnt)
      INTEGER p,pp,k,kp,r,cnt
      REAL beta(pp,pp),
     +alpha(pp,pp),ga(pp,(kp-1)*pp)
C     Computes the number of roots of the characteristic 
C     polynomials which are not eqal to one or outside
C     the unit  sirkel
C     Output: number of roots   
C     Input:   
C        p = dimension of VAR
C        k = order of VAR
C        r = rank of Pi matrix
C        alpha,beta,ga,mu=parameters of VAR
      INTEGER i,j,l,dim
 
      REAL amat(k*p,k*p),lamrl(k*p),lamim(k*p),
     +pi(p,p),aux1,aux2,tol
      tol=0.0001
      dim=k*p
      do i=1,dim
        do j=1,dim
          amat(i,j)=0.0
        enddo  
      enddo
      do i=1,p
        do j=1,p
          pi(i,j)=0.0
          do l=1,r
            pi(i,j)=pi(i,j)+alpha(i,l)*beta(j,l)
          enddo         
        enddo  
      enddo
      do i=1,p
        if(k.GT.1) then
          do j=1,p
            amat(i,j)=ga(i,j)+pi(i,j)
          enddo
        else
          do j=1,p
            amat(i,j)=pi(i,j)
          enddo
        endif 
        amat(i,i)=amat(i,i)+1.0
        do j=1,k-2
          do l=1,p
            amat(i,j*p+l)=ga(i,j*p+l)-ga(i,(j-1)*p+l)
          enddo
        enddo
        if(k.GT.1) then 
          do j=1,p
            amat(i,(k-1)*p+j)=-ga(i,(k-2)*p+j)
          enddo
        endif
      enddo
      do i=1,(k-1)*p
          amat(p+i,i)=1.0
      enddo         
      do i=1,dim
C        write(*,*)(amat(i,j),j=1,dim)
      enddo


      call balanc(amat,dim,dim)
      call elmhes(amat,dim,dim)
      call hqr(amat,dim,dim,lamrl,lamim)
      cnt=0
      do i=1,dim
        aux1= lamrl(i)*lamrl(i)+lamim(i)*lamim(i)
        aux1=sqrt(aux1)-1.0
        aux2= (lamrl(i)-1)*(lamrl(i)-1)+lamim(i)*lamim(i)
        aux2=sqrt(aux2) 
        if(aux1.GT.-tol.AND.aux2.GT.tol)cnt=cnt+1
      enddo       


      return

      END  

     







