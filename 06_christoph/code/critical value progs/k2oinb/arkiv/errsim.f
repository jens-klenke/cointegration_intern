      SUBROUTINE errsim(n,np,k,kp,p,pp,eri,pr,const,idum,err,coverr)
      INTEGER p,pp,k,kp,eri
      REAL err(pp,np),pr,const,x,coverr(pp,pp),chdia(p),a(p,p),
     +eps(p),cov(p,p)
C     Computes the px(n-k) errors used for
C     simulating a VAR model.
C     Output: pxn matrix of errors   
C     Input: eri: Indicates type of distribution
C     1= standard normal
C     2= standard normal with pr, st.sev= const with (1-pr)  
      INTEGER i,j,l,m
      do i=1,p
        do j=1,p
          cov(i,j)=coverr(i,j)
        enddo
      enddo 
      call choldc(cov,p,p,chdia)
      do i=1,p
        do j=1,i-1
          a(i,j)=cov(i,j)
        enddo
      a(i,i)=chdia(i)  
        do j=i+1,p
          a(i,j)=0.0
        enddo
      enddo
      if (eri.eq.1) then
        do i=1,n
          do j=1,p
            eps(j)=gasdev(idum)
          enddo
          do j=1,p
            err(j,i)=0.0
            do l=1,p
              err(j,i)= err(j,i)+a(j,l)*eps(l)
            enddo
          enddo         
        enddo
      elseif (eri.eq.2) then
        do i=1,n
          do j=1,p
            eps(j)=gasdev(idum)
          enddo
          do j=1,p
            err(j,i)=0.0
            do l=1,p
              err(j,i)= err(j,i)+a(j,l)*eps(l)
            enddo
          enddo 
          x=ran2(idum)
          if (x.gt.pr) then
            do j=1,p 
            err(j,i)=const*err(j,i)
            enddo
          endif
        enddo         
      else
      endif
       

      return

      END  

     







