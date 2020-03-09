      SUBROUTINE berrsim(res,n,np,k,kp,p,pp,idum,bootres)
      INTEGER p,pp,k,kp,idum,n,np
      REAL res(pp,np),bootres(pp,np)
C     Samples from the  pxn matrix res  
C     to get pxn errors for bootstrapping
C     Input p(n-k) matrix of residuals from fit  
C     Output: pxn matrix of errors for bootstrapping
      INTEGER i,j

      do i=1,n
        sim=ran2(idum) 
        sim=float(n-k)*sim+1.0
        isim=INT(sim)
        do j=1,p
           bootres(j,i)=res(j,isim)
        enddo
      enddo

      return

      END  

     







