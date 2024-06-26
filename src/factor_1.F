c      subroutine factor(N)
c **********************************************************************
c *                                                                    *
c *     CALCULATION OF FACTORIALS                                      *
c *     FAC(N) = 1*2*3*...*(N-1)                                       *
c *     DAC(N) = 1*3*5*...*(2*(N-1)-1)                                 *
c *                                                                    *
c **********************************************************************
c     implicit real*8 (a-h,o-z)
c      complex*16 cim, ccim
c      common/b1/cim,ccim,c(16,16),fac(13),dac(13),ngln(10),
c     1alamda,twopi,fpl30m,sws,vol,soo,swsl,fpilam,rmax,gmax,icnvrg,
c     2numvr,numvg,nq,nro,nl2,nl2p,klmo,nl,nlm,lmq,nowrt,idim


      subroutine factor(n)
c **********************************************************************
c *                                                                    *
c *     CALCULATION OF FACTORIALS                                      *
c *     FAC(N) = 1*2*3*...*(N-1)                                       *
c *     DAC(N) = 1*3*5*...*(2*(N-1)-1)                                 *
c *                                                                    *
c **********************************************************************
      IMPLICIT NONE
      INTEGER :: n, i, np
      real*8 fac(13), dac(13)
c      COMMON /b1/ fac, dac
c      real*8 fac(13), dac(13) 
      n = n + 1
      fac(1) = 1.d0
      fac(2) = 1.d0
      do i = 3, np
         fac(i) = fac(i-1) * (i-1)
      end do

      dac(1) = 1.d0
      dac(2) = 1.d0
      do i = 3, np
         dac(i) = dac(i-1) * (2 * i - 3)
      end do

      end

      
c      IMPLICIT NONE                                                                                                                                                                                      c     INTEGER :: n, i,np                                                                                                                                                                                  
c      DOUBLE PRECISION fac(13), dac(13)                                                                                                                                                                 
c      COMMON /b1/ fac, dac    

c      np=n+1
c      fac(1)=1.d0
c      fac(2)=1.d0
c      do 10 i=3,np
c   10 fac(i)=fac(i-1)*(i-1)
c      dac(1)=1.d0
c      dac(2)=1.d0
c      do 11 i=3,np
c   11 dac(i)=dac(i-1)*(2*i-3)
c      end





      
