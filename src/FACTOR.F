      subroutine factor(N)
c **********************************************************************
c *                                                                    *
c *     CALCULATION OF FACTORIALS                                      *
c *     FAC(N) = 1*2*3*...*(N-1)                                       *
c *     DAC(N) = 1*3*5*...*(2*(N-1)-1)                                 *
c *                                                                    *
c **********************************************************************
      implicit real*8 (a-h,o-z)
      complex*16 cim, ccim
      common/b1/cim,ccim,c(16,16),fac(13),dac(13),ngln(10),
     1alamda,twopi,fpl30m,sws,vol,soo,swsl,fpilam,rmax,gmax,icnvrg,
     2numvr,numvg,nq,nro,nl2,nl2p,klmo,nl,nlm,lmq,nowrt,idim
c
      np=n+1
      fac(1)=1.d0
      fac(2)=1.d0
      do 10 i=3,np
   10 fac(i)=fac(i-1)*(i-1)
      dac(1)=1.d0
      dac(2)=1.d0
      do 11 i=3,np
   11 dac(i)=dac(i-1)*(2*i-3)
      return
      end
