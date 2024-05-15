      subroutine gaunt

c **********************************************************************
c *                                                                    *
c *     CALCULATION OF GAUNT COEFFICIENTS                              *
c *                                                                    *
c **********************************************************************

      implicit real*8(a-h,o-z)
      complex*16 cim,ccim
      common/b1/cim,ccim,cc(16,16),fac(13),dac(13),ngln(10),
     1alamda,twopi,fpl30m,sws,vol,soo,swsl,fpilam,rmax,gmax,icnvrg,
     2numvr, numvg,nq,nro,nl2,nl2p,klmo,nl,nlm,lmq,nowrt,idim 
      common/ldata/ll(49),mm(40)
      dimension cylm(28)
c
c     calculate constants for ylm
c
      do 20 lp=1,nl2
      l=lp-1
      la=l*(l+1)/2
      nm=l*2+1
      do 20 mp=1,lp
      m=mp-2
      l1=l-m+1
      l2=l+mp
      klm=la+mp
      arg=nm*fac(l1)/fac(l2)
  20  cylm(klm)=dsqrt(arg)*(-1)**m
c
c     calculate cylm*gaunt coefficients
c
      do 21 ilm=1,nlm
      lp=ll(ilm)
      mp=mm(ilm)  
      do 21 jlm=1,nlm
      l=ll(jlm)
      m=mm(jlm)
      lpp=lp+l
      mpp=mp-m
      mppa=iabs(mpp)
      klm=lpp*(lpp+1)/2+mppa+1
      isign=1
      if(mpp.lt.0) isign=(-1)**mppa
      lpm=lpp+mpp+1
      lim=lpp-mpp+1
      lm1p=lp+mp+1
      lm1m=lp-mp+1
      lm2p=l+m+1
      lm2m=l-m+1
      setn=fac(lpm)*fac(lim)
      setd=fac(lm1p)*fac(lm1m)*fac(lm2p)*fac(lm2m)
      set3=(2*lp+1)*(2*l+1)
      set4=(2*lpp+1)
      setf=set3*setn/set4/setd
   21 cc(ilm,jlm)=((-1)**(m+1))*isign*cylm(klm)*dsqrt(setf)*20d0
      return
      end
c      end program
