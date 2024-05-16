      subroutine ylmrk(xx,yy,zz,nlmax)
c **********************************************************************
c *                                                                    *
c *     CALCULATES SPHERICAL HARMONICS                                 *
c *     YLM=CST*DSQRT(4PI)*COMPLEX CONJUGATE(YLM)                      *
c *     CST IS CALCULATED IN <GAUNT>                                   *
c *                                                                    *
c **********************************************************************
      implicit real*8(a-h,o-z)
      complex*16 ylm,dcmplx,cim,ccim
c the code had cc(16,16) but i dont think that is right
      common/b1/cim,ccim,c(16,16),fac(13),dac(13),ngln(10),
     1alamda,twopi,fpl30m,sws,vol,soo,swsl,fpilam,rmax,gmax,icnvrg,
     2numvr,numvg,nq,nro,nl2,nl2p,klmo,nl,nlm,lmq,nowrt,idim
      common/ylmd/ylm(49)
      common/ldata/ll(49),mm(49)
      dimension plm(28),cosmp(7),sinmp(7)
c
c     calculate legendre polynomials by recursion
c
c      PRINT *, 'Before calculation, YLM:', YLM

      PRINT *, 'Before calculation, fac:', fac


      nlp=nlmax
      p=dsqrt(xx*xx+yy*yy)
      x=zz
      y=p
      xa=dabs(x)

      PRINT *, 'nlmax (nlp) =', nlp
      PRINT *, 'p (distance in XY plane) =', p
      PRINT *, 'x (z-coordinate) =', x
      PRINT *, 'y (radial distance in XY plane) =', y
      PRINT *, 'xa (absolute value of x) =', xa
      
      if(xa.gt.1.d-06) go to 10
c      
c     abs(x)=0
c      
      do 11 lp=1,nlp
      l=lp-1
      la=l*(l+1)/2+1
      ta=2.d0**l
      do 11 mp=1,lp
      m=mp-1   
      k=l+m
      if(k-2*(k/2).eq.0) go to 12
      j=la+m
      plm(j)=0.d0
      go to 11
   12 ia=k+1
      ib=k/2+1
      jc=(l-m)/2
      ic=jc+1
      j=la+m
      plm(j)=(((-1)**jc)*fac(ia))/(ta*fac(ib)*fac(ic))
   11 continue
      go to 32
   10 if(xa.lt.0.000000d0) go to 20
c       
c     abs(x)=1

      plm(1)=1.d0
      plm(2)=x
      do 13 lp=3,nlp
      l=lp-1
      j=l*(l+1)/2+1
      l2=2*l-1
      k=j-l
      m=j-l2
   13 plm(j)=(l2*x*plm(k)-(l-1)*plm(m))/l
      do 14 lp=2,nlp
      l=lp-1
      la=l*(l+1)/2
      do 14 mp=2,lp
      j=la+mp
   14 plm(k)=0.d0
      go to 32
c
c     0<abs(x)<1
c
   20 plm(1)=1.d0
      plm(2)=x
      plm(3)=y
      plm(5)=3.d0*y*x
      do 21 lp=3,nlp
      l=lp-1
      j=l*(l+1)/2+1
      l2=2*l-1
      k=j-l
      m=j-l2
   21 plm(j)=(l2*x*plm(k)-(l-1)*plm(m))/l
      do 22 lp=4,nlp
      l=lp-1
      j=l*(l+1)/2+2 
      l2=2*l-1
      k=j-l
      m=j-l2
   22 plm(j) = (l2*x*plm(k) - l*plm(m)) / (l-1) 
      do 23 lp=3, nlp
      l=lp-1
      nm=l*2+1
      la=l*(l+1)/2
      do 23 mp=3, lp
      m=mp-1
      j=la+mp
      k=j-1
      n=k-1
      a=(m-1)*2.d0*x/y
      b=(l-m-1)*(l-m-2)
   23 plm(j) = a*plm(k) - b*plm(n)
      PRINT *, 'LP =', LP, 'MP =', MP, 'J =', J
      PRINT *, 'A =', A, 'PLM(K) =', PLM(K)
      PRINT *, 'B =', B, 'PLM(N) =', PLM(N)
      PRINT *, 'PLM(J) =', PLM(J)
   32 continue
c
c     form spherical harmonics
c
      if (p.gt.1.d-06) go to 34
      cosphi=1.d0
      sinphi=0.d0
      PRINT *, 'cosphi =', cosphi, 'sinphi =', sinphi
      go to 35
   34 cosphi=xx/p
      sinphi=yy/p
      PRINT *, 'cosphi =', cosphi, 'sinphi =', sinphi
   35 cosmp(1)=1.d0
      sinmp(1)=0.d0
      do 33 mp=2, nlp
      cosmp(mp)=cosmp(mp-1)*cosphi-sinmp(mp-1)*sinphi
   33 sinmp(mp)=sinmp(mp-1)*cosphi+cosmp(mp-1)*sinphi
      PRINT *, 'MP =', MP, 'cosmp(MP) =', cosmp(mp)
      PRINT *, 'MP =', MP, 'sinmp(MP) =', sinmp(mp)

      klm=0
      do 36 lp=1,nlp
      l=lp-1
      nm=l*2+1
      do 36 mp=1, nm
      klm=klm+1
      m=mm(klm)
      ma=iabs(m)+1
      mb=ma-1
      lb=l*(l+1)/2+ma
      if (m.le.0) go to 37
      ylm(klm)=plm(lb)*dcmplx(cosmp(ma), -sinmp(ma))
      PRINT *, 'ylm =',ylm(klm)
      go to 36
   37 ylm(klm)=plm(lb)*dcmplx(cosmp(ma),sinmp(ma))
      PRINT *, 'ylm2 =',ylm(klm)
   36 continue
      PRINT *, 'ylm2 =',ylm(klm)
c    PRINT *, 'After calculation, YLM:', YLM


      return
      end



