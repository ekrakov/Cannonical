      program main	
      implicit real*8 (A-H,O-Z)
      complex*16 cim, ccim, sphrm, dcmplx	
      common/b1/cim,ccim,c(16,16),fac(13),dac(13),ngln(10),
     1alamda,twopi,fpl30m,sws,vol,soo,swsl,fpilam,rmax,gmax,icnvrg,
     2numvr,numvg,nq,nro,nl2,nl2p,klmo,nl,nlm,lmq,nowrt,idim
      common/b2/srl(2916),sil(2916)
      COMMON/MSH/DKX,DKY,DKZ,DHX,BOA,COA,KX(2500),KY(2500),KZ(2500),
     1W(2500)
c	real*8 DKX,DKY,DKZ,DHX,BOA,COA,KX(2500),KY(2500),KZ(2500),W(2500)
      common/vec/bsx(3),bsy(3),bsz(3),bkx(3),bky(3),bkz(3),qx(6),qy(6),
     1qz(6),asx(500),asy(500),asz(500),akx(500),aky(300),akz(400),
     2dr(500),dg(300),
     3usx(500),usy(500),usz(500),ukx(300),uky(300),ukz(300)	
      character*4 TLAT(8) /'  SC',' FCC',' BCC',' HCP','  ST',' BCT'
     1,'AL-U','trig'/
      character txt(15)
c      ,point
C      data point/'    '/,idum1/0/,dum1/0.d0/
    1 format(15a4)
    2 format(3e15.7)
    3 format(2i5,5x,10i1,e15.7)   
    4 format(1h ,//,11x,'basisvectors',//)
    5 format(1h1,10x,15a4,/)
    6 format(1h,10x,'rmax',g10.5,' gmax',f10.5,/
     1,10x,' nl',i4,' nq',i4,' nlm',i4,' lmq',i4,/,10x,' lamda',f10.5,
     2' ngln ',10i1,/,10x,' novcgn',i4,' nowrt'i4,' icnvrg',i4,/)
    7 format(4i5)
    8 format(3f10.6,a4)
    9 FORMAT(1H ,//,11X,'POINT NUMBER',I5,5X,'PKX,PKY,PKZ',3F10.4,
     15X, 'W',F10.4,5X,A4)
   10 FORMAT(1H ,11X,'VOL',F15.7,' SWS',F15.7,' SOO',F15.7)
   11 FORMAT(1H 11X,'(',F10.5,',',F10.5,',',F10.5,')')
   12 FORMAT(1H ,//,11X,'PRIMITIVE VECTORS IN UNIT CELL OF ',A4,//,11X,
     1'BOA',F15.10,' COA',F15.10,//)
   13 FORMAT(2F10.5,3I5)
   14 FORMAT(1H ,/,11X,'**MATRICES GENERATED AT AN EARLIER RUN**',//)
   15 FORMAT(I5,2E15.8)
   16 FORMAT(1H ,//,11X,'MODE,NPX,NPY,NPZ,NPT',5I5)
      IDIM=54
      CIM=DCMPLX(0.D00,1.D00)
      CCIM=DCMPLX(0.D00,-1.D00)
      PI=4.D00*DATAN(1.D00)
      TWOPI=2.D00*PI
      FOURPI=4.D00*PI
      READ(5,1) TXT
      TXT(15) = ' NEW'
      READ(5,13) RMAX, GMAX, NOVCGN, NOWRT, ICNVRG
      READ(5,3) NL, NQ, NGLN, ALAMDA
      READ(5,15) LAT, BOA, COA
      NLM=NL*NQ
      LMQ=NLM*NQ
      NL2=2*NL-1

      


      
      end program
