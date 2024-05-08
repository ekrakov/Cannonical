c      program main	
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
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      
      character txt(15)
      CHARACTER*4 POINT      
      INTEGER IDUM1            
      DOUBLE PRECISION DUM1
      DATA POINT/'    '/, IDUM1/0/, DUM1/0.D0/
CCCCCCCCCCCCCCCCCCCCC WORKING AROUND ORIGINAL ISSUE CODE COMMENTED BELOW  CCC    
      
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
     
c      IF LPP.GT.5 ONLY REAL SPACE SUMMATIONS IN <SMTRX>
      
      NL2P=MIN0(5,NL2)
      L=NL2-1
      KLMO=L*(L+1)+L+1

C      CALL FACTOR(4*(NL-1))
C      CALL GAUNT
C      CALL LIMAO(NGLN(1),SOO)

      WRITE(6,5) TXT
      WRITE(6,6) RMAX, GMAX, NL, NQ, NLM, LMQ, ALAMDA, NGLN, NOVCGN,
     1NOWRT, ICNVRG
C      READ PRIMITIVE VECTORS BSX, BSY, BSZ OF REAL SPACE
      DO 20 I=1,3
   20 READ(5,2) BSX(I), BSY(I), BSZ(I)
           
      DO 21 I=1,NQ
   21 READ(5,2) QX(I), QY(I), QZ(I)
      WRITE(6,12) TLAT(LAT), BOA, COA
      DO 30 I=1,3
   30 WRITE(6,11) BSX(I), BSY(I), BSZ(I)
      DO 25 I=1,NQ
   25 WRITE(6,11) QX(I), QY(I), QZ(I)
      
C      GENERATE THE NUMVR AND NUMVG SHORTEST VECTORS
c      OF REAL AND RECIPROCAL SPACE
      
c      CALL VECGEN(NUMR, NUMG, NOVCGN, LAT)
      SWS = (3.D00*VOL/FOURPI)**(1.D00/3.D00)
      SWSL = SWS*ALAMDA
      FPL30M = FOURPI*ALAMDA**3.D00/VOL
      FPLAM = FOURPI*ALAMDA
      WRITE(6,10) VOL, SWS, SOO
      READ(5,7) MODE, NPX, NPY, NPZ
      IF (MODE.EQ.1) GO TO 28
C     GENERATE MESH IN K-SPACE
c      CALL MESH(LAT, NPX, NPY, NPZ, NPT)

      WRITE(6,16) MODE, NPX, NPY, NPZ, NPT
      WRITE(9) TXT, NL, NQ, NLM, LMQ, SWS
      WRITE(19) NPX, NPY, NPZ, LAT, NPT, BOA, COA, DUM1
      
      DO 23 NP=1,NPT
      KKX=KX(NP)
      KKY=KY(NP)
      KKZ=KZ(NP)
c      CALL KTRNSF(LAT, KKX, KKY, KKZ, PKX, PKY, PKZ)
      IF (NOWRT.GT.0) WRITE(6,9) NP, PKX, PKY, PKZ, W(NP), POINT
      W=W(NP)
      WRITE(9) NP, PKX, PKY, PKZ, W, POINT, KKX, KKY, KKZ
      PKXP=PKX*PI
      PKYP=PKY*PI
      PKZP=PKZ*PI
      
c      CALL SMTRX(PKXP, PKYP, PKZP)
      
      J1=-IDIM
      DO 24 J=J1, LMQ
      J1=J1+IDIM
      I=J1+JDIM
      I2=J1+LMQ
      WRITE(9) (SRL(LIN), LIN=I1, I2)
   24 WRITE(9) (SIL(LIN), LIN=I1, I2)
   23 CONTINUE
c sigh this does not work     CLOSE (9,DISP=CRUNCH) tyrying below?
      CLOSE (UNIT=9) 

      STOP

   28 WW = 0.D0
      NPOINT = NPX
      WRITE(6, 16) MODE, NPX, NPY, NPZ, NPOINT
      WRITE(9) TXT, NL, NQ, NLM, LMQ, SWS
      KKY = 0
      KKX = 0
      KKZ = 0
      FX = PI
      FY = PI / BOA
      FZ = PI / COA
      FHX = 0.D0
      IF (LAT .NE. 4) GO TO 27
      FX = 2.D0 / 3.D0 * PI
      FY = PI / DSQRT(3.D0)
      FHX = -0.5D0 * FX

            
   27 DO 29 I=1,NPOINT
      READ(5,8) OKX, OKY, OKZ, POINT
      WRITE(6,9) I, OKX, OKY, OKZ, W, POINT
      WRITE(9) OKX, OKY, OKZ, W, POINT, KKX, KKY, KKZ
      PKXP = OKX*FX*FHX*OKY
      PKYP = OKY*FY
      PKZP = OKZ*FZ
   29 CONTINUE
      
c      CALL SMTRX(PKXP, PKYP, PKZP)
      
      J1 = -IDIM
      DO 26 J=J1, LMQ
      J1 = J1 + IDIM
      I = J1 + JDIM
      I2 = J1 + LMQ
      WRITE(9) (SRL(LIN), LIN=I1, I2)
      WRITE(9) (SIL(LIN), LIN=I1, I2)
   26 CONTINUE
c      CLOSE(9, DISP=CRUNCH)
      CLOSE(unit=9)
      STOP
      END
      BLOCK DATA
C **********************************************************************
C*                                                                     *     
C*    INSERT LM-DATA                                                   *
C*                                                                     *
C **********************************************************************
      COMMON /LDATA/ LL(49), MM(49)
      DATA LL/0,3*1,5*2,7*3,9*4,11*5,13*6/
      DATA MM/0,-1,0,1,-2,-1,0,1,2,-3,-2,-1,0,1,2,3,-4,-3,-2,-1,0,1,2,3,
     14,-5,-4,-3,-2,-1,0,1,2,3,4,5,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6/
      END

      SUBROUTINE VECGEN(NUMR, NUMG, NOVCGN, LAT)
C **********************************************************************
C *                                                                    *
C *                                                                    *
C *                                                                    *
C * GENERATE VECTORS OF DIRECT AND RECIPROCAL SPACE FROM               *
C * BASIC TRANSLATION VECTORS (BSX, BSY, BSZ)                          *
C *                                                                    *
C **********************************************************************
      IMPLICIT REAL*8 (A-H, O-Z)
      COMPLEX*16 CIM, CCIM
      COMMON /B1/CIM, CCIM, CC(16, 16), FAC(13), DAC(13), NGLN(10),
     1ALAMDA,TWOPI,FPL30M,SWS,VOL,SOO, SWSL, FPLAM, RMAX, GMAX, ICNVRG,
     2NUMVR, NUMVG, NQ, NR0,NL2, NL2P, KLM0, NL, NLM, LMQ, NOWRT, IDIM
      COMMON /MSH/DKX,DKY,DKZ,DHX,BOA,COA,KX(2000),KY(2000),KZ(2000),
     1W(2000)
      COMMON/VEC/BSX(3),BSY(3),BSZ(3),BKX(3),BKY(3),BKZ(3),QX(6),QY(6),
     1QZ(6),ASX(500),ASY(500),ASZ(500),AKX(500),AKY(300),AKZ(300),
     2DR(500),DG(300),
     3USX(500),USY(500),USZ(500),UKX(300),UKY(300),UKZ(300)
      DIMENSION CSX(500),CSY(500),CSZ(500),D(500),DD(3),DK(3)
      INTEGER NSR(50),NSG(50)

c     end program
