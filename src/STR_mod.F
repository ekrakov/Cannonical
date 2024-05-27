c     program NPmain	AANP
      implicit real*8 (A-H,O-Z)
      complex*16 cim, ccim, sphrm,dcmplx	
      common/b1/cim,ccim,c(16,16),fac(13),dac(13),ngln(10),
     1alamda,twopi,fpl3om,sws,vol,soo,swsl,fpilam,rmax,gmax,icnvrg,
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
      real*8 DUM1
      DATA POINT/'    '/, IDUM1/0/, DUM1/0.D0/
CCCCCCCCCCCCCCCCCCCCC WORKING AROUND ORIGINAL ISSUE CODE COMMENTED BELOW  CCC    
      
c      ,point
C      data point/'    '/,idum1/0/,dum1/0.d0/
    1 format(15a4)
    2 format(3e15.7)
    3 format(2i5,5x,10i1,e15.7)   
    4 format(1h ,//,11x,'basis vectors',//)
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
      print*,'NL:',NL
      PRINT *,'NQ:',NQ
      READ(5,15) LAT, BOA, COA
      NLM=NL*NL
      LMQ=NLM*NQ
      NL2=2*NL-1
     
c      IF LPP.GT.5 ONLY REAL SPACE SUMMATIONS IN <SMTRX>
      
      NL2P=MIN0(5,NL2)
      L=NL2-1
      KLMO=L*(L+1)+L+1

      CALL FACTOR(4*(NL-1))
      CALL GAUNT
      CALL LIMAO(NGLN(1),SOO)

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
      WRITE(6,4)   
      DO 25 I=1,NQ
   25 WRITE(6,11) QX(I), QY(I), QZ(I)
      
C      GENERATE THE NUMVR AND NUMVG SHORTEST VECTORS
c      OF REAL AND RECIPROCAL SPACE


      CALL VECGEN(NUMR, NUMG, NOVCGN, LAT)
     

      SWS = (3.D00*VOL/FOURPI/NQ)**(1.D00/3.D00)
      SWSL = SWS*ALAMDA
      FPL3OM = FOURPI*ALAMDA**3.D00/VOL
      FPILAM = FOURPI*ALAMDA
      WRITE(6,10) VOL, SWS, SOO
      READ(5,7) MODE, NPX, NPY, NPZ
      IF (MODE.EQ.1) GO TO 28
C     GENERATE MESH IN K-SPACE

      CALL MESH(LAT, NPX, NPY, NPZ, NPT)
      WRITE(6,16) MODE, NPX, NPY, NPZ, NPT
      WRITE(9) TXT, NL, NQ, NLM, LMQ, SWS
      WRITE(19) NPX, NPY, NPZ, LAT, NPT, BOA, COA, DUM1

      DO 23 NP=1,NPT
      KKX=KX(NP)
      KKY=KY(NP)
      KKZ=KZ(NP)
      CALL KTRNSF(LAT, KKX, KKY, KKZ, PKX, PKY, PKZ)
CCCCCCCCCCCCCCCCCCCCCCCCTHIS IS TRUE GT =1 IN THIS PROGRAM CCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCC GT =0 MEANS       
      IF (NOWRT.GT.0) WRITE(6,9) NP, PKX, PKY, PKZ, W(NP), POINT
      WW=W(NP)
      WRITE(9) NP, PKX, PKY, PKZ, WW, POINT, KKX, KKY, KKZ
      PKXP=PKX*PI
      PKYP=PKY*PI
      PKZP=PKZ*PI
      CALL SMTRX(PKXP, PKYP, PKZP)
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
      FY = PI/BOA
      FZ = PI/COA
      FHX = 0.D0
      IF (LAT .NE. 4) GO TO 27
      FX = 2.D0/3.D0 * PI
      FY = PI/DSQRT(3.D0)
      FHX = -0.5D0*FX

            
   27 DO 29 I=1,NPOINT
      READ(5,8) OKX, OKY, OKZ, POINT
      WRITE(6,9) I, OKX, OKY, OKZ, WW, POINT
      WRITE(9) I OKX, OKY, OKZ, W, POINT, KKX, KKY, KKZ
      PKXP = OKX*FX+FHX*OKY
      PKYP = OKY*FY
      PKZP = OKZ*FZ
   29 CONTINUE
      CALL SMTRX(PKXP, PKYP, PKZP)
      
      J1 = -IDIM
      DO 26 J=J1, LMQ
      J1 = J1 + IDIM
      I1=J1+J
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
      COMMON/LDATA/LL(49),MM(49)
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
     1ALAMDA,TWOPI,FPL3OM,SWS,VOL,SOO, SWSL, FPILAM, RMAX, GMAX, ICNVRG,
     2NUMVR, NUMVG, NQ, NRO,NL2, NL2P, KLMO, NL, NLM, LMQ, NOWRT, IDIM
      COMMON /MSH/DKX,DKY,DKZ,DHX,BOA,COA,KX(2000),KY(2000),KZ(2000),
     1W(2000)
      COMMON/VEC/BSX(3),BSY(3),BSZ(3),BKX(3),BKY(3),BKZ(3),QX(6),QY(6),
     1QZ(6),ASX(500),ASY(500),ASZ(500),AKX(500),AKY(300),AKZ(300),
     2DR(500),DG(300),
     3USX(500),USY(500),USZ(500),UKX(300),UKY(300),UKZ(300)
      DIMENSION CSX(500),CSY(500),CSZ(500),D(500),DD(3),DK(3)
      INTEGER NSR(50),NSG(50)
      integer NR    
    1 FORMAT(1H,//, 15X,'SHELL NUMBER',I5,' WITH',I5,' POINTS',//)
    2 FORMAT(1H 10X,I5,4F10.6)
    3 FORMAT(1H ,10X,'(',F10.5,',',F10.5,',',F10.5,')',/)
    4 FORMAT(1H1 ,14X,'RESULT FROM VECGEN FOR REAL SPACE VECTORS',//,
     114X,'NO',5X,'SX',8X,'SY',8X,'SZ',8X,'D',/)
    5 FORMAT(1H1 ,10X,'RESULT FROM VECGEN FOR RECIPROCAL SPACE VECTORS',
     1 //,'NO',5X,'KX',8X,'KY',8X,'KZ',8X,'D',/)
    6 FORMAT(1H,//,11X,'PRIMITIVE VECTORS OF RECIPROCAL SPACE',/)
    7 FORMAT(1H,//,11X,'R1',F10.4,'RA',F10.4,/,11X,'G1',F10.4,
     1' GA',F10.4,/) 
    8 FORMAT(1H ,'//',11X,'VALUES GENERATED AT PREVIOUS RUN OF VECGEN , (R
     1IMAX,NUMVR,GMAX,NUMVG) =',//,11X,F10.5,I5,5X,F10.5,I5,//)
    9 FORMAT(1H ,'*** NR EXCEEDS DIMENSIONS OF ASX, ASY, ASZ',I5,3F10.4)
   10 FORMAT(1H ,'*** NG EXCEEDS DIMENSIONS OF AKX, AKY, AKZ',I5,3F10.4)
   11 FORMAT(1H ,10X,'NUMR',I4,' NUMG',I4,//)


     
      IF (NOVCGN.EQ.0) GO TO 32
      WRITE(6,6)
      
C      PRIMITIVE VECTORS (BKX, BKY, BKZ) OF RECIPROCAL SPACE
      
      DO 20 I=1,3
      I1=1+MOD(I,3)
      I2=1+MOD(I1,3)
      CALL CROSS(BSX(I1),BSY(I1),BSZ(I1),BSX(I2),BSY(I2),BSZ(I2),
     1BKX(I),BKY(I),BKZ(I))
   20 CONTINUE
      VOL=DABS(BSX(1)*BKX(1)+BSY(1)*BKY(1)+BSZ(1)*BKZ(1))
      
      DO 21 I=1,3
      BKX(I)=BKX(I)/VOL*TWOPI
      BKY(I)=BKY(I)/VOL*TWOPI
      BKZ(I)=BKZ(I)/VOL*TWOPI
   21 WRITE(6,3) BKX(I),BKY(I),BKZ(I)

C
C CALCULATE RADII RA, GA OF SPHERES HOLDING ALL VECTORS USED IN LAT-
C TICE SUMS. R1 IS LONGEST BASIS VECTOR. G1 IS LONGEST VECTOR IN
C BRILLOUIN ZONE. MUST BE RECONSIDERED IN ANY NEW APPLICATIONS
C
      R1 = 1.E-06
      DO 35 IQ=1,NQ
      PQX = QX(IQ)
      PQY = QY(IQ)
      PQZ = QZ(IQ)
      DO 35 JQ=1,NQ
      X = PQX - QX(JQ)
      Y = PQY - QY(JQ)
      Z = PQZ - QZ(JQ)
      DQ = DSQRT(X*X + Y*Y + Z*Z)
      IF (DQ .LT. R1) GO TO 35
      R1 = DQ
   35 CONTINUE
      R1 = R1 * 1.001
      RA = RMAX + R1
      B1 = BKX(1) + BKX(2) + BKX(3)
      B2 = BKY(1) + BKY(2) + BKY(3)
      B3 = BKZ(1) + BKZ(2) + BKZ(3)
      G1 = DSQRT(B1*B1 + B2*B2 + B3*B3) / 2.D00
      GA = GMAX + G1
      IF (LAT.EQ.3) GA = GMAX + 6.2832
      IF (LAT.EQ.8) GA = GMAX + 6.52
      IF (LAT.EQ.6) GA = GMAX + DSQRT(0.5 + (1./COA)**2) * 6.29
      WRITE(6,7) R1, RA, G1, GA


      DO 36 I=1,3
      DD(I)=DSQRT(BSX(I)**2+BSY(I)**2+BSZ(I)**2)
   36 DK(I)=DSQRT(BKX(I)**2+BKY(I)**2+BKZ(I)**2)
      DDM=DMAX1(DD(1),DD(2),DD(3))
      DKM=DMAX1(DK(1),DK(2),DK(3))
      DDM=TWOPI/DDM
      DKM=TWOPI/DKM
      NUMR=2*(INT(RA/DKM)+1)+1
      NUMG=2*(INT(GA/DDM)+1)+1
      NUMRH=NUMR/2+1
      NUMGH=NUMG/2+1
      WRITE(6,11) NUMR, NUMG
C
C REAL SPACE
C
      WRITE(6,4)
      NR=0
      NRO=0
      DO 22 L=1,NUMR   
      A=L-NUMRH
      DO 22 M=1,NUMR
      B=M-NUMRH
      DO 22 N=1,NUMR
      C=N-NUMRH
      SX=A*BSX(1)+B*BSX(2)+C*BSX(3)
      SY=A*BSY(1)+B*BSY(2)+C*BSY(3)
      SZ=A*BSZ(1)+B*BSZ(2)+C*BSZ(3)
      DX=DSQRT(SX*SX+SY*SY+SZ*SZ)
      IF (DX.GT.RA) GO TO 22
      IF (DX.LE.RMAX) NRO=NRO+1
      NR=NR+1
      IF (NR .GT. 500) GO TO 33
      D(NR)=DX
      CSX(NR)=SX
      CSY(NR)=SY
      CSZ(NR)=SZ
   22 CONTINUE
C
C SORT VECTORS IN ORDER OF INCREASING D
C

   
      DA=1.D-06
      NSH=0
      NSHL=-1
      DO 23 K=1,NR
      AMIN=1000.
      DO 24 N=1,NR
      IF (D(N)-AMIN) 25, 24, 24
   25 AMIN=D(N)
      N1=N
   24 CONTINUE
      NSHL=NSHL+1
      ASX(K)=CSX(N1)
      ASY(K)=CSY(N1)
      ASZ(K)=CSZ(N1)
      DB=D(N1)
      DR(K)=DB
      IF (DB .GT. DA+1.D-06) GO TO 26
      WRITE(6,2)K, ASX(K), ASY(K), ASZ(K), DB
      GO TO 23
   26 NSH=NSH+1
      WRITE(6,1) NSH, NSHL
      NSR(NSH)=NSHL
      WRITE(6,2)K, ASX(K), ASY(K), ASZ(K), DB
      NSHL=0
      DA=DB
   23 D(N1)=1000.
      NSH=NSH+1
      NSHL=NSHL+1
      NSR(NSH)=NSHL
      WRITE(6,1) NSH, NSHL
      NUMVR=NR
      NSHLR=NSH

c      
C RECIPROCAL SPACE
C

      WRITE(6,5)
      NG=0
      DO 27 L=1,NUMG
      A=L-NUMGH
      DO 27 M=1,NUMG
      B=M-NUMGH
      DO 27 N=1,NUMG
      C=N-NUMGH
      GX=A*BKX(1)+B*BKX(2)+C*BKX(3)
      GY=A*BKY(1)+B*BKY(2)+C*BKY(3)
      GZ=A*BKZ(1)+B*BKZ(2)+C*BKZ(3)
      DX=DSQRT(GX*GX+GY*GY+GZ*GZ)
      IF (DX .GT. GA) GO TO 27
      NG=NG+1
      IF (NG .GT. 300) GO TO 34
      D(NG)=DX
      CSX(NG)=GX
      CSY(NG)=GY
      CSZ(NG)=GZ
   27 CONTINUE
C

C SORT VECTORS IN ORDER OF INCREASING D
C

      DA=1.E-06
      NSH=0
      NSHL=-1
      DO 28 K=1,NG
      AMIN=1000.
      DO 29 N=1,NG
      IF (D(N)-AMIN) 30, 29, 29 
   30 AMIN=D(N)
      N1=N
   29 CONTINUE


      NSHL=NSHL+1
      AKX(K)=CSX(N1)
      AKY(K)=CSY(N1)
      AKZ(K)=CSZ(N1)
      DB=D(N1)
      IF (DB.GT.DA+1.D-07) GO TO 31
      WRITE(6,2) K, AKX(K), AKY(K), AKZ(K), DB
      GO TO 28
   31 NSH=NSH+1
      WRITE(6,1) NSH, NSHL
      WRITE(6,2) K, AKX(K), AKY(K), AKZ(K), DB
      NSG(NSH)=NSHL
      NSHL=0
      DA=DB
   28 D(N1)=1000.
      NSH=NSH+1
      NSHL=NSHL+1
      NSG(NSH)=NSHL
      WRITE(6,1) NSH, NSHL
      NUMVG=NG
      NSHLG=NSH
C
C       STORE VECTORS ON FILE 1
C      

      WRITE(1) BKX,BKY,BKZ,VOL,RMAX,GMAX,NR,NG,NRO
      WRITE(1)(ASX(N),N=1,NR)
      WRITE(1)(ASY(N),N=1,NR)
      WRITE(1)(ASZ(N),N=1,NR)
      WRITE(1)(DR(N), N=1,NR)
      WRITE(1)(AKX(N),N=1,NG)
      WRITE(1)(AKY(N),N=1,NG)
      WRITE(1)(AKZ(N),N=1,NG)
      RETURN
   32 continue
c     READ VECTORS GENERATED IN AN EARLIER RUN
      READ(1) BKX,BKY,BKZ,VOL,RMAX,GMAX,NUMVR,NUMVG,NRO
      READ(1)(ASX(N),N=1,NUMVR)
      READ(1)(ASY(N),N=1,NUMVR)
      READ(1)(ASZ(N),N=1,NUMVR)
      READ(1)(DR(N),N=1,NUMVR)
      READ(1)(AKX(N),N=1,NUMVG)
      READ(1)(AKY(N),N=1,NUMVG)
      READ(1)(AKZ(N),N=1,NUMVG)
      WRITE(6,8) RMAX,NUMVR,GMAX,NUMVG
      RETURN
   33 WRITE(6,9) NR,DX,RA,RMAX
      STOP
   34 WRITE(6,10)NG,DX,GA,GMAX
      STOP
      END


      SUBROUTINE SMTRX(PKXP,PKYP,PKZP)
C**********************************************************************
C     
C CALCULATION OF THE STRUCTURE CONSTANT MATRIX
C
C **********************************************************************
      IMPLICIT REAL*8 (A-H, O-Z)
      COMPLEX*16 SUMR, SUMG, SR, CIM, CCIM,CIL,DCMPLX,YLM,CDEXP,CXPRDK,
     1SRC
      REAL*8 DREAL, DIMAG

      COMMON /B1/CIM, CCIM, C(256), FAC(13), DAC(13), NGLN(10),
     1ALAMDA,TWOPI,FPL3OM,SWS,VOL,SOO, SWSL, FPILAM, RMAX, GMAX, ICNVRG,
     2NUMVR, NUMVG, NQ, NRO,NL2, NL2P, KLMO, NL, NLM, LMQ, NOWRT, IDIM


      
C      Common /B1/ CIM, CCIM, C(256), FAC(13), DAC(13), NGLN(10),
C     1ALAMDA, TWOPI, FPL3OM, SWS,VOL,SOO,SWSL,FPILAM RMAX,GMAX,ICNVRG,
C     2NUMVR, NUMVG, NQ, NRO, NL2, NL2P, KLMO, NL, NLM, LMQ, NOWRT, IDIM
      COMMON /B2/ SRKLP(2916), SIKLP(2916)
      COMMON /LDATA/LL(49), MM(49)
      COMMON /VEC/BSX(3),BSY(3),BSZ(3),BKX(3),BKY(3),BKZ(3),QX(6),QY(6),
     1QZ(6), ASX(500), ASY(500), ASZ(500), AKX(500), AKY(300), AKZ(300),
     2DR(500),DG(300),
     3USX(500),USY(500),USZ(500),UKX(300),UKY(300),UKZ(300)
      COMMON/YLMD/YLM(49)
      DIMENSION GLNS(500,7),FBLN(500,7),BDOTQ(500),DRQ(500),
     1SR(49), CXPRDK(500), SUMR(49), SUMG(49)
  100 FORMAT(1H ,9(1PE13.5))
  101 FORMAT(1H ,///)
  102 FORMAT(1HO,10X,'QQ DIAGONAL BLOCK',/)
  103 FORMAT(1HO,10X,1PE11.3)
  104 FORMAT(1HO,10X,'QP Q OFF DIAGONAL BLOCK ', 'QP, Q =',2I5,/)
  105 FORMAT(1HO,25X,'CONVERGENCE TEST')
  106 FORMAT(1HO,'RECIPROCAL SPACE SUMS',/)
  107 FORMAT(1HO,'REAL SPACE SUMS',/)

C     QQ DIAGONAL BLOCKS
      

      PKX=PKXP
      PKY=PKYP
      PKZ=PKZP
      CXPRDK(1)=DCMPLX(1.D00,0.D00)

      DO 22 I=2,NRO
      DRI=DR(I)

      USX(I)=ASX(I)/DRI
      USY(I)=ASY(I)/DRI
      USZ(I)=ASZ(I)/DRI
      PDR=PKX*ASX(I)+PKY*ASY(I)+PKZ*ASZ(I)
      CXPRDK(I)=CDEXP(PDR*CIM)
      ALFA=ALAMDA*DRI

C    CONVERGENCE FUNCTION, REAL SPACE


      DO 25 LP=1, NL2P 
C      DO 25 LP=1, 3            

      L=LP-1
      NG=NGLN(LP)
   25 GLNS(I,LP)=GLN(ALFA,L,NG)

   22 CONTINUE
      DO 35 I=l,NUMVG
      X=AKX(I)+PKX
      Y=AKY(I)+PKY
      Z=AKZ(I)+PKZ
      D=DSQRT(X*X+Y*Y+Z*Z)

      IF(D.GT.GMAX) GO TO 35
      BETA=D/ALAMDA
      UKX(I)=X/D
      UKY(I)=Y/D
      UKZ(I)=Z/D

c     CONVERGENCE FUNCTION, RECIP. SPACE

      DO 33 LP=1,NL2P
      L=LP-1

      NG=NGLN(LP)
   33 FBLN(I,LP)=BETA**(L-2)/(1.D00+BETA*BETA)**(NG+1)
   35 DG(I)=D

c     EVALUATION OF LATTICE SUMS FOR QPP .EQ. 0

      DO 44 KLM=1,KLMO
      SUMR(KLM)=CMPLX(0.D00,0.D00)
   44 SUMG(KLM)=CMPLX(0.D00,0.D00)

c     real space

      DO 49 I=2,NRO
      CALL YLMRK(USX(I),USY(I),USZ(I),NL2)
      CALL YLMRK(PI/2,PI/2,PI/2,5)                                                        

      PRINT *, 'USX(', I, ') = ', USX(I)
      PRINT *, 'USY(', I, ') = ', USY(I)
      PRINT *, 'USZ(', I, ') = ', USZ(I)
      PRINT *, 'NL2 = ', NL2
      PRINT *, 'YLM = ',YLM
      KLM=0
      DO 26 LP=1,NL2P
      L=LP-1
      NM=L*2+1
      DO 26 M=1,NM
C     PRINT *,"MADE IT TO 553 NL2",NL2
      KLM=KLM+1
   26 SUMR(KLM)=SUMR(KLM)+GLNS(I,LP)*YLM(KLM)*CXPRDK(I)   
C 26   YLM(1)=YLM(1)
C    26 KLM=1  


      IF(NL2.LE.5) GO TO 49
      
c     IF L IS LARGER THAN 4, REAL SPACE SUMMATIONS ONLY

      DO 50 LP=6,NL2
      NM=LP*2-1
      SORO=(SWS/DR(I))**LP
      DO 50 M=1,NM
      KLM=KLM+1
   50 SUMR(KLM)=SUMR(KLM)+SORO*YLM(KLM)*CXPRDK(I)
   49 CONTINUE

c     recip space
      
      DO 34 I=1,NUMVG
      IF(DG(I).GT.GMAX) GO TO 34
      CALL YLMRK(UKX(I),UKY(I),UKZ(I),NL2P)
      KLM=0
      DO 28 LP=1,NL2P
      L=LP-1
      NM=L*2+1
      DO 28 M=1,NP
      KLM=KLM+1
   28 SUMG(KLM)=SUMG(KLM)+FBLN(I,LP)*YLM(KLM)
   34 CONTINUE 

C     SET UP STRUCTURE CONSTANT MATRIX FROM LATTICE SUMS

      CIL=CIM
      KLM=0
      DO 51 LP=1,NL2
      L=LP-1
      NM=L*2+1
      CIL=CIL*CCIM
      IF(LP.GT.5) GO TO 52
      AR=SWSL**LP
      AG=FPL3OM*AR/DAC(LP)
      DO 45 M=1,NM
      KLM=KLM+1
   45 SR(KLM)=AR*SUMR(KLM)*CIL+AG*SUMG(KLM)
      GO TO 51
   52 DO 53 M=1,NM
      KLM=KLM+1
   53 SR(KLM)=SUMR(KLM)*CIL
   51 CONTINUE  
      SR(1)=SR(1)-SOO*SWSL
C     INSERT INTO LOWER TRIANGLE

      J1=-IDIM
      J2=-16
      DO 27 JLM=1, NLM
      J2=J2+16
      J1=J1+IDIM
      L=LL(JLM)
      M=MM(JLM)
      DO 27 ILM=JLM, NLM
      LIN=J1+ILM
      LINC=J2+ILM
      LPP=L+LL(ILM)
      MPP=MM(ILM)-M
      KLM=LPP*LPP+LPP+MPP+1
      SRC=SR(KLM)*C(LINC)
      SRKLP(LIN)=REAL(SRC)
      SIKLP(LIN)=AIMAG(SRC)
   27 CONTINUE

      IF(ICNVRG.EQ.0) GO TO 58

c     PRINT INDIVIDUAL CONTRIBUTIONS TO THE LATTICE SUMS. USED
c     FOR TEST OF THE CHOICE OF SPLITTING PARAMETERS <LAMDA)
c     AND <NGLN>
      
C      WRITE(6, 105)
      KLM=0
      DO 59 LP=1, NL2P
      NM=LP*2-1
      AR=SWSL**LP
      AG=FPL3OM*AR/DAC(LP)
      DO 59 M=1, NM
      KLM=KLM+1
      SUMR(KLM)=SUMR(KLM)*AR
      SUMG(KLM)=SUMG(KLM)*AG
   59 CONTINUE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCC THESE SHOULD NOT BE COMMENTED CCCCCCCCCCCC      
C      WRITE(6, 107)
C      WRITE(6,103) (SUMR(K),K=1,KLM)
C      WRITE(6, 106)
C      WRITE(6, 103) (SUMG(K), K=1, KLM)
   58 IF (NOWRT.EQ.0) GO TO 30


      
C PRINT STRUCTURE CONSTANT MATRIX
C     WRITE(6, 102)
      I1=-IDIM
      DO 42 ILM=1, NLM
      I1=I1+IDIM
      I2=I1+ILM
   42 WRITE(6, 100) (SRKLP(LIN), LIN=ILM, I2, IDIM)
C      WRITE(6, 101)
      I1=-IDIM
      DO 43 ILM=1, NLM
      I1=I1+IDIM
      I2=I1+ILM
   43 WRITE(6, 100) (SIKLP(LIN), LIN=ILM, I2, IDIM)
C      WRITE(6, 101)
      
   30 CONTINUE
      IF(NQ.EQ.1) RETURN 
C
C     REPEAT FIRST DIAGONAL BLOCK
c
      JQQ=0
      DO 24 JQ=2, NQ
      JQQ=JQQ+(1+IDIM)*NLM
      J1=IDIM
      DO 24 JLM=1, NLM
      J1=J1+IDIM
      DO 24 ILM=JLM, NLM
      LIN=J1+ILM

      LNIQ=JQQ+LIN
      SRKLP(LINQ)=SRKLP(LIN)
C ERROR MUST HAVE POINTER VALRIABLE 
   24 SIKLP(LINQ)=SIKLP(LIN)
      NROP=NRO+1
C     PG 141
      
      DO 29 I=1,NROP,NUMVR
      PRD=PKX*ASX(I)+PKY*ASY(I)+PKZ*ASZ(I)
  29  CXPRDK(I)=CDEXP(CIM*PRD)
      
c      QQP OFF DIAGONAL BLOCKS (PG 141)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC CHECK THIS CHECK THIS 5/22/24
      NQM=NQ-1
      JQQ=-NLM*IDIM
      DO 20 JQ=1,NQM
      JQP=JQ+1
      JQQ=JQQ+NLM*IDIM
      JQP=JQ+1
      JQ0Q=JQO+NLM*IDIM
      QPX=QX(JQ)
      QPY=QY(JQ)
      QPZ=QZ(JQ)

      DO 20 IQ=JQP,NQ
      IQQ=(IQ-1)*NLM*JQQ
      QPPX=QX(IQ)-QPX
      QPPY=QY(IQ)-QPY
      QPPZ=QZ(IQ)-QPZ
      DO 36 I=1,NUMVR
      X=ASX(I)-QPPX
      Y=ASY(I)-QPPY
      Z=ASZ(I)-QPPZ
      D=DSQRT(X*X+Y*Y+Z*Z)
      DRQ(I)=D
      IF (D.GT.RMAX) GO TO 36
      ALFA=ALAMDA*D
      USX(I)=X/D
      USY(I)=Y/D
      USZ(I)=Z/D
 
C
C     CONVERGENCE FUNCTION, REAL SPACE
C

      DO 37 LP=1,NL2P
      L=LP-1
      NG=NGLN(LP)
 37   GLNS(I,LP)=GLN(ALFA,L,NG)
  36  CONTINUE
      DO 38 I=1,NUMVG
      IF (D.GT.GMAX) GO TO 38
      BDOTQ(I)=(UKX(I)*QPPX+UKY(I)*QPPY+UKZ(I)*QPPZ)*D
      BETA=D/ALAMDA
C
C     CONVERGENCE FUNCTION, RECIP. SPACE
C
      DO 39 LP=1,NL2P
      L=LP-1
      NG=NGLN(LP)
  39  FBLN(I,LP)=BETA**((L-2)/(1.D00+BETA*BETA))*(NG+1)
  38  CONTINUE

ccccccccccccccccccccccccccccccccccccccccccccccccc
C
C     EVALUATION OF LATTICE SUMS FOR QPP .NE. 0
C
      DO 46 KLM=1,KLM0
      SUMR(KLM)=CMPLX(0.D00,0.D00)
   46 SUMG(KLM)=CMPLX(0.D00,0.D00)
   
C
C     REAL SPACE
C
      DO 40 I=1,NUMVR
      IF(DRQ(I).GT.RMAX) GO TO 40
      CALL YLMRK(USX(I),USY(I),USZ(I),NL2)
      KLM=0
      DO 21 LP=1,NL2P
      L=LP-1
      NM=L*2+1
      DO 21 M=1,NM
      KLM=KLM+1
   21 SUMR(KLM)=SUMR(KLM)+GLNS(I,LP)*YLM(KLM)*CXPRDK(I)
      IF(NL2. LE.5) GO TO 40
C
C     IF L IS LARGER THAN 4, REAL SPACE SUMMATIONS ONLY
C
      DO 54 LP=6,NL2
      SORO=(SWS/DRQ(I))**LP
      NM=LP*2-1
      DO 54 M=1,NM
      KLM=KLM+1
   54 SUMR(KLM)=SUMR(KLM)+SORO*YLM(KLM)*CXPRDK(I)
   40 CONTINUE
C
C     RECIP. SPACE
C
      DO 41 I=1,NUMVG
      IF(DG(I).GT.GMAX) GO TO 41
      CALL YLMRK(UKX(I),UKY(I),UKZ(I),NL2P)
      KLM=0
      DO 47 LP=1,NL2P
      L=LP-1
      NM=L*2+1
      DO 47 M=1,NM
      KLM=KLM+1
   47 SUMG(KLM)=SUMG(KLM)+FBLN(I,LP)*YLM(KLM)*CDEXP(CIM*BDOTQ(I))
   41 CONTINUE
C

C
C     SET UP STRUCTURE CONSTANT MATRIX FROM THE LATTICE SUMS
C
      CIL=CIM
      KLM=0
      DO 56 LP=1,NL2
      L=LP-1
      NM=L*2+1
      CLI=CIL+CCIM
      IF(LP.GT.5) GO TO 55
      AR=SWSL*LP
      AG=PL30M*AR/DAC(LP)
      DO 48 M=1,NM
      KLM=KLM+1
   48 SR(KLM)=AR*SUMR(KLM)*CIL+AG*SUMG(KLM)
      GO TO 56
   55 DO 57 M=1,NM
      KLM=KLM+1
   57 SR(KLM)=SUMR(KLM)*CIL
   56 CONTINUE
C
C     INSERT INTO LOWER TRIANGLE
C

C
      J1E=IDIM
      J2=16
      DO 23 JLM=1,NLM
      J1=J1+IDIM
      J2=J2+16
      L=LL(JLM)
      DO 23 ILM=1,NLM
      LIN=J1+ILM
      LINC=J2+ILM
      LINQ=JQ+LIN
      LPP=L+LL(ILM)
      MPP=MM(ILM)-M
      KLM=LPP*LPP+LPP+MPP+1
      SRC=SR(KLM)+C(LINC)
      SRKLP(LINQ)=REAL(SRC)
   23 SIKLP(LINQ)=AIMAG(SRC)
 
C
      IF(ICNVG.EQ.0) GO TO 60
C
C     TEST OF SPLITTING PARAMETERS
C

      KLM=0
      DO 61 LP=1,NL2P
      AR=SWSL*LP
      NM=LP*2+1
      AG=PL30M*AR/DAC(LP)
      DO 61 M=1,NM
      KLM=KLM+1
      SUMR(KLM)=SUMR(KLM)*AR
   61 SUMG(KLM)=SUMG(KLM)*AG
      WRITE(6,107)
      WRITE(6,103) (SUMR(K),K=1,KLM)
      WRITE(6,106)
      WRITE(6,103) (SUMG(K),K=1,KLM)
   60 IF(NOWRT.EQ.0) GO TO 20
C
C     PRINT STRUCTURE CONSTANT MATRIX
C     PG 143

      WRITE(6,104) IQ,JQ
      DO 31 ILM=1,NLM
      J1=IQQ+ILM
      J2=J1+(NLM-1)*IDIM
   31 WRITE(6,100) (SRKLP(LINQ),LINQ=J1,J2,IDIM)
      WRITE(6,101)
      DO 32 ILM=1,NLM
      J1=IQQ+ILM
      J2=J1+(NLM-1)*IDIM
   32 WRITE(6,100) (SIKLP(LINQ),LINQ=J1,J2,IDIM)
   20 CONTINUE
      RETURN
      END

      SUBROUTINE CROSS(AX,AY,AZ,BX,BY,BZ,CX,CY,CZ)

C
C     *****************************************************************
C     * CROSS PRODUCT (CX,CY,CZ) = (AX,AY,AZ) * (BX,BY,BZ)           *
C     *****************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
      CX=AY*BZ-BY*AZ
      CY=BX*AZ-AX*BZ
      CZ=AX*BY-BX*AY
      RETURN
      END

      FUNCTION GLN(ALFA,LN,NN)
C
C     *****************************************************************
C     * CALCULATION OF CONVERGENCE FUNCTIONS                          *
C     *****************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 CIM,CCIM
      COMMON /B1/CIM,CCIM,C(16,16),FAC(13),DAC(13),NGLN(10),
     1ALAMDA,TWOPI,FPL3OM,SWS,VOL,SOO,SWSL,FPILAM,RMAX,GMAX,ICNVRG,
     2NUMVR,NUMVG,NQ,NRO,NL2,NL2P,KLMO,NL,NLM,LMQ,NOWRT,IDIM
C
      ALFASF=ALFA*ALFA/4.D00
      EXA=DEXP(-ALFA)
      HOO=EXA/ALFA
      HO1=0.5D00*EXA
      H1=HOO
      H2=HO1

C
C     FOR L = 0
C
      GLN=H1+H2
      DO 20 NP=2,NN
      N=NP-1
      H3=(ALFASF/N*H1-(0.5D00-N)*H2)/NP
      GLN=GLN+H3
      H1=H2
   20 H2=H3
      IF(LN.EQ.0) RETURN
C
C     FOR L = 1
C
      H11=0.5D00*EXA
      H1O=(1.D00+ALFA)/ALFA/ALFA*EXA
      H1=H1O
      H2=H11
      GLN=H1+H2
      DO 21 NP=2,NN
      N=NP-1
      H3=(ALFASF/N*H1-(1.5D00-N)*H2)/NP
      GLN=GLN+H3
      H1=H2
   21 H2=H3
      IF(LN.EQ.1) RETURN
C
C     FOR L .GT. 1
C
      L=1
   23 L=L+1

      LP=L+1
      LM=L-1
      H21=H1O+(2*LM-1)*H11/ALFA
      H2O=H1O+(2*LM+1)*H1O/ALFA
      H1=H2O
      H2=H21
      GLN=H1+H2
      DO 22 NP=2,NN
      N=NP-1
      XLN=L-N
      XLN=XLN+0.500
      H3=(ALFASF/N*H1-XLMN*H2)/NP
      GLN=GLN+H3
      H1=H2
  22  H2=H3

      GLN=GLN/DAC(LP)

      IF(L.EQ.LN) RETURN
      HOO=H1O
      H1O=H2O
      HO1=H11
      H11=H21
      GO TO 23
      END

c      program main
c      implicit real*8(a-h,o-z)
c      call gaunt
c      call factor
c      end program main

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

      SUBROUTINE LIMAO(NN,FALFB)

C     *****************************************************************
C     * CALCULATES ALPHA**-1-GLN(ALPHA) FOR L=0 AND ALPHA=0          *
C     *****************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
      NUM=0
      ALFA=1.E-03
      FALFA=10.
   10 ALFAS=ALFA*ALFA/4.D00
      EXA=DEXP(-ALFA)
      H1=EXA/ALFA
      H2=0.5D00*EXA
      GON=H1+H2
      DO 11 NP=2,NN
      N=NP-1
      H3=((ALFAS/N)*H1-(0.5-N)*H2)/NP                                                                                                                                                                     
      GON=GON+H3
      H1=H2
   11 H2=H3
      FALFB=(1.-ALFA*GON)/ALFA
      DIF=FALFA-FALFB
      NUM=NUM+1
      IF(DIF.LT.1.E-05) RETURN
      IF(NUM.GT.10) WRITE(6,1)
    1 FORMAT(1H1,10X,'**TOO MANY ITERATIONS IN LIMA0**',/)
      FALFA=FALFB
      ALFA=ALFA-0.1E-04
      GO TO 10
      END
      

      SUBROUTINE YLMRK(XX,YY,ZZ,NLMAX)
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
      common/b1/cim,ccim,cc(16,16),fac(13),dac(13),ngln(10),
     1alamda,twopi,fpl3om,sws,vol,soo,swsl,fpilam,rmax,gmax,icnvrg,
     2numvr,numvg,nq,nro,nl2,nl2p,klmo,nl,nlm,lmq,nowrt,idim
      common/ylmd/ylm(49)
      common/ldata/ll(49),mm(49)
      dimension plm(28),cosmp(7),sinmp(7)
c
c     calculate legendre polynomials by recursion
c
      nlp=nlmax
      p=dsqrt(xx*xx+yy*yy)
      x=zz
      y=p
      xa=dabs(x)
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
   10 if(xa.lt.0.999999d0) go to 20
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
   14 plm(j)=0.d0
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
      b=(l+m-1)*(l-m+2)
   23 plm(j) = a*plm(k) - b*plm(n)
   32 continue
c
c     form spherical harmonics
c
      if (p.gt.1.d-06) go to 34
      cosphi=1.d0
      sinphi=0.d0
      go to 35
   34 cosphi=xx/p
      sinphi=yy/p
   35 cosmp(1)=1.d0
      sinmp(1)=0.d0

      do 33 mp=2, nlp
      cosmp(mp)=cosmp(mp-1)*cosphi-sinmp(mp-1)*sinphi
   33 sinmp(mp)=sinmp(mp-1)*cosphi+cosmp(mp-1)*sinphi
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
      PRINT *, 'KLM',KLM
      PRINT *, 'YLM', YLM(KLM)
      go to 36
   37 ylm(klm)=plm(lb)*dcmplx(cosmp(ma),sinmp(ma))
   36 continue
      return
      end



      subroutine gaunt

c **********************************************************************
c *                                                                    *
c *     CALCULATION OF GAUNT COEFFICIENTS                              *
c *                                                                    *
c **********************************************************************

      implicit real*8(a-h,o-z)
      complex*16 cim,ccim
      common/b1/cim,ccim,cc(16,16),fac(13),dac(13),ngln(10),
     1alamda,twopi,fpl3Om,sws,vol,soo,swsl,fpilam,rmax,gmax,icnvrg,
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
      m=mp-1
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
   21 cc(ilm,jlm)=((-1)**(m+1))*isign*cylm(klm)*dsqrt(setf)*2.0d0
      return
      end

      SUBROUTINE KTRNSF(LAT,KKX,KKY,KKZ,PKX,PKY,PKZ)
C     *****************************************************************
C     * TRANSFORMATION FROM SYMMETRY TO RECTANGULAR COORDINATES       *
C     *****************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/MSH/DKX,DKY,DKZ,DHX,BOA,COA,KX(2500),KY(2500),KZ(2500),
     1W(2500)
C
      GO TO (21,21,21,21,21,21,21,22), LAT
   21 PKX=DKX*KKX+DHX*KKY
      PKY=DKY*KKY
      PKZ=DKZ*KKZ
      DK=PKX*PKX+PKY*PKY+PKZ*PKZ
      IF (DK.LT.1.D-08) GO TO 20
      RETURN
   22 QKX=DKX*KKX+DHX*KKY 
      QKY=DKY*KKY 
      QKZ=DKZ*KKZ-1.5D0/DSQRT(3.D0)
      DK=QKX*QKX+QKY*QKY+QKZ*QKZ
      IF(DK.LT.1.D-08) GO TO 20
      SQ2=DSQRT(2.D0)
      SQ3=DSQRT(3.D0)
      SQ6=DSQRT(6.D0)
      PXX=QKX/SQ2+QKY/SQ6+QKZ/SQ3
      PXY=-QKX/SQ2+QKY/SQ6+QKZ/SQ3
      PXZ=-2.D0*QKY/SQ6+QKZ/SQ3
      RETURN
C
C     AVOID (0,0,0)
C    
   20 PKX=0.D0
      PKY=0.D0
      PKZ=0.01D0
      RETURN
      END 



      SUBROUTINE MESH(LAT,NPX,NPY,NPZ,NPT)
C
C     *****************************************************************
C     * CONSTRUCTION OF MESH IN K-SPACE                                *
C     * LAT=1 SIMPLE CUBIC                                            *
C     * LAT=2 FACE CENTRED CUBIC                                      *
C     * LAT=3 BODY CENTRED CUBIC                                      *
C     * LAT=4 HEXAGONAL CLOSE PACKED                                  *
C     * LAT=5 SIMPLE TETRAGONAL                                        *
C     * LAT=6 BODY CENTRED TETRAGONAL                                 *
C     * LAT=7 ALPHA URANIUM                                           *
C     * LAT=8 TRIGONAL ZONE FOR ANTIFERROMAGNETIC OXIDES              *
C     *****************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER X,Y,Z
      COMMON /MSH/DKX,DKY,DKZ,DHX,BOA,COA,KX(2500),KY(2500),KZ(2500),
     1WW(2500)
    1 FORMAT(1H ,10X,'NUMBER OF POINTS ON MESH EXCEEDS DIMENSION OF ARRA
     1YS',2I5)
C
      NDIM=2500
      NPXM=NPX-1
      NPYM=NPY-1
      NPZM=NPZ-1

      GO TO (20,21,22,23,24,25,26,27), LAT
   20 NP=0
C
C     SIMPLE CUBIC IRREDUCIBLE ZONE DEFINED BY
C       .LE. KZ .LE. KX .LE. KY .LE. PI/A
C
      DKX=1.D0/NPYM
      DKY=DKX
      DKZ=DKX
      DHX=0.D0
      DO 30 J=1,NPY
      Y=J-1
      DO 30 I=1,J
      X=I-1
      DO 30 K=1,I
      Z=K-1
      NP=NP+1
      IF (NP.GT.NDIM) GO TO 999
      KX(NP)=X
      KY(NP)=Y
      KZ(NP)=Z
      W=48.D0
      IF (X.EQ.0.OR.Y.EQ.0.OR.Z.EQ.0) W=W/2.D0
      IF (X.EQ.Y.AND.Y.EQ.Z) W=W/3.D0
      IF (X.EQ.0) W=W/2.D0
      IF (Y.EQ.0) W=W/2.D0
      IF (Z.EQ.0) W=W/2.D0
      IF (X.EQ.NPYM) W=W/2.D0
      IF (Y.EQ.NPYM) W=W/2.D0
      IF (Z.EQ.NPYM) W=W/2.D0
      IF (X+Y+Z.EQ.3*NPYM) W=W/1.D0
      IF (X+Y+Z.EQ.0) W=W/1.D0
      WW(NP)=W
   30 CONTINUE
      NPT=NP
      RETURN
   21 NP=0

C
C     FCC IRREDUCIBLE ZONE DEFINED BY
C       .LE. KZ .LE. KX .LE. KY .LE. 2PI/A
C       +X +Z .LE. 3PI/A
C
      IF(NPYM.NE.4*(NPYM/4)) GO TO 998
      DKX=2.D0/NPYM
      DKY=DKX
      DKZ=DKX
      DHX=0.D0
      NPX=NPY
      NPH=NPY/2+1
      NPTH=(NPYM/2)**3
      DO 31 I=1,NPX
      Y=I-1
      M1=NPX-I+NPH
      NQY=MIN0(I,M1)
      DO 31 J=1,NQY
      X=J-1
      M5=M1-J+1
      NPZ=MIN0(J,M5)
      DO 31 K=1,NPZ
      Z=K-1
      NP=NP+1
      IF(NP.GT.NDIM) GO TO 999
      KX(NP)=X
      KY(NP)=Y
      KZ(NP)=Z
      W=48.D0
      IF(X.EQ.Y.OR.X.EQ.Z.OR.Y.EQ.Z) W=W/2.D0
      IF(X.EQ.Y.AND.Y.EQ.Z) W=W/3.D0
      IF(Z.EQ.0) W=W/2.D0
      IF(X.EQ.0) W=W/2.D0
      IF(Y.EQ.NPYM) W=W/2.D0
      IF(X+Y+Z.EQ.0) W=1.D0
      WW(NP)=W
   31 CONTINUE  
      NPT=NP
      RETURN


      
   22 NP=0
C     BCC IRREDUCIBLE ZONE DEFINED BY
      PRINT *,'DOING BCC CELL MESH'
      IF(NPYM.NE.2*(NPYM/2)) GO TO 998
      DKX=2.D0/NPYM
      DKY=DKX
      DKZ=DKX
      DHX=0.D0
      NPH=NPY/2+1
      NPYH=NPYM/2
      DO 32 I=1,NPY
      X=I-1
      JM=MIN0(I,NPY-I+1)
      DO 32 J=1,JM
      Y=J-1
      DO 32 K=1,J
      Z=K-1
      NP=NP+1
      IF(NP.GT.NDIM) GO TO 999
      KX(NP)=X
      KY(NP)=Y
      KZ(NP)=Z
      W=48.D0
      IF(X.EQ.Y.OR.X.EQ.Z.OR.Y.EQ.Z) W=W/2.D0
      IF(X.EQ.Y.AND.Y.EQ.Z) W=W/3.D0
      IF(Z.EQ.0) W=W/2.D0
      IF(Y.EQ.0) W=W/2.D0
      IF (X+Y.EQ.NPYM) W=W/2.D0
c      IF(Y.EQ.NPYM) W=W/2.D0
      IF(Y .EQ.Z.AND.X+Y .EQ.NPYM) W=8.D0
      IF(Z.EQ.NPYH) W=2.D0
      IF(X.EQ.NPYM) W=1.D0 
      IF(X+Y+Z.EQ.0) W=1.D0
      WW(NP)=W
   32 CONTINUE  
      NPT=NP
      RETURN
   23 NP=0



C HEXAGONAL CLOSE PACKED IRREDUCIBLE ZONE DEFINED BY
C 0 .LE. 2KY .LE. KX .LE. 4/3 PI/A
C 0 .LE. KZ .LE. A/C PI/A
C

C HEXAGONAL CLOSE PACKED IRREDUCIBLE ZONE DEFINED BY
C 0 .LE. 2KY .LE. KX .LE. 4/3 PI/A
C 0 .LE. KZ .LE. A/C PI/A
C
      IF(NPYM.NE.2*(NPYM/2)) GO TO 998
      NPYH=NPYM/2
      NPH=NPY/2+1
      NPZM=NPZ-1
      AOC=1.D0/COA
      DKX=4.D0/3.D0/NPYM
      DKY=DSQRT(3.D0)/2.D0*DKX
      DKZ=AOC/NPZM
      DHX=-0.5D0*DKX
      DO 33 K=1,NPZ
      Z=K-1
      DO 33 J=1,NPH
      Y=J-1
      IM=2*Y+1
      DO 33 I=IM,NPY
      X=I-1
      NP=NP+1
      IF(NP.GT.NDIM) GO TO 999
      KX(NP)=X
      KY(NP)=Y
      KZ(NP)=Z
      W=24.D0
      IF(Z.EQ.0) W=W/2.D0
      IF(Z.EQ.0) W=W/2.D0
      IF(Y.EQ.0) W=W/2.D0
      IF(X.EQ.0) W=W/3.D0
      IF(X.EQ.2*Y) W=W/2.D0
      IF(Z.EQ.NPZM) W=W/2.D0
      IF(X.EQ.NPYM.AND.Y.EQ.0) W=W/3.D0
      IF(X.EQ.NPYM.AND.Y.NE.0) W=W/2.D0
      IF(X+Y.EQ.0.AND.Z.EQ.NPZM) W=1.D0
      WW(NP)=W
 33   CONTINUE
      NPT=NP
      RETURN
 24   NP=0



C
C SIMPLE TETRAGONAL IRREDUCIBLE ZONE DEFINED BY
C 0 .LE. KY .LE. KX .LE. PI/A
C 0 .LE. KZ .LE. A/C PI/A
C
      AOC=1.D0/COA
      DKX=1.D0/NPYM
      DKY=DKX
      DKZ=AOC/NPZM
      DHX=0.D0
      DO 34 I=1,NPY
      X=I-1
      DO 34 J=1,I
      Y=J-1
      DO 34 K=1,NPZ
      Z=K-1
      NP=NP+1
      IF(NP.GT.NDIM) GO TO 999
      KX(NP)=X
      KY(NP)=Y
      KZ(NP)=Z
      WW(NP)=0.
 34   CONTINUE
      NPT=NP
      RETURN
 25   NP=0
C
C BODY CENTRED TETRAGONAL IRREDUCIBLE ZONE DEFINED BY
C 0 .LE. KY .LE. KX .LE. PI/A
C 0 .LE. KZ .LE. 2*A/C PI/A
C

      AOC=1.D0/COA
      DKX=1.D0/NPYM
      DKY=DKX
      DKZ=2.D0*AOC/NPZM
      DHX=0.D0
      DO 35 I=1,NPY
      X=I-1
      DO 35 J=1,I
      Y=J-1
      DO 35 K=1,NPZ
      Z=K-1
      NP=NP+1
      IF(NP.GT.NDIM) GO TO 999
      KX(NP)=X
      KY(NP)=Y
      KZ(NP)=Z
      W=16.D0
      IF(X.EQ.0) W=W/2.D0
      IF(Y.EQ.0) W=W/2.D0
      IF(Z.EQ.0) W=W/2.D0
      IF(X.EQ.Y) W=W/2.D0
      IF(X.EQ.NPYM) W=W/2.D0
      IF(Y.EQ.NPYM) W=W/2.D0
      IF(Z.EQ.NPZM) W=W/2.D0
      WW(NP)=W
 35   CONTINUE
      NPT=NP
      RETURN
 26   NP=0
C
C ALPHA URANIUM
C
      AOC=1.D0/COA
      DKX=1.D0/NPXM
      DKY=2.D0/BOA/NPYM
      DKZ=AOC/NPZM
      DHX=0.D0
      DO 36 I=1,NPX
      X=I-1
      DO 36 J=1,NPY
      Y=J-1
      DO 36 K=1,NPZ
      Z=K-1
      NP=NP+1
      IF(NP.GT.NDIM) GO TO 999
      KX(NP)=X
      KY(NP)=Y
      KZ(NP)=Z
      WW(NP)=0.D0
 36   CONTINUE
      NPT=NP
      RETURN

C TRIGONAL ZONE FOR ANTIFERROMAGNETIC OXIDES
C
C 0. < KX < 4.*SQRT(2.)/3.
C 0. < KY < SQRT(2./3.) * (PI/A)
C -1.5/SQRT(3.) < KZ < 1.5/SQRT(3.)
C
 27   NP=0
      IF(NPYM.NE.2*(NPYM/2)) GO TO 998
      NPYH=NPYM/2
      NPH=NPY/2+1
      NPZM=NPZ-1
      DKX=4.D0*DSQRT(2.D0)/3.D0/NPYM
      DKY=DSQRT(3.D0)/2.D0*DKX
      DKZ=3.D0/DSQRT(3.D0)/NPZM
      DHX=-0.5D0*DKX
      DO 37 K=1,NPZ
      Z=K-1
      DO 37 J=1,NPH
      Y=J-1
      IM=2*Y+1
      DO 37 I=IM,NPY
      X=I-1
      NP=NP+1
      IF(NP.GT.NDIM) GO TO 999
      KX(NP)=X
      KY(NP)=Y
      KZ(NP)=Z
      WW(NP)=0
 37   CONTINUE
      NPT=NP
      RETURN
 998  WRITE(6,2) NPY
 2    FORMAT(1H , '*** WRONG NPY =, NPY=',I4,' ***')
      STOP

      

 999  WRITE(6,1) NP,NDIM
      STOP
      END

      
      


      









      
      
C     FIRST PART OF PAGE 140
      
      
    
c     end program
      
