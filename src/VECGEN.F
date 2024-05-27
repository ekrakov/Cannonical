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

    1 FORMAT(1H,//, 15X,'SHELL NUMBER',I5,' WITH',I5,' POINTS',//)
    2 FORMAT(1H 10X,I5,4F10.6)
    3 FORMAT(1H ,10X,'(',F10.5,',',F10.5,',',F10.5,')',/)
    4 FORMAT(1H1 ,14X,'RESULT FROM VECGEN FOR REAL SPACE VECTORS',//,
     114X,'NO',5X,'SX',8X,'SY',8X,'SZ',8X,'D',/)
    5 FORMAT(1H ,10X,'RESULT FROM VECGEN FOR RECIPROCAL SPACE VECTORS',
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
      I2=1+MOD(I,3)
      CALL CROSS(BSX(I1),BSY(I1),BSZ(I1),BSX(I2),BSY(I2),BSZ(I2),
     1BKX(I),BKY(I),BKZ(I))
      VOL=DABS(BSX(1)*BKY(1)+BSY(1)*BKZ(1))
   20 CONTINUE
      
      DO 21 I=1,3
      BKX(I)=BKX(I)*VOL*TWOPI
      BKY(I)=BKY(I)*VOL*TWOPI
      BKZ(I)=BKZ(I)*VOL*TWOPI
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
      PRINT * 'NRVAL NL:', NR
      IF (NR.GT.500) GO TO 33
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
