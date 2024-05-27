      SUBROUTINE LIMA0(NN,FALFB)

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
