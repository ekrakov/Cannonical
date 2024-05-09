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
      Y=J
      DO 30 I=1,J
      X=I
      DO 30 K=1,I
      Z=K
      NP=NP+1
      IF (NP.GT.NDIM) GO TO 999
      KX(NP)=X
      KY(NP)=Y
      KZ(NP)=Z
      W=4.8D0
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


      IF(NPYM.NE.2*(NPYM/2)) GO TO 998
      DKX=2.D0/NPYM
      DKY=DKX
      DKZ=DKX
      DHX=0.D0
      NPH=NPY/2+1
      NPYH=NPYM/2

      DO 32 I=1,NPX
      X=I-1
      M1=NPX-I+NPH
      NQY=MIN0(I,M1)
      DO 31 J=1,NQY
      X=J-1
C     check this
      


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
      IF(X+Y+Z.EQ.O) W=1.D0
      WW(NP)=W
   32 CONTINUE  
      NPT=NP
      RETURN
   23 NP=0

CCC STILL MISSING PAGE 152-153
      STOP
  999 WRITE(6,1) NP,NDIM
      STOP
C      END
C      
