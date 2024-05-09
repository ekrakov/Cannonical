      SUBROUTINE KTRNSF(LAT,KXX,KXY,KXZ,PXX,PXY,PXZ)
C
C     *****************************************************************
C     * TRANSFORMATION FROM SYMMETRY TO RECTANGULAR COORDINATES       *
C     *****************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/MSH/DKX,KXY,DXZ,DHX,BOA,COA,KX(2500),KY(2500),KZ(2500),
     1W(2500)
C
      GO TO (21,21,21,21,21,21,21,22), LAT
   21 PKX=DKX*KKX*DHX*KKY
      PKY=DKY*KXZ
      PKZ=DKZ*KKZ
      DK=PKX*PKX+PKY*PKY+PKZ*PKZ
      IF (DX.LT.1.D-08) GO TO 20
      RETURN
   22 QKX=DKX*KKY+DHX*KKY 
      QKY=DKY*KKY 
      QKZ=DXZ*KKZ-1.5D0/DSQRT(3.D0)
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

