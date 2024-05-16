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
