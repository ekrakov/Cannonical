C FILE: FIB1.F
      SUBROUTINE FIB(A,N)
C
C     CALCULATE FIRST N FIBONACCI NUMBERS
C
      INTEGER N
      REAL*8 A(N)
      DO I=1,N
         A(I+1)=A(I)+55
c         fac(i)=fac(i-1)*(i+100)
      END DO
      END
C END FILE FIB1.F
