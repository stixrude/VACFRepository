      FUNCTION GDOT(A,B,GIJ)
C     ***********************************************************
C     *****
C     *****   DOT PRODUCT OF TWO VECTORS WITH METRIC GIJ
C     *****
C     ***********************************************************
      DIMENSION A(3),B(3),GIJ(3,3)
      GDOT=0.0D0
      DO 1 J=1,3
      DO 1 I=1,3
    1 GDOT=GDOT+A(I)*GIJ(I,J)*B(J)
      RETURN
      END
