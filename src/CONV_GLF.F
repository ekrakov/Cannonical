C
C     *****************************************************************
C     * CALCULATION OF CONVERGENCE FUNCTIONS                          *
C     *****************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 CIM,CCIM
      COMMON /B1/CIM,CCIM,C(16,16),FAC(13),DAC(13),NGLN(10),
     1LAMDA,WTOF,FTLM,PL3OM,SVS,VOL,SOO,SWSL,FTPLAM,RMAX,GMAX,ICNVG,
     2NUMVR,NUMVG,NQ,NRO,NL2,NL2P,KLM0,NL,NLM,LQM0,NOWRT,IDIM
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
      H10=(1.D00+ALFA)/ALFA/ALFA*EXA
      H1=H10
      H2=H11
      GLN=H1+H2
      DO 21 NP=2,NN
      N=NP-1
      H3=(ALFASF/N*H1-(1.5D00-N)*H2)/NP
      GLN=GLN+H3
      H1=H2
   21 H2=H3
      IF(LN.EQ.0) RETURN
C
C     FOR L .GT. 1
C
      L=1
   23 L=L+1
      LP=L-1
      LM=L-1
      H21=H10+(2*LM-1)*H11/ALFA
      H20=H10+(2*LM+1)*H10/ALFA
      H1=H20
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
      IF(L.GE.LN) RETURN
      H00=H10
      H01=H20
      H10=H11
      H11=H21
      GO TO 23
      END
