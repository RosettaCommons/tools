
      PROGRAM TEST

	implicit none
      INTEGER    IP(9),IP2312(4),I,J,K,L,M1,M,IER,MODE,N
      REAL      W(7),X(3,7),Y(3,7),U(3,3),T(3),RMS,SIGMA
	REAL*16 H,G,D, SQRTH,CTH,STH, SS(6), A(3,3),
     1		SS1, SS2, SS3, SS4, SS5, SS6, P
      REAL*8     R(3,3),XC(3),YC(3),WC,B(3,3),E0,
     1 E(3),E1,E2,E3,SPUR,DET,COF,TOL,
     2 RR(6),RR1,RR2,RR3,RR4,RR5,RR6,
     3 ZERO,ONE,TWO,THREE,SQRT3
      EQUIVALENCE (RR1,RR(1)),(RR2,RR(2)),(RR3,RR(3)),
     1        (RR4,RR(4)),(RR5,RR(5)),(RR6,RR(6)),
     2        (SS1,SS(1)),(SS2,SS(2)),(SS3,SS(3)),
     3        (SS4,SS(4)),(SS5,SS(5)),(SS6,SS(6)),
     4        (E1,E(1)),(E2,E(2)),(E3,E(3))
      DATA SQRT3,TOL/1.73205080756888D+00, 1.0D-2/
      DATA ZERO,ONE,TWO,THREE/0.0D+00, 1.0D+00, 2.0D+00, 3.0D+00/
      DATA IP/1,2,4,  2,3,5,  4,5,6/
      DATA IP2312/2,3,1,2/
     
      WC=ZERO
      RMS=0.0
      E0=ZERO
      DO 1 I=1,3
      XC(I)=ZERO
      YC(I)=ZERO
      T(I)=0.0
      DO 1 J=1,3
      D=ZERO
      IF (I.EQ.J)D=ONE
      U(I,J)=D
      A(I,J)=D
1     R(I,J)=ZERO
      IER=-1
      IF (N.LT.1)RETURN
C**** DETERMINE CENTROIDS OF BOTH VECTOR SETS X AND Y
      IER=-2
      DO 2 M=1,N
      IF (W(M).LT.0.0)RETURN
      WC=WC+W(M)
      DO 2 I=1,3
      XC(I)=XC(I)+W(M)*X(I,M)
2     YC(I)=YC(I)+W(M)*Y(I,M)
      IF (WC.LE.ZERO)RETURN
      DO 3 I=1,3
      XC(I)=XC(I)/WC
3     YC(I)=YC(I)/WC
C**** DETERMINE CORRELATION MATRIX R BETWEEN VECTOR SETS Y AND X
      DO 4 M=1,N
      DO 4 I=1,3
      E0=E0+W(M)*((X(I,M)-XC(I))**2+(Y(I,M)-YC(I))**2)
      D=W(M)*(Y(I,M)-YC(I))
      DO 4 J=1,3
4     R(I,J)=R(I,J)+D*(X(J,M)-XC(J))
C**** CALCULATE DETERMINANT OF R(I,J)
      DET=R(1,1)*(R(2,2)*R(3,3)-R(2,3)*R(3,2))
     1   -R(1,2)*(R(2,1)*R(3,3)-R(2,3)*R(3,1))
     2   +R(1,3)*(R(2,1)*R(3,2)-R(2,2)*R(3,1))
      SIGMA=DET
C**** FORM UPPER TRIANGLE OF TRANSPOSED(R)*R
      M=0
      DO 5 J=1,3
      DO 5 I=1,J
      M=M+1
5     RR(M)=R(1,I)*R(1,J)+R(2,I)*R(2,J)+R(3,I)*R(3,J)
C***************** EIGENVALUES *****************************************
C**** FORM CHARACTERISTIC CUBIC  X**3-3*SPUR*X**2+3*COF*X-DET=0
      SPUR=(RR1+RR3+RR6)/THREE
      COF=(RR3*RR6-RR5*RR5+RR1*RR6-RR4*RR4+RR1*RR3-RR2*RR2)/THREE
      DET=DET*DET
      DO 6 I=1,3
6     E(I)=SPUR
      IF (SPUR.LE.ZERO)GO TO 40
C**** REDUCE CUBIC TO STANDARD FORM Y**3-3HY+2G=0 BY PUTTING X=Y+SPUR
      D=SPUR*SPUR
      H=D-COF
      G=(SPUR*COF-DET)/TWO-SPUR*H
C**** SOLVE CUBIC. ROOTS ARE E1,E2,E3 IN DECREASING ORDER
      IF (H.LE.ZERO)GO TO 8
      SQRTH=QSQRT(H)
7      D=H*H*H-G*G
9      IF (D.LT.ZERO)D=ZERO
      D=QATAN2(QSQRT(D),-G)/THREE
      CTH=SQRTH*QCOS(D)
      STH=SQRTH*SQRT3*QSIN(D)
      E1=SPUR+CTH+CTH
      E2=SPUR-CTH+STH
      E3=SPUR-CTH-STH
      IF (MODE)10,50,10
C.....HANDLE SPECIAL CASE OF 3 IDENTICAL ROOTS
8     IF (MODE)30,50,30
C**************** EIGENVECTORS *****************************************
10    DO 15 L=1,3,2
      D=E(L)
      SS1=(D-RR3)*(D-RR6)-RR5*RR5
      SS2=(D-RR6)*RR2+RR4*RR5
      SS3=(D-RR1)*(D-RR6)-RR4*RR4
      SS4=(D-RR3)*RR4+RR2*RR5
      SS5=(D-RR1)*RR5+RR2*RR4
      SS6=(D-RR1)*(D-RR3)-RR2*RR2
      J=1
      IF (QABS(SS1).GE.QABS(SS3))GO TO 12
      J=2
      IF (QABS(SS3).GE.QABS(SS6))GO TO 13
11    J=3
      GO TO 13
12    IF (QABS(SS1).LT.QABS(SS6))GO TO 11
13    D=ZERO
      J=3*(J-1)
      DO 14 I=1,3
      K=IP(I+J)
      A(I,L)=SS(K)
14    D=D+SS(K)*SS(K)
      IF (D.GT.ZERO)D=ONE/QSQRT(D)
      DO 15 I=1,3
15    A(I,L)=A(I,L)*D
      D=A(1,1)*A(1,3)+A(2,1)*A(2,3)+A(3,1)*A(3,3)
      M1=3
      M=1
      IF ((E1-E2).GT.(E2-E3))GO TO 16
      M1=1
      M=3
16    P=ZERO
      DO 17 I=1,3
      A(I,M1)=A(I,M1)-D*A(I,M)
17    P=P+A(I,M1)**2
      IF (P.LE.TOL)GO TO 19
      P=ONE/QSQRT(P)
      DO 18 I=1,3
18    A(I,M1)=A(I,M1)*P
      GO TO 21
19    P=ONE
      DO 20 I=1,3
      IF (P.LT.QABS(A(I,M)))GO TO 20
      P=QABS(A(I,M))
      J=I
20    CONTINUE
      K=IP2312(J)
      L=IP2312(J+1)
      P=QSQRT(A(K,M)**2+A(L,M)**2)
      IF (P.LE.TOL)GO TO 40
      A(J,M1)=ZERO
      A(K,M1)=-A(L,M)/P
      A(L,M1)= A(K,M)/P
21    A(1,2)=A(2,3)*A(3,1)-A(2,1)*A(3,3)
      A(2,2)=A(3,3)*A(1,1)-A(3,1)*A(1,3)
      A(3,2)=A(1,3)*A(2,1)-A(1,1)*A(2,3)
C****************** ROTATION MATRIX ************************************
30    DO 32 L=1,2
      D=ZERO
      DO 31 I=1,3
      B(I,L)=R(I,1)*A(1,L)+R(I,2)*A(2,L)+R(I,3)*A(3,L)
31    D=D+B(I,L)**2
      IF (D.GT.ZERO)D=ONE/QSQRT(D)
      DO 32 I=1,3
32    B(I,L)=B(I,L)*D
      D=B(1,1)*B(1,2)+B(2,1)*B(2,2)+B(3,1)*B(3,2)
      P=ZERO
      DO 33 I=1,3
      B(I,2)=B(I,2)-D*B(I,1)
33    P=P+B(I,2)**2
      IF (P.LE.TOL)GO TO 35
      P=ONE/QSQRT(P)
      DO 34 I=1,3
34    B(I,2)=B(I,2)*P
      GO TO 37
35    P=ONE
      DO 36 I=1,3
      IF (P.LT.DABS(B(I,1)))GO TO 36
      P=DABS(B(I,1))
      J=I
36    CONTINUE
      K=IP2312(J)                                                       
      L=IP2312(J+1)
      P=DSQRT(B(K,1)**2+B(L,1)**2)
      IF (P.LE.TOL)GO TO 40
      B(J,2)=ZERO
      B(K,2)=-B(L,1)/P
      B(L,2)= B(K,1)/P
37    B(1,3)=B(2,1)*B(3,2)-B(2,2)*B(3,1)
      B(2,3)=B(3,1)*B(1,2)-B(3,2)*B(1,1)
      B(3,3)=B(1,1)*B(2,2)-B(1,2)*B(2,1)
      DO 39 I=1,3
      DO 39 J=1,3
39    U(I,J)=B(I,1)*A(J,1)+B(I,2)*A(J,2)+B(I,3)*A(J,3)
C****************** TRANSLATION VECTOR *********************************
40    DO 41 I=1,3
41    T(I)=YC(I)-U(I,1)*XC(1)-U(I,2)*XC(2)-U(I,3)*XC(3)
C********************** RMS ERROR **************************************
50    DO 51 I=1,3
      IF (E(I).LT.ZERO)E(I)=ZERO
51    E(I)=DSQRT(E(I))
      IER=0
      IF (E2.LE.(E1*1.0D-05))IER=-1
      D=E3
      IF (SIGMA.GE.0.0)GO TO 52
      D=-D
      IF ((E2-E3).LE.(E1*1.0D-05))IER=-1
52    D=D+E2+E1
      RMS=E0-D-D
      IF (RMS.LT.0.0)RMS=0.0
     
      RETURN
      END
C.....END U3B...........................................................

      END
