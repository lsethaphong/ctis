CAnd substituting the value for R x given by eq. (28) and
      SUBROUTINE GINV2 (A,U,AFLAG,ATEMP,MR,NR,NC)
C
C THIS ROUTINE CALCULATES THF GENERALIZED INVERSE OF A
C AND STORES THE TRANSPOSE OF IT IN A.
C MR#FIRST DIMENSION Nno OF A.
C NR # NO. OF POI;S IN A
C NC # NO, C~F COLUMN� Itl A
C U IS THE BOOKKEEPING ~ATPIX.
C AFLAG AND ATEMP ARE TEMPORARY WORKING STORAGE.
C
      DIMENSION A(MR,NC),U(NC,NC),AFLAG(NC),ATEMP(r,NC)
      DO IO I # I,NC
      DO 5 J # I,NC
5     U(I,J) # D.q
      I0 U(I,I}#1.O
      FAC # DOT(MP,NR+A~I.I}
      FAC# I . O/SORT ( FAC )
      DO I 5 I # I ,NR
15    A(I,II#AII,Ii+~AC
      DO 20 I # I ,NC
20    U(I,I}#U( I,I)~cAC
      AFLAG( I }#I.~
C
C DEPENDENT COL TOLERANCE FOR N BIT FLOATING POINT FRACTION
C
      N # 2T
      TOL # (10. * O.~**N)**~
      DO ION J # 2,pNC
      DOT I # DOTIMR,NR,A,J~J)
      JMI#J-I
      DO 5n L#I ~7
      DO 9r~ K#I,JMI
30    ATEMP(K} # DOT(r~R,NR+A,,I,K)
      DO 45 K#I ,JMI
      DO 35 I # I ,NP
35    A( I ,J)#A( I ,J)-ATEMP(K)*~( I ,K)*AFLAG(K)
      DO 40 I # 1 ,NC
40    U( I ,J )#U( I ,JI-ATEMP(KI*U( I ~K )
45    CONTINUE
50    CONTINUE
      DOT2 # DOT(MR,NR,A,J~J)
      IF((DOT2/D~TI) - TOL) m5.55~70
55    DO 60 I#I,JMI
      ATEMP ( I }#3.q
      DO 60 K#I ,I
60    ATEMP( I }#ATEMP( I )+U(K, I)*U(K,J}
      DO 65 I # I ,NR
      A{ I.J)#O�F]
      DO 65 K#I,JM[
65    A( I ,JI#A( I ,JI-A( I ,KI*ATEMP(K)*AFLAG(K)
      AFLAG(J ) #I].[~
      FAC # DOT(NC,NC,U,J,J)
      FAC# I . D/SORT ( FAC }
      GO TO 75
7 r] AFLAGfJ)#I.n
FAC#1 .O/SORT(DOT2)
75    DO 80 I # I',NR
80    A( I,J)#A( I,J)*FAC
      DO 85 I # I ,NC
85    U( I ~J)#U( I,J)*FAC
10O   CONTINUE
      DO 130 J#l +NC
      DO 130 I#l ,NR
FAC#n.q
DO 12Q K#J,NC
120 FAC#FAC+A( I ,K}*U{J,K1
1'3D AiI*J)#FAC
RETURN
FNO
FUNCTION DmTIMP~NP+A,JC~KC)
COMPUTES THE INNPP RPODUCT OF COLUMNS JC AND KC
OF MATPlX A.
DIMENSION A(MR~I)
DOT#P.n
DO 5 I # I,NR
DOT # DOT + A(I,JC)*A(I,KC}
RFTUPN
END
