      SUBROUTINE SGG_TTBG(P1,ANS)
C  ifort -O0 -check all -c gg_ttbg.f -I/home/schulze/lib/HELAS-3.0/
C FUNCTION GENERATED BY MADGRAPH
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C AND HELICITIES
C FOR THE POINT IN PHASE SPACE P(0:3,NEXTERNAL)
C
C FOR PROCESS : g g -> t t~ g
C
C Crossing   1 is g g -> t t~ g
      IMPLICIT NONE
C
C CONSTANTS
C
      INTEGER    NEXTERNAL,   NCOMB,     NCROSS
      PARAMETER (NEXTERNAL=5, NCOMB= 32, NCROSS=  1)
      INTEGER    THEL
      PARAMETER (THEL=NCOMB*NCROSS)
C
C ARGUMENTS
C
      REAL*8 P1(0:3,NEXTERNAL),ANS(NCROSS)
C
C LOCAL VARIABLES
C
      INTEGER NHEL(NEXTERNAL,NCOMB),NTRY
      REAL*8 T, P(0:3,NEXTERNAL)
      REAL*8 GG_TTBG
      INTEGER IHEL,IDEN(NCROSS),IC(NEXTERNAL,NCROSS)
      INTEGER IPROC,JC(NEXTERNAL)
      LOGICAL GOODHEL(NCOMB,NCROSS)
      DATA GOODHEL/THEL*.FALSE./
      DATA NTRY/0/
      DATA (NHEL(IHEL,  1),IHEL=1,5) / -1, -1, -1, -1, -1/
      DATA (NHEL(IHEL,  2),IHEL=1,5) / -1, -1, -1, -1,  1/
      DATA (NHEL(IHEL,  3),IHEL=1,5) / -1, -1, -1,  1, -1/
      DATA (NHEL(IHEL,  4),IHEL=1,5) / -1, -1, -1,  1,  1/
      DATA (NHEL(IHEL,  5),IHEL=1,5) / -1, -1,  1, -1, -1/
      DATA (NHEL(IHEL,  6),IHEL=1,5) / -1, -1,  1, -1,  1/
      DATA (NHEL(IHEL,  7),IHEL=1,5) / -1, -1,  1,  1, -1/
      DATA (NHEL(IHEL,  8),IHEL=1,5) / -1, -1,  1,  1,  1/
      DATA (NHEL(IHEL,  9),IHEL=1,5) / -1,  1, -1, -1, -1/
      DATA (NHEL(IHEL, 10),IHEL=1,5) / -1,  1, -1, -1,  1/
      DATA (NHEL(IHEL, 11),IHEL=1,5) / -1,  1, -1,  1, -1/
      DATA (NHEL(IHEL, 12),IHEL=1,5) / -1,  1, -1,  1,  1/
      DATA (NHEL(IHEL, 13),IHEL=1,5) / -1,  1,  1, -1, -1/
      DATA (NHEL(IHEL, 14),IHEL=1,5) / -1,  1,  1, -1,  1/
      DATA (NHEL(IHEL, 15),IHEL=1,5) / -1,  1,  1,  1, -1/
      DATA (NHEL(IHEL, 16),IHEL=1,5) / -1,  1,  1,  1,  1/
      DATA (NHEL(IHEL, 17),IHEL=1,5) /  1, -1, -1, -1, -1/
      DATA (NHEL(IHEL, 18),IHEL=1,5) /  1, -1, -1, -1,  1/
      DATA (NHEL(IHEL, 19),IHEL=1,5) /  1, -1, -1,  1, -1/
      DATA (NHEL(IHEL, 20),IHEL=1,5) /  1, -1, -1,  1,  1/
      DATA (NHEL(IHEL, 21),IHEL=1,5) /  1, -1,  1, -1, -1/
      DATA (NHEL(IHEL, 22),IHEL=1,5) /  1, -1,  1, -1,  1/
      DATA (NHEL(IHEL, 23),IHEL=1,5) /  1, -1,  1,  1, -1/
      DATA (NHEL(IHEL, 24),IHEL=1,5) /  1, -1,  1,  1,  1/
      DATA (NHEL(IHEL, 25),IHEL=1,5) /  1,  1, -1, -1, -1/
      DATA (NHEL(IHEL, 26),IHEL=1,5) /  1,  1, -1, -1,  1/
      DATA (NHEL(IHEL, 27),IHEL=1,5) /  1,  1, -1,  1, -1/
      DATA (NHEL(IHEL, 28),IHEL=1,5) /  1,  1, -1,  1,  1/
      DATA (NHEL(IHEL, 29),IHEL=1,5) /  1,  1,  1, -1, -1/
      DATA (NHEL(IHEL, 30),IHEL=1,5) /  1,  1,  1, -1,  1/
      DATA (NHEL(IHEL, 31),IHEL=1,5) /  1,  1,  1,  1, -1/
      DATA (NHEL(IHEL, 32),IHEL=1,5) /  1,  1,  1,  1,  1/
      DATA (  IC(IHEL,  1),IHEL=1,5) /  1,  2,  3,  4,  5/
      DATA (IDEN(IHEL),IHEL=  1,  1) / 256/
C ----------
C BEGIN CODE
C ----------
      NTRY=NTRY+1
      DO IPROC=1,NCROSS
      CALL SWITCHMOM(P1,P,IC(1,IPROC),JC,NEXTERNAL)
      DO IHEL=1,NEXTERNAL
         JC(IHEL) = +1
      ENDDO
      ANS(IPROC) = 0D0
      DO IHEL=1,NCOMB
          IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
             T=GG_TTBG(P ,NHEL(1,IHEL),JC(1))
!              print *, "HEL",IHEL,":",T
             ANS(IPROC)=ANS(IPROC)+T
              IF (T .GT. 0D0 .AND. .NOT. GOODHEL(IHEL,IPROC)) THEN
                  GOODHEL(IHEL,IPROC)=.TRUE.
C             WRITE(*,*) IHEL,T
              ENDIF
          ENDIF
      ENDDO
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END






      REAL*8 FUNCTION GG_TTBG(P,NHEL,IC)
C
C FUNCTION GENERATED BY MADGRAPH
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C
C FOR PROCESS : g g -> t t~ g
C
      IMPLICIT NONE
C
C CONSTANTS
C
      INTEGER    NGRAPHS,    NEIGEN,    NEXTERNAL
      PARAMETER (NGRAPHS=  18,NEIGEN=  6,NEXTERNAL=5)
      INTEGER    NWAVEFUNCS     , NCOLOR
      PARAMETER (NWAVEFUNCS=  23, NCOLOR=   6)
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
C
C ARGUMENTS
C
      REAL*8 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL), IC(NEXTERNAL)
C
C LOCAL VARIABLES
C
      INTEGER I,J
      COMPLEX*16 ZTEMP
      REAL*8 DENOM(NCOLOR), CF(NCOLOR,NCOLOR)
      COMPLEX*16 AMP(NGRAPHS), JAMP(NCOLOR)
      COMPLEX*16 W(6,NWAVEFUNCS)
C
C GLOBAL VARIABLES
C
      INCLUDE 'coupl.inc'
C
C COLOR DATA
C
      DATA Denom(1  )/            9/
      DATA (CF(i,1  ),i=1  ,6  ) /    64,   -8,    1,    1,   10,   -8/
C               T[3,4,2,1,5]
      DATA Denom(2  )/            9/
      DATA (CF(i,2  ),i=1  ,6  ) /    -8,   64,   -8,   10,    1,    1/
C               T[3,4,2,5,1]
      DATA Denom(3  )/            9/
      DATA (CF(i,3  ),i=1  ,6  ) /     1,   -8,   64,    1,   -8,   10/
C               T[3,4,5,2,1]
      DATA Denom(4  )/            9/
      DATA (CF(i,4  ),i=1  ,6  ) /     1,   10,    1,   64,   -8,   -8/
C               T[3,4,1,5,2]
      DATA Denom(5  )/            9/
      DATA (CF(i,5  ),i=1  ,6  ) /    10,    1,   -8,   -8,   64,    1/
C               T[3,4,5,1,2]
      DATA Denom(6  )/            9/
      DATA (CF(i,6  ),i=1  ,6  ) /    -8,    1,   10,   -8,    1,   64/
C               T[3,4,1,2,5]
C ----------
C BEGIN CODE
C ----------
      CALL VXXXXX(P(0,1   ),ZERO ,NHEL(1   ),-1*IC(1   ),W(1,1   ))
      CALL VXXXXX(P(0,2   ),ZERO ,NHEL(2   ),-1*IC(2   ),W(1,2   ))
      CALL OXXXXX(P(0,3   ),TMASS ,NHEL(3   ),+1*IC(3   ),W(1,3   ))
      CALL IXXXXX(P(0,4   ),TMASS ,NHEL(4   ),-1*IC(4   ),W(1,4   ))
      CALL VXXXXX(P(0,5   ),ZERO ,NHEL(5   ),+1*IC(5   ),W(1,5   ))
      CALL FVOXXX(W(1,3   ),W(1,2   ),GG ,TMASS   ,TWIDTH  ,W(1,6   ))
      CALL JIOXXX(W(1,4   ),W(1,6   ),GG ,ZERO    ,ZERO    ,W(1,7   ))
      CALL VVVXXX(W(1,5   ),W(1,1   ),W(1,7   ),G ,AMP(1   ))
      CALL FVIXXX(W(1,4   ),W(1,1   ),GG ,TMASS   ,TWIDTH  ,W(1,8   ))
      CALL JVVXXX(W(1,5   ),W(1,2   ),G ,ZERO    ,ZERO    ,W(1,9   ))
      CALL IOVXXX(W(1,8   ),W(1,3   ),W(1,9   ),GG ,AMP(2   ))
      CALL FVOXXX(W(1,3   ),W(1,5   ),GG ,TMASS   ,TWIDTH  ,W(1,10  ))
      CALL IOVXXX(W(1,8   ),W(1,10  ),W(1,2   ),GG ,AMP(3   ))
      CALL FVOXXX(W(1,6   ),W(1,1   ),GG ,TMASS   ,TWIDTH  ,W(1,11  ))
      CALL IOVXXX(W(1,4   ),W(1,11  ),W(1,5   ),GG ,AMP(4   ))
      CALL IOVXXX(W(1,8   ),W(1,6   ),W(1,5   ),GG ,AMP(5   ))
      CALL FVIXXX(W(1,4   ),W(1,2   ),GG ,TMASS   ,TWIDTH  ,W(1,12  ))
      CALL JVVXXX(W(1,5   ),W(1,1   ),G ,ZERO    ,ZERO    ,W(1,13  ))
      CALL IOVXXX(W(1,12  ),W(1,3   ),W(1,13  ),GG ,AMP(6   ))
      CALL FVOXXX(W(1,3   ),W(1,1   ),GG ,TMASS   ,TWIDTH  ,W(1,14  ))
      CALL JIOXXX(W(1,4   ),W(1,14  ),GG ,ZERO    ,ZERO    ,W(1,15  ))
      CALL VVVXXX(W(1,5   ),W(1,2   ),W(1,15  ),G ,AMP(7   ))
      CALL IOVXXX(W(1,12  ),W(1,10  ),W(1,1   ),GG ,AMP(8   ))
      CALL FVOXXX(W(1,14  ),W(1,2   ),GG ,TMASS   ,TWIDTH  ,W(1,16  ))
      CALL IOVXXX(W(1,4   ),W(1,16  ),W(1,5   ),GG ,AMP(9   ))
      CALL IOVXXX(W(1,12  ),W(1,14  ),W(1,5   ),GG ,AMP(10  ))
      CALL JIOXXX(W(1,4   ),W(1,3   ),GG ,ZERO    ,ZERO    ,W(1,17  ))
      CALL VVVXXX(W(1,13  ),W(1,2   ),W(1,17  ),G ,AMP(11  ))
      CALL VVVXXX(W(1,1   ),W(1,9   ),W(1,17  ),G ,AMP(12  ))
      CALL JVVXXX(W(1,1   ),W(1,2   ),G ,ZERO    ,ZERO    ,W(1,18  ))
      CALL FVIXXX(W(1,4   ),W(1,18  ),GG ,TMASS   ,TWIDTH  ,W(1,19  ))
      CALL IOVXXX(W(1,19  ),W(1,3   ),W(1,5   ),GG ,AMP(13  ))
      CALL FVOXXX(W(1,3   ),W(1,18  ),GG ,TMASS   ,TWIDTH  ,W(1,20  ))
      CALL IOVXXX(W(1,4   ),W(1,20  ),W(1,5   ),GG ,AMP(14  ))
      CALL VVVXXX(W(1,5   ),W(1,17  ),W(1,18  ),G ,AMP(15  ))
      CALL JGGGXX(W(1,1   ),W(1,2   ),W(1,5   ),G ,W(1,21  ))
      CALL IOVXXX(W(1,4   ),W(1,3   ),W(1,21  ),GG ,AMP(16  ))
      CALL JGGGXX(W(1,5   ),W(1,1   ),W(1,2   ),G ,W(1,22  ))
      CALL IOVXXX(W(1,4   ),W(1,3   ),W(1,22  ),GG ,AMP(17  ))
      CALL JGGGXX(W(1,2   ),W(1,5   ),W(1,1   ),G ,W(1,23  ))
      CALL IOVXXX(W(1,4   ),W(1,3   ),W(1,23  ),GG ,AMP(18  ))
      JAMP(   1) = +AMP(   1)-AMP(   4)-AMP(  11)+AMP(  14)+AMP(  15)
     &             -AMP(  16)+AMP(  17)
      JAMP(   2) = -AMP(   1)+AMP(   2)-AMP(   5)+AMP(  11)-AMP(  12)
     &             -AMP(  17)+AMP(  18)
      JAMP(   3) = -AMP(   2)-AMP(   3)+AMP(  12)+AMP(  13)-AMP(  15)
     &             +AMP(  16)-AMP(  18)
      JAMP(   4) = +AMP(   6)-AMP(   7)-AMP(  10)+AMP(  11)-AMP(  12)
     &             -AMP(  17)+AMP(  18)
      JAMP(   5) = -AMP(   6)-AMP(   8)-AMP(  11)-AMP(  13)+AMP(  15)
     &             -AMP(  16)+AMP(  17)
      JAMP(   6) = +AMP(   7)-AMP(   9)+AMP(  12)-AMP(  14)-AMP(  15)
     &             +AMP(  16)-AMP(  18)
      GG_TTBG = 0.D0
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          GG_TTBG =GG_TTBG+ZTEMP*DCONJG(JAMP(I))/DENOM(I)
      ENDDO
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN)
      END
