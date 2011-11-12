      SUBROUTINE SWP_EPVE(P1,ANS)
C  
C FUNCTION GENERATED BY MADGRAPH
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C AND HELICITIES
C FOR THE POINT IN PHASE SPACE P(0:3,NEXTERNAL)
C  
C FOR PROCESS : w+ -> e+ ve  
C  
C Crossing   1 is w+ e- -> ve  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NEXTERNAL,   NCOMB,     NCROSS         
      PARAMETER (NEXTERNAL=3, NCOMB= 12, NCROSS=  1)
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
      REAL*8 WP_EPVE
      INTEGER IHEL,IDEN(NCROSS),IC(NEXTERNAL,NCROSS)
      INTEGER IPROC,JC(NEXTERNAL)
      LOGICAL GOODHEL(NCOMB,NCROSS)
      DATA GOODHEL/THEL*.FALSE./
      DATA NTRY/0/
      DATA (NHEL(IHEL,  1),IHEL=1,3) / -1, -1, -1/
      DATA (NHEL(IHEL,  2),IHEL=1,3) / -1, -1,  1/
      DATA (NHEL(IHEL,  3),IHEL=1,3) / -1,  1, -1/
      DATA (NHEL(IHEL,  4),IHEL=1,3) / -1,  1,  1/
      DATA (NHEL(IHEL,  5),IHEL=1,3) /  0, -1, -1/
      DATA (NHEL(IHEL,  6),IHEL=1,3) /  0, -1,  1/
      DATA (NHEL(IHEL,  7),IHEL=1,3) /  0,  1, -1/
      DATA (NHEL(IHEL,  8),IHEL=1,3) /  0,  1,  1/
      DATA (NHEL(IHEL,  9),IHEL=1,3) /  1, -1, -1/
      DATA (NHEL(IHEL, 10),IHEL=1,3) /  1, -1,  1/
      DATA (NHEL(IHEL, 11),IHEL=1,3) /  1,  1, -1/
      DATA (NHEL(IHEL, 12),IHEL=1,3) /  1,  1,  1/
      DATA (  IC(IHEL,  1),IHEL=1,3) /  1,  2,  3/
      DATA (IDEN(IHEL),IHEL=  1,  1) /   3/
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
             T=WP_EPVE(P ,NHEL(1,IHEL),JC(1))            
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
       
       
      REAL*8 FUNCTION WP_EPVE(P,NHEL,IC)
C  
C FUNCTION GENERATED BY MADGRAPH
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : w+ -> e+ ve  
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS,    NEIGEN,    NEXTERNAL       
      PARAMETER (NGRAPHS=   1,NEIGEN=  1,NEXTERNAL=3)   
      INTEGER    NWAVEFUNCS     , NCOLOR
      PARAMETER (NWAVEFUNCS=   3, NCOLOR=   1) 
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
      DATA Denom(1  )/            1/                                       
      DATA (CF(i,1  ),i=1  ,1  ) /     1/                                  
C                                                                          
C ----------
C BEGIN CODE
C ----------
      CALL VXXXXX(P(0,1   ),WMASS ,NHEL(1   ),-1*IC(1   ),W(1,1   ))       
      CALL IXXXXX(P(0,2   ),ZERO ,NHEL(2   ),-1*IC(2   ),W(1,2   ))        
      CALL OXXXXX(P(0,3   ),ZERO ,NHEL(3   ),+1*IC(3   ),W(1,3   ))        
      CALL IOVXXX(W(1,2   ),W(1,3   ),W(1,1   ),GWF ,AMP(1   ))            
      JAMP(   1) = +AMP(   1)
      WP_EPVE = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          WP_EPVE =WP_EPVE+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END
