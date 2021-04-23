!                               一维非稳态导热通用程序(不变部分)
! This is a general purpose program to solve 1-D diffusion.
! problem in the form of:

!      ρcdt/dz=1/a(x)d/dx(a(x)Γdt/dx)+s

!******************2003.7 revised************************
!...................Define Variables.........................
MODULE VARIABLES
INTEGER,PARAMETER::L1=130
REAL,DIMENSION(L1):: X,XF,XM,XP,R,RF,AP
REAL,DIMENSION(L1)::AE,AW,CN,T,TA
REAL,DIMENSION(L1)::TG,GM,RC
INTEGER:: K=1,KM=1,KP=1,OM=1
INTEGER:: JB,JE,KE,KI,KF,KN,KR,KT,LS,MD,M1,M2,NF    
REAL:: AEC,AI,BE,BI,DF,DS,DT,EP,EX,PW,TU,TM,XE,XI
REAL:: A1,A2,T2,TC,SC,SP,RO,TE,DN,LM
END

!..................Main Program...............
PROGRAM MAIN
USE VARIABLES
IMPLICIT NONE
INTEGER I

OPEN(1,FILE="q.dat")
OPEN(2,FILE="temp.dat")
     NF=1
     KN=1
     TU=0
150  KT=1
     CALL Speci   !First to specify the problem
     CALL Grid    !Set up grid points
200  CALL Difsor  !Specify the diff-coeff and source term
220  CALL InterOutput  !Output intermediate results
     CALL Coeff   !Set up coefficients of discretization equation
     CALL TDMA    !Solve the algebraic equation by TDMA
     IF(LS.EQ.2.OR.LS.EQ.4) THEN
        IF(DF.GT.EP) THEN
          DO I=1,M1
            TA(I)=TA(I)+OM*(T(I)-TA(I))
          END DO
          DF=0
          KT=KT+1
          GOTO 200
         END IF
     END IF
	  CALL GPRINT
	  IF(LS.EQ.3.OR.LS.EQ.4) THEN
	    IF(TU.LT.TM) THEN
		  DO I=1,M1
		    TG(I)=T(I)
		  END DO
          KT=1
		  IF(LS.EQ.3) THEN
		    CALL SPRINT
            GOTO 200
          ELSE
		    GOTO 220
		  END IF
		END IF
	  END IF
!.............special results print out, if not, just leave it open..........
	  
	  IF(NF.NE.KM) THEN
	    NF=NF+1
		GOTO 150
	  END IF
CLOSE(2)
CLOSE(1)
END 


!...................Subroutine.................
SUBROUTINE SETUP

USE VARIABLES
REAL,DIMENSION(L1)::P,Q
ENTRY COEFF
!................coefficients of boundary points..........
IF(KI.LE.1) THEN
  AP(1)=1
  AE(1)=0
  AW(1)=0
  CN(1)=AI
ELSE
  AE(1)=GM(1)/XM(2)
  AP(1)=AE(1)+BI
  AW(1)=0
  CN(1)=AI
END IF
IF(KE.LE.1) THEN
  AP(M1)=1
  AE(M1)=0
  AW(M1)=0
  CN(M1)=AEC  
ELSE
  AW(M1)=GM(M1)/XP(M2)
  AP(M1)=AW(M1)+BE
  AE(M1)=0
  CN(M1)=AEC  
END IF
!...................coefficients of internal points........
IF(LS.NE.3.OR.TU.LT.0.5*DT) THEN
  EX=1
  IF(MD.EQ.3)  EX=2
  AW(2)=GM(2)/XM(2)*RF(2)**EX
  AE(M2)=GM(M2)/XP(M2)*RF(M1)**EX
END IF
DO I=2,M2-1
  AE(I)=RF(I+1)**EX/(XP(I)/GM(I)+XM(I+1)/GM(I+1))
  AW(I+1)=AE(I)
END DO
DO I=2,M2
  AP(I)=AE(I)+AW(I)-AP(I)*(XF(I+1)-XF(I))*R(I)**EX
  CN(I)=CN(I)*(XF(I+1)-XF(I))*R(I)**EX
  IF(LS.EQ.3.OR.LS.EQ.4) THEN
    AP(I)=AP(I)+RC(I)*(XF(I+1)-XF(I))*R(I)**EX/DT
    CN(I)=CN(I)+RC(I)*(XF(I+1)-XF(I))*R(I)**EX/DT*TG(I)
  END IF  
END DO
RETURN

ENTRY TDMA
!....................elimination.............
P(1)=AE(1)/AP(1)
Q(1)=CN(1)/AP(1)
DO I=2,M1
  P(I)=AE(I)/(AP(I)-AW(I)*P(I-1))
  Q(I)=(CN(I)+AW(I)*Q(I-1))/(AP(I)-AW(I)*P(I-1))
END DO
!..................back substitution..........
T(M1)=Q(M1)
DO I=M2,1,-1    
  T(I)=P(I)*T(I+1)+Q(I)
END DO
IF(LS.EQ.2.OR.LS.EQ.4) THEN
  DO I=1,M1
    DS=ABS(T(I)-TA(I))
    IF(T(I).GT.1.E-20) DS=DS/T(I)
	IF(DF.LT.DS) DF=DS
  END DO
END IF
RETURN

ENTRY GPRINT
IF(LS.EQ.3.OR.LS.EQ.4) THEN
  M=(TU+0.5*DT)/(K*DT)
  IF(M.NE.KN) THEN
    TU=TU+DT
	RETURN
  END IF
END IF
!...................surface flux calculation..........
SELECT CASE(KI)
  CASE(1) 
    QI=GM(1)*(T(1)-T(2))/XM(2)
  CASE(2)
    QI=AI
  CASE(3)
    QI=AI-BI*T(1)
END SELECT
SELECT CASE(KE)
  CASE(1)
    QE=GM(M1)*(T(M1)-T(M2))/XP(M2)
  CASE(2)
    QE=AEC
  CASE(3)
    QE=AEC-BE*T(M1)
END SELECT
IF(LS.EQ.1.OR.LS.EQ.2) THEN
  S=0
  DO I=2,M2
    QE=QE*RF(M1)**EX
	AP(I)=AE(I)+AW(I)-AP(I)
	S=S+CN(I)+AP(I)*T(I)
  END DO
END IF
!................now it is ready to print out.........
KN=KN+1
IF(KP.NE.2) THEN
  WRITE(*,*) "Dependent Variables Distribution"
  SELECT CASE(MD)
    CASE(1)
      WRITE(*,*) "Cartisian Coordinates"
    CASE(2)
      WRITE(*,*) "Cylindrical Coordinates"
    CASE(3)
      WRITE(*,*) "Spherical Coordinates"
    CASE(4)
      WRITE(*,*) "Nonuniform cross section"
  END SELECT
  SELECT CASE(LS)
    CASE(1)
      WRITE(*,*) "Linear Steady Problem"
    CASE(2)
      WRITE(*,*) "Nonlinear Steady Problem"
      WRITE(*,*) "Iterative Times=", KT
    CASE(3)
      WRITE(*,*) "Linear Unsteady Problem"
    CASE(4)
      WRITE(*,*) "Nonlinear Unsteady Problem"
  	  WRITE(*,*) "Iterative Times=", KT
	  WRITE(*,*) "At Time=",T
  END SELECT
  JE=0
  DO WHILE(JE.LT.M1)
    JB=JE+1
    JE=JE+4
    IF(JE.GT.M1) JE=M1
    WRITE(*,*) "J"
    DO J=JB,JE
      WRITE(*,*) J
    END DO
    IF(MD.EQ.2.OR.MD.EQ.3) THEN
      WRITE(*,*) "R"
    ELSE
      WRITE(*,*) "X"
    END IF
    DO J=JB,JE
      WRITE(*,*) T(J)
    END DO
  END DO 
END IF
QI=QI*RF(2)**EX
QE=QE*RF(M1)**EX
IF(KF.EQ.1) THEN
  WRITE(*,*) "Total Heat Flow At Int.Surface Qi=",QI
  WRITE(*,*) "Total Heat Flow At Ext.Surface Qe=",QE
  IF(LS.EQ.1.OR.LS.EQ.2) THEN
    WRITE(*,*) "Total heat Input Form Source Term S=",S
	WRITE(*,*) "Heat Balance:Qi+Qe+ys=",QI+QE+S
  END IF
END IF
TU=TU+DT
RETURN
END

