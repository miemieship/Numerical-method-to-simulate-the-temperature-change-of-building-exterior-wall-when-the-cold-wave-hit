!*. .................EXAMPLE 2-7 ONE-D UNSTEADY HEAT CONDUCTION..........
!*...................THROUGH A SLAB WITH THIRD KIND OF B.C AT TWO SURFACES........
	SUBROUTINE USER
 	USE VARIABLES
	IMPLICIT NONE
	INTEGER I
	REAL T1

!*----------------SPECI-------FOR SPECIFYING PROBLEM----------
	ENTRY SPECI
	LS=3		!非稳态线性问题        
	MD=1        !直角坐标系
	KI=3        !内墙，第三类边界条件            
	KE=3        !外墙，第三类边界条件  
	A1=6        !内墙表面换热系数   
	A2=35       !外墙表面换热系数    

	T1=20       !内墙温度初值     
	T2=-10      !外墙温度初值     

	AEC=T2*A2   !外墙边界源项中的Sc     
	BE=A2		!外墙边界源项中的Sp 

    AI=T1*A1	!内墙边界源项中的Sc
	BI=A1		!内墙边界源项中的Sp

	DT=10		!时间步长
	K=1E10		!打印控制变量，每隔K*DT打印一次
	TM=30000    !最后时刻
	KF=2		!不打印边界通量
	KP=1		!输出温度场
	RETURN

!*----------------GRID-----FOR GRID SPACING--------
	ENTRY GRID
	XI=0		!内墙表面边界坐标
	XE=.3		!外墙表面边界坐标
	PW=1		!均分网格
	M1=64		!节点数
	M2=M1-1
	DO I=2,M1
	  XF(I)=(XE-XI)*(REAL(I-2)/(M1-2))**PW+XI		!控制容积的界面坐标
      RF(I)=1										!控制容积界面半径，对直角坐标RF=1
	ENDDO
	X(1)=XI
	X(M1)=XE
	DO I=2,M2
	  X(I)=0.5*(XF(I+1)+XF(I))						!节点坐标
	  XP(I)=XF(I+1)-X(I)							!节点与右界面的距离
	  XM(I)=X(I)-XF(I)								!节点与左界面的距离
	  R(I)=1
	ENDDO
!*................SPECIFY THE INITIAL FIELD........
	IF(TU.LE.2*DT)THEN 
      DO I=1,M1
	    TG(I)=15.0-REAL(30.0*X(I))/0.85
	  END DO
	END IF
	RETURN

!*----------- GAMSOR----
	ENTRY DIFSOR
	LM=0.85           !导热系数    
	SC=0              !内部节点源项的Sc     
	SP=0              !内部节点源项的Sp 
	RO=1.05E+06                
	DO I=2,M2                   
	  GM(I)=LM
	  CN(I)=SC
	  AP(I)=SP
	  RC(I)=RO
	ENDDO
	GM(1)=LM
	GM(M1)=LM
	RETURN

!*------SUBROUTINE OF INTEROUTPUT-----------
    ENTRY  INTEROUTPUT
	WRITE(*,"(1X,A5,F7.0,A7,F8.4)") "Time=",TU, "T(1)=",T(1)
	RETURN	

!*............SUBROUTINE SPRINT..........
	ENTRY SPRINT
    OPEN(2,FILE='T.DAT',access='append')
	DO I=1,1
	WRITE(2,*) TU,X(I),T(I)
	END DO
	CLOSE(2)
	RETURN
	END
