!*. .................EXAMPLE 2-7 ONE-D UNSTEADY HEAT CONDUCTION..........
!*...................THROUGH A SLAB WITH THIRD KIND OF B.C AT TWO SURFACES........
	SUBROUTINE USER
 	USE VARIABLES
	IMPLICIT NONE
	INTEGER I
	REAL T1

!*----------------SPECI-------FOR SPECIFYING PROBLEM----------
	ENTRY SPECI
	LS=3		!����̬��������        
	MD=1        !ֱ������ϵ
	KI=3        !��ǽ��������߽�����            
	KE=3        !��ǽ��������߽�����  
	A1=6        !��ǽ���滻��ϵ��   
	A2=35       !��ǽ���滻��ϵ��    

	T1=20       !��ǽ�¶ȳ�ֵ     
	T2=-10      !��ǽ�¶ȳ�ֵ     

	AEC=T2*A2   !��ǽ�߽�Դ���е�Sc     
	BE=A2		!��ǽ�߽�Դ���е�Sp 

    AI=T1*A1	!��ǽ�߽�Դ���е�Sc
	BI=A1		!��ǽ�߽�Դ���е�Sp

	DT=10		!ʱ�䲽��
	K=1E10		!��ӡ���Ʊ�����ÿ��K*DT��ӡһ��
	TM=30000    !���ʱ��
	KF=2		!����ӡ�߽�ͨ��
	KP=1		!����¶ȳ�
	RETURN

!*----------------GRID-----FOR GRID SPACING--------
	ENTRY GRID
	XI=0		!��ǽ����߽�����
	XE=.3		!��ǽ����߽�����
	PW=1		!��������
	M1=64		!�ڵ���
	M2=M1-1
	DO I=2,M1
	  XF(I)=(XE-XI)*(REAL(I-2)/(M1-2))**PW+XI		!�����ݻ��Ľ�������
      RF(I)=1										!�����ݻ�����뾶����ֱ������RF=1
	ENDDO
	X(1)=XI
	X(M1)=XE
	DO I=2,M2
	  X(I)=0.5*(XF(I+1)+XF(I))						!�ڵ�����
	  XP(I)=XF(I+1)-X(I)							!�ڵ����ҽ���ľ���
	  XM(I)=X(I)-XF(I)								!�ڵ��������ľ���
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
	LM=0.85           !����ϵ��    
	SC=0              !�ڲ��ڵ�Դ���Sc     
	SP=0              !�ڲ��ڵ�Դ���Sp 
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
