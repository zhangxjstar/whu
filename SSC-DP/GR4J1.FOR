       SUBROUTINE GR4J_MOD1(n,n1,n2,X0,consi,consit,p2,e,Q,npar,ncon,v
     &,Tave)   
	 use numerical_libraries

c >>>>>>>>>>>>>>>>GR4J_MOD(参数X(1)~X(4),状态V,降雨,蒸发,模拟流量)
C      X: Vector of model parameters:      模型参数
C       - X(1): Capacity of the production store (mm) (>=0) 
C       - X(2): Water exchange coefficient (mm) 
C       - X(3): Capacity of the routing store (mm) (>=0) 
C       - X(4): Time base of the unit hydrograph (d) (>=0.5) 
C       - X(NPX+1) to X(NPX+NH): Ordinates of UH1 (calculated in Subroutine UH1)   UH1单位线的横坐标值
C       - X(NPX+NH+1) to X(NPX+3*NH): Ordinates of UH2 (calculated in Subroutine UH2)  UH2单位线的横坐标值
C******************************************************************************************************
C        v: Vector of model states:       模型状态变量

      PARAMETER (NPX=5,NH=10)
!	COMMON /area/area 
      DIMENSION V(3*NH+2),X(NPX+3*NH),X0(NPAR),consi(ncon),consit(ncon)
	DIMENSION P(n2),E(n2),EP(n2),ES(n2),Q(n2),Qobs(n2),Tave(n2),P1(n2)
     &,melt(n2),p2(n2)

	INTEGER n1,n2,I,npar,ncon,n
	REAL PSQ,PSP,QB
C    - V(1): Level in production store 
C    - V(2): level in routing store 
C    - V(3) to V(2+3*NH): Storages for UH1 and UH2 initially set to zero 
C***************************************
C    - X(5) to X(NPX+3*NH): Ordinates of UH1 & UH2



!	X0(1)=   350.5524 
!	X0(2)= 1.136275
!	X0(3)= 20.59231    
!	X0(4)= 3.678994 
!	X0(5)= 0.4004774    
!




	

      DO 3 I=1,NPX
	  X(I)=X0(I)

!	write(*,*) x(i)
    3	CONTINUE



		!测试GR4J参数



       CKE=X0(5)
		DDF=X0(6)


!zxj.由于分段V不能重置

	V(1)=CONSI(1)        !要改为读par文件里常数
	V(2)=CONSI(2)

	snowpack=CONSI(3)

c	DO I=3,2+3*NH 
c         V(I)=0. 
c      ENDDO




C       QB=0.9*PSQ/X(3) 
C 	 V(2)=1.-1./(4.*(1.+QB)**4-3.+(4.*QB)**0.2+(8.*QB)**0.4) 
C        V(2)=V(2)*X(3) 
C        V(1)=X(1)*MIN(1.,SQRT(PSQ/PSP)) 


C Production store: 
      DO 100  I=n1,n2





	CALL UH1(X,X(4))     !得到X(5)~X(14),UH1的横坐标值
      CALL UH2(X,X(4))     !得到X(15)~X(34),UH2的横坐标值

		
		P1(I)=P2(I)
	P(I)=P2(I)
		AA=-0.0
	melt(I)=0.0

	IF (Tave(I).gt.AA) then
		melt(I)=DDF*(Tave(I)-AA)		!T>0 melt
		cc=melt(I);
		IF (melt(I).gt.snowpack) melt(I)=snowpack
		P(I)=P1(I)+melt(I)		!T>0 rain=rain+melt
		snowpack=snowpack-melt(I)
		cc=melt(I);
	else
		P(I)=0		!T<0 snow
		snowpack=snowpack+P1(I)
	endif



	EP(I)=E(I)*CKE*37.586/15.54

      IF(P(I).GE.EP(I))THEN                   
	ES(I)=0.    
	WS=(P(I)-EP(I))/X(1)   
	IF(WS.GT.13)WS=13     
	PS=X(1)*(1.-(V(1)/X(1))**2.)*tanh(WS)/(1.+V(1)/X(1)*tanh(WS))                  
	PR=P(I)-EP(I)-PS   !不需要加perc???                
	ELSE    
      WS=(EP(I)-P(I))/X(1)          
	IF(WS.GT.13)WS=13     
	ES(I)=V(1)*(2.-V(1)/X(1))*tanh(WS)/(1.+(1.-V(1)/X(1))*tanh(WS))                   
	PS=0.               
	PR=0.               
	ENDIF   
      V(1)=V(1)-ES(I)+PS    !不需要减perc???
C Percolation:  
      S2=V(1)/(1+(V(1)/2.25/X(1))**4.)**(0.25)      !公式、V(1)与perc 
	PERC=V(1)-S2                                  !计算的顺序与书上不同 
	V(1)=S2                                      !先产流再渗流？
      PR=PR+PERC 
C UH1: 
      DO 1 K=1,NH-1 
	V(2+K)=V(3+K)+X(NPX+K)*PR    
    1 CONTINUE    
      V(2+NH)=X(NPX+NH)*PR 
C UH2: 
      DO 2 K=1,2*NH-1           
	V(2+NH+K)=V(3+NH+K)+X(NPX+NH+K)*PR        
    2 CONTINUE    
      V(2+3*NH)=X(NPX+3*NH)*PR 
C Water exchange:    
      ECH=X(2)*(V(2)/X(3))**3.5   
C QR calculation (routing store):    
      V(2)=max(0.,V(2)+V(3)*0.9+ECH)       
	R2=V(2)/(1+(V(2)/X(3))**4.)**(0.25)           
	QR=V(2)-R2       
	V(2)=R2 
C QD calculation:    
      QD=MAX(0.,V(3+NH)*0.1+ECH)   
C Total streamflow:    
      Q(I)=QR+QD                
! 	Q(I)=Q(I)*1.0E-3*area*1.0E6/(24*3600)
	
100    CONTINUE 
 
c    updating the initial conditions at the end of computing 
c    time
c
        consit(1)=V(1)
        consit(2)=V(2)
	consit(3)=snowpack
      
        return
        END

 
 
C**************************************************************************************** 
      SUBROUTINE UH1(X,C) 
      PARAMETER (NPX=4,NH=10) 
      DIMENSION X(*) 
      DO 1 I=1,NH 
      X(NPX+I)=SS1(I,C)-SS1(I-1,C) 
    1 CONTINUE 
      END 
 
C**************************************************************************************** 
      FUNCTION SS1(I,C) 
      FI=I 
      IF(FI.LE.0.)THEN 
      SS1=0. 
      RETURN 
      ENDIF 
      IF(FI.LT.C)THEN 
      SS1=(FI/C)**2.5 
      RETURN 
      ENDIF 
      SS1=1. 
      END 
 

 
C**************************************************************************************** 
      SUBROUTINE UH2(X,C) 
      PARAMETER (NPX=4,NH=10) 
      DIMENSION X(*) 
      DO 1 I =1,2*NH 
      X(NPX+NH+I)=SS2(I,C)-SS2(I-1,C) 
    1 CONTINUE 
      END 
 
C**************************************************************************************** 
      FUNCTION SS2(I,C) 
      FI=I 
      IF(FI.LE.0.)THEN 
      SS2=0. 
      RETURN 
      ENDIF 
      IF(FI.LE.C)THEN 
      SS2=0.5*(FI/C)**2.5 
      RETURN 
      ENDIF 
      IF(FI.LT.2.*C)THEN 
      SS2=1.-0.5*(2.-FI/C)**2.5 
      RETURN 
      ENDIF 
      SS2=1. 
      END 
