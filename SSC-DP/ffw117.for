c
	 subroutine ffw117(p,ei,qt,n1,n2,par,consi,qi,qg,area,delta)
c
c======= The Xinanjiang(2) Model====================
c
c    >>>>inputs<<<<
c    P(I):    Precipitation                           (mm/delta)
c             (delta: time interval of simulation) 
c    EI(I):   Pan evaporation                         (mm/delta)
c    N1:      Starting position for calculation
c    N2:      End position for calculation
c
c    >>>>outputs<<<<
c    QT(I):   Inflow to channel system                (mm/delta)
c
c    >>>>parameters<<<<
c    >>>>parameters<<<<
c    PAR(1)= WM:     Areal mean tension water capacity     (mm)  
c    PAR(2)= X:      Ratio of the WUM to WM            (0 to 1)
c                    WUM: Average basin storage capacity 
c                    of the upper layer
c    PAR(3)= Y:      Ratio of the WLM to (1-X)WM       (0 to 1)
c                    WLM: Average basin storage capacity       
c                    of the lower layer
c                    (WM = WUM + WLM + WDM)
c    PAR(4)= CKE:    Ratio of potential Evp. to pan Evp.
c    PAR(5)= B:      Exponential number of storage capacity 
c                    distribution curve 
c    PAR(6)= SM:     Areal mean free water storage capacity    (mm)
c    PAR(7)= EX:     A parameter in the distribution of free water
c                    storage capacity
c    PAR(8)= CI:     A coefficient relating RI, a contribution to
c                    interflow, to free water storage (areal mean)
c                                                         (1/delta)
c    PAR(9)= CG:     A coefficient relating RG, a contribution to
c                    groundwater, to free water storage (areal mean)
c                                                         (1/delta)
c    PAR(10)= CIMP:  Proportion of impermeable area to the total 
c                    area.
c    PAR(11)= C:     Evapotranspiration coefficient from deep layer
c    PAR(12)= CKI    The interflow recession coefficient  (0 to 1)
c    PAR(13)=CKG:    The groundwater recession coefficient (0 to 1)
c    PAR(14)= CN:    Number of cascade linear reservoir for 
c                    runoff routing
c    PAR(15)= CNK:   Scale parameter of cascade linear     (delta)
c                    reservoir  
c
c    >>>>initial conditions<<<<
c    CONSI(1)= W1    Basin storage at the beginning of the  (mm)
c                    computing time
c    CONSI(2)= WU1:  Upper layer storage at the beginning of (mm)
c                    the computing time
c    CONSI(3)= WL1:  Lower layer storage at the beginning of (mm)
c                    the computing time
c                   (W1 = WU1 + WL1 + WD1) 
c    CONSI(4)= S1:   Free water storage at the begining of   (mm)
c                     the computing time
c    CONSI(5)= QI1:  Interflow discharge at the beginning of(mm/delta)
c                    the computing time   
c    CONSI(6)=QG1:   Groundwater discharge at the beginning  (mm/delta)
c                    of the computing time
c    CONSI(7)=nwarm: length of warm up period
c    CONSI(8)=m:     memory length for UH of the surface runoff    
c
	  implicit none
	  integer NPAR,i,n1,n2
	  PARAMETER(NPAR=3)
        real P(n2), EI(n2+1),QT(n2+1),par(npar),consi(6)
	  real qi(n2+1),qg(n2+1)
!        COMMON /BNDLOW/PARMIN(NPAR),/BNDUP/PARMAX(NPAR)
	  
	  real WM,WUM,WLM,WDM,CKE,B,SM,EX,CI,CG,CIMP,C,CKI,CKG
	  real W1,WU1,WL1,WD1,S1,S2,WWMM,SSM
	  real WU10,WL2,WU2,WD2,W0,WD20,WU20,WL20
	  real EM,PE,A,R,EU,EL,ED,E,W,FR,FR0,FR1,BU,ST,RS,RI,RG,RS1,RS2
	  real area,delta,temp
    
        WM=155.271
	  !WUM=WM*PAR(2) 
		WUM=WM* 0.450
	!write(*,*) WUM,WM,PAR(2)
        !WLM=WM*(1-PAR(2))*PAR(3) 
		WLM=WM*(1- 0.450 )*   0.340
        WDM=WM-WUM-WLM 
!        X=PAR0(2)
!        Y=PAR0(3)
        CKE=par(1)
        B=par(2)
        SM=44.340
	EX=0.582 
	CI=   0.899
	CG=0.100  
        CIMP=0.104
        C=0.118
	CKI=par(3)
	CKG=  0.100
c
	IF((CI+CG).GT.1.0) CG=1.0 - CI
	IF(WDM.LT.0) WDM=0.0
c
        W1=CONSI(1)
        WU1=CONSI(2)
        WL1=CONSI(3)
	S1=CONSI(4)
	QI(n1)=CONSI(5)
c	QG(n2)=CONSI(6)
	QG(n1)=CONSI(6)
c
!        WUM=X*WM
!        WLM=Y*(WM-WUM)
!        WDM=WM-WUM-WLM        
c
c
c       check the initial conditions
c
        
        if(w1.gt.wm) w1=wm
        if(wu1.gt.wum) wu1=wum
        if(wl1.gt.wlm) wl1=wlm
        wd1=w1-wu1-wl1
        if(wd1.gt.wdm) wd1=wdm
	if(s1.gt.sm) s1=sm
c
c       The maximum storage (point) capacity on the basin.        
c	(relations assumed on the permeable area)
c
        WWMM=WM*(1+B)
c
c	The maximum point free water storage capacity on the basin
c
	SSM = SM*(1.0+EX)
c
        DO I=N1,N2
            QT(I)=0.0
            EM=CKE*EI(I)
            PE=P(I)-EM
c
c	Runoff producing area FR
c
c
            IF(PE.GT.0.0) THEN
c
c       Calculating total runoff R
c 
              IF(W1.LT.WM) THEN
                  A=WWMM*(1-(1-W1/WM)**(1.0/(1.+B)))
               ELSE
                  A=WWMM
               ENDIF 
               IF((PE+A).LT.WWMM) THEN
                  R=PE-((WM-W1)-WM*(1-(PE+A)/WWMM)**(1+B))
               ELSE  
                  R=PE-(WM-W1)
               ENDIF
            ELSE  
	       R=0.0
            ENDIF  
c      
c     Start to calculate evaporation
c     EU, EL and ED represent upper, lower and deep layer evaporation
c     respectively.
c  
           EU=0.0
           EL=0.0
           ED=0.0
           WU10=WU1+P(I) 
           IF(WU10.GT.EM) THEN
              EU=EM
           ELSE   
              EU=WU10
              EL=(EM-EU)*WL1/WLM
              IF(EL.LT.(C*(EM-EU))) EL=C*(EM-EU)
              IF(WL1.LT.(C*(EM-EU))) THEN  
                  EL=WL1
                  ED=C*(EM-EU)-EL
              IF(ED.GT.WD1) ED=WD1
              ENDIF
           ENDIF 
c
c     E:  actual evaporation for the time interval
c
           E=EU+EL+ED
c
c      Start to calculate soil moisture storage (tension water):
c      WU2, WL2 and WD2 represent the upper, lower and deep layer
c      soil moisture storage respectively.
c 
          WU20=WU1+P(I)-R-EU
          IF(WU20.GE.WUM) THEN
             WU2=WUM
          ELSE  
             WU2=WU20
          ENDIF  
             TEMP=WU20-WUM
          IF(TEMP.LT.0.0) TEMP=0.0
          WL20=WL1-EL+TEMP
          IF(WL20.GE.WLM) THEN  
             WL2=WLM
          ELSE  
             WL2=WL20
          ENDIF  
             TEMP=WL20-WLM
          IF(TEMP.LT.0.0) TEMP=0.0   
             WD20=WD1-ED+TEMP
          IF(WD20.LT.WDM) THEN  
             WD2=WD20
          ELSE  
c 
             WD2=WDM
          ENDIF  
c
c     W:  soil moisture storage at the end of time interval
c
	  W0=W1
          W=WU2+WL2+WD2
          W1=W
          WU1=WU2
          WL1=WL2
          WD1=WD2
c
c     Runoff producing area at the beginning and end of time interval
c
	  IF(W0.LT.WM) THEN
             FR0 = 1.0 - (1.0-W0/WM)**(B/(1.0+B))
	  ELSE
	     FR0=1.0
	  ENDIF
	  IF(W1.LT.WM) THEN
             FR1 = 1.0 - (1.0-W1/WM)**(B/(1.0+B))
	  ELSE
	     FR1=1.0
	  ENDIF
          FR=(FR0+FR1)/2.0
c
c     Runoff R, seperated into RS, RI and RG
c
          IF(PE.GT.0.0) THEN
	       IF(S1.LT.SM) THEN
	          BU=SSM*(1-(1-S1/SM)**(1/(1+EX)))
	       ELSE
	          BU=SSM
	       ENDIF
	       IF((PE+BU).LT.SSM)  THEN
	          ST=SM-SM*(1-(PE+BU)/SSM)**(1+EX)
	          RS=(PE-ST+S1)*FR
	          RI = ST*CI*FR
	          RG = ST*CG*FR
	          S2=(1-CI-CG)*ST
               ELSE
	          RS=(PE-SM+S1)*FR
	          RI=SM*CI*FR
	          RG=SM*CG*FR
	          S2=(1-CI-CG)*SM
	       ENDIF
            ELSE  
               RS=0.0
	       RI=S1*CI*FR
	       RG=S1*CG*FR
	       S2=(1-CI-CG)*S1
            ENDIF  
c
c     Adjustment for runoff volumes
cc
c     Theoretically, the total runoff R calculated above should 
c     equal to,   rs+ri+rg+(s2-s1)*fr
c     However, due to RS, RI, and RG are functions of runoff producing
c     area FR, and FR is a function of time, therefore, these
c     two values are not equal in general, due to large discrete time
c     steps used.  One way to elliminate
c     these errors is to decrease computing time steps, however, this may
c     increase computing time considerablly. Here, approximate center 
c     difference scheme is adopted, the ratio of 
c        DR = R- (RS+RI+RG)+(S2-S1)*FR
c     to R is within 2% at the Bird creek catchment. The remainder DR 
c     is added to the S2 for the further modification of the runoff volume.
c
	IF(FR.EQ.0.0) THEN
            S2=0.0
        ELSE
            S2=(R-RS-RI-RG+S1*FR)/FR
	ENDIF
	IF(S2.LT.0.0) S2=0.0
	IF(S2.GT.SM) S2=SM
c
	  S1=S2
c
c   Adjustment in considering the impermeable area
c
          IF(PE.GT.0.0) THEN
             RS1=CIMP*PE
             RS2=RS*(1-CIMP)
             RS=RS1+RS2
	     RI=RI*(1-CIMP)
             RG=RG*(1-CIMP)
          ENDIF
          E=E*(1-CIMP)
c
c	RI and RG are routed through corresponding linear reservoirs
c	to obtain QI and QG.  These two flows then added to the RS to 
c 	get total inflow to the channel system, QT

	if(i==365) then
!	write(*,*) i
	endif

c
	QI(i+1)=QI(i)*CKI+RI*(1-CKI)*area/(3.6*delta)
	QG(i+1)=QG(i)*CKG+RG*(1-CKG)*area/(3.6*delta)
	QT(I)=RS
!
c 
	enddo        
c
c       updating the initial conditions at the end of computing 
c       time
c

        consi(1)=W1
        consi(2)=WU1
        consi(3)=WL1
	consi(4)=S1
	consi(5)=qi(n2+1)
	consi(6)=qg(n2+1)
c
        return
        END
c        
