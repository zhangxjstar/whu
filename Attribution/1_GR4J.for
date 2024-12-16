c       Program UCGXIN2
c      ===============
c      Calibration and verification of the Xinanjiang-2 model 
c
c      Version no.2:  the difference between V.2 and V.1 lies
c                     in that in this version, initial conditions
c                     kept constant for both model optimisation
c                     and model simulation processes, while	   
c                     in V.1, the initial conditions for the 
c                     next run are the updated values from
c                     the previous run. On most of the catchments,
c                     this change will not make much difference,
c                     however, in some cases, optimisations are
c                     very sensitive to initial conditions, and
c                     as a consequence, the V.1 might not converge.
c
c      Note: in the identification process, the simulated discharges
c            will not be filed out.  If you want to get the simulated
c            dischargs, as well as MSE for different zones, then
c            you have to run the program again with the optimum
c            model parameters without any further optimisation
c            required.
c
c      Calling subroutines:
c      FFW1  FFW2  FFW17  FFW42  FFW117  FFW118  FFW119  FFW120 XINOB2 
c
       parameter (idnx=80000,npar=7,ncon=4,nosub=438)
       dimension rain(idnx),evap(idnx),qobs(idnx),eobs(idnx),qhat(idnx),
     &           ehat(idnx),year(idnx),month(idnx),day(idnx),Tave(idnx),
     &		   tmin(idnx),tmax(idnx),pres(idnx),radi(idnx),xLAI(idnx)
     &,slope(idnx),corr(17),posi(17),a1(17,idnx),a2(idnx)
       dimension dpar(npar)
	 dimension xnsec(2000),xnsev(2000),SRM011(2000),SRM11(2000)
     &,SRM021(2000),SRM21(2000),SRM22(2000),SRM022(2000)
       dimension par(npar),parmin(npar),parmax(npar),qi(idnx),qg(idnx),
     &       ititle(16),consi(ncon),h(8000),consi0(ncon),data_sum(4000),
     &factor_sum(4000),a1min(1),a1max(1),xnseh(2000)
       dimension X1(idnx),X2(idnx),X3(idnx),X4(idnx),CKE(idnx)
	 dimension rsqr11(idnx),rsqr12(idnx),rsqr21(idnx),rsqr22(idnx),
     &           ratio1(idnx),ratio2(idnx),RR1(idnx),RR2(idnx),Xc(6),
     &           RIOA11(idnx),RIOA12(idnx),RIOA21(idnx),RIOA22(idnx),
     &		   rsqr31(idnx),rsqr32(idnx),ratio3(idnx),RR3(idnx),
     &		   RIOA31(idnx),RIOA32(idnx),qsimc(idnx)
       character*80 catchm,fyle,filename,fyout,fyresult,fydata,
     &          fyle_dc,fyle_wbi,fyle_ioa,fyle_par,fydata2,fydata3,
     &filename1,filename2,filename3,filename4,filename5,filename6,
     &    fyf1,fyf2,fyf3,fyf4,fyf5,fyresult1,fyf6

       integer i_file,i_x,i_x1,iyea_in(idnx),iyea_fa(idnx)
	integer ifa1,ifg1,ite,itest1,nx11,nx21,nx31,isn,ism1,ism2

       common /syspar/iuniti,iunito,idebug
       common /rain/rain,/evap/evap,/disc/qobs,/actu/eobs,
     &    /qhat/qhat,/ehat/ehat,/consi/consi,/uh/h
       common /number/nc1,nc2,nwarm,m, /switch/isw,/area/area
       common /bndlow/pmin(npar),/bndup/pmax(npar),/Tave/Tave
	common /posi/posi
	common /ifa1/ifa1,/iyea_in/iyea_in,/a1/a1,/Xc/Xc
	common /dpar/dpar


       external xinob2


       data iuniti/5/,iunito/6/,idebug/0/
		character*7 filenamex(3000),x

		character(len=80):: basin



c
c  Inputs to the Xinanjiang(2) model
c  ------------------------------
c  RAIN:   contains rainfall data time series on entry and unchanged 
c          on exit.
c  EVAP:   contains evaporation data time sries on entry and unchanged
c          on exit.
c  QOBS:   contains observed discharge time series on entry and 
c          unchanged on exit.
c  EOBS:  contains actual evaporation data time sries on entry and unchanged
c          on exit.
c  QHAT:   contains estimated discharges on exit.
c  PAR:    contains initial Xinanjiang(2) model parameters on entry; 
c          on exit it will contain optimized model parameters.
c  PARMIN: the lower bound of the  model parameters on entry and 
c          unchanged on exit.
c  PARMAX: the upper bound of the  model parameters on entry and
c          unchanged on exit.
c  CONSI:  contains initial conditions of the catchment as detailed 
c          later
c 
c    >>>>parameters<<<<
c    PAR(1)= WM:     Areal mean tension water capacity     (mm)  
c    PAR(2)= X:      Ratio of the WUM to WM            (0 to 1)
c                    WUM: Average basin storage capacity 
c                    of the upper layer
c    PAR(3)= Y:      Ratio of the WLM to (1-X)WM       (0 to 1)
c                    WLM: Average basin storage capacity       
c                    of the lower layer
c                    (WM = WUM + WLM + WDM)
c    PAR(4)= KE:     Ratio of potential Evp. to pan Evp.
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
c    CONSI(1)= W1    Arel tension water storage at the beginning (mm)
c                    of the computing time
c    CONSI(2)= WU1:  Upper layer storage at the beginning of (mm)
c                    the computing time
c    CONSI(3)= WL1:  Lower layer storage at the beginning of (mm)
c                    the computing time
c                   (W1 = WU1 + WL1 + WD1) 
c    CONSI(4)= S1:   Areal free water storage at the beginning (mm)
c                    of the computing time
c    CONSI(5)= QI1:  Interflow discharge at the beginning (mm/delta)
c                    of the computing time   
c    CONSI(6)=QG1:   Groundwater discharge at the beginning  (mm/delta)
c                    of the computing time
c    CONSI(7)=nwarm: length of warm up period
c    CONSI(8)=m:     memory length for UH of the surface runoff    
c
c
c        number of stations or catchments.
c	   nosub=11

	open(unit=103,file="LIST.txt"
     &	,status='old')

	do i_file=1,361

	read(103,*) filenamex(i_file)

c	write(*,*) filenamex(i_file)


	enddo







       write(iunito,150)
150    format(//4x,'Name of the catchment :? ',t50,'>',$)       
c       read(iuniti,160)catchm
      
160      format(a20)
         write(iunito,170)
170      format(/4x,'Name of the general control file '/
     &    4x,'(Example is given in UCGXIN2.PAR) :? ',t50,'>',$)
c        read(iuniti,160)fyle



	!遍历不同影响因素




!	open(unit=105,file="nse.txt",status='old')


	   fyle='input.txt'
         
!         	do i_file=1,1767

!	read(105,*) xnsec(i_file),xnsev(i_file),xnseh(i_file)


!	enddo

	             
	open(unit=416,file="factor_d.txt",status='old')

	open(unit=417,file="para_d.txt",status='old')

	open(unit=418,file="runoff_d.txt")

		open(unit=419,file="rep.txt")

       do i_file=1,361

	ite=0
	itest1=0
	write (fyf6,'("../sim2/",(a7),"sim.txt")') filenamex(i_file)

	open(unit=104,file=fyf6,status='old')
		
      read(104,*) area,NX

!	do i=1,NX
	
!	read(104,*) ix,xx,qsimc(i)
!	enddo

	close(104)


      

	write (fydata,'((a7), "srm.txt")') filenamex(i_file)
	 write (fyout,'((i4), "output.txt")') i_file
	 write (fyresult,'((i4), "qet.dat")') i_file
		write (fydata2,'((a7), "m_srm.txt")') filenamex(i_file)
		write (fydata3,'((a7), "recal.txt")') filenamex(i_file)

		write (fyresult1,'((a7), "rep.txt")') filenamex(i_file)


	write (fyf1,'("../fac2/",(a7), "factor1.txt")') filenamex(i_file)
	write (fyf2,'("../fac2/",(a7), "factor2.txt")') filenamex(i_file)
	write (fyf3,'("../fac2/",(a7), "factor3.txt")') filenamex(i_file)
	write (fyf4,'("../fac2/",(a7), "factor4.txt")') filenamex(i_file)

!		open(8,file=fydata,status='unknown')
!		open(33,file=fydata2)


		open(371,file=fyresult1)

!	write (fyf5,'("../i1/",(a7),"par.txt")') filenamex(i_file)


!	open(unit=315,file=fyf5,status='old')

!	do i=1,6

!	read(315,*) Xc(i)

!	enddo
	
!	write(*,*) Xc(1:6)


	write (filename5,'("../co2/",(a7), "p2.txt")') filenamex(i_file)
	write (filename6,'("../co2/",(a7), "s2.txt")') filenamex(i_file)

	open(unit=310,file=filename6,status='old')

	read(310,*) corr(1:17)

	close(310)

	open(unit=311,file=filename5,status='old')

	read(311,*) posi(1:17)

	close(311)	


	read(416,*) ifa1

	read(417,*) dpar(1:7)


		       


!	do ifa1=1,17

!	if (abs(corr(ifa1))/=0) then
	
!	itest1=itest1+1

	a1=0


	!根据因子选择不同的时期的降雨蒸发径流资料


	write(filename1,'("../inp2/",(a7), "input1.txt")') filenamex(i_file)
	write(filename2,'("../inp2/",(a7), "input2.txt")') filenamex(i_file)
	write(filename3,'("../inp2/",(a7), "input3.txt")') filenamex(i_file)
	write(filename4,'("../inp2/",(a7), "input4.txt")') filenamex(i_file)
	
		if ((ifa1>=1).and.(ifa1<=9)) then

	open(unit=306,file=filename1,status='old')
		
		i=0

	do while(.true.)

	i=i+1

	READ(306,*,end=301) iyea_in(i), rain(i),evap(i),tmin(I),tmax(I)
     &,qobs(I),pres(I),radi(I),Tave(I),xLAI(I),slope(I),qsimc(i)

!	write(*,*) iyea_in(i), rain(i),evap(i),tmin(I),tmax(I)
!     &,qobs(I),pres(I),radi(I),Tave(I),xLAI(I),slope(I)


	enddo



301	Close (unit=306)


	data_sum(i_file) = i-1

	write(*,*) i_file

	open(35,file=fyf1)	

		i=0

	do while(.true.)

	i=i+1

	READ(35,*,end=305) iyea_fa(i), a1(1:4,iyea_fa(i))
     &,a2(i),a1(5:9,iyea_fa(i))

!	write(*,*) iyea_fa(i), a1(1:4,iyea_fa(i))
!     &,a2(i),a1(5:9,iyea_fa(i))

	enddo

305	Close (unit=35)	

	factor_sum(i_file) = i-1


	endif



	if ((ifa1>=10).and.(ifa1<=12)) then

	open(unit=307,file=filename2,status='old')
		
		i=0

	do while(.true.)

	i=i+1

	READ(307,*,end=302)  iyea_in(i), rain(i),evap(i),tmin(I),tmax(I)
     &,qobs(I),pres(I),radi(I),Tave(I),xLAI(I),slope(I),qsimc(i)


	enddo



302	Close (unit=307)

	data_sum(i_file) = i-1

	write(*,*) data_sum(i_file)

	open(36,file=fyf2)
		
		i=0

	do while(.true.)

	i=i+1

	READ(36,*,end=306) iyea_fa(i), a1(10:12,iyea_fa(i))

	enddo


306	Close (unit=36)



	factor_sum(i_file) = i-1




	endif







	if ((ifa1>=13).and.(ifa1<=15)) then

	open(unit=308,file=filename3,status='old')
		
		i=0

	do while(.true.)

	i=i+1

	READ(308,*,end=303)  iyea_in(i), 	rain(i),evap(i),tmin(I),tmax(I)
     &,qobs(I),pres(I),radi(I),Tave(I),xLAI(I),slope(I),qsimc(i)


	enddo



303	Close (unit=308)

	data_sum(i_file) = i-1

	write(*,*) data_sum(i_file)

	open(37,file=fyf3)
		
		i=0

	do while(.true.)

	i=i+1

	READ(37,*,end=307) iyea_fa(i), a1(13:15,iyea_fa(i))

	enddo


307	Close (unit=37)



	factor_sum(i_file) = i-1



	endif




	if ((ifa1>=16).and.(ifa1<=17)) then

	open(unit=309,file=filename4,status='old')
		
		i=0

	do while(.true.)

	i=i+1

	READ(309,*,end=304)  iyea_in(i), rain(i),evap(i),tmin(I),tmax(I)
     &,qobs(I),pres(I),radi(I),Tave(I),xLAI(I),slope(I),qsimc(i)


	enddo



304	Close (unit=309)

	data_sum(i_file) = i-1

	write(*,*) data_sum(i_file)

	open(38,file=fyf4)
		
		i=0

	do while(.true.)

	i=i+1

	READ(38,*,end=308) iyea_fa(i), a1(16:17,iyea_fa(i))

	enddo


308	Close (unit=38)




	factor_sum(i_file) = i-1


	endif





!	do ifg1=1,2

!	do ifg2=1,2

!	write(*,*) ifa1,ifg1





         open(unit=20,file=fyle,status='old')
         read(20,175)ititle
175      format(16a4)
c         read(20,*)area
c         read(20,160)fyler
c         read(20,*)nc1
c         read(20,*)nc2
c        read(20,*)nv1
c         read(20,*)nv2
         read(20,175)zzz
         read(20,175)zzz
       do 1 i=1,npar
         read(20,*)par(i),parmin(i),parmax(i)
         pmin(i)=parmin(i)
         pmax(i)=parmax(i)
         if(parmin(i).eq.parmax(i)) par(i)=parmax(i)
1      continue
         read(20,175)zzz
         read(20,175)zzz
         read(20,175)zzz
         read(20,*)(consi(i),i=1,ncon)
         nwarm=consi(3)
c         m=consi(8)

!xiugai3:parameter scale

!		parmin(3)=0.0
!		parmax(3)=(0.0+ifg1*1.0/2.0)*parmax(3)

!	 pmin(3)=parmin(3)
!         pmax(3)=parmax(3)
	   	
!		parmin(7)=0.0
!		parmax(7)=(0.0+ifg2*1.0/2.0)*parmax(7)

!	 pmin(7)=parmin(7)
 !        pmax(7)=parmax(7)

c
c      Read parameters for genetic optimization
c
        read(20,175)zzz
        read(20,175)zzz
        read(20,*)nvar,npop,nstr1,neva,pm,cmax,iseed
c
c      Read information for Rosenbrock's optimisation
c
        read(20,175)zzz
        read(20,175)zzz
        read(20,*)ioprsb,tolx,tolf,itmxrs,iprt
c
c       Read information for Simplex optimisation
c
        read(20,175)zzz
        read(20,175)zzz
        read(20,*)iopsmp,atol,tolf,itmxsp,iprt
       close(unit=20)
c
       do 2 i=1,ncon
          consi0(i)=consi(i)
2      continue

c
	Close (unit=20)
c

c       i_file=1



	 do 9 i=1,npar
	   par(i)=(parmax(i)-parmin(i))*rand()+parmin(i)
9      continue
       X1(i_file)=0.0
       X2(i_file)=0.0
       X3(i_file)=0.0
       X4(i_file)=0.0
       CKE(i_file)=0.0


c      read precipitation,evaporation,discharge and actual evaporation data
c
c       ier=0
c      call ffw1(rain,evap,qobs,eobs,nx,area,idnx,filename,ititle,iyear0,
c     &      imon0,iday0,ih0,imin0,isec0,idt,l29,icode,iunit,iprint,ier)
c	Subroutine ffw1(x1,x2,x3,x4,nx,area,idimx,fyle,ititle,iyear,imon,iday,ihour,
c    1      imin,isec,idt,l29,icode,iunit,nprint,ier)
      






c      DO 250 i=1,NX
c	   READ(30,400)year(i),month(i),day(i),rain(i),evap(i),
c     &             qobs(i),eobs(i)
c400      format(F7.0,i4,F4.0,f6.2,3f8.4)
c400        format(F7.0,i6,F7.0,4f10.4)
c
c250   CONTINUE

c	Close (unit=30)

	
        nc1=1
	 nc2=data_sum(i_file)
	 nv1=data_sum(i_file)*2/3+1
	 nv2=data_sum(i_file)
       nn=nc2
       if(nn.lt.nv2) nn=nv2


c       if(ier.ne.0) then
c         write(iunito,180)ier    
c180       format(/4x,'Error detected in FFW1: ier=',i2)
c         stop
c       endif
c       if(nn.gt.nx) then
c       write(iunito,185)nx,nn
c185     format(/4x,'Observed data series nx= ',i4, ' but the '
c     &       'required data length nn= ',i4, ' STOP!')
c       stop
c       endif
c
c  
       idt=3600
	 delta=idt/3600
       write(iunito,190)   i_file
190    format(/4x,'File for general model output :? ',i3)
c       read(iuniti,160)fyout
c
       open(14,file=fyout,status='unknown')

       write(14,5) catchm
5      format(//4x,'Xinanjiang(2) model for catchment ',t50,a20//)
       write(14,6) nc1,nc2
6      format(4x,'Calibration period starts at: ',i4,
     &    ' and finishes at: ',i4)
       write(14,7) nv1,nv2
7      format(4x,'Verification period starts at: ',i4,
     &      ' and finishes at: ',i4)
       write(14,8) nwarm
8      format(4x,'Number of days of warm up period: ',i4)
c
c      Changing the discharge unit in to mms/delta, here delta= 1 day
c
c       if(iunitq.eq.2) then
c          do 25 i=1,nn
c             qobs(i)=qobs(i)*3.6/76264.8.
c25        continue
c        endif
c
c
      qbar=0.0
      do 35 i=nc1,nc2
         qbar=qbar+qobs(i)
35    continue
      qbar=qbar/(nc2-nc1+1)
      var0=0.0
      do 36 i=nc1,nc2
         var0=var0+(qobs(i)-qbar)**2
36    continue
      var0=var0/(nc2-nc1+1)
c  
       write(14,40)
40     format(//4x,'The starting parameters for Optimisation')
       write(14,100) (par(i),i=1,npar)
c       write(14,101)(par(i),i=9,15)
c

       write(14,45)
45     format(/4x,'The initial conditions used in calculations')
       write(14,46)(consi(i),i=1,ncon)
46     format(/4x,9x,'V(1)',8x,'V(2)',8X,'nwarm',/4X,3F10.3)
c
      isw=0
      it=0
      if(neva.eq.0) goto 1000
       call ffw120(par,fb,xinob2,nvar,npop,nstr1,neva,pm,cmax,
     1                parmin,parmax,iseed,it)

c
       write(14,80)
80     format(//4x,'Results of optimization by the genetic algorithm:')
1000   write(14,85)it
85     format(/4x,'Total number of function evaluations: ',i5)
       write(14,90)
90     format(/4x,'Values of parameters:')
       write(14,100) (par(i),i=1,npar)
100    format(/4x,7X,'X(1)',8X,'X(2)',8X,'X(3)',7X,'X(4)',7X,'CKE',
     &         /4X,5F9.3)
c       write(14,101)(par(i),i=9,15)
c101    format(4x,7x,'CG',6x,'IMP',8x,'C',7x,'KI',7x,'KG',8X,'N',7X,'NK'
c     &        /4x,7f9.3/)
c
       nc10=nc1
       nc20=nc2
       n1=nc1
       if(n1.gt.nv1) n1=nv1
       n2=nv2
       if(n2.lt.nc2) n2=nc2
c
c      Using XINOB2 to calculate qhat over the whole period
c
c      Initializing starting condition for model simulation
c
       do 140 i=1,ncon
140    consi(i) = consi0(i)
c
       nc1=n1
       nc2=n2
      isw=1
C 用优化得到的一套参数run模型率定期+检验期，再分别输出率定期和检验期的计算指标
          FF=XINOB2(PAR)     
      isw=0
       nc1=nc10
       nc2=nc20
c
c       write(14,103)m
c103    format(/4x,'Unit hydrograph (surface runoff, m=',i3,'):')
c       write(14,105)(h(i),i=1,m)
c105    format(4x,10f7.4)

c
       qm=0
	 em=0
       do 30 i=nc1,nc2
          qm=qm+qobs(i)
	    em=em+eobs(i)
30     continue
       qm=qm/(nc2-nc1+1)
	 em=em/(nc2-nc1+1)
c
       write(14,110) 
110    format(///4x,'Results in the calibration period'/)
c
        n1=nc1+nwarm
        itemp=iunito
        iunito=14
       call ffw17(qobs,eobs,n1,nc2,qhat,ehat,qm,em,f1,fhat1,fagr1,
     & fagr11,f2,fhat2,fagr2,fagr22,f3,f4,rsqr11(i_file),rsqr12(i_file),
     & 1,0,ratio1(i_file),RR1(i_file),RIOA11(i_file),RIOA12(i_file))
c	  subroutine ffw17(x,y,n1,n2,xhat,yhat,xbar,ybar,f1,fhat1,f2,fhat2,
C     &	       f3,f4,rsqr1,rsqr2,nprint,ier,ratio,RR)
c
       write(14,120)
120    format(///4x,'Results in the verification period'/)        
       n1=nv1
       if(n1.le.(nc1+nwarm)) n1=nc1+nwarm
       call ffw17(qobs,eobs,n1,nv2,qhat,ehat,qm,em,f1,fhat1,fagr1,
     & fagr11,f2,fhat2,fagr2,fagr22,f3,f4,rsqr21(i_file),rsqr22(i_file),
     & 1,0,ratio2(i_file),RR2(i_file),RIOA21(i_file),RIOA22(i_file))
       
       write(14,130)
130    format(///4x,'Results in the whole period'/) 
	 
	 
	 n1=nc1+nwarm  

	 call ffw17(qobs,eobs,n1,nv1-1,qhat,ehat,qm,em,f1,fhat1,fagr1,
     & fagr11,f2,fhat2,fagr2,fagr22,f3,f4,rsqr31(i_file),rsqr32(i_file),
     & 1,0,ratio3(i_file),RR3(i_file),RIOA31(i_file),RIOA32(i_file))
	n1=nv1
       if(n1.le.(nc1+nwarm)) n1=nc1+nwarm
c 
	iunito=itemp

       write(14,123)
123    format(/4x,'The final conditions at the end of calculations')
       write(14,46)(consi(i),i=1,ncon)
c
c      Initializing starting conditions for optimisation
c
       do 124 i=1,ncon
         consi(i)=consi0(i)
124    continue
c
       if((ioprsb.eq.1).or.(iopsmp.eq.1)) then
c
          if(ioprsb.eq.1) then
             write(14,220)
220          format(//4x,'Results of further tunning by Rosenbrock'
     1          ' optimisation')
             tolf=tolf*var0
             it=0
             call ffw119(nvar,par,parmin,parmax,xinob2,tolx,tolf,
     1               itmxrs,iprt,it)
             tolf=tolf/var0
             ioprs0=ioprsb
             ioprsb=0
             goto 1000
          endif
          if(iopsmp.eq.1) then
             write(14,225)
225          format(//4x,'Results of further tunning by Simplex '
     1              'optimisation')
             tolf=tolf*var0
             it=0
             a=1.0
             call ffw118(par,parmin,parmax,f0,xinob2,nvar,a,atol,
     1            tolf,it,itmxsp,iprt)
             tolf=tolf/var0
             iopsm0=iopsmp
             iopsmp=0
             goto 1000
          endif
       endif

!保存最后优化的一套参数
      X1(i_file)= par(1)
      X2(i_file)= par(2)
      X3(i_file)= par(3)



    

c	open(8,file=fydata,status='unknown')
c	do i=1,n2
c	write(8,*)i,qobs(i),qhat(i),eobs(i),ehat(i)
c	enddo
c	Close (unit=8) 

	sum5=0.0
	sum6=0.0
	sum7=0.0
	sum8=0.0

	isnx=INT(data_sum(i_file)/365)



	ss1=0.0
	ss2=0.0
	ss3=0.0

	ss01=0.0
	ss02=0.0
	ss03=0.0

	sumr1=0.0
	sumr2=0.0



	do isn=2,isnx

		sum1=0.0
		sum2=0.0

		sum3=0.0
		sum4=0.0

		sum9=0.0
		sum10=0.0
		sum11=0.0
		sum12=0.0
		
		sum13=0.0
		sum14=0.0
		sum15=0.0
		sum16=0.0

		sum17=0.0
		sum18=0.0
		sum19=0.0

		sum20=0.0
		sum21=0.0
		sum22=0.0


	ism1=0
	ism2=0


	do i =366, data_sum(i_file)

		if (i<(isn-1)*365.or.i>isn*365) then

		sum1=sum1+(qobs(i)-qhat(i))**2
		sum9=sum9+qobs(i)
		sum17=sum17+(qobs(i)-qsimc(i))**2

	ism1=ism1+1

		else

		sum3=sum3+(qobs(i)-qhat(i))**2
		sum10=sum10+qobs(i)
		sum18=sum18+(qobs(i)-qsimc(i))**2

	ism2=ism2+1

		endif

		sum13=sum13+(qobs(i)-qhat(i))**2
		sum14=sum14+qobs(i)
		sum19=sum19+(qobs(i)-qsimc(i))**2

	

	enddo

	sum9=sum9/ism1

	sum10=sum10/ism2

	sum14=sum14/(data_sum(i_file)-366+1)

	do i =366, data_sum(i_file)

		if (i<(isn-1)*365.or.i>isn*365) then
	
		sum11=sum11+(qobs(i)-sum9)**2

		else

		sum12=sum12+(qobs(i)-sum10)**2

		endif
		
		sum15=sum15+(qobs(i)-sum14)**2

	enddo

		rsqr31(i_file)=1.0-sum1/sum11
		rsqr21(i_file)=1.0-sum3/sum12
		rsqr11(i_file)=1.0-sum13/sum15

	ss1=ss1+rsqr31(i_file)
	ss2=ss2+rsqr21(i_file)
	ss3=ss3+rsqr11(i_file)

		xnsec(i_file)=1.0-sum17/sum11
		xnsev(i_file)=1.0-sum18/sum12
		xnseh(i_file)=1.0-sum19/sum15

	ss01=ss01+xnsec(i_file)
	ss02=ss02+xnsev(i_file)
	ss03=ss03+xnseh(i_file)


!	write(*,*) 	isnx,isn,xnseh(i_file)

	enddo


		do i =366, data_sum(i_file)

		sumr1=sumr1+qhat(i)
		sumr2=sumr2+qsimc(i)

		enddo

	RE001=(sumr1-sumr2)/factor_sum(i_file)

	write(*,*) factor_sum(i_file)

	write(418,*) RE001

	ss1=ss1/(isnx-1)
	ss2=ss2/(isnx-1)
	ss3=ss3/(isnx-1)

	ss01=ss01/(isnx-1)
	ss02=ss02/(isnx-1)
	ss03=ss03/(isnx-1)

!	write(*,*) ss1,ss2,ss3,ss01,ss02,ss03
		
		nx11=data_sum(i_file)*2/3
		nx21=data_sum(i_file)-nx11
		nx31=data_sum(i_file)

!		write(*,*)  nx11,nx21


	srm011(i_file)= (1-ss01)*1.0/(1.0-
     &(2.0/nx11-2.0/nx11*log(2.0/nx11)+log(nx11*1.0)/(2.0*nx11))**0.5)
	srm22(i_file) = (1-ss3)*1.0/(1.0-
     &(3.0/nx31-3.0/nx31*log(3.0/nx31)+log(nx31*1.0)/(2.0*nx31))**0.5)

	srm021(i_file)= (1-ss02)*1.0/(1.0-
     &(2.0/nx21-2.0/nx21*log(2.0/nx21)+log(nx21*1.0)/(2.0*nx21))**0.5)
	srm21(i_file) = (1-ss2*1.0)/(1.0-
     &(3.0/nx21-3.0/nx21*log(3.0/nx21)+log(nx21*1.0)/(2.0*nx21))**0.5)

!	Qave_c=sum(qobs(nc1:nc2))/(nc2-nc1+1)
!	Qave_v=sum(qobs(nv1:nv2))/(nv2-nv1+1)
!	Qave_w=sum(qobs(nc1:nv2))/(nv2-nc1+1)


!	do i=nc1,nc2



	srm11(i_file)=(1-ss1)*1.0/(1.0-
     &(3.0/nx11-3.0/nx11*log(3.0/nx11)+log(nx11*1.0)/(2.0*nx11))**0.5)
	srm022(i_file)=(1-ss03)*1.0/(1.0-
     &(2.0/nx31-2.0/nx31*log(2.0/nx31)+log(nx31*1.0)/(2.0*nx31))**0.5)

	write(419,'(2i,19f)') i_file,ifa1,ss1
     &,ss2,ss3,ss01,ss02,ss03,
     &srm11(i_file),srm21(i_file),srm22(i_file),srm011(i_file),
     &srm021(i_file),srm022(i_file), par(1:7)

!	write(*,*)  xnsec(i_file),rsqr11(i_file)
!     &,xnsev(i_file),rsqr21(i_file)

!	write(*,*)  srm011(i_file),srm021(i_file)
!     &,srm11(i_file),srm21(i_file)


!	if(itest1==0) then

!	temp_srm=srm022(i_file)

!	write(*,*) temp_srm

!	endif

!	if(srm11(i_file)<srm011(i_file)
 !    &.and.srm21(i_file)<srm021(i_file)) then

!	write(8,'(3i,12f)') i_file,ifa1,ifg1,rsqr31(i_file),rsqr21(i_file)
!     &,rsqr11(i_file) ,srm11(i_file),srm21(i_file),srm22(i_file),
!     &srm011(i_file),srm021(i_file),srm022(i_file), par(1),par(2),par(3)







!	if(srm22(i_file)<temp_srm) then

!	ite=ite+1

!	temp_srm = srm22(i_file)

!		write(*,*) temp_srm
	
		
!	i_ifa1 = ifa1
!	i_ifg1 = ifg1

!	t_rsqr11=rsqr31(i_file)
!	t_rsqr21=rsqr21(i_file)
!	t_rsqr31=rsqr11(i_file)
 !     t_srm11=srm11(i_file)
!	t_srm21=srm21(i_file)
!	t_srm22=srm22(i_file)
!	t_srm011=srm011(i_file)
 !     t_srm021=srm021(i_file)
!	t_srm022=srm022(i_file)
!	t_par1=par(1)
!	t_par2=par(2)
!	t_par3=par(3)


!	endif

!	write(*,*) "temp", itest1,ite,temp_srm

!	if(srm22(i_file)==temp_srm) then
!
!	write(33,'(3i,10i)') i_file,ifa1,ifg1,rsqr11(i_file),rsqr21(i_file)
!     & , srm11(i_file),srm21(i_file),srm22(i_file),srm011(i_file)
!     &,srm021(i_file),par(1),par(2),par(3)	

!	endif

!	endif


!	enddo
!	enddo

!	enddo


!	if(ite>0) then
!	write(33,'(3i,12f)')i_file,i_ifa1,i_ifg1,t_rsqr11,t_rsqr21,
 !    & t_rsqr31, t_srm11,t_srm21,t_srm22,t_srm011
 !    &,t_srm021,t_srm022,t_par1,t_par2,t_par3

!		endif






!	do i=1,data_sum(i_file)
!	write(8,*)i,qobs(i),qhat(i)
!	enddo
!	Close (unit=8)




 
!	endif


!	enddo

	enddo
	
	Close (unit=66)
	Close (unit=77)
      Close (unit=88)
      Close (unit=99)

	   
      stop
       end
c
c
        
c
c       to calculate objective function of Xinanjiang(2) model
c
         function xinob2(par)

        parameter(nmax=19000,npar=7,ncon=4,maxw=2500)
c       maxw > nacc*(m+2)  (nacc=24)
        dimension r(nmax),e(nmax),qobs(nmax),eobs(nmax),qt(nmax),
     &	        qi(n2+1),qg(n2+1),Tave(nmax)
        dimension par(npar),qhat(nmax),ehat(nmax),posi(17),
     &	        work1(maxw),work2(maxw),Xc(6),a1(17,nmax)
        Dimension h(200),consi(ncon),par0(npar),consit(ncon)
        common /rain/r,/evap/e,/disc/qobs,/actu/eobs,/qhat/qhat,
     &/ehat/ehat,/consi/consi,/uh/h,/Tave/Tave,/iyea_in/iyea_in,
     &/a1/a1,/Xc/Xc
        common /number/n1,n2,nwarm,m, /switch/isw,/area/area
        common /bndlow/parmin(npar),/bndup/parmax(npar)
	common /ifa1/ifa1,/posi/posi
c
c  if isw.eq.1, then initial conditions will be updated
c
         do 1 i=1,ncon
1        consit(i)=consi(i)
c
         do 2 i=1,npar
           par0(i)=par(i)
c           if(par0(i).lt.parmin(i)) par0(i)=parmin(i)
c           if(par0(i).gt.parmax(i)) par0(i)=parmax(i)
2       continue
c
c        cn=par0(14)
c        ck=par0(15)/cn
        nwarm=consi(3)
c        m=consi(8)

c
        do 5 i=n1,n2
            ehat(i)=0.0
5	   continue 
	 
c       call runoff generation subroutine
c
c	  call ffw117(r,e,ehat,qt,n1,n2,par0,consi,qi,qg)
c	  subroutine ffw117(p,ei,et,qt,n1,n2,par,consi,qi,qg)

!	write(*,*) a1(ifa1,1971:1990)

        call GR4J_MOD(n1,n2,par,consi,r,e,ehat,qobs,qhat,
     &Tave,iyea_in,Xc,ifa1,a1,posi)

		!write(*,*) ifa1,a1(ifa1,iyea_in(i))

c
c       SUBROUTINE GR4J_MOD(n1,n2,X0,consi,P,E,ES,Qobs,Q) 
c
c       Unit hydrograph of Gamma function
c        
c        d=1.0
c        nopt=2
c        nacc=24
c        maxs=nacc*(m+2)
c        nprint=0
c        ier=0
c
c        call   ffw42(h,d,m,cn,ck,nopt,nacc,work1,work2,maxs,nprint,ier)
c       Subroutine ffw42(h,d,m,an,ak,nopt,nacc,z,s,maxs,nprint,ier)
c
c       Discharge estimation
c																																									
        n=n1+nwarm
c
c       Estimated discharges qhat(i)
c
c         do 20 i=n1+m-1,n2
c            qhat(i)=0.0
c           do 30 j=1,m
c              qhat(i)=qhat(i)+h(j)*qt(i-j+1)                 
c30         continue
c		    qhat(i)=qhat(i)*area/86.4
c20      continue
c        
c	do 2111 i= n1+m-1,n2
c	qhat(i)=qhat(i)+qi(i)+qg(i)
c2111    continue
c       Sum of square errors

        f1=0.
        qbar=0.0
        qesbar=0.0
	  f3=0.
        ebar=0.0
        eesbar=0.0
        do 40 i=n,n2	      !n=n1+nwarm
	
	 
         qbar=qbar+qobs(i)
         qesbar=qesbar+qhat(i)
c	   f1=f1+(qobs(i)-qhat(i))**2
	   f1=f1+((qobs(i)-qhat(i))*86.4/area)**2

c	   把流量单位换成mm计算目标函数

  	     f3=f3+(eobs(i)-ehat(i))**2
	     ebar=ebar+eobs(i)
           eesbar=eesbar+ehat(i)
40      continue
        
	  
        f1=f1/(n2-n+1)
	  f3=f3/(n2-n+1)								   
        f2=0.0
        do 50 i=1,npar
            if(par(i).lt.parmin(i)) f2 = f2+abs(par(i)-parmin(i))
            if(par(i).gt.parmax(i)) f2 = f2+abs(par(i)-parmax(i))


50      continue
!	write(*,*) parmin(6),parmax(6)
c        CI=PAR(8),  CG=PAR(9),  CKI=PAR(12),  CKG=PAR(13)

c         if((par(8)+par(9)).gt.1.0) f2=f2+abs(par(8)+par(9)-1)
c	   if(par(12).gt.par(13)) f2=f2+abs(par(12)-par(13))
	  

c         主体为F1 F2为惩罚项

 !		xinob2=f1*(1+10000*f2 )
       xinob2=f1*(1+10000*f2 )
c     &	    +f3*(1+abs(ebar-eesbar)/ebar)
         
	  
c
         if(isw.eq.0) then
            do 100 i=1,ncon
100         consi(i)=consit(i)
         endif

c
         return
         end

c	   This subroutine reads the data from the standard data files used
c       in the flood forecasting workshop, 1985.
c
c

	Subroutine ffw1(x1,x2,x3,x4,nx,area,idimx,fyle,ititle,iyear,imon,
     1  iday,ihour,imin,isec,idt,l29,icode,iunit,nprint,ier)

c
	DIMENSION X1(IDIMX),X2(IDIMX),X3(IDIMX),X4(IDIMX),ITITLE(16),
     &	          MONS(12)
	Character*20 fyle,uname(5)
	COMMON /SYSPAR/IUNITI,IUNITO,IDEBUG
	DATA MONS/31,28,31,30,31,30,31,31,30,31,30,31/
	DATA UNAME / 'Not given','mm.','Cumecs.','inches','Cusecs.'/
c
c       Check the debug information
c
 	IF (IDEBUG.NE.0) THEN
	write(iunito,801)
801     Format(/,15x,'Debug information')
	 write(iunito,802)
802     Format(15x,30('-'))
	
	write(iunito,802)
	Endif
c
	OPEN(UNIT=30,FILE=FYLE,STATUS='unknown')
c
c       Read the first line. It should explain the contents.
c
	READ(30,100)(ITITLE(I),I=1,16)
100     FORMAT(16A4)
	if(nprint.eq.1)WRITE(IUNITO,200)FYLE,(ITITLE(I),I=1,16)
200     FORMAT(/,5x,' File:',a20,/,5x,' contains : ',16A4)
c
c       Read the second line. This should be the file-type 1.
c
	READ(30,300)I,IFREE
300     FORMAT(2(1X,I1))
c	IF (I.EQ.1) THEN
c       Time series file.
c       Read the third line. It should have the misc. information.
c
	IER=0
	READ(30,*)NX,IYEAR,IMON,IDAY,IHOUR,IMIN,
     1              ISEC,IDT,L29,ICODE,IUNIT,area
400     FORMAT(1X,I6,1X,I4,1X,5(I2,1X),I9,1X,I4,2(1X,I1),1x,F8.3)
	IF (nprint.eq.1) THEN
	WRITE(IUNITO,900)NX,IDAY,IMON,IYEAR,IHOUR,
     1                   IMIN,ISEC,IDT,L29
900     FormAT(/,5x,' This file contains ',I5,' values'/,5x,
     1              ' Starting from ',I2,'/',I2,'/',I4/,5x,
     2              '  at time ',I2,':',I2,':',I2/,5x,
     3              ' The time-interval is ',I9,' seconds.'/,5x,
     4              ' L29 is ',I4)
			IF (IUNIT.GE.0.AND.IUNIT.LE.4) THEN
	WRITE(IUNITO,901)UNAME(IUNIT+1)
901     FORMAT(/,5x,' Unit of data is : ',a20/,5x)
	ENDIF
	IF (ICODE.EQ.1) THEN
	WRITE(IUNITO,910)
910     FORMAT(6x,'These are averaged values.')
	ELSE
	WRITE(IUNITO,920)
920     FORMAT(' These are sampled values.')
	ENDIF
	ENDIF

	IF (NX.GT.IDIMX) THEN
c       Too many input values.
	if(idebug.ne.0)WRITE(IUNITO,500)NX,IDIMX,IDIMX
500     FORMAT(/,20x,' Warning - ffw1'/,20x,
     1               ' Data file contains ',I5,' values.'/,20x,
     2               ' This package is dimensioned for a max. of ',
     3               I5,/,20x,' Proceeding with ',I5,
     4               ' time-series values only')
	NX=IDIMX
	ier = 1
c
c       Non fatal error.
c
	ENDIF
c
c       Check that information is legal.
c
	IF(IUNIT.LT.0.OR.IUNIT.GT.4)THEN
	if(idebug.ne.0)write(IUNITO,902)IUNIT
902     FORMAT(/,20x,' Warning - ffw1, Illegal units '/,20x,
     1               ' Unit code ',I6,' is not defined.')
	IER=7
	ENDIF
	IF (IYEAR.LT.1800.OR.IYEAR.GT.2100) THEN
	if(idebug.ne.0)write(IUNITO,810)IYEAR
810     FORMAT(/,20x,' Warning - ffw1, Illegal year ',I4)
	IER=2
	ENDIF
	IF(IMON.LT.1.OR.IMON.GT.12) THEN
	if(idebug.ne.0)write(IUNITO,820)IMON
820     Format(/,20x,'Warning - ffw1, Illegal year',i4)
	IER=3
	ELSE
	NMON=MONS(IMON)
	IF(IMON.EQ.2) THEN
	IF(IYEAR.EQ.IYEAR/4*4) THEN
	NMON=29
	ENDIF
	ENDIF
	IF(IDAY.LT.1.OR.IDAY.GT.NMON) THEN
	if(idebug.ne.0)write(IUNITO,930)IDAY,IMON,IYEAR
930     FORMAT(/,20x,' Warning - ffw1 ',
     1               ' Illegal date ',
     2               I2,'/',I2,'/',I4)
	IER=4
	ENDIF
	ENDIF
	TIME=IHOUR+(IMIN+ISEC/60.)/60.
	IF(TIME.LT.0.0.OR.TIME.GT.24.) THEN
	if(idebug.ne.0)write(IUNITO,840)IHOUR,IMIN,ISEC,TIME
840     FORMAT(/,20x,' Warning - ffw1   Illegal time ',
     1            I2,':',I2,'.',I2,5X,'(',F10.5,')')
	IER=5
	ENDIF
c
	IF(IDT.LE.0) THEN
	if(idebug.ne.0)write(IUNITO,850)IDT
850     FORMAT(/,20x,' Warning - ffw1  Illegal time',
     1               ' increment ',I10)
	IER=6
	ENDIF
c
c       Now read the data.
c
c	IF(IFREE.LE.0) THEN

	DO 250 I=1,NX
	   READ(30,*)X1(I),X2(I),X3(I),X4(I)
250   CONTINUE

c600     FORMAT(2X,6E13.5)
c
c	ELSE
c	READ(30,*)(X(I),I=1,NX)
c	ENDIF
c       Sumarise the data as a check.
c
c	XMEAN=X(1)
c	XVAR=0.0
c	XMIN=XMEAN
c	XMAX=XMIN
c	DO 10 I=2,NX
c	TEMP=X(I)
c	IF (TEMP.GT.XMAX) THEN
c	XMAX=TEMP
c	ELSEIF(TEMP.LT.XMIN) THEN
c	XMIN=TEMP
c	ENDIF
c	XVAR=((I-2)*XVAR+(TEMP-XMEAN)**2*(I-1)/I)/(I-1)
c	XMEAN=((I-1)*XMEAN+TEMP)/I
c10      CONTINUE
c	if(nprint.eq.1)WRITE(IUNITO,1000)XMEAN,XVAR,XMIN,XMAX,NX
c1000    FORMAT(/,5x,' Short Summary of Data'/,5x,
c     1              ' Mean of Series ',t25,E13.6/,5x,
c     2              ' Variance of Series ',t25,E13.6/,5x,
c     3              ' Minimum value',t25,E13.6/,5x,
c     4              ' Maximum value',t25,E13.6/,5x,
c     5         /,5x,' No. of values used',t25,I6)
c	ELSE
c
c       Wrong fyle type.
c
c700        Format(/,30x,'Fatal Error - ffw1'/,20x,
c     1          ' This program expected a Time-series format data file'
c     2          ,/,20x,' (File-Type No. 1)'/,20x,
c     3          ' but instead found a file type no. ',I2)
c	IER=999
c	ENDIF
	Close (unit=30)
	RETURN
	END
c
	 subroutine ffw117(p,ei,et,qt,n1,n2,par,consi,qi,qg)
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
	PARAMETER(NPAR=15)
        DIMENSION P(n2), EI(n2+1),QT(n2+1),ET(n2+1),par(npar),
     &	        consi(8),par0(npar)
	  dimension qi(n2+1),qg(n2+1)
        COMMON /BNDLOW/PARMIN(NPAR),/BNDUP/PARMAX(NPAR)
	  COMMON /area/area
    
        DO 2 I=1,13
           PAR0(I)=PAR(I)
           IF(PAR0(I).LT.PARMIN(I)) PAR0(I)=PARMIN(I)
           IF(PAR0(I).GT.PARMAX(I)) PAR0(I)=PARMAX(I)
2       CONTINUE
c
        WM=PAR0(1)
        X=PAR0(2)
        Y=PAR0(3)
        CKE=PAR0(4)
        B=PAR0(5)
        SM=PAR0(6)
	EX=PAR0(7)
	CI=PAR0(8)
	CG=PAR0(9)
        CIMP=PAR0(10)
        C=PAR0(11)
	CKI=PAR0(12)
	CKG=PAR0(13)
c
	IF((CI+CG).GT.1.0) CG=1.0 - CI
c
        W1=CONSI(1)
        WU1=CONSI(2)
        WL1=CONSI(3)
	S1=CONSI(4)
	QI(n1)=CONSI(5)
	QG(n2)=CONSI(6)
c
        WUM=X*WM
        WLM=Y*(WM-WUM)
        WDM=WM-WUM-WLM        
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
        DO 100 I=N1,N2
            QT(I)=0.0
	      ET(I)=0.0
            EP=CKE*EI(I)
            PE=P(I)-EP
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
           IF(WU10.GT.EP) THEN
              EU=EP
           ELSE   
              EU=WU10
              EL=(EP-EU)*WL1/WLM
              IF(EL.LT.(C*(EP-EU))) EL=C*(EP-EU)
              IF(WL1.LT.(C*(EP-EU))) THEN  
                  EL=WL1
                  ED=C*(EP-EU)-EL
              IF(ED.GT.WD1) ED=WD1
              ENDIF
           ENDIF 
c
c     E:  actual evaporation for the time interval
c
           ET(I)=EU+EL+ED
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
          ET(I)=ET(I)*(1-CIMP)
c
c	RI and RG are routed through corresponding linear reservoirs
c	to obtain QI and QG.  These two flows then added to the RS to 
c 	get total inflow to the channel system, QT
c
	QI(i+1)=QI(i)*CKI+RI*(1-CKI)*area/86.4
	QG(i+1)=QG(i)*CKG+RG*(1-CKG)*area/86.4
	QT(I)=RS
c 
100     CONTINUE        
c
c       updating the initial conditions at the end of computing 
c       time
c
        consi(1)=W1
        consi(2)=WU1
        consi(3)=WL1
	consi(4)=S1
	consi(5)=qi(n2)
	consi(6)=qg(n2)
c
        return
        END
c        
c	     The multiparameter optimisation by the Simplex method
c
      SUBROUTINE FFW118(X,XMIN,XMAX,Y0,FUNC,NP,A,ATOL,TOLF,
     1                ITER,ITMAX,IPRT) 
      PARAMETER (NMAX=20,ALPHA=1.0,BETA=0.5,GAMMA=2.0,MMAX=NMAX+1)
      DIMENSION X(NP),XMIN(NP),XMAX(NP),P(MMAX,NMAX),Y(NMAX),
     1          PR(NMAX),PRR(NMAX),PBAR(NMAX),S(NMAX)
      COMMON/SYSPAR/IUNITI,IUNITO,IDEBUG
c
c  Calling arguments:
c
c  Simplex search method
c
c  X    -- Input as intial values and output as optimum
c          values of parameters
c  XMIN -- The lower bound of the parameters
c  XMAX -- The upper bound of the parameters
c  Y0   -- Output as value of objective function value at 
c          the optimum point
c  FUNC -- A user supplied external function, FUNC(X), which is to 
c          be optimised
c  NP   -- Number of parameters to be optimized
c  A    -- Input as intial size of standardized simplex
c  ATOL -- Stop condition for the distance between best and worst 
c          points
c  TOLF -- Stop condition for difference in function values at
c          the best and worst points
c  ITER -- total number of function evaluations
c  ITMAX-- Maximum allowed number of function evaluations
c  IPRT -- Print control.  If IPRT=0, no intermediate results printed.
c
c  Forming initial Simplex
c
           iter=0
           d1=a*(sqrt(np+1.0)+np-1.0)/(np*sqrt(2.0))
           d2=a*(sqrt(np+1.0)-1.0)/(np*sqrt(2.0))
C
           do 2 i=1,np
              s(i)=xmax(i)-xmin(i)
2          continue
           do 100 i=1,np+1
              do 200 j=1,np
                 if(S(j).lt.0) S(j)=1
                 if(abs(s(j)).le.atol) S(j)=0
                 if(i.eq.1) then
                    P(i,j)=X(j)
                 else
                    if(j.eq.i-1) then
                       P(i,j)=X(j)+d1*S(j)
                    else
                       P(i,j)=X(j)+d2*S(j)
                    end if
                 end if
200           continue
100        continue
c
           do 300 i=1,np+1
              do 310 j=1,np
                 X(j)=P(i,j)
310           continue
              Y(i)=func(X)
              iter=iter+1
300        continue

c
c  
      MPTS=NP+1
      ITT=0
c
c  First we determine which point is the highest (worst), next highest, 
c  and lowest (best).
c
1     ILO=1
      IF(Y(1).GT.Y(2))THEN
        IHI=1
        INHI=2
      ELSE
        IHI=2
        INHI=1
      ENDIF
      DO 11 I=1,MPTS
        IF(Y(I).LT.Y(ILO)) ILO=I
        IF(Y(I).GT.Y(IHI))THEN
          INHI=IHI
          IHI=I
        ELSE IF(Y(I).GT.Y(INHI))THEN
          IF(I.NE.IHI) INHI=I
        ENDIF
11    CONTINUE
c
c  Write intermediate results if required.
c  Check if the stop tolenrance ATOL is satisfied.
c
      IF(IPRT.NE.0) THEN
         WRITE(IUNITO,60)ITT,ITER
60       FORMAT(//4X,'Iteration number: ',i4,/
     1          4x,'Totle number of function evaluations: ',i4)
        WRITE(IUNITO,65)
65      FORMAT(/4X,'Objective function and corresponding vertex'/
     1        12x,'F',12x,'X')
        DO 70 I=1,NP+1
           WRITE(IUNITO,75)Y(I),(P(I,J),J=1,NP)
70      CONTINUE
75      FORMAT(4X,E13.6,2X,5(E13.6,2X))
      ENDIF
c
      DD=0.0
      DO 50 I=1,NP
        DD = DD+(P(IHI,I)-P(ILO,I))*(P(IHI,I)-P(ILO,I))
50    CONTINUE
      DD=SQRT(DD)
      IF((DD.LE.ATOL).OR.(ABS(Y(IHI)-Y(ILO)).LE.TOLF)
     1          .OR.ITER.GT.ITMAX) THEN
         Y0=Y(ILO)
         DO 55 I=1,NP
55       X(I)=P(ILO,I)
         WRITE(IUNITO,60)ITT,ITER
         WRITE(IUNITO,65)
         WRITE(IUNITO,75)Y0,(X(I),I=1,NP)
         RETURN
      ENDIF
c
      ITT=ITT+1
      DO 12 J=1,NP
        PBAR(J)=0.
12    CONTINUE
c
c  Begin a new iteration.  Compute the vector average of all points 
c  except the highest, i.e., the center of the "face" of the Simplex 
c  across from the high point.  We will subsequently explore along
c  the ray from the high point through that center
c
      DO 14 I=1,MPTS
        IF(I.NE.IHI)THEN
          DO 13 J=1,NP
            PBAR(J)=PBAR(J)+P(I,J)
13        CONTINUE
        ENDIF
14    CONTINUE
c
c  Extrapolate by a factor ALPHA through the face, i.e. reflect the
c  Simplex from the high point.
c
      DO 15 J=1,NP
        PBAR(J)=PBAR(J)/NP
        PR(J)=(1.+ALPHA)*PBAR(J)-ALPHA*P(IHI,J)
15    CONTINUE
c
c  Evaluate the function at the reflected point.  If gives a result
c  better than the best point, then try an additional extrapolation 
c  by a factor GAMMA.
c
      YPR=FUNC(PR)
      iter=iter+1
      IF(YPR.LT.Y(ILO))THEN
        DO 16 J=1,NP
          PRR(J)=GAMMA*PR(J)+(1.-GAMMA)*PBAR(J)
16      CONTINUE
        YPRR=FUNC(PRR)
        iter=iter+1
c
c  Check out the function there.  If the additional extrapolation
c  succeeded, then replaces the high point.
c 
       IF(YPRR.LT.Y(ILO))THEN
          DO 17 J=1,NP
            P(IHI,J)=PRR(J)
17        CONTINUE
          Y(IHI)=YPRR
        ELSE
c
c  The additional extrapolation failed, but we can still use 
c  reflected point.
c
          DO 18 J=1,NP
            P(IHI,J)=PR(J)
18        CONTINUE
          Y(IHI)=YPR
        ENDIF
      ELSE IF(YPR.GE.Y(INHI))THEN
c
c  The reflected point is worse than the second-highest.
c
        IF(YPR.LT.Y(IHI))THEN
c
c  If it is better than the highest, then replace the highest.
c
          DO 19 J=1,NP
            P(IHI,J)=PR(J)
19        CONTINUE
          Y(IHI)=YPR
        ENDIF
c
c  Look for an intermediate lower point, in other words, perform
c  a contraction of the Simplex along one dimension. Then evaluate
c  the function.
c
        DO 21 J=1,NP
          PRR(J)=BETA*P(IHI,J)+(1.-BETA)*PBAR(J)
21      CONTINUE
        YPRR=FUNC(PRR)
        iter=iter+1
        IF(YPRR.LT.Y(IHI))THEN
c
c  Contraction gives an improvement, so accept it.
c
          DO 22 J=1,NP
            P(IHI,J)=PRR(J)
22        CONTINUE
          Y(IHI)=YPRR
        ELSE
c
c  Can't seem to get rid of that point.  Better contract around
c  the lowest (best) point.
c
          DO 24 I=1,MPTS
            IF(I.NE.ILO)THEN
              DO 23 J=1,NP
                PR(J)=0.5*(P(I,J)+P(ILO,J))
                P(I,J)=PR(J)
23            CONTINUE
              Y(I)=FUNC(PR)
              iter=iter+1
            ENDIF
24        CONTINUE
        ENDIF
      ELSE
c
c  We arrive here if the original reflection gives a middling point.
c  Replace the old high point and continue.
c
        DO 25 J=1,NP
          P(IHI,J)=PR(J)
25      CONTINUE
        Y(IHI)=YPR
       ENDIF
      GO TO 1

c  For the test of doneness and the next iteration
       
	   
	   return
	   end

c
c
        SUBROUTINE FFW119(NP,X,XMIN,XMAX,FUNC,TOLX,TOLF,ITMAX,
     1                  IPRT,IT)
        PARAMETER(MP=20)
        DIMENSION X(NP),XMIN(NP),XMAX(NP),BETA(MP),V(MP,MP),
     1            S(MP,MP),W(MP)
        COMMON /XSCALE/XTEMP(MP)
        COMMON/SYSPAR/IUNITI,IUNITO,IDEBUG
c
c
c       Rosenbrock's method for multiple parameters optimisation
c       Reference "Computation methods in Hydrology" by G.C. Liang
c       (the subroutine follows lecture notes exactly, i.e, use
c       uniaxial search LINMIN at each stage)
c
c       Date:  9 March, 1992
c

C  Calling arguments
C  _________________________________________________________________
C  NP:    number of parameters to be optimised
C  X:     contains initial values of parameters on entry,
C         on exit, it contains optimised values
C  FUNC:  an user supplied external objective function, FUNC(X).
C         ***  has to be declared as an external function in the
C              calling program.
C  TOLX:  contains search stop criterion for search steps
C  TOLF:  comtains search stop criterion for the difference between
C         objective function values at two neighbour stages.
C  ITMAX: specifying maximum number of function evaluations
C  IPRT:  print control.  If IPRT not equal to zero, then the 
C         intermediate results are printed.
C  ----------------------------------------------------------------
C
c       IRTMAX:  Maximum number of rotations
c
        EXTERNAL FUNC
	  IRTMAX=100
c
        IRT=0
        IT=0
        DO 10 I=1,NP
           DO 20 J=1,NP
              V(I,J)=0.0
20         CONTINUE
           V(I,I)=1.0
10      CONTINUE
C
        FN=FUNC(X)
c
c       scalling the parameters
c
        do 401 j=1,np
          xtemp(j)=xmax(j)
          if(abs(xtemp(j)).lt.abs(xmin(j))) xtemp(j)=xmin(j)
          if(xtemp(j).eq.0.0) x(j)=0.0
          if(xtemp(j).ne.0.0) x(j)=x(j)/xtemp(j)
401     continue
c
        IT=IT+1
C
        SCAL=0.1
1000    CONTINUE
        FOLD=FN 
        BETAMX=0.0
c
c       Here SCALEV used to contain initial step for the next search
c
        SCALEV=0.0
        NP0=0
        DO 50 I=1,NP
           BETA(I)=0.0
           IF(XMIN(I).EQ.XMAX(I)) GOTO 50
           NP0=NP0+1
           DO 60 J=1,NP
              W(J)=V(I,J)
60         CONTINUE
c
c       calling subroutine LINMIN for one dimensional optimisation
c
        CALL LINMIN(X,W,NP,FRET,PLAMDA,SCAL,IT,FUNC)
        BETA(I)=PLAMDA
        SCALEV=SCALEV+ABS(PLAMDA)
        IF(BETAMX.LT.ABS(PLAMDA)) BETAMX=ABS(PLAMDA)
        FN=FRET
        IF(IT.GT.ITMAX) GOTO 301
50      CONTINUE
        SCAL=SCALEV/NP0
C
C       Print values of stage, function, variables
C
301     IF(IPRT.NE.0) THEN
           WRITE(IUNITO,300) IRT
300        FORMAT(//4X,'Stage number =',i4)
           WRITE(IUNITO,310)IT,FN
310        FORMAT(4X,'Number of function evaluations upto to now =',
     1     I4,/4X,'Value of the objective function = ',e16.6)
           WRITE(IUNITO,320)(X(I),I=1,NP)
320        FORMAT(4X,'Values of the parameters :',/4x,5(2x,e14.6))
        ENDIF
C
        IF(((BETAMX.LE.TOLX).OR.(ABS(FN-FOLD).LE.TOLF)).
     1      OR.(IT.GE.ITMAX)) GOTO 1100

C
C       STARTING ROTATING AXES
C
        DO 75 I=1,NP
        DO 75 J=1,NP
75      S(I,J)=0.0
C
C       Forming S vectors
C
        DO 100 I=1,NP
        DO 100 J=1,NP
           DO 110 K=I,NP
              IF(ABS(BETA(K)).LT.TOLX) GOTO 110
              S(I,J)=BETA(K)*V(K,J)+S(I,J)
110        CONTINUE
100     CONTINUE
C
C       Forming the first search direction
C
        IF(ABS(BETA(1)).GE.TOLX) THEN
           WNORM=0.0
           DO 150 J=1,NP
              W(J)=S(1,J)
              WNORM=WNORM+S(1,J)*S(1,J)
150        CONTINUE
           WNORM=SQRT(WNORM)
           DO 160 J=1,NP
              V(1,J)=W(J)/WNORM
160        CONTINUE
        ENDIF
C
C       Forming the 2 to n search direactions
C
        DO 205 I=2,NP
         IF(ABS(BETA(I)).GE.TOLX) THEN
           II=I-1
           DO 210 J=1,NP
              SUMV=0.0
              DO 220 K=1,II
                 SUMD=0.0
                 DO 230 KK= 1,NP
                    SUMD=SUMD+S(I,KK)*V(K,KK)
230              CONTINUE
                 SUMV=SUMD*V(K,J)+SUMV
220           CONTINUE
              W(J)=S(I,J)-SUMV
210        CONTINUE
           WNORM=0.0
           DO 250 J=1,NP
              WNORM=WNORM+W(J)*W(J)
250        CONTINUE
           WNORM=SQRT(WNORM)
           DO 260 J=1,NP
              V(I,J)=W(J)/WNORM
260        CONTINUE
         ENDIF
205     CONTINUE
C
        IRT=IRT+1
        IF(IRT.LT.IRTMAX) GOTO 1000
C
1100    WRITE(IUNITO,500)IRT,IT
500     FORMAT(//4X,'Total number of rotations :',i3/
     1           4x,'Total number of function evaluations: ',i4)
!        WRITE(IUNITO,550)(BETA(I),I=1,NP)
!550     FORMAT(/4X,'Final search step for each parameter: '/
!     1   4x,5(e13.6,2x))
c
c       scalling back to oroginal
c
        do 402 j=1,np
           x(j)=x(j)*xtemp(j)
402     continue
c
!        WRITE(IUNITO,560)FN,(X(I),I=1,NP)
!560     FORMAT(/4X,'The final objective function: ',e13.6/
!     1  4x,'Corresponding independent variables:'/
!     1  4x,5(e13.6,2x))
C
        RETURN
        END
C
C
        
      SUBROUTINE LINMIN(X,P,N,FRET,PLAMDA,S,IT,FUNOBJ)
c
c  To optimise PLAMDA in FUNC(X+PLAMDA*P), where
c  given an N dimensional point X and and N dimensional direction P,
c  moves and resets X to where the function FUNC(X) takes on a minimum 
c  along the direction P from X, and replaces P by the actual vector
c  displacement that X was moved. Also returns as FRET the value of
c  FUNC at the returned location X.  This is actually all acoomplished
c  by calling the routines MNBRAK and BRENT. 
c
c  S:  the initial guess of the search parameter, default value is one.
c  IT: for counting number of function evaluations.
c  FUNOBJ:  An external function which will be called in F1DIM for
c           calculating function values which is to be optimised.
c 
      PARAMETER (NMAX=50,TOL=1.E-2)
      EXTERNAL F1DIM
c   
c     F1DIM must accompany LINMIN
c
      DIMENSION X(N),P(N)
      COMMON /F1COM/ NCOM,XCOM(NMAX),PCOM(NMAX)
c
c     Set up the common block
c
      NCOM=N
      DO 11 J=1,N
        XCOM(J)=X(J)
        PCOM(J)=P(J)
11    CONTINUE
c
c   Initial guess for brackets.
c
      IF(S.EQ.0.0) S=1.0
      AX=0.
      XX=1.0*S
      BX=2.0*S
      CALL MNBRAK(AX,XX,BX,FA,FX,FB,F1DIM,IT,FUNOBJ)
      FRET=BRENT(AX,XX,BX,F1DIM,TOL,PLAMDA,IT,FUNOBJ)
c
c   Construct the vector results to return
c
      DO 12 J=1,N
        P(J)=PLAMDA*P(J)
        X(J)=X(J)+P(J)
12    CONTINUE
      RETURN
      END
c
c
        FUNCTION F1DIM(X,FUNOBJ)

c    Must accompany LINMIN
c
        PARAMETER(NMAX=50)
c
c    Maximum number of parameters
c
        COMMON/F1COM/NOCOM,XCOM(NMAX),PCOM(NMAX)
        DIMENSION XT(NMAX)
        common/xscale/xtemp(20)
c
        DO 11 J=1,NMAX
           XT(J)=XCOM(J)+X*PCOM(J)
11      CONTINUE
        do 20 j=1,7
           xt(j)=xt(j)*xtemp(j)
20      continue
c
        F1DIM=FUNOBJ(XT)
c
        RETURN
        END
c

      SUBROUTINE MNBRAK(AX,BX,CX,FA,FB,FC,FUNC,IT,FUNOBJ)
c
c   The first parameter is the default ratio by which successive 
c   intervals are magnified; the second is the maximum magnification 
c   allowed for a parabolic-fit step.
c 
      PARAMETER (GOLD=1.618034, GLIMIT=100., TINY=1.E-10)
c
c   Given a function FUNC, and ginven distinct initial points AX 
c   BX, this routine searches in the downhill direction (defined by 
c   the function as evaluated at the initial points) and returns new
c   points AX, BX, CX which bracket a minimum of the function.  Also
c   returned are the function values at the three points, FA, FB 
c   and FC.
c
      FA=FUNC(AX,FUNOBJ)
      FB=FUNC(BX,FUNOBJ)
      IT=IT+2
      IF(FB.GT.FA)THEN
        DUM=AX
        AX=BX
        BX=DUM
        DUM=FB
        FB=FA
        FA=DUM
      ENDIF
      CX=BX+GOLD*(BX-AX)
      FC=FUNC(CX,FUNOBJ)
      IT=IT+1
c
c     modified by Liang, 20 March, 1992
c     originally it is IF(FB.GE.FC) THEN
c     the scheme will not converge if FB=FC
c
1     IF(FB.GT.FC)THEN
        R=(BX-AX)*(FB-FC)
        Q=(BX-CX)*(FB-FA)
        U=BX-((BX-CX)*Q-(BX-AX)*R)/(2.*SIGN(MAX(ABS(Q-R),TINY),Q-R))
        ULIM=BX+GLIMIT*(CX-BX)
        IF((BX-U)*(U-CX).GT.0.)THEN
          FU=FUNC(U,FUNOBJ)
          IT=IT+1
          IF(FU.LT.FC)THEN
            AX=BX
            FA=FB
            BX=U
            FB=FU
            GO TO 1
          ELSE IF(FU.GT.FB)THEN
            CX=U
            FC=FU
            GO TO 1
          ENDIF
          U=CX+GOLD*(CX-BX)
          FU=FUNC(U,FUNOBJ)
          IT=IT+1
        ELSE IF((CX-U)*(U-ULIM).GT.0.)THEN
          FU=FUNC(U,FUNOBJ)
          IT=IT+1
          IF(FU.LT.FC)THEN
            BX=CX
            CX=U
            U=CX+GOLD*(CX-BX)
            FB=FC
            FC=FU
            FU=FUNC(U,FUNOBJ)
            IT=IT+1
          ENDIF
        ELSE IF((U-ULIM)*(ULIM-CX).GE.0.)THEN
          U=ULIM
          FU=FUNC(U,FUNOBJ)
          IT=IT+1
        ELSE
          U=CX+GOLD*(CX-BX)
          FU=FUNC(U,FUNOBJ)
          IT=IT+1
        ENDIF
        AX=BX
        BX=CX
        CX=U
        FA=FB
        FB=FC
        FC=FU
        GO TO 1
      ENDIF
      RETURN
      END
c
      FUNCTION BRENT(AX,BX,CX,F,TOL,XMIN,IT,FUNOBJ)
c
c  Given a function FUNC, and given a bracketing triplet of abscissas
c  AX, BX, CX (such that BX is between AX and CX, and F(BX) is less
c  than both F(AX) and F(CX), this routine isolates the minimum to
c  a fractional precision of about TOL using Brent's method.  The
c  abscissa of the minimum is returned as XMIN, and the minimum
c  function is returned as BRENT, the returned function value.
c
      PARAMETER (ITMAX=100,CGOLD=.3819660,ZEPS=1.0E-10)
      A=MIN(AX,CX)
      B=MAX(AX,CX)
      V=BX
      W=V
      X=V
      E=0.
      FX=F(X,FUNOBJ)
      IT=IT+1
      FV=FX
      FW=FX
      DO 11 ITER=1,ITMAX
        XM=0.5*(A+B)
        TOL1=TOL*ABS(X)+ZEPS
        TOL2=2.*TOL1
        IF(ABS(X-XM).LE.(TOL2-.5*(B-A))) GOTO 3
        IF(ABS(E).GT.TOL1) THEN
          R=(X-W)*(FX-FV)
          Q=(X-V)*(FX-FW)
          P=(X-V)*Q-(X-W)*R
          Q=2.*(Q-R)
          IF(Q.GT.0.) P=-P
          Q=ABS(Q)
          ETEMP=E
          E=D
          IF(ABS(P).GE.ABS(.5*Q*ETEMP).OR.P.LE.Q*(A-X).OR. 
     *        P.GE.Q*(B-X)) GOTO 1
          D=P/Q
          U=X+D
          IF(U-A.LT.TOL2 .OR. B-U.LT.TOL2) D=SIGN(TOL1,XM-X)
          GOTO 2
        ENDIF
1       IF(X.GE.XM) THEN
          E=A-X
        ELSE
          E=B-X
        ENDIF
        D=CGOLD*E
2       IF(ABS(D).GE.TOL1) THEN
          U=X+D
        ELSE
          U=X+SIGN(TOL1,D)
        ENDIF
        FU=F(U,FUNOBJ)
        IT=IT+1
        IF(FU.LE.FX) THEN
          IF(U.GE.X) THEN
            A=X
          ELSE
            B=X
          ENDIF
          V=W
          FV=FW
          W=X
          FW=FX
          X=U
          FX=FU
        ELSE
          IF(U.LT.X) THEN
            A=U
          ELSE
            B=U
          ENDIF
          IF(FU.LE.FW .OR. W.EQ.X) THEN
            V=W
            FV=FW
            W=U
            FW=FU
          ELSE IF(FU.LE.FV .OR. V.EQ.X .OR. V.EQ.W) THEN
            V=U
            FV=FU
          ENDIF
        ENDIF
11    CONTINUE
      PAUSE 'Brent exceed maximum iterations.'
3     XMIN=X
      BRENT=FX
      RETURN
      END

c
       Subroutine ffw120(x,fb,func,nvar,npop,nstr1,neva,pm,cmax,                  
     1                xmin,xmax,iseed,itga)
                                                                                     !!!  func=xinbob2
c      Author :  Q.J. Wang
c                Department of Engineering Hydrology
c                University College Galway
c                Ireland

c      Date   :  29 May, 1991
c
c      Modified by:   G.C.Liang, on 6th June, 1991
c
c      Genetic algorithm for function minimization
c      Ranking procedure is used to assign selection probabilities

c      Variable descriptions:
c      x     -- An array, on entry contains initial guessed parameters 
c               and return as the best point
c      fb    -- Return as the objective function value of the best point
c      func  -- Input as the external function, func(x), to be minimized            !!!!!!!!!!
c      nvar  -- Input as the number of variables
c      npop  -- Input as the number of population
c      nstr1 -- Input as the lenth of binary code for each variable
c      neva  -- Input as the number of evaluations of the objective function
c      pm    -- Input as the probability of mutation
c      cmax  -- Input as the expected value of the number of copies of 
c               the best point of each generation
c      xmin  -- An array, input as the lower bound of x region, will not
c               be changed on exit.
c      xmax  -- An array, input as the upper bound of x region
c               If xmax(i)=xmin(i), then the ith parameter will not
c               be optimized by the subroutine. This array will not
c               be changed on exit.
c      iseed -- Input as the seed for initializing the generator of random
c               numbers distributed uniformly between 0 and 1; Set iseed
c               to any negative integer value

       parameter(npopmx=150,nstrmx=80000,npmax=50)
c
c      nstrmx > nvar*nstr1
c
       integer pop(npopmx,nstrmx),pop1(npopmx,nstrmx)
       integer who(nstrmx),hero(nstrmx)
       integer note(nstrmx)
       dimension p(npopmx),pcum(npopmx),f(npopmx),m(npopmx)
       dimension x(nvar),xmin(nvar),xmax(nvar),xtemp(nstrmx)
       dimension xmin0(npmax),xmax0(npmax)

c      Total length of the binary codes of all the variables
c
       do 1 i=1,nvar
          xmin0(i)=xmin(i)
          xmax0(i)=xmax(i)
1      continue
c
       nvar0=0
       do 2 i=1,nvar
         note(i)=1
         if(xmin0(i).eq.xmax0(i)) note(i)=0
         nvar0=nvar0+note(i)
         xtemp(i)=xmin0(i)
2      continue
c
       nstr=nvar0*nstr1

c      Number of generations

       ng=neva/npop

c      Assign the selection probabilities to the ranked points.
c      Also find the cumulative probabilities.

       p(1)=(2.0-cmax)/npop
       p(npop)=cmax/npop
       deltap=(p(npop)-p(1))/(npop-1)
       do 10 i=2,npop-1
          p(i)=p(i-1)+deltap
10     continue
       pcum(1)=p(1)
       pcum(npop)=1
       do 20 i=2,npop-1
          pcum(i)=pcum(i-1)+p(i)
20     continue

c      Initialize the (0,1) random number generator
       idum=iseed
       idum0=-iseed
       ir=irbit(idum0)

c      Create the initial population

       do 30 i=1,npop
          do 40 k=1,nstr
             ir=irbit(idum0)
                pop(i,k)=ir
40        continue
30     continue
c
c      calculating the objective function for the initial guessed 
c      parameters
c
c      Rearrange the parameters based on whether it will be optimized
c      or not.
c
       i1=0
       do 3 i=1, nvar
       if(note(i).eq.0) goto 3
       i1=i1+note(i)
       x(i1)=x(i)
       xmin0(i1)=xmin0(i)
       xmax0(i1)=xmax0(i)
3      continue
c
c    1.  converting x values to binary codes
c
       ntemp=2**nstr1-1
       j1=0
       do 5 i=1,nvar0
          ii=(x(i)-xmin0(i))/(xmax0(i)-xmin0(i))*ntemp+0.5
          call decbin(ii,who,nstr1)
          do 7 j=1,nstr1
             j1=j1+1
             hero(j1)=who(j)
7         continue
5      continue
c
c  2. Converts binary codes to decimal form
c
             
             call pheno(x,hero,nvar0,nstr1,xmin0,xmax0)
c
             ntemp=nvar0+1
             do 11 j=1,nvar
             j1=nvar-j+1
               ntemp=ntemp-note(j1)
               if(note(j1).eq.1) then
                  x(j1)=x(ntemp)              
               else
                  x(j1)=xtemp(j1)
               endif
11           continue
c
c  3. calculates the objective function
c
             fb=func(x)
             itga=itga+1
c              
       do 45 ig=1,ng

c         Calculate the objective function value for each point

          do 50 i=1,npop
             do 60 k=1,nstr
                who(k)=pop(i,k)                                                                                                        
60           continue
             call pheno(x,who,nvar0,nstr1,xmin0,xmax0)
c
             ntemp=nvar0+1
             do 55 j=1,nvar
             j1=nvar-j+1
               ntemp=ntemp-note(j1)
               if(note(j1).eq.1) then
                  x(j1)=x(ntemp)              
               else
                  x(j1)=xtemp(j1)
               endif
55           continue
c
             f(i)=func(x)
             itga=itga+1
50        continue

c         Rank the points so that their objective function values
c         are in a descending order, i.e., f(i).ge.f(i+1). 

          do 70 i=1,npop
             m(i)=i
70        continue
          do 80 i=1,npop-1
             do 90 j=i+1,npop
                if(f(i).lt.f(j)) then
                   m0=m(i)
                   f0=f(i)
                   m(i)=m(j)
                   f(i)=f(j)
                   m(j)=m0
                   f(j)=f0
                endif
90           continue
80        continue

          do 100 i=1,npop
             do 110 k=1,nstr
                pop1(i,k)=pop(m(i),k)
110          continue
100       continue

c         Find the best point so far and insert the best point of the 
c         last population into the middle rank if it is better than
c         the present best point
c
          if(fb.gt.f(npop)) then
                fb=f(npop)
                do 116 k=1,nstr
                   hero(k)=pop1(npop,k)
116             continue
          else
            do 117 i=1, (npop+1)/2-1
              do 117 k=1, nstr
                pop1(i,k)=pop1(i+1,k)
117         continue
            mid=(npop+1)/2
            do 118 k=1,nstr
118         pop1(mid,k)=hero(k)
c
          endif

c         Reproduction

          do 120 i=1,npop

c            Find a pair of points

             r=ran1(idum)
             i1=0
130          continue
                i1=i1+1
             if(r.gt.pcum(i1)) goto 130

             r=ran1(idum)
             i2=0
140          continue
                i2=i2+1
             if(r.gt.pcum(i2)) goto 140

c            Crossover and mutation

c            Find two crossover points around the circle

             r=ran1(idum)
             k1=int(r*nstr)+1

145          continue
                r=ran1(idum)
                k2=int(r*nstr)+1
             if(k2.eq.k1) goto 145

c           Crossover

             if(k1.gt.k2) then
                k=k1
                k1=k2
                k2=k
             endif

             do 150 k=k1,k2-1
                pop(i,k)=pop1(i1,k)
150          continue
             do 160 k=k2,nstr
                pop(i,k)=pop1(i2,k)
160          continue
             do 170 k=1,k1-1
                pop(i,k)=pop1(i2,k)
170          continue

c            Mutation

             do 180 k=1,nstr
                r=ran1(idum)
                if(r.le.pm) then
                   pop(i,k)=abs(pop(i,k)-1)
                endif
180          continue

120       continue

45     continue

       call pheno(x,hero,nvar0,nstr1,xmin0,xmax0)
c
             ntemp=nvar0+1
             do 46 j=1,nvar
             j1=nvar-j+1
               ntemp=ntemp-note(j1)
               if(note(j1).eq.1) then
                  x(j1)=x(ntemp)              
               else
                  x(j1)=xtemp(j1)
               endif
46           continue
c
       return
       end


       subroutine pheno(x,who,nvar,nstr1,xmin,xmax)

       parameter(nstrmx=80000)
       integer who(nstrmx),me(nstrmx)
       real x(nvar),xmin(nvar),xmax(nvar)

       ww=2**nstr1-1

       do 20 j=1,nvar
          j0=(j-1)*nstr1
          do 30 k=1,nstr1
             me(k)=who(j0+k)
30        continue
          call bindec(i,me,nstr1)
          x(j)=xmin(j)+(xmax(j)-xmin(j))*i/ww
20     continue

       return
       end


       subroutine bindec(i,k,l)

c      Convert an integer from binary form to decimal form
c      Knowing k(l) to calculate i
c
       dimension k(l)

       i=k(l)
       i1=1
       do 100 j=2,l
          j1=l-j+1
          i1=i1*2
          i=i+k(j1)*i1
100    continue
        
       return
       end
c
c
       subroutine decbin(i,k,l)
c
c     convert an integer from decimal form to binary form
c     knowing i to calculate k(l)
c
      dimension k(l)
c
      i0=i
      do 10 j=1,l-1
        j1=l-j
        k(j)=0 
        jj=2**j1
        if(i0.ge.jj) then
           k(j) = 1
           i0=i0-jj
        endif
10    continue
      k(l)=i0
      return
      end
c
c

      FUNCTION RAN1(IDUM)
      DIMENSION R(97)
      PARAMETER (M1=259200,IA1=7141,IC1=54773,RM1=3.8580247E-6)
      PARAMETER (M2=134456,IA2=8121,IC2=28411,RM2=7.4373773E-6)
      PARAMETER (M3=243000,IA3=4561,IC3=51349)
      DATA IFF /0/
      IF (IDUM.LT.0.OR.IFF.EQ.0) THEN
        IFF=1
        IX1=MOD(IC1-IDUM,M1)
        IX1=MOD(IA1*IX1+IC1,M1)
        IX2=MOD(IX1,M2)
        IX1=MOD(IA1*IX1+IC1,M1)
        IX3=MOD(IX1,M3)
        DO 11 J=1,97
          IX1=MOD(IA1*IX1+IC1,M1)
          IX2=MOD(IA2*IX2+IC2,M2)
          R(J)=(FLOAT(IX1)+FLOAT(IX2)*RM2)*RM1
11      CONTINUE
        IDUM=1
      ENDIF
      IX1=MOD(IA1*IX1+IC1,M1)
      IX2=MOD(IA2*IX2+IC2,M2)
      IX3=MOD(IA3*IX3+IC3,M3)
      J=1+(97*IX3)/M3
      IF(J.GT.97.OR.J.LT.1)PAUSE
      RAN1=R(J)
      R(J)=(FLOAT(IX1)+FLOAT(IX2)*RM2)*RM1
      RETURN
      END
c
c
    
C
C
      FUNCTION IRBIT(ISEED)
c
c       ISEED must be any integer value on entry.       
c       Returns as an integer a random bit, based on the low
c       -significance bits in ISEED (which is modified for the 
c       next call)
      LOGICAL NEWBIT
      PARAMETER (IB1=1,IB4=8,IB6=32,IB30=536870910)
      NEWBIT=IAND(ISEED,IB30).NE.0
      IF(IAND(ISEED,IB6).NE.0) NEWBIT=.NOT.NEWBIT
      IF(IAND(ISEED,IB4).NE.0) NEWBIT=.NOT.NEWBIT
      IF(IAND(ISEED,IB1).NE.0) NEWBIT=.NOT.NEWBIT
      IRBIT=0
      ISEED=IAND(ISHFT(ISEED,1),NOT(IB1))
      IF(NEWBIT) THEN
         IRBIT=1
         ISEED=IOR(ISEED,IB1)
      ENDIF
      RETURN
      END
c
c


       SUBROUTINE GR4J_MOD(n1,n2,X0,consi,P,E,ES,Qobs,Q,Tave
     &,iyea_in,Xc,ifa1,a1,posi)   
	! use numerical_libraries

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

      PARAMETER (NPX=7,NPAR=7,NH=10,ncon=4)

      DIMENSION V(3*NH+2),X(7+3*NH),X0(NPAR),consi(ncon),Xc(6)
	DIMENSION P(n2),E(n2),EP(n2),ES(n2),Q(n2),Qobs(n2),posi(17)
     &	,melt(n2),Tave(n2),P1(n2),P2(n2),iyea_in(n2),a1(17,n2)
	dimension dpar(NPAR)
	INTEGER n1,n2,I,ifa1,ilength,nlength
	REAL PSQ,PSP,QB

	COMMON /area/area,/dpar/dpar
C    - V(1): Level in production store 
C    - V(2): level in routing store 
C    - V(3) to V(2+3*NH): Storages for UH1 and UH2 initially set to zero 
C***************************************
C    - X(5) to X(NPX+3*NH): Ordinates of UH1 & UH2
!

   
!	X0(1)=229.8586
!	X0(2)=-2.209420
!	X0(3)=67.14509
!	X0(4)=4.344584
!	X0(5)=0.1051328
!	X0(6)=3.059196





      DO 3 I=1,7
	  X0(I)=dpar(I)
	X(I)= X0(I)

!	write(*,*) dpar(I)
    3	CONTINUE

	


	
	

	
       CKE=X(5)
	DDF=X(6)

	snowpack=CONSI(4)

	CALL UH1(X,X(4))     !得到X(5)~X(14),UH1的横坐标值
      CALL UH2(X,X(4))     !得到X(15)~X(34),UH2的横坐标值

C       PSP=SUM(P(n1:n2))/(n2-n1+1)
C 	PSQ=SUM(Qobs(n1:n2))/(n2-n1+1)
C  	PSQ=PSQ*24*3.6/area
	V(1)=CONSI(1)        !要改为读par文件里常数
	V(2)=CONSI(2)


	DO I=3,2+3*NH 
         V(I)=0. 
      ENDDO
C       QB=0.9*PSQ/X(3) 
C 	 V(2)=1.-1./(4.*(1.+QB)**4-3.+(4.*QB)**0.2+(8.*QB)**0.4) 
C        V(2)=V(2)*X(3) 
C        V(1)=X(1)*MIN(1.,SQRT(PSQ/PSP)) 


C Production store: 
   
	ilength=1
	nlength=0   

 !  	write(*,*) i,X0(1),X0(3)

 	if(X0(3)+X0(7)<1.0) then
          
          X(3)=1.0

	X0(3)=1.0-X0(7)

	endif

	if(X0(3)+X0(7)>200.0) then
	
          X(3)=200.0 
	X0(3)=200.0-X0(7)

	endif


      DO 100  I=n1,n2
	


	IF(I>1)	THEN
	IF(iyea_in(i)>iyea_in(i-1)) THEN

	ilength=ilength+1


	endif
	endif

!	X(2)=X0(1)

!	write(*,*) i,posi(ifa1),ilength,iyea_in(i)
	
	if(ilength<=posi(ifa1)) then

	nlength=nlength+1

      else 
!xiugai2: parameter value
          !ceshi
		X(3)=X0(3)+X0(7)*a1(ifa1,iyea_in(i-nlength))






	endif

!	write(*,*) X(1:7)






!	write(*,*) X(2),X(3),ifa1,iyea_in(i),a1(ifa1,iyea_in(i))

		P1(I)=P(I)

!	write(*,*) p(i),E(i),Tave(i)
!		P(I)=P2(I)
!		P1(I)=P2(I)
	melt(I)=0.0
		AA=-0.0

	IF (Tave(I).gt.AA) then
		melt(I)=DDF*(Tave(I)-AA)		!T>0 melt
		cc=melt(I);
		IF (melt(I).gt.snowpack) melt(I)=snowpack
		P2(I)=P1(I)+melt(I)		!T>0 rain=rain+melt
		snowpack=snowpack-melt(I)
		cc=melt(I);
	else
		P2(I)=0		!T<0 snow
		snowpack=snowpack+P1(I)
	endif

	EP(I)=E(I)*CKE*37.586/15.54

      IF(P2(I).GE.EP(I))THEN                   
	ES(I)=0.    
	WS=(P2(I)-EP(I))/X(1)   
	IF(WS.GT.13)WS=13     
	PS=X(1)*(1.-(V(1)/X(1))**2.)*tanh(WS)/(1.+V(1)/X(1)*tanh(WS))                  
	PR=P2(I)-EP(I)-PS   !不需要加perc???                
	ELSE    
      WS=(EP(I)-P2(I))/X(1)          
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
 	Q(I)=Q(I)*1.0E-3*area*1.0E6/(24*3600)

	IF(I==365) THEN

!	WRITE(*,*) "OK"

	ENDIF



	
100    CONTINUE 
 
c    updating the initial conditions at the end of computing 
c    time
c
        consi(1)=V(1)
        consi(2)=V(2)

	consi(4)=snowpack
      
        return
        END

 
 
C**************************************************************************************** 
      SUBROUTINE UH1(X,C) 
      PARAMETER (NPX=6,NH=10) 
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
      PARAMETER (NPX=6,NH=10) 
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



c       This subroutine calculates the sum of squares of differences
c       between the the observed series X and the estimated series Xhat.
c       It also calcultes the model efficiency R square.
c       The difference from ffw17 is that this routine also provides
c       the mean square errors in different flow zones. At present,
c       the flow zones are limited to a number of three and only        
c       the upper limits of zone-1 and zone-2 are to be entered.
c
	subroutine ffw17(x,y,n1,n2,xhat,yhat,xbar,ybar,f1,fhat1,fagr1,
     &	 fagr11,f2,fhat2,fagr2,fagr22,f3,f4,rsqr1,rsqr2,nprint,ier,
     &     ratio,RR,RIOA1,RIOA2)
c
	Dimension x(n2),xhat(n2),y(n2),yhat(n2),y11(10),y22(10),fmse(10)
	common    /syspar/iuniti,iunito,idebug
c
c       Print debugging information if required.
c
	if(idebug.ne.0)then
	Write(iunito,1000)
1000    Format(/,15x,'Debug information',/,15x,30('-'))
	write(iunito,1001)x(n1),x(n2),xhat(n1),xhat(n2),n1,n2,xbar
1001    Format(/,15x,'Subroutine ffw17 called',/,15x,'x(n1) = ',t30,
     1  e13.6, /,15x,'x(n2) = ',t30,e13.6,
     2         /,15x,'xhat(n1) = ',t30,e13.6,
     3         /,15x,'xhat(n2) = ',t30,e13.6,
     4         /,15x,'n1 = ',t30,i4,/,15x,'n2 = ',t30,i4,/,15x,
     5                'xbar = ',t30,e13.6,/,15x,30('-'))
	endif

c
c       Check if the information is legal.
c
	if(n1.le.0)then
	if(idebug.ne.0)then
	write(iunito,998)
998     Format(/,20x,'Illegal parameter n1 in ffw17')
	endif
	n1 = 1
	ier = 1
	endif
c
	if(n1.gt.n2)then
	write(iunito,995)n2
995     Format(/,20x,'Illegal parameter in ffw17 - n2 = ',i3)
	ier = 999
	go to 999
	endif
c
c       calculate the mean of x series and xhat series.均值
c
	n = n2-n1+1
	sxhat = 0.0
	sx = 0.0
	syhat = 0.0
	sy = 0.0
	do 15 j = n1,n2
	sx = sx +x(j)
	sxhat = sxhat + xhat(j)
	sy = sy +y(j)
	syhat = syhat + yhat(j)
15      continue
	sx = sx/n
	sxhat = sxhat/n
	sy = sy/n
	syhat = syhat/n
c
c       calculate the initial variance and the residual variance after
c       fitting the model expressed /day.标准差和均方差
c
	f1 = 0.0
	fhat1 = 0.0
      fagr1 = 0.0
      fagr11 = 0.0
	f2 = 0.0
	fhat2 = 0.0
      fagr2 = 0.0
      fagr22 = 0.0
	f3 = 0.0
	f4 = 0.0
	do 20 j = n1,n2
	f1 = f1+(x(j)-xbar)**2
	fhat1 = fhat1+(x(j)-xhat(j))**2
      fagr1 = fagr1+abs(x(j)-xhat(j))
      fagr11 = fagr11+(abs(xhat(j)-xbar)+abs(x(j)-xbar))
	f2 = f2+(y(j)-ybar)**2
	fhat2 = fhat2+(y(j)-yhat(j))**2
      fagr2 = fagr2+abs(y(j)-yhat(j))
      fagr22 = fagr22+(abs(yhat(j)-ybar)+abs(y(j)-ybar))
	f3 = f3+(yhat(j)-syhat)**2
	f4 = f4+(y(j)-ybar)*(yhat(j)-syhat)
20      continue
	f1 = f1/n
	fhat1 = fhat1/n
      fagr1 = fagr1/n
      fagr11 = fagr11/n
	f2 = f2/n
	fhat2 = fhat2/n
      fagr2 = fagr2/n
      fagr22 = fagr22/n
	f3 = f3/n
	f4 = f4/n
c
c       Calculate the index rsqr.   rsqr确定性系数 ratio水量平衡系数,RR蒸发相关系数，一致性系数IOA
c
	rsqr1 = (1-(fhat1/f1))*100.0
	rsqr2 = (1-(fhat2/f2))*100.0
	RR = (f4**2)/(f2*f3)
      RIOA1 = (1-(fagr1/fagr11))*100.0
      RIOA2 = (1-(fagr2/fagr22))*100.0
c
c       check if printing of the results are required.
c
	if(nprint.eq.1)then
	Write(iunito,500)
500     Format(///,10x,'TABLE',/,10x,15('-'),//,10x,'Catchment',/,10x,
     1  'Model',/,10x,'Calibration/verification period',/,10x,
     2  'Design/updating mode',//)

	ratio = sxhat/sx

	write(iunito,505)n,n1,n2,xbar,sx,sxhat,ratio,f1,fhat1,rsqr1,rsqr2
505     Format(/,10x,'Short summary of the results',/,10x,29('-'),
     1  //,10x,'(for',i5,' values, from',i5,' to',i5,' )',/,
     2  10x,60('-'),/,10x,'1. Mean of the outflow in calibration'
     3  ,t50,'=',e13.6,//,10x,'2. Mean of the observed series',t50,
     4  '=',e13.6,/,10x,'3. Mean of the estimated series',t50,'=',
     5  e13.6,/,10x,'4. Ratio of the estimated to the ',
     6  /,10x,'   observed mean of the outflow',
     7  t50,'=',f10.4,//,10x,'5. The initial S.O.S per unit time',t50,
     8  '=',e13.6,/,10x,'6. The final S.O.S per unit time',t50,'=',
     9  e13.6,/,10x,'7. The performance index (R sqr. %)',t50,'=',
     1  f10.2,//,10x,'8. The performance index (E sqr. %)',t50,'=',
     2	  f10.2,/,10x,60('-'))
	endif
c
	If(nprint .eq. 2)then

	nband = 3
	write(iunito,510)
510     format(//2x,'Enter the band limits (2 values):',t60,'=>')
	read(iuniti,*)y22(1),y22(2)
	y22(3) = 1.0e10
	y11(1) = 0.0    
	do 600 i=2,nband
	y11(i) = y22(i-1)                       
600     continue

	do 720 j=1,nband
	fmse(j) = 0.0
	num = 0 
	do 710 i=n1,n2
	if(x(i) .ge. y11(j) .and. x(i).lt.y22(j))then   
	fmse(j) = fmse(j) + (x(i)-xhat(i))**2
	num = num + 1
	endif
710     continue
	if(num .gt. 0)then
	fmse(j) = fmse(j)/num
	else
	fmse(j) = -9.9
	endif
	write(6,640)j,num,j,fmse(j)
640     format(/2x,'Number of flows in zone-',i1,' = ',i5,
     $         /2x,'MSE of zone-',i1,' flows = ',e13.6)
720     continue

	Endif
999     return
	end

c
c       This subroutine writes the data in a standard time series
c       data files used in the flood forecasting workshop, 1985.
c
c

        Subroutine ffw2(x,nx,idimx,fyle,ititle,iyear,imon,iday,ihour,
     1                 imin,isec,idt,l29,icode,iunit,nprint,ier)
c
c
c
        DIMENSION X(IDIMX),ITITLE(16),MONS(12)
        Character*20,fyle,uname(5)
        COMMON /SYSPAR/IUNITI,IUNITO,IDEBUG
        DATA MONS/31,28,31,30,31,30,31,31,30,31,30,31/
        DATA UNAME / 'Not given','mm.','Cumecs.','inches','Cusecs.'/
c
c        Check the debug information
c
        IF (IDEBUG.NE.0) THEN
        write(iunito,801)
801        Format(/,15x,'Debug information')
        write(iunito,802)
802        Format(15x,30('-'))
                WRITE(IUNITO,800)FYLE,IDIMX
800                FORMAT(/,15x,'Subroutine ffw2 called to'
     1           /,15x,'write file :          ',a20,
     2           /,15x,'idimx =              ',i6)
        write(iunito,802)
        ENDIF
c
        open(unit=20,file=fyle,access='append',status='old',err=805)
	go to 806
805	open(unit=20,access='sequential',file=fyle,status = 'new')
806	continue
c
c        write the first line. It should explain the contents.
c
        write(20,100)(ITITLE(I),I=1,16)
100        FORMAT(16A4)
        if(nprint.eq.1)WRITE(IUNITO,200)FYLE,(ITITLE(I),I=1,16)
200        FORMAT(/,5x,' File:',a20,/,5x,' contains : ',16A4)
c
c        write the second line. This should be the file-type 1.
c
        itp = 1
        itp1 = 0
        write(20,300)itp,itp1
300        FORMAT(2(1x,i1))
c        write the third line. It should have the misc. information.
c
                write(20,400)NX,IYEAR,IMON,IDAY,IHOUR,IMIN,
     1           ISEC,IDT,L29,ICODE,IUNIT
400                FORMAT(1X,I6,1X,I4,1X,5(I2,1X),I9,1X,I4,2(1X,I1))
                IF (nprint.eq.1) THEN
                          WRITE(IUNITO,900)NX,IDAY,IMON,IYEAR,IHOUR,
     1                        IMIN,ISEC,IDT,L29
900        FormAT(/,5x,' This file contains ',I5,' values'/,5x,
     1                ' Starting from ',I2,'/',I2,'/',I4/,5x,
     2                '  at time ',I2,':',I2,':',I2/,5x,
     3              ' The time-interval is ',I9,' seconds.'/,5x,
     4               ' L29 is ',I4)
                        IF (IUNIT.GE.0.AND.IUNIT.LE.4) THEN
                WRITE(IUNITO,901)UNAME(IUNIT+1)
901                FORMAT(/,5x,' Unit of data is : ',a20/,5x)
                        ENDIF
                        IF (ICODE.EQ.1) THEN
                WRITE(IUNITO,910)
910                      FORMAT(6x,'These are averaged values.')
                        ELSE
                WRITE(IUNITO,920)
920                       FORMAT(' These are sampled values.')
                        ENDIF
                endif
                IF (NX.GT.IDIMX) THEN
c        Too many input values.
        if(idebug.ne.0)write(IUNITO,500)NX,IDIMX,IDIMX
500                        FORMAT(/,20x,' Warning - ffw2'/,20x,
     1               ' Data file contains ',I5,' values.'/,20x,
     2                ' This package is dimensioned for a max. of ',
     3                I5,/,20x,' Proceeding with ',I5,
     4               ' time-series values only')
                        NX=IDIMX
                        ier = 1
c
c        Non fatal error.
c
                ENDIF
c
c        Check that information is legal.
c
                IF(IUNIT.LT.0.OR.IUNIT.GT.4)THEN
        if(idebug.ne.0)write(IUNITO,902)IUNIT
902                FORMAT(/,20x,' Warning - ffw2, Illegal units '/,20x,
     1                        ' Unit code ',I6,' is not defined.')
                        IER=7
                ENDIF
                IF (IYEAR.LT.1800.OR.IYEAR.GT.2100) THEN
                if(idebug.ne.0)write(IUNITO,810)IYEAR
810                FORMAT(/,20x,' Warning - ffw2, Illegal year ',I4)
                        IER=2
                ENDIF
                IF(IMON.LT.1.OR.IMON.GT.12) THEN
                if(idebug.ne.0)write(IUNITO,820)IMON
820                FORMAT(/,20x,' Warning - ffw2, Illegal month',I4)
                        IER=3
                ELSE
                        NMON=MONS(IMON)
                        IF(IMON.EQ.2) THEN
                                IF(IYEAR.EQ.IYEAR/4*4) THEN
                                        NMON=29
                                ENDIF
                        ENDIF
                        IF(IDAY.LT.1.OR.IDAY.GT.NMON) THEN
                if(idebug.ne.0)write(IUNITO,930)IDAY,IMON,IYEAR
930                                FORMAT(/,20x,' Warning - ffw2 ',
     1                                 ' Illegal date ',
     2                                 I2,'/',I2,'/',I4)
                                IER=4
                        ENDIF
                ENDIF
                TIME=IHOUR+(IMIN+ISEC/60.)/60.
                IF(TIME.LT.0.0.OR.TIME.GT.24.) THEN
        if(idebug.ne.0)write(IUNITO,840)IHOUR,IMIN,ISEC,TIME
840              FORMAT(/,20x,' Warning - ffw2   Illegal time ',
     1                         I2,':',I2,'.',I2,5X,'(',F10.5,')')
                        IER=5
                ENDIF
c
                IF(IDT.LE.0) THEN
                if(idebug.ne.0)write(IUNITO,850)IDT
850                        FORMAT(/,20x,' Warning - ffw2  Illegal time',
     1                         ' increment ',I10)
                        IER=6
                ENDIF
c
c        Now write the data.
c
                        write(20,600)(X(I),I=1,NX)
600                        FORMAT(2X,6E13.5)
c
c        Sumarise the data as a check.
c
                XMEAN=X(1)
                XVAR=0.0
                XMIN=XMEAN
                XMAX=XMIN
                DO 10 I=2,NX
                TEMP=X(I)
                IF (TEMP.GT.XMAX) THEN
                        XMAX=TEMP
                ELSEIF(TEMP.LT.XMIN) THEN
                        XMIN=TEMP
                ENDIF
                XVAR=((I-2)*XVAR+(TEMP-XMEAN)**2*(I-1)/I)/(I-1)
                XMEAN=((I-1)*XMEAN+TEMP)/I
10                CONTINUE
        if(nprint.eq.1)   WRITE(IUNITO,1000)XMEAN,XVAR,XMIN,XMAX,NX
1000                FORMAT(/,5x,' Short Summary of Data'/,5x,
     1                ' Mean of Series ',t25,E13.6/,5x,
     2                ' Variance of Series ',t25,E13.6/,5x,
     3                ' Minimum value',t25,E13.6/,5x,
     4                ' Maximum value',t25,E13.6/,5x,
     5          /,5x,' No. of values used',t25,I6)
c
c
        RETURN
        END

c       This subroutine calculates the Gamma function impulse response,
c       sampled pulse response of duration D (Nash model) or the
c       averaged pulse response of duration D averaged over the
c	duration D (Nash model).
c
c
        Subroutine ffw42(h,d,m,an,ak,nopt,nacc,z,s,maxs,nprint,ier)
c
c        Dimension h(m),z(maxs),s(maxs)
	  Dimension h(8000),z(8000),s(8000)

        Common/syspar/iuniti,iunito,idebug
c
        if(nopt.lt.1.or.nopt.gt.3)nopt = 1
        if(idebug.ne.0)then
        Write(iunito,1000)
1000    Format(/,15x,'Debug information',/,15x,30('-'))
        write(iunito,1002)d,m,an,ak,nopt,nacc,maxs
1002    Format(/,15x,'Subroutine ffw42 called',/,15x,'D',t21,'=',
     1  e13.4,/,15x,'m',t21,'=',i4,/,15x,'an',t21,'=',e13.4,/,15x,
     2  'ak',t21,'=',e13.4,/,15x,'nopt',t21,'=',i1,/,15x,
     3  'nacc',t21,'=',i3,/,15x,'maxs',t21,'=',i5)
        endif
c
        do 2 i=1,m
         h(i)=0.0
2       continue
        do 3 i=1,maxs
          z(i)=0.0
          s(i)=0.0
3       continue
c
        if(nacc.lt.1)nacc=1
        if(d.lt.0.0)d=1.0
        delt=d/nacc
        nsamp = nacc*m
        if((nsamp+nacc+1).gt.maxs)then
        write(iunito,998)maxs
998     Format(/,20x,'Dimensional error in ffw42 - maxs = ',i5)
        ier = 999
        go to 999
        endif
        t = 0
        tbyk = delt/ak
c
	call gamma(an,x,ifail)
        if (ifail.ne.0) then
        write(iunito,1003) ifail
1003    format(5x,'Error detected in gamma - ier =',i3)
        ier = 998
        goto 999
        endif
c
        z(1) = 0.0
        const = 1.0/(x*ak)
        if(nopt.eq.1)then
        if(abs(an-1.0).le.0.00005)z(1) = const
        do 5 j = 2,nsamp
        t = t+tbyk
        z(j) = const*exp(-t)*t**(an-1.0)
5       continue
        go to 100
        endif
c
        ntotal = nsamp+nacc+1
        do 10 j = 2,ntotal
        t = t+tbyk
        iere = 0
        call incgam(t,an,x,ss,iere)
        if(iere.ne.0)then
        write(iunito,996)iere
996     Format(/,20x,'Error detected in incgam - ier =',i3)
        ier = 997
        go to 999
        endif
        s(j) = ss
10      continue
        do 15 j = 2,nacc+1
        z(j) = s(j)/d
15      continue
        do 16 j = nacc+2,nsamp
        z(j) = (s(j)-s(j-nacc))/d
16      continue
        if(nopt.eq.2)go to 100
        h(1) = 0.0
        do 50 j = 1,nacc+1
50      h(1) = h(1)+z(j)
        h(1) = h(1)/nacc
        nt = nacc+1
        do 20 j = 2,m
        h(j) = 0.0
        do 25 k = nt+1,nt+nacc
        h(j) = h(j)+z(k)
25      continue
        nt = nt+nacc
        h(j) = h(j)/nacc
20      Continue
        go to 60
100     h(1) = z(1)
        k = nacc+1
        do 30 j = 2,m
        if(k.gt.nsamp)go to 60
        h(j) = z(k)
        k = k+nacc
30      Continue
60      if(nprint.ne.0)then
        write(iunito,520)d,delt,an,ak,m
520     Format(/,10x,'TABLE',/,10x,20('-'),//
     1  /,10x,'Method : NASH MODEL',
     2  /,10x,'Data : Pulsed / Sampled output and pulsed input',
     3  /,10x,'Duration',t30,'=',f10.2,
     4  /,10x,'Sampling interval',t30,'=',f10.2,
     5  /,10x,'Model parameters',/,15x,'n',t30,'=',f10.2,
     6  /,15x,'k',t30,'=',f10.2,
     7  /,10x,'No. of ordinates',t30,'=',i3,//,
     8  /,10x,'Catchment : ',//)
        if(nopt.eq.1)write(iunito,500)
        if(nopt.eq.2)write(iunito,505)
        if(nopt.eq.3)write(iunito,510)
500     Format(10x,'The sampled impulse response function',
     1  /,10x,38('-'))
505     Format(10x,'The sampled pulse response function',
     1  /,10x,36('-'))
510     Format(10x,'The averaged pulse response function',
     1  /,10x,37('-'))
        write(iunito,525)(h(i),i=1,m)
525     Format(10x,6f10.3)
        endif
999     return
        end
c
c
C================================================================
      SUBROUTINE INCGAM(X,P,GP,D6,IFAULT)
C
C     ADAPTATION OF MOORE'S ALG AS 187 IN APPL STAT, FOR COMPUTING
C     INCOMPLETE GAMMA INTEGRAL & ITS DERIVATIVES. THIS ADAPTATION
C     ONLY USES THOSE PORTIONS PERTAINING TO THE INCOMPLETE GAMMA
C     INTEGRAL AND DOES NOT INCLUDE THE DERIVATIVES.
C
C     X - ARGUEMENT OF INC GAMMA FN
C     P - SHAPE FACTOR OF INC GAMMA FN
C     GP- GAMMA FN OF P
C     D6- INC GAMMA FN OF X, FOR SHAPE FACTOR P
C     IFAULT- ERROR CODE
C
C
C     THE FOLLOWING CHANGES HAVE BEEN MADE IN THIS ADAPTATION:
C     1) AS 187 RESULTS WERE STORED IN VECTOR D(6); RESULTS HERE
C        FOR INCOMPLETE GAMMA INTEGRAL ARE STORED IN SCALAR D6;
C     2) EXPRESSIONS FOR VARIABLE "F" HAVE BEEN REPLACED WITH
C        EXACT EQUIVALENT EXPRESSIONS;
C     3)
C
C==============================================================
      DIMENSION PN(6)
      DATA E,OFLO,TMAX,ZERO/1.0E-6,1.0E30,100.0,1.0E-30/
      IFAULT=0
      PM1=P-1.0
C     XLOG=ALOG(X)
      IF (X.GT.1.0.AND.X.GE.P) GOTO 30
C     F=EXP(P*XLOG-GP1LOG-X)
      F=(X**P)*EXP(-X)/(P*GP)
      TMAXP=TMAX-P
      C=1.0
      S=1.0
      A=P
    1 A=A+1.0
      C=C*X/A
      S=S+C
      IF (A.GT.TMAXP) GOTO 1001
      IF (C.GT.E*S) GOTO 1
      D6=S*F
      RETURN
C
C        CONTINUED FRACTION EXPANSION
C
   30 F=(X**P)*EXP(-X)/GP
C
      A=PM1
      B=X+1.0-A
      TERM=0.0
      PN(1)=1.0
      PN(2)=X
      PN(3)=X+1.0
      PN(4)=X*B
      S0=PN(3)/PN(4)
   32 A=A-1.0
      B=B+2.0
      TERM=TERM+1.0
      AN=A*TERM
      PN(5)=B*PN(3)+AN*PN(1)
      PN(6)=B*PN(4)+AN*PN(2)
C
      IF (ABS(PN(6)).LT.ZERO) GOTO 35
      S=PN(5)/PN(6)
      C=ABS(S-S0)
      IF (C*P.GT.E) GOTO 34
      IF (C.LE.E*S) GOTO 42
C
   34 S0=S
   35 DO 36 I=1,4
      I2=I+2
      PN(I)=PN(I2)
   36 CONTINUE
C
      IF (TERM.GT.TMAX) GOTO 1001
      IF (ABS(PN(5)).LT.OFLO) GOTO 32
      DO 41 I=1,4
      PN(I)=PN(I)/OFLO
   41 CONTINUE
      GOTO 32
C
   42 D6=1.0-F*S
      RETURN
C
C         SET FAULT INDICATOR
C
 1001 IFAULT=1
      RETURN
      END
c
c
	SUBROUTINE GAMMA(X,GAM,IER)
C- Calculates the Gamma function for positive arguments.
C- M.Bruen, Hydrology Dept., UCG, 1983
C- Cf. Hart,J.F. et al., "Computer Approximations", Wiley, 1968
C-
C- The calling arguments are
C-
C- 	X	The argument of the Gamma function.
C-		 It must be a positive number less than 34.8
C-	GAM	The resulting value of the Gamma function.
C-	IER	An error parameter which indicates if this
C-		 routine is successful, IER=0, or fails,
C-		 IER > 0
C-
	IER=999
	IF(X.LT.0.0)RETURN
	IF(X.GT.34.8)RETURN
C- THIS GIVES A GAMMA WHICH IS JUST LESS THAN THE MAXIMUM
C- WHICH CAN BE STORED IN THE DEC-20
	IER=0
	IF(X.LE.20.0) GOTO 10
C- STIRLINGS FORM
	Y=1.0/(X*X)
	P=(0.77783067E-3*Y-0.277765545E-2)*Y+0.8333333309E-1
	P=P/X
	GAM=(X-0.5)*ALOG(X)-X+0.9189385+P
	GAM=EXP(GAM)
	RETURN
10	Y=AINT(X)
	N=Y-2.
	Y=X-Y
	GAM=(((.1082985985E-1*Y-.3427052255E-2)*Y+.77549276E-1)*Y)
	GAM=(((GAM+.8017824769E-1)*Y+.4121029027)*Y+.4227663678)*Y
	GAM=GAM+1.000000199
	T1=1.0
	YP2=Y+2.0
	IF(N)40,70,60
40	CONTINUE
C- NEGATIVE N
	N=IABS(N)
	DO 45 I=1,N
45	T1=T1*(YP2-I)
	T1=1.0/T1
	GOTO 70
60	CONTINUE
C- POSITIVE N
	N=N-1
	DO 65 I=0,N
65	T1=T1*(YP2+I)
70	GAM=GAM*T1
	RETURN
	END


	 