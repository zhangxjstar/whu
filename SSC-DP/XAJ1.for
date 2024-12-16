      subroutine xaj1(r,e,Qhat,consi,par,consit,
     $                 nmax,npar,ncon,area,delta,ni
     &,qout,swobs,sw,skr,nobs,vobs)
	implicit none

	  integer nmax,npar,ncon,m
        parameter(m=60)
c       maxw > nacc*(m+2)  (nacc=24)
        real r(nmax),e(nmax),qt(nmax),qi(nmax),qg(nmax)
        real par(npar),qhat(nmax)
        real h(200),consi(ncon),consit(ncon)
	  real area,delta
		real  qout(80000),sw(80000),skr(nmax),swobs(80000),swhat(nmax)
	real  v1(nmax),v(nmax),vxx(nmax),vobs(80000),nobs(1000)
	  real cn,ck,d,g1,g2,kr,sc,gama,hw(nmax+1),hw0
	real se(nmax),sqs(nmax),sqg(nmax),sq(nmax),rain(nmax),panev(nmax)
	  integer i,j,nopt,nacc,maxs,nprint,ier,ni
	integer am,jv,inti,n2
!       common /rain/r,/evap/e,/disc/qobs,/qhat/qhat,
!    &	  /consi/consi,/uh/h  !,QWE(nmax)
!       common /number/n1,n2,nwarm,m, /switch/isw
!       common /bndlow/parmin(npar),/bndup/parmax(npar)
!	  common /area/area,delta

         do 1 i=1,ncon
1        consit(i)=consi(i)


c
       
	

	 g1=par(1)
      g2=par(2)
      kr=par(3)
      sc=par(4)
      gama=par(5)
      hw(1)=consit(1)

      do i=1,nmax
	rain(i)=r(i)
	panev(i)=e(i)

 ! simulated evaporation
        se(i)=g1*e(i)*TANH(r(i)/e(i))
   ! simulated discharge
        sq(i)=( hw(i)+(r(i)-se(i)) )*TANH( (hw(i)+(r(i)-se(i)))/g2 )
   ! simulated water depth
       IF(i<nmax)  then
	     hw(i+1)=hw(i)+(r(i)-(sq(i)+se(i)))
       ENDIF
	   
	IF(hw(i+1)<EPSILON(1.)) hw(i+1)=EPSILON(1.) 


			      ! simulated evaporation
!     hw0=hw(i)+rain(i)
!     se(i)=panev(i)*((hw(i)/sc)**gama)
      ! simulated discharge
!      sqs(i)=g1*((hw(i)/sc)**g2)*rain(i)
!      sqg(i)=kr*(hw(i)+hw0)/2
!      sq(i)=sqs(i)+sqg(i)
      ! simulated water depth
!     IF(i<nmax)  hw(i+1)=hw0-(sq(i)+se(i))
!      IF(hw(i+1)<EPSILON(1.)) hw(i+1)=EPSILON(1.)
      	
		Qhat(i)=sq(i)
      enddo

!	write(*,*) ni,1,qhat(1),swhat(1),swobs(1+720*(ni-1)),vx(1)
c
         return
         end