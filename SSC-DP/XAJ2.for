      subroutine xaj(r,e,Qhat,consi,par,consit,
     $               nmax,npar,ncon,area,delta)


	implicit none

	  integer nmax,npar,ncon,m
        parameter(m=60)
c       maxw > nacc*(m+2)  (nacc=24)
        real r(nmax),e(nmax),qt(nmax),qi(nmax),qg(nmax),q(nmax)
        real par(npar),qhat(nmax),qhat0(m+nmax),qhat00(nmax)
        real h(200),consi(ncon),consit(ncon),par0(npar)
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

!         do 1 i=1,ncon
!1        consit(i)=consi(i)


c
       
c	do 2 i=1,npar
c           par0(i)=par(i)
c           if(par0(i).lt.parmin(i)) par0(i)=parmin(i)
c           if(par0(i).gt.parmax(i)) par0(i)=parmax(i)
c2       continue



        call GR4J_MOD(1,nmax,par,consi,consit,r,e,q,npar,ncon)




		do  20 i=1,nmax
	qhat(i)=q(i)*area/(3.6*delta*24)
20      continue





	


	

		

c
         return
         end