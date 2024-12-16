		! nsplit----number of segementations
		! nfeasible----number of parameters set

module define		
		integer nsplit,nfeasible,nparameter,nobjective,nistate,gni,ns
		integer ndata,ndatamax,nyear0,i_file
		real area,delta
		character*7 filenamex(3000)
		character*20 catchm,fyle,filename,fyout,fyresult,fyle_dc,fyle_wbi,fydata,fydata1,fydata2,fydata3


!		parameter(ndatamax=80000,area=10860.0,delta=1.0)  !lhsk
!		parameter(nsplit=11,nfeasible=10000,nparameter=5,nobjective=1,nistate=1,length=120) !lhsk&XAJ

		parameter(ndatamax=80000,delta=1.0,ns=54)  !lhsk
		parameter(nfeasible=5000,nparameter=1,nobjective=1,nistate=3) !lhsk&XAJ
		
		real para(ns,nfeasible,nparameter)
		real obj(ns,nfeasible,nobjective)
		real paramin(nparameter),paramax(nparameter)

		!area,nsplit,lengthĞèÒª¶ÁÈ¡

		integer nbegin(ns),nend(ns),length(ns)


		real rain(ndatamax),panev(ndatamax),Qobs(ndatamax),Qsim(ndatamax),qhat0(ndatamax),qhat00(ndatamax),Tave(ndatamax)
		real istate(ns+1,nistate),bestpara(ns,nparameter),bestobj(nobjective),bep(ns,nparameter)
		real qout(ndatamax),swobs(ndatamax)
		real sw(61),skr(61),nobs(74),vobs(74)

		data paramin/1.0/
		data paramax/1000.0/



end module define