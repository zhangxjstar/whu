		SUBROUTINE ENKF(xp, xa, yob, n, m,  l, funcc)

		USE numerical_libraries
		IMPLICIT NONE

		! THIS PROGRAM PERFORM ENSEMBLE KALMAN FILTERING FOR MODELED SURFACE SOIL MOISTURE AND
		! SATELLITE-BASED REMOTELY SENSED SURFACE SOIL MOISTURE. 
		!
		! THE VARIABLES USED ARE DESCRIBED FOR CLARITY
		! 
		! xf     : forecast of state varibales 
		! xa     : anaylsis of state varibales 
		! n      : number of states
		! m      : number of ensemble members
		! yob    : observation of Y
		! obs_cov: error std in xob
		! Pa     : analysis/posterior error covaiance of state variables 
		! Pf     : forecast/prior error covaiance of state variables
		!
		!
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! Change log
		!
		!
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!
		! PAN LIU (12-21-2013)
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		!
		! Declare variables (in and out variables)

		INTEGER, INTENT(IN)::             n, m, l
		REAL*4, INTENT(IN)::              xp(n,m)  !状态转移方程得到的下一时段状态
		REAL*4, INTENT(IN)::              yob(l)
		REAL*4, INTENT(OUT)::             xa(n,m)  !状态修正值

!		REAL*4, INTENT(OUT)::             Pa(n,n) Pf(n,n),
!		REAL*4, INTENT(OUT)::             obs_cov

		!Declare local variables
		REAL*4                            R(l,l),err(l),siga(m,l),HI(l,m),D(l,m)
		REAL*4                            temp(n,l),temp1(n,m),RSIG(l,l)
		REAL*4                            xmean(n),ymean(l),XK(n,m),YK(l,m),PH(n,l),HPH(l,l),A(l,l),AINV(l,l)
		INTEGER*4                         ii, jj, kk
!		external funcc
		integer :: IRANK, LDR, LDRSIG, NOUT, NR,K

		
 		err=0.00 !1.0E-8

		R=0.0
		do ii=1,l
			R(ii,ii)=yob(ii)*err(ii)
 		enddo
! 
		CALL UMACH (2, NOUT)
		NR =m
		K = l
		LDRSIG =l
		LDR =m
		! Obtain the Cholesky factorization.
		CALL CHFAC (K,R,l,0.001,irank,RSIG,LDRSIG)
		CALL RNMVN (NR,K,RSIG, LDRSIG, siga, LDR)  !生成多维正态分布伪随机数


		do ii = 1, m
			call funcc(xp(:,ii),HI(:,ii),n,l)
			D(:,ii)=yob(:)-HI(:,ii)+siga(ii,:)
		enddo

		do ii=1,n
			xmean(ii)=sum(xp(ii,1:m))/m
		enddo

		do ii=1,l
			ymean(ii)=sum(HI(ii,1:m))/m
		enddo


		do ii=1,m	
			XK(:,ii)=xp(:,ii)-xmean(:)
			YK(:,ii)=HI(:,ii)-ymean(:)
		enddo


		PH = (MATMUL (XK, TRANSPOSE(YK)))/(m-1)
		HPH= (MATMUL (YK, TRANSPOSE(YK)))/(m-1)

		A=(HPH + R)
		CALL LINRG (l, A, l, AINV, l)   !矩阵求逆
		
		temp=MATMUL(PH,AINV)
		temp1=MATMUL(temp,D)

		xa=xp+temp1


		END SUBROUTINE ENKF




