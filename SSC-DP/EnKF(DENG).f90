		SUBROUTINE ENKF(xp, xa, yob, n, m,  l, sq)

!		use global
!		USE numerical_libraries
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
		REAL, INTENT(IN)::              xp(n,m)  !状态转移方程得到的下一时段状态
		REAL, INTENT(IN)::              yob(l)
		REAL, INTENT(IN)::              sq(l,m)
		REAL, INTENT(OUT)::             xa(n,m)  !状态修正值

!		REAL*4, INTENT(OUT)::             Pa(n,n) Pf(n,n),
!		REAL*4, INTENT(OUT)::             obs_cov

		!Declare local variables
		REAL                            R(l,l),err(l),siga(m,l),HI(l,m),D(l,m)
		REAL                            temp(n,l),temp1(n,m),RSIG(l,l)
		REAL                            xmean(n),ymean(l),XK(n,m),YK(l,m),PH(n,l),HPH(l,l),A(l,l),AINV(l,l)
		INTEGER                         ii, jj, kk
		integer :: IRANK, LDR, LDRSIG, NOUT, NR,K

!--------------------------		
		HI=sq
!--------------------------
		err=3.0E-1

		R=0.0
		do ii=1,l
			R(ii,ii)=yob(ii)*err(ii)
		enddo

		CALL UMACH (2, NOUT)
		NR =m
		K = l
		LDRSIG =l
		LDR =m
		! Obtain the Cholesky factorization.
		CALL CHFAC (K,R,l,0.001,irank,RSIG,LDRSIG)
		CALL RNMVN (NR,K,RSIG, LDRSIG, siga, LDR)  !生成多维正态分布伪随机数


		do ii = 1, m
!			call funcc(xp(:,ii),HI(:,ii),n,l)
!			D(:,ii)=yob(:)+siga(ii,:)-HI(:,ii)     !original
			D(:,ii)=yob(:)+0.05*yob(:)*rnnof()-HI(:,ii)      !加入流量观测值采用取乘以观测值的百分数形式
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
		
		temp=MATMUL(PH,AINV)    ! 增益矩阵
		temp1=MATMUL(temp,D)

		xa=xp+temp1    ! updated states
!----------------------------------------参数范围限制
! 	do jj=1,m
! 	  if (xa(1,jj)<0.0) then
! 	      xa(1,jj)=0.0
!  	  elseif (xa(1,jj)>xmax(1)) then
!  	      xa(1,jj)=xmax(1)
! 	  endif
! 
! 	  if (xa(2,jj)<0.0) then
! 	      xa(2,jj)=0.0
!  	  elseif (xa(2,jj)>xmax(2)) then
!  	      xa(2,jj)=xmax(2)
! 	  endif
! 
! 	  if (xa(3,jj)<0.0) then
! 	      xa(3,jj)=0.0
! 	  elseif (xa(3,jj)>xa(2,jj)) then
! 	      xa(3,jj)=xa(2,jj)
! 	  endif
!      enddo
!-----------------------------
!    do jj=1,m
!        if ( xa(1,jj)<xmin(1) .or. xa(1,jj)>xmax(1) .or. xa(2,jj)<xmin(2) .or. xa(2,jj)>xmax(2) ) then
! 	     
! 
!    enddo

		END SUBROUTINE ENKF


