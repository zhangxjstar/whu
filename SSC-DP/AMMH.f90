!!*****************************
!!Monte Carlo Markov Chain
!!*****************************
	subroutine AMMH(func,npara,nobj,nistat,itermax,xmin,xmax,seqout,objout,bestparax,stateq,parx)
	USE numerical_libraries
	implicit none
	


	integer npara,itermax,nobj,nistat
!   	integer,parameter :: itermax=50000  !原来采用参数设置该值，发现可能出现数组出界但不报错。
	real :: seqout(itermax,npara),objout(itermax,nobj)
	
	integer :: L  
    real :: z
	REAL::xmin(npara),xmax(npara)
	real :: cov(npara,npara),R(1,npara),ave(npara),RSIG(npara,npara),parx(6)
	real :: ofsp(npara), seq(npara,itermax), seqlike(itermax),bestparax(npara),stateq(itermax,nistat),stateq1(nistat),seqsum(npara),seqsum1(npara,npara),seqsum2(npara)

	real :: funoff  
	integer :: i,j,k,item,iter
	integer :: IRANK, LDR, LDRSIG, NOUT, NR
	real :: rand, temp, obj(nobj) !, RSIG(m2,m2)
	real :: sd
	real,parameter :: omiga = 1.0E-5
	logical :: flag


!	call random_number(rand)		!生成0.----1.的随机数 rand(nPara)
	CALL SEED (12358)

	
	do i=1, 1
	!	seq(i,1) = xmin(i) + rand()*(xmax(i)-xmin(i)) 

seq(i,1)=bestparax(i)

!write(*,*) seq(i,1)
	enddo


	do j=1,npara
	seqsum(j)=0.0

	seqsum2(j)=0.0
	enddo

	seqsum1=0.0



	call func(seq(:,1),obj,stateq1,parx)	!调用目标函数
	seqlike(1) = obj(1)

	sd=2.4**2.0/(npara)
	CALL ERSET (3, 0, 0)  !设置不输出警告信息
	do iter=1,itermax-1
		!write(*,*) iter



		if(iter > 100)then
			!**************************
			!计算协方差矩阵
			!**************************	
	do j=1,npara
				ave(j) = sum(seq(j,1:iter))/iter
				
			enddo

			cov=0
			do j=1, npara
				do k=1, npara
					if(j /= k)then
						temp = 0
					else
						temp = 1
					endif

			     cov(k,j) = sum((seq(j,1:iter) - ave(j))*(seq(k,1:iter)-ave(k)))

				cov(k,j) = sd*cov(k,j)/iter + sd*omiga*temp
				enddo			
			enddo						
		else
			do j=1, npara
				do k=1, npara
			      cov(k,j) = rand()*(xmax(j)-xmin(j))*(xmax(k)-xmin(k))*1.0 !1.0E   !!!change 
				enddo			
			enddo	
		endif


		
		do while(.true.)
			do while(.true.)       !利用多元正态分布生成满足决策变量可行域的的新样本
				CALL UMACH (2, NOUT)
				NR =1
				K = npara
				LDRSIG =npara
				LDR =1
				! Obtain the Cholesky factorization.
				CALL CHFAC (K,COV,npara,0.001,irank,RSIG,LDRSIG)

				CALL RNMVN (NR,K,RSIG, LDRSIG, R, LDR)
				ofsp=seq(:,iter)+r(1,:)

				!write(*,*) iter,ofsp(1:2)

				flag = .true.
				do i=1,npara
					if(ofsp(i) > xmax(i) .or. ofsp(i) < xmin(i)) then
						flag = .false.
						exit
					endif
				enddo
				
	!			if(iter < 100 .or. flag) exit
				if(flag) exit

			!	exit
			enddo

			call func(ofsp,obj,stateq1,parx)
			funoff = obj(1)
		!	if(obj(4) <= 0.0) exit
			exit
		enddo
					
		temp = min(funoff,seqlike(iter)) !为负值需要处理
		if(temp < 0.0) then
			if(funoff > seqlike(iter))then
				z = 1.0
			elseif(funoff * seqlike(iter) > 0 )then
				z=seqlike(iter)/funoff
			else
				z = 0.0
			endif
		else
			z=funoff/seqlike(iter)  !liki(item2,item)
		endif


!write(*,*) iter,z,funoff,seqlike(iter)

if(z.ge.1) then
			seq(:,iter+1)=ofsp	
			seqlike(iter+1) = funoff
!				seq(l,iter+1)=ofsp(l)	
		else
			if(z.ge.(rand()+0.2)) then
			seq(:,iter+1)=ofsp	
			seqlike(iter+1) = funoff
		else

			seq(:,iter+1)=seq(:,iter)	
			seqlike(iter+1) = seqlike(iter)
			endif
!				seq(l,iter+1)=seq(l,iter)	
		end if
!		end do

!		call func(seq(:,iter+1),obj)
!		seqlike(iter+1) = obj(1)
!		write(100001,"(1x,200f18.2)")1*log(seqlike(iter+1))+800,seq(:,iter+1)

!		write(*,*)iter+1,seqlike(iter+1)
	end do
   				!	seq(1,1)=233.1052
				!	seq(2,1)=0.2796 
				!	seq(3,1)=0.3206
				!	seq(4,1)=0.6888 
				!	seq(5,1)=4.8905
				!	seq(6,1)=25.5421
				!	seq(7,1)=1.8386 
				!	seq(8,1)=0.1226
				!	seq(9,1)=0.3344
				!	seq(10,1)=0.1037 
				!	seq(11,1)=0.3309 
				!	seq(12,1)=0.9484
				!	seq(13,1)=0.9883
				!	seq(14,1)=5.178
				!	seq(15,1)=8.8354

!	open(100,file='AMMH2.dat')
	seqout=transpose(seq)
	
!	open(12,file="W.txt")
	do i=1,itermax
		call func(seq(:,i),objout(i,:),stateq(i,:),parx)
		!write(12,'(I8, 1X, 2f16.3)') i,seq(1,i),objout(i,1)
!		write(100,"(1x,200f12.3)")seqlike(i),seq(:,i)
	enddo
!	write(*,*) "ok"
	close(100)
				
!				
				
!	return
end subroutine AMMH