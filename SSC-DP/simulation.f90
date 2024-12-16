		module simulation
		use define
		use imsl

		implicit none

		contains

			subroutine hydromodel(ni,np1,np2,listate,lpar,laststate,v,parx)
			implicit none

			integer i,np1,np2,ni
			real listate(nistate),lpar(nparameter),laststate(nistate)
			real qhat0(ndatamax),qhat00(ndatamax),v(32),parx(6)
			real c,sc,hw,s
			
			goto 50
!////////////////////////two parameter water banlance model///////////////////////////
			c=lpar(1)
			sc=lpar(2)
			hw=listate(1)			

			do i=np1,np2
				call TWmodel(rain(i),panev(i),c,sc,hw,Qsim(i),s) 
				hw=s
			enddo
			laststate(1)=hw

			return

	50		continue
!write(*,*) "listate1",listate

			
			call xaj(ni,rain(np1:np2),panev(np1:np2),Qsim(np1:np2),listate,lpar,laststate,np2-np1+1,nparameter,nistate,area,delta,v,parx,Tave(np1:np2)) 

			      
!write(*,*) "listate2",listate
!write(*,*) "laststate",laststate(1)

			end subroutine hydromodel



			subroutine hydromodel1(ni,np1,np2,listate,lpar,laststate,v,parx)
			implicit none

			integer i,np1,np2,ni
			real listate(nistate),lpar(nparameter),laststate(nistate),v(32),parx(6)
			real c,sc,hw,s
			
			goto 50
!////////////////////////two parameter water banlance model///////////////////////////
			c=lpar(1)
			sc=lpar(2)
			hw=listate(1)			

			do i=np1,np2
				call TWmodel(rain(i),panev(i),c,sc,hw,Qsim(i),s) 
				hw=s
			enddo
			laststate(1)=hw

			return

	50		continue
!write(*,*) "listate1",listate


			call xaj(ni,rain(np1:np2),panev(np1:np2),Qsim(np1:np2),listate,lpar,laststate,np2-np1+1,nparameter,nistate,area,delta,v,parx,Tave(np1:np2)) 
!write(*,*) "listate2",listate
!write(*,*) "laststate",laststate(1)

			end subroutine hydromodel1


!////////////////////////////////////////////////////////////
			subroutine runDDP()
			implicit none

			integer i,j,k
			integer ia(nsplit,nfeasible),nbest(nsplit)
			real temp,uu
			real mindis(nfeasible),tempdis(nfeasible)

			
			
			do i = 1, nsplit
				call generation(i,para(i,:,:),obj(i,:,:))  !nfeasible,nobjective,
			!write(*,*) i
			!	write(*,*) i, obj(i,1:100,1)
			!	write(*,*) para(i,1:100,:)
			enddo
		
		!	write(*,*) obj(1:20,1,1)
			temp=0.0
			mindis = 0.0
			do i = 1, nsplit - 1
				tempdis = 1.0E15
				do k = 1,nfeasible
					do j = 1,nfeasible
					
							temp = mindis(j)+50 - (obj(i,j,1) + obj(i+1,k,1))*1 + dist(para(i,j,:),para(i+1,k,:))*0.01
						
						if(temp < tempdis(k)) then
							tempdis(k) = temp
							ia(i,k)=j
						endif
					enddo

				enddo
				mindis = tempdis
			enddo
			
			nbest(nsplit) = 1
			temp = mindis(1)
			do i = 2,nfeasible
				if(mindis(i) < temp) then
					nbest(nsplit) = i
					temp = mindis(i)
				endif
			enddo

			do i = nsplit-1, 1, -1
				k = nbest(i + 1)
				nbest(i)=ia(i,k)
			enddo

			do i=1,nsplit
				bestpara(i,:) = para(i,nbest(i),:)
				write(*,"(7f9.4)") bestpara(i,:)
			enddo

	!		call simulate(1,nsplit)


			end subroutine runDDP
				

!/////////////////////////////////////////////////////////////////////////////////
			subroutine generation(ni,parai,obji)  !nfea,nobj,
			implicit none

			integer i,j,ni,ntest,h   !,nfea,nobj
			parameter(ntest=nfeasible*1)
			real parai(nfeasible,nparameter),obji(nfeasible,nobjective),parx(6)
			real aparai(nsplit,nfeasible,nparameter),paratex(nfeasible,nparameter)
			real objtemp(ntest,nobjective),paratemp(ntest,nparameter),laststate(nistate),iistaten(ntest,nistate)
			real istaten(ntest,nistate),tempn(nistate),aistate(nsplit,nfeasible,nistate),aiistate(nsplit,nfeasible,nistate)
			real tempp(nparameter),tempo(nobjective),tempm(nistate),bestpara0(nparameter),bestpara01(nsplit,nparameter)
			real rand
	
	gni=ni
				
		open(12,file=fydata2)

do i=1,6

				read(12,*) parx(i)



			!	write(*,*)  parx(i)
	enddo
       close(12)

do i=1,nparameter

bestpara0(i)=parx(i)

enddo

		

!	X0(1)=59.425
!	X0(2)=-6.338
!	X0(3)=98.425
!	X0(4)=2.055
!	X0(5)=0.181

				!write(*,*) bestpara0
			if (gni==1) then
				call AMMH(simulate,nparameter,nobjective,nistate,ntest,paramin,paramax,paratemp,objtemp,bestpara0,istaten,parx)

			end if

			if (gni>1) then
				!call AMMH(simulate,nparameter,nobjective,nistate,ntest,paramin,paramax,paratemp,objtemp,bestpara(gni-1,:),istaten)
				call AMMH(simulate,nparameter,nobjective,nistate,ntest,paramin,paramax,paratemp,objtemp,bestpara0,istaten,parx)
			end if
			
			

			do i=1,nfeasible
				do j=i+1,ntest
					if (objtemp(i,1) < objtemp(j,1)) then
						tempp=paratemp(i,:)
						tempo=objtemp(i,:)
						tempn=istaten(i,:)
						!tempm=iistaten(i,:)

						paratemp(i,:)=paratemp(j,:)
						objtemp(i,:)=objtemp(j,:)
						istaten(i,:)=istaten(j,:)
						!iistaten(i,:)=iistaten(j,:)

						
						paratemp(j,:)=tempp
						objtemp(j,:)=tempo
						istaten(j,:)=tempn
						!iistaten(i,:)=tempm

					endif
				enddo
			enddo			

			parai(1:nfeasible,:)=paratemp(1:nfeasible,:)
			aparai(ni,1:nfeasible,:)=paratemp(1:nfeasible,:)
			obji(1:nfeasible,:)=objtemp(1:nfeasible,:)
			bestpara(gni,:)=parai(1,:)
			!aistate(ni,1:nfeasible,:)=istaten(1:nfeasible,:)
			!aiistate(ni,1:nfeasible,:)=iistaten(1:nfeasible,:)
			
		!	do i=1,nistate
		!	istate(ni+1,i)=istaten(1,i)
		!	end do
				!write(12,*) gni
				write(*,*) gni
	!	write(*,*) "objtemp",objtemp(1:nfeasible,1)
			!write(*,*) "iistaten",iistaten(1:nfeasible,1)
			!write(*,"(6f9.2)") istate

			

!			do i=1,nfeasible
!			write(12,"(15f9.4)") paratemp(i,1)
!			enddo
			
			!write(12,*) "objtemp",objtemp(1:100,1)			
!		write(*,*) ni,obji(1:10,1)

		!call hydromodel(nbegin(ni),nend(ni),istate(ni,:),paratemp(1,:),istate(ni+1,:))
		
			end subroutine generation			



!/////////////////////////////////////////////////////////////////////////////////
			subroutine simulate(paratemp,objtemp,laststate,parx)
			implicit none

			integer ni
			real paratemp(nparameter),objtemp(nobjective),laststate(nistate),v(32),parx(6)			
			
			ni=gni

			v=0
			call hydromodel(ni,nbegin(ni),nend(ni),istate(ni,:),paratemp,laststate,v,parx)	
			
		!	write(*,*)			istate(ni,:)
			call sta(nbegin(ni),nend(ni),objtemp)		
		
			end subroutine simulate		

!/////////////////////////////////////////////////////////////////////////////////
			subroutine sta(np1,np2,objective)
			implicit none

			integer np1,np2,i
			real nse,balance,sum1,sum2,avg,avg1,nse2,r,a,b,sum3
			real objective(nobjective)
			sum1=0.0

			sum2=0.0
			sum3=0.0
			!评价指标改动
			avg=sum(Qsim(np1:np2))/(np2-np1+1)

			avg1=sum(Qsim(np1:np2))/(np2-np1+1)

			do i=np1,np2
				sum1=sum1+(Qsim(i)-Qobs(i))**2
				sum2=sum2+(Qobs(i)-avg)**2
	
			enddo
			
			objective(1)=1-sum1/sum2
			
			sum1=0.0

			sum2=0.0


			do i=np1,np2
				sum1=sum1+abs(Qobs(i)-Qsim(i))
				sum2=sum2+abs(Qobs(i)-avg)
	
			enddo
			
			objective(2)=1-sum1/sum2


			sum1=0.0

			sum2=0.0


			do i=np1,np2
			if(Qobs(i)>0.and.(Qsim(i)>0)) then
				sum1=sum1+(log(Qsim(i))-log(Qobs(i)))**2
				sum2=sum2+(log(Qobs(i))-log(avg))**2
				endif
	
			enddo
			
			objective(3)=1-sum1/sum2

		
			end subroutine sta

!/////////////////////////////////////////////////////////////////////////////////


!/////////////////////////////////////////////////////////////////////////////////
			subroutine sta1(np1,np2,objective)
			implicit none

			integer np1,np2,i,j
			real nse,balance,sum1,sum2,avg,avg1,nse2,r,a,b,sum3
			real objective(nobjective)
		sum1=0.0

!文件路径1
open(227,file="./qqh.txt")

			sum2=0.0
			sum3=0.0
			!评价指标改动
			avg=sum(Qobs(np1:np2))/(np2-np1+1)
			avg1=sum(Qsim(np1:np2))/(np2-np1+1)

		

do i=1,nsplit-1

do j=nend(i)+1,nend(i)+length(i+1)

!write(*,*) j

if(j/=nend(i)+1) then
sum1=sum1+(Qsim(j)-Qobs(j))**2
sum2=sum2+(Qobs(j)-avg)**2
write(227,*) j,Qobs(j),Qsim(j),avg
endif
			
	
	enddo
	enddo

	!		close(227)
			
			objective(1)=1-sum1/sum2
			
			sum1=0.0

			sum2=0.0


			do i=np1,np2
				sum1=sum1+abs(Qobs(i)-Qsim(i))
				sum2=sum2+abs(Qobs(i)-avg)
	
			enddo
			
			objective(2)=1-sum1/sum2



			sum1=0.0

			sum2=0.0


			do i=np1,np2
			if(Qobs(i)>0.and.(Qsim(i)>0)) then
				sum1=sum1+(log(Qsim(i))-log(Qobs(i)))**2
				sum2=sum2+(log(Qobs(i))-log(avg))**2
				endif
	
			enddo
			
			objective(3)=1-sum1/sum2



		
			end subroutine sta1

!/////////////////////////////////////////////////////////////////////////////////
			subroutine evaluate()
			implicit none
			integer i
			real laststate(nistate),v(32),parx(6)
		
		open(12,file=fydata2)

do i=1,6

				read(12,*) parx(i)



				!write(*,*)  bestpara0(i)
	enddo
       close(12)

	!	open(59,file="Wz.txt")

		v=0
			do i=1,nsplit
				call hydromodel(i,nbegin(i),nend(i),istate(i,:),bestpara(i,:),istate(i+1,:),v,parx)

				
		!write(59,*) i,istate(i,1)
			!	write(*,"(I5,6f9.2)") i,istate(i,:)
			!	write(*,"(I5,6f9.2)") i,istate(i+1,:)
			enddo
			call sta(nbegin(2),ndata,bestobj)

			write(*,*) bestobj

			end subroutine evaluate	

!/////////////////////////////////////////////////////////////////////////////////
			subroutine evaluate1()
			implicit none
			integer i
			real laststate(nistate),laststate1(nistate),v(32),parx(6)

open(12,file=fydata2)
do i=1,6

				read(12,*) parx(i)



				!write(*,*)  bestpara0(i)
	enddo
       close(12)
		v=0
			do i=1,nsplit
				call hydromodel(i,nbegin(i),nend(i),istate(i,:),bestpara(i,:),istate(i+1,:),v,parx)

				
	write(*,*) i,istate(i,:)
			enddo
			call sta1(nbegin(2),ndata,bestobj)

			write(*,*) bestobj

			end subroutine evaluate1
	
!/////////////////////////////////////////////////////////////////////////////////
			subroutine initial()
			implicit none
			integer i

			
			call input()

!提取分段数目及每段长度




	open(unit=39,file=fydata1,status='old')

	i=0
		do while(.true.)
			i=i+1
	   READ(39,*,end=400) length(i)

if(i==1) then
	   nbegin(i)=1
	   nend(i)=length(i)
	   endif

if(i>1) then
	   nbegin(i)=nend(i-1)+1
	   nend(i)=nend(i-1)+length(i)
	   endif

	   !	write(*,*) Qobs(i),Qsim(i)

	enddo
400 	Close (unit=39)

nsplit=i-1

write(*,*) nsplit

ndata=nend(nsplit)

write(*,*) ndata

	do i=1,nsplit
				istate(i,1)=25.0
				istate(i,2)=10.0
				istate(i,3)=1.0
			enddo

            end subroutine initial

!/////////////////////////////////////////////////////////////////////////////////
			subroutine input()
			implicit none
			real rand

			integer :: dmleap(12)=(/ 31,29,31,30,31,30,31,31,30,31,30,31 /)
			integer :: ndmleap(12)=(/ 31,28,31,30,31,30,31,31,30,31,30,31 /)
			integer i,j,ii,ny0,NX
			real cf,dis(ndatamax)


			goto 200

!///////////////////////monthly data///////////////////////////////////
			
			
			
			open(10,file="./input.txt")
			i=0
			do while(.true.)
				i=i+1
				read(10,*,end=100) rain(i),panev(i),Qobs(i)
				!write(*,*) rain(i)
				!rain(i)=rain(i)+rain(i)*(rand()*0.6-0.03)
					!write(*,*) rain(i)
				!panev(i)=panev(i)+panev(i)*(rand()*0.6-0.03)
				!Qobs(i)=Qobs(i)+Qobs(i)*(rand()*0.6-0.03)
			enddo
			100  close(10)
			ndata=i-1

			return

			cf=area/86.4      ! scaling factor
			!---------------------------------------
			! convert the discharge to monthly water depth 
			do i=1,ndata
				j=(i-1)/12
				ny0=nyear0+j
				ii=i-12*j
				IF( (MOD(ny0,4)==0 .and. MOD(ny0,100)/=0) .or. (MOD(ny0,400)==0) ) then
					  Qobs(i)=dis(i)*dmleap(ii)/cf
				else
					  Qobs(i)=dis(i)*ndmleap(ii)/cf
				endif
			enddo

!///////////////////////hourly data for XAJ///////////////////////////////////



200			open(10,file=filename)

			read(10,*) NX,area

			i=0
			do while(.true.)
				i=i+1
				read(10,*,end=300) rain(i),panev(i),Qobs(i),Tave(i)

	!			write(*,*) rain(i),panev(i),Qobs(i),Tave(i)
			enddo
300         close(10)
		!	ndata=i-1









         close(10)
	


			end subroutine input


!/////////////////////////////////////////////////////////////////////////////////
			subroutine output()
			implicit none

			integer i
			real distt

			open(10,access='append',file="./output.txt")
				open(34,file="./output1.txt")
			write(10,*)bestobj
			
			distt=0.0

			do i=1,nsplit
				write(10,'(1x,18f10.4)') bestpara(i,:)
				if(i < nsplit) then
				distt=distt+dist(bestpara(i,:),bestpara(i+1,:))
				end if
				!write(10,*) distt
			enddo

			write(34,*)bestobj,distt
			!do i=1,ndata
			!	write(34,*) rain(i),panev(i),Qobs(i),Qsim(i)
			!enddo


			close(10)

			end subroutine output
!/////////////////////////////////////////////////////////////////////////////////
			subroutine output1()
			implicit none

			integer i

			open(17,access='append',file="./output0.txt")
				open(347,access='append',file="./output01.txt")
			write(17,*)bestobj
			do i=1,nsplit
				write(17,'(1x,18f10.4)') bestpara(i,:)
			enddo

			write(347,*)bestobj
			!do i=1,ndata
			!	write(34,*) rain(i),panev(i),Qobs(i),Qsim(i)
			!enddo


			close(10)

			end subroutine output1
!/////////////////////////////////////////////////////////////////////////////////
			function dist(paraa,parab)
			implicit none

			real dist,paraa(nparameter),parab(nparameter)
			integer i

			dist=0.0
			do i=1,nparameter
				dist=dist+(abs(paraa(i)-parab(i))/(paramax(i)-paramin(i)))
!				if (abs(paraa(i)-parab(i))/(paramax(i)-paramin(i))>0.01) then
!				dist=dist+500*abs(paraa(i)-parab(i))/(paramax(i)-paramin(i))
!				end if
			enddo
			
			end function dist

!/////////////////////////////////////////////////////////////////////////////////


		end module simulation