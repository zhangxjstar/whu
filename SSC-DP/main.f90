	Program DDP
	use simulation

	implicit none
	integer im,i,j,in
	real uzi,besto(200),bestnse

bestnse=-100.0


	
		open(unit=103,file=".\LIST.txt",status='old')

	do i_file=1,1262
	read(103,*) filenamex(i_file)
!	write(*,*) filenamex(i_file)
	enddo


	do i_file=1,1

		write (filename,'("../i/",(a7),".txt")') filenamex(i_file)	
write(*,*) filename
	
		write (fydata1,'("../i/",(a7),"sp.txt")') filenamex(i_file)
		write(*,*) fydata1
	
		write (fydata2,'("../i/",(a7),"par.txt")') filenamex(i_file)
		write(*,*) fydata2

		write (fydata3,'("../i/",(a7),"re.txt")') filenamex(i_file)
		write(*,*) fydata3


	call initial()


	
	call runDDP()
	call evaluate1()
	call output()



	
	do im=1,15
	
	call runDDP()	
	call evaluate1()
	call output()




besto(im)=bestobj(1)



if(besto(im)>bestnse) then

bestnse=besto(im)

bep=bestpara



endif






	end do

!Êä³öÎÄ¼ş
open(1031,file=fydata3)
write(1031,*) bestnse

do in=1,nsplit

write(1031,*) bep(in,1)

enddo





enddo 




	end program DDP