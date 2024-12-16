	subroutine insert(arg_n,arg_x,arg_y,arg_ax,arg_ay)
	implicit none

	integer arg_n,imax,imin,imid,i
	real arg_x(arg_n),arg_y(arg_n),arg_ax,arg_ay

	if (arg_ax >= arg_x(arg_n)) then
		arg_ay = arg_y(arg_n) + (arg_ax - arg_x(arg_n))*(arg_y(arg_n) - arg_y(arg_n - 1))/(arg_x(arg_n) - arg_x(arg_n - 1))
		return
	elseif(arg_ax <= arg_x(1)) then
		arg_ay = arg_y(1) + (arg_ax - arg_x(1))*(arg_y(1) - arg_y(2))/(arg_x(1) - arg_x(2))
		return
	end if
	
	imax=arg_n
	imin=1

10	imid = (imax+imin)/2
	if(arg_ax > arg_x(imid) )then
		imin = imid
	elseif(arg_ax < arg_x(imid))then
		imax = imid
	else
		arg_ay = arg_y(imid)
		return
	endif
	
	if(imax - imin /= 1) goto 10 

	arg_ay = arg_y(imin) + (arg_ax - arg_x(imin)) &
	  *(arg_y(imax) - arg_y(imin))/(arg_x(imax) - arg_x(imin))

    end subroutine insert


!	subroutine insert(arg_n,arg_x,arg_y,arg_ax,arg_ay)
!	implicit none
!
!	integer arg_n,i
!	real arg_x(arg_n),arg_y(arg_n),arg_ax,arg_ay
!
!	if (arg_ax >= arg_x(arg_n)) then
!		arg_ay = arg_y(arg_n) + (arg_ax - arg_x(arg_n))*(arg_y(arg_n) - arg_y(arg_n - 1))/(arg_x(arg_n) - arg_x(arg_n - 1))
!		return
!	end if
!
!	do i = 1,arg_n - 1
!		if (arg_ax < arg_x(i+1)) then
!			arg_ay = arg_y(i) + (arg_ax - arg_x(i)) &
!			  *(arg_y(i+1) - arg_y(i))/(arg_x(i+1) - arg_x(i))
!			return
!		end if
!	enddo
!
!    end subroutine insert