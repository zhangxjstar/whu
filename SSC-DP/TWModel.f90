subroutine TWmodel(r,e,c,sc,hw,sq,s)  !���ݸ��������ꡢ����������������C��SC������,SΪʱ��ĩ��������ˮ������hw(i+1)
!   use global
!   use numerical_libraries
   implicit none
 !  integer :: ndata
   real,parameter :: dt=1.0
   real :: c,sc
   real :: r,e,se,sq,hw,s
 !  integer :: ii

     

   ! simulated evaporation
        se=c*e*TANH(r/e)*dt
   ! simulated discharge
        sq=( hw+(r-se)*dt )*TANH( (hw+(r-se)*dt)/sc )
   ! simulated water depth
 !      IF(ii<ndata)  then
	     s=hw+(r-(sq+se))*dt
  !     ENDIF
	   
	   IF(s<EPSILON(1.)) s=EPSILON(1.)

end subroutine TWmodel