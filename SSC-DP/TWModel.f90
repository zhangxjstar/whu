subroutine TWmodel(r,e,c,sc,hw,sq,s)  !数据个数、降雨、蒸发皿蒸发、参数C、SC，流量,S为时段末土壤净含水量，即hw(i+1)
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