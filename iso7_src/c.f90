      subroutine invert_helm_ed_c_interface(d, T, e, p, abar, zbar, dp_dd, dp_dt, de_dt, failure ) 
      include 'implno.dek'
      include 'const.dek'
      include 'vector_eos.dek'

      double precision, intent(in) :: d, e, abar, zbar
      double precision, intent(inout) :: T, p
      double precision, intent(out) :: dp_dd, dp_dt, de_dt
      integer, intent(out) :: failure
   
      etot_row(1) = e / d
      den_row(1) = d
      abar_row(1) = abar
      zbar_row(1) = zbar
      temp_row(1) = T
      jlo_eos = 1 
      jhi_eos = 1
   !   write(*,*) d, T, e, p, abar, zbar
      call invert_helm_ed
      T = temp_row(1)
      dp_dd = dpd_row(1)
      dp_dt = dpt_row(1)
      de_dt = det_row(1)
      p = ptot_row(1)
      if( eosfail ) then
        failure = 1
      else 
        failure = 0
      endif
  
      

      end subroutine



      subroutine helmeos_c_interface(d, T, e, p, abar, zbar, dp_dd, dp_dt, de_dt, failure ) 
      include 'implno.dek'
      include 'const.dek'
      include 'vector_eos.dek'

      double precision, intent(in) :: d, T, abar, zbar
      double precision, intent(inout) :: e, p
      double precision, intent(out) :: dp_dd, dp_dt, de_dt
      integer, intent(out) :: failure
    
      den_row(1) = d
      abar_row(1) = abar
      zbar_row(1) = zbar
      temp_row(1) = T
      jlo_eos = 1 
      jhi_eos = 1
   !   write(*,*) d, T, e, p, abar, zbar
      call helmeos
      e = etot_row(1)*d
      p = ptot_row(1)
      dp_dd = dpd_row(1)
      dp_dt = dpt_row(1)
      de_dt = det_row(1)
    if( eosfail ) then
        failure = 1
      else 
        failure = 0
      endif
         

      end subroutine











!---------------------------------------------------------------------
      subroutine burn_dt(din,ein,tin,xin,dt)
      include 'implno.dek'
      include 'burn_common.dek'
      include 'network.dek'

     double precision nse_temp_switch
      common /nsetsw/  nse_temp_switch
	
       double precision, intent(in) :: ein, din, tin, xin(ionmax)
       double precision, intent(out) :: dt
       integer :: i
       double precision y(ionmax+3),dydx(ionmax+3), yy(ionmax), abar, zbar, wbar, ye, xcess
       call azbar(xin,aion,zion,wion,ionmax, &
                 yy,abar,zbar,wbar,ye,xcess)
	!write(*,*) xin(1), xin(2), xin(3)
        nse_temp_switch = 10.0d9
        y(1:ionmax) = yy(1:ionmax)
       y(itemp) = tin
       y(iener) = ein / din
       y(iden) = din
       do i=1,ionmax
        y(i) = min(1.0d0, max(y(i),1.0d-30))
       enddo
       call iso7(0.0d0,y,dydx)
       dt = abs(y(iener)/dydx(iener))

      end



