 


 program fs

  use constants
  use fs
  real(kind = 8) :: fpp, eta_end, delta, h

  call general(fpp, eta_end)
  
! call delta_calculation(fpp, eta_end, delta)

! if(stf.eq.1.d0) then
!  h = 4.d0 * delta / dble(jmax-1)
! else        
!  h = (4.d0 * delta) * (stf - 1.d0) / ( stf**(jmax-1) - 1.d0 )
! endif

  call baseflow_fs(fpp) 

 end program fs     
